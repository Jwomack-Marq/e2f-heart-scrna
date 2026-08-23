# ---- Figure Studio handoff ---------------------------------------------------
# "Figure Studio" is a separate companion app (own repo/container) for
# publication-grade restyling: exact journal column widths, per-element fonts,
# per-level colours, TIFF export. Rather than duplicate ~30 plot builders there,
# this app hands over the *built ggplot object*: the plot reactive's result is
# environment-stripped, wrapped in a versioned "figspec" list, saved as
# HANDOFF_DIR/<token>.rds on a volume both containers mount, and the studio is
# opened at FIGURE_STUDIO_BASE?fig=<token> in a new tab. The figspec contract
# (fields, TTL, token grammar, ggplot2 version pin) is documented in the
# figure-studio repo's README, which is the contract of record.
#
# Everything here is inert unless both env vars point somewhere real, so local
# runs, rsconnect, and the test suite see exactly the app they saw before.

STUDIO_BASE <- Sys.getenv("FIGURE_STUDIO_BASE", "")
if (nzchar(STUDIO_BASE)) STUDIO_BASE <- sub("/*$", "/", STUDIO_BASE)
HANDOFF_DIR <- Sys.getenv("HANDOFF_DIR", "")
STUDIO_ON   <- nzchar(STUDIO_BASE) && nzchar(HANDOFF_DIR) && dir.exists(HANDOFF_DIR)

STUDIO_WARN_BYTES   <- 25e6   # notify: unusually fat figure, studio will be slow
STUDIO_REFUSE_BYTES <- 100e6  # refuse: something (an unstripped env) went wrong
STUDIO_TTL_SECS     <- 24 * 3600

# A ggplot built inside a reactive references the reactive's frame -- and through
# it the whole server frame (input/output/session, every reactive) -- via
# plot_env and the quosures in each mapping. saveRDS would serialize all of that
# by value (only globalenv and package namespaces are stored by reference).
# Every aes() in this app uses literal column names on small frames embedded in
# the plot, so those environments are dead weight: rebind them to globalenv and
# prove the plot still builds. Layers and facets are ggproto (environments), so
# edits to them touch the caller's object too -- on a failed build we restore
# the originals and hand back the plot untouched; the size guard downstream is
# the backstop.
strip_plot_env <- function(p) {
  fixq <- function(m) {
    for (nm in names(m)) if (rlang::is_quosure(m[[nm]]))
      m[[nm]] <- rlang::quo_set_env(m[[nm]], globalenv())
    m
  }
  # A ggproto INSTANCE built inside user code (a layer, a facet, an added
  # scale) carries a `super` closure whose environment holds the constructor
  # call's promises -- and through them the whole calling frame, quosure fix
  # or not. Swap it for an equivalent closure over the already-resolved parent:
  # dispatch behaves identically, the frames drop. In place on purpose.
  slim_ggproto <- function(obj) {
    ok <- inherits(obj, "ggproto") && is.environment(obj) &&
      exists("super", envir = obj, inherits = FALSE)
    if (!ok) return(invisible(NULL))
    parent <- tryCatch(get("super", envir = obj)(), error = function(e) NULL)
    if (inherits(parent, "ggproto")) {
      f <- function() parent
      environment(f) <- list2env(list(parent = parent), parent = emptyenv())
      assign("super", f, envir = obj)
    }
    invisible(NULL)
  }
  old_layer  <- lapply(p$layers, function(l) l$mapping)
  old_facet  <- tryCatch(p$facet$params, error = function(e) NULL)
  q <- p                                    # list copy; layers/facet stay shared
  ok <- tryCatch({
    q$plot_env <- rlang::new_environment(parent = globalenv())
    q$mapping  <- fixq(q$mapping)
    for (l in q$layers) { l$mapping <- fixq(l$mapping); slim_ggproto(l) }
    if (!is.null(old_facet)) {
      fp <- q$facet$params
      for (s in intersect(c("facets", "rows", "cols"), names(fp)))
        fp[[s]] <- fixq(fp[[s]])
      q$facet$params <- fp                                    # in place (ggproto)
    }
    # sweep every ggproto component the plot carries (facet, coordinates,
    # layout, scales, guides, ...): ggplot2 4.x stores them as S7 properties
    # (= attributes), 3.x as list slots. One unslimmed instance is enough to
    # keep the constructor frames alive, so hit them all, one level deep.
    slim_all <- function(x) {
      slim_ggproto(x)
      if (is.list(x)) for (el in x) slim_ggproto(el)
      for (fld in c("scales", "guides")) {
        sub <- tryCatch(x[[fld]], error = function(e) NULL)
        if (is.list(sub)) for (el in sub) slim_ggproto(el)
      }
    }
    for (a in c(attributes(q), if (is.list(q)) unclass(q))) slim_all(a)
    is.list(ggplot2::ggplot_build(q)$data)
  }, error = function(e) FALSE)
  if (ok) return(q)
  for (i in seq_along(p$layers)) p$layers[[i]]$mapping <- old_layer[[i]]
  if (!is.null(old_facet)) p$facet$params <- old_facet
  p
}

# showNotification needs a real session; under shiny::testServer (or any
# headless call) fall back to message() so the tryCatch in studio_handoff
# reports the real failure instead of the notifier's.
studio_notify <- function(msg, type = "warning")
  tryCatch(showNotification(msg, type = type), error = function(e) message(msg))

# Timestamp prefix so `ls` of the handoff dir reads chronologically; random
# suffix so a token can't be guessed from another. Auth itself is the lab
# login + same-origin, not the token.
studio_token <- function()
  paste0(format(Sys.time(), "%Y%m%d%H%M%S"),
         paste(sample(c(letters, 0:9), 12, replace = TRUE), collapse = ""))

# Handoffs are throwaway by design (TTL 24 h). Called on every write here and on
# every studio session start, so neither app needs a cron.
prune_handoff <- function(dir = HANDOFF_DIR, max_age = STUDIO_TTL_SECS) {
  fs  <- list.files(dir, pattern = "\\.rds(\\.tmp)?$", full.names = TRUE)
  old <- fs[difftime(Sys.time(), file.mtime(fs), units = "secs") > max_age]
  unlink(old)
  invisible(length(old))
}

# Build -> apply the in-app figure options -> strip -> figspec -> atomic write.
# Returns the token, or NULL after notifying the user. The .tmp + rename dance
# keeps the studio from ever reading a half-written file off the shared volume.
studio_handoff <- function(prefix, plot_reactive, input, opts_prefix = prefix) {
  tryCatch({
    g   <- function(s) input[[paste0(opts_prefix, "_", s)]]
    spec <- list(
      version = 1L,
      plot    = strip_plot_env(apply_fig_opts(plot_reactive(), opts_prefix, input)),
      meta    = list(app = "e2f-heart-scrna", prefix = prefix,
                     palette = g("palette"),
                     w_in = g("w") %||% 7, h_in = g("h") %||% 5, dpi = g("dpi") %||% 300,
                     ggplot2 = as.character(utils::packageVersion("ggplot2")),
                     created = Sys.time()))
    tok <- studio_token()
    tmp <- file.path(HANDOFF_DIR, paste0(tok, ".rds.tmp"))
    saveRDS(spec, tmp)
    sz <- file.size(tmp)
    if (sz > STUDIO_REFUSE_BYTES) {
      unlink(tmp)
      stop(sprintf("this figure serializes to %.0f MB, which means environment stripping failed",
                   sz / 1e6))
    }
    if (sz > STUDIO_WARN_BYTES)
      studio_notify(sprintf("Large figure (%.0f MB) — the studio may take a moment to open it.",
                            sz / 1e6))
    file.rename(tmp, file.path(HANDOFF_DIR, paste0(tok, ".rds")))
    prune_handoff()
    tok
  }, error = function(e) {
    studio_notify(paste("Couldn't send this figure to Figure Studio:", conditionMessage(e)),
                  type = "error")
    NULL
  })
}

# The button only exists when the studio is configured, so an unconfigured app
# renders byte-identical UI to before this feature existed.
studio_btn <- function(prefix) {
  if (!STUDIO_ON) return(NULL)
  actionButton(paste0(prefix, "_studio"), "Figure Studio",
               class = "btn-sm btn-outline-primary",
               title = "Open this figure in the Figure Studio editor (new tab)")
}

# Server-side twin of studio_btn; register_fig calls it so every registered
# figure gets the handoff without per-tab wiring.
register_studio <- function(prefix, plot_reactive, input, opts_prefix = prefix,
                            session = shiny::getDefaultReactiveDomain()) {
  if (!STUDIO_ON) return(invisible(NULL))
  observeEvent(input[[paste0(prefix, "_studio")]], {
    tok <- studio_handoff(prefix, plot_reactive, input, opts_prefix)
    if (!is.null(tok))
      session$sendCustomMessage("open_studio", paste0(STUDIO_BASE, "?fig=", tok))
  }, ignoreInit = TRUE)
  invisible(NULL)
}

# window.open from an async custom-message handler is what pop-up blockers
# block; the fallback input drives a modal with a plain link (see server).
studio_js <- function() tags$script(HTML(
  "Shiny.addCustomMessageHandler('open_studio', function(url) {
     var w = window.open(url, '_blank');
     if (!w) Shiny.setInputValue('studio_popup_blocked', url, {priority: 'event'});
   });"))
