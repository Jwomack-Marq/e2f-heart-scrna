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
#
# Most aes() here name data columns, so their environments are dead weight. But
# not all: commun_heat_gg() builds aes(text = paste0(..., metric, ...)) where
# `metric` is one of its ARGUMENTS. Blindly rebinding that quosure to globalenv
# makes `metric` unresolvable, the verification build fails, and the old code
# then handed back the untouched plot -- 1.26 GB with the whole app hanging off
# it, which saveRDS would then spend minutes gzipping. So instead of discarding
# each quosure's environment, copy the few small values it actually needs into a
# fresh one whose parent is globalenv.
#
# Returns list(plot, ok, why). ok = FALSE means "do not hand this off": an
# unstripped plot is both enormous and useless at the far end, since it carries
# environments the studio cannot reconstruct.
strip_plot_env <- function(p) {
  # plotly-only aesthetics never render in a static export -- ggplot2 ignores
  # `text` with a warning -- and they are exactly where the awkward closures
  # live. Dropping them removes the problem rather than working around it.
  PLOTLY_AES <- c("text", "customdata", "key", "frame", "ids", "hovertext")

  # Values a quosure needs that are NOT data columns: copy them if they are
  # small and self-contained. Anything big or exotic (a matrix, a Seurat object,
  # an environment) is left behind deliberately -- carrying it is the bug.
  small_enough <- function(v)
    (is.atomic(v) || is.factor(v)) && !is.matrix(v) && length(v) <= 1000
  slim_quo <- function(q, cols) {
    env <- rlang::quo_get_env(q)
    if (identical(env, globalenv()) || rlang::is_namespace(env)) return(q)
    need <- setdiff(all.vars(rlang::quo_get_expr(q)), cols)
    vals <- list()
    for (v in need) {
      val <- tryCatch(get(v, envir = env), error = function(e) NULL)
      # functions resolve through the search path at build time; don't copy them
      if (!is.null(val) && !is.function(val) && small_enough(val)) vals[[v]] <- val
    }
    rlang::quo_set_env(q, rlang::new_environment(vals, parent = globalenv()))
  }
  fixq <- function(m, cols) {
    m <- m[setdiff(names(m), PLOTLY_AES)]
    for (nm in names(m)) if (rlang::is_quosure(m[[nm]])) m[[nm]] <- slim_quo(m[[nm]], cols)
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
    cols <- unique(c(names(p$data),
                     unlist(lapply(p$layers, function(l) names(l$data)))))
    q$plot_env <- rlang::new_environment(parent = globalenv())
    q$mapping  <- fixq(q$mapping, cols)
    for (l in q$layers) { l$mapping <- fixq(l$mapping, cols); slim_ggproto(l) }
    if (!is.null(old_facet)) {
      fp <- q$facet$params
      for (s in intersect(c("facets", "rows", "cols"), names(fp)))
        fp[[s]] <- fixq(fp[[s]], cols)
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
    is.list(suppressWarnings(ggplot2::ggplot_build(q))$data)
  }, error = function(e) conditionMessage(e))
  if (isTRUE(ok)) return(list(plot = q, ok = TRUE, why = NULL))
  for (i in seq_along(p$layers)) p$layers[[i]]$mapping <- old_layer[[i]]
  if (!is.null(old_facet)) p$facet$params <- old_facet
  list(plot = p, ok = FALSE,
       why = if (is.character(ok)) ok else "the detached plot would not rebuild")
}

# Serialized size WITHOUT touching the disk. The point is to refuse an oversized
# figure before saveRDS spends minutes gzipping it -- the failure mode that read
# as "the app crashed": the old code wrote first and measured afterwards.
plot_nbytes <- function(x) tryCatch(length(serialize(x, NULL)), error = function(e) NA_real_)

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
    g <- function(s) input[[paste0(opts_prefix, "_", s)]]
    st <- strip_plot_env(apply_fig_opts(plot_reactive(), opts_prefix, input))
    # Refuse rather than ship an attached plot: it would carry the app's data by
    # value (gigabytes) and still not rebuild at the far end.
    if (!isTRUE(st$ok))
      stop(sprintf("this figure could not be detached from the app (%s)", st$why))
    # Measure in memory, BEFORE writing. saveRDS on an oversized object is the
    # slow, session-killing step, so it must never be how we find out.
    sz <- plot_nbytes(st$plot)
    if (is.na(sz) || sz > STUDIO_REFUSE_BYTES)
      stop(sprintf("this figure is too large to hand off (%s)",
                   if (is.na(sz)) "size could not be measured" else sprintf("%.0f MB", sz / 1e6)))
    if (sz > STUDIO_WARN_BYTES)
      studio_notify(sprintf("Large figure (%.0f MB) — the studio may take a moment to open it.",
                            sz / 1e6))
    spec <- list(
      version = 1L,
      plot    = st$plot,
      meta    = list(app = "e2f-heart-scrna", prefix = prefix,
                     palette = g("palette"),
                     w_in = g("w") %||% 7, h_in = g("h") %||% 5, dpi = g("dpi") %||% 300,
                     ggplot2 = as.character(utils::packageVersion("ggplot2")),
                     created = Sys.time()))
    tok <- studio_token()
    tmp <- file.path(HANDOFF_DIR, paste0(tok, ".rds.tmp"))
    saveRDS(spec, tmp)
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
