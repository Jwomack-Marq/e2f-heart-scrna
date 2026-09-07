#!/usr/bin/env python3
"""Assert the methods book (docs/) covers the app and references only files that exist.

Three failures this catches, all of which have happened to documentation before:

  1. A tab gets added to the app and no chapter documents it. The book still renders
     and still looks complete, which is worse than an obvious gap.
  2. A chapter claims a tab that no longer exists (renamed, merged, removed), so a
     reader goes looking for a panel that isn't there.
  3. A chapter references a figure or generated table that docs/export_assets.R never
     produced. Quarto renders a broken image without failing.

Same job as tools/check_download_coverage.R does for the app's download handlers, and
written in Python for the same reason the render image has no R: this must run on a
machine with neither R nor the data bundle.

    python3 tools/check_docs_coverage.py            # exit 1 on any failure
    python3 tools/check_docs_coverage.py --list     # print what it found and exit 0
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
APP = REPO / "shiny_app" / "app.R"
DOCS = REPO / "docs"

# A top-level navbar tab is a nav_panel( at exactly two spaces of indentation inside
# page_navbar(); sub-tabs live inside navset_card_tab() at 6+ spaces. Anchoring on the
# indentation is what keeps this from claiming all ~60 nav_panels need a chapter.
TOP_TAB = re.compile(r'^  nav_panel\(\s*"([^"]+)"')
FRONT_TABS = re.compile(r"^tabs:\s*$")
LIST_ITEM = re.compile(r'^\s*-\s*"?([^"\n]+?)"?\s*$')
ASSET_REF = re.compile(r"(?:\]\(|src=\"|include\s+)((?:assets|_generated)/[A-Za-z0-9._/-]+)")


def app_tabs() -> list[str]:
    if not APP.exists():
        sys.exit(f"error: {APP} not found")
    out = []
    for line in APP.read_text(encoding="utf-8").splitlines():
        m = TOP_TAB.match(line)
        if m:
            out.append(m.group(1))
    return out


def chapter_front_matter(path: Path) -> tuple[dict, str]:
    """Return (front-matter-ish dict, body). Only the keys this check needs are parsed --
    a real YAML parser is not worth a dependency for `title` and a list of tabs."""
    text = path.read_text(encoding="utf-8")
    if not text.startswith("---"):
        return {}, text
    end = text.find("\n---", 3)
    if end == -1:
        return {}, text
    front, body = text[3:end], text[end + 4 :]
    meta: dict = {}
    in_tabs = False
    for line in front.splitlines():
        if FRONT_TABS.match(line):
            in_tabs, meta["tabs"] = True, []
            continue
        if in_tabs:
            m = LIST_ITEM.match(line)
            if m and line.startswith((" ", "-")):
                meta["tabs"].append(m.group(1).strip())
                continue
            in_tabs = False
        if ":" in line and not line.startswith((" ", "-")):
            k, v = line.split(":", 1)
            meta[k.strip()] = v.strip().strip('"')
    return meta, body


def main() -> int:
    listing = "--list" in sys.argv
    tabs = app_tabs()
    chapters = sorted(p for p in DOCS.glob("*.qmd"))
    if not chapters:
        sys.exit("error: no .qmd chapters in docs/")

    claimed: dict[str, list[str]] = {}
    problems: list[str] = []
    referenced: set[str] = set()

    for ch in chapters:
        meta, body = chapter_front_matter(ch)
        for t in meta.get("tabs", []):
            claimed.setdefault(t, []).append(ch.name)
        for ref in ASSET_REF.findall(body):
            referenced.add(ref)
            if not (DOCS / ref).exists():
                problems.append(f"{ch.name}: references {ref}, which does not exist")
        if "engine" in meta and meta["engine"] != "markdown":
            problems.append(f"{ch.name}: engine is {meta['engine']!r}, expected 'markdown' "
                            "(the render image has no R)")
        if "engine" not in meta:
            problems.append(f"{ch.name}: no `engine: markdown` in the front matter")

    for t in tabs:
        who = claimed.get(t, [])
        if not who:
            problems.append(f"app tab {t!r} is not documented by any chapter")
        elif len(who) > 1:
            problems.append(f"app tab {t!r} is claimed by {', '.join(who)} -- pick one")
    for t, who in claimed.items():
        if t not in tabs:
            problems.append(f"{who[0]}: claims tab {t!r}, which is not a top-level tab in app.R")

    if listing:
        print(f"{len(tabs)} top-level app tabs, {len(chapters)} chapters, "
              f"{len(referenced)} asset references")
        for t in tabs:
            print(f"  {'OK ' if t in claimed else 'GAP'}  {t:<34} {', '.join(claimed.get(t, []))}")

    on_disk = {f"assets/{p.name}" for p in (DOCS / "assets").glob("*")}
    on_disk |= {f"_generated/{p.name}" for p in (DOCS / "_generated").glob("*.md")}
    unused = sorted(on_disk - referenced - {"_generated/_manifest.md"})
    if unused:
        print(f"note: {len(unused)} exported file(s) no chapter references: "
              f"{', '.join(unused[:6])}{' ...' if len(unused) > 6 else ''}")

    if problems:
        print(f"\n{len(problems)} problem(s):", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1
    print(f"ok: {len(tabs)} app tabs all documented; "
          f"{len(referenced)} referenced assets all present")
    return 0


if __name__ == "__main__":
    sys.exit(main())
