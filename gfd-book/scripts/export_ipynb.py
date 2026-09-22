#!/usr/bin/env python3
"""Export marimo notebooks as downloadable Jupyter notebooks.

Wraps `marimo export ipynb` -- marimo's own supported conversion -- and applies
the edits a *downloaded* notebook needs and the in-repo one does not:

1. **A dependency cell.** This book pins through `requirements.txt` and a `gfdlib`
   wheel rather than a PEP 723 header, so the install list is assembled from its
   own declarations: `marimo` exactly as pinned in `requirements.txt`, and the
   science libraries from `[project] dependencies` in `pyproject.toml`. That
   becomes a `%pip install` cell, letting the download run from a bare kernel.

2. **Bootstrap surgery.** A chapter resolves `gfdlib` either by micropip-installing
   a bundled wheel (in the browser) or by adding the repo root to `sys.path`
   (locally). Neither works for a file someone downloaded, so that block is removed
   and `gfdlib` is pip-installed from the public repository instead.

3. **Public data assets over HTTPS.** A chapter that ships a precomputed field reads
   it from `mo.notebook_location() / "public" / ...`, which resolves against the
   served bundle and means nothing on a laptop. Those paths are rewritten to
   `_gfd_public(...)`, a helper that downloads the same file from the published
   site. Currently only chapter 16 needs it; the rewrite is general, so a future
   chapter that ships data gets a working download without touching this script.

A header markdown cell is prepended explaining where the notebook came from and
the one thing that genuinely differs from the live version: marimo's sliders are
inert in Jupyter, so a parameter is changed by editing `value=` and re-running.

Usage:
    python3 scripts/export_ipynb.py notebooks/chNN.py -o site/static/ipynb/chNN.ipynb
    python3 scripts/export_ipynb.py notebooks/*.py --outdir site/static/ipynb
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path

#: Installed from the public repo, because a downloaded notebook has no checkout.
GFDLIB_SPEC = (
    '"gfdlib @ git+https://github.com/anigfd/anigfd.github.io'
    '#subdirectory=gfd-book"'
)

#: Where the published site serves the files in `notebooks/public/`.
PUBLIC_BASE_URL = "https://anigfd.github.io/nb/public"

#: `mo.notebook_location() / "public" / "x.npz"` and the notebook_dir() variant.
PUBLIC_ASSET_RE = re.compile(
    r"mo\.notebook_(?:location|dir)\(\)\s*/\s*[\"']public[\"']\s*/\s*[\"']([^\"']+)[\"']"
)


def read_dependencies(pyproject: Path, requirements: Path) -> list[str]:
    """What the download must install, from the book's own declarations.

    Two sources, mirroring how the book pins: `marimo` exactly out of
    requirements.txt -- the notebook is still a marimo notebook, its `mo.md(...)`
    cells need it, and that pin is what fixes the published runtime -- and the
    science libraries out of `[project] dependencies` in pyproject.toml, left
    unpinned so pip resolves a build that suits the reader's Python.
    """
    deps: list[str] = []

    if requirements.is_file():
        for line in requirements.read_text().splitlines():
            spec = line.split("#", 1)[0].strip()
            if re.match(r"marimo\s*[=<>~]", spec):
                deps.append(spec)
                break

    if pyproject.is_file():
        data = tomllib.loads(pyproject.read_text())
        deps += list(data.get("project", {}).get("dependencies", []))

    return deps


def _indent_width(line: str) -> int:
    return len(line) - len(line.lstrip())


def _reindent(block: list[str], indent: str) -> list[str]:
    """Dedent a captured suite to its own left margin, then set `indent`."""
    meaningful = [line for line in block if line.strip()]
    if not meaningful:
        return []
    base = min(_indent_width(line) for line in meaningful)
    return [indent + line[base:] if line.strip() else "" for line in block]


def desugar_emscripten(source: str) -> tuple[str, bool, bool]:
    """Resolve a cell's `sys.platform == "emscripten"` branch for a plain kernel.

    Returns ``(new_source, dropped_bootstrap, kept_else)``.

    Two shapes occur, and they need opposite treatment:

    * An **import cell**, whose `else` only manipulates `sys.path`. Both branches
      are repo-specific, so the whole if/else goes and the surviving imports
      (marimo, numpy, gfdlib, ...) are left exactly as the chapter wrote them.
    * A **data cell**, whose `else` does the real work of loading a field. Here the
      emscripten branch is discarded and the `else` body is kept, dedented into
      place. Deleting it -- as a bootstrap-only rule would -- silently drops the
      chapter's data.

    Line-based rather than AST-based on purpose: the cell must come out looking
    like hand-written code, and an AST round-trip would discard the comments that
    explain it.
    """
    if "emscripten" not in source:
        return source, False, False

    lines = source.splitlines()
    out: list[str] = []
    dropped_bootstrap = False
    kept_else = False
    i = 0
    while i < len(lines):
        line = lines[i]

        # The `import sys` (or `import sys as _sys`) exists only to serve the branch.
        if re.fullmatch(r"\s*import\s+sys(\s+as\s+_\w+)?\s*", line):
            i += 1
            continue

        match = re.match(
            r"(\s*)if\s+_?\w*sys\.platform\s*==\s*[\"']emscripten[\"']\s*:", line
        )
        if not match:
            out.append(line)
            i += 1
            continue

        indent = match.group(1)
        i += 1
        # Skip the emscripten suite: everything blank or indented past `if`.
        while i < len(lines) and (
            not lines[i].strip() or _indent_width(lines[i]) > len(indent)
        ):
            i += 1

        else_body: list[str] = []
        if i < len(lines) and re.match(rf"{indent}else\s*:\s*$", lines[i]):
            i += 1
            while i < len(lines) and (
                not lines[i].strip() or _indent_width(lines[i]) > len(indent)
            ):
                else_body.append(lines[i])
                i += 1

        meaningful = [
            l for l in else_body if l.strip() and not l.strip().startswith("#")
        ]
        if meaningful and all("sys.path" in l for l in meaningful):
            dropped_bootstrap = True  # import cell: both branches are repo-specific
            continue

        out.extend(_reindent(else_body, indent))
        kept_else = True

    # Collapse the run of blank lines the removal leaves behind.
    cleaned: list[str] = []
    for line in out:
        if not line.strip() and cleaned and not cleaned[-1].strip():
            continue
        cleaned.append(line)
    return "\n".join(cleaned).strip() + "\n", dropped_bootstrap, kept_else


def rewrite_public_assets(source: str) -> tuple[str, list[str]]:
    """Point `.../public/x` reads at the published site instead of the bundle."""
    names = PUBLIC_ASSET_RE.findall(source)
    if not names:
        return source, []
    rewritten = PUBLIC_ASSET_RE.sub(lambda m: f'_gfd_public("{m.group(1)}")', source)
    return rewritten, names


def strip_marimo_stop(source: str) -> bool | str:
    """Remove `mo.stop(...)` gates so a downloaded notebook computes eagerly.

    Chapters gate their expensive cells behind a Run button:
    `mo.stop(not run_btn.value, mo.callout(...))`. Under marimo that halts the
    cell until the reader clicks. Under Jupyter the button is a static widget
    whose value is False, so the gate always fires -- and `mo.stop` raises
    `MarimoStopError` outside marimo's runtime, so the download dies on the first
    gated cell instead of merely skipping it.

    Dropping the call makes the cell run when the reader runs it, which is what a
    downloaded notebook should do. The gated cell is also the cell that computes
    the result, so removing the gate is what lets the cells below it find one.

    Returns the new source, or False if the cell had no gate.
    """
    if "mo.stop(" not in source:
        return False

    out = source
    while (start := out.find("mo.stop(")) != -1:
        i = start + len("mo.stop(")
        depth, quote = 1, ""
        while i < len(out) and depth:
            ch = out[i]
            if quote:
                if ch == "\\":
                    i += 1
                elif ch == quote:
                    quote = ""
            elif ch in "\"'":
                quote = ch
            elif ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            i += 1
        if depth:                      # unbalanced: leave the cell alone
            return False
        out = out[:start] + out[i:]

    # Drop the lines the removal blanked out, and collapse the gap it leaves.
    lines = [line for line in out.splitlines() if line.strip()]
    return "\n".join(lines) + "\n"


def hoist_marimo_import(cells: list[dict]) -> bool:
    """Move the `import marimo as mo` cell to the front of the notebook.

    marimo keeps that cell at the *bottom* of the .py -- its own scheduler is a
    dependency graph, so file order is irrelevant there. Jupyter runs cells top to
    bottom, so left where it is, every `mo.md(...)` above it raises NameError and
    the whole download fails on the first markdown cell. Reordering is safe for
    the same reason it was safe to put it last: nothing else in the cell.
    """
    for index, cell in enumerate(cells):
        if cell["cell_type"] != "code":
            continue
        body = [
            line.strip()
            for line in "".join(cell["source"]).splitlines()
            if line.strip() and not line.strip().startswith("#")
        ]
        if body and all(line in ("import marimo", "import marimo as mo") for line in body):
            cells.insert(0, cells.pop(index))
            return True
    return False


def build_install_cell(
    deps: list[str], needs_gfdlib: bool, needs_public: bool
) -> list[str]:
    """The `%pip install` cell placed at the top of the download."""
    lines = [
        "# Dependencies for this notebook. Safe to skip if you already have them.\n",
        "# marimo is pinned to the version the published chapter runs on; the science\n",
        "# libraries are left to pip, matching how the book itself declares them.\n",
    ]
    if deps:
        lines.append("%pip install -q " + " ".join(f'"{d}"' for d in deps) + "\n")
    if needs_gfdlib:
        lines += [
            "\n",
            "# gfdlib holds the book's shared numerics. Installed from the public\n",
            "# repository, since a downloaded notebook has no checkout to import from.\n",
            f"%pip install -q {GFDLIB_SPEC}\n",
        ]
    if needs_public:
        lines += [
            "\n",
            "# This chapter reads a precomputed field that the published version serves\n",
            "# alongside itself. Fetch it from there, and cache it so a re-run is offline.\n",
            "def _gfd_public(name):\n",
            "    import pathlib, tempfile, urllib.request\n",
            "    dest = pathlib.Path(tempfile.gettempdir()) / name\n",
            "    if not dest.exists():\n",
            f'        urllib.request.urlretrieve(f"{PUBLIC_BASE_URL}/{{name}}", dest)\n',
            "    return dest\n",
        ]
    return lines


def build_header_cell(
    title: str, source_name: str, page_url: str | None, had_gates: bool
) -> list[str]:
    """Provenance and the one behavioural difference worth warning about."""
    lines = [
        f"# {title}\n",
        "\n",
        "*A Jupyter conversion of an interactive [marimo](https://marimo.io) notebook.*\n",
        "\n",
    ]
    if page_url:
        lines += [f"Run it live, in the browser, at <{page_url}>.\n", "\n"]
    lines += [
        "**One thing behaves differently here.** In the live version the parameters are\n",
        "sliders you drag, and every figure redraws as you move them. Jupyter cannot drive\n",
        "marimo's widgets, so they render as *static* controls showing their default\n",
        "values. To change a parameter, edit the `value=` argument where the control is\n",
        "defined and re-run the cells below it.\n",
        "\n",
        "Cells are in top-down reading order, matching the published chapter. The numerics\n",
        "are unchanged.\n",
    ]
    if had_gates:
        lines += [
            "\n",
            "The published chapter waits for a **Run** button before its expensive cells\n",
            "execute. That gate is removed here, so a cell computes when you run it.\n",
        ]
    lines += [
        "\n",
        f"Generated from `{source_name}` by `scripts/export_ipynb.py`.\n",
    ]
    return lines


def notebook_title(src: Path) -> str:
    """The chapter's own H1, so the download is titled like the published page.

    Falls back to a prettified filename. Reading the heading out of the source
    beats deriving it from the stem, which yields things like "Ch06 Geostrophic".
    """
    text = src.read_text()
    # Search only after the first `mo.md(`, so ordinary Python comments cannot be
    # mistaken for the title. The first markdown H1 after that point is the heading.
    start = text.find("mo.md(")
    if start != -1:
        for match in re.finditer(r"^\s*#\s+(\S.*?)\s*$", text[start:], re.M):
            heading = match.group(1).strip()
            if heading.startswith(("/", "-", "=", "#")):
                continue
            return heading
    return src.stem.replace("_", " ").replace("-", " ").title()


def chapter_url(src: Path, site_root: Path, base_url: str | None) -> str | None:
    """Locate the chapter's page in the Hugo content tree and build its URL.

    The part a chapter lives in is recorded only by its directory, so it is read
    from the filesystem rather than assumed -- hard-coding one base URL gets the
    other parts wrong.
    """
    if not base_url:
        return None
    matches = sorted(site_root.glob(f"content/part*/{src.stem}.md"))
    if not matches:
        return None
    part = matches[0].parent.name
    return f"{base_url.rstrip('/')}/{part}/{src.stem}/"


def convert(src: Path, dest: Path, page_url: str | None, deps: list[str]) -> None:
    """Export one marimo notebook to a downloadable .ipynb."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        raw = Path(tmp) / "raw.ipynb"
        # --sort top-down keeps the chapter's reading order; the default
        # topological order rearranges cells and makes the prose jump around.
        subprocess.run(
            [
                sys.executable, "-m", "marimo", "export", "ipynb",
                "--sort", "top-down", str(src), "-o", str(raw),
            ],
            check=True,
            capture_output=True,
        )
        nb = json.loads(raw.read_text())

    needs_gfdlib = False
    had_gates = False
    assets: list[str] = []
    for cell in nb["cells"]:
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        source, dropped, _kept = desugar_emscripten(source)
        source, names = rewrite_public_assets(source)
        ungated = strip_marimo_stop(source)
        if ungated is not False:
            source, gates = ungated, True
        else:
            gates = False
        assets += names
        if dropped or names or gates:
            cell["source"] = source.splitlines(keepends=True)
        needs_gfdlib = needs_gfdlib or dropped
        had_gates = had_gates or gates

    hoist_marimo_import(nb["cells"])

    title = notebook_title(src)
    nb["cells"] = [
        {"cell_type": "markdown", "id": "generated-header", "metadata": {},
         "source": build_header_cell(title, src.name, page_url, had_gates)},
        {"cell_type": "code", "id": "generated-install", "execution_count": None,
         "metadata": {}, "outputs": [],
         "source": build_install_cell(deps, needs_gfdlib, bool(assets))},
    ] + nb["cells"]

    # nbformat >= 4.5 requires a cell id; marimo's export omits them, which makes
    # every downstream tool emit a MissingIDFieldWarning.
    for index, cell in enumerate(nb["cells"]):
        cell.setdefault("id", f"cell-{index:03d}")

    dest.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")
    extra = f", fetches {', '.join(sorted(set(assets)))}" if assets else ""
    print(f"  {src.name} -> {dest}  ({len(nb['cells'])} cells{extra})")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument("--outdir", type=Path)
    parser.add_argument(
        "--base-url",
        help="Site root for the 'read it live' link, e.g. https://anigfd.github.io",
    )
    parser.add_argument(
        "--site-root", type=Path, default=Path("site"),
        help="Hugo site directory, used to find which part each chapter lives in",
    )
    parser.add_argument(
        "--pyproject", type=Path, default=Path("pyproject.toml"),
        help="Where to read the science-library dependencies from",
    )
    parser.add_argument(
        "--requirements", type=Path, default=Path("requirements.txt"),
        help="Where to read the exact marimo pin from",
    )
    args = parser.parse_args()

    if args.output and len(args.sources) != 1:
        parser.error("-o takes exactly one source; use --outdir for several")
    if not args.output and not args.outdir:
        parser.error("give -o or --outdir")

    deps = read_dependencies(args.pyproject, args.requirements)
    for src in args.sources:
        dest = args.output or (args.outdir / f"{src.stem}.ipynb")
        convert(src, dest, chapter_url(src, args.site_root, args.base_url), deps)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
