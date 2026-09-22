# Build pipeline: marimo notebooks -> ONE shared WASM bundle in the Hugo static tree.
NB_DIR   := notebooks
OUT_DIR  := site/static/nb
SITE_DIR := site
NOTEBOOKS := $(wildcard $(NB_DIR)/ch*.py)
# gfdlib wheel, shipped in notebooks/public/ so every WASM export carries it
# (notebooks micropip-install it when running under Pyodide).
WHEEL    := $(NB_DIR)/public/gfdlib-0.1.0-py3-none-any.whl

.PHONY: notebooks nb-one wheel serve test clean

wheel: $(WHEEL)
$(WHEEL): gfdlib/*.py pyproject.toml
	pip wheel --no-deps -q -w $(NB_DIR)/public .

# Every chapter exports into the SAME directory, so marimo emits one shared
# assets/ and one shared public/ (holding the wheel and ch16's precomputed
# life-cycle) rather than a ~27 MB self-contained copy per chapter. The
# alternative form, `-o $(OUT_DIR)/$$name/`, costs ~650 MB across 24 chapters
# and makes the reader re-download Pyodide on every chapter.
#
# `clean` runs first because marimo merges into an existing assets/ rather than
# replacing it (shutil.copytree(..., dirs_exist_ok=True)), so chunks from an
# older marimo version would otherwise accumulate forever.
notebooks: clean $(WHEEL)
	@for nb in $(NOTEBOOKS); do \
		name=$$(basename $$nb .py); \
		echo "exporting $$name"; \
		marimo export html-wasm $$nb -o $(OUT_DIR)/$$name.html --mode run || exit 1; \
	done
	@rm -f $(OUT_DIR)/CLAUDE.md
	@python3 scripts/patch_chunk_reload.py $(OUT_DIR)/*.html
	@echo "done -> $(OUT_DIR)"

# Iterate on one chapter without a full re-export:  make nb-one NB=ch16_baroclinic-instability
nb-one: $(WHEEL)
	@test -n "$(NB)" || { echo "usage: make nb-one NB=ch16_baroclinic-instability"; exit 1; }
	marimo export html-wasm $(NB_DIR)/$(NB).py -o $(OUT_DIR)/$(NB).html --mode run
	@rm -f $(OUT_DIR)/CLAUDE.md
	@python3 scripts/patch_chunk_reload.py $(OUT_DIR)/$(NB).html

serve:
	cd $(SITE_DIR) && hugo server -D

test:
	pytest -q tests/

clean:
	rm -rf $(OUT_DIR)
