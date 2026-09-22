# Build pipeline: marimo notebooks -> ONE shared WASM bundle in the Hugo static tree.
NB_DIR    := notebooks
OUT_DIR   := site/static/nb
IPYNB_DIR := site/static/ipynb
SITE_DIR  := site
BASE_URL  := https://anigfd.github.io
NOTEBOOKS := $(wildcard $(NB_DIR)/ch*.py)
# gfdlib wheel, shipped in notebooks/public/ so every WASM export carries it
# (notebooks micropip-install it when running under Pyodide).
WHEEL     := $(NB_DIR)/public/gfdlib-0.1.0-py3-none-any.whl

# Chapters that also get an editable export, for the cards on aneeshcs.com/gfd/.
# A named list rather than every chapter: --mode edit re-executes the notebook, so
# doing all 24 would roughly double every deploy for a feature two pages use.
# Adding a chapter here is all it takes to showcase it.
SHOWCASE  := ch06_geostrophic-adjustment ch16_baroclinic-instability

.PHONY: notebooks nb-one ipynb wheel serve test clean

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
	@for name in $(SHOWCASE); do \
		echo "exporting $$name (editable)"; \
		marimo export html-wasm $(NB_DIR)/$$name.py -o $(OUT_DIR)/$${name}_edit.html \
			--mode edit || exit 1; \
	done
	@rm -f $(OUT_DIR)/CLAUDE.md
	@python3 scripts/patch_chunk_reload.py $(OUT_DIR)/*.html
	@$(MAKE) --no-print-directory ipynb
	@echo "done -> $(OUT_DIR)"

# Downloadable Jupyter versions of every chapter, linked from each chapter page.
# These are marimo's own `export ipynb` plus a dependency cell and a header --
# see scripts/export_ipynb.py. Cheap (no execution), so it runs as part of
# `notebooks`.
ipynb:
	@mkdir -p $(IPYNB_DIR)
	@python3 scripts/export_ipynb.py $(NOTEBOOKS) \
		--outdir $(IPYNB_DIR) --site-root $(SITE_DIR) --base-url $(BASE_URL)

# Iterate on one chapter without a full re-export:  make nb-one NB=ch16_baroclinic-instability
nb-one: $(WHEEL)
	@test -n "$(NB)" || { echo "usage: make nb-one NB=ch16_baroclinic-instability"; exit 1; }
	marimo export html-wasm $(NB_DIR)/$(NB).py -o $(OUT_DIR)/$(NB).html --mode run
	@case " $(SHOWCASE) " in *" $(NB) "*) \
		marimo export html-wasm $(NB_DIR)/$(NB).py -o $(OUT_DIR)/$(NB)_edit.html --mode edit ;; \
	esac
	@rm -f $(OUT_DIR)/CLAUDE.md
	@python3 scripts/patch_chunk_reload.py $(OUT_DIR)/$(NB)*.html
	@python3 scripts/export_ipynb.py $(NB_DIR)/$(NB).py \
		--outdir $(IPYNB_DIR) --site-root $(SITE_DIR) --base-url $(BASE_URL)

serve:
	cd $(SITE_DIR) && hugo server -D

test:
	pytest -q tests/

clean:
	rm -rf $(OUT_DIR) $(IPYNB_DIR)
