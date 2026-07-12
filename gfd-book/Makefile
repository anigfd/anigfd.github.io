# Build pipeline: marimo notebooks -> WASM bundles in the Hugo static tree.
NB_DIR   := notebooks
OUT_DIR  := site/static/nb
SITE_DIR := site
NOTEBOOKS := $(wildcard $(NB_DIR)/ch*.py)
# gfdlib wheel, shipped in notebooks/public/ so every WASM export carries it
# (notebooks micropip-install it when running under Pyodide).
WHEEL    := $(NB_DIR)/public/gfdlib-0.1.0-py3-none-any.whl

.PHONY: notebooks wheel serve test clean

wheel: $(WHEEL)
$(WHEEL): gfdlib/*.py pyproject.toml
	pip wheel --no-deps -q -w $(NB_DIR)/public .

notebooks: $(WHEEL) $(NOTEBOOKS)
	@for nb in $(NOTEBOOKS); do \
		name=$$(basename $$nb .py); \
		echo "exporting $$name"; \
		marimo export html-wasm $$nb -o $(OUT_DIR)/$$name --mode run || exit 1; \
		if [ -f $(NB_DIR)/data/$$name.npz ]; then \
			cp $(NB_DIR)/data/$$name.npz $(OUT_DIR)/$$name/public/; \
		fi; \
	done
	@echo "done -> $(OUT_DIR)"

serve:
	cd $(SITE_DIR) && hugo server -D

test:
	pytest -q tests/

clean:
	rm -rf $(OUT_DIR)/*
