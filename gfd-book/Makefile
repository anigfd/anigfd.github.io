# Build pipeline: marimo notebooks -> WASM bundles in the Hugo static tree.
NB_DIR   := notebooks
OUT_DIR  := site/static/nb
SITE_DIR := site
NOTEBOOKS := $(wildcard $(NB_DIR)/ch*.py)

.PHONY: notebooks serve test clean
notebooks: $(NOTEBOOKS)
	@for nb in $(NOTEBOOKS); do \
		name=$$(basename $$nb .py); \
		echo "exporting $$name"; \
		marimo export html-wasm $$nb -o $(OUT_DIR)/$$name --mode run || exit 1; \
	done
	@echo "done -> $(OUT_DIR)"

serve:
	cd $(SITE_DIR) && hugo server -D

test:
	pytest -q tests/

clean:
	rm -rf $(OUT_DIR)/*
