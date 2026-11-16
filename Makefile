PYTHON ?= python3

.PHONY: help install test test-fast lint check clean

help:
	@echo "install  install the package and development dependencies"
	@echo "test     run the test suite"
	@echo "test-fast run the test suite without the end-to-end integration tests"
	@echo "lint     run static checks"
	@echo "check    run lint and tests"
	@echo "synthetic generate the artificial demonstration dataset"
	@echo "demo     generate synthetic data and run the pipeline on it"
	@echo "demo-check run the demo twice and check its outputs are byte-identical"
	@echo "clean-clone-check  verify the committed branch from a fresh clone and venv"
	@echo "clean    remove caches and build artefacts (never touches data/)"

install:
	pip install -e ".[dev]"

test:
	pytest -q

test-fast:
	pytest -q -m "not integration"

lint:
	ruff check src scripts tests

check: lint test

clean:
	rm -rf .pytest_cache .ruff_cache build dist *.egg-info
	find . -name __pycache__ -type d -prune -exec rm -rf {} +

.PHONY: synthetic demo demo-check clean-clone-check

synthetic:
	python scripts/make_synthetic_data.py --out data/synthetic --n 400 --seed 20260101

demo: synthetic
	rm -rf results/demo
	MEFV_RAW_AUG=data/synthetic/raw/synthetic_method.csv \
	MEFV_RAW_VARIANTS=data/synthetic/raw/synthetic_variants.csv \
	MEFV_RAW_PATIENTS=data/synthetic/raw/synthetic_patients.csv \
	MEFV_PROCESSED=results/demo/processed \
	MEFV_SOURCE_ANNOTATIONS=data/synthetic/external/source_annotations.csv \
	MEFV_POPULATION_FREQUENCIES=data/synthetic/external/population_frequencies.csv \
	MEFV_ANNOTATION_OVERLAY=data/synthetic/reviewed/variant_annotations_reviewed.csv \
	MEFV_ANNOTATION_DECISION_LOG=data/synthetic/reviewed/annotation_decision_log.csv \
	MEFV_STUDY_PARAMETERS=data/synthetic/study_parameters.json \
	MEFV_ANNOTATIONS_GENERATED=results/demo/external/generated.csv \
	MEFV_DICT=results/demo/external/variant_dictionary.csv \
	MEFV_TABLES=results/demo/tables \
	MEFV_FIGURES=results/demo/figures \
	MEFV_ACMG=results/demo/classification \
	MEFV_RUN_MANIFEST=results/demo/run_manifest.json \
	python scripts/run_pipeline.py --clean
	@echo
	@echo "Demonstration output is in results/demo/. These are NOT study findings."

demo-check:
	$(MAKE) demo
	cp results/demo/run_manifest.json results/demo-first-manifest.json
	$(MAKE) demo
	python scripts/compare_manifests.py results/demo-first-manifest.json results/demo/run_manifest.json

clean-clone-check:
	scripts/check_clean_clone.sh $(PYTHON)
