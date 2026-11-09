.PHONY: test quickstart

test:
	python -m unittest discover -s tests -v

quickstart:
	python scripts/run_pipeline.py --config examples/minimal/config.yaml --workdir out/minimal
