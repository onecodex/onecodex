test:
	py.test --ignore tests/test_speed.py --ignore tests/test_reports.py tests/
	@echo "Successfully passed all tests (one environment only)."

lint:
	pre-commit run --all-files
	@echo "Successfully linted all files."

coverage:
	py.test --cov-report=term-missing --cov=onecodex --ignore tests/test_speed.py --ignore tests/test_reports.py tests/

coveragehtml:
	py.test --cov-report=html --cov=onecodex --ignore tests/test_speed.py --ignore tests/test_reports.py tests/

install:
	uv sync --all-extras --locked
