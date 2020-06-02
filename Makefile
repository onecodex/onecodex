test:
	py.test -W ignore::DeprecationWarning --ignore tests/test_speed.py --ignore tests/test_reports.py tests/
	@echo "Successfully passed all tests (one environment only, use tox to full suite)."

lint:
	pre-commit run --all-files
	@echo "Successfully linted all files."

coverage:
	py.test -W ignore::DeprecationWarning --cov-report=term-missing --cov=onecodex --ignore tests/test_speed.py --ignore tests/test_reports.py tests/

coveragehtml:
	py.test -W ignore::DeprecationWarning --cov-report=html --cov=onecodex --ignore tests/test_speed.py --ignore tests/test_reports.py tests/

install:
	python setup.py install

format:
	black -l 100 --exclude vendored/* onecodex/ tests/
