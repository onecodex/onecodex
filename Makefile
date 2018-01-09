test:
	py.test --ignore tests/test_minimal.py --ignore tests/test_simplejson.py tests/
	@echo "Successfully passed all tests (one environment only, use tox to full suite)."

lint:
	flake8 --ignore E501 onecodex/
	flake8 --ignore E501 tests/
	@echo "Successfully linted all files."

coverage:
	py.test --cov-report=term-missing --cov=onecodex --ignore tests/test_minimal.py --ignore tests/test_simplejson.py tests/

coveragehtml:
	py.test --cov-report=html --cov=onecodex --ignore tests/test_minimal.py --ignore tests/test_simplejson.py tests/

install:
	python setup.py install
