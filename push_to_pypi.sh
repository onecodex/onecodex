set -e

source venv/bin/activate
make lint
make test
echo "Tests successful. Pushing to PyPI..."
python setup.py sdist
twine upload dist/*
