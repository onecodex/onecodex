set -e

uv run make lint
uv run make test
echo "Tests successful. Pushing to PyPI..."
rm -rf dist/
uv build
uv publish
