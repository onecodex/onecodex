set -e

uv run make lint
uv run make test
echo "Tests successful. Pushing to PyPI..."
uv build
uv publish
