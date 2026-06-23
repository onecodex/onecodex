set -e

if [ -e dist/ ]; then
  echo "--- dist/ already exists. Remove or move before running"
  exit 1
fi

uv run make lint
uv run make test
echo "Tests successful. Pushing to PyPI..."
uv build
uv publish
