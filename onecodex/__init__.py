import logging
import sys

__all__ = ["Api"]

formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
stream_handler = logging.StreamHandler(sys.stderr)
stream_handler.setFormatter(formatter)

log = logging.getLogger("onecodex")
log.setLevel(logging.INFO)
log.addHandler(stream_handler)

from onecodex.api import Api  # noqa
from onecodex.cli import onecodex as Cli  # noqa
