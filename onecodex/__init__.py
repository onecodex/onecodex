"""
__init__.py
author: @mbiokyle29
"""
import logging

__all__ = ["Api"]

log = logging.getLogger(__name__)
log.setLevel(logging.WARN)
log_formatter = logging.Formatter('%(asctime)s {%(levelname)s}: %(message)s')

from onecodex.api import Api  # noqa
from onecodex.cli import onecodex as Cli  # noqa
