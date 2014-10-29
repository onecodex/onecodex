from __future__ import print_function
import sys


def stderr(*objs):
    print(*objs, file=sys.stderr)
