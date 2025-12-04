from collections import deque
from typing import Union
import time
from hashlib import sha256


class RateLimiter:
    """Simple rate limiter: at most `max_actions` per `period` seconds."""

    def __init__(self, max_actions: int, period: float):
        self.max_actions = max_actions
        self.period = period
        self.actions = deque()  # timestamps of recent actions

    def acquire(self):
        now = time.time()
        # Remove timestamps older than `period`
        while self.actions and self.actions[0] <= now - self.period:
            self.actions.popleft()

        if len(self.actions) >= self.max_actions:
            # Sleep until the oldest action falls out of the window
            sleep_for = self.period - (now - self.actions[0])
            if sleep_for > 0:
                time.sleep(sleep_for)

            # Clean up
            now = time.time()
            while self.actions and self.actions[0] <= now - self.period:
                self.actions.popleft()

        # Record this action
        self.actions.append(time.time())


def hash_to_hex(value: Union[str, tuple]) -> str:
    """Hash a string or a tuple of strings into a 128-bit hex string."""

    h = sha256()

    if isinstance(value, str):
        h.update(b"S")
        h.update(value.encode("utf-8"))
    else:
        h.update(b"T")
        for part in value:
            part_bytes = part.encode("utf-8")
            # Length-prefix each part to avoid ambiguity
            h.update(len(part_bytes).to_bytes(4, "big"))
            h.update(part_bytes)

    full_digest = h.digest()  # 32 bytes (256 bits)
    digest_128 = full_digest[:16]
    return digest_128.hex()  # 16 bytes = 32 hex chars
