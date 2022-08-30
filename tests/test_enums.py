import pytest

from onecodex.lib.enums import Rank


def test_rank_level():
    levels = []
    for rank in Rank:
        if rank == Rank.Auto:
            with pytest.raises(ValueError):
                rank.level
        else:
            level = rank.level
            levels.append(level)
            assert isinstance(level, int)

    assert len(levels) == len(set(levels))
