import pytest

from onecodex.lib.enums import BaseEnum, Rank


def test_base_enum_format_compat():
    class MyEnum(BaseEnum):
        Foo = "foo"

    assert f"{MyEnum.Foo}" == "foo"


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
