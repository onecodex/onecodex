import pytest

from onecodex.lib.enums import BaseEnum, Rank, Metric


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


def test_metric_is_abundance_metric():
    assert Metric.Abundance.is_abundance_metric
    assert Metric.AbundanceWChildren.is_abundance_metric
    assert not Metric.Readcount.is_abundance_metric
    assert not Metric.NormalizedReadcount.is_abundance_metric
    assert not Metric.PropClassified.is_abundance_metric


def test_metric_includes_children():
    for metric in Metric:
        assert metric.includes_children == ("children" in metric.value)


def test_metric_dtype():
    """this just makes sure that we have defined a dtype mapping for all metrics"""
    for metric in set(Metric) - {Metric.Auto}:
        assert metric.dtype is not None


def test_metric_results_key():
    """this just makes sure that we have defined a results key mapping for all metrics"""
    for metric in set(Metric) - {Metric.Auto}:
        assert metric.results_key


def test_metric_is_normalized():
    """this just makes sure that we have defined an is_normalized mapping for all metrics"""
    for metric in set(Metric) - {Metric.Auto}:
        assert metric.is_normalized is not None
