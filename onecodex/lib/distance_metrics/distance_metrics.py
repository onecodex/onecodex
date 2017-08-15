import skbio.diversity

from onecodex.lib.distance.helpers import alpha_counts, beta_counts, ACCEPTABLE_FIELDS
from onecodex.lib.taxonomy import generate_skbio_tree


def alpha_diversity(classification, distance_metric, ids=None,
                    field='readcount_w_children', rank='species', **kwargs):
    """Caculate the diversity within a community"""
    assert field in ACCEPTABLE_FIELDS

    counts, ids = alpha_counts(classification, field=field, rank=rank)
    return skbio.diversity.alpha_diversity(distance_metric, counts, ids, **kwargs)


def beta_diversity(classifications, distance_metric,
                   field='readcount_w_children', rank='species', **kwargs):
    """Calculate the diversity between 2 communities"""
    assert field in ACCEPTABLE_FIELDS

    counts, tax_ids, ids = beta_counts(classifications, field=field, rank=rank)
    return skbio.diversity.beta_diversity(distance_metric, counts, ids, otu_ids=tax_ids)


def simpson(classification, field='readcount_w_children', rank='species'):
    """
    An alpha diversity metric that takes into account the
    number of species present and their abundances.
    """
    assert field in ACCEPTABLE_FIELDS

    counts, ids = alpha_counts(classification, field=field, rank=rank)
    return skbio.diversity.alpha.simpson(counts, ids)


def chao1(classification, bias_corrected=True,
          field='readcount_w_children', rank='species'):
    assert field in ACCEPTABLE_FIELDS

    counts, ids = alpha_counts(classification, field=field, rank=rank)
    return skbio.diversity.alpha.chao1(counts, ids=ids, bias_corrected=bias_corrected)


def unifrac(classifications, weighted=True,
            field='readcount_w_children', rank='species'):
    """
    A beta diversity metric that takes into account the relative relatedness of community members.
    Weighted UniFrac looks at abundances, unweighted UniFrac looks at presence
    """
    assert field in ACCEPTABLE_FIELDS

    counts, tax_ids, ids = beta_counts(classifications, field=field, rank=rank)

    tree = None
    for c in classifications:
        tree = generate_skbio_tree(c, existing_tree=tree)

    if weighted:
        return skbio.diversity.beta_diversity("weighted_unifrac", counts, ids,
                                              tree=tree, otu_ids=tax_ids)
    else:
        return skbio.diversity.beta_diversity("unweighted_unifrac", counts, ids,
                                              tree=tree, otu_ids=tax_ids)


def jaccard_dissimilarity(classifications, field='readcount_w_children', rank='species'):
    """Compute the Jaccard dissimilarity between two classifications."""
    assert field in ACCEPTABLE_FIELDS
    counts, tax_ids, ids = beta_counts(classifications, field=field, rank=rank)
    return skbio.diversity.beta_diversity('jaccard', counts, ids, otu_ids=tax_ids)


def bray_curtis(classifications, field='readcount_w_children', rank='species'):
    """Compute the Bray-Curtis dissimilarity between two classifications."""
    assert field in ACCEPTABLE_FIELDS

    counts, tax_ids, ids = beta_counts(classifications, field=field, rank=rank)
    return skbio.diversity.beta_diversity('braycurtis', counts, ids, otu_ids=tax_ids)


def cityblock(classifications, field='readcount_w_children', rank='species'):
    """Compute the Manhattan distance between two classifications."""
    assert field in ACCEPTABLE_FIELDS

    counts, tax_ids, ids = beta_counts(classifications, field=field, rank=rank)
    return skbio.diversity.beta_diversity('cityblock', counts, ids, otu_ids=tax_ids)
