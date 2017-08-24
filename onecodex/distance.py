import skbio.diversity

from skbio.tree import TreeNode

from onecodex.exceptions import OneCodexException
from onecodex.helpers import normalize_classifications, collate_classification_results
from onecodex.taxonomy import generate_skbio_tree, prune_to_rank


__all__ = ['alpha_diversity', 'simpson', 'chao1',
           'beta_diversity', 'unifrac', 'jaccard', 'braycurtis', 'cityblock']


ACCEPTABLE_FIELDS = ['abundance', 'readcount_w_children', 'readcount']


def alpha_counts(classification, field='readcount_w_children', rank='species'):
    counts = [t[field] for t in classification.results()['table']
              if t['rank'] == rank]
    name = classification.sample.metadata.name \
        if classification.sample.metadata.name else classification.sample.filename
    ids = [name]
    return (counts, ids)


def beta_counts(classifications, field='readcount_w_children', rank='species'):
    normed_classifications, _ = normalize_classifications(classifications)
    df, tax_info = collate_classification_results(normed_classifications, field=field, rank=rank)

    tax_ids = df.columns.values.tolist()
    vectors = df.values.tolist()
    vectors = [[int(i) for i in row] for row in vectors]
    ids = df.index.values
    return (vectors, tax_ids, ids)


def alpha_diversity(classification, metric, ids=None,
                    field='readcount_w_children', rank='species', **kwargs):
    """Caculate the diversity within a community"""
    assert field in ACCEPTABLE_FIELDS
    counts, ids = alpha_counts(classification, field=field, rank=rank)
    return skbio.diversity.alpha_diversity(metric, counts, ids, **kwargs).values[0]


def beta_diversity(classifications, metric,
                   field='readcount_w_children', rank='species', **kwargs):
    """Calculate the diversity between 2 communities"""
    assert field in ACCEPTABLE_FIELDS
    counts, tax_ids, ids = beta_counts(classifications, field=field, rank=rank)
    return skbio.diversity.beta_diversity(metric, counts, ids, **kwargs)


def simpson(classification, field='readcount_w_children', rank='species'):
    """
    An alpha diversity metric that takes into account the
    number of species present and their abundances.
    """
    assert field in ACCEPTABLE_FIELDS
    counts, ids = alpha_counts(classification, field=field, rank=rank)
    return skbio.diversity.alpha_diversity('simpson', counts, ids).values[0]


def chao1(classification, bias_corrected=True,
          field='readcount_w_children', rank='species'):
    assert field in ACCEPTABLE_FIELDS
    counts, ids = alpha_counts(classification, field=field, rank=rank)
    return skbio.diversity \
                .alpha_diversity('chao1', counts, ids, bias_corrected=bias_corrected) \
                .values[0]


def unifrac(classifications, weighted=True,
            field='readcount_w_children', rank='species', strict=False):
    """
    A beta diversity metric that takes into account the relative relatedness of community members.
    Weighted UniFrac looks at abundances, unweighted UniFrac looks at presence
    """
    assert field in ACCEPTABLE_FIELDS
    counts, tax_ids, ids = beta_counts(classifications, field=field, rank=rank)

    tree = None
    for c in classifications:
        if strict and c.job.id != classifications[0].job.id:
            raise OneCodexException('All Classifications must have the same Job for Unifrac')
        tree = generate_skbio_tree(c, existing_tree=tree)

    # there's a bug (?) in skbio where it expects the root to only have
    # one child, so we do a little faking here
    new_tree = TreeNode(name='fake root')
    new_tree.rank = 'no rank'
    new_tree.append(tree)

    # prune low-level nodes off the tree so the tips are what we're comparing
    prune_to_rank(new_tree, rank=rank)

    if weighted:
        return skbio.diversity.beta_diversity('weighted_unifrac', counts, ids,
                                              tree=new_tree, otu_ids=tax_ids)
    else:
        return skbio.diversity.beta_diversity('unweighted_unifrac', counts, ids,
                                              tree=new_tree, otu_ids=tax_ids)


def jaccard(classifications, field='readcount_w_children', rank='species'):
    """Compute the Jaccard dissimilarity between two classifications."""
    assert field in ACCEPTABLE_FIELDS
    counts, tax_ids, ids = beta_counts(classifications, field=field, rank=rank)
    return skbio.diversity.beta_diversity('jaccard', counts, ids)


def braycurtis(classifications, field='readcount_w_children', rank='species'):
    """Compute the Bray-Curtis dissimilarity between two classifications."""
    assert field in ACCEPTABLE_FIELDS
    counts, _, ids = beta_counts(classifications, field=field, rank=rank)
    return skbio.diversity.beta_diversity('braycurtis', counts, ids)


def cityblock(classifications, field='readcount_w_children', rank='species'):
    """Compute the Manhattan distance between two classifications."""
    assert field in ACCEPTABLE_FIELDS
    counts, tax_ids, ids = beta_counts(classifications, field=field, rank=rank)
    return skbio.diversity.beta_diversity('cityblock', counts, ids)
