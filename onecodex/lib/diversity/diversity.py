import skbio.diversity

from onecodex.lib.distance.helpers import alpha_counts, beta_counts


ACCEPTABLE_METRICS = ['abundance', 'readcount_w_children']


def alpha_diversity(classification, diversity_index, ids=None,
                    metric='readcount_w_children', rank='species', **kwargs):
    """
    Caculate the diversity within a community
    """
    assert metric in ACCEPTABLE_METRICS

    counts = alpha_counts(classification, calculate_by=metric, rank=rank)
    return skbio.diversity.alpha_diversity(diversity_index, counts, ids, **kwargs)


def beta_diversity(classification1, classification2, diversity_index,
                   metric='readcount_w_children', rank='species', **kwargs):
    """
    Calculate the diversity between 2 communities
    """
    assert metric in ACCEPTABLE_METRICS

    tax_ids, uv_counts, _, _ = beta_counts(classification1, classification2,
                                           calculate_by=metric, rank=rank)
    return skbio.diversity.beta.pw_distances(uv_counts, diversity_index, tax_ids, **kwargs)


def simpson(classification, metric='readcount_w_children', rank='species'):
    """
    An alpha diversity metric that takes into account the
    number of species present and their abundances.
    """
    assert metric in ACCEPTABLE_METRICS

    counts = alpha_counts(classification, calculate_by=metric, rank=rank)
    return skbio.diversity.alpha.simpson(counts)


def unifrac(classification1, classification2, weighted=True,
            metric='readcount_w_children', rank='species'):
    """
    A beta diversity metric that takes into account the relative relatedness of community members.
    Weighted UniFrac looks at abundances, unweighted UniFrac looks at presence
    """
    assert metric in ACCEPTABLE_METRICS

    counts = {}
    for row in classification1.results()['table']:
        counts[row['tax_id']] = [row[metric], 0]

    for row in classification2.results()['table']:
        if counts[row['tax_id']][1] is 0:
            counts[row['tax_id']][1] = row[metric]
        else:
            counts[row['tax_id']] = [0, row[metric]]

    tax_ids, _, u_counts, v_counts = beta_counts(classification1, classification2,
                                                 calculate_by=metric, rank=rank)

    # FIXME: get tree
    # if weighted:
    #     return skbio.diversity.beta.weighted_unifrac(u_counts, v_counts, tax_ids, tree)
    # else:
    #     return skbio.diversity.beta.unweighted_unifrac(u_counts, v_counts, tax_ids, tree)


def chao1(classification, bias_corrected=True,
          metric='readcount_w_children', rank='species'):
    assert metric in ACCEPTABLE_METRICS

    counts = alpha_counts(classification, calculate_by=metric, rank=rank)
    return skbio.diversity.alpha.chao1(counts, bias_corrected=bias_corrected)
