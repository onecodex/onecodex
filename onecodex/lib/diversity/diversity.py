import skbio.diversity

from onecodex.distance.helpers import alpha_counts, beta_counts


def alpha_diversity(classification, metric, ids=None,
                    calculate_by='readcount_w_children', rank='species', **kwargs):
    # calculate_by: either 'abundance' or 'readcount_w_children'
    """
    Caculate the diversity within a community
    """
    assert calculate_by == 'abundance' or calculate_by == 'readcount_w_children'

    counts = alpha_counts(classification, calculate_by=calculate_by, rank=rank)
    return skbio.diversity.alpha_diversity(metric, counts, ids, **kwargs)


def beta_diversity(classification1, classification2, metric,
                   calculate_by='readcount_w_children', rank='species', **kwargs):
    # calculate_by: either 'abundance' or 'readcount_w_children'
    """
    Calculate the diversity between 2 communities
    """
    assert calculate_by == 'abundance' or calculate_by == 'readcount_w_children'

    tax_ids, uv_counts, _, _ = beta_counts(classification1, classification2,
                                           calculate_by=calculate_by, rank=rank)
    return skbio.diversity.beta.pw_distances(uv_counts, metric, tax_ids, **kwargs)


def simpson(classification, calculate_by='readcount_w_children', rank='species'):
    """
    An alpha diversity metric that takes into account the
    number of species present and their abundances.
    """
    assert calculate_by == 'abundance' or calculate_by == 'readcount_w_children'

    counts = alpha_counts(classification, calculate_by=calculate_by, rank=rank)
    return skbio.diversity.alpha.simpson(counts)


def unifrac(classification1, classification2, weighted=True,
            calculate_by='readcount_w_children', rank='species'):
    """
    A beta diversity metric that takes into account the relative relatedness of community members.
    Weighted UniFrac looks at abundances, unweighted UniFrac looks at presence
    """
    assert calculate_by == 'abundance' or calculate_by == 'readcount_w_children'

    counts = {}
    for row in classification1.results()['table']:
        counts[row['tax_id']] = [row[calculate_by], 0]

    for row in classification2.results()['table']:
        if counts[row['tax_id']][1] is 0:
            counts[row['tax_id']][1] = row[calculate_by]
        else:
            counts[row['tax_id']] = [0, row[calculate_by]]

    tax_ids, _, u_counts, v_counts = beta_counts(classification1, classification2,
                                                 calculate_by=calculate_by, rank=rank)

    # if weighted:
    #     return skbio.diversity.beta.weighted_unifrac(u_counts, v_counts, tax_ids, tree)
    # else:
    #     return skbio.diversity.beta.unweighted_unifrac(u_counts, v_counts, tax_ids, tree)


def chao1(classification, bias_corrected=True,
          calculate_by='readcount_w_children', rank='species'):
    counts = alpha_counts(classification, calculate_by=calculate_by, rank=rank)
    return skbio.diversity.alpha.chao1(counts, bias_corrected=bias_corrected)
