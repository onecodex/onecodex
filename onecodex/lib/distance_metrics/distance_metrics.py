import skbio.diversity

from onecodex.lib.distance.helpers import alpha_counts, beta_counts, ACCEPTABLE_FIELDS


def alpha_diversity(classification, distance_metric, ids=None,
                    field='readcount_w_children', rank='species', **kwargs):
    """Caculate the diversity within a community"""
    assert field in ACCEPTABLE_FIELDS

    counts = alpha_counts(classification, field=field, rank=rank)
    return skbio.diversity.alpha_diversity(distance_metric, counts, ids, **kwargs)


def beta_diversity(classification1, classification2, distance_metric,
                   field='readcount_w_children', rank='species', **kwargs):
    """Calculate the diversity between 2 communities"""
    assert field in ACCEPTABLE_FIELDS

    tax_ids, uv_counts, _, _ = beta_counts(classification1, classification2,
                                           field=field, rank=rank)
    # bc_dm = beta_diversity("braycurtis", data, ids)
    return skbio.diversity.beta_diversity(distance_metric, uv_counts, tax_ids)


def simpson(classification, field='readcount_w_children', rank='species'):
    """
    An alpha diversity metric that takes into account the
    number of species present and their abundances.
    """
    assert field in ACCEPTABLE_FIELDS

    counts = alpha_counts(classification, field=field, rank=rank)
    return skbio.diversity.alpha.simpson(counts)


def chao1(classification, bias_corrected=True,
          field='readcount_w_children', rank='species'):
    assert field in ACCEPTABLE_FIELDS

    counts = alpha_counts(classification, field=field, rank=rank)
    return skbio.diversity.alpha.chao1(counts, bias_corrected=bias_corrected)


def unifrac(classification1, classification2, weighted=True,
            field='readcount_w_children', rank='species'):
    """
    A beta diversity metric that takes into account the relative relatedness of community members.
    Weighted UniFrac looks at abundances, unweighted UniFrac looks at presence
    """
    assert field in ACCEPTABLE_FIELDS

    counts = {}
    for row in classification1.results()['table']:
        counts[row['tax_id']] = [row[field], 0]

    for row in classification2.results()['table']:
        if counts[row['tax_id']][1] is 0:
            counts[row['tax_id']][1] = row[field]
        else:
            counts[row['tax_id']] = [0, row[field]]

    tax_ids, _, u_counts, v_counts = beta_counts(classification1, classification2,
                                                 field=field, rank=rank)

    # FIXME: get tree
    # if weighted:
    #     return skbio.diversity.beta.weighted_unifrac(u_counts, v_counts, tax_ids, tree)
    # else:
    #     return skbio.diversity.beta.unweighted_unifrac(u_counts, v_counts, tax_ids, tree)


def jaccard_dissimilarity(classification1, classification2,
                          field='readcount_w_children', rank='species'):
    """Compute the Jaccard dissimilarity between two classifications."""
    assert field in ACCEPTABLE_FIELDS

    _, _, u_counts, v_counts = beta_counts(classification1, classification2,
                                           field=field, rank=rank)
    n_intersection = float(len(set(u_counts) & set(v_counts)))
    n_union = float(len(set(u_counts) | set(v_counts)))
    return 1 - (n_intersection / n_union)


def bray_curtis(classification1, classification2, field='readcount_w_children', rank='species'):
    """Compute the Bray-Curtis dissimilarity between two classifications."""
    assert field in ACCEPTABLE_FIELDS

    tax_ids, uv_counts, _, _ = beta_counts(classification1, classification2, field=field, rank=rank)
    return skbio.diversity.beta_diversity('braycurtis', uv_counts, tax_ids)


def cityblock(classification1, classification2, field='readcount_w_children', rank='species'):
    """Compute the Manhattan distance between two classifications."""
    assert field in ACCEPTABLE_FIELDS

    tax_ids, uv_counts, u_counts, v_counts = beta_counts(classification1, classification2,
                                                         field=field, rank=rank)
    return skbio.diversity.beta_diversity('cityblock', uv_counts, tax_ids)
