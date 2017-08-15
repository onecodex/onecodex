ACCEPTABLE_FIELDS = ['abundance', 'readcount_w_children', 'readcount']


def alpha_counts(classification, field='readcount_w_children', rank='species'):
    counts = [t[field] for t in classification.results()['table']
              if t['rank'] == rank]
    name = classification.results()['table']['name'] \
        if classification.results()['table']['name'] \
        else classification.results()['table']['file_name']
    ids = [name]
    return (counts, ids)


def beta_counts(classifications, field='readcount_w_children', rank='species'):
    ids = []

    counts = {}
    for c, idx in enumerate(classifications):
        name = c.results()['table']['name'] if c.results()['table']['name'] \
                                            else c.results()['table']['file_name']
        ids.append(name)
        for row in c.results()['table']:
            counts[row['tax_id']][idx] = row[field]

    for k in counts:
        counts[k] = [counts[k][i] if counts[k][i] else 0 for i in counts[k]]

    counts = counts.values()
    tax_ids = counts.keys()
    return (counts, tax_ids, ids)
