ACCEPTABLE_FIELDS = ['abundance', 'readcount_w_children', 'readcount']


def alpha_counts(classification, field='readcount_w_children', rank='species'):
    return [t[field] for t in classification.results()['table']
            if t['rank'] == rank]


def beta_counts(classification1, classification2,
                field='readcount_w_children', rank='species'):
        counts = {}
        for row in classification1.results()['table']:
            counts[row['tax_id']] = [row[field], 0]

        for row in classification2.results()['table']:
            if counts[row['tax_id']][1] is 0:
                counts[row['tax_id']][1] = row[field]
            else:
                counts[row['tax_id']] = [0, row[field]]

        tax_ids = counts.keys()
        uv_counts = counts.values()
        u_counts = [count[0] for count in uv_counts]
        v_counts = [count[1] for count in uv_counts]
        return (tax_ids, uv_counts, u_counts, v_counts)
