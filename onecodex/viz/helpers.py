import pandas as pd

from onecodex.exceptions import MissingDataException
from onecodex.models import Analyses, Samples


def normalize_analyses(analyses):
    normed_analyses = []
    metadata = []

    for a in analyses:
        if isinstance(a, Samples):
            normed_analyses.append(a.primary_classification)
            metadata.append(a.metadata)
        elif isinstance(a, Analyses):
            normed_analyses.append(a)
            metadata.append(a.sample.metadata)

    return normed_analyses, metadata


def collate_analysis_results(analyses, metric='abundance', skip_missing=True):
    """For a set of analyses, return the results as a Pandas DataFrame."""
    assert metric in ['abundance', 'readcount', 'readcount_w_children']
    
    # Keep track of all of the microbial abundances
    dat = {}
    # Keep track of information for each tax_id
    tax_id_info = {}
    
    # Get results for each of the Sample objects that are passed in
    for a in analyses:
        if a.success is False:
            if not skip_missing:
                raise MissingDataException('{} was not successful'.format(a.id))
            continue
            
        # Get the results in table format
        result = a.results()['table']
        
        # Record the information (name, parent) for each organism by  its tax ID
        for d in result:
            if d['tax_id'] not in tax_id_info:
                tax_id_info[d['tax_id']] = {k: d[k] for k in ['name', 'rank', 'parent_tax_id']}

        # Reformat detection infromation as dict of {taxid: value}
        result = {d['tax_id']: d[metric] for d in result if metric in d}

        # Remove entries without the specified metric
        result = {taxid: value for taxid, value in result.items() if value is not None}

        # Remove about 0-values
        result = {taxid: value for taxid, value in result.items() if value > 0}
        
        # Catch any samples that don't have the specified metric
        if len(result) == 0:
            if not skip_missing:
                raise MissingDataException('{} has no entries for {}, skipping.'.format(a.id,
                                                                                        metric))
            continue

        # Save the set of microbial abundances
        dat[str(a.id)] = result

    if len(dat) == 0:
        return None

    # Format as a Pandas DataFrame
    df = pd.DataFrame(dat).fillna(0)

    # add an index with the tax ids name
    df['tax_name'] = df['tax_id'].map(lambda tid: tax_id_info[tid]['name'])
    df.set_index(['tax_name'], append=True)

    # Remove columns (tax_ids) with no values that are > 0
    df = df.T.loc[:, df.sum() > 0]
    
    return df
