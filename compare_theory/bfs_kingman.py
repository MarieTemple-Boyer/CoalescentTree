"""
Return the theoric branch frequency spectrum for a haploid Kingman coalescent
"""

import numpy as np


def bfs_kingman(sample_size, normalized=False):
    """ Return the theoric branch frequency (bfs) spectrum for a haploid Kingman coalescent
        
        Reminder of the definition of the branch frequency spectrum of a coalescent tree:
        bfs is a table of length sample_size-1 where bfs[i] correspond to the total that have exactly (i+1) descendant in the coalescent tree

        - sample_size is the number of individuals sampled
        - if normalized then the bfs is normalized (by default the bfs is not normalized)
        - return the bfs as a numpy array of size sample_size-1

    >>> sample_size=6
    >>> bfs_kingman(sample_size)
    array([1.        , 0.5       , 0.33333333, 0.25      , 0.2       ])
    >>> bfs_kingman(sample_size, normalized=True)
    array([0.8265843 , 0.41329215, 0.2755281 , 0.20664607, 0.16531686])
    """
    abscissa = np.arange(1, sample_size)
    kingman_bfs = 1/abscissa

    if normalized:
        kingman_bfs = kingman_bfs / np.linalg.norm(kingman_bfs)

    return kingman_bfs


if __name__ == "__main__":
    import doctest
    doctest.testmod()
