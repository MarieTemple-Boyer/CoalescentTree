"""
Compare a normalised unfolded bfs with normalised theoric bfs (bfs = branch frequency spectrum)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from compare_theory.bfs_kingman import bfs_kingman
from compare_theory.bfs_beta import bfs_beta


def compare_bfs(bfs, kingman=True, beta=False, alpha=0, bfs_label='', title=''):
    """ Compare a normalised unfolded bfs with normalised theoric bfs by ploting them and returninh the L2 error.
        - if kingman then the kingman coalescent theoric unfolded bfs is ploted
        - if beta then the beta coalescent of parameter alpha theoric unfolded bfs is ploted
        - bfs_label correspond the the label of bfs on the plot
        - title correspond to the title of the plot
        - return a tuple of floats of the L2 error when comparing bfs normalised to a theoric unfolded normalised bfs
        (the first one is compared to a kingman coalesent and the second is compard to a beta coalescent of parameter alpha)
        
        Reminder of the definition of the branch frequency spectrum of a coalescent tree:
        bfs is a table of length sample_size-1 where bfs[i] correspond to the total that have exactly (i+1) descendant in the coalescent tree
        
        /!\ When beta is true, not all the values of alpha and of sample_size=len(bfs)+1 are accepted.
            - alpha has to be among these: [1.1, 1.3, 1.5, 1.9, 2]
            - sample_size has to be among these: KNOWN_SAMPLE_SIZE_VALUES = [10]
    
        Don't forget the the beta coalescent of parameter alpha=2 is a Kingman coalescent.

    >>> sample_size=10
    >>> bfs = 1 / np.arange(1, sample_size)
    >>> compare_bfs(bfs)
    (0.0, None)
    >>> compare_bfs(bfs, beta=True, alpha=1.5)
    (0.0, 0.2091883481380384)
    >>> compare_bfs(bfs, beta=True, alpha=1.7)
    Traceback (most recent call last):
        ...
    Exception: The theoric values of the branch frequency spectrum of a beta coalescent for alpha=1.7 are not stored.
    Please choose a value of alpha among these : [1.1, 1.3, 1.5, 1.9, 2].
    >>> compare_bfs(bfs, beta=True, alpha=2)
    (0.0, 0.0)
    """
    
    sample_size = len(bfs)+1
    abscissa = np.arange(1, sample_size)

    bfs_norm = bfs / np.linalg.norm(bfs)
    kingman_bfs_norm = None
    beta_bfs_norm = None

    err_kingman = None
    err_beta = None

    if beta:
        if alpha==2:
            kingman=True
        else:
            beta_bfs_norm = bfs_beta(sample_size, alpha, normalized=True)
            err_beta = np.linalg.norm(beta_bfs_norm-bfs_norm)

    if kingman:
        kingman_bfs_norm = bfs_kingman(sample_size, normalized=True)
        err_kingman = np.linalg.norm(kingman_bfs_norm-bfs_norm)
    
    if alpha==2:
        err_beta = err_kingman

    # plots
    plotdata = pd.DataFrame({
        'Kingman coalescent': kingman_bfs_norm,
        bfs_label: bfs_norm,
        f'Beta coalescent of parameter alpha={alpha}': beta_bfs_norm},
        index=abscissa)

    plotdata.plot(kind='bar', figsize=(15, 8))

    plt.title(title)

    return err_kingman, err_beta


if __name__ == "__main__":
    import doctest
    doctest.testmod()
