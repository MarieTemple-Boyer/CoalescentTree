"""
Return the theoric branch frequency spectrum for a haploid Kingman coalescent
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from compare_theory.bfs_kingman import bfs_kingman

KNOWN_ALPHA_VALUES = [1.1, 1.3, 1.5, 1.9, 2]
KNOWN_SAMPLE_SIZE_VALUES = [10]

BETA_BFS_VALUES = {
    '1.1': [0.580175, 0.119103, 0.066440, 0.047197, 0.038166, 0.033879, 0.032796, 0.035382, 0.046863],
    '1.3': [0.521296, 0.137166, 0.078487, 0.056070, 0.045115, 0.039481, 0.037258, 0.038479, 0.046649],
    '1.5': [0.467491, 0.152216, 0.090245, 0.065103, 0.052216, 0.045067, 0.041436, 0.040898, 0.045330],
    '1.9': [0.374086, 0.173264, 0.112565, 0.083644, 0.066914, 0.056165, 0.048856, 0.043826, 0.040681]}


def bfs_beta(sample_size, alpha, normalized=False):
    """ Return the theoric branch frequency (bfs) spectrum for a haploid beta coalescent of parameter alpha
        
        Reminder of the definition of the branch frequency spectrum of a coalescent tree:
        bfs is a table of length sample_size-1 where bfs[i] correspond to the total that have exactly (i+1) descendant in the coalescent tree

        - sample_size is the number of individuals sampled
        - if normalized then the bfs is normalized (by default the bfs is not normalized)
        - return the bfs as a numpy array of size sample_size-1

        /!\ When beta is true, not all the values of alpha and of sample_size=len(bfs)+1 are accepted.
            - alpha=2 is always accepted
            - else alpha has to be among these: [1.1, 1.3, 1.5, 1.9, 2]
               and sample_size has to be among these: [10]
    
        Don't forget the the beta coalescent of parameter alpha=2 is a Kingman coalescent.
    
    >>> bfs_beta(sample_size=10, alpha=1.3)
    array([0.521296, 0.137166, 0.078487, 0.05607 , 0.045115, 0.039481,
           0.037258, 0.038479, 0.046649])
    >>> bfs_beta(sample_size=10, alpha=1.7)
    Traceback (most recent call last):
        ...
    Exception: The theoric values of the branch frequency spectrum of a beta coalescent for alpha=1.7 are not stored.
    Please choose a value of alpha among these : [1.1, 1.3, 1.5, 1.9, 2].
    >>> bfs_beta(sample_size=6, alpha=2, normalized=True)
    array([0.8265843 , 0.41329215, 0.2755281 , 0.20664607, 0.16531686])
    """
    abscissa = np.arange(1, sample_size)

    if alpha == 2:
        # the beta coalescent for alpha=2 is a Kingman coalescent
        beta_bfs = 1 / abscissa
    elif alpha not in KNOWN_ALPHA_VALUES:
        raise Exception(
            f'The theoric values of the branch frequency spectrum of a beta coalescent for alpha={alpha} are not stored.\nPlease choose a value of alpha among these : {KNOWN_ALPHA_VALUES}.')

    elif sample_size not in KNOWN_SAMPLE_SIZE_VALUES:
        raise Exception(
            f'The theoric values of the branch frequency spectrum of a beta coalescent for a sample size {sample_size} are not stored.\nPlease choose a value of the sample size among these : {KNOWN_SAMPLE_SIZE_VALUES}.')

    else:
        beta_bfs = np.array(BETA_BFS_VALUES[f'{alpha}'])

    if normalized:
        beta_bfs = beta_bfs / np.linalg.norm(beta_bfs)

    return beta_bfs


if __name__ == "__main__":
    import doctest
    doctest.testmod()
