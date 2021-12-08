"""
Compare a normalised unfolded bfs with normalised theoric bfs (bfs = branch frequency spectrum)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

KNOWN_ALPHA_VALUES = [1.1, 1.3, 1.5, 1.9, 2]
KNOWN_SAMPLE_SIZE_VALUES = [10]

BETA_BFS_VALUES = {
    '1.1': [0.580175, 0.119103, 0.066440, 0.047197, 0.038166, 0.033879, 0.032796, 0.035382, 0.046863],
    '1.3': [0.521296, 0.137166, 0.078487, 0.056070, 0.045115, 0.039481, 0.037258, 0.038479, 0.046649],
    '1.5': [0.467491, 0.152216, 0.090245, 0.065103, 0.052216, 0.045067, 0.041436, 0.040898, 0.045330],
    '1.9': [0.374086, 0.173264, 0.112565, 0.083644, 0.066914, 0.056165, 0.048856, 0.043826, 0.040681]}


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
    >>> compare_bfs(bfs, beta=True, alpha=1.55)
    Traceback (most recent call last):
        ...
    Exception: The theoric values of the branch frequency spectrum of a beta coalescent for alpha=1.55 are not stored.
    Please choose a value of alpha among these : [1.1, 1.3, 1.5, 1.9, 2].
    """
    sample_size = len(bfs)+1
    abscissa = np.arange(1, sample_size)

    bfs_norm = bfs / np.linalg.norm(bfs)
    kingman_bfs_norm = None
    beta_bfs_norm = None

    err_kingman = None
    err_beta = None

    if beta:
        if alpha == 2:
            # the beta coalescent for alpha=2 is a Kingman coalescent
            kingman = True
        elif alpha not in KNOWN_ALPHA_VALUES:
            raise Exception(
                f'The theoric values of the branch frequency spectrum of a beta coalescent for alpha={alpha} are not stored.\nPlease choose a value of alpha among these : {KNOWN_ALPHA_VALUES}.')

        elif sample_size not in KNOWN_SAMPLE_SIZE_VALUES:
            raise Exception(
                f'The theoric values of the branch frequency spectrum of a beta coalescent for a sample size {sample_size} are not stored.\nPlease choose a value of the sample size among these : {KNOWN_SAMPLE_SIZE_VALUES}.')

        else:
            beta_bfs = BETA_BFS_VALUES[f'{alpha}']
            beta_bfs_norm = beta_bfs / np.linalg.norm(beta_bfs)
            err_beta = np.linalg.norm(beta_bfs_norm-bfs_norm)

    if kingman:
        kingman_bfs = 1/abscissa
        kingman_bfs_norm = kingman_bfs / np.linalg.norm(kingman_bfs)
        err_kingman = np.linalg.norm(kingman_bfs_norm-bfs_norm)

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
