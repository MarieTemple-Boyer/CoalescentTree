""" Show the site frequency spectrum and the branch frequency spectrum of a tree obtained with SLiM """

import numpy as np
import matplotlib.pyplot as plt
import tskit


def plot_sfs_and_bfs(tree_seq, title = None, knowing_ancestral_state=True):
    """ Show a the site frequency spectrum and the branch frequency spectrum of the tskit coalescent tree tree_seq
        - title is the title of the plot (file_name by defautl)
        - knowing_ancestral_state determine if we consider that the ancestral is known
    """

    sfs = tree_seq.allele_frequency_spectrum(
        mode='site', span_normalise=False, polarised=knowing_ancestral_state)

    bfs = tree_seq.allele_frequency_spectrum(
        mode='branch', span_normalise=False, polarised=knowing_ancestral_state)

    plt.subplot(1, 2, 1)

    abs_sfs = np.arange(len(sfs))
    plt.bar(abs_sfs, sfs, width=0.9)
    plt.title('Site frequency spectrum', y=1.03)

    plt.subplot(1, 2, 2)

    abs_bfs = np.arange(len(bfs))
    plt.bar(abs_bfs, bfs, width=0.9)
    plt.title('Branch frequency spectrum', y=1.03)

    if title is not None:
        plt.suptitle(title, y=0.05)

    plt.show()


if __name__ == "__main__":
    import doctest
    doctest.testmod()
