""" Sample individuals for a tskit coalescent tree """

import random
import tskit


def sampling(tree_seq, sample_size):
    """ Sample individuals from a tskit coalescent tree tree_seq.
        - conserve only sample of size sample size of the tree tree_seq
        - return the new sampled tree

    >>> tree = tskit.load('examples/tree_coalesced.trees')
    >>> tree_sample = sampling(tree, 5)
    >>> tree_sample.num_samples
    5
    """

    ind = list(range(tree_seq.num_samples))
    sample = random.sample(ind, sample_size)
    tree_seq_sampled = tree_seq.simplify(sample)

    return tree_seq_sampled


if __name__ == "__main__":
    import doctest
    doctest.testmod()
