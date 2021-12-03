""" Handle a tskit haploid tree sequence obtained with SLiM """

import msprime
import tskit


def handle_haploid(tree_seq, mutation_rate=None, check_coalescence=False, random_seed=None):
    """ Handle a tskit  haploid tree sequence named tree_seq.
        - remove the nodes that are NULL (trick for modeling haploid individuals with SLiM)
        - checks that the trees have coalesced if mutation_rate is not None
        - add neutral mutations with a rate mutation_rate according to the random seed random_seed
        - return the new tree sequence

    >>> tree_coalesced = tskit.load('examples/tree_coalesced.trees')
    >>> tree_coalesced = handle_haploid(tree_coalesced, check_coalescence=True)
    >>> tree_coalesced.num_samples
    20
    >>> tree_not_coalesced = tskit.load('examples/tree_not_coalesced.trees')
    >>> handle_haploid(tree_not_coalesced, check_coalescence=True)
    Traceback (most recent call last):
        ...
    Exception: Tree not coalesced !
    """

    # Removing chromosomes that are NULL

    # load the tree on a table format in order to modify it
    table = tree_seq.dump_tables()

    # computing the list of the id of the node that are not NULL
    node_non_null = []

    for (i, node) in enumerate(table.nodes):
        if not node.metadata['is_null']:
            node_non_null.append(i)

    # remove all the node of the table that are NULL
    table.subset(node_non_null)

    # load the tree in a tree format
    tree_seq = table.tree_sequence()

    if mutation_rate is not None or check_coalescence:
        # Check that the tree has coalesced
        for tree in tree_seq.trees():
            if tree.has_multiple_roots:
                raise Exception('Tree not coalesced !')

        # Add neutral mutation to the coalescent tree
        tree_seq = msprime.mutate(
            tree_seq, rate=mutation_rate, random_seed=random_seed, keep=True)

    # return the new tree
    return tree_seq


if __name__ == "__main__":
    import doctest
    doctest.testmod()
