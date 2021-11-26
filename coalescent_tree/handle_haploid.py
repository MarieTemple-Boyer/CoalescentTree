""" Handle a tskit haploid tree sequence obtained with SLiM """

import random
import msprime
import tskit


def handle_haploid(tree_seq, mutation_rate=None, sample_size=None, check_coalescence=False, upload=False, new_name=None, random_seed=None):
    """ Handle a tskit  haploid tree sequence named tree_seq.
        - remove the nodes that are NULL (trick for modeling haploid individuals with SLiM)
        - checks that the trees have coalesced if mutation_rate is not None
        - add neutral mutations with a rate mutation_rate according to the random seed random_seed
        - conserve only a sample of size sample_size of the current individuals
            (keep the whole tree if sample_size is None)
        - upload the new tree sequence in a file named new_name (named file_name + '.trees' if new_name is None) if upload is True
        - return the new tree sequence
    
    >>> tree_coalesced = tskit.load('examples/tree_coalesced.trees')
    >>> tree_coalesced_handled = handle_haploid(tree_coalesced, sample_size=5, check_coalescence=True)
    >>> len(tree_coalesced_handled.individuals())
    5
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
            if tree.num_roots != 1:
                raise Exception('Tree not coalesced !')

        # Add neutral mutation to the coalescent tree
        tree_seq = msprime.mutate(
            tree_seq, rate=mutation_rate, random_seed=random_seed, keep=True)

    # Keep only a sample of the total tree
    if sample_size is not None:
        ind = list(range(tree_seq.num_samples))
        sample = random.sample(ind, sample_size)
        tree_seq = tree_seq.simplify(sample)

    # load the new tree
    if upload:
        if new_name is None:
            new_name = file_name + '_handled' 
        tree_seq.dump('./' + new_name + '.trees')

    # return the new tree
    return tree_seq


if __name__ == "__main__":
    import doctest
    doctest.testmod()
