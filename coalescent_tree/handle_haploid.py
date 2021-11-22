""" Handle an haploid tree sequence obtained with SLiM """

import random
import msprime
import tskit


def handle_haploid(file_name, mutation_rate=None, sample_size=None, upload=False, check_coalescence=False, random_seed=None):
    """ Handle an haploid tree sequence named file_name + '.trees'.
        - remove the nodes that are NULL (trick for modeling haploid individuals with SLiM)
        - checks that the trees have coalesced if mutation_rate is not None
        - add neutral mutations with a rate mutation_rate according to the random seed random_seed
        - conserve only a sample of size sample_size of the current individuals
            (keep the whole tree if sample_size is None)
        - upload the new tree sequence in a file named file_name + '_mutated.trees' if upload is True
        - return the new tree sequence
    """

    tree_seq = tskit.load(file_name + '.trees')
    tree_seq = tree_seq.simplify()

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
            assert tree.num_roots == 1, (
                f'not coalesced! on segment {tree.interval[0]} to {tree.interval[1]}')

        # Add neutral mutation to the coalescent tree
        tree_seq = msprime.mutate(
            tree_seq, rate=mutation_rate, random_seed=random_seed, keep=True)

    # Keep only a sample of the total tree
    if sample_size is not None:
        ind = list(range(tree_seq.num_samples))
        sample = random.sample(ind, sample_size)
        tree_seq = tree_seq.simplify(sample)

    # return and load the new mutated tree
    tree_seq.dump('./' + file_name + '_mutated.trees')
    return tree_seq
