""" Give a visualisation of a coalescent tree obtained with SLiM """

import tskit
from cairosvg import svg2png


def visualisation(file_name, print_tree=True, load_tree=False, output_name=None):
    """ Give a visualisation of a tree_sequence named file_name + '.trees'.
        - print the tree as a string if print_tree
        - load an image of the tree sequnce named output_name + '.png'
            (or file_name + '.png' if output_name is None)
    """

    tree_seq = tskit.load(file_name + '.trees')
    tree_seq = tree_seq.simplify()

    if print_tree:
        tree_string = tree_seq.draw_text()
        print(tree_string)

    if load_tree:
        if output_name is None:
            output_name = file_name
        svg_size = (800, 250)
        svg_string = tree_seq.draw_svg(
            size=svg_size,
            y_axis=True, y_label=" ",  # optional: show a time scale on the left
            # Match the axis coordinate systems to the text view
            time_scale="time", x_scale="treewise",
        )
        svg2png(bytestring=svg_string, write_to=output_name+'.png')