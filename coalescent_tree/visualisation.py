""" Give a visualisation of a coalescent tree obtained with SLiM """

import tskit
from cairosvg import svg2png


def visualisation(tree_seq, print_tree=True, load_tree=False, output_name=None, svg_size=None, time_scale='rank'):
    """ Give a visualisation of a tskit tree_sequence named tree_seq.
        - print the tree as a string if print_tree
        - load an image of of the tree sequence named output_name + '.png'
            (or 'tree.png' if output_name is None) if load_tree
        - svg_size is a couple of intergers that corresponds to the size of the image 
        - time_scale define the time scale considered
    
    >>> small_tree = tskit.load('examples/small_tree.trees')
    >>> visualisation(small_tree)
    4409.00┊   6     ┊
           ┊ ┏━┻━┓   ┊
    206.00 ┊ ┃   5   ┊
           ┊ ┃  ┏┻━┓ ┊
    90.00  ┊ ┃  4  ┃ ┊
           ┊ ┃ ┏┻┓ ┃ ┊
    0.00   ┊ 0 1 2 3 ┊
           0         1
    <BLANKLINE>
    """
        
    if print_tree:
        tree_string = tree_seq.draw_text()
        print(tree_string)

    if load_tree:
        if output_name is None:
            output_name = 'tree'
        if svg_size is None:
            svg_size = (800, 250)
        svg_string = tree_seq.draw_svg(
            size=svg_size,
            y_axis=True, y_label=" ",  # optional: show a time scale on the left
            # Match the axis coordinate systems to the text view
            time_scale=time_scale, x_scale="treewise",
        )
        svg2png(bytestring=svg_string, write_to=output_name+'.png')



if __name__ == "__main__":
    import doctest
    doctest.testmod()
