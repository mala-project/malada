"""WILL BE DELETED: Transformation matrices for supercell creation."""

# TODO: I know this can be replaced by a mathematical expression. This is just
# for the first tests.
# Also it is not really correct in its current form.


two_atoms_per_cell = {2: [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
                        4: [[2, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
                      16: [[2, 0, 0],
                            [0, 2, 0],
                            [0, 0, 2]],
                      32: [[4, 0, 0],
                            [0, 2, 0],
                            [0, 0, 2]],
                      64: [[4, 0, 0],
                            [0, 4, 0],
                            [0, 0, 2]],
                       128: [[4, 0, 0],
                             [0, 4, 0],
                             [0, 0, 4]],
                       256: [[8, 0, 0],
                             [0, 4, 0],
                             [0, 0, 4]]}
structure_to_transformation = {"bcc": two_atoms_per_cell,
                               "hcp": two_atoms_per_cell}
