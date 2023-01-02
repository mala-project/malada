"""Initial guesses for convergence calculations."""
# TODO: Find a smarter way to do this.

cutoff_guesses_qe = {"Fe": [40, 50, 60, 70, 80, 90, 100],
                     "Be": [40, 50, 60, 70],
                     "Al": [20, 30, 40, 50, 60, 70],
                     "H":  [26, 36, 46, 56, 66, 76, 86, 96, 106]}

cutoff_guesses_vasp = {"Fe": [268, 368, 468, 568, 668, 768],
                       "Al": [240, 340, 440, 540],
                       "Be": [248, 258, 268, 278]}

kpoints_guesses = {"Fe": [2, 3, 4],
                   "Be": [2, 3, 4, 5, 6],
                   "Al": [2, 3, 4]}
