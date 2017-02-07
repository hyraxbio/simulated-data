import numpy
import re


def get_motif_indices(sequence, motif):
    """
    Find indices of substring motif in sequence.

    Args:
        sequence: a string
        motif: a substring

    Returns: a list of left indices
    """
    indices = [m.start() for m in re.finditer(re.escape(motif), sequence)]
    return indices

if __name__ == '__main__':
    pass
