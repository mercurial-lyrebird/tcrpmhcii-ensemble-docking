
from Bio import Align


BLOSUM62 = Align.substitution_matrices.load("BLOSUM62")
ALPHABET = BLOSUM62.alphabet
