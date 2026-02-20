"""File for amino acid definitions
"""

from collections import Counter

class AminoAcids():
    def __init__(self):
        self.codes321 = {
            'ALA': 'A',
            'ARG': 'R',
            'ASN': 'N',           
            'ASP': 'D',
            'ASH': 'D',
            'CYS': 'C',
            'GLN': 'Q',
            'GLU': 'E',
            'GLH': 'E',
            'GLY': 'G',
            'HIS': 'H',
            'HIE': 'H',
            'HID': 'H',
            'HIP': 'H',
            'ILE': 'I',
            'LEU': 'L',
            'LYS': 'K',
            'MET': 'M',
            'PHE': 'F',
            'PRO': 'P',
            # 'PYL': 'O',
            'SER': 'S',
            # 'SEC': 'U',
            'THR': 'T',
            'TRP': 'W', 
            'TYR': 'Y',
            'VAL': 'V',
            '-': '-',
            # 'NH2': 'Z'
        }
        self.codes123 = {
            'A': 'ALA',
            'R': 'ARG',
            'N': 'ASN',           
            'D': 'ASP',
            'C': 'CYS',
            'Q': 'GLN',
            'E': 'GLU',
            'G': 'GLY',
            'H': 'HIS',
            'I': 'ILE',
            'L': 'LEU',
            'K': 'LYS',
            'M': 'MET',
            'F': 'PHE',
            'P': 'PRO',
            # 'O': 'PYL',
            'S': 'SER',
            # 'U': 'SEC',
            'T': 'THR',
            'W': 'TRP', 
            'Y': 'TYR',
            'V': 'VAL',
            '-': '-',
            # 'Z': 'NH2'
        }
        self.sideChains = {
            'ALA': Counter( 
                N = 1,
                O = 1,
                C = 3
                ),
            'ARG': Counter(
                N = 4,
                O = 1,
                C = 6,
                ),
            'ASN': Counter(
                N = 2,
                O = 2,
                C = 4,
            ),
            'ASP': Counter(
                N = 1, 
                O = 3,
                C = 4
            ),
            'ASH': Counter(
                N = 1,
                O = 3,
                C = 4,
            ),
            'CYS': Counter(
                N = 1,
                O = 1,
                C = 3,
                S = 1
            ),
            'GLN': Counter(
                N = 2,
                O = 2,
                C = 5
            ),
            'GLU': Counter(
                N = 1,
                O = 3,
                C = 5
            ),
            'GLH': Counter(
                N = 1,
                O = 3,
                C = 5
            ),
            'GLY': Counter(
                N = 1,
                O = 1,
                C = 2
            ),
            'HIS': Counter(
                N = 3,
                O = 1,
                C = 6
             ),
            'HIE': Counter(
                N = 3,
                O = 1,
                C = 6
             ),
            'HID': Counter(
                N = 3,
                O = 1,
                C = 6
             ),
            'HIP': Counter(
                N = 3,
                O = 1,
                C = 6
             ),
            'ILE': Counter(
                N = 1,
                O = 1,
                C = 6
            ),
            'LEU': Counter(
                N = 1,
                O = 1,
                C = 6
            ),
            'LYS': Counter(
                N = 2,
                O = 1,
                C = 6
            ),
            'MET': Counter(
                N = 1,
                O = 1,
                C = 5,
                S = 1
            ),
            'PHE': Counter(
                N = 1,
                O = 1, 
                C = 9
            ),
            'PRO': Counter(
                N = 1,
                O = 1,
                C = 5
            ),
            'SER': Counter(
                N = 1,
                O = 2,
                C = 3
            ),
            'THR': Counter(
                N = 1,
                O = 2,
                C = 4
            ),
            'TRP': Counter(
                N = 2,
                O = 1,
                C = 11
            ),
            'TYR': Counter(
                N = 1,
                O = 2,
                C = 9
            ),
            'VAL': Counter(
                N = 1, 
                O = 1,
                C = 5
            ),
            #'NH2': Counter(
            #    N = 1
            #),
        }
