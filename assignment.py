from itertools import product, combinations

translation_table = {
    "TTT": "F",
    "TCT": "S",
    "TAT": "Y",
    "TGT": "C",
    "TTC": "F",
    "TCC": "S",
    "TAC": "Y",
    "TGC": "C",
    "TTA": "L",
    "TCA": "S",
    "TAA": "*",
    "TGA": "*",
    "TTG": "L",
    "TCG": "S",
    "TAG": "*",
    "TGG": "W",
    "CTT": "L",
    "CCT": "P",
    "CAT": "H",
    "CGT": "R",
    "CTC": "L",
    "CCC": "P",
    "CAC": "H",
    "CGC": "R",
    "CTA": "L",
    "CCA": "P",
    "CAA": "Q",
    "CGA": "R",
    "CTG": "L",
    "CCG": "P",
    "CAG": "Q",
    "CGG": "R",
    "ATT": "I",
    "ACT": "T",
    "AAT": "N",
    "AGT": "S",
    "ATC": "I",
    "ACC": "T",
    "AAC": "N",
    "AGC": "S",
    "ATA": "I",
    "ACA": "T",
    "AAA": "K",
    "AGA": "R",
    "ATG": "M",
    "ACG": "T",
    "AAG": "K",
    "AGG": "R",
    "GTT": "V",
    "GCT": "A",
    "GAT": "D",
    "GGT": "G",
    "GTC": "V",
    "GCC": "A",
    "GAC": "D",
    "GGC": "G",
    "GTA": "V",
    "GCA": "A",
    "GAA": "E",
    "GGA": "G",
    "GTG": "V",
    "GCG": "A",
    "GAG": "E",
    "GGG": "G",
}

# nomenclature for degenerate codons
expanded_code = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "W": ["A", "T"],
    "S": ["C", "G"],
    "M": ["A", "C"],
    "K": ["G", "T"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"],
}

# helpful for validating input
valid_nucleotides = "ACGTWSMKRYBDHVN"
valid_aa = "GAVLIMFWPSTCYNQDEKRH*"

# WE WILL BE PRE-COMPUTING AN EXPANDED TRANSLATION TABLE, SO IT CAN BE REUSED FOR REPEAT CALLS OF THE FUNCTION
# There are only 3375 possible permutations of the degenerate nucleotide, which is not too large
possible_nucleotide_permutations = set(
    [
        "".join(degenerate_codon)
        for degenerate_codon in product(
            valid_nucleotides, valid_nucleotides, valid_nucleotides
        )
    ]
)

# generate expanded nucleotide table e.g. {'RYA':['GCA', 'ACA', 'ATA', 'GTA']}
expanded_nucleotides_table = {
    degenerate_codon: [
        "".join(codon)
        for codon in product(
            expanded_code[degenerate_codon[0]],
            expanded_code[degenerate_codon[1]],
            expanded_code[degenerate_codon[2]],
        )
    ]
    for degenerate_codon in possible_nucleotide_permutations
}

# generate expanded translation table e.g. {'RYA':['A', 'I', 'V', 'T']}
expanded_translation_table = {
    degenerate_codon: [translation_table[codon] for codon in codons]
    for degenerate_codon, codons in expanded_nucleotides_table.items()
}


def validate_nucleotides(nucleotides):
    """
    :param nucleotides: string
        the string of nucleotides we want validate, i.e. 'RYA'
    :rtype: bool
        returns True if input string contains only valid nucleotides. Otherwise, returns False.
    """
    return set(valid_nucleotides).issuperset(nucleotides)


def validate_amino_acids(amino_acids):
    """
    :param amino_acids: set
        the set of amino acids we want validate, i.e. {'A','I','V'}
    :rtype: bool
        returns True if input set list contains only valid amino acids. Otherwise, returns False.
    """
    return set(valid_aa).issuperset(amino_acids)


def find_all_codons_for_amino_acids(amino_acids):
    """
    :param amino_acids: set
            the set of valid amino acids we want to code for, i.e. {'A','I','V'}
    :rtype: set
            returns the set of all degenerate codons that code for all amino acids in the set e.g. {'RYA', 'RYH', 'RYC'}
    """
    assert validate_amino_acids(
        amino_acids
    ), f"{amino_acids} contains invalid amino acid"

    if len(amino_acids) == 0:
        return amino_acids

    return set(
        [
            degenerate_codon
            for degenerate_codon, amino_acid in expanded_translation_table.items()
            if set(amino_acid).issuperset(amino_acids)
        ]
    )


def generate_codon_efficiency_table(amino_acids, codons):
    """
    :param amino_acids: set
            the set of valid amino acids we want to code for, i.e. {'A','I','V'}
    :param codons: set
            the set of all degenerate codons that code for all amino acids in the first set e.g. {'RYA', 'RYH', 'RYC'}
    :rtype: dict
            returns dictionary with key as degenerate codon and value as achieved efficiency e.g. {'RYA':0.75}
    """
    assert validate_amino_acids(
        amino_acids
    ), f"{amino_acids} contains invalid amino acid"

    assert validate_nucleotides(
        "".join(codons)
    ), f"Codons contain invalid nucleotide"

    if len(amino_acids) == 0:
        return {}

    # codon efficiency = length of the list only containing the amino acids of interest
    # divided by the length of the list that contains all amino acids the codon codes for
    return {
        codon: len(
            [
                amino_acid
                for amino_acid in expanded_translation_table[codon]
                if amino_acid in amino_acids
            ]
        )
        / len(expanded_translation_table[codon])
        for codon in codons
    }


def get_codon_for_amino_acids(amino_acids):
    """
    :param amino_acids: set
        the valid amino acids we want to code for, i.e. {'A','I','V'}
    :rtype: set, float
        returns two values the set of most efficient codons for the input set list,
        e.g. {'RYA', 'RYH', 'RYC', 'RYW', 'RYM', 'RYY', 'RYT'} and the achieved efficiency e.g. 0.75
    """
    assert validate_amino_acids(
        amino_acids
    ), f"{amino_acids} contains invalid amino acid"

    if len(amino_acids) == 0:
        return amino_acids, 0

    # generate dict with key as degenerate codon and value as efficiency e.g. {'RYA':0.75}
    codon_efficiency_table = generate_codon_efficiency_table(
        amino_acids, find_all_codons_for_amino_acids(amino_acids)
    )

    # calculate codon efficiency and return results
    return (
        set(
            [
                degenerate_codon
                for degenerate_codon, efficiency in codon_efficiency_table.items()
                if efficiency == max(codon_efficiency_table.values())
            ]
        ),
        max(codon_efficiency_table.values()),
    )


def truncate_list_of_amino_acids(amino_acids):
    """
    :param amino_acids: set
        the valid amino acids we want to code for, i.e. {'A','I','V'}
    :rtype: set
        the set of sets of amino acids that can be coded with 100% efficiency,
        i.e. {frozenset({'V', 'A'}), frozenset({'V', 'I'})}
    """
    assert validate_amino_acids(
        amino_acids
    ), f"{amino_acids} contains invalid amino acid"

    if len(amino_acids) == 0:
        return amino_acids

    amino_acid_list_size = len(amino_acids)
    codon_efficiency_table = {}

    # remove one amino acid from the list at a time until there are none left
    while amino_acid_list_size > 0:

        # find all possible combinations of amino acids after removal of 0, 1, 2, etc. amino acids
        for amino_acids_shorter in set(combinations(amino_acids, amino_acid_list_size)):

            # codon_efficiency_table stores all max efficiencies after removal of amino acids
            codon_for_amino_acids = get_codon_for_amino_acids(set(amino_acids_shorter))
            codon_efficiency_table[frozenset(amino_acids_shorter)] = codon_for_amino_acids[1]

        # if amino acids can be encoded with 100% efficiency, store set of amino acids in a list
        if max(codon_efficiency_table.values()) == 1:
            truncate_list = [
                amino_list
                for amino_list, efficiency in codon_efficiency_table.items()
                if efficiency == 1
            ]

            # if the list contains <=1 set of amino acids, return a set; else, return a set of sets
            return set(truncate_list[0]) if len(truncate_list) <= 1 else set(truncate_list)

        else:
            amino_acid_list_size -= 1
            codon_efficiency_table = {}

    # return empty set if task not possible
    return set()


if __name__ == "__main__":
    # using sets instead of lists throughout the code since the order doesn't matter and all items should be unique
    assert get_codon_for_amino_acids({"A", "I", "V"}) == (
        {"RYA", "RYH", "RYC", "RYW", "RYM", "RYY", "RYT"},
        0.75,
    )
    assert get_codon_for_amino_acids({"M", "F"}) == ({"WTS", "WTK", "WTB"}, 0.5)

    # "frozenset" here since this seems to be the only way to get a set of sets
    # - see https://stackoverflow.com/questions/5931291/how-can-i-create-a-set-of-sets-in-python
    assert truncate_list_of_amino_acids({"A", "V", "I"}) == {
        frozenset({"V", "A"}),
        frozenset({"V", "I"}),
    }
