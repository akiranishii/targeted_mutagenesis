import unittest
from assignment import (
    validate_nucleotides,
    validate_amino_acids,
    find_all_codons_for_amino_acids,
    get_codon_for_amino_acids,
    truncate_list_of_amino_acids,
)


class Assignment(unittest.TestCase):
    """
    Tests for helper functions, get_codon_for_amino_acids, and truncate_list_of_amino_acids
    """

    def test_validate_nucleotides(self):

        # (1) test for empty string
        nucleotides = ""
        self.assertEqual(validate_nucleotides(nucleotides), True)

        # (2) test for valid amino_acid
        nucleotides = "RYA"
        self.assertEqual(validate_nucleotides(nucleotides), True)

        # (3) test for one invalid amino_acid
        nucleotides = "ROA"
        self.assertEqual(validate_nucleotides(nucleotides), False)

    def test_validate_amino_acids(self):

        # (1) test for empty set
        amino_acids = set()
        self.assertEqual(validate_amino_acids(amino_acids), True)

        # (2) test for valid amino_acid
        amino_acids = {"A", "V", "I"}
        self.assertEqual(validate_amino_acids(amino_acids), True)

        # (3) test for one invalid amino_acid
        amino_acids = {"A", "O", "I"}
        self.assertEqual(validate_amino_acids(amino_acids), False)

    def test_find_all_codons_for_amino_acids(self):

        # (1) test for empty set
        amino_acids = set()
        expected_result = set()
        self.assertEqual(find_all_codons_for_amino_acids(amino_acids), expected_result)

        # (2) test for valid amino_acid (expected result extracted using expanded translation table)
        amino_acids = {
            "K",
            "F",
            "I",
            "V",
            "R",
            "E",
            "D",
            "S",
            "Y",
            "C",
            "N",
            "T",
            "A",
            "M",
            "W",
            "G",
            "L",
        }
        expected_result = {
            "DNK",
            "DNV",
            "NNB",
            "NNV",
            "NNN",
            "DNS",
            "DNN",
            "DND",
            "NNK",
            "DNB",
            "NNS",
            "NND",
        }
        self.assertEqual(find_all_codons_for_amino_acids(amino_acids), expected_result)

    def test_get_codon_for_amino_acids(self):

        # (1) test for empty set
        amino_acids = set()
        expected_result = (set(), 0)
        self.assertEqual(get_codon_for_amino_acids(amino_acids), expected_result)

        # (2) test for valid amino_acids
        amino_acids = {"A", "I", "V"}
        expected_result = ({"RYA", "RYH", "RYC", "RYW", "RYM", "RYY", "RYT"}, 0.75)
        self.assertEqual(get_codon_for_amino_acids(amino_acids), expected_result)

        # (3) test for valid amino_acids
        amino_acids = {"M", "F"}
        expected_result = ({"WTS", "WTK", "WTB"}, 0.5)
        self.assertEqual(get_codon_for_amino_acids(amino_acids), expected_result)

    def test_truncate_list_of_amino_acids(self):

        # (1) test for empty set
        amino_acids = set()
        expected_result = set()
        self.assertEqual(truncate_list_of_amino_acids(amino_acids), expected_result)

        # (2) test for valid amino_acids
        amino_acids = {"A", "V", "I"}
        expected_result = {frozenset({"V", "A"}), frozenset({"V", "I"})}
        self.assertEqual(truncate_list_of_amino_acids(amino_acids), expected_result)

    """
    End test
    """


if __name__ == "__main__":
    unittest.main()
