from assignment import get_codon_for_amino_acids, truncate_list_of_amino_acids, valid_aa
import time
import random

# USE THIS SCRIPT TO COMPUTE TIME TO RUN ALGORITHM >100 TIMES

st = time.time()

for i in range(100):
    amino_acids = set([random.choice(valid_aa) for x in range(3)])
    codons = get_codon_for_amino_acids(amino_acids)
    truncate_list = truncate_list_of_amino_acids(amino_acids)

et = time.time()
elapsed_time = et - st
print("Execution time:", elapsed_time, "seconds")
