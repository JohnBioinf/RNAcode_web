#!/usr/bin/python3
"""Select an optimal set of sequence to use for RNAcode."""
import sys
import os

from Bio import SeqIO

import SeqSelection
from SeqSelection import collect_blast_res
from SeqSelection import list_candidate_seq
from SeqSelection import check_candidates
from SeqSelection import write_candidates
from SeqSelection import eprint

# paths
REFERENCE_SPECIES_FILE_PATH = "reference_species.txt"

REFERENCE_SPECIES = {}


def get_arguments():
    """Get arguments from system."""
    if len(sys.argv) != 5:
        eprint("4 Arguments needed " + str(len(sys.argv) - 1) + " given")
        sys.exit(1)

    with open(SeqSelection.SELECTED_SPECIES_FILE_PATH, "r", encoding="UTF-8") as file_handle:
        try:
            SeqSelection.SELECTED_SPECIES = [int(line) for line in file_handle.read()]
        except ValueError:
            SeqSelection.ELECTED_SPECIES = []

    SeqSelection.MIN_PAIR_DIST = float(sys.argv[1])
    SeqSelection.MAX_PAIR_DIST = float(sys.argv[2])
    SeqSelection.ITERATION = int(sys.argv[3])
    SeqSelection.BLAST_DB = sys.argv[4]

    input_seq = list(SeqIO.parse(SeqSelection.INPUT_FILE_PATH, "fasta"))[0].seq
    SeqSelection.INPUT_SEQ = input_seq

    selected_sequences_path = SeqSelection.SELECTED_SEQUENCES_PATH.format(SeqSelection.ITERATION)

    SeqSelection.FINAL_CANDIDATES = {}
    if not os.path.isfile(selected_sequences_path):
        SeqSelection.FINAL_CANDIDATES["Target"] = str(input_seq)
    else:
        with open(selected_sequences_path, "r", encoding="UTF-8") as file_handle:
            for line in file_handle:
                if line[0] == ">":
                    name = line[1:-1]
                    SeqSelection.FINAL_CANDIDATES[name] = ""
                else:
                    SeqSelection.FINAL_CANDIDATES[name] += line[:-1]

    SeqSelection.NCBI = SeqSelection.load_ncbi()


def main():
    """Run main function."""
    global REFERENCE_SPECIES
    get_arguments()

    print("Collect blast results")
    blast_results = collect_blast_res()
    if len(blast_results) == 0:
        print("No candidates to add. No blast results.")
        sys.exit()

    # set reference species
    if not os.path.isfile(REFERENCE_SPECIES_FILE_PATH):
        with open(REFERENCE_SPECIES_FILE_PATH, "w", encoding="UTF-8") as file_handle:
            file_handle.write(str(blast_results[0].staxid))

    print("Get full sequences for candidates")
    list_candidate_seq(blast_results)

    print("Select candidates")
    check_candidates(blast_results)

    # save newly selected species
    if SeqSelection.SELECTED_SPECIES != []:
        with open(SeqSelection.SELECTED_SPECIES_FILE_PATH, "w", encoding="UTF-8") as file_handle:
            file_handle.write("\n".join([str(species) for species in SeqSelection.SELECTED_SPECIES]))

    print("Write candidates")
    write_candidates()

    print("Finished!")


if __name__ == "__main__":
    main()
