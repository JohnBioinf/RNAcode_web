"""Select an optimal set of sequence to use for RNAcode based on a custom blast db."""

import sys

from Bio import SeqIO

import SeqSelection
from SeqSelection import eprint
from SeqSelection import collect_blast_res
from SeqSelection import list_candidate_seq
from SeqSelection import check_candidates
from SeqSelection import write_candidates


SeqSelection.VERBOSE = True


def get_arguments():
    """Get arguments from system."""
    if len(sys.argv) != 4:
        eprint("3 Arguments needed " + str(len(sys.argv) - 1) + " given")
        sys.exit(1)

    min_pair_dist = float(sys.argv[1])
    max_pair_dist = float(sys.argv[2])
    blast_db = sys.argv[3]

    SeqSelection.MIN_PAIR_DIST = min_pair_dist
    SeqSelection.MAX_PAIR_DIST = max_pair_dist

    input_seq = str(list(SeqIO.parse(SeqSelection.INPUT_FILE_PATH, "fasta"))[0].seq)
    SeqSelection.INPUT_SEQ = input_seq

    final_candidates = {}
    final_candidates["Target"] = input_seq
    SeqSelection.FINAL_CANDIDATES = final_candidates

    SeqSelection.BLAST_DB = blast_db
    SeqSelection.ITERATION = 1
    SeqSelection.NCBI = SeqSelection.load_ncbi()


def main():
    """Run main function."""
    get_arguments()

    print("Collect blast results")
    blast_results = collect_blast_res(select_genomic_seq=False)
    if len(blast_results) == 0:
        print("No candidates to add. No blast results.")
        sys.exit()

    print("Get full sequences for candidates")
    list_candidate_seq(blast_results)

    print("Select candidates")
    check_candidates(blast_results, select_species=False)

    print("Write candidates")
    write_candidates()

    print("Finished!")


if __name__ == "__main__":
    main()
