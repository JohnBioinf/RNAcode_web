"""Select an optimal set of sequence to build blastdb."""

import sys

from Bio import SeqIO

import SeqSelection
from SeqSelection import vprint
from SeqSelection import DB_SIZE
from SeqSelection import collect_sequences

SeqSelection.VERBOSE = True

TAXIDMAPFILE_PATH = "TaxIDMapFile"
# Factor for the input sequence which expanse the target sequence to both sides.
EXPANSION_FACTOR = 5


INPUT_SEQ = None
ITERATION = None
SELECTED_SEQUENCES_PATH = None
SELECTED_SPECIES = None
FOUND_REGIONS = None


def get_arguments():
    """Get arguments from system."""
    global INPUT_SEQ, ITERATION, SELECTED_SEQUENCES_PATH, SELECTED_SPECIES

    # Iteration is last blast call.
    ITERATION = int(sys.argv[1])
    SeqSelection.BLAST_DB = sys.argv[2]

    SeqSelection.ITERATION = ITERATION
    SeqSelection.NCBI = SeqSelection.load_ncbi()

    INPUT_SEQ = list(SeqIO.parse(SeqSelection.INPUT_FILE_PATH, "fasta"))[0].seq
    SELECTED_SEQUENCES_PATH = SeqSelection.SELECTED_SEQUENCES_PATH.format(ITERATION)

    with open(SeqSelection.SELECTED_SPECIES_FILE_PATH, "r", encoding="UTF-8") as file_handle:
        try:
            SELECTED_SPECIES = [int(line.strip()) for line in file_handle.readlines()]
        except ValueError:
            SELECTED_SPECIES = []


def select_regions(blast_results):
    """Make a list of postion for which sequences can be retrieved by blastdbcmd."""
    position_list = []
    taxidmapfile = []
    for blast_result in blast_results:
        accession = blast_result.sacc
        if accession in FOUND_REGIONS:
            continue
        species_taxid = SeqSelection.get_species_taxid(blast_result.staxid)
        if species_taxid in SELECTED_SPECIES:
            continue
        # species taxid can be None, do not add None to list
        if species_taxid:
            SELECTED_SPECIES.append(species_taxid)
        taxidmapfile.append(f"{accession} {blast_result.staxid}")

        if blast_result.sframe == 1:
            full_sstart = blast_result.sstart - len(INPUT_SEQ) * EXPANSION_FACTOR
            full_send = blast_result.send + len(INPUT_SEQ) * EXPANSION_FACTOR
            strand = "plus"
        else:
            full_sstart = blast_result.send - len(INPUT_SEQ) * EXPANSION_FACTOR
            full_send = blast_result.sstart + len(INPUT_SEQ) * EXPANSION_FACTOR
            strand = "minus"

        position_list.append(f"{accession} {full_sstart}-{full_send} {strand}")
        FOUND_REGIONS[accession] = ""
        if len(FOUND_REGIONS) >= DB_SIZE:
            break

    if len(position_list) == 0:
        print("No possible regions found.")
        sys.exit()

    with open(TAXIDMAPFILE_PATH, "a", encoding="UTF-8") as f_handle:
        f_handle.write("\n".join(taxidmapfile) + "\n")

    SeqSelection.call_blastdbcmd(position_list, SELECTED_SEQUENCES_PATH)


def main():
    """Run main function."""
    global FOUND_REGIONS
    get_arguments()
    FOUND_REGIONS = collect_sequences(ITERATION - 1)
    vprint(f"Iteration: {ITERATION}")
    vprint(f"{len(FOUND_REGIONS)} regions already found.")

    vprint("Collect blast results")
    blast_results = SeqSelection.collect_blast_res()
    if len(blast_results) == 0:
        print("No blast results.")
        sys.exit()

    vprint("Get full sequences for candidate regions")
    select_regions(blast_results)
    vprint(f"{len(FOUND_REGIONS)} regions found.")

    with open(SeqSelection.SELECTED_SPECIES_FILE_PATH, "w", encoding="UTF-8") as file_handle:
        file_handle.write("\n".join([str(species) for species in SELECTED_SPECIES]))

    vprint("Finished!")


if __name__ == "__main__":
    main()
