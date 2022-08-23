"""Select an optimal set of sequence to use for RNAcode."""
import sys
import re
import subprocess
from time import sleep
import random
import sqlite3
import os

from ete3 import NCBITaxa
from Bio import SeqIO, Align


# paths
SELECTED_SPECIES_FILE_PATH = "selected_species.txt"
BLAST_RESULT_PATH = "{}_blastn.result"
CANDIDATES_FASTA_PATH = "candidates.fasta"
SELECTED_SEQUENCES_PATH = "{}_selected_sequences.fasta"
POSITION_LIST_PATH = "positions_list_blastdbcmd.txt"
INPUT_FILE_PATH = "input.fasta"

DB_SIZE = 300

# Variables for analysis
MAX_NUM_SEQ = 25
MIN_NUM_SEQ = 2
FINAL_CANDIDATES = {}
ITERATION = None
EXPAND_TARGET = None
INPUT_SEQ = None
SELECTED_SPECIES = []
BLAST_DB = None

# User defined variables
MIN_PAIR_DIST = None
MAX_PAIR_DIST = None

# Must be loaded with load_ncbi() like so
# SeqSelection.NCBI = SeqSelection.load_ncbi()
NCBI = None

VERBOSE = True

ALIGNER = Align.PairwiseAligner()
ALIGNER.mode = "local"
ALIGNER.match_score = 1
ALIGNER.mismatch_score = -1
ALIGNER.open_gap_score = -1
ALIGNER.extend_gap_score = -1


def load_ncbi():
    """Load NCBI class from ete3 and manage sql error.

    The loop and try, except are needed if another NCBITaxa instance is updating
    the database. This could be avoided if a "no update" argument would exist
    for the initiation of the class NCBITaxa. I made a pull request if it will
    get accepted this function would be obsolete.
    """
    for i in range(10):
        try:
            ncbi = NCBITaxa()
            break
        except sqlite3.OperationalError:
            if i == 9:
                raise
            sleep(60)
    return ncbi


def vprint(*a, **k):
    """Print to verbose."""
    if VERBOSE:
        print(*a, **k, flush=True)


def eprint(*a, **k):
    """Print to stderr."""
    print(*a, file=sys.stderr, flush=True, **k)


class BlastResult:
    """Represents blast result.

    Expect blast result in specific format:
    outfmt = "6 evalue bitscore pident sseqid slen sstart send sframe qstart qend qframe qseq sseq staxid"
    """

    def __init__(self, result_line):
        """One line as input."""
        lines = result_line.split("\t")
        self.evalue = float(lines[0])
        self.bitscore = float(lines[1])
        self.pident = float(lines[2])
        self.sseqid = lines[3]
        self.slen = int(lines[4])
        self.sstart = int(lines[5])
        self.send = int(lines[6])
        self.sframe = int(lines[7])
        self.qstart = int(lines[8])
        self.qend = int(lines[9])
        self.qframe = int(lines[10])
        self.qseq = lines[11]
        self.sseq = lines[12]
        self.staxid = int(lines[13].strip())

        try:
            self.sacc = self.sseqid.split("|")[3]
        except IndexError:
            try:
                self.sacc = self.sseqid.split("|")[1]
            except IndexError:
                self.sacc = self.sseqid
        try:
            self.sgi = self.sseqid.split("|")[1]
        except IndexError:
            self.sgi = None

    def __str__(self):
        """Show some general results of the hit."""
        return (
            f"Blast result {self.sseqid} from {self.staxid} txid with "
            "{self.evalue} e-val, {self.bitscore} bitscore"
        )

    def fasta_header(self):
        """Nice format for fasta."""
        species_name = NCBI.get_taxid_translator([self.staxid])[self.staxid]
        species_name = "_".join(species_name.split(" ")[0:2])
        return (
            f"{species_name}-{self.sacc}-{self.sstart}_{self.send} "
            f"taxid:{self.staxid} e-value:{self.evalue}"
        )


def get_mol_type(gi_list, retry=True):
    """Get the molecular type for a accession list."""
    vprint(f"Get mol type for {len(gi_list)} gis.")
    is_dna_dic = {}

    if len(gi_list) > 5000:
        vprint("To big -> recursive.")
        is_dna_dic.update(get_mol_type(gi_list[5000:]))
        gi_list = gi_list[:5000]

    call_str = f"esummary -db nuccore -id {','.join(gi_list)} | xtract -pattern DocumentSummary -element Gi Biomol"

    try:
        completed_process = subprocess.run(
            call_str, capture_output=True, text=True, check=True, shell=True
        )
    except subprocess.CalledProcessError as exc:
        eprint("Entrez search failed!")
        eprint(call_str)
        eprint(exc.stdout)
        eprint(exc.stderr)
        raise

    for line in completed_process.stdout.split("\n")[:-1]:
        line = line.split()
        is_dna_dic[line[0]] = line[1] == "genomic"

    not_found_gi = [gi for gi in gi_list if gi not in is_dna_dic.keys()]
    if len(not_found_gi) != 0 and retry:
        vprint("Error entrez search. {len(not_found_gi)} are missing. Retry")
        for i in range(4):
            vprint(f"{i}. retry")
            is_dna_dic.update(get_mol_type(not_found_gi, retry=False))
            not_found_gi = [gi for gi in gi_list if gi not in is_dna_dic.keys()]
            if len(not_found_gi) == 0:
                break
            sleep(random.choice(range(80, 200)))

        if len(not_found_gi) != 0:
            eprint("Entrez search failed!")
            eprint("The following Gis could not be found by entrez:")
            eprint("\n".join(not_found_gi))
            raise subprocess.CalledProcessError

    return is_dna_dic


def collect_blast_res(select_genomic_seq=True):
    """Build list of  BlastResult objects.

    expect blast result in specific format:
    outfmt = "6 evalue bitscore pident sallseqid sstart send qstart qend qseq sseq staxid"
    """
    blast_results = []
    gi_list = []
    with open(
        BLAST_RESULT_PATH.format(ITERATION), "r", encoding="UTF-8"
    ) as file_handle:
        for line in file_handle:
            blast_result = BlastResult(line)
            gi_list.append(blast_result.sgi)
            blast_results.append(blast_result)

    # If the nt data base is used, each blast hit should be checked if the
    # sequence is a genomic region.
    if select_genomic_seq and BLAST_DB == "nt":
        vprint("Entrez call")
        is_dna_dic = get_mol_type(gi_list)

        blast_results = [result for result in blast_results if is_dna_dic[result.sgi]]
    return blast_results


def list_candidate_seq(blast_results):
    """Filter blast result based on lineage and sequence similarity.

    Sorts by species and bitscore. From this list to calls all with blastdbcmd
    to get full sequence length.
    """
    position_list = []
    len_input_seq = len(INPUT_SEQ)
    for blast_result in blast_results:
        if blast_result.qframe != 1:
            eprint("Query negative strand!")
            eprint(blast_result)
            continue
        if blast_result.staxid in SELECTED_SPECIES:
            continue
        len_query = blast_result.qend - blast_result.qstart + 1
        # Check if the query covers the input completely and distance is in
        # not in boundaries
        if len_query == len_input_seq and (
            1 - blast_result.pident < MIN_PAIR_DIST
            or 1 - blast_result.pident > MAX_PAIR_DIST
        ):
            continue
        accession = blast_result.sacc

        if blast_result.sframe == 1:
            full_sstart = blast_result.sstart - blast_result.qstart
            full_send = blast_result.send + (len_input_seq - blast_result.qend)
            # reduce or expand because of loss or gain through gaps
            gap_loss = (full_send - full_sstart) - len_input_seq
            full_sstart = full_sstart + int(gap_loss / 2)
            full_send = full_send - int(gap_loss / 2)
            full_sstart = max(full_sstart, 1)
            if full_send > blast_result.slen:
                full_send = blast_result.slen
            position = f"{accession} {full_sstart}-{full_send} plus"
        else:
            full_sstart = blast_result.send - blast_result.qstart
            full_send = blast_result.sstart + (len_input_seq - blast_result.qend)
            # reduce or expand because of loss or gain through gaps
            gap_loss = (full_send - full_sstart) - len_input_seq
            full_sstart = full_sstart + int(gap_loss / 2)
            full_send = full_send - int(gap_loss / 2)
            full_sstart = max(full_sstart, 1)
            if full_send > blast_result.slen:
                full_send = blast_result.slen
            position = f"{accession} {full_sstart}-{full_send} minus"
        position_list.append(position)

    call_blastdbcmd(position_list, CANDIDATES_FASTA_PATH)


def call_blastdbcmd(position_list, out_file):
    """Call blastdbcmd."""
    with open(POSITION_LIST_PATH, "w", encoding="UTF-8") as file_handle:
        for position in position_list:
            file_handle.write(position + "\n")

    call_str = (
        f'blastdbcmd -entry_batch {POSITION_LIST_PATH} -out {out_file} -db "{BLAST_DB}"'
    )
    vprint(call_str)
    try:
        subprocess.run(
            call_str, capture_output=True, text=True, check=True, shell=True
        )
    except subprocess.CalledProcessError as exc:
        eprint("Blastdbcmd failed!")
        eprint(exc.stdout)
        eprint(exc.stderr)
        raise


def print_alignment(alignment, as_error=False):
    """Print alignment generated by biopython function."""
    alignment = str(alignment).split("\n")
    line_size = 150
    for i in range(0, len(alignment[0]), line_size):
        if as_error:
            eprint(alignment[0][i : i + line_size])
            eprint(alignment[1][i : i + line_size])
            eprint(alignment[2][i : i + line_size])
        else:
            vprint(alignment[0][i : i + line_size])
            vprint(alignment[1][i : i + line_size])
            vprint(alignment[2][i : i + line_size])


def get_dist_alg(alignment, len_target):
    """Calculate percentage of gaps in alignment."""
    alignment_case = str(alignment).split("\n")[1]
    matches = alignment_case.count("|")
    return (1 - round(matches / len_target, 2)) * 100


def get_dist_seq(seq_a, seq_b):
    """Return sequence distance between two sequences."""
    alignment = ALIGNER.align(seq_a, seq_b)[0]
    len_seq = len(seq_a) if len(seq_a) < len(seq_b) else len(seq_b)
    return get_dist_alg(alignment, len_seq)


def build_dist_seq_dic(candidates):
    """Build distance matrix for candidates.

    As distance is symmetric only compute ones. And keys are ordered list of to
    integer so accessing is easy.
    """
    dist_seq_dic = {}
    for i, candidate_i in enumerate(candidates):
        for j, candidate_j in enumerate(candidates):
            key = frozenset([j, i])
            if key in dist_seq_dic:
                continue
            if j == i:
                dist = 0
            else:
                dist = get_dist_seq(candidate_i[1], candidate_j[1])
            dist_seq_dic[key] = dist
    return dist_seq_dic


def reduce_cand_min_dist(candidates):
    """Reduce candidates.

    Such that no two sequence have a distance smaller than MIN_PAIR_DIST. Build
    clusters based on MIN_PAIR_DIST. So that any clusters contains maximal many
    candidates.
    """
    # Build pairwise distance matrix
    dist_seq_dic = build_dist_seq_dic(candidates)

    # result list of indices to be deleted
    del_cand = []
    # label for skip: Already clustered
    clustered_candidates = []
    # Do clustering as long until all candidates are considered
    while len(clustered_candidates) < len(candidates):
        for i in range(0, len(candidates)):
            if i in clustered_candidates:
                continue
            # build cluster all in distance
            neighbourhood = [
                j
                for j in range(0, len(candidates))
                if dist_seq_dic[frozenset([i, j])] <= MIN_PAIR_DIST
                and j not in clustered_candidates
            ]
            # if cluster only contains itself continue
            if len(neighbourhood) == 1:
                clustered_candidates.append(i)
                continue
            # see if a bigger cluster could be build
            centroid = i
            for j in neighbourhood:
                if j == i:
                    continue
                new_neighbourhood = [
                    k
                    for k in range(0, len(candidates))
                    if dist_seq_dic[frozenset([j, k])] <= MIN_PAIR_DIST
                    and k not in clustered_candidates
                ]
                if len(new_neighbourhood) > len(neighbourhood):
                    centroid = j
                    neighbourhood = new_neighbourhood
            clustered_candidates += neighbourhood
            # N is no biggest cluster within min pair dist
            # remove all candidates from cluster but keep centroid
            neighbourhood.remove(centroid)
            del_cand += neighbourhood
    for i in sorted(del_cand, reverse=True):
        del candidates[i]
    return candidates


def calculate_mean_pairwise_dis(sequences):
    """Calculate the mean pairwise distance for a self of sequences."""
    the_sum = 0
    count = 0
    for i, seq_i in enumerate(sequences):
        for j, seq_j in enumerate(sequences):
            if j == i:
                continue
            count += 1
            the_sum += get_dist_seq(seq_i, seq_j)
    return the_sum / count


def get_species_taxid(taxid):
    """Get the species taxid for a specific taxid.

    If the taxid is above species level return None.
    """
    try:
        rank = NCBI.get_rank([taxid])[taxid]
    except KeyError:
        return None
    if rank == "species":
        species_taxid = taxid
    else:
        species_taxid = [
            txd
            for txd, rank in NCBI.get_rank(NCBI.get_lineage(taxid)).items()
            if rank == "species"
        ]
        try:
            species_taxid = species_taxid[0]
        except IndexError:
            species_taxid = None
    return species_taxid


def add_candidates(candidates, select_species=True):
    """Add candidates sorted by blast evalue."""
    vprint("Add candidates")
    candidates = sorted(candidates, key=lambda x: float(x[0].evalue))
    if select_species:
        candidates = [
            cand
            for cand in candidates
            if get_species_taxid(cand[0].staxid) not in SELECTED_SPECIES
        ]
    while len(candidates) > 0 and len(FINAL_CANDIDATES) < MAX_NUM_SEQ:
        name = candidates[0][0].fasta_header()
        if select_species:
            species_taxid = get_species_taxid(candidates[0][0].staxid)
        else:
            species_taxid = None
        vprint(f"Found best candidate {name}")
        FINAL_CANDIDATES[name] = candidates[0][1]
        del candidates[0]
        if species_taxid:
            SELECTED_SPECIES.append(species_taxid)
            candidates = [
                cand
                for cand in candidates
                if get_species_taxid(cand[0].staxid) not in SELECTED_SPECIES
            ]


def get_min_max_dist_to_final_set(candidate_seq):
    """Calculate distances to final set.

    Return the minimal distance to any sequence in the final set and the
    distance to the target.
    """
    max_dist_target = get_dist_seq(str(candidate_seq.seq), FINAL_CANDIDATES["Target"])
    min_dist = float("inf")
    for _final_candidate_name, final_candidate_seq in FINAL_CANDIDATES.items():
        dist = get_dist_seq(str(candidate_seq.seq), final_candidate_seq)
        min_dist = min(min_dist, dist)
    return min_dist, max_dist_target


def check_candidates(blast_results, select_species=True):
    """Map full sequence from blasdbcmd to candidates.

    If expand target option aligners with target seq and checks if target is
    present. Removes candidates based on max dist. Reduces set based on min
    dist. Adds candidates to global variable FINAL_CANDIDATES.
    """
    candidates = []
    candidate = None
    vprint("Collecting candidates")
    for candidate_seq in SeqIO.parse(CANDIDATES_FASTA_PATH, "fasta"):
        for candidate in blast_results:
            if candidate.sacc in candidate_seq.description:
                break
        else:
            print("Can not map candidate sequence to candidate")
            print(candidate_seq.description)
            continue

        if select_species and get_species_taxid(candidate.staxid) in SELECTED_SPECIES:
            continue

        # Check if candidate sequence is allowed by min max distance.
        min_dist, max_dist_to_target = get_min_max_dist_to_final_set(candidate_seq)

        if max_dist_to_target > MAX_PAIR_DIST or min_dist < MIN_PAIR_DIST:
            continue
        candidates.append([candidate, str(candidate_seq.seq)])

        # Check if already enough candidates had been added
        if len(candidates) > 200:
            vprint("Check if enough candidates were collected")
            candidates = reduce_cand_min_dist(candidates)
            vprint(f"{len(candidates)} candidates after reduction")

            if len(candidates) == 0:
                vprint("No candidates to add. After reduction.")
                candidates = []
                continue
            add_candidates(candidates, select_species=select_species)
            if len(FINAL_CANDIDATES) == MAX_NUM_SEQ:
                vprint("Enough candidates found no need to add more")
                return
            candidates = []

    vprint(f"{len(candidates)} candidates collected")

    if len(candidates) == 0:
        print("No candidates to add. Before reduction.")
        return
    # First reduce list so that each candidate only has minimum distance to any
    # other candidate
    vprint("Reducing candidates")
    candidates = reduce_cand_min_dist(candidates)
    vprint(f"{len(candidates)} candidates after reduction")

    if len(candidates) == 0:
        print("No candidates to add. After reduction.")
        return
    # Second add candidates until max is reached or candidates are exhausted
    vprint("Add candidates")
    add_candidates(candidates, select_species=select_species)


def write_candidates():
    """Write the final multiple fasta."""
    line_length = 60
    with open(
        SELECTED_SEQUENCES_PATH.format(ITERATION), "w", encoding="UTF-8"
    ) as file_handle:
        file_handle.write(">Target\n")
        seq = FINAL_CANDIDATES["Target"]
        print_seq = "\n".join(re.findall(f".{{1,{line_length}}}", seq)) + "\n"
        file_handle.write(print_seq)
        for name, seq in FINAL_CANDIDATES.items():
            if name == "Target":
                continue
            print_seq = "\n".join(re.findall(f".{{1,{line_length}}}", seq)) + "\n"
            file_handle.write(f">{name}\n")
            file_handle.write(print_seq)


def collect_sequences(iteration, current_work_dir="."):
    """Collect all previously selected sequences."""
    sequences = {}
    for i in range(1, iteration + 1):
        selected_sequences_path = (
            f"{current_work_dir}/{SELECTED_SEQUENCES_PATH.format(i)}"
        )
        if not os.path.isfile(selected_sequences_path):
            continue
        key = None
        current_sequences = {}
        with open(selected_sequences_path, "r", encoding="UTF-8") as file_handle:
            for line in file_handle:
                if line[0] == ">":
                    key = re.sub(r":c*[0-9]+-[0-9]+", "", line[1:].split()[0])
                    current_sequences[key] = ""
                else:
                    current_sequences[key] += line.strip()
        sequences.update(current_sequences)
    return sequences
