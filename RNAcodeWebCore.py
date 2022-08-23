"""
Functions and general parameters for back end service for RNAcode web.

Handles all computation etc.

Shared between:
RNAcodeWeb_build_db
...

"""

import json
import sys
import os
import urllib.request
import subprocess
import re
import time
import traceback
from glob import glob
from collections import defaultdict

from ete3 import Tree

import jinja2

import SeqSelection
from SeqSelection import MIN_NUM_SEQ
from SeqSelection import collect_sequences

with open("./parameters_backend_local.json", "r", encoding="UTF-8") as file_handle:
    PARAMETERS_BACKEND = json.load(file_handle)

with open("./parameters_frontend_local.json", "r", encoding="UTF-8") as file_handle:
    PARAMETERS_FRONTEND = json.load(file_handle)

PORT_FRONTEND = PARAMETERS_FRONTEND["port"]
IP_FRONTEND = PARAMETERS_FRONTEND["ip"]

PATH = ":".join(
    [
        PARAMETERS_BACKEND["blast_bin"],
        PARAMETERS_BACKEND["RNAcode_web_repo_backend"],
        PARAMETERS_BACKEND["clustal_bin"],
        PARAMETERS_BACKEND["RNAcode_bin"],
        PARAMETERS_BACKEND["entrez_bin"],
        PARAMETERS_BACKEND["PATH"],
    ]
)

PYTHON_ENV_SEQSEL_PATH = PARAMETERS_BACKEND["python_env_seqSel_path"]
BLAST_DB_PATH = PARAMETERS_BACKEND["blast_db"]
WORK_DIR = PARAMETERS_BACKEND["work_dir"]

# Depending on number of idle cores (range) allow n cores. For blast
CORE_USE_DIC = {
    range(0, 11): 1,
    range(11, 21): 2,
    range(21, 31): 3,
    range(31, 41): 4,
    range(41, 161): 5,
}

P_THRESHOLD = 0.05

# Should be set in sub modules in get_arguments()
JOB_ID = None
GENOME_START = None
DB_TYPE = None
BLAST_DB = None
VERBOSE = True
REFERENCE_SPECIES = None
ANCESTORS = []
NCBI = None
ALL_TAXIDS_BLAST_DB = None

# general structure
STDOUT_FILE_PATH = "./logs/{}.out"
STDERR_FILE_PATH = "./logs/{}.err"

CURRENT_WORK_DIR_TEMPLATE = WORK_DIR + "/{}"
SLURM_OUTPUT_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/{}.out"
SLURM_ERROR_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/{}.err"
JOB_STATUS_FILE_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/job_status.json"
README_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/README.md"
RUN_ALL_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/run_all.sh"
PYTHON_REQ_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/python_requirements.txt"
REFERENCE_SPECIES_FILE_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/reference_species.txt"
SELECTED_SPECIES_FILE_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/selected_species.txt"
# blastn
BLAST_RESULT_PATH_TEMPLATE = (
    CURRENT_WORK_DIR_TEMPLATE + "/" + SeqSelection.BLAST_RESULT_PATH
)
INPUT_FILE_PATH_TEMPLATE = (
    CURRENT_WORK_DIR_TEMPLATE + "/" + SeqSelection.INPUT_FILE_PATH
)
TAXIDS_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/{}_taxid_list.txt"
BLAST_WRAPPER_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/{}_blastn.sh"
GI_LIST_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/{}_gilist.txt"
# sequence selection
SEQ_SELECTION_SCRIPT_PATH_TEMPLATE = (
    CURRENT_WORK_DIR_TEMPLATE + "/SeqSelection_pipeline.py"
)
SEQ_SELECTION_MODULE_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/SeqSelection.py"
SEQSEL_WRAPPER_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/{}_seqSel.sh"
SELECTED_SEQUENCES_PATH_TEMPLATE = (
    CURRENT_WORK_DIR_TEMPLATE + "/" + SeqSelection.SELECTED_SEQUENCES_PATH
)
SELECTED_SPECIES_FILE_PATH_TEMPLATE = (
    CURRENT_WORK_DIR_TEMPLATE + "/" + SeqSelection.SELECTED_SPECIES_FILE_PATH
)
# alignment
ALIGN_SCRIPT_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/alignment.sh"
ALIGN_FASTA_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/alignment.fasta"
ALIGN_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/alignment.aln"
ALIGN_MAF_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/alignment.maf"
ALIGN_TREE_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/alignment.tree"
ALIGN_PLOT_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/alignment.svg"
# RNAcode
RNACODE_SCRIPT_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/RNAcode.sh"
RNACODE_RESULT_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/rnacode_result.tsv"
EPS_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/eps"
# Build DB
# This is only path that is relative not absolute
CUSTOM_DB_PATH = "blast_db"
BLAST_DB_FASTA_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/blast_db.fasta"
BUILD_DB_WRAPPER_PATH_TEMPLATE = CURRENT_WORK_DIR_TEMPLATE + "/build_db.sh"

OUTFMT = (
    "6 evalue bitscore pident sseqid slen sstart send sframe qstart "
    "qend qframe qseq sseq staxid"
)

SEQSEL_WRAPPER_TEMPLATE = (
    "#!/bin/bash\n\nset -e\n\n"
    "source $PYTHON_ENV/bin/activate\n\n"
    "python3 seq_selection.py {} {} {} {} {} {}\n"
)

CLUSTALO_WRAPPER_TEMPLATE = (
    "#!/bin/bash\n\n"
    "set -e\n\n"
    "clustalo-1.2.4-Ubuntu-x86_64 -i {} -o {} "
    "--output-order=input-order --outfmt=clu --force "
    "--guidetree-out={}\n"
)

RNACODE_SCRIPT_TEMPLATE = (
    "#!/bin/bash\n\n"
    f"set -e\n\nRNAcode -o {{}} -t {{}} -p {P_THRESHOLD} -e --eps-dir {{}}\n"
)


class PipelineFailed(Exception):
    """Exception if pipeline failed before finishing completly."""


class PipelineError(Exception):
    """Exception if pipeline failed before finishing completly."""


def get_number_idles_cpus():
    """Return number of idle cores."""
    call_str = "sinfo -a --format='%C'"
    try:
        completed_process = subprocess.run(call_str.split(), capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        eprint("Error determining free cores")
        eprint(call_str)
        eprint(exc.stdout)
        eprint(exc.stderr)
        raise
    return int(completed_process.stdout.split("\n")[1].split("/")[1])


def start_pipeline(main):
    """Start pipeline and handle errors."""
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
    except PipelineFailed:
        eprint("Pipeline failed.")
        notify_frontend(JOB_ID)
        sys.exit(0)
    except PipelineError:
        error = traceback.format_exc()
        eprint("Pipeline error.")
        write_status_file(f"fullJob.{JOB_ID}", ["E", "pipe_broken"])
        eprint(error)
        notify_frontend(JOB_ID)
        sys.exit(1)
    except Exception:
        error = traceback.format_exc()
        eprint("Pipeline unexpected error.")
        eprint(error)
        write_status_file(f"fullJob.{JOB_ID}", ["E", "pipe_broken"])
        notify_frontend(JOB_ID)
        sys.exit(1)


def notify_frontend(job_id):
    """Notify frontend that job finished."""
    if PORT_FRONTEND == "None":
        url = f"http://{IP_FRONTEND}/submission/{job_id}"
    else:
        url = f"http://{IP_FRONTEND}:{PORT_FRONTEND}/submission/{job_id}"
    for i in range(1, 4):
        try:
            urllib.request.urlopen(url)
            break
        except urllib.error.URLError:
            if i == 3:
                eprint("Error notifiying frontend")
                eprint(url)
                raise
            time.sleep(4)
            continue


def eprint(*a, **k):
    """Print to stderr."""
    if not STDERR_FILE_PATH:
        print(*a, file=sys.stderr, flush=True, **k)
    else:
        with open(STDERR_FILE_PATH, "a", encoding="UTF-8") as f_handle:
            print(*a, file=f_handle, flush=True, **k)


def check_process(job_name):
    """Check how process ended."""
    job_type = job_name.split(".")[0]
    vprint(f"Job {job_type} ended")
    error_file = SLURM_ERROR_TEMPLATE.format(JOB_ID, job_type)
    normal_jobs = ["blastn", "alignment", "seqSel", "RNAcode", "buildDB"]

    with open(error_file, "r", encoding="UTF-8") as f_handle:
        error_content = f_handle.read()
        # weird slurm error catch seem to be have no influence
        error_content = "\n".join(
            [
                line
                for line in error_content.split("\n")
                if "error: _is_a_lwp: " not in line
            ]
        )

    if job_type == "RNAcode":
        rnacode_result_path = RNACODE_RESULT_PATH_TEMPLATE.format(JOB_ID)
        if not os.path.isfile(rnacode_result_path):
            write_status_file(job_name, ["F", 1])
            raise PipelineError

        with open(rnacode_result_path, "r", encoding="UTF-8") as f_handle:
            content = f_handle.read()
        if content != "":
            write_status_file(job_name, ["CD", 0])
            return True
        write_status_file(job_name, ["F", "no_res"])
        write_status_file(f"fullJob.{JOB_ID}", ["F", "no_res"])
        raise PipelineFailed

    if any(job in job_type for job in normal_jobs):
        if error_content == "":
            write_status_file(job_name, ["CD", 0])
            return True

        vprint(error_content)
        if "DUE TO TIME LIMIT" in error_content:
            write_status_file(job_name, ["F", "time"])
            raise PipelineError
        # this should catch the oom killer, but this is fairly generic. Yet
        # SLURM does not give any better option at the moment.
        if "Killed" in error_content:
            write_status_file(job_name, ["F", "mem"])
            raise PipelineError

        write_status_file(job_name, ["F", 1])
        raise PipelineError

    eprint("Unknown jobtype: " + job_type)
    eprint("check_process()")
    raise PipelineError


def slurm_listen(job_name):
    """Listen to SLURM batch and updating status_file.

    Afterwards check output of script with check_process().
    """
    call_str = f"squeue --name={job_name} -h"
    while True:
        try:
            completed_process = subprocess.run(
                call_str.split(), capture_output=True, text=True, check=True
            )
        except subprocess.CalledProcessError as exc:
            eprint("Error listening to slurm job")
            eprint(call_str)
            eprint(exc.stdout)
            eprint(exc.stderr)
            write_status_file(job_name, ["F", 1])
            raise
        if completed_process.stdout == "":
            break
        status = re.split(" +", completed_process.stdout)[5]
        write_status_file(job_name, [status, 0])
        time.sleep(30)

    time.sleep(4)
    return check_process(job_name)


def slurm_batch(script_path, job_type, cores=1):
    """Submit SLURM job, listens to job with slurm_listen().

    Checks if job finished normally with check_process().
    """
    job_name = job_type + "." + JOB_ID
    current_work_dir = CURRENT_WORK_DIR_TEMPLATE.format(JOB_ID)
    output = SLURM_OUTPUT_TEMPLATE.format(JOB_ID, job_type)
    error = SLURM_ERROR_TEMPLATE.format(JOB_ID, job_type)
    jobs = ["blastn", "seqSel", "alignment", "RNAcode", "buildDB"]
    if any(job in job_type for job in jobs):
        max_time = str(10 * 60 * 60)
        partition = "main"
    else:
        eprint("Unknown jobtype: " + job_type)
        eprint("slurm_batch()")
        raise PipelineError

    call_str = (
        f"sbatch --job-name={job_name} --output={output} --error={error} "
        f"--time={max_time} --partition={partition} --chdir={current_work_dir} "
        f"--cpus-per-task={cores} "
        f"--export=ALL,PATH={PATH},BLASTDB={BLAST_DB_PATH},PYTHON_ENV={PYTHON_ENV_SEQSEL_PATH} "
        f"{script_path}"
    )
    vprint(call_str)

    for i in range(4):
        try:
            completed_process = subprocess.run(call_str.split(), capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as exc:
            eprint("error submitting slurm job")
            eprint(exc.stdout)
            eprint(exc.stderr)
            raise
        vprint(completed_process.stdout)
        time.sleep(3)
        try:
            return slurm_listen(job_name)
        # Sometimes slurm does not start the job
        except FileNotFoundError:
            vprint("Slurm did not execute the job. Try again.")
            if i == 3:
                raise
            time.sleep(3)
            continue


def write_status_file(job_name, status):
    """Write status file for job."""
    job_type = job_name.split(".")[0]
    job_status_file = JOB_STATUS_FILE_TEMPLATE.format(JOB_ID)
    try:
        with open(job_status_file, "r", encoding="UTF-8") as f_handle:
            content = f_handle.read()
            job_status = json.loads(content)
    except json.decoder.JSONDecodeError:
        eprint("JSONDDecoderError!")
        eprint(content)
        eprint(traceback.format_exc())
        raise PipelineError

    job_status[job_type] = status
    with open(job_status_file, "w", encoding="UTF-8") as f_handle:
        json.dump(job_status, f_handle, indent=4)


def vprint(*a, **k):
    """Print verbose."""
    if VERBOSE and not STDOUT_FILE_PATH:
        print(*a, **k, flush=True)
    elif VERBOSE and STDOUT_FILE_PATH:
        with open(STDOUT_FILE_PATH, "a", encoding="UTF-8") as f_handle:
            print(*a, file=f_handle, flush=True, **k)


def check_num_seq(iteration):
    """Check if enough sequence were found."""
    current_work_dir = CURRENT_WORK_DIR_TEMPLATE.format(JOB_ID)
    num_seqs = len(collect_sequences(iteration, current_work_dir=current_work_dir))
    if num_seqs - 1 >= MIN_NUM_SEQ:
        return

    vprint("Not enough sequences!")
    len_blast_result = 0
    for blast_result_files in glob(current_work_dir + "/*_blastn.result"):
        with open(blast_result_files, "r", encoding="UTF-8") as f_handle:
            len_blast_result += sum(1 for i in f_handle)

    if len_blast_result - 1 < MIN_NUM_SEQ:
        write_status_file(f"fullJob.{JOB_ID}", ["F", "blast"])
    else:
        write_status_file(f"fullJob.{JOB_ID}", ["F", "selection"])
    raise PipelineFailed


def concat_sequences_fasta(iteration):
    """Concatinate all selected sequences in one."""
    line_length = 60
    sequences = collect_sequences(iteration, current_work_dir=CURRENT_WORK_DIR_TEMPLATE.format(JOB_ID))
    align_fasta = ALIGN_FASTA_TEMPLATE.format(JOB_ID)
    with open(align_fasta, "w", encoding="UTF-8") as db_handle:
        for accession, sequence in sequences.items():
            print_seq = "\n".join(re.findall(f".{{1,{line_length}}}", sequence))
            db_handle.write(f">{accession}\n{print_seq}\n")


def clustalo():
    """3. Step in analysis. Align sequence previously selected."""
    job_type = "alignment"
    align_script_path = ALIGN_SCRIPT_PATH_TEMPLATE.format(JOB_ID)
    error_file = SLURM_ERROR_TEMPLATE.format(JOB_ID, job_type)
    align_fasta = ALIGN_FASTA_TEMPLATE.format(JOB_ID)
    align_path = ALIGN_PATH_TEMPLATE.format(JOB_ID)
    align_tree_path = ALIGN_TREE_PATH_TEMPLATE.format(JOB_ID)

    align_wrapper = CLUSTALO_WRAPPER_TEMPLATE.format(
        align_fasta.split("/")[-1],
        align_path.split("/")[-1],
        align_tree_path.split("/")[-1],
    )

    with open(align_script_path, "w", encoding="UTF-8") as f_handle:
        f_handle.write(align_wrapper)

    if not slurm_batch(align_script_path, job_type):
        eprint("Alignment exited with error!")
        with open(error_file, "r", encoding="UTF-8") as f_handle:
            eprint(f_handle.read())
        raise PipelineError


def convert_to_maf():
    """Convert the clustal alignment into a .maf file.

    Import to hold positional information of the parent.
    """
    align_path = ALIGN_PATH_TEMPLATE.format(JOB_ID)
    align_maf_path = ALIGN_MAF_PATH_TEMPLATE.format(JOB_ID)
    alignment_dic = defaultdict(str)
    maf_block = [["a", "score=0"]]

    with open(align_path, "r", encoding="UTF-8") as f_handle:
        for line in f_handle:
            if line[0:7] == "CLUSTAL":
                continue
            if line[0] in [" ", "", "\n"]:
                continue
            src = line.split()[0]
            text = line.split()[1]
            alignment_dic[src] += text

    for src, text in alignment_dic.items():
        if src == "Target":
            start = str(GENOME_START)
        else:
            start = "1"
        length_seq = str(len(text.replace("-", "")))
        maf_block.append(["s", src, start, length_seq, "+", length_seq, text])

    with open(align_maf_path, "w", encoding="UTF-8") as f_handle:
        f_handle.write("\n".join([" ".join(entry) for entry in maf_block]) + "\n")


def plot_alignment():
    """Plot alignment tree."""
    align_tree_path = ALIGN_TREE_PATH_TEMPLATE.format(JOB_ID)
    align_plot_path = ALIGN_PLOT_PATH_TEMPLATE.format(JOB_ID)
    with open(align_tree_path, "r", encoding="UTF-8") as f_handle:
        tree = Tree(f_handle.read())

    for node in tree:
        node.name = node.name.split("-")[0]

    tree.render(align_plot_path)


def rnacode():
    """Last step of analysis. Execute RNAcode on aligment."""
    job_type = "RNAcode"
    rnacode_script_path = RNACODE_SCRIPT_PATH_TEMPLATE.format(JOB_ID)
    error_file = SLURM_ERROR_TEMPLATE.format(JOB_ID, job_type)
    rnacode_result_path = RNACODE_RESULT_PATH_TEMPLATE.format(JOB_ID)
    align_path = ALIGN_MAF_PATH_TEMPLATE.format(JOB_ID)
    eps_path = EPS_PATH_TEMPLATE.format(JOB_ID)

    rnacode_script = RNACODE_SCRIPT_TEMPLATE.format(
        rnacode_result_path.split("/")[-1],
        align_path.split("/")[-1],
        eps_path.split("/")[-1],
    )

    with open(rnacode_script_path, "w", encoding="UTF-8") as f_handle:
        f_handle.write(rnacode_script)

    for _i in range(2):
        if slurm_batch(rnacode_script_path, job_type):
            error_content = ""
            break
        with open(error_file, "r", encoding="UTF-8") as f_handle:
            error_content = f_handle.read()
        if "Segmentation fault" not in error_content:
            break

    if error_content != "":
        eprint("RNAcode exited with error!")
        raise PipelineError


def render_jinja(file_name, context, template_dir="./"):
    """General function to render template file with jinja."""
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_dir + "/"))
    return env.get_template(file_name).render(context)


def get_version(program):
    """Get version for programm."""
    if program == "blastn":
        call_str = "blastn -version"
    elif program == "clustalo":
        call_str = "clustalo-1.2.4-Ubuntu-x86_64 --version"
    elif program == "RNAcode":
        call_str = "RNAcode --version"
    elif program == "blastdb":
        call_str = f"stat -c '%y' {BLAST_DB_PATH}/nt.00.nhd"
    else:
        eprint(f"Unknown program {program}")
        raise PipelineError

    current_env = os.environ.copy()
    current_env["PATH"] = PATH
    try:
        completed_process = subprocess.run(
            call_str.split(), capture_output=True, text=True, check=True, env=current_env
        )
    except subprocess.CalledProcessError as exc:
        eprint("get_version() ended with an error!")
        eprint(call_str)
        eprint(exc.stdout)
        eprint(exc.stderr)
        raise

    out = completed_process.stdout
    if program == "blastn":
        version = out.split("\n", maxsplit=1)[0].split(" ")[1]
    if program == "clustalo":
        version = out
    if program == "RNAcode":
        version = out.split()[2]
    if program == "blastdb":
        version = out.split()[0]
    return version


def get_kingdom(taxid):
    """Get the name of the kindom for a given taxid."""
    ancestors = NCBI.get_lineage(taxid)
    ancestors_rank = NCBI.get_rank(ancestors)
    try:
        taxid_kingdom = [t for t, g in ancestors_rank.items() if g == "superkingdom"][0]
    except IndexError:
        return None
    return NCBI.get_taxid_translator([taxid_kingdom])[taxid_kingdom]


def _set_reference_species(iteration):
    global REFERENCE_SPECIES, BLAST_DB, ALL_TAXIDS_BLAST_DB
    if REFERENCE_SPECIES is not None:
        return

    blast_result_path = BLAST_RESULT_PATH_TEMPLATE.format(JOB_ID, iteration - 1)
    with open(blast_result_path, "r", encoding="UTF-8") as f_handle:
        try:
            REFERENCE_SPECIES = [
                int(li.split()[13].strip()) for li in f_handle
            ][0]
        except IndexError:
            REFERENCE_SPECIES = None

    if REFERENCE_SPECIES is None:
        return

    # if the refseq DBs are used the exact one can be now specified
    if DB_TYPE == "refseq":
        kingdom = get_kingdom(REFERENCE_SPECIES)
        if kingdom == "Viruses":
            BLAST_DB = "ref_viroids_rep_genomes ref_viruses_rep_genomes"
        elif kingdom in ["Archaea", "Bacteria"]:
            BLAST_DB = "ref_prok_rep_genomes"
        elif kingdom == "Eukaryota":
            BLAST_DB = "ref_euk_rep_genomes"
        else:
            eprint(f"Unknown kingdom: {kingdom}")
            raise PipelineError

    for db in BLAST_DB.split(" "):
        with open(f"{BLAST_DB_PATH}/{db}.taxidlist", "r", encoding="UTF-8") as f_handle:
            ALL_TAXIDS_BLAST_DB = set(line.strip() for line in f_handle)


def _build_positive_taxid_list(iteration):
    """Build a list of all taxids which should be included in the next blast search.

    The taxonomic restricition is based on the iteration. The higher the
    iteration the more species should be included. All previously found species
    will be excluded from the next search, hence removed from the list. The
    function will return False if no rank ca be assigned to the reference.
    Further if the DB_TYPE is "refseq" each species must be removed which does
    not belong to the same kingdom as "REFERENCE_SPECIES"
    """
    # The taxonomic rank to which the search should be restricted. Depending on
    # the iteration.
    rank_level_dic = {
        2: "family",
        2.5: "suborder",
        3: "order",
        3.5: "subclass",
        4: "class",
        4.5: "subphylum",
        5: "phylum",
        5.5: "subkingdom",
        6: "kingdom",
    }
    taxids_path = TAXIDS_PATH_TEMPLATE.format(JOB_ID, iteration)
    selected_species_file_path = SELECTED_SPECIES_FILE_PATH_TEMPLATE.format(JOB_ID)

    ancestors = NCBI.get_lineage(REFERENCE_SPECIES)

    # try to find the taxid corresponding to the rank in rank_level_dic. If no
    # taxid can be found try to find a taxonomic higher.
    rank_add = 0
    while iteration + rank_add <= 6:
        rank = rank_level_dic[iteration]
        ancestors_rank = NCBI.get_rank(ancestors)
        try:
            ancestor_threshold = [t for t, g in ancestors_rank.items() if g == rank][0]
        except IndexError:
            rank_add += 0.5
        else:
            break
    else:
        return False

    call_str = f"{PARAMETERS_BACKEND['blast_bin']}/get_species_taxids.sh -t {ancestor_threshold} > {taxids_path}"
    try:
        subprocess.run(
            call_str, capture_output=True, text=True, shell=True, check=True
        )
    except subprocess.CalledProcessError as exc:
        eprint("get_species_taxids.sh ended with an error!")
        eprint(call_str)
        eprint(exc.stdout)
        eprint(exc.stderr)
        write_status_file(f"blastn_{iteration}", ["F", 1])
        raise

    vprint(call_str)

    # Remove from the taxid list all previously selected taxids
    with open(selected_species_file_path, "r", encoding="UTF-8") as f_handle:
        try:
            selected_species = set(line.strip() for line in f_handle)
        except ValueError:
            # if no species had been selected
            return True

    all_taxids = ALL_TAXIDS_BLAST_DB - selected_species

    with open(taxids_path, "r", encoding="UTF-8") as f_handle:
        taxid_list = set(
            line.strip()
            for line in f_handle
        )

    taxid_list = taxid_list.intersection(all_taxids)

    with open(taxids_path, "w", encoding="UTF-8") as f_handle:
        f_handle.write("\n".join(taxid_list) + "\n")
    return True


def build_taxid_list(iteration):
    """Generate a list of taxids that should be included in the blast search.

    The list contains all sequence that belong to a certain taxonomic rank. But
    excludes all previously found sequences and all taxids which are not present
    in the current data base.
    """
    taxids_path = TAXIDS_PATH_TEMPLATE.format(JOB_ID, iteration)
    with open(taxids_path, "w", encoding="UTF-8") as f_handle:
        f_handle.write("")
    # First iteration is always with empty restriction
    if iteration == 1:
        return
    # First check if a reference was set if not try to set
    vprint("Set reference species")
    _set_reference_species(iteration)
    # If no reference was set skip
    if REFERENCE_SPECIES is None:
        return

    # Build list for species taxid below a certain rank and remove from this
    # list all previously found taxids
    vprint("Build positive taxid list")
    if not _build_positive_taxid_list(iteration):
        return
