"""Back end service for RNAcode web based on a custom blast db.

Handles all computation etc.
"""
import sys
import os
import re
import json
from shutil import rmtree, copyfile

import RNAcodeWebCore
from RNAcodeWebCore import start_pipeline
from RNAcodeWebCore import eprint
from RNAcodeWebCore import vprint
from RNAcodeWebCore import write_status_file
from RNAcodeWebCore import get_version
from RNAcodeWebCore import render_jinja
from RNAcodeWebCore import notify_frontend
from RNAcodeWebCore import check_num_seq
from RNAcodeWebCore import clustalo
from RNAcodeWebCore import convert_to_maf
from RNAcodeWebCore import rnacode
from RNAcodeWebCore import plot_alignment
from RNAcodeWebCore import build_taxid_list
from RNAcodeWebCore import concat_sequences_fasta

import SeqSelection
from SeqSelection import collect_sequences

RNAcodeWebCore.VERBOSE = True

MAX_TARGET_SEQS_BLASTN = 10000
E_VALUE_CUTOFF = 0.01

WORD_SIZE_DIC = {1: 14, 2: 11, 3: 9, 4: 8, 5: 7}

# Will be set in get_arguments()
JOB_ID = None
MIN_PAIR_DIST = None
MAX_PAIR_DIST = None
INPUT_SEQ_NUC = None

BLAST_WRAPPER_TEMPLATE = (
    "#!/bin/bash\n\n"
    "set -e\n\n"
    f'blastn -db "{{}}" -query {{}} -max_target_seqs {MAX_TARGET_SEQS_BLASTN} '
    f'-taxidlist {{}} -outfmt "{RNAcodeWebCore.OUTFMT}" -out {{}} -evalue {E_VALUE_CUTOFF} '
    "-num_threads {} -max_hsps 1 -task blastn -word_size {}\n"
)

SEQSEL_WRAPPER_TEMPLATE = (
    "#!/bin/bash\n\n"
    "set -e\n\n"
    "source $PYTHON_ENV/bin/activate\n\n"
    'python3 SeqSelection_pipeline.py {} {} {} "{}"\n'
)


def make_readme():
    """Build README for workdir."""
    readme_path = RNAcodeWebCore.README_PATH_TEMPLATE.format(JOB_ID)
    context = {}
    context["blast_version"] = get_version("blastn")

    context["name_aligner"] = "Clustal Omega"
    context["aligner_version"] = get_version("clustalo")

    context["rnacode_version"] = get_version("RNAcode")

    with open(readme_path, "w", encoding="UTF-8") as file_handle:
        file_handle.write(render_jinja("templates/pipeline/README_NCBI_db.md", context))


def init_work_dir():
    """Initialize work directory for analysis."""
    line_length = 60
    current_work_dir = RNAcodeWebCore.CURRENT_WORK_DIR_TEMPLATE.format(JOB_ID)
    job_status_file = RNAcodeWebCore.JOB_STATUS_FILE_TEMPLATE.format(JOB_ID)
    python_req_path = RNAcodeWebCore.PYTHON_REQ_PATH_TEMPLATE.format(JOB_ID)
    reference_species_file = RNAcodeWebCore.REFERENCE_SPECIES_FILE_TEMPLATE.format(
        JOB_ID
    )
    selected_species_file = RNAcodeWebCore.SELECTED_SPECIES_FILE_TEMPLATE.format(JOB_ID)
    seq_selection_module_path = (
        RNAcodeWebCore.SEQ_SELECTION_MODULE_PATH_TEMPLATE.format(JOB_ID)
    )
    seq_selection_script_path = (
        RNAcodeWebCore.SEQ_SELECTION_SCRIPT_PATH_TEMPLATE.format(JOB_ID)
    )
    input_file_path = RNAcodeWebCore.INPUT_FILE_PATH_TEMPLATE.format(JOB_ID)
    eps_path = RNAcodeWebCore.EPS_PATH_TEMPLATE.format(JOB_ID)

    if not os.path.isdir(current_work_dir):
        os.mkdir(current_work_dir)
        os.mkdir(eps_path)
    else:
        rmtree(current_work_dir)
        os.mkdir(current_work_dir)
        os.mkdir(eps_path)

    with open(input_file_path, "w", encoding="UTF-8") as file_handle:
        file_handle.write(">Target\n")
        seq = "\n".join(re.findall(f".{{1,{line_length}}}", INPUT_SEQ_NUC)) + "\n"
        file_handle.write(seq)

    job_status = {
        "blastn_1": ["R", 0],
        "seqSel_1": ["NS", 0],
        "alignment": ["NS", 0],
        "RNAcode": ["NS", 0],
        "fullJob": ["R", 0],
    }

    with open(job_status_file, "w", encoding="UTF-8") as file_handle:
        json.dump(job_status, file_handle, indent=4)

    with open(reference_species_file, "w", encoding="UTF-8") as file_handle:
        file_handle.write("None")

    with open(selected_species_file, "w", encoding="UTF-8") as file_handle:
        file_handle.write("None")

    make_readme()

    copyfile("./templates/pipeline/python_requirements.txt", python_req_path)
    copyfile("./SeqSelection.py", seq_selection_module_path)
    copyfile("./SeqSelection_NCBI_DB.py", seq_selection_script_path)


def get_arguments():
    """Get arguments from sys."""
    global JOB_ID, MIN_PAIR_DIST, MAX_PAIR_DIST, INPUT_SEQ_NUC

    JOB_ID = sys.argv[1]
    RNAcodeWebCore.JOB_ID = JOB_ID

    RNAcodeWebCore.STDOUT_FILE_PATH = RNAcodeWebCore.STDOUT_FILE_PATH.format(JOB_ID)
    RNAcodeWebCore.STDERR_FILE_PATH = RNAcodeWebCore.STDERR_FILE_PATH.format(JOB_ID)

    with open(RNAcodeWebCore.STDOUT_FILE_PATH, "w", encoding="UTF-8") as file_handle:
        file_handle.write(" ".join(sys.argv) + "\n")
    with open(RNAcodeWebCore.STDERR_FILE_PATH, "w", encoding="UTF-8") as file_handle:
        file_handle.write("")

    # Uncomment will print to stdout/ stderr
    # RNAcodeWebCore.STDOUT_FILE_PATH = None
    # RNAcodeWebCore.STDERR_FILE_PATH = None

    MIN_PAIR_DIST = float(sys.argv[2])
    MAX_PAIR_DIST = float(sys.argv[3])
    RNAcodeWebCore.GENOME_START = int(sys.argv[4])
    RNAcodeWebCore.DB_TYPE = sys.argv[5]
    INPUT_SEQ_NUC = sys.argv[6]

    if RNAcodeWebCore.DB_TYPE == "nt":
        RNAcodeWebCore.BLAST_DB = "nt"
    elif RNAcodeWebCore.DB_TYPE == "refseq":
        RNAcodeWebCore.BLAST_DB = "ref_euk_rep_genomes ref_prok_rep_genomes ref_viroids_rep_genomes ref_viruses_rep_genomes"

    RNAcodeWebCore.NCBI = SeqSelection.load_ncbi()


def blastn(iteration):
    """Call blast."""
    job_type = f"blastn_{iteration}"
    input_file_path = RNAcodeWebCore.INPUT_FILE_PATH_TEMPLATE.format(JOB_ID)
    blast_result_path = RNAcodeWebCore.BLAST_RESULT_PATH_TEMPLATE.format(
        JOB_ID, iteration
    )
    blast_wrapper_path = RNAcodeWebCore.BLAST_WRAPPER_PATH_TEMPLATE.format(
        JOB_ID, iteration
    )
    taxids_path = RNAcodeWebCore.TAXIDS_PATH_TEMPLATE.format(JOB_ID, iteration)

    error_file = RNAcodeWebCore.SLURM_ERROR_TEMPLATE.format(JOB_ID, job_type)
    idle_cores = RNAcodeWebCore.get_number_idles_cpus()
    cores = [v for k, v in RNAcodeWebCore.CORE_USE_DIC.items() if idle_cores in k][0]

    build_taxid_list(iteration)

    blast_wrapper = BLAST_WRAPPER_TEMPLATE.format(
        RNAcodeWebCore.BLAST_DB,
        input_file_path.split("/")[-1],
        taxids_path.split("/")[-1],
        blast_result_path.split("/")[-1],
        cores,
        WORD_SIZE_DIC[iteration],
    )

    with open(blast_wrapper_path, "w", encoding="UTF-8") as file_handle:
        file_handle.write(blast_wrapper)

    if not RNAcodeWebCore.slurm_batch(blast_wrapper_path, job_type, cores=cores):
        eprint("Blast exited with error!")
        with open(error_file, "r", encoding="UTF-8") as file_handle:
            eprint(file_handle.read())
        raise RNAcodeWebCore.PipelineError


def seq_selection(i):
    """2. Step of analysis. Select candidates from blast output."""
    job_type = f"seqSel_{i}"
    seqsel_wrapper_path = RNAcodeWebCore.SEQSEL_WRAPPER_PATH_TEMPLATE.format(JOB_ID, i)
    error_file = RNAcodeWebCore.SLURM_ERROR_TEMPLATE.format(JOB_ID, job_type)

    seqsel_wrapper = SEQSEL_WRAPPER_TEMPLATE.format(
        MIN_PAIR_DIST,
        MAX_PAIR_DIST,
        i,
        RNAcodeWebCore.BLAST_DB,
    )

    with open(seqsel_wrapper_path, "w", encoding="UTF-8") as file_handle:
        file_handle.write(seqsel_wrapper)

    if not RNAcodeWebCore.slurm_batch(seqsel_wrapper_path, job_type):
        eprint("Sequence selection exited with error!")
        with open(error_file, "r", encoding="UTF-8") as file_handle:
            eprint(file_handle.read())
        raise RNAcodeWebCore.PipelineError


def main():
    """Back end service for RNAcode web. Handles all computation etc."""
    get_arguments()

    vprint("Initializing work directory")
    init_work_dir()
    for iteration in range(1, 6):
        write_status_file(f"seqSel_{iteration}.{JOB_ID}", ["NS", "0"])
        write_status_file(f"blastn_{iteration}.{JOB_ID}", ["NS", "0"])
        vprint(f"{iteration}. Iteration of selection")
        vprint("Start blast")
        blastn(iteration)
        vprint("Select sequences")
        seq_selection(iteration)
        num_seqs = len(
            collect_sequences(
                iteration,
                current_work_dir=RNAcodeWebCore.CURRENT_WORK_DIR_TEMPLATE.format(
                    JOB_ID
                ),
            )
        )
        vprint(f"{num_seqs} sequences found.")
        if num_seqs >= SeqSelection.MAX_NUM_SEQ:
            break

    # checks if enough candidates had been found
    check_num_seq(iteration)
    concat_sequences_fasta(iteration)
    vprint("Aligne sequences")
    clustalo()
    vprint("Plot tree")
    plot_alignment()
    vprint("Make MAF")
    convert_to_maf()
    vprint("Run RNAcode")
    rnacode()
    write_status_file(f"fullJob.{JOB_ID}", ["CD", 0])
    notify_frontend(JOB_ID)
    vprint("Finished")


if __name__ == "__main__":
    start_pipeline(main)
