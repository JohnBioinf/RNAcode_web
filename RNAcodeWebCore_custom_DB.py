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
from RNAcodeWebCore import concat_sequences_fasta

from SeqSelection import DB_SIZE

RNAcodeWebCore.VERBOSE = True

E_VALUE_CUTOFF = 0.01

# Will be set in get_arguments()
JOB_ID = None
MIN_PAIR_DIST = None
MAX_PAIR_DIST = None
INPUT_SEQ_NUC = None
CUSTOM_DB_PATH = None

BLAST_WRAPPER_TEMPLATE = (
    "#!/bin/bash\n\n"
    "set -e\n\n"
    f"blastn -db {{}} -query {{}} -max_target_seqs {DB_SIZE} "
    f'-outfmt "{RNAcodeWebCore.OUTFMT}" -out {{}} -evalue {E_VALUE_CUTOFF} '
    "-num_threads {} -max_hsps 1 -task blastn -word_size 5\n"
)

SEQSEL_WRAPPER_TEMPLATE = (
    "#!/bin/bash\n\n"
    "set -e\n\n"
    "source $PYTHON_ENV/bin/activate\n\n"
    "python3 SeqSelection_pipeline.py {} {} {}\n"
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
        file_handle.write(render_jinja("templates/pipeline/README_custom_db.md", context))


def init_work_dir():
    """Initialize work directory for analysis."""
    line_length = 60
    current_work_dir = RNAcodeWebCore.CURRENT_WORK_DIR_TEMPLATE.format(JOB_ID)
    job_status_file = RNAcodeWebCore.JOB_STATUS_FILE_TEMPLATE.format(JOB_ID)
    python_req_path = RNAcodeWebCore.PYTHON_REQ_PATH_TEMPLATE.format(JOB_ID)
    seq_selection_module_path = RNAcodeWebCore.SEQ_SELECTION_MODULE_PATH_TEMPLATE.format(JOB_ID)
    seq_selection_script_path = RNAcodeWebCore.SEQ_SELECTION_SCRIPT_PATH_TEMPLATE.format(JOB_ID)
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

    make_readme()

    copyfile("./templates/pipeline/python_requirements.txt", python_req_path)
    copyfile("./SeqSelection.py", seq_selection_module_path)
    copyfile("./SeqSelection_custom_DB.py", seq_selection_script_path)


def get_arguments():
    """Get arguments from sys."""
    global JOB_ID, MIN_PAIR_DIST, MAX_PAIR_DIST, INPUT_SEQ_NUC, CUSTOM_DB_PATH

    JOB_ID = sys.argv[1]
    RNAcodeWebCore.JOB_ID = JOB_ID
    RNAcodeWebCore.STDOUT_FILE_PATH = RNAcodeWebCore.STDOUT_FILE_PATH.format(JOB_ID)
    RNAcodeWebCore.STDERR_FILE_PATH = RNAcodeWebCore.STDERR_FILE_PATH.format(JOB_ID)

    with open(RNAcodeWebCore.STDOUT_FILE_PATH, "w", encoding="UTF-8") as file_handle:
        file_handle.write(" ".join(sys.argv) + "\n")
    with open(RNAcodeWebCore.STDERR_FILE_PATH, "w", encoding="UTF-8") as file_handle:
        file_handle.write("")

    # Uncomment to print to stdout/ stderr
    # RNAcodeWebCore.STDOUT_FILE_PATH = None
    # RNAcodeWebCore.STDERR_FILE_PATH = None

    MIN_PAIR_DIST = float(sys.argv[2])
    MAX_PAIR_DIST = float(sys.argv[3])
    RNAcodeWebCore.GENOME_START = int(sys.argv[4])
    INPUT_SEQ_NUC = sys.argv[5]

    parent_job_id = "-".join(JOB_ID.split("-")[:-1])
    CUSTOM_DB_PATH = f"../{parent_job_id}/{RNAcodeWebCore.CUSTOM_DB_PATH}/{parent_job_id}"


def blastn():
    """Call blast."""
    job_type = "blastn_1"
    input_file_path = RNAcodeWebCore.INPUT_FILE_PATH_TEMPLATE.format(JOB_ID)
    blast_result_path = RNAcodeWebCore.BLAST_RESULT_PATH_TEMPLATE.format(JOB_ID, 1)
    blast_wrapper_path = RNAcodeWebCore.BLAST_WRAPPER_PATH_TEMPLATE.format(JOB_ID, 1)

    error_file = RNAcodeWebCore.SLURM_ERROR_TEMPLATE.format(JOB_ID, job_type)
    idle_cores = RNAcodeWebCore.get_number_idles_cpus()
    cores = [v for k, v in RNAcodeWebCore.CORE_USE_DIC.items() if idle_cores in k][0]

    blast_wrapper = BLAST_WRAPPER_TEMPLATE.format(
        CUSTOM_DB_PATH,
        input_file_path.split("/")[-1],
        blast_result_path.split("/")[-1],
        cores,
    )

    with open(blast_wrapper_path, "w", encoding="UTF-8") as file_handle:
        file_handle.write(blast_wrapper)

    if not RNAcodeWebCore.slurm_batch(blast_wrapper_path, job_type, cores=cores):
        eprint("Blast exited with error!")
        with open(error_file, "r", encoding="UTF-8") as file_handle:
            eprint(file_handle.read())
        raise RNAcodeWebCore.PipelineError


def seq_selection():
    """2. Step of analysis. Select candidates from blast output."""
    job_type = "seqSel_1"
    seqsel_wrapper_path = RNAcodeWebCore.SEQSEL_WRAPPER_PATH_TEMPLATE.format(JOB_ID, 1)
    error_file = RNAcodeWebCore.SLURM_ERROR_TEMPLATE.format(JOB_ID, job_type)

    with open(seqsel_wrapper_path, "w", encoding="UTF-8") as file_handle:
        file_handle.write(SEQSEL_WRAPPER_TEMPLATE.format(MIN_PAIR_DIST, MAX_PAIR_DIST, CUSTOM_DB_PATH))

    if not RNAcodeWebCore.slurm_batch(seqsel_wrapper_path, job_type):
        eprint("Region selection exited with error!")
        with open(error_file, "r", encoding="UTF-8") as file_handle:
            eprint(file_handle.read())
        raise RNAcodeWebCore.PipelineError


def main():
    """Back end service for RNAcode web. Handles all computation etc."""
    get_arguments()

    vprint("Initializing work directory")
    init_work_dir()
    vprint("Blast")
    blastn()
    vprint("Sequence selection")
    seq_selection()
    # checks if enough candidates had been found
    check_num_seq(1)
    concat_sequences_fasta(1)
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
