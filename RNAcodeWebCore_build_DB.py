"""Build custom blast DB for parent jobs."""
import sys
import os
import re
import json
from shutil import rmtree, copyfile
from time import sleep

import RNAcodeWebCore
from RNAcodeWebCore import start_pipeline
from RNAcodeWebCore import eprint
from RNAcodeWebCore import vprint
from RNAcodeWebCore import write_status_file
from RNAcodeWebCore import notify_frontend
from RNAcodeWebCore import build_taxid_list
from RNAcodeWebCore import get_version
from RNAcodeWebCore import render_jinja

import SeqSelection
from SeqSelection import DB_SIZE
from SeqSelection import collect_sequences
from SeqSelection_build_DB import TAXIDMAPFILE_PATH

RNAcodeWebCore.VERBOSE = True

MAX_TARGET_SEQS_BLASTN = DB_SIZE * 1000
E_VALUE_CUTOFF = 0.01

# Word size for blast in each iteration.
WORD_SIZE_DIC = {1: 18, 2: 14, 3: 10, 4: 8, 5: 7}

# will be set in get_arguments()
JOB_ID = None
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
    'python3 SeqSelection_pipeline.py {} "{}"\n'
)

BUILD_DB_WRAPPER_TEMPLATE = (
    "#!/bin/bash\n\n"
    "set -e\n\n"
    "makeblastdb -in {0} -parse_seqids -dbtype nucl -taxid_map {3} -title {1} -out {2}/{1}\n"
)


def get_arguments():
    """Get arguments from sys."""
    global JOB_ID, INPUT_SEQ_NUC
    JOB_ID = sys.argv[1]
    RNAcodeWebCore.JOB_ID = JOB_ID

    RNAcodeWebCore.STDOUT_FILE_PATH = RNAcodeWebCore.STDOUT_FILE_PATH.format(JOB_ID)
    RNAcodeWebCore.STDERR_FILE_PATH = RNAcodeWebCore.STDERR_FILE_PATH.format(JOB_ID)

    with open(RNAcodeWebCore.STDOUT_FILE_PATH, "w", encoding="UTF-8") as file_handle:
        file_handle.write(" ".join(sys.argv) + "\n")
    with open(RNAcodeWebCore.STDERR_FILE_PATH, "w", encoding="UTF-8") as file_handle:
        file_handle.write("")
    # Uncomment to print to stdout/stderr
    # RNAcodeWebCore.STDOUT_FILE_PATH = None
    # RNAcodeWebCore.STDERR_FILE_PATH = None

    INPUT_SEQ_NUC = sys.argv[2]
    RNAcodeWebCore.DB_TYPE = sys.argv[3]

    if RNAcodeWebCore.DB_TYPE == "nt":
        RNAcodeWebCore.BLAST_DB = "nt"
    elif RNAcodeWebCore.DB_TYPE == "refseq":
        RNAcodeWebCore.BLAST_DB = "ref_euk_rep_genomes ref_prok_rep_genomes ref_viroids_rep_genomes ref_viruses_rep_genomes"

    RNAcodeWebCore.NCBI = SeqSelection.load_ncbi()


def make_readme():
    """Build README for workdir."""
    readme_path = RNAcodeWebCore.README_PATH_TEMPLATE.format(JOB_ID)
    context = {}
    context["blast_version"] = get_version("blastn")

    context["name_aligner"] = "Clustal Omega"

    with open(readme_path, "w", encoding="UTF-8") as file_handle:
        file_handle.write(render_jinja("templates/pipeline/README_build_db.md", context))


def init_work_dir():
    """Initialize work directory for analysis."""
    line_length = 60
    current_work_dir = RNAcodeWebCore.CURRENT_WORK_DIR_TEMPLATE.format(JOB_ID)
    job_status_file = RNAcodeWebCore.JOB_STATUS_FILE_TEMPLATE.format(JOB_ID)
    input_file_path = RNAcodeWebCore.INPUT_FILE_PATH_TEMPLATE.format(JOB_ID)
    seq_selection_script_path = RNAcodeWebCore.SEQ_SELECTION_SCRIPT_PATH_TEMPLATE.format(JOB_ID)
    seq_selection_module_path = RNAcodeWebCore.SEQ_SELECTION_MODULE_PATH_TEMPLATE.format(JOB_ID)
    python_req_path = RNAcodeWebCore.PYTHON_REQ_PATH_TEMPLATE.format(JOB_ID)
    selected_species_file = current_work_dir + "/" + SeqSelection.SELECTED_SPECIES_FILE_PATH
    taxidmapfile_path = current_work_dir + "/" + TAXIDMAPFILE_PATH
    custom_db_path = current_work_dir + "/" + RNAcodeWebCore.CUSTOM_DB_PATH

    if not os.path.isdir(current_work_dir):
        os.mkdir(current_work_dir)
        os.mkdir(custom_db_path)
    else:
        rmtree(current_work_dir)
        os.mkdir(current_work_dir)
        os.mkdir(custom_db_path)

    with open(input_file_path, "w", encoding="UTF-8") as file_handle:
        file_handle.write(">Target\n")
        seq = "\n".join(re.findall(f".{{1,{line_length}}}", INPUT_SEQ_NUC)) + "\n"
        file_handle.write(seq)

    job_status = {
        "blastn_1": ["R", 0],
        "seqSel_1": ["NS", 0],
        "buildDB": ["NS", 0],
        "fullJob": ["R", 0],
    }

    with open(job_status_file, "w", encoding="UTF-8") as file_handle:
        json.dump(job_status, file_handle, indent=4)

    with open(selected_species_file, "w", encoding="UTF-8") as file_handle:
        file_handle.write("None")

    with open(taxidmapfile_path, "w", encoding="UTF-8") as file_handle:
        file_handle.write("")

    copyfile("./templates/pipeline/python_requirements.txt", python_req_path)
    copyfile("./SeqSelection.py", seq_selection_module_path)
    copyfile("./SeqSelection_build_DB.py", seq_selection_script_path)
    make_readme()


def blastn(iteration):
    """Call blast."""
    job_type = f"blastn_{iteration}"
    input_file_path = RNAcodeWebCore.INPUT_FILE_PATH_TEMPLATE.format(JOB_ID)
    blast_result_path = RNAcodeWebCore.BLAST_RESULT_PATH_TEMPLATE.format(JOB_ID, iteration)
    blast_wrapper_path = RNAcodeWebCore.BLAST_WRAPPER_PATH_TEMPLATE.format(JOB_ID, iteration)
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


def build_db():
    """Build blast db."""
    job_type = "buildDB"
    blast_db_fasta_path = RNAcodeWebCore.BLAST_DB_FASTA_PATH_TEMPLATE.format(JOB_ID)
    custom_db_path = RNAcodeWebCore.CUSTOM_DB_PATH
    build_db_wrapper_path = RNAcodeWebCore.BUILD_DB_WRAPPER_PATH_TEMPLATE.format(JOB_ID)

    error_file = RNAcodeWebCore.SLURM_ERROR_TEMPLATE.format(JOB_ID, job_type)

    build_db_wrapper = BUILD_DB_WRAPPER_TEMPLATE.format(
        blast_db_fasta_path.split("/")[-1],
        JOB_ID,
        custom_db_path,
        TAXIDMAPFILE_PATH,
    )

    with open(build_db_wrapper_path, "w", encoding="UTF-8") as file_handle:
        file_handle.write(build_db_wrapper)

    if not RNAcodeWebCore.slurm_batch(build_db_wrapper_path, job_type):
        eprint("makeblastdb exited with error!")
        with open(error_file, "r", encoding="UTF-8") as file_handle:
            eprint(file_handle.read())
        raise RNAcodeWebCore.PipelineError


def seq_selection(iteration):
    """2. Step of analysis. Select candidates from blast output."""
    job_type = f"seqSel_{iteration}"
    seqsel_wrapper_path = RNAcodeWebCore.SEQSEL_WRAPPER_PATH_TEMPLATE.format(JOB_ID, iteration)
    error_file = RNAcodeWebCore.SLURM_ERROR_TEMPLATE.format(JOB_ID, job_type)

    with open(seqsel_wrapper_path, "w", encoding="UTF-8") as file_handle:
        file_handle.write(SEQSEL_WRAPPER_TEMPLATE.format(iteration, RNAcodeWebCore.BLAST_DB))

    if not RNAcodeWebCore.slurm_batch(seqsel_wrapper_path, job_type):
        eprint("Region selection exited with error!")
        with open(error_file, "r", encoding="UTF-8") as file_handle:
            eprint(file_handle.read())
        raise RNAcodeWebCore.PipelineError


def count_last_blast_results(iteration):
    """Count the number of hits found by blastn in the last search."""
    blast_result_path = RNAcodeWebCore.BLAST_RESULT_PATH_TEMPLATE.format(JOB_ID, iteration)
    with open(blast_result_path, "r", encoding="UTF-8") as file_handle:
        return len(list(file_handle.readlines()))


def concat_sequences_fasta(iteration):
    """Concatinate all selected sequences in one."""
    blast_db_fasta_path = RNAcodeWebCore.BLAST_DB_FASTA_PATH_TEMPLATE.format(JOB_ID)
    with open(blast_db_fasta_path, "w", encoding="UTF-8") as db_handle:
        for i in range(1, iteration + 1):
            selected_sequences_path = f"{RNAcodeWebCore.SELECTED_SEQUENCES_PATH_TEMPLATE.format(JOB_ID, i)}"
            if not os.path.isfile(selected_sequences_path):
                continue
            with open(selected_sequences_path, "r", encoding="UTF-8") as seq_handle:
                db_handle.write(
                    re.sub(r":c*[0-9]+-[0-9]+", "", seq_handle.read())
                )


def main():
    """Build custom blast DB for parent jobs."""
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
        num_seqs = len(collect_sequences(iteration + 1, current_work_dir=RNAcodeWebCore.CURRENT_WORK_DIR_TEMPLATE.format(JOB_ID)))
        vprint(f"{num_seqs} sequences found.")
        if num_seqs >= DB_SIZE:
            break

    if num_seqs < RNAcodeWebCore.MIN_NUM_SEQ:
        write_status_file(f"fullJob.{JOB_ID}", ["F", num_seqs])
        raise RNAcodeWebCore.PipelineFailed

    concat_sequences_fasta(iteration)
    build_db()
    # Wait sometime. For reason beyond my comprehension some files seem to be
    # not written to disk when the process ended. Might cause error in child
    # processes.
    sleep(4)
    write_status_file(f"fullJob.{JOB_ID}", ["CD", num_seqs])
    notify_frontend(JOB_ID)


if __name__ == "__main__":
    start_pipeline(main)
