#!/usr/bin/python3
"""Module for the class job submission.

This module gives all functionality so that the frontend server can handle the
job. This includes all communication with the backend. Further each instance can
be locked so that only one instance can write. To use this you can use the "with
... as ...:" syntax.
"""

import os
import json
import datetime
from time import sleep
import re
import subprocess
import getpass

from filelock import FileLock
from coolname import generate_slug
from Bio import SeqIO
from Bio.Seq import Seq
from validate_email import validate_email

from RNAcodeWebCore import JOB_STATUS_FILE_TEMPLATE
from RNAcodeWebCore import RNACODE_RESULT_PATH_TEMPLATE
from RNAcodeWebCore import EPS_PATH_TEMPLATE
from RNAcodeWebCore import ALIGN_PATH_TEMPLATE
from RNAcodeWebCore import ALIGN_PLOT_PATH_TEMPLATE
from RNAcodeWebCore import WORK_DIR
from RNAcodeWebCore import CURRENT_WORK_DIR_TEMPLATE


def vprint(msg, verbose_level=0):
    """Print to verbose."""
    check_verbose_level(verbose_level)
    msg = datetime.datetime.today().strftime("[%d/%b/%Y %H:%M:%S] - - ") + str(msg)
    if verbose_level <= VERBOSE_LEVEL:
        with open(f"{LOG_DIR}/{os.getpid()}.log", "a", encoding="UTF-8") as f_handle:
            print(msg, file=f_handle, flush=True)


# Parameters
with open("./parameters_frontend_local.json", "r", encoding="UTF-8") as file_handle:
    server_parameters_frontend = json.load(file_handle)

with open("./parameters_backend_local.json", "r", encoding="UTF-8") as file_handle:
    server_parameters_backend = json.load(file_handle)

RNACODE_WEB_REPO_BACKEND = server_parameters_backend["RNAcode_web_repo_backend"]
R_BIN = server_parameters_frontend["R_bin"]

# setup directory's
if getpass.getuser() == server_parameters_frontend["user_dev"]:
    LOG_DIR = "./logs"
    RESULT_DIR = "static/results"
    JOB_SUBMISSION_PATH_TEMPLATE = "./job_submissions/{}.json"
elif getpass.getuser() == server_parameters_frontend["user_production"]:
    PRODUCTION_DIR = server_parameters_frontend["production_dir"]
    LOG_DIR = PRODUCTION_DIR + "/logs"
    RESULT_DIR = PRODUCTION_DIR + "/results"
    JOB_SUBMISSION_DIR = PRODUCTION_DIR + "/job_submissions"
    JOB_SUBMISSION_PATH_TEMPLATE = JOB_SUBMISSION_DIR + "/{}.json"
else:
    LOG_DIR = "./logs"
    vprint(f"Unknown user: {os.getlogin()}")
    raise ValueError

# Result files
# ! Dependencies in ssh-scp-api/ bash script
ALIG_FILE_PATH_TEMPLATE = f"{RESULT_DIR}/{{}}_aln.html"
RNACODE_RESULT_FRONTEND_PATH_TEMPLATE = f"{RESULT_DIR}/{{}}_rnacode_result.tsv"
FULL_RESULTS_PATH_TEMPLATE = f"{RESULT_DIR}/{{}}_full_results.tgz"
HSS_TABLE_R_PATH_TEMPLATE = f"{RESULT_DIR}/{{}}_HSS_segment_R.csv"
HSS_PLOT_FRONTEND_PATH_TEMPLATE = f"{RESULT_DIR}/{{}}_HSS_segment.svg"
ALIG_FASTA_PATH_TEMPLATE = f"{RESULT_DIR}/{{}}_aln.fasta"
ALIG_INFO_CHILD_PATH_TEMPLATE = f"{RESULT_DIR}/{{}}_alignment_info_children.tsv"
HSS_TABLE_PATH_TEMPLATE = f"{RESULT_DIR}/{{}}_HSS_table.tsv"
PROTEIN_LIST_PATH_TEMPLATE = f"{RESULT_DIR}/{{}}_Protein_list.csv"

# Default parameters for the jobs
MIN_PAIR_DIST_DEFAULT = "10.0"
MAX_PAIR_DIST_DEFAULT = "60.0"

DATE_FORMATION = "%B %d, %Y"
# Level of how much stuff should be logged. Either 0, 1, 2, 3.
VERBOSE_LEVEL = 3

IUPAC_ALPHABET = set("ACGTRYSWKMBDHVN")
# Window size into which a child process will be cut from its parent job.
CHILD_LENGTH = 1000

# Every job with an input sequence greater than MAX_LEN_NO_SPLIT will be be
# split up into multiple child jobs.
MAX_LEN_NO_SPLIT = CHILD_LENGTH * 2

# Length below a job will use the standard NCBI DB. Every length above will
# build its custom database.
MIN_INPUT_LEN = 30
LENGTH_NO_CUSTOM = CHILD_LENGTH * 4
MAX_INPUT_LEN = 200 * CHILD_LENGTH
MAX_JOB_NAME_LENGTH = 80

# Bulk Submission
MAX_JOBS_BULK = 10
# keys of parameters in json for bulk submission
JSON_KEYS = set(
    [
        "input_seq_nuc",
        "job_id",
        "db_type",
        "email",
        "min_pair_dist",
        "max_pair_dist",
    ]
)


def blank_job_submission():
    """Return the attributes of an empty default job submission."""
    while True:
        job_id = generate_slug()
        path = JOB_SUBMISSION_PATH_TEMPLATE.format(job_id)
        if not os.path.isfile(path):
            break

    return {
        "date": datetime.datetime.today().strftime(DATE_FORMATION),
        "job_id": job_id,
        "invalid_job_id": "",
        "path": path,
        "lock_path": f"{path}.lock",
        "job_hierarchy": "orphan",
        "input_seq_nuc": "",
        "genome_start": None,
        "genome_stop": None,
        "min_pair_dist": MIN_PAIR_DIST_DEFAULT,
        "max_pair_dist": MAX_PAIR_DIST_DEFAULT,
        "db_type": "refseq",
        "custom_db": False,
        "email": "Nobody",
        "submitted": False,
        "custom_db_submitted": False,
        "email_notification": True,
        "submission_email_send": False,
        "finished": False,
        "got_full_results": False,
        "status_dic": {
            "blastn_1": ["NS", 0],
            "seqSel_1": ["NS", 0],
            "alignment": ["NS", 0],
            "RNAcode": ["NS", 0],
            "fullJob": ["R", 0],
            "jobStatus": ["R", 0],
        },
        "error_out": "",
        "input_check_dic": {
            "input_seq_nuc": "",
            "input_fasta": "",
            "email": "",
            "job_id": "",
            "fasta": "",
            "min_pair_dist": "",
            "max_pair_dist": "",
            "no_error": True,
        },
        "hss_table": [],
        "protein_list": [],
        "alig_info_table": [],
        "alignment_print": "",
        "last_plot_par": [],
    }


def check_bulk_json(json_file_path):
    """Check if the given json file is valid.

    Only checks the basic validity of the json fields fields etc. Not if the
    submission do have correct values.
    """
    bulk_check_dic = {
        "error": False,
        "no_file": "",
        "num_job": "",
        "json_formatting": "",
        "dic": [],
        "keys": [],
        "values": [],
        "job_ids": "",
    }
    bulk_check_dic["error"] = False

    with open(json_file_path, "r", encoding="UTF-8") as f_handle:
        try:
            parameter_list = json.load(f_handle)
        except UnicodeDecodeError:
            bulk_check_dic["json_formatting"] = "File must be UTF-8 decoded."
            bulk_check_dic["error"] = True
        except json.decoder.JSONDecodeError as err:
            bulk_check_dic["json_formatting"] = str(err)
            bulk_check_dic["error"] = True
        else:
            num_jobs = len(parameter_list)
            if num_jobs > MAX_JOBS_BULK + 1:
                bulk_check_dic["num_jobs"] = (
                    f"{MAX_JOBS_BULK}, the "
                    "array of the json file "
                    f"contained {num_jobs} fields."
                )
                bulk_check_dic["error"] = True
            else:
                job_ids = []
                bulk_check_dic["dic"] = []
                bulk_check_dic["keys"] = []
                bulk_check_dic["values"] = []
                for i, parameters in enumerate(parameter_list):
                    if not isinstance(parameters, dict):
                        err = f"Field {i+1} is not an associative arrays."
                        bulk_check_dic["dic"].append(err)
                        bulk_check_dic["error"] = True
                        continue
                    if set(parameters.keys()) != JSON_KEYS:
                        err = f"Field {i+1} has not the correct keys."
                        bulk_check_dic["keys"].append(err)
                        bulk_check_dic["error"] = True
                        continue
                    if not all((isinstance(v, str) for k, v in parameters.items())):
                        err = f"Field {i+1} has none string values."
                        bulk_check_dic["values"].append(err)
                        bulk_check_dic["error"] = True
                        continue
                    job_ids.append(parameters["job_id"])

            if (
                bulk_check_dic["dic"] == []
                and bulk_check_dic["keys"] == []
                and bulk_check_dic["values"] == []
            ):
                if len(job_ids) != len(set(job_ids)):
                    bulk_check_dic["job_ids"] = "The job ids must be unique"
                    bulk_check_dic["error"] = True
    return bulk_check_dic


def read_bulk_json(f_handle):
    """Read a json file for bulk submission.

    Returns a list of job_submissions and a dictionary that contains checks for
    each submission.
    """
    vprint("Read bulk json", verbose_level=1)
    parameter_list = json.load(f_handle)
    job_list = []

    for parameters in parameter_list:
        if parameters["email"] != "":
            parameters["email_notification"] = True

        with JobSubmission(form=parameters) as job_submission:
            job_list.append(job_submission)

    return job_list


def _read_from_json(job_submission_path):
    for _i in range(0, 3):
        with open(job_submission_path, "r", encoding="UTF-8") as f_handle:
            try:
                attributes = json.load(f_handle)
            except json.decoder.JSONDecodeError:
                sleep(2)
                continue
            else:
                break
    else:
        vprint(f"Persistent JSONDecodeError with file {job_submission_path}")
        raise JobNotFoundException
    return attributes


def generate_bulk_submission_example():
    """Return a json str for an example bulk submission."""
    job_1 = {}
    job_1["input_seq_nuc"] = "ACGT"
    job_1["job_id"] = blank_job_submission()["job_id"]
    job_1["email"] = "john_doe@nasa.gov"
    job_1["min_pair_dist"] = MIN_PAIR_DIST_DEFAULT
    job_1["max_pair_dist"] = MAX_PAIR_DIST_DEFAULT
    job_1["blast_db"] = "refseq"

    job_2 = {}
    job_2["input_seq_nuc"] = "TGCA"
    job_2["job_id"] = blank_job_submission()["job_id"]
    job_2["email"] = ""
    job_2["min_pair_dist"] = MIN_PAIR_DIST_DEFAULT
    job_2["max_pair_dist"] = MAX_PAIR_DIST_DEFAULT
    job_1["blast_db"] = "nt"

    return json.dumps([job_1, job_2], indent=4)


def check_verbose_level(verbose_level):
    """Check if verbose levl is in correct form."""
    if not isinstance(verbose_level, int):
        raise ValueError("verbose level must be an integer")
    if 0 > verbose_level > 3:
        raise ValueError("verbose level must be between 0 and 3")


def _ssh_call(call_str):
    for i in range(1, 4):
        try:
            completed_process = subprocess.run(
                call_str.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
            )
            break
        except subprocess.CalledProcessError as exc:
            # Some weird error catch.
            if "kex_exchange_identification" in str(exc.stderr.decode('UTF8')) and i < 3:
                sleep(2)
                continue
            vprint("error ssh call")
            vprint(f"Command\n{call_str}")
            vprint(f"stdout:\n{exc.stdout.decode('UTF8')}")
            vprint(f"stderr:\n{exc.stderr.decode('UTF8')}")
            raise
    return completed_process.stdout.decode('UTF8')


def merge_hss_table(hss_table):
    """Merge segments."""
    vprint("Merge hss table", verbose_level=2)

    def do_overlap(seg_a, seg_b):
        return (
            seg_a[1] == seg_b[1]
            # frame
            and seg_a[2] == seg_b[2]
            # overlap
            and int(seg_b[7]) <= int(seg_a[8])
            and int(seg_a[7]) <= int(seg_b[8])
        )

    def merge(seg_a, seg_b):
        new_id = seg_a[0] + "-" + seg_b[0]
        plot_list = seg_a[11] + seg_b[11]
        new_start = str(min((int(seg_a[7]), int(seg_b[7]))))
        new_end = str(max((int(seg_a[8]), int(seg_b[8]))))
        if float(seg_a[10]) > float(seg_b[10]):
            seg_a = seg_b
        seg_a[0] = new_id
        seg_a[7] = new_start
        seg_a[8] = new_end
        seg_a[11] = plot_list
        return seg_a

    segments = []
    # Make list of list so I can merge in place
    for line in hss_table:
        for segment in segments:
            if do_overlap(segment[0], line):
                segment[0] = merge(segment[0], line)
                break
        else:
            segments.append([line])
    # Unpack list of lists
    return [seg[0] for seg in segments]


def _build_r_table(segments):
    r_table_list = [
        [
            "start",
            "end",
            "height_arrow",
            "height_text",
            "color",
            "arrow_end",
            "label",
        ]
    ]
    for segment in segments:
        start = int(segment[7]) - 1
        end = int(segment[8])
        frame = int(segment[2])
        hss_id = segment[0]
        if segment[1] == "+":
            strand = 1
        else:
            strand = -1

        height_arrow = (0.2 * strand * frame) - (0.06 * strand)
        height_text = height_arrow + 0.1 * strand
        arrow_end = "last" if strand == 1 else "first"
        r_table_list.append(
            [
                start,
                end,
                height_text,
                height_arrow,
                strand * frame,
                arrow_end,
                hss_id,
            ]
        )
    return r_table_list


PUBLIC_ATTR = blank_job_submission().keys()


class NoSaveException(Exception):
    """Error is used to not save a JobSubmission.

    Throw in method if you wish that the JobSubmission should not be saved.
    Will be cought in __exit__ and prevents the method self.save().
    Important for mehtod self.delete().
    """


class JobNotFoundException(Exception):
    """Exception if a job can not be found."""


class AlreadySubmittedException(Exception):
    """Exception raised if a job that hast already been submitted is tried to submit again."""


class JobSubmission:
    """Represent a job submitted by a user."""

    def __init__(
        self,
        job_id=None,
        form=None,
        attributes=None,
        fasta_file=None,
    ):
        """Init empty class."""
        if job_id:
            vprint(f"Load submission {job_id}", verbose_level=2)
            job_submission_path = JOB_SUBMISSION_PATH_TEMPLATE.format(job_id)
            if len(re.findall("[^a-zA-Z0-9-_]+", job_id)) != 0:
                vprint("Invalid chars", verbose_level=2)
                raise JobNotFoundException
            if not os.path.isfile(job_submission_path):
                vprint("Job does not exist", verbose_level=2)
                raise JobNotFoundException
            attributes = _read_from_json(job_submission_path)
            is_new = False
        elif attributes:
            vprint("Set from attributes", verbose_level=2)
            is_new = True
        else:
            vprint("Make blank job", verbose_level=2)
            attributes = blank_job_submission()
            is_new = True

        # Will be set on enter
        self.lock = None

        for key, value in attributes.items():
            setattr(self, key, value)

        if form:
            vprint("Set attributes from form")
            self.set_attributes(form, fasta_file=fasta_file)

        if is_new:
            self.read_only = False
            self.save()
        self.read_only = True

    def __enter__(self):
        """Make lock on job submission file."""
        self.lock = FileLock(self.lock_path, timeout=600)
        vprint(
            f"Try to enter {self.job_id}.",
            verbose_level=3,
        )
        self.lock.acquire()
        vprint("Lock acquired", verbose_level=3)
        attributes = _read_from_json(self.path)
        for key, value in attributes.items():
            setattr(self, key, value)
        vprint(f"Enter {self.job_id}", verbose_level=3)
        self.read_only = False
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        """Close lock on job submission file."""
        if exc_value is None:
            self.save()
        self.lock.release()
        vprint(f"Exit {self.job_id}", verbose_level=3)
        if exc_type is NoSaveException:
            return True
        return exc_type is None

    def get_attributes(self):
        """Return all public attributes."""
        return {attr: getattr(self, attr) for attr in PUBLIC_ATTR}

    def save(self):
        """Safe job submission to file."""
        vprint(f"Save submission {self.job_id}", verbose_level=2)
        attributes = self.get_attributes()
        with open(self.path, "w", encoding="UTF-8") as f_handle:
            json.dump(attributes, f_handle, indent=4)

    def delete(self, backend=True, recursive=True):
        """Remove job submission file."""
        vprint(f"Delete job {self.job_id}", verbose_level=2)
        vprint(self.job_hierarchy, verbose_level=2)

        if self.job_hierarchy == "parent" and recursive and self.submitted:
            vprint("Delete all children", verbose_level=3)
            for child_job_id, _child_status in self.status_dic.items():
                if child_job_id in ["jobStatus", "customDB"]:
                    continue
                try:
                    with JobSubmission(child_job_id) as child_job:
                        child_job.delete()
                # The clean up process might delete child before parent
                except JobNotFoundException:
                    continue

        os.remove(self.path)
        if backend and self.submitted and self.job_hierarchy != "parent":
            current_work_dir = CURRENT_WORK_DIR_TEMPLATE.format(self.job_id)
            out_file = RNACODE_WEB_REPO_BACKEND + f"logs/{self.job_id}.out"
            err_file = RNACODE_WEB_REPO_BACKEND + f"logs/{self.job_id}.err"

            call_str = (
                f"ssh-scp-api/delete_job.sh {current_work_dir} {out_file} {err_file}"
            )
            _ssh_call(call_str)
        raise NoSaveException

    def set_attributes(self, form, fasta_file=None):
        """Set attributes from html request."""
        self._set_input_seq_nuc(form, fasta_file)
        self._set_db_type(form)
        self._set_dist(form, "max_pair_dist")
        self._set_dist(form, "min_pair_dist")
        self._check_dist()
        self._set_email(form)
        self._set_job_id(form)
        self.input_check_dic["no_error"] = all(
            (v == "" for k, v in self.input_check_dic.items() if k != "no_error")
        )
        if self.input_check_dic["no_error"]:
            vprint("Form was vaild", verbose_level=2)
        else:
            vprint("Form was invalid", verbose_level=2)
            for attribute, error in self.input_check_dic.items():
                if attribute == "no_error":
                    continue
                if error != "":
                    vprint(
                        f"Attribute {attribute}: {getattr(self,attribute)} has error: {error}.",
                        verbose_level=2,
                    )

    def _set_input_seq_nuc(self, form, fasta_file):
        if fasta_file:
            input_seq_nuc = list(SeqIO.parse(fasta_file, "fasta"))
            if len(input_seq_nuc) == 0:
                self.input_check_dic["input_fasta"] = "Invalid fasta"
                self.input_seq_nuc = ""
                return
            if len(input_seq_nuc) != 1:
                self.input_check_dic["input_fasta"] = "Multi fasta"
                self.input_seq_nuc = ""
                return
            input_seq_nuc = (
                str(input_seq_nuc[0].seq).upper().replace("\n", "").replace("\r", "")
            )
        else:
            input_seq_nuc = (
                form["input_seq_nuc"].upper().replace("\n", "").replace("\r", "")
            )

        self.input_seq_nuc = input_seq_nuc

        if input_seq_nuc == "":
            self.input_check_dic["input_seq_nuc"] = "Empty"

        if len(set(input_seq_nuc).difference(IUPAC_ALPHABET)) != 0:
            bad_chars = ", ".join(set(input_seq_nuc).difference(IUPAC_ALPHABET))
            self.input_check_dic["input_seq_nuc"] = f"Bad chars:{bad_chars}"

        if len(input_seq_nuc) > MAX_INPUT_LEN:
            self.input_check_dic["input_seq_nuc"] = "To long"

        if len(input_seq_nuc) < MIN_INPUT_LEN:
            self.input_check_dic["input_seq_nuc"] = "To short"

        self.genome_start = 0
        self.genome_stop = len(input_seq_nuc)

    def _set_db_type(self, form):
        if form["db_type"] not in ["nt", "refseq"]:
            self.input_check_dic["db_type"] = "Wrong type"
        else:
            self.input_check_dic["db_type"] = ""
        self.db_type = form["db_type"]

    def _set_dist(self, form, name):
        distance = form[name].replace("%", "")
        try:
            distance = float(distance)
            setattr(self, name, distance)
        except ValueError:
            setattr(self, name, distance)
            self.input_check_dic[name] = "NaN"
        else:
            if distance < 0 or distance > 100:
                self.input_check_dic[name] = "Range"

    def _check_dist(self):
        if (
            self.input_check_dic["min_pair_dist"] != ""
            or self.input_check_dic["max_pair_dist"] != ""
        ):
            return

        if self.min_pair_dist >= self.max_pair_dist:
            self.input_check_dic["min_pair_dist"] = "Range"
            self.input_check_dic["max_pair_dist"] = "Range"

    def _set_email(self, form):
        if "email_notification" in form:
            self.email = form["email"]
            if not validate_email(self.email):
                self.input_check_dic["email"] = "Form"
            self.email_notification = True
        else:
            self.email_notification = False
            self.email = "Nobody"

    def _set_job_id(self, form):
        self.invalid_job_id = form["job_id"]
        job_submission_path = JOB_SUBMISSION_PATH_TEMPLATE.format(self.invalid_job_id)

        if len(self.invalid_job_id) > MAX_JOB_NAME_LENGTH:
            self.input_check_dic["job_id"] = "Length"
            return
        if len(re.findall("[^a-zA-Z0-9-_]+", self.invalid_job_id)) != 0:
            self.input_check_dic["job_id"] = "Char"
            return
        if os.path.isfile(job_submission_path):
            self.input_check_dic["job_id"] = "Already exists"
            return

        self.path = job_submission_path
        self.lock_path = self.path + ".lock"
        self.job_id = self.invalid_job_id

        if len(self.input_seq_nuc) > MAX_LEN_NO_SPLIT:
            self.job_hierarchy = "parent"
            self.status_dic = {}
            self.status_dic["jobStatus"] = ["R", 0]
            for iterator, start in enumerate(
                range(0, len(self.input_seq_nuc), int(CHILD_LENGTH / 2))
            ):
                self.status_dic[f"{self.job_id}-child_{iterator+1}"] = ["NS", 0]
                if start + CHILD_LENGTH >= len(self.input_seq_nuc):
                    break

        if len(self.input_seq_nuc) > LENGTH_NO_CUSTOM:
            self.custom_db = True
            self.status_dic["customDB"] = {
                "blastn_1": ["NS", 0],
                "seqSel_1": ["NS", 0],
                "buildDB": ["NS", 0],
                "fullJob": ["NS", 0],
            }
        else:
            self.custom_db = False

    def submit(self):
        """Submit job to backend.

        This method needs that the job submission is locked, hence should be
        used with 'with ... as ...' syntax.
        """
        vprint(f"Submit job {self.job_id}", verbose_level=2)
        if not self.input_check_dic["no_error"]:
            raise AttributeError(f"the input for job {self.job_id} is not valid.")
        if self.read_only:
            raise AttributeError("job submission has no write privilege")
        # First see if custom db must be build
        if self.custom_db and not self.custom_db_submitted:
            vprint("Submit build of custom db to backend.", verbose_level=2)
            call_str = f"ssh-scp-api/submit_build_db.sh {self.job_id} {self.input_seq_nuc} {self.db_type}"
            vprint(call_str, verbose_level=3)
            out = _ssh_call(call_str)
            vprint(out, verbose_level=3)
            self.custom_db_submitted = True
            self.status_dic["customDB"] = {
                "blastn_1": ["R", 0],
                "seqSel_1": ["NS", 0],
                "buildDB": ["NS", 0],
                "fullJob": ["NS", 0],
            }
            return

        if self.custom_db and self.status_dic["customDB"]["fullJob"][0] != "CD":
            vprint("Check if custom DB is finished.", verbose_level=3)
            self.check_status()
            if self.status_dic["customDB"]["fullJob"][0] != "CD":
                vprint("Building custom db is still running or failed", verbose_level=3)
                return

        # Now submit real job
        if self.submitted:
            vprint("The job was already submitted.", verbose_level=2)
            raise AlreadySubmittedException

        if self.job_hierarchy == "parent":
            vprint(f"Job {self.job_id} is long will be split up", verbose_level=2)
            len_parent_seq = len(self.input_seq_nuc)
            for iterator, start in enumerate(
                range(0, len_parent_seq, int(CHILD_LENGTH / 2))
            ):
                self._spawn_child(iterator, start)
                if start + CHILD_LENGTH >= len_parent_seq:
                    break
            vprint("All children submitted", verbose_level=3)
        else:
            if self.custom_db:
                call_str = (
                    "ssh-scp-api/submit_job_custom_db.sh "
                    f"{self.job_id} "
                    f"{self.min_pair_dist} "
                    f"{self.max_pair_dist} "
                    f"{self.genome_start} "
                    f"{self.input_seq_nuc}"
                )
            else:
                call_str = (
                    "ssh-scp-api/submit_job_ncbi_db.sh "
                    f"{self.job_id} "
                    f"{self.min_pair_dist} "
                    f"{self.max_pair_dist} "
                    f"{self.genome_start} "
                    f"{self.db_type} "
                    f"{self.input_seq_nuc}"
                )

            vprint(call_str, verbose_level=3)
            out = _ssh_call(call_str)
            vprint(out, verbose_level=3)

        self.submitted = True

    def _spawn_child(self, iterator, start):
        attributes = self.get_attributes()
        len_parent_seq = len(self.input_seq_nuc)

        attributes["job_hierarchy"] = "child"
        attributes["email_notification"] = False
        attributes["email"] = "Nobody"

        attributes["job_id"] = self.job_id + f"-child_{iterator+1}"
        attributes["status_dic"] = blank_job_submission()["status_dic"]

        # Important. When child calls submit method make job submission and not
        # start building custom db
        if attributes["custom_db"]:
            attributes["status_dic"]["customDB"] = self.status_dic["customDB"]

        attributes["path"] = JOB_SUBMISSION_PATH_TEMPLATE.format(attributes["job_id"])
        attributes["lock_path"] = attributes["path"] + ".lock"
        attributes["submitted"] = False
        attributes["finished"] = False
        attributes["got_full_results"] = False
        attributes["error_out"] = ""
        attributes["hss_table"] = []
        attributes["protein_list"] = []
        attributes["alig_info_table"] = []
        attributes["alignment_print"] = ""

        stop = start + CHILD_LENGTH
        stop = min(stop, len_parent_seq)

        attributes["genome_start"] = start
        attributes["genome_stop"] = stop
        attributes["input_seq_nuc"] = self.input_seq_nuc[start:stop]
        with JobSubmission(attributes=attributes) as child_job:
            child_job.submit()

    def _get_results(self):
        """Copy RNAcode results into results folder."""
        vprint(f"Get results for {self.job_id}", verbose_level=2)
        if self.status_dic["jobStatus"][1] in ["selection", "blast"]:
            return
        if self.job_hierarchy != "parent":
            call_str = (
                "ssh-scp-api/get_results.sh "
                f"{self.job_id} "
                f"{RNACODE_RESULT_PATH_TEMPLATE.format(self.job_id)} "
                f"{EPS_PATH_TEMPLATE.format(self.job_id)} "
                f"{ALIGN_PATH_TEMPLATE.format(self.job_id)} "
                f"{ALIGN_PLOT_PATH_TEMPLATE.format(self.job_id)} "
                f"{RESULT_DIR}"
            )
            vprint("Download results", verbose_level=3)
            vprint(call_str, verbose_level=3)
            _ssh_call(call_str)
        if self.status_dic["jobStatus"][0] == "CD":
            self.get_hss_table()
            self.get_protein_list()

        self.get_alig_info_table()
        self.load_alignment_print()

    def _read_error_call(self):
        call_str = f"ssh-scp-api/read_error.sh {self.job_id}"
        vprint(call_str, verbose_level=3)
        self.error_out = _ssh_call(call_str)

    def read_error(self):
        """Read error output from backend.

        This method needs that the job submission is locked or finished already,
        hence should be used with 'with ... as ...' syntax.
        """
        vprint(f"Read error job {self.job_id}", verbose_level=2)
        if self.finished:
            return self.error_out
        if self.read_only:
            raise AttributeError("job submission has no write privilege")
        if self.job_hierarchy == "parent" and self.custom_db and not self.submitted:
            vprint("Collect errors for build custom DB.")
            self._read_error_call()
        elif self.job_hierarchy == "parent":
            vprint(
                f"Collect errors from child processes {self.job_id}", verbose_level=3
            )
            for child_job_id, status in self.status_dic.items():
                if child_job_id in ["jobStatus", "customDB"]:
                    continue
                if status[0] != "E":
                    continue
                with JobSubmission(job_id=child_job_id) as child_job:
                    child_job = JobSubmission(job_id=child_job_id)
                    self.error_out += f"Error {child_job_id}\n"
                    self.error_out += child_job.read_error()
        else:
            self._read_error_call()
        return self.error_out

    def _get_full_job_status(self):
        """Reduce the states of each step or each child to a state for the whole job ("CD", "F, "E") and fail_reason."""
        if self.job_hierarchy == "parent":
            job_status = "R"
            fail_reason = 0

            if self.custom_db and not self.submitted:
                job_status = self.status_dic["customDB"]["fullJob"][0]
                fail_reason = self.status_dic["customDB"]["fullJob"][1]
                if job_status == "CD":
                    job_status = "R"
            elif all(
                (
                    status[0] == "F"
                    for job_id, status in self.status_dic.items()
                    if job_id not in ["jobStatus", "customDB"]
                )
            ):
                job_status = "F"
            elif any(
                (
                    status[0] == "E"
                    for job_id, status in self.status_dic.items()
                    if job_id not in ["jobStatus", "customDB"]
                )
            ):
                job_status = "E"
            elif all(
                (
                    status[0] in ["F", "CD", "E"]
                    for job_id, status in self.status_dic.items()
                    if job_id not in ["jobStatus", "customDB"]
                )
            ):
                job_status = "CD"
        else:
            job_status = self.status_dic["fullJob"][0]
            fail_reason = self.status_dic["fullJob"][1]

        self.status_dic["jobStatus"] = [job_status, fail_reason]
        vprint(f"Job status is: {job_status}; {fail_reason}", verbose_level=3)

    def _check_status_call(self):
        # simplifies state for website
        state_map = {
            "R": "R",
            "CD": "CD",
            "CG": "R",
            "F": "F",
            "PD": "R",
            "S": "F",
            "ST": "F",
            "NS": "NS",
            "E": "E",
        }
        call_str = f"ssh-scp-api/check_status.sh {JOB_STATUS_FILE_TEMPLATE.format(self.job_id)}"
        vprint(call_str, verbose_level=3)
        out = _ssh_call(call_str)
        try:
            status_dic = {k: [state_map[v[0]], v[1]] for k, v in json.loads(out).items()}
        except json.decoder.JSONDecodeError:
            vprint("Json decoder Error!")
            vprint(out)
            raise
        vprint(f"Backend status is:\n{json.dumps(status_dic, indent=2)}", verbose_level=3)
        if self.custom_db and self.job_hierarchy == "parent":
            self.status_dic["customDB"] = status_dic
        else:
            self.status_dic.update(status_dic)

    def check_status(self):
        """Check status of job.

        This method needs that the job submission is locked or finished already,
        hence should be used with 'with ... as ...' syntax.
        """
        vprint(f"Check status {self.job_id}", verbose_level=2)
        if self.finished:
            return self.status_dic
        if self.read_only:
            raise AttributeError("job submission has no write privilege")
        # First check status of custom db
        if self.custom_db and not self.submitted:
            self._check_status_call()
        elif self.job_hierarchy == "parent":
            vprint(f"Check status of child processes {self.job_id}", verbose_level=3)
            child_status_dic = {}
            for child_job_id in self.status_dic.keys():
                if child_job_id in ["jobStatus", "customDB"]:
                    continue
                vprint(f"Child: {child_job_id}", verbose_level=3)
                child_job = JobSubmission(job_id=child_job_id)
                with JobSubmission(job_id=child_job_id) as child_job:
                    child_status_dic[child_job_id] = child_job.check_status()[
                        "jobStatus"
                    ]
            self.status_dic.update(child_status_dic)
        else:
            self._check_status_call()

        self._get_full_job_status()
        if self.status_dic["jobStatus"][0] not in ["R", "NS"] and not self.finished:
            vprint("Not longer running", verbose_level=2)
            if self.status_dic["jobStatus"][0] == "E":
                self.read_error()
                self.finished = True
            # If custom DB parent job is checked and the job was not submitted
            # the custom DB was build.
            elif (
                self.custom_db and self.job_hierarchy == "parent" and not self.submitted
            ):
                if self.status_dic["jobStatus"][0] == "F":
                    vprint("Custom DB failed!", verbose_level=3)
                    self.submitted = True
                    self.finished = True
            else:
                self._get_results()
                self.finished = True
        return self.status_dic

    def download_full_results(self):
        """Create zip file of the complete job folder.

        This method needs that the job submission is locked or finished already,
        hence should be used with 'with ... as ...' syntax.
        """
        vprint(f"Get full results {self.job_id}", verbose_level=2)
        if not self.finished:
            raise AttributeError(f"job {self.job_id} did not finished.")
        if self.got_full_results:
            vprint("Full results already generated.", verbose_level=3)
            return
        if self.read_only:
            raise AttributeError("job submission has no write privilege")

        if self.job_hierarchy == "parent":
            if self.custom_db:
                vprint("Tar directory where custom db was build", verbose_level=3)
                call_str = f"ssh-scp-api/download_full_results.sh {self.job_id} {WORK_DIR} {RESULT_DIR}"
                _ssh_call(call_str)
            rnacode_full_results_path = FULL_RESULTS_PATH_TEMPLATE.format(self.job_id)
            vprint("Tar child directorys", verbose_level=3)
            # Tar in to tmp file move afterwards. Needed for custom DB jobs.
            call_str = f"tar -cvf {rnacode_full_results_path}.tmp --directory={RESULT_DIR} "
            if self.custom_db:
                call_str += f" {rnacode_full_results_path.split('/')[-1]}"

            for child_job_id in self.status_dic:
                if child_job_id in ["jobStatus", "customDB"]:
                    continue
                with JobSubmission(child_job_id) as child_job:
                    child_job.download_full_results()
                child_rnacode_full_results_path = FULL_RESULTS_PATH_TEMPLATE.format(
                    child_job_id
                ).split("/")[-1]
                call_str += f" {child_rnacode_full_results_path}"

            vprint("Tar result files", verbose_level=3)
            try:
                subprocess.run(
                    call_str.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
                )
            except subprocess.CalledProcessError as exc:
                vprint("error creating full result tar")
                vprint(f"Command\n{call_str}")
                vprint(f"stdout:\n{exc.stdout.decode('UTF8')}")
                vprint(f"stderr:\n{exc.stderr.decode('UTF8')}")
                raise
            vprint(rnacode_full_results_path)
            os.rename(rnacode_full_results_path + ".tmp", rnacode_full_results_path)
        else:
            vprint("Tar files", verbose_level=3)
            call_str = f"ssh-scp-api/download_full_results.sh {self.job_id} {WORK_DIR} {RESULT_DIR}"
            _ssh_call(call_str)

        self.got_full_results = True

    def filter_hss_table(self, p_threshold, best_overlap):
        """Filter the hss table by p-val and select best overlapping segments."""
        vprint("Filter hss table", verbose_level=2)
        segments = []
        if not self.finished:
            return segments

        for line in self.get_hss_table():
            if float(line[10]) > p_threshold:
                continue

            if not best_overlap:
                # Make list of list so it can be merged in place
                segments.append([line])
                continue
            for segment in segments:
                if (
                    # overlap
                    int(line[7]) <= int(segment[0][8])
                    and int(segment[0][7]) <= int(line[8])
                ):
                    if float(segment[0][10]) > float(line[10]):
                        segment[0] = line
                    break
            else:
                segments.append([line])
        # Unpack list of lists
        return [seg[0] for seg in segments]

    def get_hss_table(self):
        """Return information for High Scoring Segment table.

        This method needs that the job submission is locked or finished already,
        hence should be used with 'with ... as ...' syntax.
        """
        vprint(f"Get HSS table for job {self.job_id}", verbose_level=2)
        if self.finished:
            return self.hss_table
        if self.status_dic["jobStatus"][0] != "CD":
            return []
        if self.read_only:
            raise AttributeError("job submission has no write privilege")
        self.hss_table = []
        if self.job_hierarchy == "parent":
            vprint("Collect HSS table for all children", verbose_level=3)
            hss_table = []
            for child_job_id, child_status in self.status_dic.items():
                if child_job_id in ["jobStatus", "customDB"]:
                    continue
                if child_status[0] != "CD":
                    continue
                child_iterator = child_job_id.split("-")[-1]
                with JobSubmission(child_job_id) as child_job:
                    child_hss_table = [
                        [f"{child_iterator}_hss_{line[0]}"] + line[1:12]
                        for line in child_job.get_hss_table()
                    ]
                    hss_table += child_hss_table
            self.hss_table = merge_hss_table(hss_table)
        else:
            rnacode_result_path = RNACODE_RESULT_FRONTEND_PATH_TEMPLATE.format(
                self.job_id
            )
            with open(rnacode_result_path, "r", encoding="UTF-8") as f_handle:
                for line in f_handle.read().split("\n")[:-1]:
                    line = line.split("\t")
                    start = int(line[7])
                    frame = (start % 3) + 1
                    line[2] = str(frame)
                    self.hss_table.append(
                        line + [[f"results/{self.job_id}_hss-{line[0]}.eps"]]
                    )

        return self.hss_table

    def get_alig_info_table(self):
        """Return information for alignment info table.

        This method needs that the job submission is locked or finished already,
        hence should be used with 'with ... as ...' syntax.
        """
        vprint(f"Get alignment info table for job {self.job_id}")
        if self.finished:
            return self.alig_info_table
        if self.read_only:
            raise AttributeError("job submission has no write privilege")
        if self.job_hierarchy == "parent":
            vprint("Collect alignment info table for all children")
            self.alig_info_table = []
            for child_job_id in self.status_dic.keys():
                if child_job_id in ["jobStatus", "customDB"]:
                    continue
                child_iterator = child_job_id.split("-")[-1]
                with JobSubmission(child_job_id) as child_job:
                    self.alig_info_table.append(
                        [child_iterator] + child_job.get_alig_info_table()
                    )
        else:
            job_status = self.status_dic["jobStatus"][0]
            if self.status_dic["jobStatus"][0] not in ["F", "CD"]:
                self.alig_info_table = [0, "", job_status]
                return self.alig_info_table
            if self.status_dic["jobStatus"][1] not in [0, "no_res"]:
                self.alig_info_table = [0, "", job_status]
                return self.alig_info_table
            alignment_fasta_path = ALIG_FASTA_PATH_TEMPLATE.format(self.job_id)
            if os.path.isfile(alignment_fasta_path):
                with open(alignment_fasta_path, "r", encoding="UTF-8") as f_handle:
                    species_list = [
                        li.split()[0][1:] for li in f_handle if li[0] == ">"
                    ]
                    species_list = [
                        " ".join(s.split("-")[0].split("_"))
                        for s in species_list
                        if s != "Target"
                    ]
                self.alig_info_table = [
                    len(species_list),
                    ", ".join(species_list),
                    job_status,
                ]
            else:
                self.alig_info_table = [0, "", job_status]

        return self.alig_info_table

    def get_protein_list(self):
        """Return information for High Scoring Segment table.

        This method needs that the job submission is locked or finished already,
        hence should be used with 'with ... as ...' syntax.
        """
        vprint(f"Get protein list for job {self.job_id}", verbose_level=2)
        if self.finished:
            return self.protein_list
        if self.status_dic["jobStatus"][0] != "CD":
            return self.protein_list
        if self.read_only:
            raise AttributeError("job submission has no write privilege")
        self.protein_list = []
        if self.job_hierarchy == "child":
            parent_job_id = "-".join(self.job_id.split("-")[:-1])
            parent_job = JobSubmission(parent_job_id)
            input_seq_nuc = parent_job.input_seq_nuc
        else:
            input_seq_nuc = self.input_seq_nuc

        for entry in self.get_hss_table():
            hss_id = entry[0]
            strand = entry[1]
            start = int(entry[7])
            end = int(entry[8])

            if strand == "+":
                hss_seq_nuc = Seq(input_seq_nuc[start:end+1])
            else:
                hss_seq_nuc = Seq(input_seq_nuc[start:end+1]).reverse_complement()
            protein = str(hss_seq_nuc.translate())
            self.protein_list.append([hss_id, protein])
        return self.protein_list

    def plot_segment_diagram(self, p_threshold, best_overlap):
        """Use R to draw a plot for all hss segments."""
        vprint("Plot segment diagramm", verbose_level=2)
        if not self.finished or self.status_dic["jobStatus"][0] != "CD":
            vprint("Not finished", verbose_level=3)
            return
        if self.last_plot_par == [p_threshold, best_overlap]:
            vprint("Already plotted", verbose_level=3)
            return
        if not self.read_only:
            vprint(
                "This method will open its object. This might overwrite previous changes.",
                verbose_level=1,
            )
            raise AttributeError("job submission has write privilege.")

        vprint("Finished and needs to plot", verbose_level=3)

        with self:
            self.last_plot_par = [p_threshold, best_overlap]

        hss_r_table_path = HSS_TABLE_R_PATH_TEMPLATE.format(self.job_id)
        hss_plot_path = HSS_PLOT_FRONTEND_PATH_TEMPLATE.format(self.job_id)

        if self.job_hierarchy == "child":
            parent_job_id = "-".join(self.job_id.split("-")[:-1])
            parent_job = JobSubmission(parent_job_id)
            input_seq_nuc = parent_job.input_seq_nuc
        else:
            input_seq_nuc = self.input_seq_nuc

        vprint("Merge overlapping segments", verbose_level=3)
        segments = self.filter_hss_table(p_threshold, best_overlap)

        if len(segments) == 0:
            vprint("No hss to plot", verbose_level=2)
            if os.path.isfile(hss_plot_path):
                os.remove(hss_plot_path)
            return

        vprint("Build table for R", verbose_level=3)
        r_table_list = _build_r_table(segments)

        with open(hss_r_table_path, "w", encoding="UTF-8") as f_handle:
            for line in r_table_list:
                f_handle.write(",".join([repr(e) for e in line]) + "\n")

        if self.job_hierarchy == "parent":
            max_len_no_split = CHILD_LENGTH
        else:
            max_len_no_split = len(input_seq_nuc)

        call_str = (
            f"{R_BIN}/Rscript ./plot_hss_segments.R {hss_r_table_path} "
            f"{hss_plot_path} {len(input_seq_nuc)} "
            f"{max_len_no_split}"
        )
        try:
            subprocess.run(
                call_str.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
            )
        except subprocess.CalledProcessError as exc:
            vprint("error R script")
            vprint(f"Command\n{call_str}")
            vprint(f"stdout:\n{exc.stdout.decode('UTF8')}")
            vprint(f"stderr:\n{exc.stderr.decode('UTF8')}")
            raise

    def load_alignment_print(self):
        """Read html segment generated by mview and returns string for jinja context.

        This method needs that the job submission is locked or finished already,
        hence should be used with 'with ... as ...' syntax.
        """
        vprint("Load html for alingment", verbose_level=2)
        if self.job_hierarchy == "parent":
            return ""
        if self.finished:
            return self.alignment_print
        if self.status_dic["jobStatus"][0] not in ["F", "CD"]:
            self.alignment_print = ""
            return self.alignment_print
        if self.status_dic["jobStatus"][1] not in [0, "no_res"]:
            self.alignment_print = ""
            return self.alignment_print
        if self.read_only:
            raise AttributeError("job submission has no write privilege")
        alig_file_path = ALIG_FILE_PATH_TEMPLATE.format(self.job_id)
        with open(alig_file_path, "r", encoding="UTF-8") as f_handle:
            alignment_print = f_handle.read()

        self.alignment_print = alignment_print.replace("#ffff00", "#f04258")
        return self.alignment_print
