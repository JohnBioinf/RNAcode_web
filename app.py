#!/usr/bin/python3
"""Flask frontend web service for RNAcode_web."""

# Basic python modules
import json
import os
import random
import string
from datetime import datetime
from ipaddress import ip_network, ip_address
import multiprocessing
import traceback
import functools

# email
from email.message import EmailMessage
import smtplib

# Flask
from flask import Flask, render_template, request, redirect, url_for, Response
from flask import send_file, send_from_directory

# Limiting access of users
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

# Background cleaning
from apscheduler.schedulers.background import BackgroundScheduler
from apscheduler.triggers.interval import IntervalTrigger

from JobSubmission import vprint
from JobSubmission import LOG_DIR
from JobSubmission import RESULT_DIR
from JobSubmission import JobSubmission
from JobSubmission import JobNotFoundException
from JobSubmission import AlreadySubmittedException
from JobSubmission import JOB_SUBMISSION_PATH_TEMPLATE
from JobSubmission import FULL_RESULTS_PATH_TEMPLATE
from JobSubmission import check_bulk_json
from JobSubmission import read_bulk_json
from JobSubmission import generate_bulk_submission_example

from RNAcodeWebCore import P_THRESHOLD


# Parameters
with open("./parameters_frontend_local.json", "r", encoding="UTF-8") as file_handle:
    server_parameters_frontend = json.load(file_handle)

with open("./parameters_backend_local.json", "r", encoding="UTF-8") as file_handle:
    server_parameters_backend = json.load(file_handle)

SMTP_SERVER = server_parameters_frontend["smtp_server"]
RNACODE_EMAIL_ADDRESS = server_parameters_frontend["rnacode_email_address"]
ADMIN_EMAIL = server_parameters_frontend["admin_email_address"]
RNACODE_WEB_REPO_BACKEND = server_parameters_backend["RNAcode_web_repo_backend"]
IP_FRONTEND = server_parameters_frontend["ip"]
PORT_FRONTEND = server_parameters_frontend["port"]

# Names of jobs used in publication. Will be excluded from clean up.
with open("./static/jobs_publication.json", "r", encoding="UTF-8") as file_handle:
    STATIC_JOB_IDS = [e["job_id"] for e in json.load(file_handle)]


def handle_internal_exception(route_function):
    """Decorate a route to handle internal exceptions.

    This should be use as an decorator for all functions that could fail.
    Especially all routes with functionality and might throw an error.
    Send an email to the admin or print the log to console if in debug mode.
    If an internal error is handled handle_exception() must be called from the
    same process that caused the error. Or else not the correct log will be
    found.
    """
    @functools.wraps(route_function)
    def wrapper(*args, **kwargs):
        try:
            return route_function(*args, **kwargs)
        except Exception:
            vprint("Handle exception", verbose_level=2)
            error = traceback.format_exc()
            with open(f"{LOG_DIR}/{os.getpid()}.log", "r", encoding="UTF-8") as f_handle:
                logs = f_handle.read()
                logs = "\n".join(logs.split("\n")[-200:])

            content = (
                f"Unexpected error in {route_function.__name__}:\n"
                f"{error}\n"
                f"PID={os.getpid()}\n"
                f"{logs}\n"
            )

            if app.debug:
                print(content)
                raise

            subject = "RNAcode Web encounterned an error"
            send_mail(subject, ADMIN_EMAIL, content)
            return render_template("html/errors/internal_error.html")

    return wrapper


@handle_internal_exception
def clean_up():
    """Delete data from system, should be called by scheduler."""
    vprint("Start clean up", verbose_level=1)
    today = datetime.today()
    job_submission_dir = os.path.dirname(JOB_SUBMISSION_PATH_TEMPLATE.format("none"))
    for job_submission_file in os.listdir(job_submission_dir):
        if job_submission_file.split(".")[-1] != "json":
            continue
        job_id = job_submission_file.replace(".json", "")
        # Ignore jobs from publication.
        is_static_job = False
        for static_job_id in STATIC_JOB_IDS:
            # job_id can be static_job_id_child_n
            if static_job_id in job_id:
                is_static_job = True
                continue
        if is_static_job:
            continue

        with JobSubmission(job_id) as job_submission:
            submission_date = datetime.strptime(
                job_submission.date, DATE_FORMATION
            )
            # clean up jobs older than storage time
            if (today - submission_date).days > JOB_STORAGE_TIME:
                vprint(f"Delete job {job_id} is old", verbose_level=2)
                job_submission.delete()
            # clean up jobs that were not submitted and are older than 1 day
            if not job_submission.submitted and (today - submission_date).days > 1:
                vprint(f"Delete job {job_id} was not submitted", verbose_level=2)
                job_submission.delete()

    vprint("Delete left over filelocks", verbose_level=1)
    for lock_file in os.listdir(job_submission_dir):
        if lock_file.split(".")[-1] != "lock":
            continue
        job_submission_file = os.path.join(
            job_submission_dir, lock_file.replace(".lock", "")
        )
        if os.path.isfile(job_submission_file):
            continue
        vprint(f"Lock file {lock_file} is left over", verbose_level=2)
        os.remove(os.path.join(job_submission_dir, lock_file))


# Storage time of jobs in days
JOB_STORAGE_TIME = 80

# Job clean up
scheduler = BackgroundScheduler(daemon=True, timezone="Europe/Berlin")
trigger = IntervalTrigger(hours=24)
scheduler.add_job(func=clean_up, trigger=trigger)
# python3.9+
# scheduler.add_job(func=clean_up, trigger="interval", hours=24)
scheduler.start()

DATE_FORMATION = "%B %d, %Y"

app = Flask(__name__)
app.config["UPLOAD_FOLDER"] = "./tmp/"
app.config["MAX_CONTENT_LENGTH"] = 16 * 1024 * 1024  # for 16 MB max-limit
app.jinja_env.trim_blocks = True
app.jinja_env.lstrip_blocks = True

# Limiter parameters
# bulk submission
BULK_TIME_LIMIT = "1 per day"


# Currently limit is only set for bulk submissions
@app.errorhandler(429)
def rate_limit_handler(error):
    """Create limit error html."""
    return render_template("html/errors/bulk_rate_limit.html", BULK_TIME_LIMIT=BULK_TIME_LIMIT)


# Delay in seconds when a result page should be reloaded.
RELOAD_DELAY = 240
UNI_NETS = [ip_network(server_parameters_frontend["uni_net_space"])]
LIMITER = Limiter(
    app,
    key_func=get_remote_address,
    storage_uri="redis://localhost:6379",
    # default limits also limit files in static should be not the case. Maybe an
    # update would fix this.
    # default_limits=["200 per day", "50 per hour"]
)

with open(
    "./templates/emails/successful_email.txt", "r", encoding="UTF-8"
) as file_handle:
    SUCCESSFUL_EMAIL = file_handle.read()

with open("./templates/emails/failure_email.txt", "r", encoding="UTF-8") as file_handle:
    FAILURE_EMAIL = file_handle.read()

with open("./templates/emails/error_email.txt", "r", encoding="UTF-8") as file_handle:
    ERROR_EMAIL = file_handle.read()

with open(
    "./templates/emails/error_email_job_admin.txt", "r", encoding="UTF-8"
) as file_handle:
    ERROR_EMAIL_JOB_ADMIN = file_handle.read()

with open(
    "./templates/emails/error_email_internal_admin.txt", "r", encoding="UTF-8"
) as file_handle:
    ERROR_EMAIL_INTERNAL_ADMIN = file_handle.read()

with open(
    "./templates/emails/submission_email.txt", "r", encoding="UTF-8"
) as file_handle:
    SUBMISSION_EMAIL = file_handle.read()

with open(
    "./templates/emails/bulk_submission_email.txt", "r", encoding="UTF-8"
) as file_handle:
    BULK_SUBMISSION_EMAIL = file_handle.read()

# Level of how much stuff should be logged. Either 0, 1, 2, 3.
VERBOSE_LEVEL = 3


def handle_job_error(job_submission):
    """Handle notifcation if backend job had an error.

    Send an email to the admin or print the stderr from the process to console
    if in debug mode.
    """
    vprint("Handle job error", verbose_level=2)
    content = (
        f"Unexpected error in {job_submission.job_id}:\n"
        f"{job_submission.read_error()}\n"
    )

    if app.debug:
        print(content)
    else:
        subject = "RNAcode Web encounterned an error"
        send_mail(subject, ADMIN_EMAIL, content)


def visitor_has_no_restriction():
    """Check if visitor is from specified range list."""
    return any((ip_address(request.remote_addr) in net for net in UNI_NETS))


def send_mail(subject, recipient_email_address, content):
    """Send mail function."""
    msg = EmailMessage()
    msg.set_content(content)
    msg["Subject"] = subject
    msg["From"] = RNACODE_EMAIL_ADDRESS
    msg["To"] = recipient_email_address
    smtp_handle = smtplib.SMTP(SMTP_SERVER)
    smtp_handle.send_message(msg)
    smtp_handle.quit()


def bulk_notification(email, job_id_list):
    """Send notification for a bulk submission."""
    overview_jobs_url = url_for(
        "list_bulk", job_id_list=";".join(job_id_list), _external=True
    )
    content = BULK_SUBMISSION_EMAIL.format(len(job_id_list), overview_jobs_url)
    subject = "Jobs submitted to RNAcode Web"
    send_mail(subject, email, content)


def send_notification(job_submission):
    """Send notification depending on job status."""
    if job_submission.email == "Nobody":
        return

    status = job_submission.status_dic["jobStatus"][0]
    job_id = job_submission.job_id
    vprint(f"Send email to user for {job_id} with status {status}", verbose_level=2)
    url = url_for("submission_by_url", job_id=job_id, _external=True)
    if status == "CD":
        content = SUCCESSFUL_EMAIL.format(job_id, url)
        subject = f"RNAcode Web finished successfully {job_id}"
    elif status == "F":
        content = FAILURE_EMAIL.format(job_id, url)
        subject = f"RNAcode Web finished unsucessfully {job_id}"
    elif status == "E":
        handle_job_error(job_submission)
        content = ERROR_EMAIL.format(job_id, url)
        subject = f"RNAcode Web finished with an error {job_id}"
    elif status == "R":
        content = SUBMISSION_EMAIL.format(job_id, url)
        subject = "Job submitted to RNAcode Web"

    send_mail(subject, job_submission.email, content)


@handle_internal_exception
def submit_job(job_id, notification=True):
    """Start job submission in background.

    Should be started in a new thread. If the job is building a custom DB
    submit_job() becomes the check_status() function for the build until the
    child jobs are submitted.
    """
    vprint(f"Submitting job {job_id} in a new thread.", verbose_level=1)
    try:
        with JobSubmission(job_id) as job_submission:
            job_finished_before_check = job_submission.finished
            job_submission.submit()
            if job_submission.email_notification and notification and not job_submission.submission_email_send:
                send_notification(job_submission)
            job_submission.submission_email_send = True
            if job_finished_before_check != job_submission.finished:
                vprint("Job finished!", verbose_level=1)
                send_notification(job_submission)
    except AlreadySubmittedException:
        pass


@handle_internal_exception
def check_status(job_id):
    """Check status of job and if job finished send notifications.

    Should be started in a new thread.
    """
    vprint(f"Check status for job {job_id} and send notifications.", verbose_level=1)
    parent_job_id = None
    with JobSubmission(job_id=job_id) as job_submission:
        job_finished_before_check = job_submission.finished
        job_submission.check_status()
        vprint("Status checked.", verbose_level=2)

        if job_finished_before_check != job_submission.finished:
            vprint("Job finished!", verbose_level=1)
            send_notification(job_submission)
            if job_submission.job_hierarchy == "child":
                parent_job_id = "-".join(job_submission.job_id.split("-")[:-1])
    if parent_job_id is not None:
        vprint("Notify parent", verbose_level=2)
        check_status(parent_job_id)


def show_bulk(job_id_list):
    """Show progress of bulk jobs in list."""
    job_id_list = job_id_list.split(";")
    job_info_list = []
    context = {}
    for job_id in job_id_list:
        try:
            job_submission = JobSubmission(job_id=job_id)
        except JobNotFoundException:
            return render_template("html/errors/submission_not_found.html", **locals())

    for job_id in job_id_list:
        job_submission = JobSubmission(job_id=job_id)
        thread = multiprocessing.Process(
            target=check_status, args=(job_submission.job_id,)
        )
        thread.start()
        status = job_submission.status_dic["jobStatus"][0]

        job_info_list.append([job_id, status])

    context["job_info_list"] = job_info_list
    return render_template("html/results/bulk_results.html", **context)


def submit_bulk_jobs(job_id_list):
    """Submit bulk jobs and send bulk notification."""
    job_id_list = job_id_list.split(";")
    email_list = []
    for job_id in job_id_list:
        try:
            job_submission = JobSubmission(job_id=job_id)
        except JobNotFoundException:
            return
        email_list.append(job_submission.email)
        if job_submission.submitted:
            return

    # Checks if all emails are the same. If so only send one notification
    # instead of multiples.
    if len(set(email_list)) == 1 and len(email_list) != 1:
        notification = False
        send_bulk_notification = True
    else:
        notification = True
        send_bulk_notification = False

    for job_id in job_id_list:
        job_submission = JobSubmission(job_id=job_id)
        k_arg_dic = {}
        k_arg_dic["notification"] = notification
        thread = multiprocessing.Process(
            target=submit_job, args=(job_submission.job_id,), kwargs=k_arg_dic
        )
        thread.start()

    if send_bulk_notification:
        bulk_notification(email_list[0], job_id_list)


def set_context_child_status(status_dic):
    """Set context of children status, for result website."""
    context = {}
    context["num_failed"] = count_num_status_child(status_dic, "F")
    context["num_success"] = count_num_status_child(status_dic, "CD")
    context["num_running"] = count_num_status_child(status_dic, "R")
    context["num_error"] = count_num_status_child(status_dic, "E")
    context["num_waiting"] = count_num_status_child(status_dic, "NS")
    context["num_children"] = (
        context["num_error"]
        + context["num_failed"]
        + context["num_success"]
        + context["num_running"]
        + context["num_waiting"]
    )
    return context


def count_num_status_child(status_dic, status_type):
    """Count number of children that have a certain status."""
    return len(
        [
            1
            for name, status in status_dic.items()
            if name not in ["jobStatus", "customDB"] and status[0] == status_type
        ]
    )


def check_p_threshold(p_threshold):
    """Check p threshold."""
    try:
        p_threshold = float(p_threshold)
        p_threshold_error = ""
    except ValueError:
        p_threshold_error = "Not float"
        p_threshold = P_THRESHOLD
    if p_threshold > P_THRESHOLD:
        p_threshold_error = "To high"
        p_threshold = P_THRESHOLD
    if p_threshold <= 0:
        p_threshold_error = "Negativ"
        p_threshold = P_THRESHOLD
    return p_threshold, p_threshold_error


def build_result_webpage(job_submission, form):
    """Build website for result page.

    :return: html site
    :rtype: str
    """
    vprint(f"Build status website for job {job_submission.job_id}.", verbose_level=1)
    context = {}
    context.update(job_submission.get_attributes())
    if len(form) != 0:
        p_threshold, p_threshold_error = check_p_threshold(form["p_threshold"])
        best_overlap = "best_overlap" in form
    else:
        p_threshold = P_THRESHOLD
        p_threshold_error = ""
        best_overlap = False

    context["p_threshold"] = p_threshold
    context["best_overlap"] = best_overlap
    context["p_threshold_error"] = p_threshold_error
    context["P_THRESHOLD"] = P_THRESHOLD

    context["hss_table"] = job_submission.filter_hss_table(p_threshold, best_overlap)
    hss_id_list = [line[0] for line in context["hss_table"]]
    context["protein_list"] = [
        entry for entry in job_submission.get_protein_list() if entry[0] in hss_id_list
    ]

    job_submission.plot_segment_diagram(p_threshold, best_overlap)

    status_dic = job_submission.status_dic
    context["status"] = status_dic["jobStatus"][0]
    context["fail_reason"] = status_dic["jobStatus"][1]
    if context["status"] == "CD":
        context["show_results"] = True
    elif context["fail_reason"] == "no_res":
        # This only catches no results by RNAcode. Alignment results should still
        # be shown.
        context["show_results"] = True
    else:
        context["show_results"] = False

    if job_submission.job_hierarchy == "parent":
        if job_submission.custom_db:
            num_iterations = sorted(
                [k.split("_")[1] for k in status_dic["customDB"] if "_" in k]
            )[-1]
            context["num_iterations"] = int(num_iterations)
            # reduces job states of blast and seqSel to highest iteration
            job_status = {
                k.split("_")[0]: v
                for k, v in status_dic["customDB"].items()
                if f"_{num_iterations}" in k or "_" not in k
            }
            context.update(**job_status)
        context.update(set_context_child_status(status_dic))
    else:
        num_iterations = sorted([k.split("_")[1] for k in status_dic if "_" in k])[-1]
        context["num_iterations"] = int(num_iterations)
        # reduces job states of blast and seqSel to highest iteration
        job_status = {
            k.split("_")[0]: v
            for k, v in status_dic.items()
            if f"_{num_iterations}" in k or "_" not in k
        }
        context.update(**job_status)
        context["alignment_print"] = job_submission.load_alignment_print()

    if context["status"] == "R":
        context["reload_delay"] = RELOAD_DELAY
    return render_template("html/results/full_results.html", **context)


def download_file(request_file):
    """Download file and return path to file."""
    vprint("Download file", verbose_level=1)
    letters = string.ascii_lowercase
    file_name = "".join(random.choice(letters) for i in range(10))
    file_path = os.path.join(app.config["UPLOAD_FOLDER"], file_name)
    while os.path.isfile(file_path):
        file_name = "".join(random.choice(letters) for i in range(10))
        file_path = os.path.join(app.config["UPLOAD_FOLDER"], file_name)
    request_file.save(file_path)
    return file_path


@app.route("/")
@app.route("/start")
def routing():
    """Present startpage."""
    return render_template("html/content/startpage.html")


@app.route("/bulk_submission/example_bulk_submission.json")
@handle_internal_exception
def get_example_bulk_submission():
    """Return an example json file for the bulk download."""
    vprint(
        "Visiting /bulk_submission/example_bulk_submission.json routing to get_example_bulk_submission",
        verbose_level=1,
    )
    json_content = generate_bulk_submission_example()
    return Response(
        json_content,
        mimetype="text/json",
        headers={
            "Content-disposition": "attachment; filename=example_bulk_submission.json"
        },
    )


@app.route("/submission/<job_id>.tgz")
@handle_internal_exception
def get_full_results(job_id):
    """Return the full work directory for a job as compressed file."""
    vprint(
        f"Visiting /submission/{job_id}.tgz routing to get_full_results",
        verbose_level=1,
    )
    rnacode_full_results_path = FULL_RESULTS_PATH_TEMPLATE.format(job_id)
    with JobSubmission(job_id=job_id) as job_submission:
        job_submission.download_full_results()

    return send_file(rnacode_full_results_path)


@app.route("/bulk_submission", methods=["POST", "GET"])
@handle_internal_exception
def bulk_submission():
    """Generate the side for submitting in bulk."""
    vprint(
        "Visiting /bulk_submission routing to bulk_submission",
        verbose_level=1,
    )
    bulk_check_dic = {}
    context = {}
    if len(request.files) == 0:
        context["bulk_check_dic"] = bulk_check_dic
        return render_template("html/input/bulk_submission.html", **context)

    json_file = request.files["json"]
    if json_file.filename == "":
        bulk_check_dic["no_file"] = "No file provided"
        context["bulk_check_dic"] = bulk_check_dic
        return render_template("html/input/bulk_submission.html", **context)

    json_file_path = download_file(json_file)

    # This only checks basic validity of json file not if the parameters for
    # each job submission are actually valid.
    bulk_check_dic = check_bulk_json(json_file_path)

    if bulk_check_dic["error"]:
        os.remove(json_file_path)
        context["bulk_check_dic"] = bulk_check_dic
        return render_template("html/input/bulk_submission.html", **context)

    # important json file still exist after leaving through this route. Make
    # sure that it will get deleted in the next route.
    return redirect(
        url_for("confirm_bulk_submission", filename=os.path.basename(json_file_path))
    )


@app.route("/confirm_bulk_submission/<filename>", methods=["POST", "GET"])
@handle_internal_exception
def confirm_bulk_submission(filename):
    """Check bulk submission arguments.

    Show reports of invalid arguments. If all arguments are valid will ask the
    user to submit bulk jobs.
    """
    vprint(
        f"Visiting /confirm_bulk_submission/{filename} routing to confirm_bulk_submission()",
        verbose_level=1,
    )
    download_path = os.path.join(app.config["UPLOAD_FOLDER"], filename)

    if not os.path.isfile(download_path):
        return render_template("html/errors/bulk_submission_not_found.html", **locals())

    with open(download_path, "r", encoding="UTF-8") as f_handle:
        job_list = read_bulk_json(f_handle)

    os.remove(download_path)

    no_errors = all((job.input_check_dic["no_error"] for job in job_list))
    job_id_list = []
    job_dic_list = []
    for job in job_list:
        job_id_list.append(job.job_id)
        job_dic_list.append(job.get_attributes())
        if not no_errors:
            with job:
                job.delete(backend=False, recursive=False)

    job_id_list = ";".join(job_id_list)
    return render_template(
        "html/confirm/confirm_bulk_submission.html",
        job_list=job_dic_list,
        job_id_list=job_id_list,
        no_errors=no_errors,
    )


@app.route("/change_bulk_submission/<job_id_list>", methods=["POST", "GET"])
@handle_internal_exception
def change_bulk_submission(job_id_list):
    """Remove all jobs that had been checked, redirect to new submission."""
    vprint(
        f"Visiting /change_bulk_submission/{job_id_list} routing to change_bulk_submission()",
        verbose_level=1,
    )
    job_id_list = job_id_list.split(";")
    for job_id in job_id_list:
        try:
            with JobSubmission(job_id=job_id) as job_submission:
                if job_submission.submitted:
                    return render_template(
                        "html/errors/job_already_submitted.html", **locals()
                    )
                job_submission.delete()
        except JobNotFoundException:
            return render_template("html/errors/submission_not_found.html", **locals())

    return redirect(url_for("bulk_submission"))


@app.route("/list_bulk/<job_id_list>", methods=["post", "get"])
@handle_internal_exception
def list_bulk(job_id_list):
    """Will list progress of jobs in list."""
    vprint(f"Visiting /list_bulk/{job_id_list} routing to list_bulk()", verbose_level=1)
    return show_bulk(job_id_list)


@app.route("/submit_bulk/<job_id_list>", methods=["post", "get"])
@LIMITER.limit(
    BULK_TIME_LIMIT,
    exempt_when=visitor_has_no_restriction,
)
@handle_internal_exception
def submit_bulk(job_id_list):
    """Will submit all jobs if not yet submitted."""
    vprint(
        f"Visiting /submit_bulk/{job_id_list} routing to submit_bulk()", verbose_level=1
    )
    submit_bulk_jobs(job_id_list)
    return redirect(url_for("list_bulk", job_id_list=job_id_list))


@app.route(
    "/single_submission/<job_id>/protein.fasta?p_threshold=<p_threshold>,best_overlap=<best_overlap>",
    methods=["POST", "GET"],
)
@handle_internal_exception
def get_protein_fasta(job_id, p_threshold, best_overlap):
    """Return protein fasta for job."""
    vprint(
        f"Visiting /single_submission/{job_id}/protein.fasta?p_threshold={p_threshold},best_overlap={best_overlap} routing to get_protein_fasta()",
        verbose_level=1,
    )
    try:
        job_submission = JobSubmission(job_id=job_id)
    except JobNotFoundException:
        return render_template("html/errors/submission_not_found.html", **locals())

    if job_submission.status_dic["jobStatus"][0] != "CD":
        vprint(
            f"Job {job_submission.job_id} did not finished sucesfully.", verbose_level=2
        )
        return redirect(url_for("submission_by_url", job_id=job_id))

    p_threshold, p_threshold_error = check_p_threshold(p_threshold)
    if p_threshold_error != "":
        return redirect(url_for("submission_by_url", job_id=job_id))
    if best_overlap == "True":
        best_overlap = True
    elif best_overlap == "False":
        best_overlap = False
    else:
        return redirect(url_for("submission_by_url", job_id=job_id))

    hss_id_list = [
        line[0] for line in job_submission.filter_hss_table(p_threshold, best_overlap)
    ]
    protein_list = [
        entry for entry in job_submission.get_protein_list() if entry[0] in hss_id_list
    ]
    fasta_content = ""
    for entry in protein_list:
        fasta_content += f">{entry[0]}\n{entry[1]}\n"
    fasta_name = job_id
    if p_threshold != P_THRESHOLD:
        fasta_name += f"_p-{p_threshold}"

    if best_overlap:
        fasta_name += "_only-best-overlap"

    fasta_name += "_protein.fasta"
    return Response(
        fasta_content,
        mimetype="text/fasta",
        headers={"Content-disposition": f"attachment; filename={fasta_name}"},
    )


@app.route("/single_submission/<job_id>", methods=["POST", "GET"])
@handle_internal_exception
def single_submission_edit(job_id):
    """Edit job submission input."""
    vprint(
        f"Visiting /single_submission/{job_id} routing to single_submission_edit()",
        verbose_level=1,
    )
    try:
        job_submission = JobSubmission(job_id=job_id)
        if job_submission.submitted:
            vprint(
                f"Job {job_submission.job_id} was already submitted.", verbose_level=2
            )
            return render_template(
                "html/errors/job_already_submitted.html", job_id=job_id
            )
    except JobNotFoundException:
        vprint("Set new job.", verbose_level=2)
        job_submission = JobSubmission()

    return render_template(
        "html/input/single_submission.html", **job_submission.get_attributes()
    )


@app.route("/single_submission", methods=["POST", "GET"])
@handle_internal_exception
def single_submission_fresh():
    """Return on first visit a fresh input site, or test the validity of the user input."""
    vprint(
        "Visiting /single_submission routing to single_submission_fresh()",
        verbose_level=1,
    )
    context = {}
    if len(request.form) == 0:
        vprint("New submision")
        job_submission = JobSubmission()
        context.update(job_submission.get_attributes())
        return render_template("html/input/single_submission.html", **context)
    vprint("Check input of submision", verbose_level=2)
    if request.files["fasta"].filename != "":
        fasta_file_path = download_file(request.files["fasta"])
    else:
        fasta_file_path = None

    request_job_id = request.form["job_id"]
    try:
        job_submission = JobSubmission(job_id=request.form["job_id"])
        if job_submission.submitted:
            vprint(
                f"Job {job_submission.job_id} was already submitted.", verbose_level=2
            )
            if fasta_file_path:
                os.remove(fasta_file_path)

            return render_template(
                "html/errors/job_already_submitted.html", job_id=request_job_id
            )

        vprint(f"Overwrite job {job_submission.job_id}.", verbose_level=2)
        with JobSubmission(job_id=request.form["job_id"]) as job_submission:
            job_submission.delete(backend=False, recursive=False)
    except JobNotFoundException:
        vprint("Set new job.", verbose_level=2)

    job_submission = JobSubmission(
        form=request.form,
        fasta_file=fasta_file_path,
    )

    if fasta_file_path:
        os.remove(fasta_file_path)

    context.update(job_submission.get_attributes())
    if job_submission.input_check_dic["no_error"]:
        return render_template(
            "html/confirm/confirm_single_submission.html",
            **context,
        )
    return render_template("html/input/single_submission.html", **context)


@app.route("/change_submission/<job_id>", methods=["POST", "GET"])
@handle_internal_exception
def change_submission(job_id):
    """Redirect user to make a new submission that already passed all checks."""
    vprint(
        f"Visiting /change_submission/{job_id} routing to change_submission",
        verbose_level=1,
    )
    try:
        job_submission = JobSubmission(job_id=job_id)
    except JobNotFoundException:
        return render_template("html/errors/submission_not_found.html", **locals())

    if job_submission.submitted:
        return render_template("html/errors/job_already_submitted.html", **locals())

    return redirect(url_for("single_submission_edit", job_id=job_id))


@app.route("/submission", methods=["POST", "GET"])
@handle_internal_exception
def submission_by_form():
    """Redirect route if submission gets posted by check submission route."""
    vprint("Visiting /submission routing to submission_by_form", verbose_level=1)
    if "job_id_list" not in request.form:
        return redirect(url_for("check"))
    job_id_list = request.form["job_id_list"]

    if job_id_list == "":
        return redirect(url_for("check"))

    if ";" in job_id_list:
        return redirect(url_for("list_bulk", job_id_list=job_id_list))

    return redirect(url_for("submission_by_url", job_id=job_id_list))


@app.route("/submission/<job_id>", methods=["POST", "GET"])
@handle_internal_exception
def submission_by_url(job_id):
    """Submit job if not yet submitted and/ or show progress of job."""
    vprint(
        f"Visiting /submission/{job_id} routing to submission_by_url", verbose_level=1
    )
    try:
        job_submission = JobSubmission(job_id=job_id)
    except JobNotFoundException:
        return render_template("html/errors/submission_not_found.html", **locals())

    if not job_submission.input_check_dic["no_error"]:
        return redirect(url_for("single_submission_edit", job_id=job_id))

    if not job_submission.submitted:
        thread = multiprocessing.Process(
            target=submit_job, args=(job_submission.job_id,)
        )
        thread.start()
    else:
        thread = multiprocessing.Process(
            target=check_status, args=(job_submission.job_id,)
        )
        thread.start()

    return build_result_webpage(job_submission, request.form)


@app.route("/privacy")
@LIMITER.limit(
    BULK_TIME_LIMIT,
    exempt_when=visitor_has_no_restriction,
)
def privacy():
    """Return privacy notes."""
    return render_template("html/content/privacy.html")


@app.route("/check")
def check():
    """Return form that user can check her submission."""
    return render_template("html/input/check.html")


@app.route("/documentation")
def documentation():
    """Return documentation of web service."""
    return render_template("html/content/documentation.html")


@app.route("/favicon.ico")
def favicon():
    """Serve favicon."""
    return send_from_directory(
        os.path.join(app.root_path, "static"),
        "favicon.ico",
        mimetype="image/vnd.microsoft.icon",
    )


@app.route("/results/<file_name>")
@handle_internal_exception
def result_file(file_name):
    """Serve result files."""
    vprint("Visiting /results/{file_name} to result_file()", verbose_level=1)
    return send_from_directory(RESULT_DIR, file_name)


@app.route("/test_error_admin")
@handle_internal_exception
def test_error_admin():
    """Test route to check error handling."""
    raise Exception("Test Error")


@app.route("/start_clean_up_manually")
@handle_internal_exception
def start_clean_up_manually():
    """Trigger clean up not by appscheduler."""
    thread = multiprocessing.Process(
        target=clean_up
    )
    thread.start()
    return "Putzen!"


if __name__ == "__main__":
    app.run()
