#!/usr/bin/python3
"""Administration script to delete all files related to a job."""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from JobSubmission import JobSubmission  # noqa: E402
from JobSubmission import JOB_SUBMISSION_PATH_TEMPLATE  # noqa: E402
from JobSubmission import SubprocessError  # noqa: E402


def main():
    """Delete jobs frontend and backend."""
    job_id_list = sys.argv[1:]

    job_submission_dir = os.path.dirname(JOB_SUBMISSION_PATH_TEMPLATE.format("none"))
    print(f"Should {len(job_id_list)} jobs be deleted?")
    print("\n".join(job_id_list))
    while True:
        delete = input("Delete jobs? (Y/n)")
        if delete in ["y", "Y", ""]:
            break
        if delete in ["n", "N"]:
            sys.exit()
        print(f"Do not understand {delete}")

    for job_id in job_id_list:
        with JobSubmission(job_id) as job_submission:
            try:
                job_submission.delete(recursive=True)
            except SubprocessError:
                continue

    print("Delete left over filelocks")
    for lock_file in os.listdir(job_submission_dir):
        if lock_file.split(".")[-1] != "lock":
            continue
        job_submission_file = os.path.join(job_submission_dir, lock_file.replace(".lock", ""))
        if os.path.isfile(job_submission_file):
            continue
        print(f"Lock file {lock_file} is left over")
        os.remove(os.path.join(job_submission_dir, lock_file))


if __name__ == "__main__":
    main()
