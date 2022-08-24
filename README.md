
Table of Contents
=================

* [RNAcode Web](#rnacode-web)
   * [Setup](#setup)
      * [Variables in parameter_frontend:](#variables-in-parameter_frontend)
      * [Variables in parameter_backend:](#variables-in-parameter_backend)
   * [Logging](#logging)
   * [Communication Frontend Backend](#communication-frontend-backend)
   * [Front end structure](#front-end-structure)
      * [Job types](#job-types)
      * [Clean up](#clean-up)
   * [Backend](#backend)


# RNAcode Web

This repository contains all software to setup and run the web service RNAcode
web.

The web service provides a convenient way to call the coding region prediction
software RNAcode. RNAcode performs a statistical analysis on a nucleotide
alignment. The web service provides RNAcode with a nucleotide alignment based on
an input sequence and a few information to build the alignment.

The web service runs on two server. 1. The frontend web server which handles
the user input and output. 2. The backend server which performs the
computation.

The frontend server build on the flask framework.
The backend server interacts with the workload manager SLURM.

## Setup

The web service was build for a Fedora 30 (Workstation Edition) as backend and
CentOS Stream 8 as a frontend server. The frontend server that was runs the
service currently is a apache server and uses `mod_wsgi` to manage the flask
framework. For this the file `apache.conf` and `app.wsgi` are of needed. The
`appache.conf` shows the lines necessary for the web server, only the absolute
paths and the admin email must be set.

On both systems the git repository must be cloned. Then the json files
`parameters_frontend.json` and `parameters_backend.json` must be copied to
`parameters_frontend_local.json` respectively `parameters_backend_local.json`.
In both files the variables must be set accordingly. Next both scripts
`system_setup/build_env_frontend.sh` and `system_setup/build_env_backend.sh`
must be executed on the appropriate machines.

### Variables in parameter_frontend:

 * ip
   * IP address of the web server
 * user_dev
   * User name in development
 * user_production
   * User name in production
 * python_env_path
   * Path to python env
 * user_bin_path
   * Path to binaries
 * production_dir
   * Directory where the apache server will store logs, result, etc
 * smtp_server
   * address of smtp email server
 * rnacode_email_address
   * Email address from which the emails will be send
 * admin_email_address
   * Email address to which error messages will be send
 * port
   * Port for development/ set "None" in production
 * uni_net_space
   * Addressee space internally. Unlimited bulk jobs
 * R_bin
   * R binary path

### Variables in parameter_backend:

 * work_dir
   * Directory where all job submission data will be put
 * blast_db
   * Directory where current blast db can be found
 * python_env_backend_path
   * Path to python env of the backend manager
 * python_env_seqSel_path
   * Path to python env of the backend computation
 * RNAcode_web_repo_backend
   * Path to repository
 * blast_bin
   * Path to blast binaries
 * clustal_bin
   * Path to clustal binaries
 * RNAcode_bin
   * Path to RNAcode binaries
 * user
   * User name
 * machine_name
   * Name of the backend machine
 * entrez_bin
   * Path to entrez binaries
 * PATH
   * Default value of $PATH


Further the ssh keys from the frontend machine must be exchanged with the back
end. So that the command `ssh user@frontend_machine` functions without manually
typing the password.

On the backend machine a local blast db must be setup. To update or create a
new instance of the current blast db the script
`system_setup/update_blast_db.sh` should be used. Updating the blast db most be
done manually.

Both backend and frontend use the NCBI taxonomy data base. This is done with
the python module ete3. Initializing the module will look for a local version
of the database. In case it can not be found the database will be downloaded,
this might result in some downtime on the first start. Updates of the data base
are managed by the frontend on a regular base (time interval is set by the
variable `update_interval_ncbi`). The frontend will also send an update
request to the backend machine. Further the above mentioned script
`system_setup/update_blastdb.sh` will also updated the taxonomy data base.

For the deployment an apache server, which runs the module 'mod_wsgi', was used.
The config files are `RNAcode_web.wsgi` and `apache.conf`.

## Example Results

On the current startpage and the input page a few examples are given for
succesfull RNAcode Web jobs. The submission for these jobs was done as a bulk
submission using the file `static/jobs_publication.json`. If the server is setup
freshly the examples will not be available either submit new with the file (yet
this might change the outcome of the results depending how much the blast data
bases changed) or simply coppy all files from the <a
href="http://www.bioinf.uni-leipzig.de/Publications/SUPPLEMENTS/22-004/results_publication/">supplements</a>
of the publication to the workdir of the backend and use the job submision files
found in `./job_submissions`.

## Logging

Both backend and frontend will produce log files. The backend machine will use
the directory `./log`. The frontend machine will use in production mode the
given directory in `./parameters_frontend_local.json` and in development
`./logs`. For the frontend sever the verbose depth can be set with the variable
VERBOSE_LEVEL in `JobSubmission.py` and `app.py`.

## Communication Frontend Backend

The communication from frontend to backend is done via `ssh` commands. All
commands can be found in the folder `./ssh-scp-api/`. All communication between
frontend and backend is organised in the script `JobSubmission.py` the main
script of the web server `app.py` always use this module for communication.

The communication from backend to frontend is done with a url request. Every
time a job ends on the backend machine a request to the result page for the job
is made. This triggers the frontend to load the status of the job.

## Front end structure

The main logic for the flaks server is in the two scripts `app.py`, which
manages the server, handles request and renders the html sites, and
`JobSubmission.py`, where the logic for each job submission is described. The
main idea is that a user creates a job submission either by using the route
`single_submission`, filling out a form, or `bulk_submission`, providing a
`.json` file which has all parameters for the job. The web server will create
the class `JobSubmission` from the input. This is done by creating a file
`job_submissions/$job_id.json`. If the input contains error the user is informed
and in case of a bulk submission the jobs are deleted. In case of submission via
form, the job is kept and updated if a new form is submitted with corrected
values. If all values satisfy the criteria the user can submit the job. When the
job finished on the backend server and the frontend server is notified, the
results files will be downloaded to the frontend server and the user is notified
via email.

As the information for a job submission are stored in a file, simultaneous write
access can create errors or undesired behaviour. Thus the `JobSubmission` class
can be initialized with write access or without. The syntax `job =
JobSubmission(job_id)` will only allow read access to the job submission. The
syntax `with JobSubmission(job_id) as job:` will allow write access.

### Job types

There are three main job types `orphan`, `parent`, `child`. An `orphan` has an
input length below the threshold `MAX_LEN_NO_SPLIT` and will simply start one
job at the backend. For longer jobs a `parent` job will spawn several `child`
jobs with the same parameters except a shorter version of its input sequence.
If the input sequence for a `parent` job is below the threshold of
`LENGTH_NO_CUSTOM` the job will simply spawn and submit the children jobs. If
the sequence is longer the parent will first submit the construction of a custom
data base and after this is finished submit the children jobs.

Thus there are three types of pipelines can be started:

 * Using the NCBI data base.
   * Possible job types: `orphan` or `child`
   * Backend script: `RNAcodeWebCore_NCBI_DB.py`
 * Constructing a custom data base.
   * Possible job types: `parent`
   * Backend script: `RNAcodeWebCore_build_DB.py`
 * Using a custom data base.
   * Possible job types: `child`
   * Backend script: `RNAcodeWebCore_custom_DB.py`

### Clean up

The server will uses the module `BackgroundScheduler` to schedule a regular
cleaning of all job submission. This will delete all non submitted jobs and all
jobs which are older than `JOB_STORAGE_TIME`. Jobs from the publication will be
excluded.

## Backend

The backend structure is straight forward. Every time a job gets submitted one
of the scripts `RNAcodeWebCore_NCBI_DB.py`, `RNAcodeWebCore_build_DB.py` or
`RNAcodeWebCore_custom_DB.py` is started, see "Job types" for explanation.
Which will build a job directory and spawn SLURM process, which will do the
computation. The logs for each step will be written to `./logs/$job_id.err` and
`./logs/$job_id.out`.
