import sys
import os
import json

repo_path = os.path.dirname(os.path.realpath(__file__))

with open(f"{repo_path}/parameters_frontend_local.json", "r", encoding="UTF-8") as file_handle:
    server_parameters = json.load(file_handle)

sys.path.insert(0, repo_path)
sys.path.insert(0, f"{server_parameters['python_env_path']}/lib/python3.6/site-packages")

os.chdir(repo_path)

from app import app as application

# TODO Set before starting server
application.secret_key = 'some random name'
