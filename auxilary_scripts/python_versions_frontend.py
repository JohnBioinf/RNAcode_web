"""Print python versions."""

import sys

import filelock
import coolname
import Bio
import validate_email
import flask_limiter

import apscheduler

import ete3
import jinja2

print(f"\\texttt{{python (webserver)}} & version:{sys.version.split()[0]} \\\\")
python_concat = "".join(sys.version.split()[0].split("."))
python_url = f"https://www.python.org/downloads/release/python-{python_concat}"
print(f" & url:\\url{{{python_url}}} \\\\")
print(f"\\texttt{{filelock}} & version:{filelock.__version__} \\\\")
print(f"\\texttt{{coolname}} & version:{coolname.__version__} \\\\")
print(f"\\texttt{{validate_email}} & version:{1.3} \\\\")
print(f"\\texttt{{flask_limiter}} & version:{flask_limiter.__version__} \\\\")
print(f"\\texttt{{apscheduler}} & version:{apscheduler.__version__} \\\\")
print(f"\\texttt{{ete3}} & version:{ete3.__version__} \\\\")
print(f"\\texttt{{jinja2}} & version:{jinja2.__version__} \\\\")
print(f"\\texttt{{Bio}} & version:{Bio.__version__} \\\\")
print()
print(f"<dt>python (webserver)</dt>\n<dd>version:{sys.version.split()[0]}")
print(f'<a href="{python_url}">\n{python_url}</a><br>\n</dd>')
print('<div class="m-3">')
print(f"  <dt>filelock</dt>  \n<dd>version:{filelock.__version__}</dd>")
print(f"  <dt>coolname</dt>  \n<dd>version:{coolname.__version__}</dd>")
print(f"  <dt>validate_email</dt>  \n<dd>version:{1.3}</dd>")
print(f"  <dt>flask_limiter</dt>  \n<dd>version:{flask_limiter.__version__}</dd>")
print(f"  <dt>apscheduler</dt>  \n<dd>version:{apscheduler.__version__}</dd>")
print(f"  <dt>ete3</dt>  \n<dd>version:{ete3.__version__}</dd>")
print(f"  <dt>jinja2</dt>  \n<dd>version:{jinja2.__version__}</dd>")
print(f"  <dt>Bio</dt>  \n<dd>version:{Bio.__version__}</dd>")
print('<\div>')
