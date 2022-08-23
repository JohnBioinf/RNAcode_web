"""Print python versions."""

import sys

import Bio
import ete3

print(f"\\texttt{{python (cluster)}} & version:{sys.version.split()[0]} \\\\")
python_concat = "".join(sys.version.split()[0].split("."))
python_url = f"https://www.python.org/downloads/release/python-{python_concat}"
print(f" & url:\\url{{{python_url}}} \\\\")
print(f"\\texttt{{ete3}} & version:{ete3.__version__} \\\\")
print(f"\\texttt{{Bio}} & version:{Bio.__version__} \\\\")
print()
print(f"<dt>python (cluster)</dt>\n<dd>version:{sys.version.split()[0]}")
print(f'<a href="{python_url}">\n{python_url}</a><br>\n</dd>')
print('<div class="m-3">')
print(f"  <dt>ete3</dt>\n  <dd>version:{ete3.__version__}</dd>")
print(f"  <dt>Bio</dt>\n  <dd>version:{Bio.__version__}</dd>")
print('</div>')
