import os
from pathlib import Path
import subprocess
ROOT_DIR = os.path.dirname(str(Path(__file__)))
CONFIG_DIR = os.path.join(ROOT_DIR, 'config')

cwd = os.getcwd()
os.chdir(ROOT_DIR)
git_hash = subprocess.check_output(['git', 'describe', '--always']).strip()
git_date = subprocess.check_output(['git', 'log', '-1', '--format=%cd']).strip()
__version__ = git_hash.decode('utf8') + ' ' + git_date.decode('utf8')
os.chdir(cwd)
