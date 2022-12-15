#!/usr/bin/env python

from setuptools import setup
import subprocess
import sys
import os

# check that python version is 3.8 or above
python_version = sys.version_info
if python_version < (3, 7):
    sys.exit("Python < 3.7 is not supported, aborting setup")
print("Confirmed Python version {}.{}.{} >= 3.7.0".format(*python_version[:3]))


def write_version_file(version):
    """Writes a file with version information to be used at run time

    Parameters
    ----------
    version: str
        A string containing the current version information

    Returns
    -------
    version_file: str
        A path to the version file

    """
    try:
        git_log = subprocess.check_output(
            ["git", "log", "-1", "--pretty=%h %ai"]
        ).decode("utf-8")
        git_diff = (
            subprocess.check_output(["git", "diff", "."])
            + subprocess.check_output(["git", "diff", "--cached", "."])
        ).decode("utf-8")
        if git_diff == "":
            git_status = "(CLEAN) " + git_log
        else:
            git_status = "(UNCLEAN) " + git_log
    except Exception as e:
        print("Unable to obtain git version information, exception: {}".format(e))
        git_status = ""

    version_file = ".version"
    if os.path.isfile(version_file) is False:
        with open("NGSF/" + version_file, "w+") as f:
            f.write("{}: {}".format(version, git_status))

    return version_file


def get_long_description():
    """Finds the README and reads in the description"""
    here = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(here, "README.md")) as f:
        long_description = f.read()
    return long_description


def get_requirements():
    with open("requirements.txt", "r") as ff:
        requirements = ff.readlines()
    return requirements


# get version info from __init__.py
def readfile(filename):
    with open(filename) as fp:
        filecontents = fp.read()
    return filecontents


requirements = [
    "astropy",
    "extinction",
    "matplotlib",
    "numpy",
    "pandas",
    "PyAstronomy",
    "scipy",
]

VERSION = "0.0.1"
version_file = write_version_file(VERSION)
long_description = get_long_description()

setup(
    name="NGSF",
    description="Next Generation SuperFit in Python",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/samanthagoldwasser25/NGSF",
    author="Samantha Goldwasser, Ofer Yaron, Avner Sass, Ido Irani, Avishay Gal-Yam, D. Andrew Howell",
    author_email="ofer.yaron@weizmann.ac.il ",
    license="MIT",
    version=VERSION,
    packages=[
        "NGSF",
    ],
    package_dir={"NGSF": "NGSF"},
    package_data={
        "NGSF": [version_file],
    },
    python_requires=">=3.7",
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
