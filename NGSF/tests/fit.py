import json
import numpy as np
import os
import pandas as pd
import subprocess
import sys


def test_fit():

    SUPERFIT_PATH = "./"
    SUPERFIT_PARAMETERS_JSON = "./parameters.json"
    SUPERFIT_DATA_PATH = "./NGSF/tests/data"
    filebase = "SN2021urb_2021-08-06_00-00-00_Keck1_LRIS_TNS"

    params = json.loads(open(SUPERFIT_PARAMETERS_JSON).read())
    params["object_to_fit"] = f"{SUPERFIT_DATA_PATH}/{filebase}.flm"
    params["show_plot"] = 0

    JSON_FILE = f"{SUPERFIT_DATA_PATH}/{filebase}.json"
    with open(JSON_FILE, "w") as f:
        json.dump(params, f)

    subprocess.call(
        f"python run.py {SUPERFIT_DATA_PATH}/{filebase}.json",
        shell=True,
        stdout=sys.stdout,
        stderr=subprocess.STDOUT,
    )

    results_path = os.path.join(SUPERFIT_PATH, f"{filebase}.csv")
    results = pd.read_csv(results_path)
    results.sort_values(by=["CHI2/dof"], inplace=True)

    assert all([np.isclose(z, 0.127) for z in results["Z"]])
    assert all([~np.isnan(chi2) for chi2 in results["CHI2/dof"]])
