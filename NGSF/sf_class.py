import os
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.ndimage import gaussian_filter1d


from NGSF.SF_functions import Alam, all_parameter_space, remove_telluric, mask_gal_lines
from NGSF.Header_Binnings import kill_header, kill_header_and_bin
from NGSF.error_routines import linear_error, savitzky_golay
from NGSF.get_metadata import Metadata
from NGSF.params import Parameters
import NGSF_version

try:
    configfile = os.environ["NGSFCONFIG"]
except KeyError:
    configfile = os.path.join(NGSF_version.CONFIG_DIR, "parameters.json")
with open(configfile) as config_file:
    ngsf_cfg = json.load(config_file)

parameters = Parameters(ngsf_cfg)


class Superfit:
    def __init__(self, spec_file):

        ngsf_cfg["object_to_fit"] = spec_file
        parameters.object_to_fit = spec_file
        self.original_path_name = parameters.object_to_fit
        self.name = os.path.basename(self.original_path_name)
        self.name_no_extension = self.name[: self.name.rfind(".")]

        self.spectrum = kill_header(self.original_path_name)

        obj_original_res = self.spectrum[:, 0][-1] - self.spectrum[:, 0][-2]

        how_many_bins = 0

        if obj_original_res > 10:
            how_many_bins = +(self.spectrum[:,
                              0][-1] - self.spectrum[:, 0][0]) / 30

        elif obj_original_res <= 10:
            how_many_bins = +(self.spectrum[:,
                              0][-1] - self.spectrum[:, 0][0]) / 10

        # Check if spectrum is short

        if how_many_bins < 35:
            raise TypeError("This spectrum is too short to fit!")

        self.lamda, self.flux = self.spectrum[:, 0], self.spectrum[:, 1]

        self.binned_name = os.path.join(parameters.save_results_path,
                                        self.name_no_extension + "_binned.txt")
        self.results_name = os.path.join(parameters.save_results_path,
                                         self.name_no_extension)

        self.results_path = self.results_name + ".csv"

        if parameters.mask_galaxy_lines == 1 and parameters.mask_telluric == 1:
            object_spec = mask_gal_lines(self.spectrum, parameters.redshift)
            object_spec = remove_telluric(object_spec)
        if parameters.mask_galaxy_lines == 1 and parameters.mask_telluric == 0:
            object_spec = mask_gal_lines(self.spectrum, parameters.redshift)
        if parameters.mask_galaxy_lines == 0 and parameters.mask_telluric == 1:
            object_spec = np.loadtxt(self.original_path_name)
            object_spec = remove_telluric(object_spec)
        if parameters.mask_galaxy_lines == 0 and parameters.mask_telluric == 0:
            object_spec = np.loadtxt(self.original_path_name)

        object_spec[:, 1] = object_spec[:, 1] / np.nanmedian(object_spec[:, 1])

        int_obj = interpolate.interp1d(
            object_spec[:, 0], object_spec[:, 1], bounds_error=False, fill_value="nan"
        )
        self.int_obj = int_obj(parameters.lam)

        self.metadata = Metadata()

        # Make json with the used parameters
        with open(
                os.path.join(parameters.save_results_path,
                             self.name_no_extension + "_pars_used.json"),
                "w") as fp:
            json.dump(ngsf_cfg, fp)

    def plot(self):

        plt.figure(figsize=(7 * np.sqrt(2), 7))
        plt.title(str(self.name), fontsize=17)
        plt.ylabel("Flux", fontsize=16)
        plt.xlabel("Lamda", fontsize=16)
        plt.plot(self.lamda, self.flux, "k")

    def mask_telluric(self):

        masked_spectrum = remove_telluric(self.spectrum)
        plt.figure(figsize=(7 * np.sqrt(2), 7))
        plt.title("Masked Telluric for " + str(self.name), fontsize=17)
        plt.ylabel("Flux", fontsize=16)
        plt.xlabel("Lamda", fontsize=16)
        plt.plot(masked_spectrum[:, 0], masked_spectrum[:, 1], "k")

    def mask_galaxy_lines(self):

        if parameters.use_exact_z != 1:
            raise Exception(
                "Make sure to pick an exact value for z in order to mask the host lines accordingly!"
            )
        else:
            Data = mask_gal_lines(self.spectrum, z_obj=parameters.redshift[0])
            plt.figure(figsize=(7 * np.sqrt(2), 7))
            plt.ylabel("Flux arbitrary", fontsize=14)
            plt.xlabel("Lamda", fontsize=14)
            plt.title(
                "Galaxy lines masked at z=" + str(parameters.redshift[0]),
                fontsize=15,
                fontweight="bold",
            )
            plt.plot(
                self.lamda, self.flux / np.median(self.flux), "r", label=str(self.name)
            )
            plt.plot(
                Data[:, 0],
                Data[:, 1] / np.median(Data[:, 1]),
                "b",
                label="Masked object",
            )
            plt.legend(framealpha=1, frameon=True, fontsize=12)
            # plt.savefig(str(self.name) + '_masked.pdf' )

    def sg_error(self):

        Data = self.spectrum
        error = savitzky_golay(self.spectrum)[:, 1]

        plt.figure(figsize=(7 * np.sqrt(2), 7))
        plt.ylabel("Flux arbitrary", fontsize=14)
        plt.xlabel("Lamda", fontsize=14)
        plt.title("Savitzky-Golay error estimation", fontsize=15, fontweight="bold")
        plt.fill_between(
            Data[:, 0],
            Data[:, 1] / np.median(Data[:, 1]) - error,
            Data[:, 1] / np.median(Data[:, 1]) + error,
            color="#FF4500",
            label="error",
        )
        plt.plot(
            Data[:, 0], Data[:, 1] / np.median(Data[:, 1]), "k", label=str(self.name)
        )
        plt.legend(framealpha=1, frameon=True, fontsize=12)
        # plt.savefig(str(self.name) + '_sg.pdf' )

    def linear_error(self):

        Data = self.spectrum
        error = linear_error(self.spectrum)[:, 1]

        plt.figure(figsize=(7 * np.sqrt(2), 7))
        plt.ylabel("Flux arbitrary", fontsize=14)
        plt.xlabel("Lamda", fontsize=14)
        plt.title("Linear error estimation", fontsize=15, fontweight="bold")
        plt.fill_between(
            Data[:, 0],
            (Data[:, 1] - error) / np.median(Data[:, 1]),
            (Data[:, 1] + error) / np.median(Data[:, 1]),
            color="#03AC13",
            label="error",
        )
        plt.plot(
            Data[:, 0], Data[:, 1] / np.median(Data[:, 1]), "k", label=str(self.name)
        )
        plt.legend(framealpha=1, frameon=True, fontsize=12)

    def mask_gal_lines_and_telluric(self):

        Data_masked = mask_gal_lines(self.name, z_obj=parameters.redshift)
        masked_spectrum = remove_telluric(Data_masked)

        plt.figure(figsize=(7 * np.sqrt(2), 7))
        plt.title("Masked Telluric and Galaxy lines", fontsize=17)
        plt.ylabel("Flux", fontsize=16)
        plt.xlabel("Lamda", fontsize=16)
        plt.plot(
            masked_spectrum[:, 0], masked_spectrum[:, 1], "k", label=str(self.name)
        )
        plt.legend(framealpha=1, frameon=True, fontsize=12)

    def superfit(self):

        try:
            print(
                "Running optimization for spectrum file: {0} with resolution = {1} Å".format(
                    self.name_no_extension, parameters.resolution
                )
            )

            kill_header_and_bin(
                self.original_path_name,
                parameters.resolution,
                save_bin=self.binned_name,
            )

            all_parameter_space(
                self.int_obj,
                parameters.redshift,
                parameters.extconstant,
                parameters.templates_sn_trunc,
                parameters.templates_gal_trunc,
                parameters.lam,
                parameters.resolution,
                parameters.iterations,
                kind=parameters.kind,
                original=self.binned_name,
                save=self.results_name,
                show=parameters.show,
                minimum_overlap=parameters.minimum_overlap,
            )

        except Exception:

            resolution = 30
            print("NGSF failed. Retrying for resolution = {0} Å".format(resolution))

            kill_header_and_bin(
                self.original_path_name, resolution, save_bin=self.binned_name
            )

            all_parameter_space(
                self.int_obj,
                parameters.redshift,
                parameters.extconstant,
                parameters.templates_sn_trunc,
                parameters.templates_gal_trunc,
                parameters.lam,
                resolution,
                parameters.iterations,
                kind=parameters.kind,
                original=self.binned_name,
                save=self.results_name,
                show=parameters.show,
                minimum_overlap=parameters.minimum_overlap,
            )

        self.results = pd.read_csv(self.results_path)

        result_number = 0

        if parameters.n > len(self.results):

            result_number = result_number + len(self.results)

        elif len(self.results) >= parameters.n:

            result_number = result_number + parameters.n

        for j in range(result_number):

            row = self.results.iloc[j]

            hg_name = row["GALAXY"]
            short_name = row["SN"]
            bb = row["CONST_SN"]
            dd = row["CONST_GAL"]
            z = row["Z"]
            extmag = row["A_v"]
            sn_cont = row["Frac(SN)"]

            # Get all names from the dictionary
            full_names = [str(x) for x in self.metadata.shorhand_dict.keys()]
            short_names = [str(x) for x in self.metadata.shorhand_dict.values()]

            # print(full_names)

            subtype = ""
            sn_best_fullname = ""
            for i in range(0, len(short_names)):
                if str(short_names[i]) == str(short_name):
                    sn_best_fullname = full_names[i]
                    sn_short_name = short_names[i]
                    idx = sn_short_name.rfind("/")
                    subtype = sn_short_name[:idx]

            int_obj = self.int_obj

            sn_name = os.path.join(parameters.bank_dir, "binnings", "10A",
                                   "sne",  subtype,  sn_best_fullname)
            hg_name = os.path.join(parameters.bank_dir, "binnings", "10A",
                                   "gal", hg_name)

            # print(sn_name)

            nova = kill_header(sn_name)
            nova[:, 1] = nova[:, 1] / np.nanmedian(nova[:, 1])

            host = np.loadtxt(hg_name)
            host[:, 1] = host[:, 1] / np.nanmedian(host[:, 1])

            # Interpolate supernova and host galaxy
            # redshifted_nova   =  nova[:,0]*(z+1)
            # extinct_nova      =  nova[:,1]*10**(-0.4*extmag * Alam(nova[:,0]))/(1+z)

            # reshifted_host    =  host[:,0]*(z+1)
            # reshifted_hostf   =  host[:,1]/(z+1)

            redshifted_nova = nova[:, 0] * (z + 1)
            extinct_nova = (
                nova[:, 1] * 10 ** (-0.4 * extmag * Alam(nova[:, 0])) / (z + 1)
            )

            reshifted_host = host[:, 0] * (z + 1)
            reshifted_hostf = host[:, 1] / (z + 1)

            nova_int = interpolate.interp1d(
                redshifted_nova, extinct_nova, bounds_error=False, fill_value="nan"
            )
            host_int = interpolate.interp1d(
                reshifted_host, reshifted_hostf, bounds_error=False, fill_value="nan"
            )
            host_nova = bb * nova_int(parameters.lam) + dd * host_int(parameters.lam)

            sn_type = short_name[: short_name.find("/")]
            hg_name = hg_name[hg_name.rfind("/") + 1 :]
            subclass = short_name[short_name.find("/") + 1 : short_name.rfind("/")]
            phase = str(short_name[short_name.rfind(":") + 1 : -1])

            plt.figure(figsize=(8 * np.sqrt(2), 8))
            plt.plot(parameters.lam, int_obj, "r", label="Input object: " + self.name)
            plt.plot(
                parameters.lam,
                host_nova,
                "g",
                label="SN: "
                + sn_type
                + " - "
                + subclass
                + " - Phase: "
                + phase
                + "\nHost: "
                + str(hg_name)
                + "\nSN contrib: {0: .1f}%".format(100 * sn_cont),
            )
            plt.legend(framealpha=1, frameon=True, fontsize=12)
            plt.ylabel("Flux arbitrary", fontsize=14)
            plt.xlabel("Lamda", fontsize=14)
            plt.title("Best fit for z = " + str(z), fontsize=15, fontweight="bold")

            if parameters.show_plot_png:
                plt.savefig(self.results_name + "_" + str(j) + ".png")
            else:
                plt.savefig(self.results_name + "_" + str(j) + ".pdf")

            if parameters.show == 1:
                plt.show()

    def results(self):

        if os.path.isfile(self.results_path):
            results = self.results
        else:
            raise Exception("Do the superfit! <( @_@" ")> ")

        return results

    def any_result(self, j):

        row = self.results.iloc[j]

        hg_name = row["GALAXY"]
        short_name = row["SN"]
        bb = row["CONST_SN"]
        dd = row["CONST_GAL"]
        z = row["Z"]
        extmag = row["A_v"]
        sn_cont = row["Frac(SN)"]

        # Get all names from the dictionary
        full_names = [str(x) for x in self.metadata.shorhand_dict.keys()]
        short_names = [str(x) for x in self.metadata.shorhand_dict.values()]

        for i in range(0, len(short_names)):
            if str(short_names[i]) == str(short_name):
                sn_best_fullname = full_names[i]
                sn_short_name = short_names[i]
                idx = sn_short_name.rfind("/")
                subtype = sn_short_name[:idx]

        int_obj = self.int_obj

        sn_name = "bank/binnings/10A/sne/" + subtype + "/" + sn_best_fullname
        hg_name = "bank/binnings/10A/gal/" + hg_name

        nova = kill_header(sn_name)
        nova[:, 1] = nova[:, 1] / np.nanmedian(nova[:, 1])

        host = np.loadtxt(hg_name)
        host[:, 1] = host[:, 1] / np.nanmedian(host[:, 1])

        # Interpolate supernova and host galaxy
        redshifted_nova = nova[:, 0] * (z + 1)
        extinct_nova = nova[:, 1] * 10 ** (-0.4 * extmag * Alam(nova[:, 0])) / (z + 1)

        reshifted_host = host[:, 0] * (z + 1)
        reshifted_hostf = host[:, 1] / (z + 1)

        nova_int = interpolate.interp1d(
            redshifted_nova, extinct_nova, bounds_error=False, fill_value="nan"
        )
        host_int = interpolate.interp1d(
            reshifted_host, reshifted_hostf, bounds_error=False, fill_value="nan"
        )
        host_nova = bb * nova_int(parameters.lam) + dd * host_int(parameters.lam)

        sn_type = short_name[: short_name.find("/")]
        hg_name = hg_name[hg_name.rfind("/") + 1 :]
        subclass = short_name[short_name.find("/") + 1 : short_name.rfind("/")]
        phase = str(short_name[short_name.rfind(":") + 1 : -1])
        plt.figure(figsize=(8 * np.sqrt(2), 8))
        plt.plot(parameters.lam, int_obj, "r", label="Input object: " + self.name)
        plt.plot(
            parameters.lam,
            host_nova,
            "g",
            label="SN: "
            + sn_type
            + " - "
            + subclass
            + " - Phase: "
            + phase
            + "\nHost: "
            + str(hg_name)
            + "\nSN contrib: {0: .1f}%".format(100 * sn_cont),
        )
        plt.legend(framealpha=1, frameon=True, fontsize=12)
        plt.ylabel("Flux arbitrary", fontsize=14)
        plt.xlabel("Lamda", fontsize=14)
        plt.title("Best fit for z = " + str(z), fontsize=15, fontweight="bold")

        if parameters.show_plot_png:
            plt.savefig(self.results_name + "_" + str(j) + ".png")
        else:
            plt.savefig(self.results_name + "_" + str(j) + ".pdf")

        if parameters.show == 1:
            plt.show()

    def convolution(self):

        obj_res = 100

        obj_med = np.median(self.lamda)
        width = obj_med / obj_res
        sig = width / (2 * np.sqrt(2 * np.log(2)))

        filtered = gaussian_filter1d(self.flux, sig)
        plt.plot(self.lamda, filtered)
