import glob
import numpy as np
from astropy.io import ascii
import os
import csv
import json
import NGSF_version
from NGSF.params import Parameters

try:
    configfile = os.environ["NGSFCONFIG"]
except KeyError:
    configfile = os.path.join(NGSF_version.CONFIG_DIR, "parameters.json")
with open(configfile) as config_file:
    ngsf_cfg = json.load(config_file)


def list_folders(path):

    folders = []
    dirs = glob.glob(os.path.join(path, "*"))
    for d in dirs:
        if os.path.isdir(d):
            folders.append(d)

    return folders


class Metadata(object):

    def __init__(self):

        parameters = Parameters(ngsf_cfg)

        # mjd_max_brightness = glob.glob('**/mjd**')[0]
        mjd_max_brightness = os.path.join(parameters.pkg_dir,
                                          'NGSF/mjd_of_maximum_brightness.csv')

        with open(mjd_max_brightness, mode='r') as inp:
            reader = csv.reader(inp)
            band_dictionary = {rows[0]: rows[2] for rows in reader}

        with open(mjd_max_brightness, mode='r') as inp:
            reader = csv.reader(inp)
            mjd_dictionary = {rows[0]: rows[1] for rows in reader}

        sne_folder = os.path.join(parameters.bank_dir,
                                  'original_resolution', 'sne')

        folders = [os.path.join(parameters.pkg_dir, sne_folder, x) for x in
                   parameters.temp_sn_tr]
        have_wiserep = []
        no_wiserep = []
        z_dic = {}
        path_dic = {}
        dictionary_all_trunc_objects = {}
        jd_dic = {}
        coord_dic = {}
        spec_file_dic = {}
        inst_dic = {}
        obs_date_dict = {}
        shorhand_dict = {}
        type_dic = {}
        subfolders = []
        short_path_dict = {}
        for folder in folders:
            subs = list_folders(folder)
            for sub in subs:
                subpath = sub
                idx = subpath.rfind('/')
                sub = subpath[(idx+1):]
                subfolders.append(subpath)
                idx2 = subpath[0:idx].rfind('/')
                sn_type = subpath[idx2+1:idx]
                type_dic[sub] = sn_type
                if os.path.exists(subpath+'/wiserep_spectra.csv'):
                    have_wiserep.append(subpath)
                    wise = ascii.read(subpath+'/wiserep_spectra.csv')
                    path_dic[sub] = subpath
                    z_dic[sub] = wise['Redshift'][0]
                    coord_dic[sub] = np.array(list(wise['Obj. RA',
                                                        'Obj. DEC'][0]))

                    jd_dic[sub] = np.array(wise['JD'][:])
                    obs_date_dict[sub] = np.array(wise['Obs-date'][:])
                    spec_file_dic[sub] = np.array(wise['Ascii file'][:])
                    inst_dic[sub] = np.array(wise['Instrument'][:])
                    for i, spec_file in enumerate(spec_file_dic[sub]):

                        if float(mjd_dictionary[sub]) == -1:

                            phase = 'u'

                        else:

                            phase = float(wise['JD'][i]) - (
                                float(mjd_dictionary[sub]) + 2400000.5)

                            phase = round(phase, 2)

                        if parameters.epoch_high == parameters.epoch_low:

                            band = band_dictionary[sub]

                            shorhand_dict[spec_file] = \
                                sn_type + '/' + sub + '/' + wise[
                                    'Instrument'][i]+' phase-band : ' + str(
                                    phase) + str(band)

                            short_path_dict[shorhand_dict[
                                spec_file]] = spec_file

                            dictionary_all_trunc_objects[spec_file] = \
                                os.path.join(parameters.pkg_dir, sne_folder,
                                             sn_type, sub, spec_file)

                        else:

                            if phase != 'u' and parameters.epoch_low <= phase \
                                    <= parameters.epoch_high:

                                band = band_dictionary[sub]

                                shorhand_dict[spec_file] = \
                                    sn_type + '/' + sub + '/' + \
                                    wise['Instrument'][i] + ' phase-band : ' + \
                                    str(phase) + str(band)

                                short_path_dict[shorhand_dict[
                                    spec_file]] = spec_file

                                dictionary_all_trunc_objects[spec_file] = \
                                    sne_folder + sn_type + '/' + sub + '/' + \
                                    spec_file

                else:
                    no_wiserep.append(subpath)

        self.shorhand_dict = shorhand_dict
        self.no_wiserep = no_wiserep
        self.dictionary_all_trunc_objects = dictionary_all_trunc_objects
