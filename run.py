import json 
from sf_class import *

with open("parameters.json", "r") as read_file:
    data = json.load(read_file)

this_supernova = superfit_class(object_to_fit)
this_supernova.superfit()


