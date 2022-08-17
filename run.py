from supyfit.sf_class import *

supernova = Superfit()
supernova.superfit()

for i in range(Parameters.n):
    supernova.any_result(i)