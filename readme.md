# Welcome to SuPyFit!

SuPyFit (Python Superfit) is a software for the spectral classification of Supernovae 

## Install 

To install the software the user should download all the files and unzip the `30A.zip` and `20A.zip` files into a folder called "binnings".

## To run the code for an individual object

To achieve this task the files needed are: 

- `Header_binnings.py`
- `error_routines.py`
- `SF_functions.py`
- `params.py`
- `run.py`


In the `params.py` file there are three paths that the user should change.

- The "save_bin_path" to which the binned files will be saved.
- The "save_results_path" to which the results (a csv file and pdf images of the plots) will be saved.
- The "path" which is the location of the "binnings" folder. 

In the `run.py` file the user should change the "original" path to be that of the object of interest.


## Main SuPyFit function 

In the `run.py` file we find the main function which looks like this:


`all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, 
    lam, resolution, n=2, plot=1, kind=kind, original=save_bin, path=path, save=save_results_path, show=show)`
    
    
The inputs of the function are updated in the `params.py` file and are as follow: 

- `redshift:` Can be an array or an individual number. These are the redshift values over which to optimize. 
- `templates_sn_trunc:`  truncated library of supernovae, aka: which SN types to look at when optimizing.
- `templates_gal_trunc:` truncated library of host galaxies, aka: which HG types to look at when optimizing.
- `lam:` lambda array over which to perform the fit. The default is from 3000 A to 10500 A. 
- `resolution:` resolution at which to bin and perform the fit. The default is 20 A. 
- `n:` this corresponds to the number of plots to show and save as a result. 
- `plot:` either 1 or 0, to either plot or not plot. 
- `kind:` corresponds to the type of error spectrum the user prefers, the options are 'SG':Savitsky Golay, 'linear': for obtaining the error of the spectrum 
by making linear fit every 10 points, and 'included': if the user wants to use the error that comes with the object itself. The default is 'SG'


The `templates_sn_trunc` and `templates_sn_trunc` are updated by changing the `temp_gal_tr` and `temp_sn_tr` lists on the `params.py` file, to what the user is
interested in seeing (default is full library).


The rest the inputs correspond to the paths mentioned above. 
    
## Results

The results are: an astropy table that is saved as a csv file (to the specified path) and the best fit plots saved as pdf files (to the specified path)


## To run

Once the parameters have been updated in the `params.py` file the user simply needs to run the script from the `run.py` file. 
