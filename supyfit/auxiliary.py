import numpy as np




def select_templates(DATABASE, TYPES):

    
    '''
    
    
    Selects templates of a given type(s) from a template database
   
    Input: DATEBASE   list of templates
           TYPES      which types should be selected
  
    Output: array of templates of given type(s)
    
    
    '''     
    
    
    database_trunc = list([])
    
    for type in TYPES:
        database_trunc += list([x for x in DATABASE if type in x])
    
    return np.array(database_trunc)



def kill_header(file_name):


    lines = [] 
    
    file = open(file_name,'r')
    
    lines = file.readlines()

    lines = [i for i in lines if i]

    lines = [i for i in lines if i[0].isalpha() == False and i[0] != '#' and i[0] != '%']

    lines = [i for i in lines if i[0] != '\n']
    
    lines = [s.strip('\n') for s in lines] # remove empty lines
    
    lines = [s.replace('\n', '') for s in lines]  #replace with nothing
    
    
    columns = [] 
    
    for line in lines:
        ii = line.split()
        columns.append(ii)
        
    columns = np.array(columns)
    
    lam_floats  = [float(i) for i in columns[:,0]]
    flux_floats = [float(i) for i in columns[:,1]]

    spectrum = np.array([lam_floats, flux_floats]).T
    
    return spectrum


