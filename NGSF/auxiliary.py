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



