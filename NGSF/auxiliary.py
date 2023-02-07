import numpy as np


def select_templates(database, types):

    """


    Selects templates of a given type(s) from a template database

    Input: database   list of templates
           types      which types should be selected

    Output: array of templates of given type(s)


    """

    database_trunc = list([])

    for tp in types:
        database_trunc += list([x for x in database if tp in x])

    return np.array(database_trunc)
