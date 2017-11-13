import numpy as np
import pandas as pd
import re

# date conversion imports
from dateutil.relativedelta import relativedelta as rd
import datetime
from datetime import timedelta
import time

def get_tree_names(tree_string, verbose=False):
    """Extracts the names of nodes from a tree string in newick format.
    Note that internal nodes are sometimes blank, or are sometimes numbers,
    and I have no idea why. Need a more structural way of identifying tips.

    Params
    ------
    tree_string: str; a tree string in newick format

    Returns
    -------
    d_nodes: dataframe of node names, and node types.
    """
    print("WARNING: HACKY AND EXPERIMENTAL.")

    # sometimes leaves the enclosing symbols for some reason
    pattern = re.compile(r'[(),][\s\S]*?(?=:|;)')

    names_ls = re.findall(pattern, tree_string)
    for i in range(len(names_ls)):
        for char in ["(", ")", ","]:
            names_ls[i] = names_ls[i].replace(char, "")

    if verbose:
        print("no. of nodes = %s" % len(names_ls))

    for i in range(len(names_ls)):
        nm = names_ls[i]
        if (len(nm) == 0) or _is_number(nm):
            names_ls[i] = [nm, "node"]
        elif len(nm) != 0 or (_is_number(nm) == False):
            names_ls[i] = [nm, "tip"]
    d_nodes = pd.DataFrame(data=names_ls, columns=["name", "node_type"])

    return d_nodes


def _is_number(s):
    """Checks if a string s is actually a number.
    Used in get_tree_names
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


def prep_nexus_contents(tree_string):
    """Preps a nexus format around a tree string.
    """

    d_nodes = get_tree_names(tree_string, verbose=False)
    d_tips = d_nodes.loc[d_nodes["node_type"]=="tip"]

    nex_contents = ["#NEXUS\nbegin taxa;"]
    nex_contents.append("\tdimensions ntax=%s;" % d_tips.shape[0])
    nex_contents.append("\ttaxlabels")
    for index, row in d_tips.iterrows():
        row = "\t" + str(row["name"])
        nex_contents.append(row)
    nex_contents.append(";\nend;\n")

    nex_contents.append("begin trees;\n\ttree TREE1 = [&R]\t")
    nex_contents.append("\t"+tree_string)
    nex_contents.append("end;")

    return nex_contents
