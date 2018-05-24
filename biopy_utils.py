import numpy as np
import pandas as pd
import re
import itertools


# date conversion imports
from dateutil.relativedelta import relativedelta as rd
import datetime
from datetime import timedelta
import time

import xio
from Bio import Phylo


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


def tip_to_tip_distance(my_tree, tip1, tip2):
    """Computes the tip-to-mrca-to-tip distance between tip1 and tip2.

    PARAMS
    ------
    my_tree: biopython tree object.
    tip1, tip2: biopython Clade objects

    RETURNS
    -------
    tt_dist: float. tip-to-mrca-to-tip distance between tip1 and tip2.
    """

    mrca = my_tree.common_ancestor(tip1, tip2)
    t1_mrca_dist = my_tree.distance(tip1, mrca)
    t2_mrca_dist = my_tree.distance(tip2, mrca)

    tt_dist = t1_mrca_dist + t2_mrca_dist
    return tt_dist


def genetic_distance_matrix(my_tree, names_ls):
    """Returns an upper triangular similarity matrix of all pairwise branch
    distances between possible pairs of tipnames in names_ls.

    PARAMS
    ------
    my_tree: Bio.Phylo tree object
    names_ls: list of str. List of tipnames to compute genetic distance with,
    using branch lengths as a measure.

    RETURNS
    -------
    hm_data: np array of shape (len(names_ls), len(names_ls)), type float.
    """
    all_pairs = list(itertools.combinations((names_ls), 2))
    gd_ls = []
    n_seq = len(names_ls)
    for pair in all_pairs:
        x, y = pair
        gd = tip_to_tip_distance(my_tree, x, y)
        gd_ls.append(gd)

    # Construct heatmap data
    hm_data = np.zeros((n_seq, n_seq))
    idx = 0
    for i in range(n_seq):
        for j in range(i+1, n_seq):
            hm_data[i][j] = gd_ls[idx]
            idx += 1

    return hm_data


def get_clade_labels(my_tree, ref_names_ls, verbose=True):
    """Computes tip-to-tmrc-to-tip distances for each tipname in my_tree, to
    each reference name in ref_names_ls. my_tree must contain the reference
    names in ref_names_ls.

    PARAMS
    ------
    my_tree: a tree file, readable by Bio.Phylo.
    ref_names_ls: a list of reference tip names.
    verbose: verbosity parameter.

    RETURNS
    -------
    df: pandas dataframe with columns: tip_names, their distance to every given
    reference name, the nearest reference distance, and nearest reference label.

    TODO: implement thresholding by percentile.
    """
    t0 = time.time()

    # Get all tipnames
    names_ls = [x.name for x in my_tree.get_terminals()]

    # Check
    for ref_nm in ref_names_ls:
        if ref_nm not in names_ls:
            print("WARNING: %s not found in input tree!" % ref_nm)

    # Get the non-reference names
    non_ref_names_ls = list(set.difference(set(names_ls), set(ref_names_ls)))

    if verbose:
        print("No. of reference tip names = %s" % len(ref_names_ls))
        print("No. of non-reference tip names = %s" % len(non_ref_names_ls))

    # Compute all possible distances between all tips and all references
    contents = []
    for nm in non_ref_names_ls:
        line = [nm]
        for ref_nm in ref_names_ls:
            tt_dist = tip_to_tip_distance(my_tree, nm, ref_nm)
            line.append(tt_dist)

        contents.append(line)

    col_names = ["tip_name"]
    for ref_nm in ref_names_ls:
        col_names.append("dist_to_"+ref_nm)

    df = pd.DataFrame(data=contents,
                      columns=col_names)

    if verbose:
        print("Done in %.2fs" % (time.time() - t0))

    df_cols = list(df.columns)[1:]
    df["clade_label"] = df[df_cols].idxmin(axis=1)
    df["clade_label"] = df.apply(lambda row: str(row["clade_label"]).replace("dist_to_", ""), axis=1)
    df["min_dist"] = df.loc[:, df_cols].min(axis=1)

    return df


def get_consensus_seq(seq_ls, th=0.95):
    """Iterate over the columns (amino acids) of aa_array. Form a consensus seq which consists
    of the aas which are common to at least th% of the sequences (rows) in aa_array. 
    Otherwise, return '-'. 
    
    PARAMS
    ------
    seq_ls: list of str. Must all be of the same length.
    th: majority vote threshold, ranges between 0 and 1.
    
    RETURNS
    -------
    consensus_seq: consensus sequence. 
    """
    # All sequences must be of the same length!
    try: 
        row_len_ls = [len(x) for x in seq_ls]
        seq_len_arr, counts_arr = np.unique(np.array(row_len_ls), return_counts=True)
        if len(seq_len_arr) > 1:
            raise my_exception
    except my_exception as e:
        print("ERROR: not all sequences are of the same length!")
    else:
        # Convert the list of input sequences into an array of shape (n_seqs, seq_len)
        contents = [list(x) for x in seq_ls]
        contents = np.array(contents)

        idx_ls = []
        consensus_seq = ""
        for i in range(len(contents[0])):
            elem_arr, counts_arr = np.unique(contents[:, i], return_counts=True)
            counts_arr = counts_arr/len(d_temp)
            #print(elem_arr, counts_arr)
            new_aa = "-"
            new_idx = -1
            if check_value_in_list(th, counts_arr):
                new_idx = i
                new_aa = elem_arr[np.argmax(counts_arr)]

            idx_ls.append(new_idx)
            consensus_seq = consensus_seq + new_aa

    return consensus_seq


def _check_value_in_list(th, ls):
    """Given a list of numbers ls, check if there's at least 1 value in ls > th.
    So far, exclusively used by get_consensus_seq()
    """
    ls_bool = False
    for num in ls:
        if num > th:
            ls_bool = True
    return ls_bool


def thin_tree(fn_tree, interval_n, incl, outgroup_nm="", verbose=True):
    """Thins a rooted, ladderized tree by selecting/excluding every n-th sequence.
    Because of the way the thinning function works, duplicate names may mess with the usual
    application of this function in unforeseeable ways!
    
    PARAMS
    ------
    fn_tree: str; filename of tree file. Need not be rooted or aligned. 
    interval_n: int; thinning parameter, so as to include or exclude the interval_n-th sequence.
    incl: bool; whether to include (set to True) or exclude (set to False) the n-th sequence
    outgroup_nm: str; name of outgroup. If set to an empty string, midpoint root. 
    verbose: Bool; verbosity parameter.
    
    RETURNS
    -------
    names_ls2: thinned list of sequence names.
    
    NOTE: Biopython's ladderize function is different to Figtree's, so the positions of
    entire clades may be swapped, though the tree is still exactly the same; just displayed
    differently. 
    """
    
    # Read a tree to get an ordered list of names
    tree = Phylo.read(fn_tree, "newick")
    if outgroup_nm != "":
        tree.root_with_outgroup({'name': outgroup_nm})
    else:
        tree.root_at_midpoint()
    tree.ladderize(reverse=True) # Sort increasing

    # Get names ordered from top to bottom on a ladderized (increasing) tree
    names_ls = [x.name for x in tree.get_terminals()]
    
    # Check for duplicate names
    for nm in list(set(names_ls)):
        if names_ls.count(nm) > 1:
            print("WARNING: tipname %s duplicated in input tree!" % nm)

    # Filter for every n-th sequence
    incl = True
    names_ls2 = []
    for i in range(len(names_ls)):
        # Thin
        # include every n-th seq
        if incl:
            if i%interval_n == 0:
                names_ls2.append(names_ls[i])
        # exclude every n-th seq
        elif incl == False:
            if i%interval_n != 0:
                names_ls2.append(names_ls[i])
                
    # Append the outgroup, if one was specified
    if (outgroup_nm != "") and (outgroup_nm not in names_ls):
        names_ls2.append(outgroup_nm)
    
    if verbose:
        print("%s names reduced to %s" % (len(names_ls), len(names_ls2)))
        
    return names_ls2
