import numpy as np
import pandas as pd
import re
import subprocess

from Bio import SeqIO
from Bio import Phylo
from io import StringIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna

import time

"""treesubsub helper functions. Has 3 of IO classics: read_fasta(), prep_fasta_contents()
and write_fasta().
Terminology: here, following graph theory convention, all leaves and internal
nodes in a tree are called 'nodes'.
"""
def read_fasta(fn, delimiter="|", preview=0):
    """Reads a fasta into a list of lists, where each list is one record.

    Params
    ------
    fn: string; filename
    verbose: boolean; verbosity.
    delimiter: str; delimiting char in the fasta header.
    preview: int; no. of rows to preview. DNA sequence previews are truncated.

    Returns
    -------
    contents: a list of lists.
    """
    print_limit = 0
    contents = []
    for record in SeqIO.parse(fn, "fasta"):
        header = record.description
        seq = str(record.seq)
        if delimiter != "":
            row = header.split(delimiter) + [seq]
        elif delimiter == "":
            row = [header]+[seq]
        contents.append(row)

        if preview > 0:
            if print_limit < preview:
                preview_seq = row[-1][:10] + "..."
                print_row = row[:-1] + [preview_seq]
                print(print_row)
                print_limit +=1

    # strip whitespaces, which GISAID tends to add
    for i in range(len(contents)):
        for j in range(len(contents[i])):
            contents[i][j] = contents[i][j].strip()

    return contents


def prep_fasta_contents(d0_in, header_cols, seq_col='seq', delimiter="|", preview=0):
    """Reads a dataframe into a list of lists, in fasta format.

    Params
    ------
    d0_in: input dataframe
    header_cols: list of str. Whatever columns from d0_in to put in the fasta header
    seq_col: str. The name of the column which contains actual sequence data.
        default = 'seq'.
    preview: int. No. of records to preview after the function is finished.

    Returns
    -------
    contents: a list of lists
    """

    contents = []
    for index, row in d0_in.iterrows():
        header_ls = row[header_cols]
        header_ls = [str(item) for item in header_ls]
        fasta_header = ">" + delimiter.join(header_ls)
        seq = str(row[seq_col])
        contents.append(fasta_header)
        contents.append(seq)

    # display contents
    if preview > 0:
        for i in range(preview*2):
            if i%2 == 0:# for even rows
                print(contents[i])
            else:
                print(contents[i][:20]+ " ... " + contents[i][-20:])

    return contents


def write_fasta(fn_out, d_in, header_cols, seq_col="seq", delimiter="|", preview=0):
    """Wrapper for prep_fasta_contents(), and does the additional step of
    writing out to a file.

    Params
    ------
    d_in: input dataframe
    header_cols: list of str. Whatever columns from d_in to put in the fasta header
    seq_col: str. The name of the column which contains actual sequence data.
        default = 'seq'.
    preview: int. No. of records to preview after the function is finished.
    fn_out: name of file to write out to.

    Returns
    -------
    Writes out a fasta file.
    """
    contents = prep_fasta_contents(d_in, header_cols,
                                   seq_col=seq_col,
                                   delimiter=delimiter,
                                   preview=preview)
    with open(fn_out, "w") as f:
        for line in contents:
            f.write("%s\n" % line)
    print("Wrote to file %s" % fn_out)


def main_preprocess(fn_in):
    """Prepare raxml input files and baseml input files.

    PARAMS
    ------
    fn_in: str; input file of sequence data.

    OUTPUT:
    -------
    Returns nothing, but writes out the required raxml and baseml input files.
    """
    contents = read_fasta(fn_in, delimiter="")
    d0 = pd.DataFrame(data=contents, columns=["header", "seq"])

    seq_num_range = np.arange(len(d0))+1
    d0["seq_index"] = ["seq_"+str(x) for x in seq_num_range]

    # Write out baseml input alignment file
    # It's in some weird format
    # First row is num of seqs, then seq length, then "GC", apparently
    contents = [str(len(d0))+" "+str(len(d0.iloc[0].seq))+" "+"GC"]
    for index, row in d0.iterrows():
        contents.append(str(row["seq_index"])+"   ")
        contents.append(str(row["seq"]).upper())

    with open("aln_baseml_in.phy", "w") as f:
        for line in contents:
            f.write("%s\n" % line)
    print("Wrote to file aln_baseml_in.phy")

    # Write out RAXML input file
    # This is just a normal .fas, with seq_index as headers
    write_fasta("aln_raxml_input.fas", d0, ["seq_index"])

    # write out original v_name to seq_index dictionary
    d0[["header", "seq_index"]].to_csv("aln_names.csv", index=True)
    print("Wrote to file aln_names.csv")


def main_calls():
    """External subprocess call to RAXML for tree computation, and PAML baseml
    for ancestral reconstruction.
    """
    # Run raxml, delete all temp files except for best tree
    cmd = "raxml -m GTRGAMMA -T 2 -# 10 -p 12345 -s aln_raxml_input.fas -n rax"
    subprocess.run(cmd, shell=True)
    cmd = "cp RAxML_bestTree.rax best_tree"
    subprocess.run(cmd, shell=True)

    # baseml call
    prep_baseml_ctl("aln_baseml_in.phy", "best_tree_rooted")
    cmd = "baseml baseml.ctl"
    subprocess.run(cmd, shell=True)


def main_parse_rst(rst_path):
    """Public function call that procs _rst_get_tree() and
    _rst_get_ancestral_sequences(). WIP, because I'm still trying to pick a
    standard output.

    PARAMS
    ------
    rst_path: str; path to rst file.

    RETURNS
    -------
    WIP
    """
    print("WIP")


def _rst_get_tree(rst_path):
    """The `rst` file output by `baseml` has a newick string with node labels '
    but no branch lengths, and another newick string with branch lengths but
    no node labels. Combine them both, and write out a file `tree0.tre`, that
    has internal node labels and branch lengths.

    The node labels are called `node#`, or `REALNAME`, by the actual treesub
    program.

    PARAMS
    ------
    rst_path: str; path to rst file.

    RETURNS
    -------
    tree0: rst tree.
    """
    with open(rst_path) as f:
        rst_contents = f.readlines()
    rst_contents = [x.strip() for x in rst_contents]

    for i in range(len(rst_contents)):
        # tree w/o nodes but with br lengths
        if "Ancestral reconstruction by BASEML." in rst_contents[i]:
            t_string0 = rst_contents[i+2].replace(" ", "")
            tree0 = Phylo.read(StringIO(t_string0), "newick")

        # tree w nodes (node#) but no branch lengths
        if "Rod Page's TreeView" in rst_contents[i]:
            t_string1 = rst_contents[i+1].replace(" ", "")
            tree1 = Phylo.read(StringIO(t_string1), "newick")

    # Traverse both trees, and populate tree0 nodes with names from tree1
    # default = preorder traversal
    inodes_ls0 = tree0.get_nonterminals()
    inodes_ls1 = tree1.get_nonterminals()

    for i in range(len(tree0.get_nonterminals())):
        cl = inodes_ls0[i]
        cl.name = str(inodes_ls1[i].confidence)

    return tree0


def _rst_get_ancestral_sequences(rst_path):
    """Parses the 'rst' file for the internal node sequences, by searching for
    'List of extant and reconstructed sequences'.

    PARAMS
    ------
    rst: str; path to rst file.

    RETURNS
    -------
    ar_contents: list of rows of str; nt sequences of tips and internal nodes.
    Each row is in the format ['seq_index', 'nt_sequence']
    """
    with open(rst_path) as f:
        contents = f.readlines()
    contents = [x.strip() for x in contents]
    rst_contents = []
    for line in contents:
        if line != "":
            rst_contents.append(line)

    for i in range(len(rst_contents)):
        if rst_contents[i] == "List of extant and reconstructed sequences":
            idx_start = i+2

    ar_contents = []
    line = rst_contents[idx_start]
    while "Overall accuracy of the" not in line:
        line = rst_contents[idx_start]
        ar_contents.append(line)
        idx_start +=1
    ar_contents = ar_contents[:-1]

    # Some wrangling to get rid of whitespaces
    # Use '*' as a placeholder, to keep string lengths constant
    ar_contents2 = []
    for i in range(len(ar_contents)):
        line = ar_contents[i]
        for j in range(len(line)):
            if j < len(line):
                if (line[j] == " ") and (line[j+1] != " "):
                    line = str_replace_in_place(line, j, "*")
        ar_contents2.append(line)

    # Get rid of the placeholder
    for i in range(len(ar_contents2)):
        #print(ar_contents2[i])
        ar_contents2[i] = ar_contents2[i].replace("*", "")

    # Override ar_contents
    ar_contents = []
    for row in ar_contents2:
        new_line = []
        for elem in row.split(" "):
            if len(elem) > 0:
                new_line.append(elem)
        ar_contents.append(new_line)

    return ar_contents


def biopy_translate_nuc(seq_nuc):
    """Translate a nucleotide sequence. Unknown nucleotide characters will be
    translated to a codon X. First turns all gaps '-' into 'N'.

    PARAMS
    ------
    seq_nuc: str; a nucleotide sequence.

    RETURNS
    -------
    seq_aa: str; an amino acid sequence.

    EXAMPLES
    biopy_translate_nuc("AGG")
    >> 'R'
    biopy_translate_nuc("AGX")
    >> 'X'
    biopy_translate_nuc("AGGT")
    >> Deprecation warning: this will become an error in future!
    >> 'R'
    """
    seq_nuc = seq_nuc.replace("-", "N")

    my_seq = Seq(seq_nuc, generic_rna)
    seq_aa = str(my_seq.translate())
    return seq_aa


def str_replace_in_place(my_string, idx, new_char):
    """Replace the character at my_string[idx] with new_char.
    Necessary because Python doesn't allow in-place modifications for strings.

    PARAMS
    ------
    my_string: string to modify
    idx: int; index of character to modify
    new_char: string; new character to replace the old one with

    RETURNS
    -------
    new_string: new string.
    """
    str_ls = list(my_string)
    str_ls[idx] = new_char
    new_string = "".join(str_ls)
    return new_string


def get_aa_changes(target, query, output_format="nt", output_delimiter="-"):
    """Find the synonymous and nonsynonymous mutations between two nt sequences,
    target and query.
    This should supercede get_aa_diff() at some point.

    PARAMS
    ------
    target, query: str; nucleotide sequences.
    output_format: 'nt' returns each entry in the output as an amino acid triple,
    e.g. ['AAT-5-AGT']
    'aa' returns each entry in the output as the translated amino acid, e.g.
    ['N-5-S']
    output_delimiter: str. Output delimiter character.

    RETURNS
    -------
    aa_diff_ls: list of aa differences, expressed as "target" to "query", and the index,
    delimited by "-".
    Note that the index returned is the amino acid index.

    Note: if there happens to be nt changes that are right beside each other,
    aa_diff_ls will record these separately.

    EXAMPLES
    >> seq_t = "AGGTGCCGTACT"
    >> seq_q = "AGGTGTCGTAGT"
    >> get_aa_changes(seq_t, seq_q)
    ['TGC-2-TGT', 'ACT-4-AGT']
    >> get_aa_changes(seq_t, seq_q, output_format="aa")
    ['C-2-C', 'T-4-S']
    """
    # Should this break the function?
    if len(target) != len(query):
        print("WARNING: length of input strings are different!")

    # this should definitely break the function
    if len(target) % 3 != 0:
        print("ERROR: length of target string is not divisible by 3!")
    if len(query) % 3 != 0:
        print("ERROR: length of query string is not divisible by 3!")

    aa_diff_ls = []
    aa_idx = 0
    for i in range(len(target)):
        aa_idx = i//3 + 1 # add one to convert to 1-base
        if target[i] != query[i]:
            codon_coords = get_codon_window_coords(i)
            target_codon = target[codon_coords[0]:codon_coords[1]]
            query_codon = query[codon_coords[0]:codon_coords[1]]

            if output_format == "nt":
                aa_diff_row = output_delimiter.join([target_codon,
                                                     str(aa_idx),
                                                     query_codon])
            elif output_format == "aa":
                aa_diff_row = output_delimiter.join([biopy_translate_nuc(target_codon),
                                                     str(aa_idx),
                                                     biopy_translate_nuc(query_codon)])
            aa_diff_ls.append(aa_diff_row)

    return aa_diff_ls


def get_codon_window_coords(nt_idx):
    """Given a nucleotide index, compute the starting and ending nt coordinates
    of the codon which it belongs to.

    There's an off-by-one offset built into this function because of the way
    string slicing works in python: "ABCDEF"[3:6] returns "DEF".
    """
    start_idx = 0
    end_idx = 0
    if nt_idx%3 == 0:
        start_idx = nt_idx
        end_idx = nt_idx + 3
    elif nt_idx%3 == 1:
        start_idx = nt_idx - 1
        end_idx = nt_idx + 2
    elif nt_idx%3 == 2:
        start_idx = nt_idx - 2
        end_idx = nt_idx + 1

    return [start_idx, end_idx]


def get_aa_diff(target, query):
    """Finds the amino acid differences between two aa sequences, target and query.

    PARAMS
    ------
    target, query: str; amino acid sequences. These logically must be of the same length.

    RETURNS
    -------
    aa_diff_ls: list of aa differences, expressed as "target" to "query", and the index.
    """
    print("DEPRECATED: Use get_aa_changes() instead.")

    #if len(target) != len(query):
    #    print("WARNING: length of input strings are different!")

    #aa_diff_ls = []
    #for i in range(len(target)):
    #    if target[i] != query[i]:
    #        # set index to 1 because aa coords are 1-based
    #        aa_diff_row = "".join([target[i], str(i+1), query[i]])
    #        aa_diff_ls.append(aa_diff_row)
    #return aa_diff_ls


def get_parent_name(my_tree, node_name):
    """Get the immediate parent of a node, in my_tree. Returns the parent name.

    PARAMS
    ------
    my_tree: input Biopython tree object, where my_node lives
    node_name: name of a Biopython clade object, member of my_tree.

    RETURNS
    -------
    node_parent_name: the name of the immediate parent of the internal node with
    name node_name
    """
    # Check that node_name is in my_tree at all
    # This should break the function:
    nodes_ls = my_tree.get_nonterminals() + my_tree.get_terminals()
    node_nm_ls = [x.name for x in nodes_ls]
    if node_name not in node_nm_ls:
        print("ERROR: node %s not in tree!" % node_name)

    # Get the query clade object
    # Should have a test here to check that internal node names are unique
    node_parent_name = "*"
    if len(my_tree.get_path(node_name)) > 1:
        node_parent = my_tree.get_path(node_name)[-2]
        node_parent_name = node_parent.name
    else:
        print("NOTE: node %s does not have a parent. Returning '*'." % node_name)

    return node_parent_name


def prep_baseml_ctl(aln_fn, tre_fn, log_fn="baseml.out", ctl_fn="baseml.ctl"):
    """As in the name, prepare a baseml ctl for ancestral reconstruction.
    Reverse-engineered (copied) from treesub.

    PARAMS
    ------
    aln_fn: input alignment file name
    tre_fn: input tree file name
    log_fn: output file name of the baseml run.

    RETURNS
    -------
    ctl_fn: the baseml control file.
    """

    contents = ["seqfile = "+aln_fn]
    contents.append("treefile = "+tre_fn)
    contents.append("outfile = "+log_fn)

    contents.append("noisy = 1   * 0,1,2,3: how much rubbish on the screen")
    contents.append("verbose = 0   * 1: detailed output, 0: concise output")
    contents.append("runmode = 0   * 0: user tree")
    contents.append("model = 4   * 4:HKY85")
    contents.append("Mgene = 4")

    contents.append("clock = 0   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis")
    contents.append("fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below")
    contents.append("kappa = 5  * initial or fixed kappa")
    contents.append("fix_alpha = 0   * 0: estimate alpha; 1: fix alpha at value below")
    contents.append("alpha = 0.5   * initial or fixed alpha, 0:infinity (constant rate)")
    contents.append("Malpha = 0   * 1: different alpha's for genes, 0: one alpha")
    contents.append("ncatG = 5   * # of categories in the dG, AdG, or nparK models of rates")
    contents.append("nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK")

    contents.append("nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2")
    contents.append("getSE = 0   * 0: no, 1: get S.E.s of estimates")
    contents.append("RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states")

    contents.append("Small_Diff = 7e-11")
    contents.append("cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)")
    contents.append("icode = 0  * (with RateAncestor=1. try GC in data,model=4,Mgene=4)")
    contents.append("*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed")
    contents.append("method = 1  * Optimization method 0: simultaneous; 1: one branch a time")

    with open(ctl_fn, "w") as f:
        for line in contents:
            f.write("%s\n" % line)
