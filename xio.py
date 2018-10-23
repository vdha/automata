"""Library related to any kind of file reading/writing, and
some data wrangling tools.
"""

import os
import pandas as pd
import numpy as np
from math import factorial
import json

from Bio import SeqIO


def read_fasta(fn, delimiter="|", preview=0):
    """Reads a fasta into a list of lists, where each list is one record.
    Best practice: fasta file headers comprise of different metadata columns, delimited by the specified 
    delimiter. e.g.:

    >id1|Japan|2015-05-31
    agctagctagct
    >id2|Australia|2016-10-20
    agctagctagct
    >id2|NewZealand|2014-10-25
    agctagctagct

    You can then load the output of this function into a pandas dataframe, like so:

    >>contents = xio.read_fasta(input_fn.fas, delimiter="|", preview=0)
    >>fas_cols = ["id", "country", "collection_date", "sequence"]
    >>df = pd.DataFrame(data=contents, columns=fas_cols)

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


def read_beast_logfile(fn, burnin=0.1):
    """Reads a BEAST output log file into a table. Ignores all the stuff on top; goes straight
    to the table. Verified for BEAST2.5. 
    
    PARAMS
    ------
    fn: str; input file name
    burnin: float between 0 and 1; percentage burnin. 
    
    OUTPUT
    ------
    df: pandas dataframe. 
    """
    
    with open(fn) as f:
        contents = f.readlines()
    contents = [x.strip() for x in contents]

    # Retrieve only the logger contents
    contents2 = []
    log_bool = False
    for line in contents:
        if "Sample" in line:
            log_bool = True
        if log_bool:
            contents2.append(line)
    contents2 = [x.split("\t") for x in contents2]

    # Remove burnin - first 10% of rows
    n_burnin = round(len(contents2)*burnin)

    # Put it in a df
    df = pd.DataFrame(data=contents2[n_burnin:], columns=contents2[0])
    
    return df



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


def write_fasta(fn_out, d_in, header_cols, seq_col="seq", delimiter="|", preview=0, verbose=True):
    """Wrapper for prep_fasta_contents(), and does the additional step of
    writing out to a file.

    Params
    ------
    d_in: input dataframe
    header_cols: list of str. Whatever columns from d_in to put in the fasta header
    seq_col: str. The name of the column which contains actual sequence data.
        default = 'seq'.
    preview: int. No. of records to preview after the function is finished.
    fn_out: str; name of file to write out to.
    verbose: bool; verbosity

    Returns
    -------
    Writes out a fasta file.
    """
    contents = prep_fasta_contents(d_in, header_cols, seq_col=seq_col, delimiter=delimiter, preview=preview)
    with open(fn_out, "w") as f:
        for line in contents:
            f.write("%s\n" % line)
    if verbose:
    	print("Wrote to file %s" % fn_out)


def remove_duplicates(df, colname, seq_colname):
    """
    EXPERIMENTAL
    Removes duplicate fields from df, given by "colname", breaking ties by selecting one
    with the longest 'seq' string. Note that gaps are still counted. 

    Only works for dfs with the same segment. A df with multiple segments would have in line 4:
    >>df = df.loc[df.seq.str.len().groupby([df.colname, df.seg]).idxmax()]

    In order to differentiate between segments. But that shouldn't be a problem if the segment
    is already in 'colname'...should it? Needs unit testing!

    Params
    ------
    df: pandas dataframe of a fasta file, probably read by read_fasta()
    colname: column name in df to select duplicate names for.
    seq_colname: column name of the sequence data in df to select the longest sequence by. 

    Returns
    -------
    df: dataframe with duplicates removed.
    """

    #Why does this work? (answer from SO)
    n_records0 = df.shape[0]
    print("Initial no. of records = %s" % len(df))
    n_uq0 = len(set(df[colname]))
    df = df.loc[df[seq_colname].str.len().groupby([df[colname]]).idxmax()]
    n_records1 = df.shape[0]
    n_uq1 = len(set(df[colname]))
    print("n_unique name_ids (before, after) = %s, %s" % (n_uq0, n_uq1))
    print("No. of records discarded = %s" % (n_records0-n_records1))

    return df


def read_json(json_fn):
    """An almost unnecessary function to read a json file, but I keep forgetting the full syntax.

    Params
    ------
    json_fn: str; path to json file to be read into a dictionary

    Returns
    -------
    my_json_dict: dictionary
    """
    with open(json_fn) as f:
        json_str = f.read()
    my_json_dict = json.loads(json_str)
    return my_json_dict


def lookup_key(my_dict, val):
    """Look up a value, val, from a dictionary of the format key:list_of_values. Values cannot appear in
    multiple lists.
    TO DO: exception handling!

    Params
    ------
    my_dict: dictionary of the format k:list_of_values
    val: value to lookup, to find corresponding key

    Returns
    -------
    my_key: key which value belongs to
    """
    my_key = "*"
    for k in list(my_dict.keys()):
        if val in my_dict[k]:
            my_key = k
    return my_key


def pivot_raw_tbl(df):
    """Pivots a GISAID dataframe of raw data such that each segment is its own column. That is, turns this:

    name_id      | segment | seq
    -------------|---------|----
    roger        |  PB1    | agct...
    roger        |  PB2    | ggtc...
    roger        |  HA     | aatc...

    into this:

    name_id      |   PB1   |   PB2   |   HA
    -------------|---------|---------|--------
    roger        | agct... | ggtc... | aatc...

    The column "name_id" is a concatenation of the isolate name, and its Epiflu ID, delimited by '|'. This is
    necessary as some isolate names can have multiple IDs for some reason (but, thankfully, not the other way
    round.)

    Params
    ------
    df: a pandas dataframe. Must have the 'name_id' and 'segment' columns.

    Returns
    -------
    df_piv: the requested pivoted dataframe.
    """

    #d0a: the partition that will remain unpivoted. d0b: the partition to be pivoted.
    dfa_selection = []; dfb_selection = []
    # select the columns for the partition that will remain unpivoted
    for col in df.columns.values:
        if col not in ['segment', 'seq', 'segment_number']:
            dfa_selection.append(col)

    dfa = df[dfa_selection].drop_duplicates()
    dfb = df[['name_id', 'segment', 'seq']].drop_duplicates()
    dfb_T = dfb.pivot(index='name_id', columns='segment', values='seq')
    dfb_T.columns.name=None
    dfb_T.reset_index(inplace=True)

    df = pd.merge(dfa, dfb_T, on='name_id')
    return df


def location_split(loc_ls, max_loc_len=5, verbose=False):
    """Splits a list of locations, where each location is in the format:
    'continent/country/state/city/district', into an array of 5. Location
    entries with less than 5 elements are padded with a '*' to represent empty
    values.

    Params
    ------
    loc_ls: list of str; a list of locations.
    max_loc_len: int; number of levels of the location hierarchy.

    Returns
    -------
    loc_arr: str array, shape (n_locations, 5).
    """

    loc_arr = []

    counter = 0
    for line in loc_ls:
        ln = line.split("/")
        for i in range(len(ln)):
            ln[i] = ln[i].strip()
        # Pad ln with '*'
        ln = ln + ["*"]*(max_loc_len-len(ln))

        if verbose:
            if counter < 5:
                print(ln)
                counter +=1
        loc_arr.append(ln)

    return np.array(loc_arr)


def adjust_raw_loc(df):
    """Wrapper for location_split(). Formats the existing 'location' column,
    and adds the location array as 5 new columns.
    """
    loc_ls = list(df.location)
    for i in range(len(loc_ls)):
        loca = loc_ls[i].split("/")
        loca = [i.strip() for i in loca]
        loca = "/".join(loca)
        loc_ls[i] = loca
    df["location"] = loc_ls

    loc_arr = location_split(loc_ls)
    df["continent"] = loc_arr[:, 0]
    df["country"] = loc_arr[:, 1]
    df["state"] = loc_arr[:, 2]
    df["city"] = loc_arr[:, 3]
    df["district"] = loc_arr[:, 4]

    return df


def xls_to_csv(fn):
    fn_out = fn[:-4] + ".csv"
    wb = xlrd.open_workbook(fn)
    sheet = wb.sheet_by_name('Tabelle1')
    with open(fn_out, "w") as f:
        wr = csv.writer(f, quoting=csv.QUOTE_ALL)
        for rownum in range(sheet.nrows):
            wr.writerow(sheet.row_values(rownum))
