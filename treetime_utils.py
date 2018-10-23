import numpy as np
import pandas as pd
import re

# date conversion imports
from dateutil.relativedelta import relativedelta as rd
import datetime
from datetime import timedelta
import time

import xio
import date_utils as xdu

"""A bunch of utility functions associated with Richard Neher's treetime.
"""

def get_treetime_dates(fn_in, fn_out="",
                       name_col="header",
                       date_col_pos=-1,
                       header_delimiter="|",
                       colnames=["v_name", "coll_date"]):
    """Reads the input fasta file into a dataframe, and returns/writes out the
    date file that treetime wants as input.

    PARAMS
    ------
    fn_in: str; name of input fasta file.
    fn_out: str; name of output file. If left as an empty string, does not
    write out.
    name_col: str; name of the fasta header column in df.
    date_col_pos: int; index of the dates in the fasta headers.
    header_delimiter: str; delimiting character in the fasta headers.
    colnames: list of str of length 2; output column names.

    RETURNS
    -------
    d1: dataframe; output dataframe. Optionally written out to csv, depending
    on fn_out.
    """
    contents = xio.read_fasta(fn_in, delimiter="", preview=0)
    df = pd.DataFrame(data=contents, columns=["header", "seq"])
    df["cdate"] = df.apply(lambda row: str(row["header"]).split(header_delimiter)[date_col_pos],
                           axis=1)
    d1 = df[["header", "cdate"]]
    d1.columns=colnames

    if len(fn_out) > 0:
        d1.to_csv(fn_out, index=False)

    return d1


def get_treetime_input_datefile(fn_in, fn_out, date_delim="|", date_pos=-1):
    """Read a fasta file, and write out dates in a format that treetime likes.
    
    PARAMS
    ------
    fn_in: str; name of input fasta file.
    fn_out: str; name of output text file.
    date_delim: str; date delimiter in fasta headers.
    date_pos: position of date in fasta header.
    
    RETURNS
    -------
    Writes out .csv in a format that treetime can read: names in col0, dates in col1.
    """
    c = xio.read_fasta(fn_in, delimiter="", preview=0)
    df = pd.DataFrame(data=c, columns=["header", "seq"])
    df["cdate"] = df.apply(lambda row: str(row["header"]).split(date_delim)[date_pos], axis=1)
    df["ddate"] = df.apply(lambda row: xdu.to_decimal_year(str(row["cdate"])), axis=1)

    dm = df[["header", "ddate"]]
    dm.columns = ["vname", "ddate"]
    dm.to_csv(fn_out, index=False)


def get_lsd_dates(df, name_col, date_col, fn_out):
    """Writes out dates in a format that lsd understands.
    Writes out partial dates as absolute dates, thoughself.

    PARAMS
    ------
    df: input pandas DataFrame
    name_col: str; name of tree tipnames.
    date_col: str; name of decimal date column
    fn_out: str; name of output file
    """

    contents = [str(len(df))]
    for index, row in df.iterrows():
        line = str(row[name_col])+"\t"+str(row[date_col])
        contents.append(line)

    with open(fn_out, "w") as f:
        for line in contents:
            f.write("%s\n" % line)
