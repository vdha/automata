import numpy as np
import pandas as pd
import re

# date conversion imports
from dateutil.relativedelta import relativedelta as rd
import datetime
from datetime import timedelta
import time

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