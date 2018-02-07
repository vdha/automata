import numpy as np
import pandas as pd
import re

# date conversion imports
from dateutil.relativedelta import relativedelta as rd
import datetime
from datetime import timedelta
import time

"""A very, very general library of utility functions which don't have a home yet.
We should find them one, or make them one, at some point.
"""

def sample_by_col(df, col, n_sample):
    df1 = df.groupby(col, group_keys=False).apply(lambda g: g.sample(n_sample)
                                                  if len(g) > n_sample else g)
    return df1


def translate_codon(triplet):
    """Translates an amino acid triplet into a codon.
    A triple-nuc gap is explicitly handled as a codon gap.
    Otherwise, if the triplet is anything other than a valid nucleotide triplet,
    return 'X'.
    Uracil currently not handled.

    Params
    ------
    triplet: str. Any string of length at most 3. If len(triplet) < 3, return ''.

    Returns
    -------
    codon: str of length 1. The amino acid translation.
    """
    codon_dict = {'F':['TTT', 'TTC'],
                  'L':['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
                  'I':['ATT', 'ATC', 'ATA'],
                  'M':['ATG'],
                  'V':['GTT', 'GTC', 'GTA', 'GTG'],
                  'S':['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                  'P':['CCT', 'CCC', 'CCA', 'CCG'],
                  'T':['ACT', 'ACC', 'ACA', 'ACG'],
                  'A':['GCT', 'GCC', 'GCA', 'GCG'],
                  'Y':['TAT', 'TAC'],
                  'X':['TAA', 'TAG', 'TGA'], #stop; I just used X: Ochre, Amber, Opal resp.
                  'H':['CAT', 'CAC'],
                  'Q':['CAA', 'CAG'],
                  'N':['AAT', 'AAC'],
                  'K':['AAA', 'AAG'],
                  'D':['GAT', 'GAC'],
                  'E':['GAA', 'GAG'],
                  'C':['TGT', 'TGC'],
                  'W':['TGG'],
                  'R':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                  'G':['GGT', 'GGC', 'GGA', 'GGG']
                 }
    triplet = triplet.upper()

    # Init as unknown codon, X.
    # If not found in dictionary, return X
    codon = 'X'
    for k in list(codon_dict):
        if triplet in codon_dict[k]:
            codon = k

    # Explicit handling of a 3-nuc gap, and input length < 3
    if triplet == '---':
        codon = '-'
    elif len(triplet) < 3:
        codon = ''

    return codon


def translate_aa_seq(seq):
    """Wrapper function for translte_codon().
    Translates a nucleotide sequence.

    >> test_cases = ["AGT", "AZT", "YYT", "-AT",
              "AGTCATAGTCAT",
              "AGTYATAGTCAT",
              "AGT--TAGTCAT",
              "AGT-ATAGTCAT",
              "AGT---AGTCAT",
              "AGT------AGT"]
    >> for seq in test_cases:
    >>    print("%s --> %s" % (seq, my_utils.translate_aa_seq(seq)))
    AGT --> S
    AZT --> X
    YYT --> X
    -AT --> X
    AGTCATAGTCAT --> SHSH
    AGTYATAGTCAT --> SXSH
    AGT--TAGTCAT --> SXSH
    AGT-ATAGTCAT --> SXSH
    AGT---AGTCAT --> S-SH
    AGT------AGT --> S--S
    """
    aa = ''
    for i in range(0, len(seq), 3):
        aa_temp = translate_codon(seq[i:i+3])
        aa = aa + aa_temp
    return aa


def choose_string_sections(seq, coords_ls):
    """Given a sequence seq, extract the relevant portions using the coords given in coords_ls

    Params
    ------
    seq: str. A nuc sequence.
    coords_ls: list of lists, where each element list is a set of coordinates.

    Returns
    -------
    sections_ls: list of string.

    >> seq = "ACGTACGTACGT"
    >> coords_ls = [[0, 4], [1, 5], [5, 8]]
    >> choose_gene_sections(seq, coords_ls)
    ['ACGT', 'CGTA', 'CGT']
    """
    section_ls = []
    for coord in coords_ls:
        s, t = coord
        section_ls.append(seq[s:t])
    return section_ls


def format_date(dmy, input_date_fmt="%d/%m/%y"):
    """Because Excel keeps reformatting dates into dd/mm/yy format,
    this function converts it into yyyy-mm-dd for printing purposes.

    Params
    ------
    dmy: str; a date in the format dd/mm/yyyy.
    input_date_fmt: format of the input date.

    Returns
    -------
    datum: str; a date in the format yyyy-mm-dd
    """
    if dmy.count("/") == 2:
        datum = datetime.datetime.strptime(dmy, "%d/%m/%y")
        datum = str(datum)[:10]
    elif dmy.count("/") == 1:
        datum = datetime.datetime.strptime(dmy, "%m/%y")
        datum = str(datum)[:10]
    elif dmy.count("/") == 0:
        datum = datetime.datetime.strptime(dmy, "%Y")
        datum = str(datum)[:10]
    return datum


def clean_dates(dates_ls, verbose=False):
    """Standardizes a list of dates into yyyy-mm-dd, split into a list.
    The raw date formats that don't conform to yyyy-mm-dd are assumed to have
    their first 4 digits as 'yyyy', as GISAID data tends to have.

    Params
    ------
    dates_ls: list of dates. Dates are excepted as type str.
    verbose: boolean.

    Returns
    -------
    clean_dates: 2D arr shape (n_dates, 3), where each date is [yyyy, mm, dd].
    """
    pattern = r"(\d{4})-(\d{2})-(\d{2})"
    clean_dates = []

    n_bad = 0
    n_good = 0
    for date in dates_ls:
        match = re.findall(pattern, date)
        match = [item for sublist in match for item in sublist] # flatten
        if len(match) == 3: # A good date
            clean_dates.append(match)
            n_good += 1
        else:
            n_bad += 1
            clean_dates.append([date[:4], '*', '*'])

    print("n_good, n_bad = %s, %s" % (n_good, n_bad))

    return np.array(clean_dates)


def hms_time(x):
    """Converts x, in seconds, into an (h, m, s) list.
    """
    m, s = divmod(x, 60)
    h, m = divmod(m, 60)
    return [h, m, s]


def datetime_to_decimal_year(date):
    """Some SO solution I haven't really looked into yet.

    PARAMS
    ------
    date: datetime object

    returns
    -------
    decimal year: float. the decimal year version of the input.

    e.g.
    >>> T = datetime.date(2016, 11, 29)
    >>> most_recent = to_decimal_year(T)
    >>> most_recent
    2016.9098360655737
    """
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime.datetime(year=year, month=1, day=1)
    startOfNextYear = datetime.datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction


def to_decimal_year(date_str, date_delimiter="-"):
    """Wraps datetime_to_decimal_year, with the additional precursor step
    of converting the input date_str into a datetime.date object.

    Params
    -------
    date_str: str, date in the format yyyy-mm-dd, yyyy-mm, or yyyy
    date_delimiter: str, separator for the input date_str

    Returns
    -------
    decimal year: float. the decimal year version of the input.
    """
    date_ln = date_str.split(date_delimiter)
    T = "*"
    if len(date_ln) == 3:
        T = datetime.date(int(date_ln[0]), int(date_ln[1]), int(date_ln[2]))
    elif len(date_ln) == 2:
        T = datetime.date(int(date_ln[0]), int(date_ln[1]), 1)
    elif len(date_ln) == 1:
        T = datetime.date(int(date_ln[0]), 1, 1)

    return datetime_to_decimal_year(T)


def to_calendar_year(cdate, return_type="date"):
    """Converts a decimal date to a calendar date.

    Params
    ------
    cdate: float; decimal year
    return_type: "date" returns a datetime object with just y-m-d,
        'full' returns a datetime object with h:m:s as well.

    Returns
    -------
    result: datetime object, with resolution specified by the input param
        'return_type'
    """
    year = int(cdate)
    rem = cdate - year

    base = datetime.datetime(year, 1, 1)
    full_date = base + timedelta(seconds=(base.replace(year=base.year + 1) - base).total_seconds() * rem)

    if return_type == "date":
        y = full_date.year
        m = full_date.month
        d = full_date.day
        date_vals_ls = [y, m, d]
        result = datetime.datetime(*date_vals_ls)
    elif return_type == "full":
        result = full_date

    return result
