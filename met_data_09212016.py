# coding: utf-8
"""doc string"""
# Import modules
from copy import deepcopy, copy
import argparse
import os
import numpy as np
import pandas as pd

# Set constants
ROOTDIRECT = os.getcwd()
print(ROOTDIRECT)
try:
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("input_file")
    ARGS = PARSER.parse_args()
    INPUT = ARGS.input_file
    print(INPUT)
except:
    INPUT = input("Define input file: ") # Takes filename as commandline input
INPUTFILE = os.path.join(ROOTDIRECT, INPUT)
OUTPUTFILE = os.path.join(ROOTDIRECT, "storms_"+INPUT)
OUTPUTFILE2 = os.path.join(ROOTDIRECT, "storms_raw_"+INPUT)
WAVE_HEIGHT = 3             # meters
DISTINCT_EVENT_LENGTH = 24  # hours

# Load data
if int(INPUTFILE.split(".")[0][-4:]) < 1999:
    SKIPROWS = 1
    WIDTHS = [2, 3, 3, 3, 4, 5, 5, 6, 6, 6, 4, 7, 6, 6, 6, 5]
    NAMES = ['YY', 'MM', 'DD', 'hh', 'WDIR', 'WSPD', 'GST', 'WVHT', 'DPD',
             'APD', 'MWD', 'PRES', 'ATMP', 'WTMP', 'DEWP', 'VIS']
elif int(INPUTFILE.split(".")[0][-4:]) < 2000:
    SKIPROWS = 1
    WIDTHS = [4, 3, 3, 3, 4, 5, 5, 6, 6, 6, 4, 7, 6, 6, 6, 5]
    NAMES = ['YY', 'MM', 'DD', 'hh', 'WDIR', 'WSPD', 'GST', 'WVHT', 'DPD',
             'APD', 'MWD', 'PRES', 'ATMP', 'WTMP', 'DEWP', 'VIS']
elif int(INPUTFILE.split(".")[0][-4:]) < 2005:
    SKIPROWS = 1
    WIDTHS = [4, 3, 3, 3, 4, 5, 5, 6, 6, 6, 4, 7, 6, 6, 6, 5, 6]
    NAMES = ['YY', 'MM', 'DD', 'hh', 'WDIR', 'WSPD', 'GST', 'WVHT', 'DPD',
             'APD', 'MWD', 'PRES', 'ATMP', 'WTMP', 'DEWP', 'VIS', 'TIDE']
elif int(INPUTFILE.split(".")[0][-4:]) < 2007:
    SKIPROWS = 1
    WIDTHS = [4, 3, 3, 3, 3, 4, 5, 5, 6, 6, 6, 4, 7, 6, 6, 6, 5, 6]
    NAMES = ['YY', 'MM', 'DD', 'hh', 'mm', 'WDIR', 'WSPD', 'GST', 'WVHT', 'DPD',
             'APD', 'MWD', 'PRES', 'ATMP', 'WTMP', 'DEWP', 'VIS', 'TIDE']
else:
    SKIPROWS = 2
    WIDTHS = [4, 3, 3, 3, 3, 4, 5, 5, 6, 6, 6, 4, 7, 6, 6, 6, 5, 6]
    NAMES = ['YY', 'MM', 'DD', 'hh', 'mm', 'WDIR', 'WSPD', 'GST', 'WVHT', 'DPD',
             'APD', 'MWD', 'PRES', 'ATMP', 'WTMP', 'DEWP', 'VIS', 'TIDE']

DATA = pd.read_fwf(INPUTFILE,
                   widths=WIDTHS,
                   skiprows=SKIPROWS,
                   names=NAMES,
                   na_values=["99", "99.0", "99.00", "999", "999.0", "999.00"])

# Define functions
def select_storms():
    """finding significant storm events"""
    i = 0
    storm_index = 0

    storm = {"hs":float,    # Max 1-hr significant wave height (m)
             "tps":float,   # Dominant wave period at max 1-hr hs (sec)
             "a_hs":float,  # Average significant wave height (m)
             "a_tps":float, # Average dominant wave period (sec)
             "begin":int,   # Beginning of storm
             "length":int,  # Length of storm (hrs)
             "year":int,    # Year
             "month":int,   # Month
             "day":int,     # Day
             "hour":int}    # Hour

    storms = {}

    while i < len(DATA):
        if DATA["WVHT"][i] > WAVE_HEIGHT:
            storm["hs"] = DATA["WVHT"][i]
            storm["tps"] = DATA["DPD"][i]
            storm["begin"] = i
            storm["length"] = 1
            storm["year"] = DATA["YY"][i]
            storm["month"] = DATA["MM"][i]
            storm["day"] = DATA["DD"][i]
            storm["hour"] = DATA["hh"][i]
            storm["a_hs"] = DATA["WVHT"][i]
            storm["a_tps"] = DATA["DPD"][i]
            running_index = i + 1
            current_wave_height = DATA["WVHT"][running_index]
            while current_wave_height > WAVE_HEIGHT:
                if DATA["WVHT"][running_index] > storm["hs"]:
                    storm["hs"] = DATA["WVHT"][running_index]
                    storm["tps"] = DATA["DPD"][running_index]
                storm["length"] += 1
                storm["a_hs"] = ((storm["a_hs"]*(storm["length"]-1))+\
                                  DATA["WVHT"][running_index])/storm["length"]
                storm["a_tps"] = ((storm["a_tps"]*(storm["length"]-1))+\
                                   DATA["DPD"][running_index])/storm["length"]
                running_index += 1
                if running_index >= len(DATA):
                    current_wave_height = 0
                else:
                    current_wave_height = DATA["WVHT"][running_index]
            if storm["length"] > 2:
                storms[storm_index] = deepcopy(storm)
                storm_index += 1
            i = running_index
        else:
            i += 1
    return storms

def combine_close_storms(storms_raw):
    """decide if each storm event is separate or should be combined,
       based on minimum storm interval criteria"""
    storm = {"hs":float,     # Max 1-hr significant wave height (m)
             "tps":float,    # Dominant wave period at max 1-hr hs (sec)
             "a_hs":float,  # Average significant wave height (m)
             "a_tps":float, # Average dominant wave period (sec)
             "begin":int,    # Beginning of storm
             "length":int,   # Length of storm (hrs)
             "year":int,     # Year
             "month":int,    # Month
             "day":int,      # Day
             "hour":int}     # Hour

    storms = {}
    storm_index = 0
    i = 0
    if len(storms_raw) == 0:
        return storms_raw
    else:
        while i < len(storms_raw)-1:
            between_2_storms = storms_raw[i+1]["begin"]-\
                              (storms_raw[i]["begin"]+\
                               storms_raw[i]["length"])
            if  between_2_storms < DISTINCT_EVENT_LENGTH:
                hs_list = [storms_raw[i]["hs"], storms_raw[i+1]["hs"]]
                max_ind = np.argmax(hs_list)
                storm["hs"] = copy(storms_raw[i+max_ind]["hs"])
                storm["tps"] = copy(storms_raw[i+max_ind]["tps"])
                storm["begin"] = copy(storms_raw[i]["begin"])
                new_length = storms_raw[i+1]["begin"]-\
                             storms_raw[i]["begin"]+\
                             storms_raw[i+1]["length"]
                storm["length"] = copy(new_length)
                storm["year"] = copy(storms_raw[i]["year"])
                storm["month"] = copy(storms_raw[i]["month"])
                storm["day"] = copy(storms_raw[i]["day"])
                storm["hour"] = copy(storms_raw[i]["hour"])

                storm["a_hs"] = (storms_raw[i]["a_hs"]*storms_raw[i]["length"]+\
                                 storms_raw[i+1]["a_hs"]*storms_raw[i+1]["length"])/\
                                 (storms_raw[i]["length"]+storms_raw[i+1]["length"])
                storm["a_tps"] = (storms_raw[i]["a_tps"]*storms_raw[i]["length"]+\
                                 storms_raw[i+1]["a_tps"]*storms_raw[i+1]["length"])/\
                                 (storms_raw[i]["length"]+storms_raw[i+1]["length"])

                storms[storm_index] = deepcopy(storm)
                storm_index += 1
                i += 2
            else:
                storms[storm_index] = deepcopy(storms_raw[i])
                storm_index += 1
                i += 1
        if i < len(storms_raw):
            storms[storm_index] = deepcopy(storms_raw[i])
        else:
            pass
        return storms

def write_files(storms, outputfile):
    """Write list of storms to output file"""
    with open(outputfile, "w") as output:
        line = "    Year   Month     Day    Hour   begin  length    hsig     tps  a_hsig   a_tps\n"
        output.write(line)
        for storm in storms:
            year = '{:>8}'.format(storms[storm]["year"])
            month = '{:>8}'.format(storms[storm]["month"])
            day = '{:>8}'.format(storms[storm]["day"])
            hour = '{:>8}'.format(storms[storm]["hour"])
            begin = '{:>8}'.format(storms[storm]["begin"])
            length = '{:>8}'.format(storms[storm]["length"])
            hsig = '{:>8.2f}'.format(storms[storm]["hs"])
            tps = '{:>8.2f}'.format(storms[storm]["tps"])
            a_hsig = '{:>8.2f}'.format(storms[storm]["a_hs"])
            a_tps = '{:>8.2f}'.format(storms[storm]["a_tps"])
            line = year+month+day+hour+begin+length+hsig+tps+a_hsig+a_tps+"\n"
            output.write(line)

# Define main function
def main():
    """lorem ipsum"""
    storms_raw = select_storms()
    storms = combine_close_storms(storms_raw)
    storms1 = combine_close_storms(storms)
    storms2 = combine_close_storms(storms1)
    write_files(storms_raw, OUTPUTFILE2)
    write_files(storms2, OUTPUTFILE)

if __name__ == '__main__':
    main()
