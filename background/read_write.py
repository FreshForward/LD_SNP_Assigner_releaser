import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import time
import numpy as np
import sys
import os
import h5py

def openFile(filelocation, marker_column, marker_pos_data, NA_val, chunk=None):
    """
    Function opens a given file, returns a Pandas Dataframe.
    Takes CSV file or HDF5 file location with name of the marker column.
    :PARAM: (string) filelocation: Contains location of to- read file.
                                    Can be either hdf5 or csv file.
    :PARAM: (string) marker_column: Contains the name of the marker column
    :PARAM: (array) chunk: Contains the chunk of the hdf5 file to be read.
    :RETURNS: (Pandas Dataframe) containins markers, positions & sequenced SNPs
    """
    start_time = time.time()
    #print("Opening file, starting count...")

    # Marker_column car uses "Marker" column to set index
    if filelocation.endswith(".csv"):
        with open(filelocation, "r") as f:
            data = pd.read_csv(f, na_values=NA_val)
            try:
                data.set_index(marker_column, inplace=True)
                #data.replace(NA_val, np.NaN)
            except KeyError:
                print("KeyError: None of '"+marker_column+
                      "' found in CSV file '"+filelocation+"'.\nHas the file been edited in excel? Please double check your file has a comma separator.\nPrinting found columns & exiting...\n")
                print(data.columns)
                print("If there is a big blob of data above, please double check the separator and ensure it's a comma.")
                sys.exit()
    else:
        hf = h5py.File(os.path.join(filelocation), 'r')
        if chunk is None:
            corrset1 = hf['corrset']
        else:
            corrset1 = hf['corrset'][..., chunk]
        data = pd.DataFrame(corrset1, columns=None, index=None)
        hf.close()
    return data


def write_corrmatrix_hdf(corrset, corrmat_name):
    """
    Writes matrix to HDF5file. Takes matrix and name to write. 

    Args:
        corrset (pandas dataframe): Matrix containing correlation scores, markers,
                                    position, chrom position & est chrom position.
        corrmat_name     (string) : Contains the name of the to be written file.
    """
    #start_time = time.time()
    #print("Starting data writing")
    os.makedirs(os.path.dirname(corrmat_name), exist_ok=True) #mkdir if not exist
    #corrset.to_csv((corrmat_name+".csv"))
    corrset.to_hdf((corrmat_name+".h5"), key = "df", mode = "w", index=True, complevel=6)
    print("LD matrix saved as: ", corrmat_name, ".h5",  sep="")
    
def write_statmat_csv(statistics, statfile_name, decimalnums):
    """
    Writes statistics variable dataframe to CSV file.

    Args:
        statistics      (pandas df): Matrix with position, markers, est position 
                                       and total calculated scores of the results.
        statfile_name     (string) : Contains the name of the to be written file.
        decimalnums:         (int) : Number of decimals written after the comma. 
                                     Correlates to file size.         
    """
    os.makedirs(os.path.dirname(statfile_name), exist_ok=True) #mkdir if not exist
    (statistics.round(decimals=decimalnums)).to_csv(statfile_name+".csv")
    print("Wrote statistics file as: ", statfile_name,".csv",sep="")