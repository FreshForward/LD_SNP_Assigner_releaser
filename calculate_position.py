import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import background.functions_Mainfile as func
import background.read_write as rw
import yaml
import os
import sys
import time
import background.preprocessing as prep
import numpy as np
import math


def main():
    time_start = time.time()
    try:
        if sys.argv[1] and sys.argv[1].endswith(".yaml"):
            with open(sys.argv[1], 'r') as configfile:
                config = yaml.safe_load(configfile)
    except IndexError:
        with open('config.yaml', 'r') as configfile:
            config = yaml.safe_load(configfile)

    chr_col_name = config["chr_col_name"] # Name of existing predicted chromosome column
    pred_chr_col_name = config["pred_chr_col_name"] # Column name to be created for predicted chromosome
    remove_columns = config["remove_columns"] # Remove chromosome value / column to predict
    position_col_name = config["position_col_name"] 
    est_pos_col_name = config["est_pos_col_name"]
    abs_diff_pos_name = config["abs_diff_pos_name"]
    datapoints_used_for_yadav_col = config["datapoints_used_for_yadav_col"] # Column name
    calculate_all_or_different = config["calculate_all_or_different_pos"] # 'all' or 'different': calc position for all elements or just changed elements.
                                       # This is more important for time constraints.
    #partdat = config["partdata"]
    NA_val=config["NA_val"]
    run_fullpipe = config["run_fullpipe"]


    # if partdat == True:
    #     marker_column = "Marker"
    #     file = "PartResults/CorrelationMatrixFile_RoyRoyce__PartDat_top10.csv"
    #     file2 = "PartResults/CorrelationMatrixFile_RoyRoyce__PartDat_top10.hdf5"
    #     statisticsfilename = "PartResults/Statisticsfile_RoyRoyce_PartDat_top10.csv"
    #     statfileresult = "PartResults/Statistics__PartDat_top10_Position"
    #     corrmatresult = "PartResults/Corrmat__PartDat_top10_Position"

    data_type = config["data_type"]
    de_novo_or_anchor = config["de_novo_or_anchor"]
    sumchrom_calc_method = config["sumchrom_calc_method"]
    run_ID =  config["run_ID"]
    marker_column = config["marker_column"]
    file2 = config["corrmat_name"].replace("{{ data_type }}", data_type).replace("{{ sumchrom_calc_method }}", sumchrom_calc_method).replace("{{ run_ID }}", run_ID)+".h5"
    statisticsfilename = config["statfile_name"].replace("{{ data_type }}", data_type).replace("{{ sumchrom_calc_method }}", sumchrom_calc_method).replace("{{ run_ID }}", run_ID) +".csv"
    statfileresult = config["statfile_name_with_position"].replace("{{ data_type }}", data_type).replace("{{ sumchrom_calc_method }}", sumchrom_calc_method).replace("{{ run_ID }}", run_ID)
    corrmatresult = config["corrmat_name_with_position"].replace("{{ data_type }}", data_type).replace("{{ sumchrom_calc_method }}", sumchrom_calc_method).replace("{{ run_ID }}", run_ID)
    corrmat_name = config["corrmat_name"].replace("{{ data_type }}", data_type).replace("{{ sumchrom_calc_method }}", sumchrom_calc_method).replace("{{ run_ID }}", run_ID)
    filelocation = config["sourcefile_location"]
    datacolumns = config["datacolumns"] # Number of columns containing data at the start of the file
    splitlength = config ["splitlength"]
    
    data1 = rw.openFile(filelocation, marker_column, None, NA_val)

    if de_novo_or_anchor == "anchor":
        checkdata = config["checkfile_location"]
        data2 = rw.openFile(checkdata, marker_column, None, NA_val)
        marker_pos_data2 = prep.get_marker_pos_data(data2, datacolumns)
        marker_pos_data = prep.get_marker_pos_data(data1, datacolumns)
    else:
        marker_pos_data = prep.get_marker_pos_data(data1, datacolumns)
        marker_pos_data2 = marker_pos_data

    # Open files using external function  
    # if de_novo_or_anchor == "anchor":
    #     colnames = pd.read_csv(os.path.join(corrmat_name + ".column_names.txt"),sep="\t",header=None, index_col=None)
    #     colnames = np.array(colnames[0])
    #     length = len(colnames)
    # else:
    #     index = pd.read_csv(os.path.join(corrmat_name + ".index_names.txt"),sep="\t",header=None, index_col=None)
    #     index = np.array(index[0])
    #     length=len(index)
    statistics1 = rw.openFile(statisticsfilename, marker_column, None, NA_val)

    colnames = pd.read_csv(os.path.join(corrmat_name + ".column_names.txt"),sep="\t",header=None, index_col=None)
    colnames = np.array(colnames[0])
    length=len(colnames)
    print("POS_EST: Starting estimating SNP positions.")
    if splitlength < length: 
        # time_start = time.time()
        statistics1 = rw.openFile(statisticsfilename, marker_column, None, NA_val)
        print("POS_EST: Splitting data into parts of", splitlength, "SNPs to calculate positions...")
        parts_y = math.ceil(length/splitlength) #calculate how many parts data should be divided if max snps in a single part is 10000  
        splits_y = list(np.array_split(range(length),parts_y)) #get ranges of splits
        for splitnr, split in enumerate(splits_y): 
            corrset = rw.openFile(os.path.join(str(corrmat_name) + ".h5"), marker_column, marker_pos_data ,NA_val,chunk=split)
            # corrset=corrset.to_numpy()
            if de_novo_or_anchor == "de_novo":
                print("!!! WARNING !!! \n\n Filling diagonal with NaNs is not possible for split matrices. Please use the full matrix. \n\n")
            corrset = pd.DataFrame(corrset)
            index = pd.read_csv(os.path.join(corrmat_name + ".index_names.txt"),sep="\t",header=None, index_col=None)
            index = index[0].tolist()
            colnames = pd.read_csv(os.path.join(corrmat_name + ".column_names.txt"),sep="\t",header=None, index_col=None)
            colnames = np.array(colnames[0])
            corrset.index = index
            corrset.columns = colnames[split]
            data = pd.concat([marker_pos_data,corrset],axis=1, join="inner")

            # Add predicted chrom pos to corrset
            
            statistics=statistics1.loc[colnames[split]]
            data.insert(data.columns.get_loc(chr_col_name)+1, pred_chr_col_name, statistics.loc[:, pred_chr_col_name])
            difdata, notdif = func.check_difference(data, chr_col_name, pred_chr_col_name)
            
            data, statistics = func.calculate_positions(data, statistics, 
                                                    position_col_name,
                                                    marker_pos_data2.loc[colnames[split],:], 
                                                    est_pos_col_name,
                                                    datapoints_used_for_yadav_col,
                                                    calculate_all_or_different, 
                                                    difdata, notdif, abs_diff_pos_name)
            if splitnr == 0:
                statistics2 = statistics
            else:
                statistics2 = pd.concat([statistics2,statistics])
        
        #add markers that are not in corrset (due to filtering like MAF, of missing data) to statistics
        for marker in marker_pos_data2.index:
            if marker not in colnames:
                statistics2.loc[marker, chr_col_name] = marker_pos_data2.loc[marker, chr_col_name]
                statistics2.loc[marker, position_col_name] = marker_pos_data2.loc[marker, position_col_name]      

        statistics2.sort_values(by=[chr_col_name, position_col_name, marker_column], inplace=True)

        statistics2.to_csv(statfileresult+".csv")

    else:
        # time_start = time.time()
        statistics = rw.openFile(statisticsfilename, marker_column, None, NA_val)
        corrset = rw.openFile(os.path.join(str(corrmat_name) + ".h5"), marker_column, marker_pos_data ,NA_val)
        corrset=corrset.to_numpy()
        if de_novo_or_anchor == "de_novo":
            np.fill_diagonal(corrset, np.nan) 
        corrset = pd.DataFrame(corrset)
        index = pd.read_csv(os.path.join(corrmat_name + ".index_names.txt"),sep="\t",header=None, index_col=None)
        index = index[0].tolist()
        colnames = pd.read_csv(os.path.join(corrmat_name + ".column_names.txt"),sep="\t",header=None, index_col=None)
        colnames = colnames[0].tolist()
        corrset.index = index
        corrset.columns = colnames
        data = pd.concat([marker_pos_data,corrset],axis=1, join="inner")

    
        # Add predicted chrom pos to corrset
        # print(data[marker_column])
    
        data.insert(data.columns.get_loc(chr_col_name)+1, pred_chr_col_name, statistics.loc[:, pred_chr_col_name])
        difdata, notdif = func.check_difference(data, chr_col_name, pred_chr_col_name)
        
        data, statistics = func.calculate_positions(data, statistics, 
                                                position_col_name,
                                                marker_pos_data2, 
                                                est_pos_col_name,
                                                datapoints_used_for_yadav_col,
                                                calculate_all_or_different, 
                                                difdata, notdif, abs_diff_pos_name)   
        statistics.to_csv(statfileresult+".csv")



    print("POS_EST: Estimated SNP positions in", round((time.time()-time_start)/60,1), "minutes.")

    #print("Writing statistics summary file...")
    # statistics.to_csv(statfileresult+".csv")
    print("POS_EST: Wrote statistics including positions as ",statfileresult,".csv",sep="")
    #data.to_csv(corrmatresult+".csv")
    #print("Correlation matrix written.")


    if run_fullpipe:
        try:
            if sys.argv[1] and sys.argv[1].endswith(".yaml"):
                os.system("python3 summaryfile_creation.py " + sys.argv[1])
        except IndexError:
            os.system("python3 summaryfile_creation.py")

main()