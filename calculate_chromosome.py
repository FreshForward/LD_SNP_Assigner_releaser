import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import time
import numpy as np
import os
import background.functions_Mainfile as func
import background.read_write as rw
import background.preprocessing as prep
import yaml
import sys
import math

def main():

    try:
        if sys.argv[1] and sys.argv[1].endswith(".yaml"):
            with open(sys.argv[1], 'r') as configfile:
                config = yaml.safe_load(configfile)
    except IndexError:
        with open('config.yaml', 'r') as configfile:
            config = yaml.safe_load(configfile)

    datacolumns = config["datacolumns"] # Number of columns containing data at the start of the file
    decimalnums = config["decimalnums"] # Number of decimal digits after each correlation
    remove_values = config["remove_values"] # Cutoff value
    chr_col_name = config["chr_col_name"] # Name of existing predicted chromosome column
    pred_chr_col_name = config["pred_chr_col_name"] # Column name to be created for predicted chromosome
    remove_columns = config["remove_columns"] # Remove chromosome value / column to predict
    position_col_name = config["position_col_name"] 
    remove_string_in_columns = config["remove_string_in_columns"] # String list in which files are to be removed
    secondchoice_statistic_col_name = config["secondchoice_statistic_col_name"]
    suffix_count_values = config["suffix_count_values"]
    maf_threshold = config["maf_threshold"]
    splitlength = config ["splitlength"]
    run_fullpipe = config["run_fullpipe"]
    sumchrom_calc_method = config["sumchrom_calc_method"]
    yadav = config["yadav"]

    run_ID =  config["run_ID"]
    data_type = config["data_type"]

    marker_column = config["marker_column"]
    filelocation = config["sourcefile_location"]
    statfile_name = config["statfile_name"].replace("{{ data_type }}", data_type).replace("{{ sumchrom_calc_method }}", sumchrom_calc_method).replace("{{ run_ID }}", run_ID)
    corrmat_name = config["corrmat_name"].replace("{{ data_type }}", data_type).replace("{{ sumchrom_calc_method }}", sumchrom_calc_method).replace("{{ run_ID }}", run_ID)
    subgenomemax_name = config["subgenomemax_name"]
    subgenomemin_name = config["subgenomemin_name"]
    cnt_subgenomemax_name = config["cnt_subgenomemax_name"]
    cnt_subgenomemin_name = config["cnt_subgenomemin_name"]
    allowed_nan_perc = config["allowed_nan_perc"]
    corrmethod = config["corrmethod"]
    extra_columns = config["extra_columns"]
    NA_val = config["NA_val"]


    data = rw.openFile(filelocation, marker_column, None, NA_val)

    ##check duplicate markerIDs
    if data.index.duplicated().sum() > 0:
        print(data.index.duplicated().sum(), "duplicate markernames found. \nQuick fix: removing every second occurence.")
        data = data[~data.index.duplicated(keep='first')]
    
    de_novo_or_anchor = config["de_novo_or_anchor"]
    if de_novo_or_anchor == "anchor":
        checkdata = config["checkfile_location"]
        checkdata_marker_col = config["checkdata_marker_col"]
        checkdata_chr_col = config["checkfile_chr_col"]
        checkdata_pos_col = config["checkfile_pos_col"]
        checkdata_extra_cols = config["checkfile_extra_cols"]
        marker_column = config["marker_column"]
        data2 = rw.openFile(checkdata, checkdata_marker_col,None,NA_val)
                
        remove_intersect = config["remove_intersect"] # Var indicates terminate or handle.
        data2 = prep.check_intersect(data, data2, remove_intersect)
        data = prep.check_intersect(data2, data, remove_intersect) 
    
    marker_pos_data = prep.get_marker_pos_data(data, datacolumns)
    
    if de_novo_or_anchor == "anchor":
        marker_pos_data2 = prep.get_marker_pos_data(data2, datacolumns)
        # marker_pos_data2.to_csv(os.path.join(corrmat_name + ".map2.txt"),sep="\t",header=None,index=False)
    else:
        marker_pos_data2 = marker_pos_data
    
    
    # marker_pos_data.to_csv(os.path.join(corrmat_name + ".map.txt"),sep="\t",header=None,index=False)
    
    if remove_string_in_columns: # If not empty
        data = prep.drop_cols_containing_x(data, remove_string_in_columns)
        print(data.shape[1]-2,"genotypes remaining after removing specified genotypes.")
        if de_novo_or_anchor == "anchor":
            data2 = prep.drop_cols_containing_x(data2, remove_string_in_columns)

    if de_novo_or_anchor == "anchor":
        print(data.head())
        print(data2.head())
        #check if columnIDs in both files are the same and keep the first two columns as first and second column       
        common_cols = list(set(data.columns[2:]).intersection(data2.columns[2:])) # Convert set to list
        # sort common_cols
        common_cols.sort()
        data = data[[chr_col_name, position_col_name] + common_cols]
        data2 = data2[[checkdata_chr_col, checkdata_pos_col] + common_cols]
        print(data.head())  
        print(data2.head())    
    
    if not os.path.isfile(os.path.join(corrmat_name + '.' + "h5")): #check if corrset already exists
        print("LD matrix does not exist yet.")
        if de_novo_or_anchor == "anchor":
            corrset = func.corrmatrix_calc(data, data2, datacolumns, marker_pos_data, corrmat_name,
                        marker_column, decimalnums, remove_values, maf_threshold,splitlength,
                        allowed_nan_perc, corrmethod)
        else:
            corrset = func.corrmatrix_calc(data, None, datacolumns, marker_pos_data, corrmat_name,
                        marker_column, decimalnums, remove_values, maf_threshold,splitlength,
                        allowed_nan_perc, corrmethod)    
        # rw.write_corrmatrix_hdf(corrset, corrmat_name)        
    else: 
        print("LD matrix already exists.")
        
    if de_novo_or_anchor == "anchor":
        colnames = pd.read_csv(os.path.join(corrmat_name + ".column_names.txt"),sep="\t",header=None, index_col=None)
        colnames = np.array(colnames[0])
        length = len(colnames)
    else:
        index = pd.read_csv(os.path.join(corrmat_name + ".index_names.txt"),sep="\t",header=None, index_col=None)
        index = np.array(index[0])
        length=len(index)
    
    start_time = time.time()
    print("PRED_CHROM: Starting predicting best chromosome for each SNP based on LD")

    if splitlength < length:   
        parts_y = math.ceil(length/splitlength) #calculate how many parts data should be divided if max snps in a single part is 10000  
        splits_y = list(np.array_split(range(length),parts_y)) #get ranges of splits
        for splitnr, split in enumerate(splits_y):
            corrset = rw.openFile(os.path.join(str(corrmat_name) + ".h5"), marker_column, marker_pos_data ,NA_val,chunk=split)
            index = pd.read_csv(os.path.join(corrmat_name + ".index_names.txt"),sep="\t",header=None, index_col=None)
            index = index[0].tolist()
            colnames = pd.read_csv(os.path.join(corrmat_name + ".column_names.txt"),sep="\t",header=None, index_col=None)
            colnames = np.array(colnames[0])
            corrset.index = index
            corrset.columns = colnames[split]
            corrset = pd.concat([marker_pos_data,corrset],axis=1, join="inner")
            print("Reading LD matrix for chunk ",splitnr," done.")
            # print(corrset.head())


            # ↓ Required to append to prevent bug & crash
            extra_columns.append(position_col_name)

            corrset, statistics_1 = func.calc_subchrom(marker_pos_data, 
                                                    marker_pos_data2.loc[colnames[split],:],
                                                    corrset, 
                                                    chr_col_name, 
                                                    pred_chr_col_name,
                                                    secondchoice_statistic_col_name,
                                                    suffix_count_values,
                                                    remove_columns,
                                                    subgenomemax_name,
                                                    subgenomemin_name,
                                                    cnt_subgenomemax_name,
                                                    cnt_subgenomemin_name,
                                                    extra_columns,
                                                    sumchrom_calc_method,
                                                    yadav,
                                                    marker_column,chunk=split)
            if splitnr == 0:
                statistics = statistics_1
            else:
                statistics = pd.concat([statistics, statistics_1])       
        
        #add markers that are not in corrset (due to filtering like MAF, or missing data) to statistics
        for marker in marker_pos_data2.index:
            if marker not in colnames:
                statistics.loc[marker, chr_col_name] = marker_pos_data2.loc[marker, chr_col_name]
                statistics.loc[marker, position_col_name] = marker_pos_data2.loc[marker, position_col_name]        

        statistics.sort_values(by=[chr_col_name, position_col_name, marker_column], inplace=True)
        rw.write_statmat_csv(statistics, statfile_name, decimalnums)
    else:
        corrset = rw.openFile(os.path.join(str(corrmat_name) + ".h5"), marker_column, marker_pos_data ,NA_val,chunk=None).to_numpy()
        if de_novo_or_anchor == "de_novo":
            np.fill_diagonal(corrset, np.nan) 
        corrset = pd.DataFrame(corrset, columns=None, index=None)
        index = pd.read_csv(os.path.join(corrmat_name + ".index_names.txt"),sep="\t",header=None, index_col=None)
        index = index[0].tolist()
        colnames = pd.read_csv(os.path.join(corrmat_name + ".column_names.txt"),sep="\t",header=None, index_col=None)
        colnames = colnames[0].tolist()
        corrset.index = index
        corrset.columns = colnames
        corrset = pd.concat([marker_pos_data,corrset],axis=1, join="inner")
        print("Reading LD matrix done.")
        extra_columns.append(position_col_name)
        
        corrset, statistics = func.calc_subchrom(marker_pos_data, 
                                                marker_pos_data2,
                                                corrset, 
                                                chr_col_name, 
                                                pred_chr_col_name,
                                                secondchoice_statistic_col_name,
                                                suffix_count_values,
                                                remove_columns,
                                                subgenomemax_name,
                                                subgenomemin_name,
                                                cnt_subgenomemax_name,
                                                cnt_subgenomemin_name,
                                                extra_columns,
                                                sumchrom_calc_method,
                                                yadav,
                                                marker_column,chunk=None)
        

        rw.write_statmat_csv(statistics, statfile_name, decimalnums)

    print("PRED_CHROM: Chromosome predicted for each SNP based on LD in", round((time.time() - start_time)/60,1), "minutes.")

    if run_fullpipe:
        try:
            if sys.argv[1] and sys.argv[1].endswith(".yaml"):
                os.system("python3 calculate_position.py " + sys.argv[1])
        except IndexError:
            os.system("python3 calculate_position.py")

main()
