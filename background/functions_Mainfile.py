import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import time
import numpy as np
import background.read_write as rw
import concurrent.futures
import math
import h5py
import os

def af_to_maf(value):
    """
    Convert an allele frequency to a minor allele frequency.
    Args:
        value (float): The allele frequency.
    Returns:
        float: The minor allele frequency.
    """
    if value > 0.5:
        return 1 - value
    return value

def np_pearson_cor(x, y):
    """
    Calculate the Pearson correlation coefficient between two arrays. 
    Written on Apr 23, 2020 by Philip Montgomery. Taken from https://cancerdatascience.org/blog/posts/pearson-correlation/ accessed on Oct 15, 2023.

    Parameters:
    x (numpy.ndarray): The first array.
    y (numpy.ndarray): The second array.

    Returns:
    numpy.ndarray: The Pearson correlation coefficient between x and y.
    """
    xv = x - x.mean(axis=0)
    yv = y - y.mean(axis=0)
    xvss = (xv * xv).sum(axis=0)
    yvss = (yv * yv).sum(axis=0)
    result = np.matmul(xv.transpose(), yv) / np.sqrt(np.outer(xvss, yvss))
    # bound the values to -1 to 1 in the event of precision issues
    return np.maximum(np.minimum(result, 1.0), -1.0)

def cor_to_LD_round(cor, decimalnums=3):
    """
    Converts a correlation coefficient to an LD value and rounds it to the specified number of decimal places.
    Args:
        cor (float): The correlation coefficient to convert.
        decimalnums (int): The number of decimal places to round the LD value to. Default is 3.
    Returns:
        float: The LD value, rounded to the specified number of decimal places.
    """
    LD = cor ** 2
    return LD.round(decimalnums)

def corrmatrix_calc(data, data2, datacolumns, marker_pos_data,corrmat_name, 
                    marker_column, decimalnums, remove_values, maf_threshold,splitlength,
                    allowed_nan_perc = 100, corrmethod = "leave"):
    """
    Calculates a correlation matrix using a pandas dataframe and set parameters.

    Args:
        data (pandas.DataFrame): A dataframe containing markers, positions, and sequenced SNPs of anchorset.
        data2 (pandas.DataFrame): A dataframe containing markers, positions, and sequenced SNPs of SNPs that need prediction.
        datacolumns (int): The number of columns in the dataframe containing marker, position, and chromosome information.
        marker_pos_data (pandas.DataFrame): A dataframe containing marker, position, and chromosome information.
        marker_column (str): The name of the column containing marker information.
        decimalnums (int): The number of decimal places to round the correlation matrix values to.
        remove_values (float): Values in the dataframe below this set value are removed and not used in further calculations.
        maf_threshold (float): The threshold below which SNPs are discarded in the analysis.
        corrmat_name (str): The name of the correlation matrix file to be created.
        allowed_nan_perc (int, optional): The allowed percentage of NaN in markers. Defaults to 100.
        corrmethod (str, optional): Determines which correlation method to use for computing corrmat. Can be "leave" or "fill". Defaults to "leave".

    Returns:
        corrset (pandas.DataFrame): A matrix containing LD values among all SNPs, markers, positions, and chromosomal positions.
    """

    start_time = time.time()
    print("LD: Starting computing LD matrix. Caution: this can take a while.")    
    # Create transposed dataset, excluding first data columns and remove snps with low maf
    numericTp = data.T.iloc[datacolumns:, 0:].apply(pd.to_numeric, errors='coerce')
    # print(numericTp.head())
    # print(numericTp2.head())
    print("Raw data succesfully loaded. Shape:", numericTp.shape)
    if corrmethod == "fill":    
        print("LD: Computing LD matrix using fill method.")	    
        #remove low MAF snps from analysis
        allele_count = np.sum(numericTp, axis=0)
        notnan_count = np.sum(~np.isnan(numericTp), axis=0)
        af = allele_count / (notnan_count*2)
        maf = [af_to_maf(value) for value in af]
        maflist = dict(zip(numericTp.columns, maf))
        low_maf_snps = {snp for snp, value in maflist.items() if value < maf_threshold}
        numericTp = numericTp.drop(low_maf_snps,axis=1)
        print("LD: Removed SNPs with MAF lower than", maf_threshold)
        
        # ↓ Get columns where NaN vals per sample is less than var
        # Save the columns with less NaN than cutoff
        less_than_x_missing_vals = (numericTp.isnull().sum() * 100 / len(numericTp)<allowed_nan_perc) 
        numericTp = numericTp.loc[:, less_than_x_missing_vals]
        print("LD: SNPs with more than", allowed_nan_perc, "% missing values removed.")
        # ↓ Fill NaN with median of every column column for numpy calculation
        numericTp = numericTp.fillna(numericTp.median())
        print("LD: Missing values imputed by the median of the SNP.")
        x=numericTp.to_numpy()
        xnames=numericTp.columns
        print(x.shape)
        numericTp=None #clean up mem

        if data2 is not None:
            numericTp2 = data2.T.iloc[datacolumns:, 0:].apply(pd.to_numeric, errors='coerce')
            allele_count = np.sum(numericTp2, axis=0)
            notnan_count = np.sum(~np.isnan(numericTp2), axis=0)
            af = allele_count / (notnan_count*2)
            maf = [af_to_maf(value) for value in af]
            maflist = dict(zip(numericTp2.columns, maf))
            low_maf_snps = {snp for snp, value in maflist.items() if value < maf_threshold}
            numericTp2 = numericTp2.drop(low_maf_snps,axis=1)
            print("LD: Removed SNPs with MAF lower than", maf_threshold, "from checkdata.")
            
            # ↓ Get columns where NaN vals per sample is less than var
            # Save the columns with less NaN than cutoff
            less_than_x_missing_vals = (numericTp2.isnull().sum() * 100 / len(numericTp2)<allowed_nan_perc) 
            numericTp2 = numericTp2.loc[:, less_than_x_missing_vals]
            print("LD: SNPs with more than", allowed_nan_perc, "% missing values removed from checkdata.")
            print(numericTp2.head())
            # ↓ Fill NaN with median of every column column for numpy calculation
            numericTp2.fillna(numericTp2.median(), inplace=True)
            print("LD: Missing values imputed by the median of the SNP in checkdata.")
        
            # print(numericTp2.head())
            y=numericTp2.to_numpy()
            ynames=numericTp2.columns
        else: 
            y=x
            ynames=xnames

        # parts_x = math.ceil(len(xnames)/10000) #calculate how many parts data should be divided if max snps in a single part is 10000  
        # splits_x = list(np.array_split(range(len(xnames)),parts_x)) #get ranges of splits
        # print(splitsx)
        splitlength = int(splitlength)
        parts_y = math.ceil(len(ynames)/splitlength) #calculate how many parts data should be divided if max snps in a single part is 10000  
        splits_y = list(np.array_split(range(len(ynames)),parts_y)) #get ranges of splits
        # print(splits_y)

        corrset=None
        os.makedirs(os.path.dirname(corrmat_name), exist_ok=True) #mkdir if not exist
        nan_counter = 0
        with h5py.File(os.path.join(corrmat_name + ".h5"), 'a') as hf:
            for splitnr, split in enumerate(splits_y):
                # print(split)
                # print("Decimalnums:", decimalnums)
                corrset_split = cor_to_LD_round(np_pearson_cor(x,y[...,split])) ### subset x by created splits (splitsx)
                corrset_split = pd.DataFrame(corrset_split, columns=ynames[split],index=xnames,dtype='float16') #transform to pd dataframe
                corrset_split.sort_index()
                # np.fill_diagonal(corrset_split.values, np.nan) 
                corrset_split = (corrset_split.mask(corrset_split < remove_values)) #.astype("float16")#.astype(pd.SparseDtype("float16", np.nan))
                
                ## remove LD values of variants with itself. 
                ## This for loop can be optimized.
                for i in range(len(corrset_split.index)):
                    if corrset_split.index[i] in corrset_split.columns:
                        corrset_split.loc[corrset_split.index[i], corrset_split.index[i]] = np.nan
                        nan_counter += 1
                ## create or append hdf5 file (LD matrix) with corrset_split        
                if splitnr == 0:
                    hf.create_dataset('corrset', data=corrset_split.to_numpy(), compression="gzip", chunks=True, maxshape=(None,None))
                    ynames2=ynames[split]
                    ynames2=ynames2.values.tolist()
                else:
                    hf["corrset"].resize((hf["corrset"].shape[1] + corrset_split.shape[1]), axis = 1)
                    hf["corrset"][:,-corrset_split.shape[1]:] = corrset_split
                    ynames2 = ynames2 + ynames[split].values.tolist()
                print("Iteration {} and 'data' chunk has shape:{}".format(splitnr,hf['corrset'].shape))
            #hf.close()
            index_names=pd.Series(xnames)
            ## storing seperate files for column and index names because hard to store strings in hdf5
            index_names.to_csv(os.path.join(corrmat_name + ".index_names.txt"),sep="\t",header=None,index=False)
            ynames2=pd.Series(ynames2)
            ynames2.to_csv(os.path.join(corrmat_name + ".column_names.txt"),sep="\t",header=None,index=False)
            print(nan_counter, " values set to NaN")
            print("Corrset calculation done")
            # corrset = np.array(hf['corrset']) # `data` is now an ndarray.
            # if data2 is None:
            #     print("Need to remove diagonal LD values from LD matrix as these are LD values of markers with itself. \nThis is not optimised for memory so may result in a fatal error.")
            #     corrset = np.array(hf['corrset']) # `data` is now an ndarray.
            #     np.fill_diagonal(corrset, np.nan) #remove LD values of markers with itself when "de_novo" option is used. With "anhcor" option these are not present.
            #     #replace hf['corrset'] with corrset with removed diagonal values
            #     hf.create_dataset('corrset', data=corrset, compression="gzip", chunks=True, maxshape=(None,None))
            # corrset = pd.DataFrame(corrset, columns=ynames,index=index_names,dtype='float16') #transform to pd dataframe
            # print(corrset.head())
            # corrset.to_csv(os.path.join(corrmat_name + ".corrset.csv"),sep=",",header=None,index=False)
        hf.close()
    else: 
        print("LD: Computing LD matrix using leave method. \nCAUTION: this is not optimized for memory consumption and may result in a fatal error.")
        # Create cor matrix, square values, round values to decimal number
        corrset = (numericTp.corr(method="pearson")**2).round(decimals=decimalnums) 
    # Convert object dtype to fixed length string dtype
    # corrset = pd.concat([marker_pos_data,corrset],axis=1, join="inner")
    print("LD: LD-matrix (r-squared) computed in ", round((time.time() - start_time)/60,1), "minutes.")
    return corrset


def calc_subchrom(marker_pos_data, marker_pos_data2,corrset,  
                    chr_col_name, pred_chr_col_name,
                    second_choice_statistic_col_name,
                    suffix_count_values, 
                    remove_columns, subgenomemax_name,
                    subgenomemin_name, cnt_subgenomemax_name,
                    cnt_subgenomemin_name,
                    extra_columns,
                    sumchrom_calc_method, yadav, marker_column, chunk=None): 
    """
    Calculates & creates pandas dataframe containing information about the    
    calculated correlation matrix and a pandas dataframe about 
    statistics of how the new subchromosome location is calculated.

    Args:
        marker_pos_data (pandas df): markers, positions & chrom data,
        corrset (pandas dataframe) : Matrix containing correlation scores, 
        chr_col_name       (string): Name of the previous predicted chromosome,
        pred_chr_col_name  (string): Name of the to be created predicted chrom col name,
        second_choice_statistic_col_name:  
                           (string): Contains the name for the col to be used in 
                                     writing the second best choice divided by the best,
        suffix_count_values   (str): Suffix to be used for pasting to count values,
        remove_columns       (list): List of columns to remove from the dataset.
                               Columns contain chromosome values to be removed (such as dummy chromosomes).
        subgenomemax          (str): Name of subgenome Maxvalue column.
        subgenomemin          (str): Name of subgenome Minvalue column.
        cnt_subgenomemax      (str): Name of the subgenome max count column.
        cnt_subgenomemin      (str): Name of the subgenome min count column.
        extra_columns        (list): List of extra column names plus position column name.
                                     Used in reducing dataframe and preventing bug. 
        sumchrom_calc_method(_str_): Method to be used to calculate subchrom.
                                     Must be "yadav" or "all". "yadav" uses top two 
                                     correlations in corrset, "all" uses all correlations
                                     in corrset. 

    Returns:
        corrset (pandas dataframe): Matrix containing correlation scores, markers,
                                    position, chrom position & est chrom position.
        statistics (pandas dataframe): Matrix with position, markers, est position 
                                       and total calculated scores of the results.
                                       
    """


    
    
    ########################
    ##### YADAV METHOD ##### calculate yadav to compare with other methods
    ########################
    
    if yadav:
        statistics_y = corrset
        statistics_y = statistics_y.drop(columns = extra_columns, inplace = False, errors="ignore")
        #print(corrset.head())
        #Remove all correlations (in rows) of SNPs on dummy chromosomes (from "remove_columns") to keep them from influencing the analysis
        snps_to_update = statistics_y[chr_col_name].isin(remove_columns)
        statistics_y.loc[snps_to_update, statistics_y.columns != chr_col_name] = np.nan

        #print(statistics.columns)
        yadav_pred_list = []
        for column in statistics_y.columns: 
            if column == chr_col_name:
                continue
            else:
                tmp = pd.to_numeric(statistics_y[column],errors='coerce')
                tmp = tmp.nlargest(2)
                #print(tmp)
                if tmp.isnull().values.any() == False:
                    #print("values found!")
                    if statistics_y.at[tmp.index[0], chr_col_name] == statistics_y.at[tmp.index[1], chr_col_name]:
                        yadav_pred_list.append(statistics_y.at[tmp.index[0], chr_col_name])
                    else:
                        yadav_pred_list.append(pd.NA)
                else:
                    yadav_pred_list.append(pd.NA)

        # Create a dataframe like 'pred_chr_col' dataframe generated further down the line
        yadav_pred = pd.DataFrame()
        yadav_pred[""] = yadav_pred_list
        # print(statistics_y.head())
        yadav_pred.index = statistics_y.columns[1:]
        statistics_y=None
        # print(yadav_pred.head())
    # statistics = corrset.drop(columns = extra_columns, inplace = False, errors="ignore").groupby(chr_col_name).sum().T
    
    
    # Calculating the subchromosome with "top10" or "sum_all" without breaking code
    
    if sumchrom_calc_method == "top10":
        statistics = corrset.drop(columns = extra_columns, inplace = False, errors="ignore")
        
        snps_to_update = statistics[chr_col_name].isin(remove_columns)
        statistics.loc[snps_to_update, statistics.columns != chr_col_name] = np.nan
        
        pred_list = []


        for column in statistics.columns: 
            if column == chr_col_name:
                continue
            else:
                correlations = statistics[column][~np.isnan(statistics[column])]
                if correlations.size >= 10:
                    tmp = correlations.nlargest(10)
                elif correlations.size > 0 and correlations.size < 10:
                    tmp = correlations.nlargest(correlations.size)
                else:
                    tmp = pd.Series(np.nan)
                    #print("No correlations found for", column)
                
                if pd.isna(tmp).any() == False:
                    top10_list = []
                    for ind in range(0, tmp.size):                        
                        top10_list.append(statistics.at[tmp.index[ind], chr_col_name])

                    counts = {}
                    for chrom in top10_list:
                        if chrom in counts:
                            counts[chrom] += 1
                        else:
                            counts[chrom] = 1

                    most_common_chrom = max(counts, key=counts.get)
                    count = counts[most_common_chrom]
                    #print("Chr_pred:", most_common_chrom,". ",count," occurrences.")
                    if count > 0:
                        pred_list.append(most_common_chrom)
                    else:
                        pred_list.append(pd.NA)
                else:
                    pred_list.append(pd.NA)

        pred = pd.DataFrame()
        pred[""] = pred_list
        pred.index = statistics.columns[1:]

        statistics = corrset.drop(columns = extra_columns, inplace = False, errors="ignore").groupby(chr_col_name).sum().T    
    
    elif sumchrom_calc_method == "sum_all_squared":

        #same as "sum_all", but first square (**2) all LD values to minimize the impact of many low LD SNPs
        statistics = (corrset.drop(columns=extra_columns, inplace=False, errors="ignore")
                                .apply(lambda val: val ** 2 if pd.api.types.is_numeric_dtype(val) else val)
                                .groupby(chr_col_name)
                                .sum(numeric_only=True)
                                .T)

    elif sumchrom_calc_method == "sum_all":
        statistics = (corrset.drop(columns = extra_columns, inplace = False, errors="ignore")
                                .groupby(chr_col_name)
                                .sum(numeric_only=True)
                                .T)
    else:
        print("ERROR: Please specify a valid sumchrom_calc_method in the config file.")
    
    
    # Number of values pointing to a specific subchrom
    counts = corrset.drop(columns = extra_columns, inplace = False, errors="ignore").groupby(chr_col_name).count().T

    # Increase readability by removing all 0 values
    statistics = statistics.mask(statistics <= 0.0)
    counts = counts.mask(counts <= 0.0).astype("Int64") #Test: Int16 → Memory optimization
    
    # Sort chrom columns by name
    statistics = statistics.reindex(sorted(statistics.columns), axis=1)
    counts = counts.reindex(sorted(counts.columns), axis=1)

    # TODO - might be done already - Remove possibility of assigning to -1 and other unwanted columns
    statistics = statistics.drop(columns=remove_columns, axis= 1, errors="ignore")
    counts = counts.drop(columns=remove_columns, axis = 1, errors="ignore")
    
    # add extra count column to avoid error with fetching of the subgenome (np.nan is the name of the column)
    counts[np.nan] = np.nan

    # Create temp list of highest values
    temp = statistics.apply(lambda row: row.nlargest(2).values[-1],axis=1)

    # the two subgenomes that perform best
    subgenomemax = statistics.apply(lambda row: row.nlargest(2).idxmax(), axis=1)
    ## the second best subgenome, only if there is a second best subgenome
    subgenomemin = statistics.apply(lambda row: np.nan if row.nlargest(2).idxmax() == row.nlargest(2).idxmin() else row.nlargest(2).idxmin(), axis=1)

    cnt_subgenomemax, cnt_subgenomemin = [], []
    # loop through the subgenomes and fetches the subgenome position count 
    for index in range(len(subgenomemax)):
        maxs = subgenomemax[index]
        mins = subgenomemin[index]
        # print(index, maxs, mins, counts.iloc[index][maxs], counts.iloc[index][mins])
        cnt_subgenomemax.append(counts.iloc[index][maxs])
        if pd.isna(mins):
            # Handle NaN case, e.g., append a default value like 0 or np.nan
            cnt_subgenomemin.append(np.nan)  # or any other default value
        else:
            cnt_subgenomemin.append(counts.iloc[index][mins])
    
    # Add suffix to count columns after work getting best / secondbest cols
    counts = counts.add_suffix(suffix_count_values)

    # Insert estimated chromosome using Yadav or all method
    #if sumchrom_calc_method == "yadav":
        
    Yadav_chr = yadav_pred
    if sumchrom_calc_method == "top10":
        pred_chr_col = pred
    elif sumchrom_calc_method == "sum_all_squared":
        # Put all calculations seperate from insert, then insert in statistics file
        pred_chr_col = statistics.idxmax(axis=1)
    elif sumchrom_calc_method == "sum_all":
        # Put all calculations seperate from insert, then insert in statistics file
        pred_chr_col = statistics.idxmax(axis=1)
    else:
        print("ERROR: Please specify a valid sumchrom_calc_method in the config file.")

    # second_choice_statistic = statistics.max(axis=1)/temp
    second_choice_statistic = 1-(temp.fillna(0)/statistics.max(axis=1)) #fill np.nan with 0 to get an certainty of 1 if there is no second chromosome option. Now only SNPs that do not have a predicted chromosome result in np.nan for certainty. 

    # Insert secondbest choice statistic by dividing the best choice by the secondbest choice
    statistics.insert(0, pred_chr_col_name, pred_chr_col)  
    statistics.insert(1, "Yadav_chr", Yadav_chr)
    statistics.insert(2, second_choice_statistic_col_name, second_choice_statistic)

    # store subgenome counts in statistics file
    statistics.insert(3, subgenomemax_name, subgenomemax)  
    statistics.insert(4, subgenomemin_name, subgenomemin)  
    statistics.insert(5, cnt_subgenomemax_name, cnt_subgenomemax)  
    statistics.insert(6, cnt_subgenomemin_name, cnt_subgenomemin)  


    # Add Chrom, Pos + New_chr back using 'left join' 
    if marker_pos_data2 is not None:
        statistics = marker_pos_data2.join(statistics)
    else:
        statistics = marker_pos_data.join(statistics)

    # Merge sum and count data in one df
    statistics = statistics.join(counts)
    statistics.index.name = marker_column #add marker name as index name
    return corrset, statistics




def check_difference(data, chr_col_name, pred_chr_col_name):
    """Check difference between two predicted markers. Returns dataset
    with only changed markers and dataset with only unchanged markers. 
    
    Args: 
        Data            (pd df): Dataframe with at least two given columns 
                                 chr_col_name & pred_chr_col_name 
                                 containing similar data.
        chr_col_name      (str): String containing the previous known 
                                 chrom column name.
        pred_chr_col_name (str): String containing 
                                 the predicted column name.
        
    Returns:         
        difdata (pd df): Dataframe containing rows of which 
                         predicted not equal to previous known.
        notdif  (pd df): Dataframe containing rows of which 
                         predicted is equal to previous known.
    """
    # print("difdata")
    # print(data.head())

    # Get the markers for which the predicted chromosome was changed
    difdata = (data.loc[(data[chr_col_name] != data[pred_chr_col_name])])
    # print(difdata.head())
    # Remove markers for which no correlations were found and were kept NaN
    difdata = difdata.dropna(subset=pred_chr_col_name)
 
    # Get markers for which the prediction was not changed
    notdif = data.loc[~data.index.isin(difdata.index),:]

    return difdata, notdif

def calculate_positions(data, statistics, position_col_name, marker_pos_data2,
                        est_pos_col_name ,datapoints_used_for_yadav_col,
                        calculate_all_or_different, difdata, notdif, abs_diff_pos_name):
    """
    Function to estimate positions using unchanged datapoints. 

    Args:
        data                 (pd df): Dataframe containing information about
                                      correlations between different markers.
        statistics           (pd df): Dataframe containing 
                                      position_col_name, est_pos_col_name,
                                      which are used to calculate positions.              
        position_col_name 
                               (str): Name of previous known position column.
        est_pos_col_name 
                               (str): Name of predicted position column.
                                      This column is going to be new.
        datapoints_used_for_yadav_col 
                               (str): Name of the column describing 
                                	  datapoints used for yadav position 
                                      calculation technique. 
        calculate_all_or_different 
                               (str): Variable used for describing if 
                                      all positions need to be calculated
                                      or only the different previously
                                      predicted to be different. 
        difdata              (pd df): Dataframe containing rows of which 
                                      predicted not equal to previous known.
        notdif               (pd df): Dataframe containing rows of which 
                                      predicted is equal to previous known.

    Returns:
        data       (pd df): Dataframe containing information about
                            correlations between different markers
                            including estimated position column.
        statistics (pd df): Dataframe containing 
                            position_col_name, est_pos_col_name,
                            including added estimated position col.
    """        
    

    #print("Inserting columns into dataframes...")    
    # Insert column of estimated position into statistics & corrmatrix
    data.insert(loc=data.columns.get_loc(position_col_name)+1, 
                                         column=est_pos_col_name, 
                                         value=np.nan)

    statistics.insert(loc=statistics.columns.get_loc(position_col_name)+1, 
                                                     column=est_pos_col_name, 
                                                     value=np.nan)

    statistics.insert(loc=statistics.columns.get_loc(position_col_name)+2, 
                        column=datapoints_used_for_yadav_col, value=np.nan)


    # If-else structure to check if user wants to calculate all positions
    # or only the positions for a new chromosome. 
    if calculate_all_or_different == "all":
        calc_pos = statistics.index
    elif calculate_all_or_different == "different":
        calc_pos = difdata.index
    else: 
        calc_pos = difdata.index    

    # start_time = time.time()

    # print(notdif.head())
        # TODO this for- loop could be sped up, perhaps save values and insert at the same time?
    for marker in calc_pos:
        target_col = marker
        try: #TODO Analyse this except better - too big float16 num = (infinite values)
            filled_data = notdif[~notdif[target_col].isna()]
            #print(filled_data)
            #filled_data.replace([np.inf, -np.inf], 1, inplace=True) #replace infinite values that only seldom occur
        except KeyError:
            #print("KeyError has occured while attempting to retrieve unchanged data. Ignoring.")
            continue
        statistics[datapoints_used_for_yadav_col][target_col] = filled_data[target_col].count() 
        # ↑ Number of datapoints used for calculating position

        try:
                # Fill calculated sum of estimated position as an integer 
                                ##### YADAV METHOD #####
            est_pos = int(sum(filled_data[position_col_name] * filled_data[target_col])/
                                    (filled_data[[target_col]].sum()))

            data[est_pos_col_name][target_col] = est_pos
            statistics[est_pos_col_name][target_col] = est_pos

        except (ValueError):
            # ↑ Happens when not enough data is available
            data[est_pos_col_name][target_col] = np.nan
            statistics[est_pos_col_name][target_col] = np.nan
        except (OverflowError):
            #print("Overflow error; Occured when calculating below ↓:")
            #print(str(filled_data[position_col_name]) + "*" + str(filled_data[target_col])+
            #        "\n/" + str(filled_data[[target_col]].sum()))
            try:
                position = np.float64(filled_data[position_col_name])
                certainty = np.float64(filled_data[target_col])
                sum_acc = np.float64(filled_data[[target_col]].sum())
                est_pos = int((position * certainty)/sum_acc)

                data[est_pos_col_name][target_col] = est_pos
                statistics[est_pos_col_name][target_col] = est_pos
                #print("Overflow error fixed.")
            except OverflowError:
                print("Error: Could not fix the Overflow error for below marker ↓.")
                print(str(filled_data[position_col_name]) + "*" + str(filled_data[target_col])+
                    "\n/" + str(filled_data[[target_col]].sum()))
            except TypeError:
                print("Error: Could not fix the Overflow error for below marker ↓.")
                print(str(filled_data[position_col_name]) + "*" + str(filled_data[target_col])+
                    "\n/" + str(filled_data[[target_col]].sum()))

    #print("Added values to files using 'loop' function: ", time.time() - start_time, "seconds.")

    statistics[datapoints_used_for_yadav_col] = statistics[datapoints_used_for_yadav_col].convert_dtypes()
    statistics[est_pos_col_name] = statistics[est_pos_col_name].convert_dtypes()

    statistics.insert(loc=statistics.columns.get_loc(position_col_name)+2, 
                                            column=abs_diff_pos_name, 
                                            value=abs(statistics[position_col_name]-statistics[est_pos_col_name]))
    return data, statistics
