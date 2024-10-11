import pandas as pd
import sys


def get_marker_pos_data(data, datacolumns):
    """
    Takes pandas dataframe and returns the first set columns 
    containing marker, position and chromosome location.

    Args:
        data  (pandas df): df containing columns and rows of markers
        datacolumns (int): Number of marker - pos - chrom containing columns
 
    Returns:
        marker_pos_data (pandas df): markers, positions & chrom data
    """
    #TODO Create variable to take column names instead of column index
    # Marker Pos is first columns with product info
    marker_pos_data = data.iloc[0:, :datacolumns]
    return marker_pos_data

def drop_cols_containing_x(data, remove_str_in_columns):
    """
    Drops columns from a pandas dataframe containing a specific string.
    Careful with this function, as it uses a regex pattern to find
    and remove the columns. 

    Args:
        data (pandas dataframe): Pandas dataframe containing columns.
        remove_str_in_columns (List of strings) : Contains string
                                                  that is to be removed
                                                  if they are in column names.        
    """
    for str in remove_str_in_columns:
        data = data[data.columns.drop(list(data.filter(regex=str)))]
    return data

def merge_files(data, data2, 
                chr_col_name, 
                position_col_name,
                checkdata_chr_col, 
                checkdata_pos_col, 
                checkdata_extra_cols,
                marker_column,
                extra_columns,
                remove_columns):
    """
    Combines two pandas dataframes together by

    Args:
        data                        (pd df): Anchor file
        data2                       (pd df): Testing file
        chr_col_name                (_str_): Name of chrom col in anchor file
        position_col_name           (_str_): Name of pos col in anchor file
        checkdata_chr_col           (_str_): Name of chrom col in check file
        checkdata_pos_col            (_str_): Name of pos col in check file
        checkdata_extra_cols (_list of str_): Name(s) of extra cols 
                                              in check file
        marker_column               (_str_): Name of marker col in anchor file
        extra_columns        (_list of str_): Name(s) of extra cols 
                                              in anchor file
        remove_columns       (_list_of_str_): List with name(s) 
                                              of illegal chrom predictions

    Returns:
        Merged pandas dataframe with index marker, chrom col and pos col.
        Does not contain any other cols as they are removed in merging process.
    """

    # Drop extra cols; ease merging files
    data.drop(columns= extra_columns, inplace=True, errors="ignore")
    data2.drop(columns = checkdata_extra_cols, inplace=True, errors="ignore")
    
    # Rename checkfile col names to anchor col name, ease merge
    data2.rename(columns={checkdata_pos_col:position_col_name}, inplace = True)   
    data2.rename(columns={checkdata_chr_col:chr_col_name}, inplace = True)   
    data2.index.name = marker_column

    # Rename all chr to illegal end result
    data2.loc[:, chr_col_name] = remove_columns[0]

    # Use inner join to merge files; only keep columns with same name.
    joined = pd.concat([data, data2], join="inner")
    print(joined.shape[1]-2," overlapping genotypes.", sep="")
    return joined

def check_intersect(data, data2, remove_intersect=False):
    """
    Checks if there is overlap between the indexes of two pandas dataframes. 
    If there is overlap; writes overlap to a csv file and terminates script. 

    Args:
        data (_pd df_): _Pandas dataframe with at least a set index_
        data2 (_pd df_): _Pandas dataframe with at least a set index_
    """
    if set(data.index).intersection(set(data2.index)):
        if remove_intersect == False:
            print("\n\n!!! Fatal error !!!\nOverlap in index detected in both datasets.\n"+
                "Printing overlap to 'overlap_index.csv' ...")
            pd.DataFrame(set(data.index).intersection(set(data2.index))).rename(columns={0:"overlapping_marker_names"}).to_csv("overlap_index.csv", index=False)
            print("Done. Trying to mask LD values of variants with itself with NaN... \n\n")
            #sys.exit()
        else: 
            print("\n\n!!! Warning !!!\nOverlap in index detected in both datasets.\n"+
                "Printing overlap to 'overlap_index.csv' ...")
            pd.DataFrame(set(data.index).intersection(set(data2.index))).rename(columns={0:"overlapping_marker_names"}).to_csv("overlap_index.csv", index=False)
            print("Done. Removing overlapping values in checkfile...")
            
            remove_index = set(data.index).intersection(set(data2.index))
            data2 = data2.drop(index=list(remove_index))
            print("Done.")
            
    return data2
