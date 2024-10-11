import yaml
import background.functions_Mainfile as func
import background.read_write as rw
import pandas as pd
import sys
import numpy as np

def get_number_changed_chroms(summary, statistics, 
                        prev_chr_col_name, est_chr_col_name):
    """Generate a summary of how many values changed based on two columns.
    Requires an empty pandas dataframe and a statistics file with the two
    given column names. Also indexes summary file.
        

    Args:
        summary (pd df): Empty pandas dataframe
        statistics (pd df): Statistics file with at least prev & est columns.
        prev_chr_col_name (str): Name of the baseline column in statistics.
        est_chr_col_name (str): Name of the column baseline changed into.

    Returns:
        pd df:  Inedexed dataframe with baseline count of unique values, 
                counts of what baseline changed into and 
                a column of predicted count - baseline.
    """

    # Count total number of times a chromosome occurs.

    chromlist = sorted(set(statistics.loc[:, est_chr_col_name].astype(str)))
    summary["N_Prev"] = statistics[prev_chr_col_name].astype('string').value_counts().sort_index()
    temp = pd.DataFrame(statistics[est_chr_col_name].value_counts())
    summary = summary.join(temp).rename(columns={est_chr_col_name: "N_Pred"})
    summary["N_Diff"] = summary["N_Pred"] - summary["N_Prev"]

    # Compute values for chrom that occurs in chromlist, but not in statistics[prev_chr_col_name]
    missing_chroms = set(chromlist) - set(summary.index)
    for chrom in missing_chroms:
        summary.loc[chrom] = [0, statistics[est_chr_col_name].value_counts().get(chrom, 0), statistics[est_chr_col_name].value_counts().get(chrom, 0)]
    summary = summary.sort_index()
    # print(summary)
    return summary

def create_confusion_matrix(summary, statistics, 
                    prev_chr_col_name, est_chr_col_name,
                    suffix_confusion_matrix):
    """Creates a confusion matrix based on the index of the summary file. 
    Requires statistics file with prev_chr and pred_chr column names. 
    Requires indexed summary file.

    Args:
        summary (pd df): Indexed file with indexes of chrom appearing in statistics.
        statistics (pd df): Statistics file with at least prev & est columns.
        prev_chr_col_name (str): Name of the baseline column in statistics.
        est_chr_col_name (str): Name of the column baseline changed into.
        suffix_confusion_matrix (str): Suffix to be added 
                                        for cols in confusion_matrix

    Returns:
        pd df: summary file with added confusion matrix
    """    
    # Create confusion matrix: Col became X

    for chrom in summary.index:
        temp = statistics.loc[statistics[prev_chr_col_name] == chrom]
        summary[chrom + suffix_confusion_matrix] = temp[est_chr_col_name].value_counts(dropna=False).sort_index()
        #print(summary[chrom + suffix_confusion_matrix])
    return summary

def append_sum_avg_to_list_from_df(df, colname,
                                    avglist, totlist):
    """Takes pandas dataframe and column name and two lists. 
    Performs sum and total calculations, adds to lists.

    Args:
        df (pd df): Dataframe containing numeric colname.
        colname (str): String correspnds to column name in df.
        avglist (List): List with average numbers of colname.
        totlist (List): List with total numbers of colname.

    Returns:
        (List): List with average numbers of colname.
        (List): List with total numbers of colname.
    """
    # Get averge and total of absolute difference columns
    totlist.append(df[colname].sum())    
    try: 
        avglist.append(round(df[colname].mean(),2))
    except ValueError: 
        avglist.append(0)
    return avglist, totlist
    

def get_avg_total_per_chrom(summary, statistics, 
                                abs_diff_pos_name,
                                secondchoice_stat_name, 
                                yadav_data_col_name,
                                pred_chr_col_name,
                                avg_abs_diff,
                                tot_abs_diff,
                                avg_yadav_points,
                                tot_yadav_points,
                                avg_secondbest_choice,
                                tot_secondbest_choice):
    """Gets average and total numbers of three given 
    numeric columns. Requires indexed summary file
    and complete statistics file. 

    Args:
        summary (pd df): Indexed pd df. Index corresponds to 
                            chroms in statistics.
        statistics (pd df): Statistics file three given columns
                            containing numeric data.
        abs_diff_pos_name (str): Name of column in statistics file.
        secondchoice_stat_name (str): Name of column in statistics file.
        yadav_data_col_name (str): Name of column in statistics file.
        pred_chr_col_name (str): Name of column in statistics file.
        avg_abs_diff (str): Name of abs_col avg col in statistics.
        tot_abs_diff (str): Name of abs_col tot col in statistics.
        avg_yadav_points (str): Name of yadav avg col in statistics.
        tot_yadav_points (str): Name of yadav tot col in statistics.
        avg_secondbest_choice(str): Name of first/secondbest avg col in statistics.
        tot_secondbest_choice(str): Name of first/secondbest tot col in statistics.

    Returns:
        pd df: Summary file with avg and totals of three input columns. 
    """
    # Create empty lists
    avg_abs_diff_list = []
    tot_abs_diff_list = []
    avg_yadav_points_list = []
    tot_yadav_points_list = []
    avg_secondbest_choice_list = []
    tot_secondbest_choice_list = []
    # Loop through the different possible chroms
    for chrom in summary.index:
        temp = statistics.loc[statistics[pred_chr_col_name] == chrom]
        # Fill lists
        avg_abs_diff_list, tot_abs_diff_list = append_sum_avg_to_list_from_df(
                        temp, abs_diff_pos_name,
                        avg_abs_diff_list, tot_abs_diff_list
                        )

        avg_yadav_points_list, tot_yadav_points_list = append_sum_avg_to_list_from_df(
                            temp, yadav_data_col_name,
                            avg_yadav_points_list, tot_yadav_points_list
                            )

        avg_secondbest_choice_list, tot_secondbest_choice_list = append_sum_avg_to_list_from_df(
                                temp, secondchoice_stat_name,
                                avg_secondbest_choice_list, tot_secondbest_choice_list
        )
    
    # Add summary to statistics
    summary[avg_abs_diff] = avg_abs_diff_list
    summary[tot_abs_diff] = tot_abs_diff_list
    summary[avg_yadav_points] = avg_yadav_points_list
    summary[tot_yadav_points] = tot_yadav_points_list
    summary[avg_secondbest_choice] = avg_secondbest_choice_list
    summary[tot_secondbest_choice] = tot_secondbest_choice_list
    return summary


def main():
    try:
        if sys.argv[1] and sys.argv[1].endswith(".yaml"):
            with open(sys.argv[1], 'r') as configfile:
                config = yaml.safe_load(configfile)
    except IndexError:
        with open('config.yaml', 'r') as configfile:
            config = yaml.safe_load(configfile)
    run_ID =  config["run_ID"]
    data_type = config["data_type"]
    sumchrom_calc_method = config["sumchrom_calc_method"]
    run_ID =  config["run_ID"]

    statfile = config["statfile_name_with_position"].replace("{{ data_type }}", data_type).replace("{{ sumchrom_calc_method }}", sumchrom_calc_method).replace("{{ run_ID }}", run_ID)+".csv"
    marker_column = config["marker_column"]
    prev_chr_col_name = config["chr_col_name"]
    pred_chr_col_name = config["pred_chr_col_name"]
    abs_diff_pos_name = config["abs_diff_pos_name"]
    secondchoice_stat_name = config["secondchoice_statistic_col_name"]
    yadav_data_col_name = config["datapoints_used_for_yadav_col"]
    suffix_confusion_matrix = config["suffix_confusion_matrix"]

    avg_abs_diff = config["avg_abs_diff"]
    tot_abs_diff = config["tot_abs_diff"]
    avg_yadav_points = config["avg_yadav_points"]
    tot_yadav_points = config["tot_yadav_points"]
    avg_secondbest_choice = config["avg_secondbest_choice"]
    tot_secondbest_choice = config["tot_secondbest_choice"]   
    summaryfile_name = config["summaryfile_name"].replace("{{ data_type }}", data_type).replace("{{ sumchrom_calc_method }}", sumchrom_calc_method).replace("{{ run_ID }}", run_ID) 
    run_fullpipe = config["run_fullpipe"]
    NA_val=config["NA_val"]
  
    # Read statistics file
    statistics = rw.openFile(statfile, marker_column,None,NA_val)
    statistics[prev_chr_col_name] = statistics[prev_chr_col_name].fillna('Unknown')    # Get empty pandas dataframe
    summary = pd.DataFrame()

    summary = get_number_changed_chroms(summary, statistics, 
                                prev_chr_col_name, pred_chr_col_name)
    summary = get_avg_total_per_chrom(summary, statistics, 
                                abs_diff_pos_name,
                                secondchoice_stat_name, 
                                yadav_data_col_name,
                                pred_chr_col_name,
                                avg_abs_diff,
                                tot_abs_diff,
                                avg_yadav_points,
                                tot_yadav_points,
                                avg_secondbest_choice,
                                tot_secondbest_choice)

    summary = create_confusion_matrix(summary, statistics, 
                                prev_chr_col_name, pred_chr_col_name,
                                suffix_confusion_matrix)

    #print("Writing summary...")
    summary.to_csv(summaryfile_name+".csv")
    print("Wrote summary statistics file as ",summaryfile_name,".csv", sep="")
    
    if run_fullpipe:
            print("FINISHED complete analysis.") #TODO add time?

main()
