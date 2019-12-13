from pybedtools import BedTool
import pybedtools
import pandas as pd
import uuid

"""left join based on coordinates and a name


This program can do the following:

1. specify a left table and a name to match

2. specify a right table and a name to match

        in the right table, you can have many coordinates
    
3. specify on the right table, which columns to be joined.


Method
-------


add the column (and the row) if any of the specified coordinates
 and name (in the right column) match to the ones on the left


Note
----

Below is an output from pybedtools, as you can see, colname doesn't really match

   chrom      start  end      name  score    strand  thickStart  thickEnd                itemRgb  blockCount
0  chr10  131809366  131809389  TRAC_site_1   .   -1       -1           .                      .          0
1  chr10  133305098  133305121  TRAC_site_1  chr10  133305098   133305121  TRAC_site_1  GTCAGGCCTCCGGATAACTGCGG    23
2  chr10    5627726 5627749  TRAC_site_1      .   -1       -1           .                      .          0
3  chr11  101192038  101192061  TRAC_site_1   .   -1       -1           .                      .          0
4  chr11  125280830  125280853  TRAC_site_1   .   -1       -1           .                      .          0


"""
def guess_sep(x):
    with open(x) as f:
        for line in f:
            tmp1 = len(line.strip().split(","))
            tmp2 = len(line.strip().split("\t"))
            # print (tmp1,tmp2)
            if tmp1 > tmp2:
                return ","
            if tmp2 > tmp1: 
                return "\t"
            else:
                print ("Can't determine the separator. Please input manually")
                exit()
                

def to_bed(df, col_list,output_filename):
    df[col_list].dropna().astype(int,errors='ignore').to_csv(output_filename,sep="\t",header=False,index=False)

def read_table(f):
    return pd.read_csv(f,sep=guess_sep(f))

def set_bed_index(df,col_list):
    df.index = df[col_list[0]]
    for c in col_list[1:]:
        df.index = df.index.astype(str)+"."+df[c].astype(str)

def merge_values(val1, val2):
    if val1 == val2:
        return val1
    elif val1 == "" or val1 == None:
        return val2
    elif val2 == "" or val2 == None:
        return val1
    else:
        return "%s:%s"%(val1,val2)

def merge_dict(d1,d2):
    return {key: merge_values(d1.get(key), d2.get(key)) for key in set(d1).union(d2)}


## Define input ##
left_table = "CHANGEseq_GUIDEseq_classes_validation_set_20191212.csv"
right_table = "GUIDEseq_matched_subset_20191212.csv"

## Define parameters ##
left_table_bed3 = ['chr','start','end']
left_table_match_name = "sample"


right_table_bed3_list = [['chr','Site_SubstitutionsOnly.Start','Site_SubstitutionsOnly.End'],['chr','Site_GapsAllowed.Start','Site_GapsAllowed.End']]
right_table_match_name = "sample"
right_table_columns_to_join = ['Site_SubstitutionsOnly.Sequence','Site_GapsAllowed.Sequence']

## read input ##
left = read_table(left_table)
right = read_table(right_table)

## prepare to start joining, get left bedtools obj for intersection ##
left_bed_file = "/tmp/%s"%(uuid.uuid4())
to_bed(left, left_table_bed3+[left_table_match_name],left_bed_file)
left_obj = BedTool(left_bed_file)

set_bed_index(left,left_table_bed3+[left_table_match_name])

## for each new column, get the content for each row ##
for columns_to_join in right_table_columns_to_join:
    left[columns_to_join] = left.index.tolist()
    my_replace_dict = {x:"" for x in left.index.tolist()}
    my_tmp_dict = {}
    for coordinates in right_table_bed3_list:
        tmp_out = "/tmp/%s"%(uuid.uuid4())
        to_bed(right, coordinates+[right_table_match_name,columns_to_join],tmp_out)
        tmp_obj = BedTool(tmp_out)
        bed3 = left_obj.intersect(tmp_obj,wao=True)
        tmp_df = bed3.to_dataframe()
        ## overlap and name match
        overlap = tmp_df[(tmp_df['name']==tmp_df['thickEnd'])&(tmp_df['blockCount']>0)]
        ## add column "itemRgb", this is just a dummy name, infered by pybedtools
        set_bed_index(overlap,['chrom','start','end','name'])
        my_tmp_dict= merge_dict(my_tmp_dict,overlap['itemRgb'].to_dict())
    my_replace_dict = merge_dict(my_replace_dict,my_tmp_dict)
    left[columns_to_join] = left[columns_to_join].replace(my_replace_dict)

left.to_csv("joined_table.csv")















