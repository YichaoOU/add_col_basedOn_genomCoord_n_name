## Get Started


Append columns from the right table to the left table based on a common name and overlapping genomic locations


### Dependency

**pybedtools** and **pandas**


### Edit line 82 to 93 to specify your need.

```

## Define input ##
left_table = "CHANGEseq_GUIDEseq_classes_validation_set_20191212.csv"
right_table = "GUIDEseq_matched_subset_20191212.csv"

## Define parameters ##
left_table_bed3 = ['chr','start','end']
left_table_match_name = "sample"


right_table_bed3_list = [['chr','Site_SubstitutionsOnly.Start','Site_SubstitutionsOnly.End'],['chr','Site_GapsAllowed.Start','Site_GapsAllowed.End']]
right_table_match_name = "sample"
right_table_columns_to_join = ['Site_SubstitutionsOnly.Sequence','Site_GapsAllowed.Sequence']

```

Define Input in `left_table` and `right_table`

Define left table genomic location `left_table_bed3` and a name column to do string match `left_table_match_name`

#### Note that left table should have only one genomic location

Define right table genomic location lists `right_table_bed3_list`. Here you can define multipe locations to use. 

Define right table columns to be added to the left `right_table_columns_to_join`

Also, define a name column to do string match `right_table_match_name`

### Workflow

For every columns in `right_table_match_name`, values in each row will be added to the left table if: (1) its genomic location overlaps with left **AND** (2) its name matches to the left.

For condition (1), since right table allows multiple genomic locations, the values will be added if anyone of them overlaps with left.



