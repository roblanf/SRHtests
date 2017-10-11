# SRH analysis pipeline

All scripts have hard-coded input and output destinations. If you want to run them for yourself, you will need to adjust these destinations in each script as you go.

##### 1. run `sh start_over.sh`

This just deletes (with `rm -rf`, so be careful), `/data/srh/processed_data/SRH_tables/` and `/data/srh/processed_data/IQtree/`.

##### 2. run `python srh.py`

This will calculate SRH stats and IQtree input files for all datasets in `SRHtests/datasets`. It requires python 3.6.x or higher and dependencies as in the header of srh.py. Input and output files are hardcoded near the end of the script, change them if you need to. This script creates two output folders, each of which contains one folder for each dataset in `SRHtests/datasets`.

* `/data/srh/processed_data/SRH_tables/` where each folder has two files:

	* `data.csv` - raw data on pairwise MPTS, MPTIS, and MPTMS tests between all pairs of taxa in the alingment and for all CHARSETS in the alignment.

	* `table_binom.csv` - binomial tests for each CHARSET and each test (MPTS, MPTIS, MPTMS) which ask whether we observe more p<0.05 than we would expect by chance.

* `/data/srh/processed_data/IQtree/` each folder has three subfolders (`/MPTS`, `/MPTMS`, `/MPTIS`), each of which contain three more subfolders (`/All`, `/Bad`, `/Not_bad`). Each of these folders contains input files for IQtree. `/All` contains all of the CHARSETS in the alignment, `/Bad` contains the charsets with a significant binomial test result (i.e. p<0.05 for the binomial test), `Not_bad` contains the charsets with binomial tests results >=0.05. Each folder has the following files:

	* `alignment.nex` a copy of the complete alignment for that dataset 

	* `partition.nex` a description of the relevant charsets for that subfolder (i.e. all, bad, or not_bad)


##### 3. run `sh run_iqtree.sh` 

This will do three things in the following order (note there is a `threads` argument at the top of the script which you should change as appropriate):

* Run IQtree on all of the sub-folders in `/data/srh/processed_data/IQtree/` with the following command: `iqtree -s alignment.nex -spp /partition.nex -bb 1000 -redo`

* Make one `trees.nex` file for each of the tests in each of the datasets, if and only if all three analyses for that test (i.e. All, Bad, Not_bad) produced trees. This file is then copied into each of the test subfolders (`/MPTS`, `/MPTMS`, `/MPTIS`).

* Run IQtree on all of the sub-folders in `/data/srh/processed_data/IQtree/`, this time including topology tests to compare the three trees in `trees.nex`, with the following command: `iqtree -s {}"/alignment.nex" -spp {}"/partition.nex" -bb 1000 -z {}"/trees.nex" -zb 10000 -zw -au -redo -safe`

##### 4. Run 'mkdir /data/srh/tables'

This just makes a directory for the output of the following scripts which take the raw data from steps 1-3 and convert them into summary tables for analysis.


##### 5. run `Rscript tree_dist.r` 

This will measure normalised Path Distances between all three pairs of trees (All vs Bad, All vs Not_Bad, Bad vs Not_bad) for each test (MPTS, MPTIS, MPTMS) within each dataset. It requires a couple of libraries listed at the top of the script.

The file outputs the following things:

* PDFs of cophyloplots comparing each pair of trees. These are put into the dataset and test folder (e.g. `/processed_data/Anderson_2013/MPTS`) and have names like 'cophylo_all_bad.pdf'. They are there so you can see clearly (most of the time) any differences between trees.

* `/data/srh/tables/tree_distances.csv` which is a csv file that contains the pairwise tree distances for each pair of trees in each test for each dataset. I.e. there 9 comparisons total per dataset (3 for each of 3 tests). 

this creates ```processed_data/tree_distances.csv``` which is a CSV file of tree-to-tree path distances comparing all three trees generated from each of the tree tests for each dataset. 

##### 6. run `python charsets_percentage.py` 

This looks at the data in `/data/srh/processed_data/SRH_tables/` and calculates from these the percentage of bad and not_bad charsets in each dataset. 

This is then output to a csv file: `/data/srh/tables/charsets_percentage.csv`

##### 7. run `python charset_length.py`

This looks at the input datasets in `SRHtests/datasets`, and just makes a table with the names, charset names, and lengths of each charset in each dataset.

##### 8. run `python is_charset_bad.py` 

This generate a table that divides the character sets into two categories: one that rejects the null hypothesis and one that does not

9. run `trees_lengths.py` to generate a table that contains the lengths of the trees generated by IQtree from each partition (all, bad, not_bad) and symmetry test (MPTS, MPTMS, MPTIS)

10. run `topology_tests.py` to generate a table that contains all the topology tests results

11. run `figure1.r`

12. run `figure2.r`

13. run `figure3.r`