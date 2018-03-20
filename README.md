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

This will do three things in the following order (note there is a `threads` argument at the top of the script which you should change as appropriate, it also relies on GNU parallel):

* Run IQtree on all of the sub-folders in `/data/srh/processed_data/IQtree/` with the following command: `iqtree -s alignment.nex -spp partition.nex -bb 1000 -redo`

* Make one `trees.nex` file for each of the tests in each of the datasets, if and only if all three analyses for that test (i.e. All, Bad, Not_bad) produced trees. This file is then copied into each of the test subfolders (`/MPTS`, `/MPTMS`, `/MPTIS`).

* Run IQtree on all of the sub-folders in `/data/srh/processed_data/IQtree/`, this time including topology tests to compare the three trees in `trees.nex`, with the following command: `iqtree -s alignment.nex -spp partition.nex -bb 1000 -z {}"/trees.nex" -zb 10000 -zw -au -redo -safe`

##### 4. Run `mkdir /data/srh/tables`

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

This looks at the data from the binomial tests in `/data/srh/processed_data/SRH_tables/`, and just classifies each charset as to whether the binomial fails (p<0.05) or passes (p>0.05). It generates a csv file: `/data/srh/tables/is_charset_bad.csv`. The csv file is really a summary of the results of all the tests for every charset in one file.

##### 9. run `python trees_lengths.py` 

This extracts the tree lengths of every estimated tree in `/data/srh/processed_data/IQtree/` to generate a table that contains the lengths of the trees generated by IQtree from each partition (all, bad, not_bad) and symmetry test (MPTS, MPTMS, MPTIS). The csv file is: `/data/srh/tables/trees_lengths.csv`

##### 10. run `python topology_tests.py` 

This extracts the topology tests from `/data/srh/processed_data/IQtree/` to generate a table that contains all the topology tests results. The output is: `/data/srh/tables/topology_tests.csv`

##### 11. run `python summary_charsets.py`

This creates a big summary table of everything we know about every charset we have analysed. This is what we use to make figures from. The output is: `/data/srh/tables/summary_charsets.csv`

##### 12. run `python summary_trees.py`

This creates a big summary table of everything we know about every tree we have generated. This is what we use to make figures from. The output is: `/data/srh/tables/summary_trees.csv`


11. run `figure1.r`

12. run `figure2.r`

13. run `figure3.r`