# SRH analysis pipeline

1. run ```start_over.sh``` to clean out any old results

2. run ```srh.py``` to get SRH stats and IQtree input files for all datasets in raw_data

3. run ```run_iqtree.sh``` to run IQ tree on all the folders from 2 (everything in ```processed_data/IQtree```)

4. run ```tree_dist.r``` this creates ```processed_data/tree_distances.csv``` which is a CSV file of tree-to-tree path distances comparing all three trees generated from each of the tree tests for each dataset. 

5. get all the results of the topology comparisons 