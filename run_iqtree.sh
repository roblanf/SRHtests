threads=20

# run IQ tree on every dataset.
folders=$(find /data/srh/processed_data/IQtree -type f -name 'alignment.nex' -printf '%h\n' | sort -u)
echo "$folders" | parallel -P $threads iqtree -s {}"/alignment.nex" -spp {}"/partition.nex" -bb 1000

# make a tree file for each test dataset, that contains all three trees for that test and that dataset
for folder in $(find /data/srh/processed_data/IQtree/ -maxdepth 2 -mindepth 2 -type d); do

	cd $folder
	echo $folder
	cat ./All/partition.nex.treefile ./Bad/partition.nex.treefile ./Not_Bad/partition.nex.treefile > trees.nex

	cp trees.nex All/trees.nex 
	cp trees.nex Bad/trees.nex 
	cp trees.nex Not_Bad/trees.nex

done

# now run IQtree with topology tests (this re-estimates ML trees, but oh well for now)
folders=$(find /data/srh/processed_data/IQtree -type f -name 'alignment.nex' -printf '%h\n' | sort -u)
echo "$folders" | parallel -P $threads iqtree -s {}"/alignment.nex" -spp {}"/partition.nex" -bb 1000 -z {}"/trees.nex" -zb 10000 -zw -au -redo

