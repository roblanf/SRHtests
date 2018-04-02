threads=26

folders=$(find /data/srh/processed_data/IQtree/ -mindepth 5 -wholename '*Bad/alignment.nex' -printf '%h\n' | sort -u)
echo "$folders" | parallel -P $threads iqtree -s {}"/alignment.nex" -spp {}"/partition.nex"

for folder in $(find /data/srh/processed_data/IQtree/ -mindepth 4 -maxdepth 4 -wholename "*/Random*" -type d); do
cd $folder
echo $folder
cat ./All/partition.nex.treefile ./Bad/partition.nex.treefile ./Not_Bad/partition.nex.treefile > trees.nex
cp trees.nex All/trees.nex
cp trees.nex Bad/trees.nex
cp trees.nex Not_Bad/trees.nex
done

folders=$(find /data/srh/processed_data/IQtree/  -mindepth 5 -wholename '*/alignment.nex' -printf '%h\n' | sort -u)
echo "$folders" | parallel -P $threads iqtree -s {}"/alignment.nex" -spp {}"/partition.nex" -te {}"/partition.nex.treefile" -z {}"/trees.nex" -zb 10000 -zw -redo -safe
