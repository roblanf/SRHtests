# run IQ tree on a whole bunch of datasets
for analysis in $(find /data/srh/processed_data/IQtree/ -name "alignment.nex"); do

    echo $(dirname $analysis)
    cd $(dirname $analysis)
    iqtree-omp -s alignment.nex -sp partition.nex -nt AUTO -b 100

done
