# run IQ tree on a whole bunch of datasets
threads=5

for analysis in $(find /data/srh/processed_data/IQtree/ -name "alignment.nex"); do

    cd $(dirname $analysis)
    iqtree-omp -s alignment.nex -sp partition.nex -nt $threads -b 100

done