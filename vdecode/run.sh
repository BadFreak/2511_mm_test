#!/bin/bash
# binfile=`ls ./data/*.dat`
# ./build/vdecode combine $binfile> log_decode 2>&1
# resultfile=`ls result_*.root`
# #vplot $resultfile > log_draw
# root -l -b -q draw.cxx\(1,\"$resultfile\"\) > log_draw 2>&1
# python VPlot.py > log_plot 2>&1


#!/bin/bash

# datadir="./data"
option="3"
datadir="./data/data${option}"
outdir="./result/result${option}"
mkdir -p $datadir

for dat in "$datadir"/*.dat; do
    filename=$(basename "$dat")
    prefix="${filename%.dat}"
	echo "Processing $filename"
    curr_out="$outdir/$prefix"
	if [ -f "$curr_out/result_${prefix}.root" ]; then
        echo "[Skip] $curr_out 已存在，跳过处理"
        continue
    fi
    mkdir -p "$curr_out"
    
    # decode
    ./build/vdecode combine "$dat" > "$curr_out/log_decode" 2>&1

    # find ROOT output

    resultfile="$curr_out/result_${prefix}.root"
	mv result_${prefix}.root "$curr_out/"

	# echo root output file absolute path
	ABSOLUTE_PATH=$(realpath "$curr_out/result_${prefix}.root")
	echo "Output File Path: $ABSOLUTE_PATH"
    # draw（注意这里的路径！）
    (
        cd "$curr_out"
        root -l -b -q "../../../draw.cxx(1, \"result_${prefix}.root\")" > log_draw 2>&1
    )

    # python plot
    (
        cd "$curr_out"
        python ../../../VPlot.py > log_plot 2>&1
    )
done

