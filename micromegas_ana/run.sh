#!/bin/bash
# make result directory
# option=""
option="3"
result_dir="result/result${option}"
mkdir -p $result_dir
mm_data="${result_dir}/adas_track_data${option}.root"
decode_data="${result_dir}/result_decode${option}.root"
mm_data_dir="./data/mmdata/mmdata${option}"
decode_data_dir="./data/decodedata/decodedata${option}"

# if adas_track_data.root exists, then delete
if [ -f $mm_data ]; then
    echo "Micromegas data already existed, deleting ..."
	rm $mm_data
fi
hadd $mm_data $mm_data_dir/*

# if result_decode.root exists, then delete
if [ -f $decode_data ]; then
	echo "Decode data already existed, deleting ..."
	rm $decode_data
fi 
hadd $decode_data $decode_data_dir/*

root -l -b -q "ana.cxx(\"$option\")"
