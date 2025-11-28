#!/bin/bash
# make result directory
mkdir -p result

# if adas_track_data.root exists, then delete
if [ -f "result/adas_track_data.root" ]; then
    echo "Micromegas data already existed, deleting ..."
	rm result/adas_track_data.root
fi
hadd result/adas_track_data.root mmdata/*

# if result_decode.root exists, then delete
if [ -f "result/result_decode.root" ]; then
	echo "Decode data already existed, deleting ..."
	rm result/result_decode.root
fi 
hadd result/result_decode.root decodedata/*

root -l -b -q ana.cxx
