#!/bin/bash
binfile=`ls *.dat`
vdecode combine $binfile> log_decode 2>&1
resultfile=`ls result_*.root`
#vplot $resultfile > log_draw
root -l -b -q draw.cxx\(1,\"$resultfile\"\) > log_draw 2>&1
python VPlot.py > log_plot 2>&1
