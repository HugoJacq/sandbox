#!/usr/bin/env bash

# INVESTIGATING NUMERICAL STABILITY FOR THE MULTILAYER
#
# For this script to work properly, you need to have
# - defined $BASILISK to the src/ directory of Basilisk
# - if you want netcdf, the netcdf.h header
# - the source file 'no_forcing_instrumented.c'
# - the modified multilayer files (see ./modified_layered/)
# - a Makefile with Basilisk's recipes
#
# run the tests with
# ./run_test_instrumented.sh 2>&1 | tee test_instrumented.log
#
# produce the figures with
# python3 figures_report.py
#
# All results will be in './results_intrumented'
#

set -x

# ==================
# Default parameters
# ==================
N=1024       # Horizontal resolution
NL=60        # Number of layers
NT0=100.0    # Number of periods targeted
RE=400.0     # Reynolds number, nu=1/Re, u diffusion (isotropic)
L0=100.0     # Box size (not necessarily related to peak wave number)
K_=10.0      # Peak wavelength
H_=5.0       # Depth of bassin
AK=0.01      # Steepness
CFL_H=1.0    # Barotropic CFL (advection of eta)
CFL=0.5      # Non hydrostatic CFL (advection of u, tracers)
THETA_H=0.5  # Numerical parameter to dump fast barotropic modes
MAXSLOPE=1.0 # Breaking parameter
THETA=1.3    # minmod2: 1 => minmod, 2 => superbee
DT=0.1       # this value is about the value for base case with DT auto

# Folders
NAME=no_forcing_instrumented
OUT=results_instrumented
mkdir -p $OUT

# A simple function that copy output files to a folder
save_results() {
  mkdir -p $1
  cp $NAME/{out,log} $1
}

# =======================
# Base case
# =======================

# compile and run base case
# echo "Running default case"
# make -k clean
# rm $NAME/plots
# make $NAME.tst $NAME/plots
# save_results $OUT/base
# exit 0

# =======================
# Rui's case
# =======================

# ------------------------------------------
# Rui's parameters, assuming default when not specified in preprint
# echo "Rui's case"
# cd $NAME
# ./$NAME 2048 15 300. 40000 1. 1. 5. 0.05 1.0 0.5 0.55 1.0 1.3 2>log >out
# cd ..
# save_results $OUT/Ruis_thetah0.5

# ------------------------------------------
# echo 'increasing theta_h for stability'
# cd $NAME
# ./$NAME 2048 15 300. 40000 1. 1. 5. 0.05 1.0 0.5 0.51 1.0 1.3 2>log >out
# cd ..
# mv $NAME/out $OUT/Ruis_thetah0.51/out.nc
# save_results $OUT/Ruis_thetah0.51

# ------------------------------------------
# echo 'increasing theta_h for stability'
# cd $NAME
# ./$NAME 2048 15 300. 40000 1. 1. 5. 0.05 1.0 0.5 0.55 1.0 1.3 2>log >out
# cd ..
# mv $NAME/out $OUT/Ruis_thetah0.55/out.nc
# save_results $OUT/Ruis_thetah0.55

# ------------------------------------------
echo "Rui's case but with very small CFL"
cd $NAME
./$NAME 2048 15 100. 40000 1. 1. 5. 0.05 0.05 0.05 0.51 1.0 1.3 2>log >out
cd ..
mv $NAME/out $OUT/Ruis_thetah_smallCFL/out.nc
save_results $OUT/Ruis_thetah_smallCFL

# =======================
# A well resolved case
# =======================
# echo 'a well resolved test' # at least it should be !
# 100 points per wavelength
# 100 points per period
# high Re to be sure that viscous term is small and that we can keep the analytic result
# small steepness to stay in linear regime
# depth = 5m, a=ak/k~0.008m deep water regime
# theta_h > 0.5 to ensure stability
# DT = T0/100 = 0.008 (bypass CFL_H and CFL)
# cd $NAME
# ./$NAME 1024 30 100. 40000. 1. 1. 5. 0.05 1.0 0.5 0.51 1.0 1.3 0.008 2>log >out
# cd ..
# save_results $OUT/resolved_thetah0.51

# ---------------------------------------------
# echo 'reducing theta_H to its default value'
# cd $NAME
# ./$NAME 1024 30 100. 40000. 1. 1. 5. 0.05 1.0 0.5 0.5 1.0 1.3 0.008 2>log >out
# cd ..
# rm $NAME/plots
# make $NAME/plots
# mv $NAME/out $OUT/resolved_thetah0.5/out.nc
# save_results $OUT/resolved_thetah0.5

# ---------------------------------------------
# echo 'Increasing theta_H'
# cd $NAME
# ./$NAME 1024 30 100. 40000. 1. 1. 5. 0.05 1.0 0.5 0.52 1.0 1.3 0.008 2>log >out
# cd ..
# mv $NAME/out $OUT/resolved_thetah0.52/out.nc
# save_results $OUT/resolved_thetah0.52

#=================
# OTHER TESTS
#=================
# ------------------------------------------
# echo 'TESTING N'
# NDEF=$N
# list_N=(256 512 2048 4096) # 512 2048 4096
# for i in "${list_N[@]}"; do
#   echo $i
#   N=$i
#   cd $NAME
#   ./$NAME $N $NL $NT0 $RE $L0 $K_ $H_ $AK $CFL_H $CFL $THETA_H $MAXSLOPE $THETA 2>log >out
#   cd ..
#   echo $?
#   save_results $OUT/N_$i
# done
# N=$NDEF

# ------------------------------------------
# echo 'TESTING theta'
# TTDEF=$THETA
# list_theta=(1.0 1.5 2.0) #
# for i in "${list_theta[@]}"; do
#   echo $i
#   THETA=$i
#   cd $NAME
#   ./$NAME $N $NL $NT0 $RE $L0 $K_ $H_ $AK $CFL_H $CFL $THETA_H $MAXSLOPE $THETA 2>log >out
#   cd ..
#   echo $?
#   save_results $OUT/theta_$i
# done
# THETA=$TTDEF

# ------------------------------------------
# echo 'TESTING Re'
# REDEF=$RE
# list_RE=(1000. 5000. 10000. 40000.)
# for i in "${list_RE[@]}"; do
#   echo $i
#   RE=$i
#   cd $NAME
#   ./$NAME $N $NL $NT0 $RE $L0 $K_ $H_ $AK $CFL_H $CFL $THETA_H $MAXSLOPE $THETA 2>log >out
#   cd ..
#   echo $?
#   save_results $OUT/RE_$i
# done
# RE=$REDEF

# ------------------------------------------
# echo 'TESTING ak'
# AKDEF=$AK
# list_AK=(0.05 0.08 0.1 0.2 0.3) #
# for i in "${list_AK[@]}"; do
#   echo $i
#   AK=$i
#   cd $NAME
#   ./$NAME $N $NL $NT0 $RE $L0 $K_ $H_ $AK $CFL_H $CFL $THETA_H $MAXSLOPE $THETA 2>log >out
#   cd ..
#   echo $?
#   save_results $OUT/ak_$i
# done
# AK=$AKDEF

# ------------------------------------------
# echo 'TESTING theta'
# TTDEF_H=$THETA_H
# list_thetah=(0.52 0.54 0.56) #
# for i in "${list_thetah[@]}"; do
#   echo $i
#   THETA_H=$i
#   cd $NAME
#   ./$NAME $N $NL $NT0 $RE $L0 $K_ $H_ $AK $CFL_H $CFL $THETA_H $MAXSLOPE $THETA 2>log >out
#   cd ..
#   echo $?
#   save_results $OUT/theta_h_$i
# done
# THETA_H=$TTDEF_H

# ------------------------------------------
# echo 'TESTING DT'
# DTDEF=$DT
# list_dt=(0.05 0.07 0.08 0.09 0.1)
# for i in "${list_dt[@]}"; do
#   echo $i
#   DT=$i
#   cd $NAME
#   ./$NAME $N $NL $NT0 $RE $L0 $K_ $H_ $AK $CFL_H $CFL $THETA_H $MAXSLOPE $THETA $DT 2>log >out
#   cd ..
#   save_results $OUT/DT_$i
# done
# DT=$DTDEF

# ------------------------------------------
# echo 'TESTING CFL_H'
# CFLHDEF=$CFL_H
# list_CFL_H=(0.3 0.5 0.7 0.8 0.9)
# for i in "${list_CFL_H[@]}"; do
#   echo $i
#   CFL_H=$i
#   cd $NAME
#   ./$NAME $N $NL $NT0 $RE $L0 $K_ $H_ $AK $CFL_H $CFL $THETA_H $MAXSLOPE $THETA 2>log >out
#   cd ..
#   save_results $OUT/CFL_H_$i
# done
# CFL_H=$CFLHDEF

# ------------------------------------------
# echo 'TESTING CFL'
# CFLDEF=$CFL
# list_CFL=(0.4 0.6 0.7 0.8 0.9 1.0)
# for i in "${list_CFL[@]}"; do
#   echo $i
#   CFL=$i
#   cd $NAME
#   ./$NAME $N $NL $NT0 $RE $L0 $K_ $H_ $AK $CFL_H $CFL $THETA_H $MAXSLOPE $THETA 2>log >out
#   cd ..
#   save_results $OUT/CFL_$i
# done
# CFL=$CFLDEF

# ========================
# PLOTTING
# ========================
python3 figures_report.py
