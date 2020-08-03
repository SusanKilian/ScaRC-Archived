#!/bin/bash

#VERIFICATION=`pwd`
#SCRIPTS="../Lunarc"

VERIFICATION="/home/skil/GIT/Github/ScaRC/Verification"
SCRIPTS="/home/skil/Scripts"

#
# Pressure_Solver 
#
cd $VERIFICATION/Linear/2x1

rba.sh fft 2
rba.sh glmat 2
rba.sh scarc_amg_fft_agg_l2_sm8 2
rba.sh scarc_amg_fft_agg_l3_sm8 2
rba.sh scarc_amg_fft_cub_l2_sm8 2
rba.sh scarc_amg_fft_cub_l3_sm8 2
rba.sh scarc_amg_fft_l2_sm8 2
rba.sh scarc_amg_fft_l3_sm8 2
rba.sh scarc_amg_ssor_agg_l2_sm8 2
rba.sh scarc_amg_ssor_agg_l3_sm8 2
rba.sh scarc_amg_ssor_cub_l2_sm8 2
rba.sh scarc_amg_ssor_cub_l3_sm8 2
rba.sh scarc_amg_ssor_l2_sm8 2
rba.sh scarc_amg_ssor_l3_sm8 2
rba.sh scarc_cgamg_fft_agg_l2_sm8 2
rba.sh scarc_cgamg_fft_agg_l3_sm8 2
rba.sh scarc_cgamg_fft_cub_l2_sm8 2
rba.sh scarc_cgamg_fft_cub_l3_sm8 2
rba.sh scarc_cgamg_ssor_agg_l2_sm8 2
rba.sh scarc_cgamg_ssor_agg_l3_sm8 2
rba.sh scarc_cgamg_ssor_cub_l2_sm8 2
rba.sh scarc_cgamg_ssor_cub_l3_sm8 2
rba.sh scarc_cgamg_ssor_l2_sm8 2
rba.sh scarc_cgamg_ssor_l3_sm8 2
rba.sh scarc_cggmg_fft_l2_sm8 2
rba.sh scarc_cggmg_fft_l3_sm8 2
rba.sh scarc_cggmg_ssor_agg_l2_sm8 2
rba.sh scarc_cggmg_ssor_agg_l3_sm8 2
rba.sh scarc_cggmg_ssor_cub_l2_sm8 2
rba.sh scarc_cggmg_ssor_cub_l3_sm8 2
rba.sh scarc_cggmg_ssor_l2_sm8 2
rba.sh scarc_cggmg_ssor_l3_sm8 2
rba.sh scarc_gmg_fft_l2_sm8 2
rba.sh scarc_gmg_fft_l3_sm8 2
rba.sh scarc_gmg_ssor_l2_sm8 2
rba.sh scarc_gmg_ssor_l3_sm8 2
