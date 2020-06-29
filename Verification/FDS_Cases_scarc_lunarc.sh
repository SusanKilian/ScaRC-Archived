#!/bin/bash

#VERIFICATION=`pwd`
#SCRIPTS="../Lunarc"

VERIFICATION="/home/skil/GIT/Github/ScaRC/Verification"
SCRIPTS="/home/skil/Scripts"

#
# NS_Analytical_Solution 
#
cd $VERIFICATION/NS_Analytical_Solution

sbatch -J ns2d_4mesh_8_scarc  $SCRIPTS/scarc_opt4.sh ns2d_4mesh_8_scarc.fds
sbatch -J ns2d_4mesh_8_scarc  $SCRIPTS/scarc_opt4.sh ns2d_4mesh_8_scarc_nupt1.fds
sbatch -J ns2d_4mesh_16_scarc $SCRIPTS/scarc_opt4.sh ns2d_4mesh_16_scarc.fds
sbatch -J ns2d_4mesh_16_scarc $SCRIPTS/scarc_opt4.sh ns2d_4mesh_16_scarc_nupt1.fds
sbatch -J ns2d_4mesh_32_scarc $SCRIPTS/scarc_opt4.sh ns2d_4mesh_32_scarc.fds
sbatch -J ns2d_4mesh_32_scarc $SCRIPTS/scarc_opt4.sh ns2d_4mesh_32_scarc_nupt1.fds
sbatch -J ns2d_4mesh_64_scarc $SCRIPTS/scarc_opt4.sh ns2d_4mesh_64_scarc.fds
sbatch -J ns2d_4mesh_64_scarc $SCRIPTS/scarc_opt4.sh ns2d_4mesh_64_scarc_nupt1.fds

#
# Pressure_Effects
#
cd $VERIFICATION/Pressure_Effects

sbatch -J pressure_boundary $SCRIPTS/scarc_opt1.sh pressure_boundary.fds
sbatch -J pressure_boundary_scarc $SCRIPTS/scarc_opt10.sh pressure_boundary_scarc.fds

sbatch -J pressure_rise $SCRIPTS/scarc_opt3.sh pressure_rise.fds
sbatch -J pressure_rise_glmat $SCRIPTS/scarc_opt3.sh pressure_rise_glmat.fds
sbatch -J pressure_rise_scarc $SCRIPTS/scarc_opt3.sh pressure_rise_scarc.fds

sbatch -J zone_break_fast_scarc $SCRIPTS/scarc_opt4.sh zone_break_fast_scarc.fds
sbatch -J zone_break_fast $SCRIPTS/scarc_opt1.sh zone_break_fast.fds
sbatch -J zone_break_slow_scarc $SCRIPTS/scarc_opt4.sh zone_break_slow_scarc.fds
sbatch -J zone_break_slow $SCRIPTS/scarc_opt1.sh zone_break_slow.fds

sbatch -J zone_shape $SCRIPTS/scarc_opt2.sh zone_shape.fds
sbatch -J zone_shape_scarc $SCRIPTS/scarc_opt2.sh zone_shape_scarc.fds

sbatch -J zone_shape_2 $SCRIPTS/scarc_opt8.sh zone_shape_2.fds
sbatch -J zone_shape_2_uscarc $SCRIPTS/scarc_opt8.sh zone_shape_2_uscarc.fds

#
# Pressure_Solver 
#
cd $VERIFICATION/Pressure_Solver

sbatch -J dancing_eddies_default $SCRIPTS/scarc_opt4.sh dancing_eddies_default.fds
sbatch -J dancing_eddies_scarc $SCRIPTS/scarc_opt4.sh dancing_eddies_scarc.fds
sbatch -J dancing_eddies_scarc_tight $SCRIPTS/scarc_opt4.sh dancing_eddies_scarc_tight.fds
sbatch -J dancing_eddies_tight $SCRIPTS/scarc_opt4.sh dancing_eddies_tight.fds
sbatch -J dancing_eddies_tight_overlap $SCRIPTS/scarc_opt4.sh dancing_eddies_tight_overlap.fds
sbatch -J dancing_eddies_uglmat $SCRIPTS/scarc_opt4.sh dancing_eddies_uglmat.fds
sbatch -J dancing_eddies_uscarc $SCRIPTS/scarc_opt4.sh dancing_eddies_uscarc.fds
sbatch -J dancing_eddies_1mesh $SCRIPTS/scarc_opt1.sh dancing_eddies_1mesh.fds

sbatch -J duct_flow $SCRIPTS/scarc_opt8.sh duct_flow.fds
sbatch -J duct_flow_scarc $SCRIPTS/scarc_opt8.sh duct_flow_scarc.fds
sbatch -J duct_flow_uglmat $SCRIPTS/scarc_opt8.sh duct_flow_uglmat.fds
sbatch -J duct_flow_uscarc $SCRIPTS/scarc_opt8.sh duct_flow_uscarc.fds

sbatch -J hallways $SCRIPTS/scarc_opt5.sh hallways.fds
sbatch -J hallways_scarc $SCRIPTS/scarc_opt5.sh hallways_scarc.fds

sbatch -J poisson2d_fft $SCRIPTS/scarc_opt44.sh poisson2d_fft.fds
sbatch -J poisson2d_scarc $SCRIPTS/scarc_opt44.sh poisson2d_scarc.fds
sbatch -J poisson2d_uglmat $SCRIPTS/scarc_opt44.sh poisson2d_uglmat.fds
sbatch -J poisson2d_uscarc $SCRIPTS/scarc_opt44.sh poisson2d_uscarc.fds
sbatch -J poisson2d_fft_tight $SCRIPTS/scarc_opt44.sh poisson2d_fft_tight.fds
sbatch -J poisson2d_scarc_tight $SCRIPTS/scarc_opt44.sh poisson2d_scarc_tight.fds

sbatch -J poisson3d_fft $SCRIPTS/scarc_opt4.sh poisson3d_fft.fds
sbatch -J poisson3d_glmat $SCRIPTS/scarc_opt4.sh poisson3d_glmat.fds
sbatch -J poisson3d_scarc $SCRIPTS/scarc_opt4.sh poisson3d_scarc.fds
sbatch -J poisson3d_uglmat $SCRIPTS/scarc_opt4.sh poisson3d_uglmat.fds
sbatch -J poisson3d_uscarc $SCRIPTS/scarc_opt4.sh poisson3d_uscarc.fds
sbatch -J tunnel_demo $SCRIPTS/scarc_opt8.sh tunnel_demo.fds
sbatch -J tunnel_demo_glmat $SCRIPTS/scarc_opt8.sh tunnel_demo_glmat.fds
sbatch -J tunnel_demo_scarc $SCRIPTS/scarc_opt8.sh tunnel_demo_scarc.fds
sbatch -J tunnel_demo_scarc_amg $SCRIPTS/scarc_opt8.sh tunnel_demo_scarc_amg.fds

#
# Scalar_Analytical_Solution
#
cd $VERIFICATION/Scalar_Analytical_Solutiion

sbatch -J shunn3_16mesh_128_scarc $SCRIPTS/scarc_opt16.sh shunn3_16mesh_128_scarc.fds
sbatch -J shunn3_16mesh_256_scarc $SCRIPTS/scarc_opt16.sh shunn3_16mesh_256_scarc.fds
sbatch -J shunn3_16mesh_32_scarc $SCRIPTS/scarc_opt16.sh shunn3_16mesh_32_scarc.fds
sbatch -J shunn3_16mesh_512_scarc $SCRIPTS/scarc_opt16.sh shunn3_16mesh_512_scarc.fds
sbatch -J shunn3_16mesh_64_scarc $SCRIPTS/scarc_opt16.sh shunn3_16mesh_64_scarc.fds
sbatch -J shunn3_4mesh_128_scarc $SCRIPTS/scarc_opt4.sh shunn3_4mesh_128_scarc.fds
sbatch -J shunn3_4mesh_256_scarc $SCRIPTS/scarc_opt4.sh shunn3_4mesh_256_scarc.fds
sbatch -J shunn3_4mesh_32_scarc $SCRIPTS/scarc_opt4.sh shunn3_4mesh_32_scarc.fds
sbatch -J shunn3_4mesh_512_scarc $SCRIPTS/scarc_opt4.sh shunn3_4mesh_512_scarc.fds
sbatch -J shunn3_4mesh_64_scarc $SCRIPTS/scarc_opt4.sh shunn3_4mesh_64_scarc.fds

#
#Turbulence
#
cd $VERIFICATION/Turbulence

sbatch -J ribbed_channel_160_4mesh_uscarc $SCRIPTS/scarc_opt4.sh ribbed_channel_160_4mesh_uscarc.fds
sbatch -J ribbed_channel_20_4mesh_uscarc $SCRIPTS/scarc_opt4.sh ribbed_channel_20_4mesh_uscarc.fds
sbatch -J ribbed_channel_40_4mesh_uscarc $SCRIPTS/scarc_opt4.sh ribbed_channel_40_4mesh_uscarc.fds
sbatch -J ribbed_channel_80_4mesh_uscarc $SCRIPTS/scarc_opt4.sh ribbed_channel_80_4mesh_uscarc.fds
sbatch -J ribbed_channel_20 $SCRIPTS/scarc_opt1.sh ribbed_channel_20.fds
sbatch -J ribbed_channel_40 $SCRIPTS/scarc_opt1.sh ribbed_channel_40.fds
sbatch -J ribbed_channel_80 $SCRIPTS/scarc_opt1.sh ribbed_channel_80.fds
sbatch -J ribbed_channel_160 $SCRIPTS/scarc_opt1.sh ribbed_channel_160.fds

cd VERIFICATION
