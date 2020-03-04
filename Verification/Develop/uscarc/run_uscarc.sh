sbatch -J uscarc1_single $SCRIPTS/fds_opt1.sh uscarc1_single.fds
sbatch -J uscarc2_single $SCRIPTS/fds_opt2.sh uscarc2_single.fds
sbatch -J uscarc5a_single $SCRIPTS/fds_opt5.sh uscarc5a_single.fds
sbatch -J uscarc10a_single $SCRIPTS/fds_opt10.sh uscarc10a_single.fds
sbatch -J uscarc10b_single $SCRIPTS/fds_opt10.sh uscarc10b_single.fds
sbatch -J uscarc20_single $SCRIPTS/fds_opt20.sh uscarc20_single.fds

sbatch -J uscarc1_double $SCRIPTS/fds_opt1.sh uscarc1_double.fds
sbatch -J uscarc2_double $SCRIPTS/fds_opt2.sh uscarc2_double.fds
sbatch -J uscarc5a_double $SCRIPTS/fds_opt5.sh uscarc5a_double.fds
sbatch -J uscarc10a_double $SCRIPTS/fds_opt10.sh uscarc10a_double.fds
sbatch -J uscarc10b_double $SCRIPTS/fds_opt10.sh uscarc10b_double.fds
sbatch -J uscarc20_double $SCRIPTS/fds_opt20.sh uscarc20_double.fds

