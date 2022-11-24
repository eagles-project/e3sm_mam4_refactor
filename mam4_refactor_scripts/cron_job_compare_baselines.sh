#!/bin/sh

#Unlimited stack
ulimit -d unlimited
ulimit -s unlimited
ulimit -c unlimited


#===========#===========#===========#===========#===========#===========#===========#===========
# USAGE:
# 1. This script can be launched from anywhere as it doesn't depend on the local directory
# 2. Issue command:
#           $bash  cron_job_compare_baselines.sh

# 3. The command above will create a temporary directory, clone the code and launch the comparison
#    tests
#===========#===========#===========#===========#===========#===========#===========#===========

main() {


#---------------------------------------------------------------
# User input
#---------------------------------------------------------------
scratch_dir="/compyfs/$USER/cron_jobs_refactor" #scratch directory to clone the code and run the tests
baselines="mam4_org_v2_baselines_created_11_23_2022"
#---------------------------------------------------------------
# User input ENDs
#---------------------------------------------------------------

#===================================================================
#===================================================================

#cronjobs are very lightweight, so source the modules so that they are available
source /etc/profile.d/modules.sh


#form a unique temporary directory name
temp_dir="test_`date +'%m-%d-%Y__%H_%M_%S'`"

newline && time_elapsed_min

#cd into the scratch directory
cd $scratch_dir
if_dir_exists_then_exit $scratch_dir/$temp_dir # exit if temp dir already exists

#create a new temporary directory:
newline && time_elapsed_min
echo "Creating new temporary directory: $temp_dir"
mkdir $temp_dir
cd $temp_dir

#clone the code
newline && time_elapsed_min
echo 'Cloning E3SM...'
git clone git@github.com:eagles-project/e3sm_mam4_refactor.git > /dev/null

# cd into the code directory
cd e3sm_mam4_refactor

#Merge branch into the maint-2.0 branch
newline && time_elapsed_min

#update submodules
newline && time_elapsed_min
echo 'Update submodules...'
git submodule update --init --recursive > /dev/null

#repeating this command as sometime it fails due to connection issues
git submodule update --init --recursive > /dev/null

# cd into cime/scripts
cd cime/scripts

#launch tests
newline && time_elapsed_min
echo "Launch Tests:"
./create_test eagles_mam -c -b $baselines -t comp_${temp_dir}_mam4_org_v2_baselines \
    -p esmd -r $scratch_dir/$temp_dir --output-root $scratch_dir/$temp_dir  > /dev/null

newline && time_elapsed_min

echo "-------------------------------------------------------------------------------------------------"
echo "-------------------------------------------------------------------------------------------------"
echo "To view the result:"
newline
echo "cd $scratch_dir/$temp_dir"
echo "./cs.status.comp_${temp_dir}_mam4_org_v2_baselines -s"
echo "-------------------------------------------------------------------------------------------------"
echo "-------------------------------------------------------------------------------------------------"
}

#---------------------
# Function Definitions
#---------------------

#Generate a newline
newline () { echo ''; }

#Time elapsed in minutes
time_elapsed_min() {
    endtime=$(date +%s)
    total_time=$(( $endtime - $starttime ))
    minutes=$((total_time / 60))
    seconds=$((total_time % 60))
    echo "Time elapsed (min): " $minutes:$seconds
}

if_dir_exists_then_exit () {
    if [ -d $1 ]; then
        echo "$1 directory already exists. Please remove it or choose a different temporary directory to proceed"
        exit 1
    fi
}

#capture start time to compute time elapsed in minutes
starttime=$(date +%s)
echo 'Start Time:' $(date +%T)

#call the main function
main
