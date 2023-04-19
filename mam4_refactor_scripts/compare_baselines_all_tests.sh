#!/bin/sh

#===========#===========#===========#===========#===========#===========#===========#===========
# USAGE:
# 1. This script can be launched from anywhere as it doesn't depend on the local directory
# 2. Issue command:
#           $bash  compare_baselines_all_tests.sh -b <branch-to-test> -t <temporary-directory-name>

# 3. The command above will create a temporary directory, clone the code and launch the comparison
#    tests
#===========#===========#===========#===========#===========#===========#===========#===========

main() {


#---------------------------------------------------------------
# User input
#---------------------------------------------------------------
scratch_dir="/compyfs/$USER/e3sm_scratch/" #scratch directory to clone the code and run the tests
baselines="mam4_org_v2_baselines_created_11_23_2022"
#---------------------------------------------------------------
# User input ENDs
#---------------------------------------------------------------

#===================================================================
#===================================================================

newline && time_elapsed_min
echo "Testing branch:$branch_name"

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
echo "Check out and merge branch:$branch_name"
git checkout $branch_name > /dev/null
git checkout refactor-maint-2.0 > /dev/null
git merge $branch_name > /dev/null

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
    -p esmd -r $scratch_dir/$temp_dir --output-root $scratch_dir/$temp_dir > /dev/null

newline && time_elapsed_min

echo "-------------------------------------------------------------------------------------------------"
echo "-------------------------------------------------------------------------------------------------"
echo "To view the result:"
newline
echo "$scratch_dir$temp_dir/cs.status.comp_${temp_dir}_mam4_org_v2_baselines -s"
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


#parse command line args
while getopts ":b:t:" opt; do
  case $opt in
    b) branch_name="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG; please set branch name using -b command line option" >&2
    exit 1
    ;;
    t) temp_dir="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG; please set temp directory name using -t command line option" >&2
    exit 1
    ;;
  esac

  case $OPTARG in
    -*) echo "Option $opt needs a valid argument"
    exit 1
    ;;
  esac
done

if [ -z "${branch_name}" ]; then
    echo "branch name is not set, please set it using -b command line option"
    exit 1
fi

if [ -z "${temp_dir}" ]; then
    echo "temporary directory is not set, please set it using -t command line option"
    exit 1
fi


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
