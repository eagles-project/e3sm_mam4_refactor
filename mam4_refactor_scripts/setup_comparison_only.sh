#!/bin/sh

#===========#===========#===========#===========#===========#===========#===========#===========
# USAGE:
# 1. Go to the root of the E3SM codebase (i.e., e3sm_mam4_refactor directory)
# 2. Issue command:
#           $bash  mam4_refactor_scripts/setup_comparison_only.sh -t <tag>

# 3. The command above will launch a comparison build of the E3SM
#    The tag for the comparison run will be appended by "comp_"
#===========#===========#===========#===========#===========#===========#===========#===========

main() {


    #-------------------
    #user input starts
    #-------------------

    #--------------------------------------
    #Less frequently changed user input
    #--------------------------------------

    #compiler
    compiler="intel"

    #scratch directory path
    scratch_dir="/compyfs/$USER/e3sm_scratch"

    #baseline id
    baseline_id="mam4_org_v2_baselines"

    test_name="SMS_D_Ln5_P32x1"
    grid="ne4pg2_oQU480"
    compset="F2010"
    project="esmd"

    #test_id for the test is obtained from the command line arg

    #Test_id for the comparison simulation
    comparison_test_id="comp_"$test_id


    #----------------------------------------
    #user input ends
    #----------------------------------------

    #======#======#======#======#======#======#======#======#======#======#======#======#======#======#======#======#======

    #create directory names
    test_dir="$test_name.$grid.$compset.compy_$compiler" #test dir (common name)
    comparison_dir="$scratch_dir/$test_dir.C.$comparison_test_id"

    #if any of the generate or comparison directories already exist, exit
    if_dir_exists_then_exit $comparison_dir "Comparison baseline"


    script_name="create_test" #test script name
    status_file="TestStatus" #status file to check test status

    #----- Start running tests ---------

    #create_test script exists in cime/scripts directory, cd into that directory
    cd cime/scripts

    #check if we are in the right directory and the script exists
    if [ ! -f $script_name ]; then
        echo "$script_name does not exist. Have you initialized submodules? "
        echo "Current directory is:" `pwd`
        newline
        exit -1 #exit if it is wrong directory or file doesn't exist
    fi

    newline && time_elapsed_min

    #Launch the comparison simulation but DO NOT submit it
    echo 'Compare baselines build starts...'
    ./$script_name $test_name.$grid.$compset --compiler $compiler -t $comparison_test_id -c \
        -b $baseline_id -p $project > /dev/null &

    newline && time_elapsed_min

    #Find out if the comparison directory is generated:
    echo 'Check if comparison baseline directory is there ....'
    while [ ! -d $comparison_dir ]
    do
        echo -n '.'
        sleep 2
    done

    echo 'Link fast_compile_run.sh script in the comparison directory..'
    cd $comparison_dir
    /bin/ln -sf  `./xmlquery SRCROOT --value`/mam4_refactor_scripts/fast_compile_run.sh .

    newline && time_elapsed_min

    echo '-------------------------------------------------------------------'
    echo "Comparison directory is: $comparison_dir"
    echo '-------------------------------------------------------------------'


    echo 'Comparison run should be building at this time, wait for a few minutes (4-5 min) and then check the comparison directory...'


    newline && time_elapsed_min
}

#---------------------
# Function Definitions
#---------------------

#Generate a newline
newline () { echo ''; }

#Print time
#print_time() { echo -n "Time elapsed:" && date +%T ; }
#Time elapsed in minutes
time_elapsed_min() {
    endtime=$(date +%s)
    total_time=$(( $endtime - $starttime ))
    minutes=$((total_time / 60))
    seconds=$((total_time % 60))
    echo "Time elapsed (min): " $minutes:$seconds
}

#if directory exists, we cannot proceed
if_dir_exists_then_exit () {
    if [ -d $1 ]; then
        echo "$2 directory at: $1 already exists. Please remove it to proceed"
        exit -1
    fi
}

#parse command line args
while getopts ":t:" opt; do
  case $opt in
    t) test_id="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG; please set test_id using -t command line option" >&2
    exit 1
    ;;
  esac

  case $OPTARG in
    -*) echo "Option $opt needs a valid argument"
    exit 1
    ;;
  esac
done

if [ -z "${test_id}" ]; then
    echo "test_id is not set, please set it using -t command line option"
    exit 1
fi

#capture start time to compute time elapsed in minutes
starttime=$(date +%s)
echo 'Start Time:' $(date +%T)
#call the main function
main
