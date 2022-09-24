#!/bin/sh
main() {


    #-------------------
    #user input starts
    #-------------------

    #Test_Id for the test
    test_id="tag1"

    #baseline test_id (must be different from the test_id above; you can simply add "base_" infront of test_id if you like)
    baseline_test_id="base_tag1"

    #Test_Id for the comparison simulation(must be different from the test_ids above; you can simply add "comp_" infront of test_id if you like)
    comparison_test_id="comp_tag1"

    #--------------------------------------
    #Less frequently changed user input
    #--------------------------------------

    #scratch directory path
    scratch_dir="/compyfs/sing201/e3sm_scratch"

    test_name="SMS_D_Ln5_P32x1"
    grid="ne4_oQU240"
    compset="F2010"
    project="esmd"

    #-------------------
    #user input ends
    #-------------------
    #======#======#======#======#======#======#======#======#======#======#======#======#======#======#======#======#======

    #create directory names
    test_dir="$test_name.$grid.$compset.compy_intel" #test dir (common name)
    generate_dir="$scratch_dir/$test_dir.G.$test_id"
    comparison_dir="$scratch_dir/$test_dir.C.$comparison_test_id"

    #if any of the generate or comparison directories already exist, exit
    dir_exists_then_exit $generate_dir "Generate baseline"
    dir_exists_then_exit $comparison_dir "Comparison baseline"


    script_name="create_test" #test script name
    status_file="TestStatus" #status file to check test status

    #create_test script exists in cime/scripts directory, cd into that directory
    cd cime/scripts

    #check if we are in the right directory and the script exists
    if [ ! -f $script_name ]; then
        echo "$script_name does not exist. Have you initialized submodules? "
        exit -1 #exit if it is wrong directory or file doesn't exist
    fi

    newline && print_time

    #run baseline test
    echo 'Generating baselines...'
    ./$script_name $test_name.$grid.$compset --compiler intel -t $test_id -g -b $baseline_test_id -p $project -o > /dev/null &

    newline && print_time

    #Launch the comparison simulation but DO NOT submit it
    echo 'Compare baselines build starts...'
    ./$script_name $test_name.$grid.$compset --compiler intel -t $comparison_test_id -c -b $baseline_test_id -p $project --no-run > /dev/null &

    newline && print_time
    echo "Generate directory is: $generate_dir"
    echo "Comparison directory is: $comparison_dir"

    newline && print_time
    #Find out if the baseline run is complete:
    echo 'Check if baselines has been generated...'
    while [ ! -d $generate_dir ]
    do
        echo 'waiting for 2 secs...'
        sleep 2
    done

    #cd into generate baseline directory
    cd $generate_dir

    baseline_generated=0

    #loop until baseline_generated variable catched "GENERATE" and "PASS" in the $status file
    while [ $baseline_generated -eq 0 ]
    do
        #see if TestStatus file exists
        if [ ! -f "$status_file" ]; then
            echo "$status_file doesnt exists, will check back in 2 minute... "
            sleep 120
            continue
        fi

        baseline_generated=`cat $status_file |grep GENERATE |grep PASS|wc -l`
        sleep 60
        echo 'Checking if baselines have been generated:'
    done

    newline && print_time

    echo 'Baselines have been generated, submiting comparison run..'

    cd $compare_dir

    ./case.submit > /dev/null #submit run

    echo 'Comparison run should be in queue...'

    newline && print_time

    echo "Comparison directory is: $comparison_dir"

    newline && print_time
}

#---------------------
# Function Definitions
#---------------------

#Generate a newline
newline () { echo ''; }

#Print time
print_time() { echo -n "Time:" && date +%T ; }

#if directory exists, we cannot proceed
dir_exists_then_exit () {
    if [ -d $1 ]; then
        echo "$2 directory at: $1 already exists. Please remove it to proceed"
        exit -1
    fi
}

#call the main function
main
