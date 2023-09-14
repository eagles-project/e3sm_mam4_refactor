#!/bin/sh

#===========#===========#===========#===========#===========#===========#===========#===========
# USAGE:
# 1. Go to the case directory
# 2. Copy script from the source code and place in the case directory
# 3. In the case directory, issue command:
#           $bash  fast_compile_run.sh

# 3. The command above will build and submit the run
#===========#===========#===========#===========#===========#===========#===========#===========

main () {

    #executable name
    executable="e3sm.exe"

    newline && time_elapsed_min
    #check if xmlquery is present
    if_file_not_present_then_exit xmlquery " Please cd into the case directory"

    #source the enviroment to build E3SM
    source ./.env_mach_specific.sh

    #get the build directory path
    bld_dir=`./xmlquery -value EXEROOT`
    case_dir=`pwd`

    #remove executable file in the bld directory
    echo  "Removing $executable"
    rm -rf $bld_dir/$executable

    #remove TestStatus file
    echo  "Removing TestStatus file and creating an empty file"
    rm -rf TestStatus
    /bin/touch TestStatus

    #compile only the atm model
    newline && time_elapsed_min
    echo  "Compiling just the atm model"
    cd $bld_dir/cmake-bld/cmake/atm/
    make -j40
    status=`echo $?` #DO NOT MOVE THIS LINE:Status should be captured just after make

    newline && time_elapsed_min
    #if the model failed to build, exit
    if [ $status -ne 0 ]; then
        newline && time_elapsed_min
        echo "Model Failed to compile with status:$status, exiting..."
        newline
        exit 1
    fi
    newline && time_elapsed_min

    #now link to create the exe
    echo "Linking $executable"
    cd $bld_dir/cmake-bld/cmake/cpl
    cmake -E cmake_link_script CMakeFiles/e3sm.exe.dir/link.txt
    newline && time_elapsed_min

    #see if the exe has been created
    if_file_not_present_then_exit $bld_dir/$executable " (Please look into bldlogs to see what went wrong or) issue ./case.build && ./case.submit to rebuild the model"

    newline && time_elapsed_min

    #change BUILD_COMPLETE to True
    cd $case_dir
    echo "Setting BUILD_COMPLETE to TRUE"
    ./xmlchange BUILD_COMPLETE=TRUE

    #submit run
    newline && time_elapsed_min
    echo "Submit run to the compute queue"
    cd $case_dir
    ./case.submit

    newline && time_elapsed_min
}

#---------------------
# Function Definitions
#---------------------

if_file_not_present_then_exit () {
    if [ ! -f $1 ]; then
        echo "$1 is not found.$2"
        exit 1
    fi
}

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

#capture start time to compute time elapsed in minutes
starttime=$(date +%s)
echo 'Start Time:' $(date +%T)
#call the main function
main
