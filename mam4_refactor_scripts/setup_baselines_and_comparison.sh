#!/bin/sh
main() {


#-------------------
#user input starts
#-------------------

#Tag for the test
tag="tag1"

#baseline tag (must be different from the tag above; you can simply add "base_" infront of tag if you like)
baseline_tag="base_tag1"

#Tag for the comparison simulation(must be different from the tags above; you can simply add "comp_" infront of tag if you like)
comparison_tag="comp_tag1"

#--------------------------------------
#Less frequently changed user input
#--------------------------------------5B

test_name="SMS_D_Ln5_P32x1"
grid="ne4_oQU240"
compset="F2010"
project="esmd"

#-------------------
#user input ends
#-------------------

#======#======#======#======#======#======#======#======#======#======#======#======#======#======#======#======#======

status_file="TestStatus"

#check if we are in the right directory and the script exists
cd cime/scripts
if [ ! -f "create_test" ]; then
    echo "create_test does not exist. Have you initialized submodules? "
    exit -1 #exit if it is wrong directory or file doesn't exist
fi

newline
print_time
#run baseline test
echo 'Generating baselines...'
generate_output=$(./create_test $test_name.$grid.$compset --compiler intel -t $tag -g -b $baseline_tag -p $project -o  2>&1 & )

#see if any error is returned
error_check "$generate_output"


#capture baseline directory
generate_case_dir=`echo "$generate_output" |grep "Creating test directory"| cut -d " " -f 4`

newline
print_time
#Launch the comparison simulation but DO NOT submit it
echo 'Compare baselines build starts...'
compare_output=$(./create_test $test_name.$grid.$compset --compiler intel -t $comparison_tag -c -b $baseline_tag -p $project --no-run  2>&1 & )

#see if any error is returned
error_check "$compare_output"

#capture compare directory
compare_case_dir=`echo "$compare_output" |grep "Creating test directory"| cut -d " " -f 4`

newline
print_time
#Find out if the baseline run is complete:
echo 'Check if base;ines has been generated...'
cd $generate_case_dir

baseline_generated=0

while [ $baseline_generated -eq 0 ]
do
    #see if TestStatus file exists
    if [ ! -f "$status_file" ]; then
        echo '$status_file doesnt exists, will check back in 30 seconds... '
        sleep 30
        continue
    fi

    baseline_generated=`cat $status_file |grep GENERATE |grep PASS|wc -l`
    sleep 30
    echo 'Checking if baselines have been generated:'$baseline_generated
done

newline
print_time
echo 'Baselines have been generated, submiting comparison run..'

cd $compare_case_dir

compare_output=$(./ase.submit 2>&1 & )
newline
print_time


}

#functions

#generate a newline
newline () { echo ''; }

#print time
print_time() { date +%T ; }


#check errors in command outputs and exit the script
error_check() {
    error_count=`echo "$1" |grep -i "Error"|wc -l`
    if [ $error_count -ne 0 ]; then
        echo "$1"
        newline
        echo 'Error occured, exiting..'
        exit -1
    fi
}

main
