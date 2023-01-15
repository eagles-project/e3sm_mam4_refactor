#!/bin/sh


main() {
    #User input

    #module name
    module_name=module_aero_amicphys

    #subroutine name
    sub_name=mam_pcarbon_aging_1subarea

    #file path
    dir_path=components/eam/src/chemistry/modal_aero/

    #USER INPUT ENDS

    yaml_dir_path=components/eam/src/chemistry/yaml

    #ensure that script is run in the right directory (i.e. "e3sm_mam4_refactor") as
    #relative path code works only if script is run at the code root
    this_dir=${PWD##*/}

    if [ $this_dir != "e3sm_mam4_refactor" ]; then
        echo 'Please run this script at source code root in directory "e3sm_mam4_refactor"'
    fi


    #relative path
    relative_path=`realpath --relative-to="$dir_path" "$yaml_dir_path"`

    #create directory and files if not already created
    /bin/mkdir -p $yaml_dir_path/$module_name/f90_yaml

    #code to add
    echo '-------------'
    echo 'Code to add:'
    echo '-------------'

    echo '(If not already aded, add this at module level at the top)'
    echo "#include \"$relative_path/common_files/common_uses.ymlf90\""

    newline
    create_file $module_name $sub_name $yaml_dir_path $relative_path "beginning"

    newline
    create_file $module_name $sub_name $yaml_dir_path $relative_path "beginning"

    echo '(Add this at the beginning of the subroutine)'
    f_beg_name=${sub_name}_beg.ymlf90
    echo "#include \"$relative_path/$module_name/f90_yaml/$f_beg_name\""
    #create this file if it doesn't exist
    f_beg_path=$yaml_dir_path/$module_name/f90_yaml/$f_beg_name
    if test -f "$FILE"; then
        echo "$FILE exists."
    fi
    /bin/touch $yaml_dir_path/$module_name/f90_yaml/$f_beg_name

    newline
    echo '(Add this at the end of the subroutine)'
    f_end_file=${sub_name}_end.ymlf90
    echo "#include \"$relative_path/$module_name/f90_yaml/$f_end_file\""
    /bin/touch $yaml_dir_path/$module_name/f90_yaml/$f_end_name
    #create directory (if not created already) and beg and end files
}

#---------------------
# Function Definitions
#---------------------

#Generate a newline
newline () { echo ''; }
