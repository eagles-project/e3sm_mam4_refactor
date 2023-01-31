#!/bin/sh


main() {
    #-------------------------
    #User input Starts
    #-------------------------

    #subroutine name
    sub_name=clddiag

    #module name
    module_name=wetdep

    #file path
    dir_path=components/eam/src/chemistry/aerosol/

    #CPP directive to turn on file writing
    cpp_directive=YAML_WETDEP

    #-------------------------
    #USER INPUT ENDS
    #-------------------------

    #----------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------


    #ensure that script is run in the right directory (i.e. "e3sm_mam4_refactor") as
    #relative path code works only if script is run at the code root
    this_dir=${PWD##*/}

    if [ $this_dir != "e3sm_mam4_refactor" ]; then
        echo 'Please run this script at source code root in directory "e3sm_mam4_refactor"'
        exit 1
    fi

    #Yaml directory path
    yaml_dir_path=components/eam/src/chemistry/yaml

    #path to directory containing beg and end input stubs
    stub_dir=mam4_refactor_scripts



    #relative path
    relative_path=`realpath --relative-to="$dir_path" "$yaml_dir_path"`

    #create directory and files if not already created
    /bin/mkdir -p $yaml_dir_path/$module_name/f90_yaml

    #code to add
    echo '-------------'
    echo 'Code to add:'
    echo '-------------'

    echo '(If not already added, add the following line at module level at the top)'
    echo "#include \"$relative_path/common_files/common_uses.ymlf90\""

    #Subroutine beginning file
    newline
    create_file $module_name $sub_name $yaml_dir_path $relative_path "beginning"

    #subroutine end file
    newline
    create_file $module_name $sub_name $yaml_dir_path $relative_path "end"
}

#---------------------
# Function Definitions
#---------------------

#Generate a newline
newline () { echo ''; }

create_file () {

    #create_file $module_name $sub_name $yaml_dir_path $relative_path "end"
    echo "(Add the following line at the $5 of the subroutine)"

    #find beginning or end sub string
    if [ $5 == "beginning" ]; then
        sub_str=beg
    elif [ $5 == "end" ]; then
        sub_str=end
    else
        echo 'Invalid string (beginning or end)for function create_file'
        exit 1
    fi

    #full file name
    f_name=${sub_name}_${sub_str}_yml.f90

    #code to include
    echo "#include \"$4/$1/f90_yaml/$f_name\""

    #create this file if it doesn't exist
    f_path=$3/$1/f90_yaml/$f_name
    if test -f "$f_path"; then
        newline
        echo "[INFO ONLY]: $f_path already exists."
    else
        newline
        echo '[INFO ONLY]:Creating NEW file: '$f_path
        f_tmp=tmp.inp
        sed s/"!#ifdef YAML_CPP"/"#ifdef $cpp_directive"/g $stub_dir/${sub_str}.inp > $f_tmp
        sed -i s/"SUB_NAME"/"\'$2\'"/g $f_tmp
        cat $f_tmp>>$f_path
        /bin/rm $f_tmp
    fi
}

#call the main function
main
