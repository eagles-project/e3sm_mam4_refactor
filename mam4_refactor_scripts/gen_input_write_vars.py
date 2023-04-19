#!/bin/python
import argparse

#-------------------------------------------------------------------
# Purpose: Given a list of intent-ins, inouts and outs, this script
# --------
#         print calls for writing variables for the yaml/python I/O
#         generation code
#-------------------------------------------------------------------

def main():

    #parse command line arguments
    ins,outs = parse_args()

    #------------------------------------------------------
    #Provide a string of inputs as comma-separted variables
    #
    # For example:
    # ins = 'a,b,c,d,e,f'
    #
    # The script will automatically extract variables from
    # the string above
    #------------------------------------------------------
    #ins  = 'ncnst, aer'

    #------------------------------------------------------
    #Provide a string of outputs as comma-separted variables
    # (similar to inputs)
    #------------------------------------------------------
    #outs = 'total_interstitial_aer_num'

    #------------------------------------------------------
    #User input ends
    #------------------------------------------------------


    #Break the inputs string into individual variables
    in_vars = [el.strip() for el in ins.split(',')]
    out_vars = [el.strip() for el in outs.split(',')]



    #Generate write line fortran codes
    #input calls
    print()
    print('Input calls to copy-paste:')
    print(27*'-')
    write_calls(in_vars)

    #output calls
    print('')
    print('Output calls to copy-paste:')
    print(28*'-')
    write_calls(out_vars,'output')

# Rest of the functions

def write_calls(var_list, out=None):

    out_unit = ''
    #since the input file write both yaml input and python output files,
    # we need to include "unit_output" as an argument for input calls only
    if not out: #i.e. if it is a call for printing intent-ins or intent-inouts
        out_unit = 'unit_output,' # ensure to put comma at the end

    for var in var_list:
        print(f"call write_var(unit_input,{out_unit}'{var}',{var})")
    print('')

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help="A comma separated list of intent-in/inout variables", required=True)
    parser.add_argument('-o', type=str, help="A comma separated list of intent-out/inout variables", required=True)

    args = parser.parse_args()

    #return the args captured by -i and -o flags
    return args.i, args.o


main()
