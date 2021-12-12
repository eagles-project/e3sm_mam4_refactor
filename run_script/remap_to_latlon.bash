#!/bin/bash

module load nco

   map_file="/qfs/people/zender/data/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"
  data_path="/compyfs/wanh895/scidac_scratch/macmic_subcycles_compy_F2010_ne30pg2_r05_oECv3_CondiDiag_for_MichaelBrunke/run/"
 input_file="macmic_subcycles_compy_F2010_ne30pg2_r05_oECv3_CondiDiag_for_MichaelBrunke.eam.h0.0001-01-05-00000.nc"
output_file="macmic_subcycles_compy_F2010_ne30pg2_r05_oECv3_CondiDiag_for_MichaelBrunke.eam.h0.0001-01-05-00000.latlon.nc"

ncremap  -m ${map_file} -i ${data_path}${input_file} -o ${data_path}${output_file}


