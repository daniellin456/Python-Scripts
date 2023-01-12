#!/bin/bash

original_G0_factor=1.0
modified_G0_factor=1.0

ramses_output="/data/daniellin/RAMSES_RUNS/HEATING_COOLING/TR_FLOOR_10K/ISM_FFSCT_0.100_10000M/"
figure_output="/data/daniellin/RAMSES_RUNS/HEATING_COOLING/TR_FLOOR_10K/ISM_FFSCT_0.100_10000M/profiles"
source /home/ynlee/softs/envpym

# Modify the G0 factor in calculate_theoretical_temperature file and Store the nH and equilibrium_temp
sed -i '5s/'$original_G0_factor'/'$modified_G0_factor'/' /data/daniellin/PYTHON_SCRIPTS/COOLING_FUNCTION_ANALYZE/EQUILIBRIUM/equilibrium_temp.py
python3 /data/daniellin/PYTHON_SCRIPTS/COOLING_FUNCTION_ANALYZE/EQUILIBRIUM/equilibrium_temp.py

python write_pymsesrc.py $ramses_output
python produce_maps.py $ramses_output
python histogram.py -i $ramses_output -o $figure_output


