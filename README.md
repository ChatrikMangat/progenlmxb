# progenlmxb

This is a companion repository for the database and tool described in "" 

**directory_structure**: Contains scripts to create the all the pre-set directories and input files for one set of simulations (corresponding to one accretor mass). There are four directories to segregate simulations based on the input grid resolution as described in our paper. This script can be run using the following command in each of the four directories:
```
./loop.sh
```

**MESA_lmxb_src**: Contains the extras files required to compile a binary executable in MESA (v15140) which has access to Convection And Rotation Boosted (CARB) Magnetic Braking.

**load_pickle.py**: This script is an example which navigates serially through one set of simulations in the database. Any changes required in the data files can be made using this template. After making changes to internal path variables for the database, the script can be run using:
```
python3 load_pickle.py
``` 

**PEAS_tool**: Contains the python script and sample inputs for the Progenitor Extractor for Accreting Systems (PEAS) tool. After making changes to internal path variables for the database, the script can be run using: 
```
python3 query.py sample_query.txt
```
The output containing progenitor information will be stored in the same directory as the query file. The input and output paths can be edited in the query.py file if needed. 

The output is structured in the order: initial donor mass (Msol), log10(initial orbital period/days), initial accretor mass (Msol), amount of time spent as the observed system (years), total evolution time spanned by the simulation (years), donor mass at the start of MT (Msol), log10(orbital period at the start of MT/days), accretor mass when system first matches query (Msol), accretor mass when system last matches query (Msol), age of system when it first matches query (years) and age of system when it last matches query (years). 

If the tool runs without any errors and the output file contains only column headers, it indicates that no progenitor was found in the database.
