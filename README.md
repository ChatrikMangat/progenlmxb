# progenlmxb

This is a companion repository for the database and tool described in "" 

**directory_structure**: Contains scripts to create the all the pre-set directories and input files for one set of simulations (corresponding to one accretor mass). There are four directories to segregate simulations based on the input grid resolution as described in "". This can be done using the following command in each of the four directories:
```
./loop.sh
```

**MESA_lmxb_src**: Contains the extras files required to compile a binary executable in MESA which has access to Convection And Rotation Boosted (CARB) Magnetic Braking.

**tool**: Contains the python script and sample inputs for the Progenitor Extractor for Accreting Systems (PEAS) tool. After making changes to internal path variables for the database (Zenodo link), the script can be run using: 
```
python3 query.py sample_query.txt
```

**load_pickle.py**: This script is an example which navigates serially through one set of simulations in the database (Zenodo link). Any changes required in the data files can be made using this template. After making changes to internal path variables for the database, the script can be run using:
```
python3 load_pickle.py
``` 
