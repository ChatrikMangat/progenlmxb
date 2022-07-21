import numpy as np
import pickle
import time
import os
import sys


## Arrays for creating Progenitor Grid. The m_arr arrays corrrespond to the progenitor donor masses and the p_arr arrays correspond to log10(progenitor orbital period/days). The "l" and "s" signify large and small interval sizes respectively.
m_arrl = np.arange(4.0,7.1,0.1)
m_arrs = np.arange(0.95,4.05,0.05)
p_arrl = np.arange(1.65,4.05,0.05)
p_arrs = np.arange(-0.60,1.66,0.02)

## List of black hole masses, each mass corresponds to one Set of simulation data. 
bh_masses = [5.0,7.0,10.0]

## Creating progenitor grid, which is divided into four parts based on grid density. One simulation data Set (containing all simulations for one accretor mass) is also divided into four folders:
## lmlp -> large mass intervals, large period intervals
## smlp -> small mass intervals, large period intervals 
## lmsp -> large mass intervals, small period intervals 
## smsp -> small mass intervals, small period intervals
confs = {'lmlp':[m_arrl,p_arrl],'smlp':[m_arrs,p_arrl],'lmsp':[m_arrl,p_arrs],'smsp':[m_arrs,p_arrs]}

## m1 -> Donor Mass {Msol}
## m2 -> Accretor Mass {Msol}
## mt -> log10(MT Rate) {Msol/yr}
## p -> Orbital Period {days}
## teff -> log10(Donor Effective Temp) {K}
labels = ['m1','m2','mt','p','teff']

## This list stores all the errors encountered during runtime.
errors = []

## Slightly modified binary search algorithm used to search for the queried upper and lower donor mass limits in the simulation data.
def search(array, element):
    
    mid = 0
    start = 0
    end = len(array)
    step = 0

    while (start <= end):
        step = step+1
        mid = (start + end) // 2
        
        if mid == len(array):
            return len(array)-1

        if element == array[mid]:
            return mid

        if element > array[mid]:
            end = mid - 1
        else:
            start = mid + 1
    return start


## Function to retrieve query information from .txt file and return a dictionary for the query. If the query format is wrong or any property limit is out of bounds, the function appends the error to the "errors" list and returns -1.
def get_query(qpath):
    
    try:
        q = open(qpath,"r")

        ## First character of query.txt file is stored in bhns, this determines whether to search in NS simulations or BH simulations.
        bhns = (q.readline()).split(None, 1)[0]
        q.close()

        qtemp = np.loadtxt(qpath,delimiter=",",skiprows=1)

        ## Take limits from query.txt file and create a dictionary containing the bhns variable and the limits for each property. Exceptions happen when the input data is not in the correct format.
        #print(qtemp.shape)
        query = {'bhns': int(bhns)}

        for i in range(qtemp.shape[0]):
            query[labels[i]] = qtemp[i]

    except Exception as e:
        errors.append("Error in query input, please recheck format\n"+str(e))
        return -1

    ## Various checks to determine whether the input limits are within bounds or not. It is also checked that the upper limit in query is larger than the lower limit. If any problem is found, return -1.
    if (query['m1'][0] < 0.0 or query['m1'][0] > query['m1'][1]):
        errors.append("Wrong Donor Mass range")
        return -1
    if (query['m2'][0] < 0.0 or query['m2'][0] > query['m2'][1]):
        errors.append("Wrong Accretor Mass range")
        return -1
    if (query['mt'][0] < -100.0 or query['mt'][0] > query['mt'][1]):
        errors.append("Wrong MT Rate range")
        return -1
    if (query['p'][0] < 0.0 or query['p'][0] > query['p'][1]):
        errors.append("Wrong Orbital Period range")
        return -1
    if query['teff'][0] > query['teff'][1]:
        errors.append("Wrong Donor Teff range")
        return -1

    ## If query satisfies all checks, return dictionary containing query info. 
    return query


## Function to check whether a model in the simulation (data) matches all the queried properties (query). Returns 1 if simulated system fully matches query and 0 if it does not.
def match_props(data,query):
    
    try:
        if 'm1' in query:
            if (data[0] < query['m1'][0] or data[0] > query['m1'][1]):
                #print("No m1")
                return 0
        if 'm2' in query:
            if (data[1] < query['m2'][0] or data[1] > query['m2'][1]):
                #print("No m2")
                return 0
        if 'mt' in query:
            if (data[2] < query['mt'][0] or data[2] > query['mt'][1]):
                #print("No mt")
                return 0
        if 'p' in query:
            if (data[3] < query['p'][0] or data[3] > query['p'][1]):
                #print("No p")
                return 0
        if 'teff' in query:
            if (data[4] < query['teff'][0] or data[4] > query['teff'][1]):
                #print("No teff")
                return 0
    except Exception as e:
        errors.append("Error in matching properties")
        return 0

    return 1


## Function to find the array index where the simulated system first starts MT.
def find_mt_start(mt_arr):
    
    idx_start = 0

    if len(mt_arr) <= 5:
        return idx_start

    ## If there are four entries with log10(MT) > -15, it is considered that the system has started MT. 
    for i in range(len(mt_arr)-3):
        if (mt_arr[i] >= -15 and mt_arr[i+1] >= -15 and mt_arr[i+2] >= -15 and mt_arr[i+3] >= -15):
            idx_start = i
            break
    
    return idx_start


## Function to look through data files to find progenitors for the query
def get_progens(query):
    
    progens = []
    
    data_paths = []
    ## If bhns is 0, we only search through the NS simulation Set. If bhns is 1, we choose which accretor masses will be relevant for the search. Based on those accretor masses, we choose the BH simulation Sets to be searched through.
    ## For example, if the query given is [9.0,9.5] for accretor mass, the code will not go into the 10 Msol Set at all, since there will be no matches there.
    ## The second condition for choosing the relevant Set  uses the 3.5 Msol quantity. So for our example, in case of the 5 Msol Set, the code will check whether 5+3.5=8.5 is less than the lower limit of the query (9.0). If that is the case, it will not search in the 5 Msol Set.
    ## This condition is used because the maximum mass that can be accreted by a BH in our simulation setup is 0.5*7.0 (MT efficiency * max donor mass) = 3.5
    ## No system in the 5 msol Set will reach the accretor mass given in the query, so there is no need to search through the 5 Msol Set.
    ## Therefore, for our example, only the 7 Msol Set will be selected to be searched through.
    
    db_location = os.environ['DB_LOCATION']
    print(db_location)

    if query['bhns'] == 0:
        data_paths.append(db_location+'ns_data/')
    else:
        for x in bh_masses:
            if(x > query['m2'][1]):
                continue
            if(x < query['m2'][0]-3.5):
                continue
            data_paths.append(db_location+'runs'+str(int(x))+'_data/')
    print(data_paths)

    ## Loop through the selected simulation Sets
    for dpath in data_paths:
        ## Loop through the four folders (lmlp,smlp,lmsp,smsp) in each Set
        for c in confs.keys():
            ## Loops for initial mass and orbital period values, these will go through each simulation data file in the folder serially. 
            for i in confs[c][0]:
                for j in confs[c][1]:
                    ## Any code inside this loop will be applied to all simulation data files present in all the selected Sets. 

                    ## Only go through simulations where initial donor mass is more than the lower limit of donor mass in query.
                    if (i >= query['m1'][0]):
                        try:
                            fpath = dpath+c+'/m_'+f'{i:4.2f}'+'_p_'+f'{j:4.2f}'+'.data'
                            ## Check if simulation data file exists
                            if os.path.exists(fpath):
                                infile = open(fpath,'rb')
                                vals = pickle.load(infile)
                                infile.close()
                                
                                #print(fpath)
                                #print(vals.shape)
                                
                                ## Only go through simulations where final donor mass is less than upper limit of donor mass in query. Similarly, only go through simulations where final accretor mass is more than lower limit of accretor mass in query.
                                if (vals[0][-1] <= query['m1'][1] and vals[1][-1] >= query['m2'][0]):
                                    ## Variables to indicate if a match is found between the query and the simulation data. pflag indicates whether any match was found in the simulation data file or not, it is set only once (at max) per simulation file. flag is used to check every model one by one.
                                    flag = 0
                                    pflag = 0
                                    ## Initialize the variables that define the start and end indices of the Window of models in simulation data array which satisfy query donor mass limits.
                                    start_m = 0
                                    end_m = 0
                                    
                                    ## If upper limit of donor mass in query is greater than initial donor mass, start Window at first array element (first model in simulation).
                                    if vals[0][0] <= query['m1'][1]:
                                        start_m = 0
                                    else:
                                    ## Search through the simulation data and check if any model has donor mass equal to upper donor mass limit in query. If found, save array index of the model in start_m (start Window at that model). 
                                        start_m = search(vals[0],query['m1'][1])-1
                                    
                                    ## If lower limit of donor mass in query is less than final donor mass, end Window at the last array element (last model in simulation).
                                    if vals[0][-1] >= query['m1'][0]:
                                        end_m = vals.shape[1]-1
                                    else:
                                    ## Search through the simulation data and check if any model has donor mass equal to lower donor mass limit in query. If found, save array index of the model in end_m (end Window at that model).
                                        end_m = search(vals[0],query['m1'][0])
                                    
                                    ## If a valid donor mass Window is found (i.e. 0 <= start_m <= end_m <= array length), go through all models in Window one by one, and check if all model properties match the queried properties.
                                    if (start_m <= end_m):
                                        #print("Mass matched ",start_m,end_m,end_m-start_m)
                                        ## obs_time is the amount of time spent by the simulated system while it matches the queried system.
                                        ## mt_start is the array index of the model where the simulation starts significant MT.
                                        obs_time = 0.0
                                        mt_start = 0

                                        for k in range(start_m,end_m+1):
                                            if (k < 0 or k >= vals.shape[1]):
                                                break
                                            ## Check if all model properties match the queried properties, if so, set flag to 1. 
                                            flag = match_props([vals[0][k],vals[1][k],vals[2][k],vals[3][k],vals[4][k]],query)
                                            
                                            ## If the model matches the query, first check if this is the first match in the simulation data file. If so, set pflag to 1 and save age and accretor mass of the model. After that check, add the timestep (dt) of the model to obs_time (time spent by the simulation while it matches the queried system).
                                            ## macc_max and age_max will keep getting overwritten till the last model that matches the query. Hence, those will finally store the age and accretor mass of the last model that matches the query.
                                            if (flag == 1):
                                                if (pflag == 0):
                                                    macc_min = vals[1][k]
                                                    age_min = vals[5][k]-10**vals[7][k]
                                                    pflag = 1
                                                
                                                macc_max = vals[1][k] 
                                                age_max = vals[5][k]
                                                obs_time = obs_time + 10**vals[7][k]

                                        ## If a match was found in the data, findthe array index of the model where significant MT started and append all output properties to list of progenitors.
                                        ## The output properties are: initial donor mass (Msol), log10(initial orbital period/days), initial accretor mass (Msol), amount of time spent as the observed system (years), total evolution time spanned by the simulation (years), donor mass at the start of MT (Msol), log10(orbital period at the start of MT/days), accretor mass when system first matches query (Msol), accretor mass when system last matches query (Msol), age of system when it first matches query (years) and age of system when it last matches query (years).
                                        if (pflag == 1):
                                            mt_start = find_mt_start(vals[2])
                                            tot_time = vals[5][-1]-vals[5][0]
        
                                            progens.append([i,j,vals[1][0],obs_time,tot_time,vals[0][mt_start],np.log10(vals[3][mt_start]),macc_min,macc_max,age_min,age_max])
                                            print("Found Progenitor -> "+fpath)
                                
                            else:
                                #print("No path found: "+fpath)
                                continue

                        except Exception as e:
                            errors.append("Error occurred in "+fpath+"\n"+str(e)+"\n")
    return progens


## Function to store progenitor information for the given query in a .txt file
def output_progens(rpath,progens):
    
    try:
        ## Open output file and enter the progenitor information.
        outfile = open(rpath,'w')
        
        outfile.write('{:8}'.format('M_don,i')+'{:12}'.format('lg[P_orb,i]')+'{:8}'.format('M_acc,i')+'{:16}'.format('tau_obs')+'{:16}'.format('tau_tot')+'{:9}'.format('M_don,mt')+'{:13}'.format('lg[P_orb,mt]')+'{:16}'.format('M_acc,obs start')+'{:14}'.format('M_acc,obs end')+'{:16}'.format('age_obs start')+'{:15}'.format('age_obs end')+'\n') 
        outfile.write('{:8}'.format('[Msun]')+'{:12}'.format('[P in days]')+'{:8}'.format('[Msun]')+'{:16}'.format('[years]')+'{:16}'.format('[years]')+'{:9}'.format('[Msun]')+'{:13}'.format('[P in days]')+'{:16}'.format('[Msun]')+'{:14}'.format('[Msun]')+'{:16}'.format('[years]')+'{:15}'.format('[years]')+'\n')
        for i in progens:
            s = f'{i[0]:<8.2f}'+f'{i[1]:<12.3f}'+f'{i[2]:<8.2f}'+f'{i[3]:<16.3f}'+f'{i[4]:<16.3f}'+f'{i[5]:<9.4f}'+f'{i[6]:<13.4f}'+f'{i[7]:<16.4f}'+f'{i[8]:<14.4f}'+f'{i[9]:<16.3f}'+f'{i[10]:<15.3f}'+'\n'
            outfile.write(s)

        outfile.close()
        print("Progenitor properties stored in: "+rpath)
    
    except Exception as e:
        errors.append("Error occurred in writing output:\n"+str(e))
        error_exit(rpath)


## This function is used to write any fatal error into the output file before the code shuts down. It prints out all the errors stored in the "errors" list and exits the code.
def error_exit(rpath):
    
    outfile = open(rpath,'w')
    outfile.write("Errors:\n")
    for i in errors:
        outfile.write(i+'\n')

    print("Error occurred: Check progens file for details")
    exit()


## Driver Code
print("Start Search")
start_time = time.time()

## Query .txt filename is given as command line argument, change qpath according to where file is stored.
qname = sys.argv[1]
qpath = ""+qname

## Progenitor list is stored in a .txt file corresponding to the query file, change rpath according to preferences.
rpath = "progens_"+qname

query = get_query(qpath)

## If query was not generated properly, print errors into output file and exit code.
if (query == -1):
    error_exit(rpath)

print("Query input taken:")
print(query)

## If query is succesfully generated, search for progenitors of the queried system.
progens = get_progens(query)

## If no errors were encountered during the search, create output file.
output_progens(rpath,progens)

print("Search Complete, time taken: "+str(time.time() - start_time)+" seconds")
