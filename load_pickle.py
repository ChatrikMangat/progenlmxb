# This is a small example to navigate through one set of simulation data (corresponding to a single Accretor mass) and check if all files are present. 

import numpy as np
import pickle
import time
import os.path

start_time = time.time()

print("Start Test")

m_arrl = np.arange(4.0,7.1,0.1)
m_arrs = np.arange(0.95,4.05,0.05)
p_arrl = np.arange(1.65,4.05,0.05)
p_arrs = np.arange(-0.60,1.66,0.02)

confs = {'lmlp':[m_arrl,p_arrl],'smlp':[m_arrs,p_arrl],'lmsp':[m_arrl,p_arrs],'smsp':[m_arrs,p_arrs]}

for c in confs.keys():
    for i in confs[c][0]:
        for j in confs[c][1]:
                try:

                    fpath = '/Path/to/simulation/set/'+c+'/m_'+f'{i:4.2f}'+'_p_'+f'{j:4.2f}'+'.data'
                    if os.path.exists(fpath):
                        infile = open(fpath,'rb')
                        vals = pickle.load(infile)
                        infile.close()

                        print(fpath)
                        print(vals.shape)
                    else:
                        print("No path found: "+fpath)
                    
                    # Any code written here will be applied to all simulation data files present in the set.

                except Exception as e:
                    print("Error occurred in "+fpath+"\n"+str(e)+"\n")

print("Test Complete, time taken: "+str(time.time() - start_time))
