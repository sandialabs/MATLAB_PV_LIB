def rlm(y,x,inter):

    import statsmodels.api as sm
    import numpy as np
    
    if inter:
        x = sm.add_constant(x)
        
    rlmod = sm.RLM(y,x,M= sm.robust.norms.TukeyBiweight())
    res = rlmod.fit()
    return res.params

    
import numpy as np
import csv
import os

wdir = os.getenv('TEMP')

os.chdir(wdir)

f1 = open('temp_rlm_config.txt')
intercept = int(f1.readline()); # 0 if no intercept is to be used
f1.close()

dat = np.loadtxt('temp_rlm_input.csv',delimiter=',')

X = dat[:,0:-1]
Y = dat[:,-1]

pp = rlm(Y,X,intercept)

np.savetxt('temp_rlm_output.csv', pp)
