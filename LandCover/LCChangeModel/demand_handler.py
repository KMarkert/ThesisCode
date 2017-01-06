import numpy as np

def get_demand(preLC):
    
    raterngs = {0:[0,2],1:[-1,0],2:[-2,2],3:[-1,1],4:[-1,1],5:[-0.5,0]}
    
    rates = []
    
    for i in range(6):
        rng = raterngs[i]
        rate = (rng[0]-rng[1]) * np.random.random(1) + rng[1]
        rates.append(rate)
        
    demand = []
    
    for i in range(6):
        idx = np.where(preLC == i)
        cnt = len(idx[0])
        
        d = int(cnt * (1+(rates[i]/100.)))
        
        demand.append(d)
        
    return demand
    
def rate_check(predLC,demand):
    
    check = []
    
    for i in range(6):
        
        idx = np.where(predLC==i)
        
        pd = float(len(idx[0])-demand[i])/float(demand[i])
        
        if pd >= -0.05 and pd <= 0.05:
            check.append(True)
        else:
            check.append(False)
            
    return np.array(check)
        
def update_itervar(predLC,demand,itervar,check):
    
    for i in range(6):
        if check[i] == False:
            idx = np.where(predLC==i)
            pd = -1 * (float(len(idx[0])-demand[i])/float(demand[i]))
            
            itervar[i] = itervar[i] + pd
            
    return itervar