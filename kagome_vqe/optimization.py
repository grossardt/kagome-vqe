import time
import numpy as np
from datetime import datetime

"""Optimization functions for VQE."""

def next_index(indices):
    """Get next index to optimize"""
    index = indices.pop(0)
    indices.append(index)
    return index

def parameter_rotate(theta):
    """rotate parameter vector by one position of the kagome lattice"""
    new_edge_indices = (2,3,4,5,6,7,8,9,10,11,0,1,13,14,15,16,17,12)
    new_theta = []
    for i in new_edge_indices:
        new_theta.append(theta[3*i])
        new_theta.append(theta[3*i+1])
        new_theta.append(theta[3*i+2])
    return new_theta

def optimization_step(cost_function, theta, index, c0 = None):
    """Perform one optimization step with Rotosolve"""
    # cost function at old theta values if not provided
    if c0 is None:
        c0 = cost_function(theta)
    # new theta vectors with theta_i +/- pi/2
    tm = theta[:index]+[theta[index]-np.pi/2]+theta[index+1:]
    tp = theta[:index]+[theta[index]+np.pi/2]+theta[index+1:]
    # cost function at +/- pi/2
    cm, cp = cost_function(tm, tp)
    dp = c0 - cp
    dm = c0 - cm
    # update theta
    theta[index] = (theta[index] + np.arctan2(cp-cm,dp+dm) + np.pi) % (2*np.pi)
    # return new cost function value (next c0)
    new_cost = (cp+cm)/2 - np.sqrt((dp**2+dm**2)/2)
    return (new_cost, c0-new_cost)

def save_step(filename, step, index, cost, change, theta):
    """Save optimization step to file"""
    with open(filename, 'a') as f:
        f.write('STEP {}, varying index {}\n'.format(step,index))
        f.write('  COST:   {}\n'.format(cost))
        f.write('  CHANGE: {}\n'.format(change))
        f.write('  PARAMS: {}\n'.format(theta))
        f.close()
        
def optimize(cost_function, theta, max_steps, stopping_delta, outfile, step = 0, cost = None):
    """Run the optimization with Rotosolve until stopping criterion reached"""
    n_params = len(theta)
    min_steps = n_params
    last_change = 0 # last step with minimum required change
    # optimization loop:
    start_time = time.process_time()
    indices = [39,40,41,45,46,47,51,52,53,36,37,38,42,43,44,48,49,50,3,4,5,9,10,11,15,16,17,
               21,22,23,27,28,29,33,34,35,0,1,2,6,7,8,12,13,14,18,19,20,24,25,26,30,31,32]
    
    if cost is None:
        print('Calculating initial cost...')
        cost = cost_function(theta)
        print('cost value: {}'.format(cost))
    
    # make sure indices matches last value if continuing previous calculation
    for _ in range(step):
        index = next_index(indices)
    
    while step < max_steps:
        index = next_index(indices)
        cost, change = optimization_step(cost_function, theta, index, cost)
        print(datetime.now().strftime('%Y%m%d-%H%M%S: '), end='')
        print('Saving iteration {}, cost value: {}'.format(step,cost))
        save_step(outfile, step, index, cost, change, theta)
        if change < stopping_delta:
            if step >= max(min_steps - 1,
                           last_change + min_steps):
                break
        else:
            last_change = step
        step += 1
    end_time = time.process_time()
    print('Optimization finished in {} steps and {} seconds'.format(step,end_time-start_time))
    print('with cost function value {}'.format(cost))
