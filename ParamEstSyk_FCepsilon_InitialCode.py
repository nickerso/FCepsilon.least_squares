import numpy as np

from matplotlib import pyplot as plt
from scipy.optimize import curve_fit, minimize,least_squares

import OpenCOR as oc

# Experiment data -- has to be regularly sampled
expt_data = np.loadtxt('FCepsilonRI_trial_FullResolution.csv', delimiter=',')
expt_time = expt_data[...,0]
expt_pGRB2 = expt_data[...,1]
expt_pSyk = expt_data[...,2]

# The state variable in the model that the data represents
expt_state_uri_pGrb2 = 'FCepsilonRI/pGrb2'
expt_state_uri_pSyk = 'FCepsilonRI/pSyk'

# Load and initialise the simulation
simulation = oc.openSimulation('FCepsilonRI_trial.cellml')

# In case we have reloaded an open simulation
simulation.resetParameters()
simulation.clearResults()

# Reference some simulation objects
initial_data = simulation.data()
constants = initial_data.constants()
results = simulation.results()

# Simulation time points must match experiment data points
initial_data.setStartingPoint(0.0)
initial_data.setEndingPoint(3600)
initial_data.setPointInterval(1)

# Specify as two parallel lists:
# 1. Uri's of parameters to estimate
# 2. Initial values of parameters
parameter_names = list(constants.keys())
#parameter_names.remove('FCepsilonRI/k_f1')
#parameter_names.remove('FCepsilonRI/K_1')
#parameter_names.remove('FCepsilonRI/k_f2')
#parameter_names.remove('FCepsilonRI/K_2')
#parameter_names.remove('FCepsilonRI/K_3')
#parameter_names.remove('FCepsilonRI/k_f5')
#parameter_names.remove('FCepsilonRI/k_r5')
#parameter_names.remove('FCepsilonRI/k_f6')
#parameter_names.remove('FCepsilonRI/k_r6')
#parameter_names.remove('FCepsilonRI/K_7')
#parameter_names.remove('FCepsilonRI/V_7')
#parameter_names.remove('FCepsilonRI/Lyn')

parameter_names.remove('FCepsilonRI/k_f3')
parameter_names.remove('FCepsilonRI/k_f4')
parameter_names.remove('FCepsilonRI/k_r4')
initial_params = [constants[name] for name in parameter_names]

# Set bounds for parameters (optional)
parameter_bounds = [len(initial_params)*[0], len(initial_params)*[6]]
parameter_bounds[0][parameter_names.index('FCepsilonRI/k_f1')] = 0.00023
parameter_bounds[1][parameter_names.index('FCepsilonRI/k_f1')] = 70
parameter_bounds[0][parameter_names.index('FCepsilonRI/K_1')] = 0
parameter_bounds[1][parameter_names.index('FCepsilonRI/K_1')] = 4000
parameter_bounds[0][parameter_names.index('FCepsilonRI/k_f2')] = 0.00023
parameter_bounds[1][parameter_names.index('FCepsilonRI/k_f2')] = 70
parameter_bounds[0][parameter_names.index('FCepsilonRI/K_2')] = 0
parameter_bounds[1][parameter_names.index('FCepsilonRI/K_2')] = 4000
parameter_bounds[0][parameter_names.index('FCepsilonRI/K_3')] = 0
parameter_bounds[1][parameter_names.index('FCepsilonRI/K_3')] = 4000
parameter_bounds[0][parameter_names.index('FCepsilonRI/k_f5')] = 0.0009
parameter_bounds[1][parameter_names.index('FCepsilonRI/k_f5')] = 101.16
parameter_bounds[0][parameter_names.index('FCepsilonRI/k_r5')] = 0.001
parameter_bounds[1][parameter_names.index('FCepsilonRI/k_r5')] = 18
parameter_bounds[0][parameter_names.index('FCepsilonRI/k_f6')] = 0.0009
parameter_bounds[1][parameter_names.index('FCepsilonRI/k_f6')] = 101.16
parameter_bounds[0][parameter_names.index('FCepsilonRI/k_r6')] = 0.001
parameter_bounds[1][parameter_names.index('FCepsilonRI/k_r6')] = 101.16
parameter_bounds[0][parameter_names.index('FCepsilonRI/K_7')] = 0
parameter_bounds[1][parameter_names.index('FCepsilonRI/K_7')] = 4000
parameter_bounds[0][parameter_names.index('FCepsilonRI/V_7')] = 0
parameter_bounds[1][parameter_names.index('FCepsilonRI/V_7')] = 4000
parameter_bounds[0][parameter_names.index('FCepsilonRI/Lyn')] = 0
parameter_bounds[1][parameter_names.index('FCepsilonRI/Lyn')] = 100
parameter_bounds = tuple(parameter_bounds)

# Run the simulation using given parameter values and return the
# values of the state variable for which we have experiment data
def model_function_lsq(params, expt_time, expt_pSyk, expt_pGRB2, return_type, debug=False):
    if debug:
        print('Parameters:')
        print(params)
    simulation.resetParameters()   
    for n, v in enumerate(params):
        constants[parameter_names[n]] = v
    try:
        simulation.run()
    except RuntimeError:
        print("Runtime error:")
        for n, v in enumerate(params):
            print('  {}: {}'.format(parameter_names[n], v))
        raise

    if return_type == 'optimisation':
        f1 = results.states()[expt_state_uri_pSyk].values()-expt_pSyk
        f2 = results.states()[expt_state_uri_pGrb2].values()-expt_pGRB2        
        f = np.concatenate((f1,f2))
        print('SSD:')    
        print(sum(f**2))
    elif return_type == 'visualisation':
        f1 = results.states()[expt_state_uri_pSyk].values()
        f2 = results.states()[expt_state_uri_pGrb2].values()        
        f = np.vstack((f1,f2))
    return f
    
#minimize(model_function_lsq, initial_params, args=(expt_time,expt_pSyk,expt_pGRB2), bounds=np.transpose(np.array(parameter_bounds)), method='SLSQP',options={'xtol': 1e-8, 'disp': True})
# least_squares(model_function_lsq, initial_params, args=(expt_time,expt_pSyk,expt_pGRB2),
#              bounds=parameter_bounds, xtol=1e-5,verbose=1)

opt =least_squares(model_function_lsq, initial_params, args=(expt_time,expt_pSyk,expt_pGRB2, 'optimisation'),
                               bounds=parameter_bounds,xtol=1e-4,verbose=1)
print('Optimal parameters:')
for n, v in enumerate(opt.x):
    print('  {}: {:g} ({:g})'.format(parameter_names[n], v, initial_params[n]))
    
f =model_function_lsq(opt.x, expt_time, expt_pSyk, expt_pGRB2,'visualisation', debug=True)
#print(f)

fig, ax = plt.subplots()
plt.plot(expt_time, expt_pSyk, '.', label='Experiment pSyk')
plt.plot(expt_time, f[0], '*', label='Model pSyk', color='red')
pSyk_error = f[0] - expt_pSyk
print(np.mean(pSyk_error))
print(np.std(pSyk_error))
print(max(expt_pSyk))
print(max(pSyk_error))
print(min(pSyk_error))
fig.canvas.draw()
plt.show()


fig, ax = plt.subplots()
plt.plot(expt_time, expt_pGRB2, '.', label='Experiment pGRB2')
plt.plot(expt_time, f[1], '*', label='Model pGRB2', color='red')
pGRB2_error = f[1] - expt_pGRB2
print(np.mean(pGRB2_error))
print(np.std(pGRB2_error))
print(max(expt_pGRB2))
print(max(pGRB2_error))
print(min(pGRB2_error))
fig.canvas.draw()
plt.show()















