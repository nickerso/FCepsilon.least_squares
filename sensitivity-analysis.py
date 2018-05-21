from collections import OrderedDict
import numpy as np
import math
from SALib.sample import saltelli
from SALib.analyze import sobol

import OpenCOR as oc

times = np.array([0, 300, 600, 900, 1200, 1500, 1800, 2100, 2400, 2700, 3000, 3300, 3600])
pGRB2 = np.array([0.0, 1.42, 2.97, 4.51, 5.47, 5.85, 6.19, 6.27, 6.32, 6.75, 6.49, 6.54, 6.62])


class Simulation(object):
    def __init__(self):
        self.simulation = oc.simulation()
        self.simulation.data().setStartingPoint(0)
        self.simulation.data().setEndingPoint(3600)
        self.simulation.data().setPointInterval(1)
        self.constants = self.simulation.data().constants()
        self.model_constants = OrderedDict({k: self.constants[k]
                                            for k in self.constants.keys()})
        
        # default the parameter bounds to something sensible, needs to be set directly
        bounds = []
        scale = 10.0
        for c in self.constants:
            v = self.constants[c];
            bounds.append([v/scale, v*scale])
        # define our sensitivity analysis problem
        self.problem = {
                   'num_vars': len(self.constants),
                   'names': self.constants.keys(),
                   'bounds': bounds
                   }
        self.samples = saltelli.sample(self.problem, 3)

    def run_once(self, c, v):
        self.simulation.resetParameters()
        self.constants[c] = v
        self.simulation.run()
        return (self.simulation.results().points().values(),
                self.simulation.results().states()['FCepsilonRI/pGrb2'].values())

    def run(self, c, scale=2.0):
        self.simulation.clearResults()
        v = self.model_constants[c]
        base = self.run_once(c, v)[1][times]
        divergence = 0.0
        for s in [1.0/scale, scale]:
            trial = self.run_once(c, s*v)[1][times]
            divergence += math.sqrt(np.sum((base - trial)**2))
        return divergence

    def evaluate_model(self, parameter_values):
        self.simulation.clearResults()
        for i, k in enumerate(self.model_constants.keys()):
            self.constants[k] = parameter_values[i]
        #print('Parameter set: ', parameter_values)
        self.simulation.run()
        return (self.simulation.results().states()['FCepsilonRI/pGrb2'].values()[times])
        
        
    def test_run(self):
        Y = np.zeros([self.samples.shape[0],len(times)])
        #print(Y)
        #exit
        #Y = np.array()
        for i, X in enumerate(self.samples):
            np.append(Y, self.evaluate_model(X), axis=1)
            print(Y)
            #Y[i] = self.evaluate_model(X)
            
#         variation = OrderedDict()
#         for c in self.model_constants.keys():
#             variation[c] = self.run(c)
#         return variation
        return Y

    def test(self, c, v):
        trial = self.run_once(c, s*v)[1][times]
        return math.sqrt(np.sum((pGRB2 - trial)**2))

s = Simulation()

v = s.test_run()

#print({ k: d for k, d in v.items() if d > 0.001 })
Si = sobol.analyze(s.problem, v, print_to_console=True)


'''
FCepsilonRI/k_f3 9.682260032327033
FCepsilonRI/K_3 1.3755810584148058
'''
