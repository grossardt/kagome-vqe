import os
from argparse import ArgumentParser
from datetime import datetime

from .estimator import RetryEstimator, LocalEstimator
from qiskit.primitives import Estimator as FastEstimator
from .hamiltonian import kagome_hamiltonian
from .circuits import error_mitigated_ansaetze, initial_parameters_dimer
from .optimization import optimize
from .cost import error_mitigated_cost


if __name__ == "__main__":
    parser = ArgumentParser(description='KagomeVQE with error mitigation.')
    parser.add_argument('--out', metavar='DIR', type=str, default='./out', help='output directory')
    parser.add_argument('--layers', metavar='N', type=int, default=1, help='number of layers')
    parser.add_argument('--noise', action='store_true', help='perform VQE with noise (default: no noise)')
    parser.add_argument('--slow', action='store_true', help='even without noise use Aer Estimator (no effect with noise)')
    parser.add_argument('--nocnot', action='store_true', help='skip error mitigation by cnot multiplication')
    parser.add_argument('--norot', action='store_true', help='skip error mitigation by rotation')
    parser.add_argument('--maxsteps', metavar='M', type=int, default=1000, help='max. number of optimization steps')
    parser.add_argument('--delta', metavar='X', type=float, default=0.001, help='stop if the energy changes by less than X for a number of steps equal to the number of parameters')
    parser.add_argument('--continue', metavar='FILE', type=str, default='', help='continue from previous calculation, all other arguments will be ignored')
    args = vars(parser.parse_args())
    
    # Defaults:
    step = 0
    cost = None
    slow = False
    no_cnot = False
    no_rot = False
    
    if args['continue']:
        outfile = args['continue']
        with open(outfile, 'r') as f:
            line = f.readline()
            if line[13:20] == 'without':
                noise = False
            else:
                noise = True
            num_layers = int(line.split('with ')[-1].split(' layer')[0])
            if noise:
                additional_parameters = line.split('(')[-1].split(')')[0].split('cnot ')[-1].split(', rot ')
                if additional_parameters[0] == 'no':
                    no_cnot = True
                if additional_parameters[1] == 'no':
                    no_rot = True
            else:
                if line.split('(')[-1].split(')')[0] == 'slow simulation':
                    slow == True
            line = f.readline()
            stopping_criteria = line.split('max_steps = ')[-1].split(', stopping_delta = ')
            max_steps = int(stopping_criteria[0])
            stopping_delta = float(stopping_criteria[1])
            line = f.readline()
            while line:
                if line[:4] == 'STEP':
                    step = int(line.split(',')[0][5:])
                elif line[:7] == '  COST:':
                    cost = float(line[7:])
                elif line[:9] == '  PARAMS:':
                    theta = list(map(float,line[11:-2].split(',')))
                line = f.readline()
            step += 1
        
    else:
        noise = args['noise']
        slow = args['slow']
        num_layers = args['layers']
        max_steps = args['maxsteps']
        stopping_delta = args['delta']
        no_rot = args['norot']
        rot_text = ('yes','no')[args['norot']]
        no_cnot = args['nocnot']
        cnot_text = ('yes','no')[args['nocnot']]
        outfile = args['out'] + '/' + datetime.now().strftime('%Y%m%d-%H%M%S-kagome-errmit.txt')
        if not os.path.exists(args['out']):
            os.makedirs(args['out'])
        with open(outfile, 'a') as f:
            noise_text = ('without','with')[noise] + ' noise'
            if noise:
                noise_text += ' (error mitigation: cnot {}, rot {})'.format(cnot_text,rot_text)
            elif slow:
                noise_text += ' (slow simulation)'
            f.write('STARTING VQE {} with {} layer ansatz\n'.format(noise_text,num_layers))
            f.write('STOPPING CRITERION: max_steps = {}, stopping_delta = {}\n\n'.format(max_steps,stopping_delta))
            f.close()
    
    rots = ((0,3,6,9),(0,))[no_rot]
    cnot_mults = ((0,1,2),(0,))[no_cnot]
    
    print('Welcome to KagomeVQE. Preparing...')
    
    ansaetze = error_mitigated_ansaetze(num_layers, rots, cnot_mults)
    hamiltonian = kagome_hamiltonian()
    if noise or slow:
        estimator = LocalEstimator(noise=noise)
    else:
        estimator = FastEstimator()
    if not step:
        theta = initial_parameters_dimer(num_layers)
    
    def cost_function(theta1, theta2 = None):
        return error_mitigated_cost(estimator, hamiltonian, ansaetze, rots, cnot_mults, theta1, theta2)
    
    print('Starting optimization...')
    print('  layers: {}, noise: {}, slow: {}, cnot mit.: {}, rot mit.: {}'.format(num_layers,('N','Y')[noise],('N','Y')[slow],('Y','N')[no_cnot],('Y','N')[no_rot]))
    print('  maximum steps: {}, stopping delta: {}, initial step: {}, initial cost: {}'.format(max_steps,stopping_delta,step,cost))
    print('  Output will be written to ' + outfile)
        
    optimize(cost_function, theta, max_steps, stopping_delta, outfile, step, cost)
