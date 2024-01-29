import numpy as np
from qiskit.circuit import QuantumCircuit, Parameter

"""Circuits for the kagome VQE"""

num_qubits = 16

def multi_cnot(qc,q1,q2,mult):
    """Apply CNOT gate between q1 and q2 mult times"""
    for _ in range(2*mult+1):
        qc.cx(q1,q2)
        
def multi_swap(qc,q1,q2,mult):
    """Apply SWAP gate between q1 and q2 mult times"""
    for _ in range(2*mult+1):
        qc.cx(q1,q2)
        qc.cx(q2,q1)
        qc.cx(q1,q2)

def u_heis(chi, psi, omega, cnot_multiply = 0):
    """Heisenberg gate"""
    qc = QuantumCircuit(2)
    for i in (0,1):
        qc.rz(np.pi/2, i)
        qc.sx(i)
        qc.rz(np.pi/2, i)
    multi_cnot(qc,0,1,cnot_multiply)
    qc.rz(chi,1)
    multi_cnot(qc,0,1,cnot_multiply)
    qc.sx(0)
    qc.sx(1)
    multi_cnot(qc,0,1,cnot_multiply)
    qc.rz(psi,1)
    multi_cnot(qc,0,1,cnot_multiply)
    for i in (0,1):
        qc.rz(np.pi/2, i)
        qc.sx(i)
    multi_cnot(qc,0,1,cnot_multiply)
    qc.rz(omega,1)
    multi_cnot(qc,0,1,cnot_multiply)
    return qc

def qubit_pair(p, rotate = 0):
    """Qubit pair corresponding to p after a rotation by rotate"""
    qubits = (1,2,3,5,8,11,14,13,12,10,7,4)
    i = (p[0]+rotate)%len(qubits)
    j = (p[1]+rotate)%len(qubits)
    return (qubits[i], qubits[j])

def u_h(i, theta, cnot_multiply = 0, rotate = 0):
    """Circuits U_H^i for i in (1,2,3,4)"""
    i = (i-1)%4 # i is now in (0,1,2,3)
    qubit_pairs = (((1,0),(11,10),(9,8),(7,6),(5,4),(3,2)),
                   ((0,11),(10,9),(8,7),(6,5),(4,3),(2,1)),
                   ((0,11),(8,7),(4,3)),
                   ((10,9),(6,5),(2,1)))[i]
    swap_pairs = ((),(),((1,0),(9,8),(5,4)),((11,10),(7,6),(3,2)))[i]
    qc = QuantumCircuit(num_qubits)
    for pair in swap_pairs:
        p = qubit_pair(pair, rotate)
        multi_swap(qc,p[0],p[1],cnot_multiply)
    parameter_index = 0
    for pair in qubit_pairs:
        p = qubit_pair(pair, rotate)
        par = (theta[parameter_index],theta[parameter_index+1],theta[parameter_index+2])
        parameter_index += 3
        qc.compose(u_heis(par[0],par[1],par[2],cnot_multiply),qubits=p,inplace=True)
    for pair in swap_pairs:
        p = qubit_pair(pair, rotate)
        multi_swap(qc,p[0],p[1],cnot_multiply)
    return qc

def u_dimer(cnot_multiply = 0, rotate = 0):
    """Dimer groundstate of H_1"""
    qubit_pairs = ((1,0),(11,10),(9,8),(7,6),(5,4),(3,2))  
    qc = QuantumCircuit(num_qubits)
    for pair in qubit_pairs:
        p = qubit_pair(pair, rotate)
        qc.rz(np.pi/2,p[1])
        qc.sx(p[1])
        qc.rz(np.pi/2,p[1])
        multi_cnot(qc,p[1],p[0],cnot_multiply)
        qc.rz(np.pi,p[0])
        qc.x(p[1])
    return qc

def u_hva_dimer(p, cnot_multiply = 0, rotate = 0):
    """HVA ansatz for dimer ground state"""
    theta = []
    digits_of_p = len(str(p))
    for layer in range(p):
        layer_params = []
        for edge in range(18):
            layer_params.append([Parameter('theta'+str(layer).rjust(digits_of_p,'0')
                                 +chr(97+edge)+str(par)) for par in range(3)])
        theta.append(layer_params)
    qc = u_dimer(cnot_multiply,rotate)
    circuit_index = (4,3,2,1)
    circuit_edges = ((13,15,17),(12,14,16),(1,3,5,7,9,11),(0,2,4,6,8,10))
    for layer in range(p):
        for ci, edges in zip(circuit_index,circuit_edges):
            parameters = []
            for edge in edges:
                for par in theta[layer][edge]:
                    parameters.append(par)
            qc.compose(u_h(ci,parameters,cnot_multiply,rotate),inplace=True)
    return qc

def error_mitigated_ansaetze(num_layers = 1, rots = (0,3,6,9), cnot_mults = (0,1,2)):
    """Introduction of error mitigation into the HVA dimer ansatz.

    Args:
        num_layers: number of layers in the ansatz
        rots: rotations of the ansatz
        cnot_mults: number of times the CNOT gates are repeated for zero error extrapolation
    """
    ansatz = []
    for rot in rots:
        for cnot_mult in cnot_mults:
            ansatz.append(u_hva_dimer(num_layers,cnot_mult,rot))
    return ansatz

def initial_parameters_dimer(num_layers):
    """Initial parameters for the dimer ansatz"""
    initial_hamiltonian_edges = (0,2,4,6,8,10)
    theta = []
    for layer in range(num_layers):
        sqp = np.sqrt(num_layers)
        pa = 1/num_layers
        pb = (np.floor((1+layer)/sqp) - .5) / sqp / num_layers
        for edge in range(18):
            if edge in initial_hamiltonian_edges:
                p = pb
            else:
                p = pa
            for _ in range(3):
                theta.append(p)
    return theta    
