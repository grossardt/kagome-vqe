import numpy as np
import itertools as it
import networkx as nx
from qiskit.circuit import QuantumCircuit

"""Pauli string manipulation"""

def pauli_string_product(p, q):
    """
    Returns the matrix product of two Pauli strings of equal length
    input: p, q respresented by lists of integers 0...3
    output: tuple (pauli string, relative phase)
    """
    # Check input
    input_exception = Exception("Input must be lists of integers 0...3")
    if isinstance(p, list) and isinstance(q, list) and len(p) == len(q):
        for i in p:
            if i not in [0,1,2,3]:
                raise input_exception
                return
        for i in q:
            if i not in [0,1,2,3]:
                raise input_exception
                return
    else:
        raise input_exception
        return
    # List of length 1: product of two single qubit Pauli matrices
    if len(p) == 1:
        multiplication_table = {(0,0): 0, (0,1): 1, (0,2): 2, (0,3): 3,
                                (1,0): 1, (1,1): 0, (1,2): 3, (1,3): 2,
                                (2,0): 2, (2,1): 3, (2,2): 0, (2,3): 1,
                                (3,0): 3, (3,1): 2, (3,2): 1, (3,3): 0}
        phase_table = {(0,0): 0, (0,1): 0, (0,2): 0, (0,3): 0,
                       (1,0): 0, (1,1): 0, (1,2): 1, (1,3): 3,
                       (2,0): 0, (2,1): 3, (2,2): 0, (2,3): 1,
                       (3,0): 0, (3,1): 1, (3,2): 3, (3,3): 0}
        return ([multiplication_table[(p[0],q[0])]], phase_table[(p[0],q[0])])
    # Recurse for length > 1:
    (res1, phase1) = pauli_string_product(p[:-1], q[:-1])
    (res2, phase2) = pauli_string_product(p[-1:], q[-1:])
    return (res1+res2, (phase1+phase2) % 4)
    
class Pauli:
    """
    Class for Pauli Operators represented by string of the single qubit Paulis [I, X, Y, Z]
    """
    
    def __init__(self, value, phase = 0):
        """
        value can be either a list of integers 0...3 or a Pauli string
        """
        self.__matrix = None
        if isinstance(phase, int):
            self.phase = phase % 4
        else:
            raise Exception('Phase must be integer')
            return
        if isinstance(value, str):
            value_dict = {'I': 0, 'X': 1, 'Y': 2, 'Z': 3}
            self.value = []
            for c in value:
                if c not in 'IXYZ':
                    raise Exception('Not a valid Pauli string: ' + value)
                    return
                else:
                    self.value.append(value_dict[c])
        elif isinstance(value, list):
            # check list
            for i in value:
                if i not in [0,1,2,3]:
                    raise Exception('Not a valid list of integers 0...3')
                    return
            self.value = value
        else:
            raise Exception('Unknown input type')
    
    def __len__(self):
        """
        Length is number of qubits
        """
        return len(self.value)
    
    def __mul__(self, second):
        """
        Matrix multiplication of two Pauli strings returns other Pauli string
        """
        if not isinstance(second, Pauli):
            raise Exception('Pauli string can only be multiplied with another Pauli string')
            return
        if len(second.value) != len(self.value):
            raise Exception('Multiplication requires two Pauli strings of equal length')
            return
        (prod, phase) = pauli_string_product(self.value, second.value)
        return Pauli(prod, phase + self.phase + second.phase)
        
    def string(self):
        """
        Pauli string written with operators I, X, Y, Z, including phase prefactor
        """
        phase_dict = {0: '', 1: 'i*', 2: '-', 3: '-i*'}
        char_dict = {0: 'I', 1: 'X', 2: 'Y', 3: 'Z'}
        string = phase_dict[self.phase]
        for i in self.value:
            string += char_dict[i]
        return string
    
    def phase_string(self):
        """
        Pauli string written with operators I, X, Y, Z, including phase prefactor
        """
        phase_dict = {0: 1, 1: 1j, 2: -1, 3: -1j}
        char_dict = {0: 'I', 1: 'X', 2: 'Y', 3: 'Z'}
        phase = phase_dict[self.phase]
        string = ''
        for i in self.value:
            string += char_dict[i]
        return (phase,string)
    
    def matrix(self):
        """
        Matrix representation of Pauli string
        """
        if self.__matrix is None: # matrix form has not been calculated yet
            pauli_matrices = [np.matrix([[1,0],[0,1]]),
                              np.matrix([[0,1],[1,0]]),
                              np.matrix([[0,-1j],[1j,0]]),
                              np.matrix([[1,0],[0,-1]])
                             ]
            self.__matrix = np.array(1, dtype=np.complex128)
            for i in self.value:
                self.__matrix = np.kron(self.__matrix, pauli_matrices[i])
            self.__matrix *= 1j**self.phase
        return self.__matrix
    
    def circuit(self, parameter):
        """
        Implement Pauli operator as n+1 qubit circuit in star + ancilla implementation
        Returns qiskit.QuantumCircuit()
        """
        qc = QuantumCircuit(len(self.value) + 1)
        for qubit, op in enumerate(self.value):
            if op == 2:
                qc.sdg(qubit)
            if op in [1,2]:
                qc.h(qubit)
        for qubit, op in enumerate(self.value):
            if op:
                qc.cx(qubit,len(self.value))
        qc.rz(parameter, len(self.value))
        for qubit, op in reversed(list(enumerate(self.value))):
            if op:
                qc.cx(qubit,len(self.value))        
        for qubit, op in enumerate(self.value):
            if op in [1,2]:
                qc.h(qubit)
            if op == 2:
                qc.s(qubit)
        return qc
    
    def depth(self):
        """
        Circuit depth, only counting CNOT gates. Two for every X, Y, Z
        """
        res = 0
        for i in self.value:
            if i:
                res += 2
        return res
    
    def commutes_with(self, q):
        """
        Returns True if this Pauli commutes with the argument q, False otherwise
        Input: Pauli() or string
        """
        if isinstance(q, str):
            q = Pauli(q)
        if isinstance(q, Pauli):
            c = self * q
            d = q * self
            if (c.value == d.value):
                if (c.phase - d.phase) % 4:
                    return False
                else:
                    return True
            else:
                raise Exception('Product and inverse product have different Pauli strings')
        else:
            raise Exception('Input must be of Pauli() type') 
        
def gate_cancellations(p, q):
    """
    Returns the number of cancelled CNOTs between two Pauli string operators p, q
    """
    if isinstance(p, str):
        p = Pauli(p)
    if isinstance(q, str):
        q = Pauli(q)
    if isinstance(p, Pauli) and isinstance(q, Pauli) and len(p) == len(q):
        res = 0
        for i, s in enumerate(p.value):
            if s and s == q.value[i]:
                res += 2
        return res
    else:
        raise Exception('Arguments must be of Pauli() type or valid Pauli strings of equal length')
        return
    
def cnot_distance(p, q):
    """
    Returns the number of CNOTs between two Pauli string operators p, q after gate cancellation
    """
    if isinstance(p, str):
        p = Pauli(p)
    if isinstance(q, str):
        q = Pauli(q)
    if isinstance(p, Pauli) and isinstance(q, Pauli) and len(p) == len(q):
        res = 0
        for i, s in enumerate(p.value):
            if s:
                if s != q.value[i]:
                    res += 2
            else:
                if q.value[i]:
                    res += 1
        return res
    else:
        raise Exception('Arguments must be of Pauli() type or valid Pauli strings of equal length')

def total_depth(ops):
    """
    Returns the total number of CNOTs in a given product of Pauli string operators
    Input: list of Pauli() or Pauli strings of equal length
    """
    if not (isinstance(ops, list) and len(ops)):
        raise Exception('Argument must be list of Pauli strings')
    if isinstance(ops[0], str):
        ops[0] = Pauli(ops[0])
    n = len(ops[0])
    res = ops[0].depth()
    for i in range(1,len(ops)):
        if isinstance(ops[i], str):
            ops[i] = Pauli(ops[i])
        res += ops[i].depth()
        res -= gate_cancellations(ops[i-1], ops[i])
    return res
    

def matrix_to_pauli(m):
    """
    Reversely finds the Pauli string for a given matrix via brute force search.
    Returns False if no match is found.
    """
    n = 0
    while 2**n < len(m):
        n += 1
    if 2**n == len(m):
        # loop through all n qubit Pauli strings
        for i in it.product([0,1,2,3],repeat=n):
            p = Pauli(list(i))
            if np.all( p.matrix() == m ):
                return p.string()
        return False
    else:
        raise Exception('Input matrix of invalid dimension')
            
def com(p, q):
    """
    Returns HALF the commutator between two Pauli strings
    which is either 0 or +/- 2i times another Pauli string
    Input can be Pauli() objects or strings
    """
    if isinstance(p, str):
        p = Pauli(p)
    if isinstance(q, str):
        q = Pauli(q)
    if isinstance(p, Pauli) and isinstance(q, Pauli) and len(p) == len(q):
        c = p * q
        d = q * p
        if (c.value == d.value):
            phase_diff = (c.phase - d.phase) % 4
            if phase_diff == 0:
                return 0
            elif phase_diff == 2:
                phase_dict = {0: '+2 ', 1: '+2i ', 2: '-2 ', 3: '-2i '}
                return c
            else:
                raise Exception('Phase difference is neither 0 nor pi')
        else:
            raise Exception('Product and inverse product have different Pauli strings')
    else:
        raise Exception('Arguments must be of Pauli() type or valid Pauli strings of equal length')

def commutativity_graph(nodes):
    """
    Graph of Pauli strings with edges between any two that DO commute
    """
    g = nx.Graph()
    g.add_nodes_from(nodes)
    for i in range(len(nodes)):
        for j in range(i+1,len(nodes)):
            if Pauli(nodes[i]).commutes_with(Pauli(nodes[j])):
                g.add_edge(nodes[i], nodes[j])
    return g

def noncommutativity_graph(nodes):
    """
    Graph of Pauli strings with edges between any two that DO NOT commute
    """
    g = nx.Graph()
    g.add_nodes_from(nodes)
    for i in range(len(nodes)):
        for j in range(i+1,len(nodes)):
            if not Pauli(nodes[i]).commutes_with(Pauli(nodes[j])):
                g.add_edge(nodes[i], nodes[j])
    return g

def cnot_graph(nodes):
    """
    Graph of Pauli strings with edges between all of them with weights given by CNOT distances
    """
    g = nx.Graph()
    g.add_nodes_from(nodes)
    for i in range(len(nodes)):
        for j in range(i+1,len(nodes)):
            g.add_edge(nodes[i], nodes[j], weight=cnot_distance(nodes[i],nodes[j]))
    return g
