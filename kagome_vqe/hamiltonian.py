from qiskit.quantum_info import SparsePauliOp
from qiskit.opflow.primitive_ops import PauliSumOp

"""Hamiltonian for the Kagome lattice."""

num_qubits = 16

def kagome_hamiltonian() -> PauliSumOp:
    """Returns the Hamiltonian for the Kagome lattice."""
    # every edge in the Kagome lattice contributes three Pauli
    # strings with X/Y/Z on the corresponding vertices:
    kagome_edges = [(2,1),(1,4),(4,7),(7,10),(10,12),(12,13),
                    (13,14),(14,11),(11,8),(8,5),(5,3),(3,2),
                    (2,4),(4,10),(10,13),(13,11),(11,5),(5,2)]
        
    # Construct the hamiltonians for the estimator
    pauli_strings = {'X': [], 'Y': [], 'Z': []}
    for e in kagome_edges:
        for c in 'XYZ':
            s = ''
            for i in range(num_qubits):
                if i in e:
                    s = c + s
                else:
                    s = 'I' + s
            pauli_strings[c].append(s)
    return PauliSumOp(SparsePauliOp(pauli_strings['X']
                +pauli_strings['Y']+pauli_strings['Z']))
