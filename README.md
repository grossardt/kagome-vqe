# kagome-vqe
 Ground state estimation for a Kagome lattice.
 This code was originally written motivated by the IBM Open Science Prize 2023, but it is not a 
 submission to the contest. We implement the VQE algorithm to estimate the ground state energy of
 the Kagome lattice using Rotosolve for optimization and an HVA ansatz.

 One idea was to base the HVA ansatz on trimer or double-dimer elements of the lattice rather than
 the dimer building blocks used in other papers. The main goal we were following was to explore 
 orthogonal subspaces of the Hilbert space by exploiting the orthogonality of the ansatz states.
 This idea, however, turned out to be disadvantageous.

 The code also includes some first attempts to implement certain error mitigation techniques.

 A more detailed description of the code can be found in the LaTeX and PDF files in the `notes` 
 folder.

## Reuse
The code is an unfinished work in progress and should be treated as such, but you are free to use
it under the MIT license.
