"""
QAOA Example code from qiskit textbook chapter 4.
"""

import networkx as nx
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit
from qiskit import Aer, execute
from qiskit.circuit import Parameter
from scipy.optimize import minimize
from qiskit.visualization import plot_histogram



'''Cost function'''
def maxcut_obj(x, G):
    obj = 0
    for i, j in G.edges():
        if x[i] != x[j]:
            obj -= 1
            
    return obj


'''Computes average cost of all the bitstrings outputted by the circuit, multiplicity accounted for'''
def compute_expectation(counts, G):
    avg = 0
    sum_count = 0
    for bitstring, count in counts.items():
        
        obj = maxcut_obj(bitstring, G)
        avg += obj * count
        sum_count += count
        
    return avg/sum_count



'''This function starts by hadamarding everything. Then it does the following p times: performs
Rzz on pairs of connected qubits, using the given angle as a parameter, then it performs Rx on
every qubit using, again, the appropriate parameter.'''
def create_qaoa_circ(G, theta):

    nqubits = len(G.nodes())
    p = len(theta)//2  # number of alternating unitaries
    qc = QuantumCircuit(nqubits)
    
    beta = theta[:p]
    gamma = theta[p:]
    
    # initial_state
    for i in range(0, nqubits):
        qc.h(i)
    
    for irep in range(0, p):
        
        # problem unitary
        for pair in list(G.edges()):
            '''The Rzz gates here correspond to terms in the hamiltonian that select two qubits'''
            qc.rzz(2 * gamma[irep], pair[0], pair[1])

        # mixer unitary
        for i in range(0, nqubits):
            qc.rx(2 * beta[irep], i)
            
    qc.measure_all()
        
    return qc


'''Execute circuit in simulator'''
def get_expectation(G, shots=512):
    
    backend = Aer.get_backend('qasm_simulator')
    backend.shots = shots
    
    def execute_circ(theta):
        qc = create_qaoa_circ(G, theta)
        counts = backend.run(qc, seed_simulator=10, nshots=512).result().get_counts()
        return compute_expectation(counts, G)
    
    return execute_circ


G = nx.Graph()
G.add_nodes_from([0, 1, 2])
G.add_edges_from([(0, 1), (1, 2)])

adjacency = nx.adjacency_matrix(G).todense()

nqubits = 3

expectation = get_expectation(G)
# scipy.minimize allows you to give it a whole function. the function we give it, 'expectation,'
# is the execute_circ function from get_expectation. execute_circ creates the quantum circuit,
# runs it on the backend with 512 shots, and returns the average cost function of all seen outputs.
# the second input to scipy.minimize seems to be an initial input to the function.
res = minimize(expectation, [1.0, 1.0], method='COBYLA')

# now we run the quantum circuit one last time on the trained parameters, and see the correct result.
backend = Aer.get_backend('aer_simulator')
backend.shots = 512

qc_res = create_qaoa_circ(G, res.x)

counts = backend.run(qc_res, seed_simulator=10).result().get_counts()

plot_histogram(counts)
plt.show()
