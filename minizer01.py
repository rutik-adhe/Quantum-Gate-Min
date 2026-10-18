import sympy as sp
from collections import deque


# 4x4 Identity matrix (Starting state of our circuit)
I4 = sp.eye(4)

# CNOT with control qubit 0 and target qubit 1
CNOT_01 = sp.Matrix([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0]
])

# CNOT with control qubit 1 and target qubit 0
CNOT_10 = sp.Matrix([
    [1, 0, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0],
    [0, 1, 0, 0]
])

# Pauli-Y gate (Contains complex number 'i' which is sp.I in sympy)
# We tensor it with a 2x2 Identity to make it a valid 4x4 2-qubit gate
Y_gate = sp.Matrix([
    [0, -sp.I],
    [sp.I, 0]
])
I2 = sp.eye(2)
Y_on_qubit_0 = sp.kronecker_product(Y_gate, I2)

# SWAP gate (Target matrix for the test case)
SWAP_TARGET = sp.Matrix([
    [1, 0, 0, 0],
    [0, 0, 1, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1]
])

# Dictionary of available gates to use in our search programme
basis_gates = {
    "CNOT(0,1)": CNOT_01,
    "CNOT(1,0)": CNOT_10,
    "Y(0)": Y_on_qubit_0
}

# ==========================================
# 2. Algorithm Setup: Breadth-First Search
# ==========================================

# Queue stores tuples of: (Current Matrix, List of Applied Gates)
queue = deque()
queue.append((I4, []))

# We use a set to track visited matrices to avoid redundant calculations.
# Sympy matrices are mutable and not directly hashable, so we convert 
# them to tuples of their elements before storing them in the set.
visited = set()
initial_hash = tuple(I4)
visited.add(initial_hash)

target_found = False
minimum_circuit = []

print("Starting the BFS programme to minimise gate count...")


while len(queue) > 0:
    current_matrix, current_circuit = queue.popleft()

    # Check if we successfully synthesised the target matrix
    if current_matrix == SWAP_TARGET:
        target_found = True
        minimum_circuit = current_circuit
        break

    # Apply each available basis gate to the current state
    for gate_name, gate_matrix in basis_gates.items():
        
        # Quantum gates are applied backwards mathematically: New = Gate * Current
        new_matrix = gate_matrix * current_matrix

        # Create a hashable representation to check our visited set
        matrix_hash = tuple(new_matrix)

        if matrix_hash not in visited:
            # Mark as visited
            visited.add(matrix_hash)
            
            # Record the sequence of gates
            new_circuit = list(current_circuit)
            new_circuit.append(gate_name)

            # Add to the back of the queue for later exploration
            queue.append((new_matrix, new_circuit))


if target_found:
    print("\nTarget matrix synthesised successfully!")
    print("Minimum number of gates required:", len(minimum_circuit))
    print("Optimal circuit sequence:", " -> ".join(minimum_circuit))
else:
    print("\nTarget matrix could not be synthesised with the given basis gates.")