import sympy as sp
from collections import deque

def get_cnot_01():
    return sp.Matrix([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0]
    ])

def get_cnot_10():
    return sp.Matrix([
        [1, 0, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
        [0, 1, 0, 0]
    ])

def get_complex_y_gate():
    # Demonstrating arbitrary complex numbers handling (sp.I is imaginary 'i')
    Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    I2 = sp.eye(2)
    # Tensor product to make it a 4x4 2-qubit gate
    return sp.kronecker_product(Y, I2)

def get_swap_target():
    return sp.Matrix([
        [1, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1]
    ])

def hash_matrix(matrix):
    # Sympy matrices must be converted to tuples to be used as dictionary keys
    return tuple(matrix)

def expand_state(current_matrix, current_sequence, basis_gates, queue, visited_depths):
    # Applies all basis gates to the current state and queues valid new states
    current_depth = len(current_sequence)
    new_depth = current_depth + 1
    
    for gate_name, gate_matrix in basis_gates.items():
        # Apply gate (New = Gate * Current)
        new_matrix = gate_matrix * current_matrix
        state_hash = hash_matrix(new_matrix)
        
        # Check if we have seen this matrix via a shorter or equal path
        is_visited = state_hash in visited_depths
        is_redundant = is_visited and visited_depths[state_hash] < new_depth
        
        if is_redundant:
            continue
            
        # Record the new shortest depth for this specific matrix
        visited_depths[state_hash] = new_depth
        
        # Create new sequence and enqueue
        new_sequence = list(current_sequence)
        new_sequence.append(gate_name)
        queue.append((new_matrix, new_sequence))


def find_all_minimum_sequences(target_matrix, basis_gates, max_search_depth=5):
    # Setup initial state (Identity matrix)
    dimensions = target_matrix.shape[0]
    identity = sp.eye(dimensions)
    
    queue = deque([(identity, [])])
    visited_depths = {hash_matrix(identity): 0}
    
    solutions = []
    minimum_depth_found = max_search_depth + 1
    
    # Flat loop processing the queue
    while queue:
        current_matrix, current_sequence = queue.popleft()
        current_depth = len(current_sequence)
        
        # Stop searching if we exceed the depth of our found optimal solutions
        if current_depth > minimum_depth_found:
            break
            
        # Check if we found a match
        if current_matrix == target_matrix:
            minimum_depth_found = current_depth
            solutions.append(current_sequence)
            continue
            
        # Prevent expansion beyond the optimal depth once a solution is found
        if current_depth == minimum_depth_found:
            continue
            
        # Expand current state by applying all available gates
        expand_state(current_matrix, current_sequence, basis_gates, queue, visited_depths)
        
    return solutions

print("Starting generalised quantum circuit synthesiser...\n")

target = get_swap_target()

basis = {
    "CNOT_01": get_cnot_01(),
    "CNOT_10": get_cnot_10(),
    "Y_0": get_complex_y_gate()
}

optimal_circuits = find_all_minimum_sequences(target, basis, max_search_depth=5)

print("Target Matrix: SWAP")
print(f"Total optimal sequences found: {len(optimal_circuits)}")

# Print all found minimum sequences without nested loops
def print_circuit(circuit, index):
    formatted_path = " -> ".join(circuit)
    print(f"Solution {index}: {formatted_path}")

for i, circuit in enumerate(optimal_circuits, 1):
    print_circuit(circuit, i)