import sympy as sp
import heapq
import itertools


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
    Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    I2 = sp.eye(2)
    return sp.kronecker_product(Y, I2)

def get_swap_target():
    return sp.Matrix([
        [1, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1]
    ])

# ==========================================
# 2. Heuristic and Distance Mathematics
# ==========================================

def calculate_frobenius_distance(mat1, mat2):
    # Calculate geometric distance between two matrices
    difference = mat1 - mat2
    
    # Extract numerical magnitude for complex/symbolic entries to sort in heapq
    sum_of_squares = sum(float(abs(val))**2 for val in difference)
    
    return sum_of_squares ** 0.5

def calculate_max_step_size(basis_gates, dimension=4):
    # Find the maximum distance a single gate can jump in the matrix space
    identity = sp.eye(dimension)
    max_distance = 0.0
    
    for gate in basis_gates.values():
        dist = calculate_frobenius_distance(gate, identity)
        if dist > max_distance:
            max_distance = dist
            
    return max_distance

# ==========================================
# 3. Flat Helper Functions
# ==========================================

def hash_matrix(matrix):
    return tuple(matrix)

def expand_astar_state(current, target, basis, queue, visited, step_limit, tie_breaker):
    current_matrix = current["matrix"]
    current_depth = current["depth"]
    new_depth = current_depth + 1
    
    for gate_name, gate_matrix in basis.items():
        new_matrix = gate_matrix * current_matrix
        state_hash = hash_matrix(new_matrix)
        
        # Pruning: skip if we've seen this exact state via a shorter or equal path
        is_visited = state_hash in visited
        is_redundant = is_visited and visited[state_hash] <= new_depth
        
        if is_redundant:
            continue
            
        visited[state_hash] = new_depth
        
        # Calculate heuristic (Estimated remaining gates)
        distance = calculate_frobenius_distance(new_matrix, target)
        h_score = distance / step_limit if step_limit > 0 else 0
        
        # Total priority score: f(n) = g(n) (gates used) + h(n) (estimated remaining)
        f_score = new_depth + h_score
        
        new_sequence = list(current["sequence"])
        new_sequence.append(gate_name)
        
        # Push to priority queue (tie_breaker ensures sympy matrices aren't directly compared)
        new_state = (f_score, next(tie_breaker), new_depth, new_matrix, new_sequence)
        heapq.heappush(queue, new_state)

# ==========================================
# 4. A* Search Algorithm
# ==========================================

def find_minimum_sequence_astar(target_matrix, basis_gates):
    dimensions = target_matrix.shape[0]
    identity = sp.eye(dimensions)
    
    # Pre-calculate the maximum step size for our admissible heuristic
    max_step = calculate_max_step_size(basis_gates, dimensions)
    
    # Priority Queue setup
    priority_queue = []
    tie_breaker = itertools.count() # Safely breaks ties for equal f_scores
    
    # Initial state
    initial_h = calculate_frobenius_distance(identity, target_matrix) / max_step
    initial_f = 0 + initial_h
    
    # State tuple: (f_score, tie_breaker, depth, matrix, sequence)
    heapq.heappush(priority_queue, (initial_f, next(tie_breaker), 0, identity, []))
    
    visited_depths = {hash_matrix(identity): 0}
    
    # Flat loop processing the priority queue
    while priority_queue:
        f_score, _, current_depth, current_matrix, current_sequence = heapq.heappop(priority_queue)
        
        # Because we use an admissible heuristic, the FIRST exact match popped 
        # from the priority queue is mathematically guaranteed to be the minimum depth.
        if current_matrix == target_matrix:
            return current_sequence
            
        # Repackage state for the expander function
        current_state = {
            "matrix": current_matrix,
            "depth": current_depth,
            "sequence": current_sequence
        }
        
        expand_astar_state(
            current_state, target_matrix, basis_gates, 
            priority_queue, visited_depths, max_step, tie_breaker
        )
        
    return None

# ==========================================
# 5. Main Execution
# ==========================================

print("Starting Heuristic A* quantum circuit synthesiser...\n")

target = get_swap_target()

basis = {
    "CNOT(0,1)": get_cnot_01(),
    "CNOT(1,0)": get_cnot_10(),
    "Y(0)": get_complex_y_gate()
}

optimal_circuit = find_minimum_sequence_astar(target, basis)

if optimal_circuit:
    print("Target Matrix: SWAP synthesised successfully!")
    print(f"Minimum gates required: {len(optimal_circuit)}")
    print("Optimal circuit path: " + " -> ".join(optimal_circuit))
else:
    print("Could not synthesise target matrix.")