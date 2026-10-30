import os
import gc
import sympy as sp

# Threading Configuration
os.environ["OMP_NUM_THREADS"] = str(os.cpu_count())
os.environ["OPENBLAS_NUM_THREADS"] = str(os.cpu_count())
os.environ["MKL_NUM_THREADS"] = str(os.cpu_count())

import numpy as np

# Numpy / Cupy Setup
try:
    import cupy as xp
    ACTIVE_BACKEND = "GPU Acceleration (CuPy with VRAM Batching)"
    def to_cpu(tensor): return tensor.get()
    def to_gpu(array): return xp.array(array)
except ImportError:
    xp = np
    ACTIVE_BACKEND = "CPU Multithreading (NumPy OpenBLAS with RAM Batching)"
    def to_cpu(tensor): return tensor
    def to_gpu(array): return array


def synthesize_circuit_batched(target_matrix_sympy, basis_gates_list, max_depth, target_fidelity):
    print(f"Initializing synthesis engine using {ACTIVE_BACKEND}...")
    
    # Sympy to Numpy
    gate_names = [name for name, _ in basis_gates_list]
    basis_np = np.array([np.array(gate.evalf(), dtype=np.complex128) for _, gate in basis_gates_list])
    target_np = np.array(target_matrix_sympy.evalf(), dtype=np.complex128)

    num_gates = basis_np.shape[0]
    dimension = target_matrix_sympy.shape[0]

    gates_gpu = to_gpu(np.expand_dims(basis_np, axis=1))
    target_gpu = to_gpu(target_np)

    # No of Matrices in each Batch
    VRAM_BATCH_SIZE = 500_000

    # GPU Memory Limit
    RAM_LIMIT_BYTES = 4_000_000_000 # 4 GB threshold before using SSD Swap

    current_depth = 1
    start_idx = 0
    layer_size = 1
    
    current_layer_cpu = np.eye(dimension, dtype=np.complex128).reshape(1, dimension, dimension)
    next_layer_cpu = np.empty((num_gates, dimension, dimension), dtype=np.complex128)

    best_fidelity = -1.0
    best_index = -1
    best_depth = -1
    temp_files = []

    
    while current_depth <= max_depth:
        end_idx = min(start_idx + VRAM_BATCH_SIZE, layer_size)
        chunk_size = end_idx - start_idx
        
        chunk_gpu = to_gpu(current_layer_cpu[start_idx:end_idx])

        # Applying Gates to each Node of Batch
        states_expanded = xp.broadcast_to(chunk_gpu, (num_gates, chunk_size, dimension, dimension))
        next_chunk_gpu = xp.matmul(gates_gpu, states_expanded).reshape(num_gates * chunk_size, dimension, dimension)
        
        # Batched Fidelity Calculation
        target_broadcast = xp.broadcast_to(target_gpu, (num_gates * chunk_size, dimension, dimension))
        trace_values = xp.sum(xp.conj(next_chunk_gpu) * target_broadcast, axis=(1, 2))
        fidelities = xp.abs(trace_values) / dimension
        
        # Finding Best Fidelity
        chunk_best_idx = int(xp.argmax(fidelities))
        chunk_best_fid = float(fidelities[chunk_best_idx])
        global_best_idx_for_this_chunk = (start_idx * num_gates) + chunk_best_idx
        
        if chunk_best_fid > best_fidelity:
            best_fidelity = chunk_best_fid
            best_index = global_best_idx_for_this_chunk
            best_depth = current_depth
            
        # Saving results back to CPU RAM / SSD Swap
        insert_start = start_idx * num_gates
        insert_end = end_idx * num_gates
        next_layer_cpu[insert_start:insert_end] = to_cpu(next_chunk_gpu)
        
        # Advance to Next Batch
        start_idx += VRAM_BATCH_SIZE
        
        # Move to Next Depth if current layer is completed
        if start_idx >= layer_size:
            print(f"Depth {current_depth:02d} Complete | States Evaluated: {layer_size * num_gates:<12} | Best Fidelity: {best_fidelity:.6f}")
            
            if best_fidelity >= target_fidelity or current_depth == max_depth:
                break
                
            current_depth += 1
            start_idx = 0
            layer_size = next_layer_cpu.shape[0]
            
            # Free previous layer memory completely
            if hasattr(current_layer_cpu, 'filename'):
                temp_files.append(current_layer_cpu.filename)
            del current_layer_cpu
            gc.collect()
            
            current_layer_cpu = next_layer_cpu
            
            # Dynamically allocate the next layer to RAM or Memmap (Disk)
            bytes_req = layer_size * num_gates * dimension * dimension * 16
            if bytes_req > RAM_LIMIT_BYTES:
                filename = f'quantum_layer_swap_{current_depth}.dat'
                print(f"--> Memory Threshold Exceeded ({bytes_req / 1e9:.2f} GB). Allocating next layer to SSD swap file...")
                next_layer_cpu = np.memmap(filename, dtype=np.complex128, mode='w+', shape=(layer_size * num_gates, dimension, dimension))
            else:
                next_layer_cpu = np.empty((layer_size * num_gates, dimension, dimension), dtype=np.complex128)

    # Delete SSD Swap Files
    del current_layer_cpu
    del next_layer_cpu
    gc.collect()
    for file in temp_files:
        try:
            os.remove(file)
        except Exception:
            pass

    # Creating Path Sequence
    optimal_sequence = []
    residual_index = best_index
    divisor = num_gates ** (best_depth - 1)

    while divisor >= 1:
        gate_idx = residual_index // divisor
        residual_index = residual_index % divisor
        optimal_sequence.insert(0, gate_names[gate_idx])
        divisor = divisor // num_gates

    return optimal_sequence, best_fidelity, best_depth

# Main
if __name__ == "__main__":
    
    I2 = sp.eye(2)
    H = sp.Matrix([[1, 1], [1, -1]]) / sp.sqrt(2)
    S = sp.Matrix([[1, 0], [0, sp.I]])
    T = sp.Matrix([[1, 0], [0, sp.exp(sp.I * sp.pi / 4)]])

    def get_cnot_01():
        return sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])

    def get_cnot_10():
        return sp.Matrix([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]])

    target = sp.Matrix([
        [1, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1]
    ])

    basis = [
        ("H(0)", sp.kronecker_product(H, I2)),
        ("H(1)", sp.kronecker_product(I2, H)),
        ("S(0)", sp.kronecker_product(S, I2)),
        ("S(1)", sp.kronecker_product(I2, S)),
        ("T(0)", sp.kronecker_product(T, I2)),
        ("T(1)", sp.kronecker_product(I2, T)),
        ("CNOT(0,1)", get_cnot_01()),
        ("CNOT(1,0)", get_cnot_10())
    ]
    
    sequence, final_fidelity, depth = synthesize_circuit_batched(
        target_matrix_sympy=target,
        basis_gates_list=basis,
        max_depth=10, 
        target_fidelity=0.999
    )

    print("-" * 50)
    print(f"Synthesis Complete.")
    print(f"Achieved Fidelity:     {final_fidelity:.6f}")
    print(f"Required Gates:        {depth}")
    print("Optimal Gate Sequence: " + (" -> ".join(sequence) if sequence else "None"))