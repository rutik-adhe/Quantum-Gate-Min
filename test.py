import sympy as sp
import numpy as np
from sympy.physics.quantum import TensorProduct
from sympy.physics.units import h, electron_rest_mass
from sympy.physics.quantum.constants import hbar
from circuit05 import synthesize_circuit_batched
ONE = sp.Integer(1)

HBAR = ONE
MASS = sp.Rational(1, 2)

def getExponent(num, base) :
    count = 0
    isValid = True

    while num > 0 and isValid :
        reminder = num % base

        if reminder == 0 :
            num = num / base
            count = count + 1
        else :
            if num == 1 :
                num = 0
            else :
                isValid = False

    return count, isValid

SPATIAL_STEP_SIZE = ONE
TIME_STEP_SIZE = 0.01

INTERVAL_LOW = sp.Rational(-5,2)
INTERVAL_HIGH = sp.Rational(5,2)

IntervalLength = INTERVAL_HIGH - INTERVAL_LOW
LastGridPoint : sp.Integer = ((IntervalLength / SPATIAL_STEP_SIZE) - 1)

N, isPowerOf2 = getExponent(LastGridPoint, 2)

if(type(LastGridPoint) is not sp.Integer or (not isPowerOf2)) :
    print("Enter Spatial Step Size such that it create 2^N grid points.")
    exit(0)
else :
    print(f"Total Spatial Grid Points = {LastGridPoint}")

# Shifted Row Matrix Filling
# Take Size of Row = LastGridPoint + 1
# Rotate right by 1 at each iteration append to matrix

def rotateListRight(arr, len) -> None :
    temp = arr[len - 1]

    i = len - 1
    while(i > 0) :
        arr[i] = arr[i-1]
        i = i -1

    arr[0] = temp

def createKE_operator() :
    RowLength = LastGridPoint + 1 

    InitialRow = [0] * (int(RowLength))
    InitialRow[0] = sp.Integer(2)
    InitialRow[1] = -1*ONE
    InitialRow[LastGridPoint] = -1*ONE

    KineticEnergyOperator = []

    i = 1
    while (i < RowLength) :
        KineticEnergyOperator.append(InitialRow.copy()[:-1])
        rotateListRight(InitialRow, RowLength)
        i = i + 1

    M = sp.Matrix(KineticEnergyOperator)
    multiplicationFactor = (HBAR**2) / (2 * MASS)
    return multiplicationFactor * M

# PE Operator


def createPE_operator(potential_func) -> sp.Matrix :
    PotentialEneryOperator = []

    i = 1
    while(i <= LastGridPoint) :

        x_cord = INTERVAL_LOW + i*SPATIAL_STEP_SIZE
        diag_ele = potential_func(x_cord)

        row_vector = [0] * LastGridPoint
        row_vector[i-1] = diag_ele

        PotentialEneryOperator.append(row_vector.copy())
        i = i+1

    M = sp.Matrix(PotentialEneryOperator)
    return M

def getIdentitySympy(order) :
    I = sp.zeros(order, order)

    i = 0
    while i < order:
        I[i, i] = ONE
        i = i + 1

    return I

def generate_pauliBasisCoeff(TwoPowN, N, PauliBasis, Hamiltonian):
    symbols = list(PauliBasis.keys())
    M = len(symbols)
    pauliBasisCoeffs = {}

    if N <= 0 or M <= 0 :
        return {}

    reciTwoPowN = sp.Rational(1, TwoPowN)
   
    indices = [0] * N
    current = [symbols[0]] * N

    finished = False

    while not finished:
        basisString = ''.join(current)
        basisMatrixList = [PauliBasis[gate] for gate in current]
        basisMatrix = TensorProduct(*basisMatrixList)
        basisCoeff = (basisMatrix * Hamiltonian).trace()

        theta = - (ONE / HBAR) * SPATIAL_STEP_SIZE * basisCoeff
        pauliBasisCoeffs[basisString] = {"coeff" : reciTwoPowN * basisCoeff, "matrix":basisMatrix, "theta":theta}

        position = N - 1
        carry = True

        while position >= 0 and carry:
            indices[position] = indices[position] + 1

            if indices[position] == M:
                indices[position] = 0
                current[position] = symbols[0]
                position = position - 1
            else:
                current[position] = symbols[indices[position]]
                carry = False

        if carry:
            finished = True

    return pauliBasisCoeffs

PauliBasis = {"I" : sp.Matrix([ [1,0]
                            ,[0, 1] ]),
        "X" : sp.Matrix([ [0, 1],
                        [1, 0]]),
        "Y" : sp.Matrix([[0, -sp.I],
                        [sp.I, 0]]),
        "Z" : sp.Matrix([[1, 0],
                        [0, sp.I]])}

# Generate Gates based on Qubit Counts

I2 = sp.eye(2)

H = sp.Matrix([
    [1, 1],
    [1, -1]
]) / sp.sqrt(2)

S = sp.Matrix([
    [1, 0],
    [0, sp.I]
])

T = sp.Matrix([
    [1, 0],
    [0, sp.exp(sp.I * sp.pi / 4)]
])


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


BasisGates = [
    ("H(0)", sp.kronecker_product(H, I2)),
    ("H(1)", sp.kronecker_product(I2, H)),

    ("S(0)", sp.kronecker_product(S, I2)),
    ("S(1)", sp.kronecker_product(I2, S)),

    ("T(0)", sp.kronecker_product(T, I2)),
    ("T(1)", sp.kronecker_product(I2, T)),

    ("CNOT(0,1)", get_cnot_01()),
    ("CNOT(1,0)", get_cnot_10())
]


def identityTensors(max) :
    result = []

    i=1
    while i <= max :
        result.append(sp.eye(2**i))
        i = i + 1

    return result

identityTensorList = identityTensors(LastGridPoint-1) # Here N is exponent of Nq = 2^k

def generateGateTensorsSQ(gate, count) :
    results = []

    pos = 1
    while pos <= count:
        currentGate = gate

        leftSide_IdentityCount = pos - 1
        if leftSide_IdentityCount > 0 :
            currentGate = TensorProduct(identityTensorList[leftSide_IdentityCount-1], gate)

        rightSide_IdentityCount = count - pos
        if rightSide_IdentityCount > 0:
            currentGate = TensorProduct(currentGate, identityTensorList[rightSide_IdentityCount-1])

        results.append(currentGate)
        pos = pos + 1

    return results


def main() :
    x = sp.symbols('x')
    potentialFunc = sp.Lambda(x, x**2)

    PotentialEnergyOperator = createPE_operator(potentialFunc)
    KineticEnergyOperator = createKE_operator()

    print(PotentialEnergyOperator)

    print(KineticEnergyOperator)

    HamiltonianOperator = KineticEnergyOperator + PotentialEnergyOperator
    print(HamiltonianOperator)

    # Check if Hamiltonian is Hermitian
    Hamiltonian_H = HamiltonianOperator.H
    if(HamiltonianOperator != Hamiltonian_H):
        print("Cannot Perform Unitary Time Evolution because H is not Hermitian")
        return



    pauliBasisCoeffs = generate_pauliBasisCoeff(LastGridPoint,N, PauliBasis, HamiltonianOperator)
    # print(pauliBasisCoeffs)

    # we can drop II from the list because Measurements results are phase invarient
    pauliBasisCoeffs.pop("II")
    Identity_2Nx2N = getIdentitySympy(LastGridPoint)

    theta = pauliBasisCoeffs["XX"]["theta"]
    basisMatrix = pauliBasisCoeffs["XX"]["matrix"]
    componentMatrix = sp.sin(theta)*Identity_2Nx2N + sp.I*sp.cos(theta)*basisMatrix

    print("\n", componentMatrix, "\n")


    # Run the function
    sequence, final_fidelity, depth = synthesize_circuit_batched(
        target_matrix_sympy=componentMatrix,
        basis_gates_list=BasisGates,
        max_depth=8, 
        target_fidelity=1
    )

    print("-" * 50)
    print(f"Synthesis Complete.")
    print(f"Achieved Fidelity:     {final_fidelity:.6f}")
    print(f"Required Gates:        {depth}")
    print("Optimal Gate Sequence: " + (" -> ".join(sequence) if sequence else "None"))
"""
    for basisString, value in pauliBasisCoeffs.items() :
        theta = value["theta"]

        componentMatrix = sp.sin(theta)*Identity_2Nx2N + sp.I*sp.cos(theta)*value["matrix"]

        # numpyMatrix = np.array(componentMatrix.evalf(15).tolist(), dtype=complex)

        print(basisString)
        print(numpyMatrix, "\n")
"""

    # get numpy matrix

main()
    
# Apply Trrotterization get individuals matrices

# find no. of products required
# Find Circuit Construction for each with fidelity >= 0.9998 atleast
# find time step size

# compare with analytic solution and approximated one

# generate circuit upto a fidlity 0.9999

# create multiple trotterizations and compare

# fault tolereant gates constructions

# create multiple comparisons of systems

# IBM Qiskit and Actual Quantum Deployment for any possible case. Compare Analytic Results and Errors

# Logging system in SQL and JSON

# OOPs and module architecture

# prepare eigenstate of actual system with approximate and find its evolution

# what other results of matrices can be obtained using quantum computing like eigenvalues or eigenvectors etc

