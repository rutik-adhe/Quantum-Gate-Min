\usepackage{amsmath}

# Quantum Simulation of Hamiltonian Time Evolution

Project Description

---

### (1) Hamiltonian of the System

Schrodinger Time-Dependant Wave Equation

$$ i \hbar \frac{\partial \Psi }{\partial t} = - \frac{\hbar^2}{2m} \nabla^2 \Psi + V(r, t)\Psi  $$

For Time-Independent Potentials and Solving for Stationary States of System, we can write

$$\Psi(r,t) = {\psi(r)} \cdot  {\phi(t)} $$

Substitutating we have

$$ \frac{i \hbar}{\phi}  \frac{\partial \phi }{\partial t} = - \frac{\hbar^2}{2m \psi}  \nabla^2 \psi + V(r) \psi  $$

Left hand side is completely of time domain and right hand side is of space domain, both sides variables are independant of each other. Hence both sides must be constant. Let $E$ denote the constant (Energy)

$$ \left[ i \hbar \frac{\partial  }{\partial t} \right] \phi =  E \phi $$

$$ \left[ - \frac{\hbar^2}{2m} \nabla^2  + V(r) \right] \psi = E \psi $$ 

Solving the Time Dependant Part of Equation, we can write the stationary state of system $\Psi$ as follows

$$ \phi(t) = \phi(0) e^{- \frac{i}{\hbar} E t} $$

Let $\phi(0) = 1$

$$ \Psi(r, t) = e^{- \frac{i}{\hbar} E t} \  \psi(r) $$

For the case of 1D Harmonic Oscillator - Quadratic Potential

$$ V(x) = x^2 $$

Thus, 

$$ - \frac{\hbar^2}{2m} \frac{d^2\psi}{dx^2} + x^2\psi = E \psi $$ 

Here Hamiltonian of the system is

$$ \hat H = - \frac{\hbar^2}{2m} \frac{d^2}{dx^2} + x^2 $$

### (2) Discretizing the Hamiltonian

Using Finite Differences

$$  \left. \frac{d^2\psi}{dx^2} \right|_{x_j}   \approx \frac{\psi_{j+1} + \psi_{j-1} - 2 \psi_{j} }{ (\Delta x)^2 } $$

Taking Spatial Domain $ x \in [-2,2] $ with $4$ Grid Points equally spaced

$x_0=-1.5, \ x_1=-0.5, \ x_2=0.5, \ \ x_3=1.5   $

Grid Spacing $\Delta x = 1$

With Boundary Conditions for $x_{-1} = -2$ and $x_{4}=2$

$\psi_{-1} = 0$ and $\psi_{4}=0$

Then We Let Spatial Wave Function to be a column vector of 4 entries corresponding to the values at each grid point 

$$ \vec{\psi} = \begin{pmatrix}\psi_0 \\ \psi_1 \\ \psi_2 \\ \psi_3 \end{pmatrix} $$

Then Using the finite difference equations

$$  \left. \frac{d^2\psi}{dx^2} \right|_{x_0}   \approx \frac{\psi_{1} + \psi_{-1} - 2 \psi_{0} }{ (\Delta x)^2 } = \psi_{1} - 2 \psi_{0} $$

$$  \left. \frac{d^2\psi}{dx^2} \right|_{x_1}   \approx \frac{\psi_{2} + \psi_{0} - 2 \psi_{1} }{ (\Delta x)^2 } = \psi_{2} +\psi_{0} - 2 \psi_{1} $$

$$  \left. \frac{d^2\psi}{dx^2} \right|_{x_2}   \approx \frac{\psi_{3} + \psi_{1} - 2 \psi_{2} }{ (\Delta x)^2 } = \psi_{3} + \psi_{1} - 2\psi_{2} $$

$$  \left. \frac{d^2\psi}{dx^2} \right|_{x_3}   \approx \frac{\psi_{4} + \psi_{2} - 2 \psi_{3} }{ (\Delta x)^2 } = \psi_{2} - 2\psi_{3} $$

Thus we can write as

$$  \frac{d^2 }{dx^2} \vec{\psi} = \begin{bmatrix} \psi_{1} - 2 \psi_{0} \\ \psi_{2} +\psi_{0} - 2 \psi_{1} \\ \psi_{3} + \psi_{1} - 2\psi_{2} \\ \psi_{2} - 2\psi_{3} \end{bmatrix} $$

$$ \frac{d^2 }{dx^2} \vec{\psi} = \begin{bmatrix} -2 & 1 & 0 & 0 \\ 1 & -2 & 1 & 0 \\ 0 & 1 & -2 & 1 \\ 0 & 0 & 1 & -2 \end{bmatrix} \vec{\psi} $$

Hence the Discrete Matrix form of Kinetic Energy Operator is

$$ [\hat T] = \frac{\hbar^2}{2m} \begin{bmatrix} 2 & -1 & 0 & 0 \\ -1 & 2 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 2 \end{bmatrix} $$

Similarly, Discretizing the Potential we have

$$ V(x_0=-1.5)  = \frac{9}{4} $$
$$ V(x_1=-0.5)  = \frac{1}{4} $$
$$ V(x_2=0.5)  = \frac{1}{4} $$
$$ V(x_3=1.5)  = \frac{9}{4} $$


Thus we can write as

$$ [\hat V]  = \begin{bmatrix}  \frac{9}{4} & 0 & 0 & 0 \\ 0 & \frac{1}{4} & 0 & 0 \\ 0 & 0 &  \frac{1}{4} & 0 \\ 0 & 0 & 0 & \frac{9}{4} \end{bmatrix}  $$

Hence Hamiltonian can be written as :

$$ [\hat H] = [\hat T] + [\hat V] $$

Taking $\hbar = 1 \ unit$ and $2 m = 1 \ unit$ for demonstration purpose. We have

$$ [\hat H] = \begin{bmatrix}  \frac{17}{4} & -1 & 0 & 0 \\ -1 & \frac{9}{4} & 0 & 0 \\ 0 & 0 &  \frac{9}{4} & -1 \\ 0 & 0 & -1 & \frac{17}{4} \end{bmatrix} $$

### (3) Expansion of Hamiltonian in terms of Pauli's Operators

By Pauli Decomposition

The Hamiltonian $(2^n \times 2^n)$ can be expanded in terms of Basis Pauli Matrices

Let $ P_k \in \{I, X, Y, Z\}^{\otimes \ n} $

$$ \hat H = \sum_{k=1}^{4^n} c_k P_k $$

where $c_k = \frac{1}{2^n} Tr(P_k H) $

Substiuating in, we have

$$ [\hat H] = \frac{13}{4} I \otimes I \ - I \otimes X + \frac{13}{8}(1 + i) I \otimes Z $$

$$ + \frac{13}{8}(1 + i) Z \otimes I - \frac{1}{2}(1 + i) Z \otimes X + \frac{9}{8}i Z \otimes Z  $$

We rewrite in $c_{XY}$ & $P_{XY}$ Notation as

$$ [\hat H] = c_{II} P_{II} \ + c_{IX} P_{IX} + c_{IZ} P_{IZ} $$

$$ + c_{ZI} P_{ZI} + c_{ZX} P_{ZX} + c_{ZZ} P_{ZZ}  $$

### (4) Time Evolution of Hamiltonian

$$ \Psi(r, t) = e^{- \frac{i}{\hbar} E t} \  \psi(r) $$

Here

$$ \hat{H} \psi(r) = E \psi(r) $$

Using Taylor Series Expansion $e^x = 1 + \frac{x}{1!} + \frac{x^2}{2!} + ... + \infin $
One can show that

$$ e^{- \frac{i}{\hbar} E t} \psi(r) = e^{- \frac{i}{\hbar} \hat{H} t} \psi(r) $$

Hence we say the Time Evolution Operator acting on $\psi(r)$ is -

$$ \hat{U}(t) = e^{- \frac{i}{\hbar} \hat{H} t} $$

In matrix from we write

$$ [\hat{U}(t)] = e^{- \frac{i}{\hbar} [\hat{H}] t} $$

The Operator $\hat U $ is Unitary if & only if $\hat H$ is Hermitian. (which is True in our case)

Applying the Pauli Decomposition (For our case), we get


$$ [\hat{U}(t)] = \prod_{k} e^{- \frac{i}{\hbar} c_{k} P_{k} t} $$

$$ [\hat{U}(t)] =  e^{- \frac{i}{\hbar} c_{II} P_{II} t} \cdot   e^{- \frac{i}{\hbar} c_{IX} P_{IX} t} \cdot  e^{- \frac{i}{\hbar} c_{IZ} P_{IZ} t} \cdot  e^{- \frac{i}{\hbar} c_{ZI} P_{ZI} t} \cdot  e^{- \frac{i}{\hbar} c_{ZX} P_{ZX} t} \cdot  e^{- \frac{i}{\hbar} c_{ZZ} P_{ZZ} t} $$

### (5) First Order Trotterization - Approximation

$$ [\hat{U}(\Delta t)] \approx \prod_{k} e^{- \frac{i}{\hbar} c_{k} P_{k} \Delta t} $$

The term $P_{II}$ can be neglected upto a global phase factor


$$ \because  e^{i \theta A} = I cos(\theta) + iAsin(\theta) $$

$$ e^{- \frac{i}{\hbar} \frac{13}{4} I_4 \Delta t} = \left[cos\left(\frac{-13 \Delta}{4\hbar }\right) + i sin\left(\frac{-13 \Delta}{4\hbar }\right)\right] I_4 = e^{- \frac{i}{\hbar} \frac{13}{4} \Delta t} I_4 $$

So

$$ [\hat{U}(\Delta t)] \approx e^{- \frac{i}{\hbar} c_{IX} P_{IX} \Delta t} \cdot  e^{- \frac{i}{\hbar} c_{IZ} P_{IZ} \Delta t} \cdot  e^{- \frac{i}{\hbar} c_{ZI} P_{ZI} \Delta t} \cdot  e^{- \frac{i}{\hbar} c_{ZX} P_{ZX} \Delta t} \cdot  e^{- \frac{i}{\hbar} c_{ZZ} P_{ZZ} \Delta t} $$

---

Universal Quantum Gate Set 

Hadamard Gate (H) 
Phase Gate (S)  
T Gate ($\pi$/8)  
Controlled-NOT Gate (CNOT) 
