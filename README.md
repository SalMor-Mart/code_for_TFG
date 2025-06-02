# Evaluating and Programming a 3-Qubit NMR Quantum Computer: The Case of SpinQ Triangulum

This repository contains the code and notebooks developed for my Final Year Project (TFG) titled **"Evaluating and Programming a 3-Qubit NMR Quantum Computer: The Case of SpinQ Triangulum."** The work focuses on evaluating the performance of the SpinQ Triangulum, a 3-qubit Nuclear Magnetic Resonance (NMR) quantum computer, and implementing a novel quantum algorithm on it.

## Project Goals

This project addresses two main objectives:

1. **Evaluation of the SpinQ Triangulum Device**  
   Assess the proper functioning of the Triangulum by analyzing its hardware behavior and verifying the accuracy of built-in quantum algorithms such as Grover’s and Deutsch-Jozsa. Quantum entanglement is checked through Bell inequality tests. Any inconsistencies between theoretical and experimental outcomes are documented, along with possible sources of error and suggestions for improvement.

2. **Implementation of the Quantum Imaginary Time Evolution (QITE) Algorithm**  
   For the first time on this platform, the QITE algorithm is implemented to approximate ground states of Hamiltonians. Different versions of QITE are tested and evaluated based on circuit depth and result accuracy.

## Repository Structure

- `notebooks/` – Jupyter Notebooks with the full QITE implementation and analysis.
- `scripts/` – Python scripts with helper functions for matrix generation, measurement basis construction, and unitary evolution.
- `results/` – Visualizations and circuit representations generated during the algorithm execution.

## Technologies and Libraries

- Python 3.x
- Google Colab (recommended for running the code)
- [Qiskit](https://qiskit.org/)
- NumPy
- QuTiP
- Pandas
- Math

## How to Use

1. **Run on Google Colab**  
   Open the notebook file (`qite_algorithm.ipynb`) in Google Colab for best compatibility. Most dependencies are pre-installed or can be installed within the notebook.

2. **Modify Measurement Parameters**  
   The user can change the measurement bases or the initial state parameters as desired.

3. **Execute the Cells**  
   The code will compute the necessary unitaries at each step of the QITE algorithm and construct the quantum circuit accordingly.

4. **View Results**  
   The output includes the sequence of unitary matrices and a visual representation of the corresponding quantum circuit.

## Results

The main output of this project includes:
- Custom implementation of QITE adapted for SpinQ Triangulum.
- Calculation of unitary operators per QITE iteration step.
- Graphical visualization of quantum circuits based on computed unitaries.
- Analysis of different QITE strategies in terms of depth and convergence.

## Author

**Martí Sales**  
Final Year Project – Bachelor's Degree in Physical Engineering.  
Universitat Politècnica de València, Escuela Técnica Superior de Ingeniería de Telecomunicación.

## License

This repository is shared for academic and educational purposes. You may use or adapt the content under the terms of the MIT License.
