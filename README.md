# Elastic Pendulum Chain Simulation

This project simulates a chain of **N** elastic particles (spring-mass system) swinging under gravity. It uses symbolic mathematics to derive the exact equations of motion (Lagrangian formalism), generates optimized C code for numerical solving, and visualizes the system in real-time.

## Features

*   **N-Particle Chain**: Simulate any number of connected particles forming a chain.
*   **2D and 3D Modes**: Visualize the chain swinging in a 2D plane or in full 3D space.
*   **High Performance**:
    *   **SymPy** derives the Lagrangian equations automatically.
    *   **C Code Generation**: The physics equations are converted to C code on-the-fly.
    *   **RK4 Solver**: A 4th-order Runge-Kutta integrator written in C handles the time-stepping.
*   **Visualization**: Smooth animation using `matplotlib`.

## Requirements

*   **OS**: Linux/Unix (required for GCC and `.so` shared object generation).
*   **Compiler**: GCC.
*   **Python**: 3.x
*   **Dependencies**: `numpy`, `matplotlib`, `sympy`.

## Setup

1.  **Create a virtual environment** (recommended):
    ```bash
    python3 -m venv venv_matplotlib
    source venv_matplotlib/bin/activate
    ```

2.  **Install dependencies**:
    ```bash
    pip install numpy matplotlib sympy
    ```

## How to Run

Navigate to the project root directory. The main script is `run/particles.py`.

### 1. Basic 2D Chain
To simulate a single elastic pendulum (or default number) in 2D:
```bash
python run/particles.py
```

### 2. Simulation with N Particles
To simulate a chain of `N` particles (e.g., 4 particles), pass the number as an argument:
```bash
python run/particles.py 4
```

### 3. 3D Simulation
To run the simulation in 3D mode, add the `3d` argument. You can combine this with the number of particles:
```bash
# 3D simulation with 3 particles
python run/particles.py 3 3d
```
*(Note: The particles are initialized with slight offsets to induce motion in all three dimensions.)*

## Project Structure

*   **`run/`**: Python source files.
    *   `particles.py`: Main entry point. Orchestrates physics generation, compilation, and animation.
    *   `animation.py`: Handles 2D and 3D plotting/animation.
    *   `lagrangian.py`: Generates C code from SymPy Lagrangian expressions.
    *   `ccompiler.py`: Compiles C source files into shared libraries.
    *   `cprototype.py`: CTypes wrappers for communicating with the C libraries.
*   **`solver/`**: C source files and compiled libraries.
    *   `solver.c`: Implementation of the Runge-Kutta 4 (RK4) solver.
    *   `generated_lagrangian.c`: Auto-generated C code containing the specific equations of motion for your chosen N.
    *   `libsolver.so`, `liblagrangian.so`: Compiled shared libraries.

## Troubleshooting

*   **Compilation Errors**: Ensure you have `gcc` installed and accessible in your system PATH.
*   **Visualization Issues**: If running in a headless environment (like a server or some containers), the window will not appear. Ensure you have an X server running or use X11 forwarding if accessing remotely.
