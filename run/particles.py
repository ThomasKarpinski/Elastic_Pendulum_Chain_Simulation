# === IMPORTS ===
from sys import argv
import numpy as np
import sympy as sp
from sympy.physics.mechanics import dynamicsymbols
from ctypes import CFUNCTYPE, POINTER, c_float, c_int, c_size_t, c_void_p, cast, CDLL

# Local imports
import cprototype as cp
import animation as anim
from ccompiler import CSharedLibraryCompiler
from lagrangian import LagrangianToC

# === SETTINGS ===
MODE = "2D"
if "3d" in argv:
    MODE = "3D"

NUMBER_OF_PARTICLES = 1
for arg in argv:
    try:
        NUMBER_OF_PARTICLES = int(arg)
    except:
        pass

dt = 0.01
DIMENSIONS = 2 if MODE == "2D" else 3
ANCHOR = (0, 0, 0) if MODE == "3D" else (0, 0)

# === GENERATE LAGRANGIAN (N-Particle Chain) ===
print(f"--- Generating C Code for {NUMBER_OF_PARTICLES}-Particle Elastic Chain ({MODE}) ---")

# Define Symbols
m, g, k, l0 = sp.symbols('m g k l0')
t = dynamicsymbols._t

# Lists to hold coordinates and energies
coords = []
T_total = 0
V_total = 0

# Helper to get position vector
def get_pos(i):
    if i < 0:
        return sp.Matrix([0, 0, 0]) if MODE == "3D" else sp.Matrix([0, 0])
    
    if MODE == "2D":
        return sp.Matrix([dynamicsymbols(f'x_{i}'), dynamicsymbols(f'y_{i}')])
    else:
        return sp.Matrix([dynamicsymbols(f'x_{i}'), dynamicsymbols(f'y_{i}'), dynamicsymbols(f'z_{i}')])

# Build Lagrangian
for i in range(NUMBER_OF_PARTICLES):
    # Position and Velocity of current particle
    pos = get_pos(i)
    vel = pos.diff(t)
    
    # 1. Kinetic Energy
    # v^2 = dot product of velocity vector with itself
    v_sq = vel.dot(vel)
    T_total += sp.Rational(1, 2) * m * v_sq
    
    # 2. Potential Energy (Gravity)
    # 2D: y is vertical UP usually, or DOWN? 
    # In previous code: V = m*g*y (UP). Let's stick to standard Y-up for 2D, Z-up for 3D.
    if MODE == "2D":
        V_total += m * g * pos[1] 
    else:
        V_total += m * g * pos[2]

    # 3. Potential Energy (Spring connected to previous particle/anchor)
    prev_pos = get_pos(i - 1)
    diff = pos - prev_pos
    dist = sp.sqrt(diff.dot(diff))
    V_total += sp.Rational(1, 2) * k * (dist - l0)**2
    
    # Add coordinates to flat list for the generator
    if MODE == "2D":
        coords.extend([pos[0], pos[1]])
    else:
        coords.extend([pos[0], pos[1], pos[2]])

L_chain = T_total - V_total

# 4. Generate C Code
# This will generate equations for ALL particles.
# For N=4, it will generate q[0]..q[3] logic.
gen = LagrangianToC(L_chain, coords)
consts = {'m': 1.0, 'g': 9.81, 'k': 50.0, 'l0': 1.0} # Stiffer spring (50), shorter length (1.0) for chain
c_code = gen.generate_c_function("elastic_pendulum_step", mode=MODE, constants_map=consts)

# Prepend Struct Definition
if MODE == "2D":
    struct_def = """
#include <math.h>
#include <stddef.h>

typedef struct {
	float x;
	float y;
} Vector2D;
"""
else:
    struct_def = """
#include <math.h>
#include <stddef.h>

typedef struct {
	float x;
	float y;
	float z;
} Vector3D;
"""
full_code = struct_def + c_code

# Write to file
lag_file_path = "solver/generated_lagrangian.c"
with open(lag_file_path, "w") as f:
    f.write(full_code)

print(f"Generated C code written to {lag_file_path}")

# === COMPILE LIBRARIES ===
print("Compiling libsolver.so...")
ccompiler_solver = CSharedLibraryCompiler(source_file="solver/solver.c", output_dir="solver")
solver_path = ccompiler_solver.compile(output_name="libsolver.so")

print("Compiling liblagrangian.so...")
ccompiler_lag = CSharedLibraryCompiler(source_file=lag_file_path, output_dir="solver")
lag_path = ccompiler_lag.compile(output_name="liblagrangian.so")

# === LOAD LIBRARIES ===
_libsolver = cp.EOMSolver(solver_path, NUMBER_OF_PARTICLES, DIMENSIONS=DIMENSIONS)

# Load Lagrangian Lib
_liblag = CDLL(lag_path)

# Define Callback Type
calc_derivs = cast(_liblag.elastic_pendulum_step, c_void_p)

# Bind the callback to next_step
next_step = lambda c, v, nc, nv, dt, n: _libsolver.next_step(c, v, nc, nv, dt, n, calc_derivs)

# === INITIAL CONDITIONS ===
# Initialize as a chain hanging down
positions = []
velocities = []

start_x = 0.0
start_y = 0.0
start_z = 0.0

separation = 1.0 # approx l0

for i in range(NUMBER_OF_PARTICLES):
    if MODE == "2D":
        # Hang along negative Y
        p_x = start_x + 0.5 # slight offset to swing
        p_y = start_y - (i + 1) * separation
        positions.append(_libsolver.vector(x=p_x, y=p_y))
        velocities.append(_libsolver.vector(x=0.0, y=0.0))
    else:
        # Hang along negative Z with some X and Y offset to create 3D motion
        p_x = start_x + 0.5
        p_y = start_y + 0.5  # Added Y offset
        p_z = start_z - (i + 1) * separation
        positions.append(_libsolver.vector(x=p_x, y=p_y, z=p_z))
        velocities.append(_libsolver.vector(x=0.0, y=0.0, z=0.0))

# === RUN ANIMATION ===
print(f"Starting Animation ({MODE} Chain Visualization)...")
if MODE == "2D":
    ani = anim.Animation2D(vector_factory=_libsolver.vector,
                           c_arr=_libsolver.c_arr,
                           next_step=next_step,
                           positions=positions,
                           velocities=velocities,
                           dt=dt,
                           NUMBER_OF_PARTICLES=NUMBER_OF_PARTICLES,
                           DIMENSIONS=DIMENSIONS,
                           anchor_point=ANCHOR)
    # Auto-scale limits based on N
    lim = NUMBER_OF_PARTICLES * 1.5
    ani.create_canvas(xlabel='x', ylabel='y', xlim=[-lim, lim], ylim=[-lim, 1])
else:
    ani = anim.Animation3D(vector_factory=_libsolver.vector,
                           c_arr=_libsolver.c_arr,
                           next_step=next_step,
                           positions=positions,
                           velocities=velocities,
                           dt=dt,
                           NUMBER_OF_PARTICLES=NUMBER_OF_PARTICLES,
                           DIMENSIONS=DIMENSIONS,
                           anchor_point=ANCHOR)
    lim = NUMBER_OF_PARTICLES * 1.5
    ani.create_canvas(xlabel='x', ylabel='y', zlabel='z', xlim=[-lim, lim], ylim=[-lim, lim], zlim=[-lim, 1])

if "save" in argv:
    gif_name = f"simulation_{MODE.lower()}.gif"
    # Use more frames for the GIF to show a good loop
    ani.run_animation(frames=120, save_filename=gif_name)
else:
    ani.run_animation()