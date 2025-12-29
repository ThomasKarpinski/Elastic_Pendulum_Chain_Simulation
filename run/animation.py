# === IMPORTS ===
# Standard library imports
import itertools as it

# Numpy (https://numpy.org/)
# and ctypes (https://docs.python.org/3/library/ctypes.html)
import numpy as np

# Matplotlib (https://matplotlib.org/) 
# imports for plotting and animation
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.animation as animation


class Animation2D:
    def __init__(self,vector_factory=None, c_arr=None,
                 next_step=None,positions=None, velocities=None,
                 dt=0.01,NUMBER_OF_PARTICLES=1, DIMENSIONS=2, anchor_point=None):
        self.data = np.zeros((NUMBER_OF_PARTICLES, DIMENSIONS))
        self.DIMENSIONS = DIMENSIONS
        self.anchor_point = anchor_point

        self.vector = vector_factory
        self.c_arr = c_arr
        self.next_step = next_step
        self.set_positions(positions)
        self.set_velocities(velocities)
        self.colours = it.cycle(mcolors.TABLEAU_COLORS)
        self.dt = dt
        self.NUMBER_OF_PARTICLES = NUMBER_OF_PARTICLES

    def set_positions(self, positions):
        self.positions = positions
        for i, pos in enumerate(self.positions):
            pos(self.data, i)  # Initialize 'self.data' with starting positions

    def set_velocities(self, velocities):
        self.velocities = velocities

    def create_canvas(self,**kwargs):
        """
        This function sets up the Matplotlib figure and axes for the animation.
        It initializes the plot elements that will be updated in each frame.
        """
        self.fig, self.ax = plt.subplots()
        self.ax.set_aspect('equal', adjustable='box')  # Ensure x and y axes have the same scale
        # Set plot limits and labels
        if 'xlim' not in kwargs:
            kwargs['xlim']=[-2.1, 2.1]
        if 'ylim' not in kwargs:
            kwargs['ylim']=[-2.1, 2.1]
        if 'xlabel' not in kwargs:
            kwargs['xlabel']='x'
        if 'ylabel' not in kwargs:
            kwargs['ylabel']='y'
        self.ax.set(**kwargs)
        if self.NUMBER_OF_PARTICLES > 1:
            # Draw chain lines (open loop)
            self.lines = self.ax.plot(self.data[:, 0],
                                      self.data[:, 1], lw=1)[0]
        
        # Draw anchor point and spring if provided
        self.spring_line = None
        if self.anchor_point is not None:
             self.ax.scatter(*self.anchor_point[:2], color='black', marker='x', s=100, label='Anchor')
             # Connect anchor to the first particle (or all? usually 1st for pendulum)
             # Assuming single pendulum or chain starting at anchor
             self.spring_line = self.ax.plot([self.anchor_point[0], self.data[0, 0]],
                                             [self.anchor_point[1], self.data[0, 1]], 
                                             color='gray', linestyle='--', lw=1)[0]

        # 'points' is a scatter plot of the particles themselves
        self.points = self.ax.scatter(self.data[:, 0], self.data[:, 1],
                      c=[clr for clr, _ in zip(self.colours, range(self.NUMBER_OF_PARTICLES))], s=57)

    def update_frame(self, frame):
        """
        This function is called for each frame of the animation.
        It calculates the new state of the simulation and updates the plot.
        """
        global positions,velocities

        # Create empty Vector objects to hold the C function results
        new_positions = [self.vector(x=0, y=0, z=0) for i in range(self.NUMBER_OF_PARTICLES)]
        new_velocities= [self.vector(x=0, y=0, z=0) for i in range(self.NUMBER_OF_PARTICLES)]

        c_positions       = self.c_arr(*self.positions)
        c_velocities      = self.c_arr(*self.velocities)
        c_new_positions   = self.c_arr(*new_positions)
        c_new_velocities  = self.c_arr(*new_velocities)

        # 1. Calculate the new positions and velocities
        self.next_step(c_positions, c_velocities,
                       c_new_positions, c_new_velocities,
                       self.dt, self.NUMBER_OF_PARTICLES)
        self.positions  = c_positions[:]
        self.velocities = c_velocities[:]
        new_positions   = c_new_positions[:]
        new_velocities  = c_new_velocities[:]

        for i,new_position in enumerate(new_positions):
            # 2. Update the master Python lists with the new state
            self.positions[i]  = new_position
            self.velocities[i] = new_velocities[i]

            # 3. Update the NumPy plotting array
            new_position(self.data, i)

        # --- Update Matplotlib elements ---
        # Update the positions of the scattered points
        self.update_plot_elements()

    def update_plot_elements(self):
        self.points.set_offsets(self.data[:, :2])
        if self.NUMBER_OF_PARTICLES > 1:
            self.lines.set_xdata(self.data[:, 0])
            self.lines.set_ydata(self.data[:, 1])
        
        if self.spring_line is not None:
            self.spring_line.set_xdata([self.anchor_point[0], self.data[0, 0]])
            self.spring_line.set_ydata([self.anchor_point[1], self.data[0, 1]])

    def run_animation(self, frames=60, interval=30, save_filename=None):
        self.ani = animation.FuncAnimation(fig=self.fig, func=self.update_frame,
                                           frames=frames, interval=interval)
        if save_filename:
            print(f"Saving animation to {save_filename}...")
            fps = 1000 / interval
            self.ani.save(save_filename, writer='pillow', fps=fps)
            print(f"Animation saved to {save_filename}")
        else:
            plt.show()

class Animation3D(Animation2D):
    def create_canvas(self, **kwargs):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        
        if 'xlim' not in kwargs: kwargs['xlim'] = [-4, 4]
        if 'ylim' not in kwargs: kwargs['ylim'] = [-4, 4]
        if 'zlim' not in kwargs: kwargs['zlim'] = [-4, 4]
        if 'xlabel' not in kwargs: kwargs['xlabel'] = 'x'
        if 'ylabel' not in kwargs: kwargs['ylabel'] = 'y'
        if 'zlabel' not in kwargs: kwargs['zlabel'] = 'z'
        
        self.ax.set(**kwargs)
        
        # Draw anchor point and spring if provided
        self.spring_line = None
        if self.anchor_point is not None:
             self.ax.scatter(*self.anchor_point, color='black', marker='x', s=100, label='Anchor')
             # Connect anchor to the first particle
             self.spring_line = self.ax.plot([self.anchor_point[0], self.data[0, 0]],
                                             [self.anchor_point[1], self.data[0, 1]], 
                                             [self.anchor_point[2], self.data[0, 2]],
                                             color='gray', linestyle='--', lw=1)[0]
        
        clrs = [clr for clr, _ in zip(self.colours, range(self.NUMBER_OF_PARTICLES))]
        self.points = self.ax.scatter(self.data[:, 0], self.data[:, 1], self.data[:, 2], c=clrs, s=57)
        
        if self.NUMBER_OF_PARTICLES > 1:
            self.lines = self.ax.plot(self.data[:, 0],
                                      self.data[:, 1],
                                      self.data[:, 2], lw=1)[0]

    def update_plot_elements(self):
        # Update scatter points
        self.points._offsets3d = (self.data[:, 0], self.data[:, 1], self.data[:, 2])
        
        if self.spring_line is not None:
            self.spring_line.set_data([self.anchor_point[0], self.data[0, 0]],
                                      [self.anchor_point[1], self.data[0, 1]])
            self.spring_line.set_3d_properties([self.anchor_point[2], self.data[0, 2]])
        
        if self.NUMBER_OF_PARTICLES > 1:
            self.lines.set_data(self.data[:, 0],
                                self.data[:, 1])
            self.lines.set_3d_properties(self.data[:, 2])


