import os
import sys
cfg = open("../pythia8312/Makefile.inc")
lib = "../lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)

#==========================================================================

# Import the Pythia module.

import pythia8
import numpy as np
from vispy.scene import visuals
from vispy import scene, app
from vispy.app import run
import sys
import matplotlib.pyplot as plt
from particle import Particle as NamedParticle
from particle import PDGID
from vispy.io import write_png

import imageio
import glob


def plot_point_cloud(points_matrix):
    # Create a canvas with a simple view
    canvas = scene.SceneCanvas(keys='interactive', show=True)
    view = canvas.central_widget.add_view()

    # Create scatter object and fill in the data
    scatter = visuals.Markers()
    scatter.set_data(points_matrix, edge_color=None, face_color=(1, 1, 1, 0.5), size=25)

    view.add(scatter)

    # Set the camera to 'turntable' for 3D navigation
    view.camera = 'turntable'

    # Run the application
    if sys.flags.interactive != 1:
        run()


def normalize_data(points_matrix, momentum_matrix, target_max=1.0):
    # Calculate the maximum absolute value for scaling
    max_val_points = np.max(np.abs(points_matrix))
    max_val_momentum = np.max(np.abs(momentum_matrix))

    # Normalize the matrices to the target range [-target_max, target_max]
    normalized_points = (points_matrix / max_val_points) * target_max
    normalized_momentum = (momentum_matrix / max_val_momentum) * target_max

    return normalized_points, normalized_momentum


def log_normalize_data(points_matrix, momentum_matrix, target_max=1.0):
    epsilon = 1e-9  # Small constant to ensure no log(0)

    # Calculate the magnitude of each vector
    magnitudes = np.linalg.norm(momentum_matrix, axis=1, keepdims=True)

    # Shift magnitudes to ensure all values are positive before log transformation
    shifted_magnitudes = magnitudes + epsilon

    # Apply log normalization
    log_magnitudes = np.log(shifted_magnitudes)

    # Normalize the log-transformed magnitudes to the target range [-target_max, target_max]
    max_log_magnitude = np.max(np.abs(log_magnitudes))
    normalized_magnitudes = (log_magnitudes / max_log_magnitude) * target_max

    # Recreate the momentum vectors with normalized magnitudes while preserving direction
    unit_vectors = momentum_matrix / magnitudes
    normalized_momentum_matrix = unit_vectors * normalized_magnitudes

    # normalize the points matrix
    normalized_points_matrix = (points_matrix / np.max(np.abs(points_matrix))) * target_max


    return normalized_points_matrix, normalized_momentum_matrix


def shift_particles_from_origin(points_matrix, momentum_matrix, shift_magnitude=0.0001):
    # Identify particles at the origin
    origin_indices = np.where(np.all(points_matrix == 0, axis=1))[0]

    # Calculate unit norm vectors for their momentum
    momentum_norms = np.linalg.norm(momentum_matrix[origin_indices], axis=1, keepdims=True)
    unit_momentum_vectors = momentum_matrix[origin_indices] / momentum_norms

    # Shift particles at the origin along the unit norm of their momentum vector
    points_matrix[origin_indices] += unit_momentum_vectors * shift_magnitude

    return points_matrix


def generate_positions_through_momentum(position_matrix, momentum_matrix, categories, lifetimes,
                                        spawntimes, time_schedule, step_size):
    """
    Generates a list of positions by taking steps in each momentum vector.

    Parameters:
    - position_matrix: np.array, initial positions of the particles.
    - momentum_matrix: np.array, momentum vectors for the particles.
    - categories: np.array, categories of the particles.
    - lifetimes: np.array, lifetimes of the particles. mm/c
    - spawntimes: np.array, spawntimes of the particles. mm/c
    - num_steps: int, number of steps to simulate.
    - step_size: float, size of each step.

    Returns:
    - List[(np.array, np.array)], a list of position matrices after each step and the corresponding categories for
        colouring
    """



    # modify lifetimes so any particle that has a lifetime of
    # 0 is set a value that guarantees they will survive for 30 steps
    lifetimes = np.where(lifetimes == 0, 30*step_size, lifetimes)

    positions_list = [(position_matrix.copy(), categories.copy())]

    for current_time in time_schedule:
        age = current_time - spawntimes
        alive = age < lifetimes
        new_positions = positions_list[-1][0] + np.where(alive[:, None], momentum_matrix, 0) * step_size
        categories_of_current_particles = np.where(alive, categories, 0)
        positions_list.append((new_positions, categories_of_current_particles))

    positions_list = positions_list[1:]

    # # non time dependent one
    # for _ in range(num_steps):
    #     # Update positions by adding momentum * step_size
    #     new_positions = positions_list[-1] + momentum_matrix * step_size
    #     positions_list.append(new_positions)

    return positions_list

def generate_prequel_positions(positions_matrix, momentum_matrix, num_steps, step_size):
    """
    Generates a list of positions by taking steps in each momentum vector.
    This works backwards and steps up to reach the positions in the final positions matrix
    :param positions_matrix:
    :param momentum_matrix:
    :param num_steps:
    :param step_size:
    :return: - List[np.array], a list of position matrices after each step.
    """
    prequel_positions = [positions_matrix]
    for _ in range(num_steps):
        # Calculate new positions by subtracting momentum * step_size
        new_positions = prequel_positions[0] - momentum_matrix * step_size
        # Prepend new positions to the list
        prequel_positions.insert(0, new_positions)

    return prequel_positions



def categorize_particle_numerically(pdgid):
    pdgid_instance = PDGID(pdgid)

    # Assign a unique integer to each category
    if pdgid_instance.is_lepton:
        return 1  # Lepton
    elif pdgid_instance.is_quark:
        return 2  # Quark
    elif pdgid_instance.is_gauge_boson_or_higgs:
        return 3  # Gauge Boson or Higgs
    elif pdgid_instance.is_hadron:
        if pdgid_instance.is_meson:
            return 4  # Meson
        elif pdgid_instance.is_baryon:
            return 5  # Baryon
        else:
            return 6  # Other Hadron
    else:
        return 0  # Unknown


def plot_particles_with_momentum(points_matrix, momentum_matrix):

    # Create a canvas with a simple view
    canvas = scene.SceneCanvas(keys='interactive', show=True)
    view = canvas.central_widget.add_view()

    # Generate a unique color for each unique particle ID
    unique_ids = np.unique(id_array)
    colors = plt.cm.jet(np.linspace(0, 1, len(unique_ids)))  # Using matplotlib's colormap
    color_dict = {id_: color for id_, color in zip(unique_ids, colors)}

    # Apply colors based on particle IDs
    particle_colors = np.array([color_dict[id_] for id_ in id_array])

    # Plot particles
    scatter = visuals.Markers()
    scatter.set_data(points_matrix, edge_color=None, face_color=particle_colors, size=10)
    view.add(scatter)

    # Add arrows for momentum
    for start_point, momentum, id_ in zip(points_matrix, momentum_matrix, id_array):
        end_point = start_point + momentum
        arrow_color = list(color_dict[id_])  # Convert color to list to modify
        arrow_color[3] = 0.25 # Set alpha value to 0.5 for semi-transparency
        arrow = visuals.Arrow(pos=np.array([start_point, end_point]), color=arrow_color, arrow_size=5,
                              arrow_type='triangle_60', parent=view.scene)

    # Set the camera
    view.camera = 'turntable'

    # Run the application
    if __name__ == '__main__':
        run()

def animate_particle_movement(all_positions):
    """
    Animates the movement of particles.

    """
    # check if the gif_frames directory exists, if not create it
    if not os.path.exists("gif_frames"):
        os.makedirs("gif_frames")
    # clear out the gif_frames directory
    for file in glob.glob(os.path.join("gif_frames", "*")):
        os.remove(file)

    # Create a canvas with a simple view
    canvas = scene.SceneCanvas(keys='interactive', show=True, size=(1920, 1080))
    view = canvas.central_widget.add_view()

    # Map categories to colors using a colormap
    unique_categories = np.arange(7)  # 7 unique categories
    colors = plt.cm.viridis(np.linspace(0, 1, len(unique_categories)))
    color_map = {category: color for category, color in zip(unique_categories, colors)}



    # Initial plot setup with the first set of positions
    scatter = visuals.Markers()
    scatter.set_data(all_positions[0], edge_color=None, face_color=(1, 1, 1, 0.5), size=7)
    view.add(scatter)

    # Set the camera to 'turntable' for 3D navigation
    view.camera = 'turntable'
    view.camera.azimuth = 30

    # determine prequel steps and main steps
    prequel_steps = 0
    main_steps = 0
    for i, positions in enumerate(all_positions):
        if isinstance(positions, np.ndarray):
            prequel_steps += 1
        else:
            break

    main_steps = len(all_positions) - prequel_steps

    # create a prequel distance schedule where it starts at 10000 and rapidly decreases to 1e-3 near the end
    #prequel_distance_schedule = np.logspace(1, -1, prequel_steps)
    camera_elevation_schedule = np.sin(np.linspace(0, 3*np.pi/2, len(all_positions))) * 30



    # Define the update function for the animation
    def update(ev):
        nonlocal scatter, all_positions
        current_time_step = update.current_time_step
        if isinstance(all_positions[current_time_step % len(all_positions)], np.ndarray):
            # this is a prequel scene, no color necessary
            scatter.set_data(all_positions[current_time_step % len(all_positions)], edge_color=None, face_color=(1, 1, 1, 0.5), size=20)
            # set the camera to look at the mean position of the particles and slowly zoom in
            mean_position = np.mean(all_positions[current_time_step % len(all_positions)], axis=0)
            view.camera.center = mean_position
            view.camera.distance = 1
            #view.camera.distance = prequel_distance_schedule[current_time_step % len(all_positions)]


        else:
            # Generate an array of colors for each particle based on its category
            categories = all_positions[current_time_step % len(all_positions)][1]
            particle_colors = np.array([color_map[category] for category in categories])
            # make each particle a bit transparent
            particle_colors[:, 3] = 0.75
            scatter.set_data(all_positions[current_time_step % len(all_positions)][0], edge_color=None, face_color=particle_colors, size=13)
            view.camera.distance = 5e-5


        #Rotate the camera by incrementing the azimuth angle
        view.camera.azimuth += 2
        view.camera.elevation = camera_elevation_schedule[current_time_step % len(all_positions)]
        if view.camera.azimuth >= 360:
            view.camera.azimuth -= 360


        # Capture the frame
        img = canvas.render()
        write_png(os.path.join("gif_frames", f"frame_{current_time_step:04d}.png"), img)

        update.current_time_step += 1

    update.current_time_step = 0

    # Create a timer that calls the update function
    timer = app.Timer()
    timer.connect(update)
    timer.start(0.07)  # Adjust the interval for faster or slower animation
    app.run()
    frames = sorted(glob.glob(os.path.join("gif_frames", "*.png")))
    images = [imageio.imread(frame) for frame in frames]
    imageio.mimsave('animation.gif', images, fps=25)  # Adjust fps as needed

energy_scale = 1
energyFactor = 2000
angle = np.pi / 24  # glancing angle
momentum_p1 = np.array([energy_scale * np.cos(angle), energy_scale * np.sin(angle), 0]) * energyFactor
momentum_p2 = np.array([-energy_scale, 0, 0]) * energyFactor


nevents = 10
pythia = pythia8.Pythia()
pythia.readString("Beams:idA = 2212")
pythia.readString("Beams:idB = 2212")
pythia.readString("Beams:frameType = 3")
pythia.readString(f"Beams:pxA = {momentum_p1[0]}")
pythia.readString(f"Beams:pyA = {momentum_p1[1]}")
pythia.readString(f"Beams:pzA = {momentum_p1[2]}")
pythia.readString(f"Beams:pxB = {momentum_p2[0]}")
pythia.readString(f"Beams:pyB = {momentum_p2[1]}")
pythia.readString(f"Beams:pzB = {momentum_p2[2]}")

pythia.readString("HardQCD:all = on")
pythia.readString("SoftQCD:all = on")






pythia.init()
# Begin event loop. Generate event. Skip if error. List first one.
for iEvent in range(nevents):
    if not pythia.next(): continue
    numParticles = pythia.event.size()
    if numParticles != 0:
        points_matrix = np.zeros((numParticles, 3))
        momentum_matrix = np.zeros((numParticles, 3))
        categories = np.zeros(numParticles)
        lifetimes = np.zeros(numParticles)
        spawnTimes = np.zeros(numParticles)
        id_array = np.zeros(numParticles)
        storedParticles = 0
        for i in range(numParticles):
            particle = pythia.event[i]
            print("Particle", i, "has PID = ", particle.id())
            id = particle.id()
            try:
                named_particle = NamedParticle.from_pdgid(id)
                particle_name = named_particle.name
            except Exception as e:
                # didnt find it oh well, name isnt critical
                pass

            p = particle.p()
            mass = particle.m()
            vProd = particle.vProd()
            tau = particle.tau()
            survived = particle.isFinal()
            decaySpot = particle.vDec()
            # store the starting point of the particle
            points_matrix[i] = np.array([vProd[0], vProd[1], vProd[2]])
            # store the momentum of the particle
            momentum_matrix[i] = np.array([p[0], p[1], p[2]])
            id_array[i] = id
            categories[i] = categorize_particle_numerically(id)
            lifetimes[i] = tau
            spawnTimes[i] = vProd[3]
            y = particle.y()
            velocity = momentum_matrix[i] / (y*mass)
            j=1

        collisionParticles = np.array([1, 2])

        points_matrix_normed, momentum_matrix_normed = log_normalize_data(points_matrix, momentum_matrix)
        points_matrix_normed = shift_particles_from_origin(points_matrix_normed, momentum_matrix_normed,
                                                           shift_magnitude=0.0000001)
        #plot_particles_with_momentum(points_matrix_normed, momentum_matrix_normed)
        forward_steps = 100
        # create a time schedule for the overall simulation
        time_schedule = np.linspace(0, np.max(spawnTimes) * 1.1, forward_steps)
        step_size = time_schedule[1] - time_schedule[0]

        prequel_positions = generate_prequel_positions(points_matrix[collisionParticles],
                                                       momentum_matrix[collisionParticles],
                                                       num_steps=50, step_size=step_size)
        all_positions = generate_positions_through_momentum(points_matrix_normed, momentum_matrix_normed, categories,
                                                            lifetimes, spawnTimes, time_schedule, step_size)
        all_positions = prequel_positions + all_positions
        animate_particle_movement(all_positions)

# End of event loop. Statistics. Histogram. Done
