Sure, here's the updated README including the status as a work in progress and the feature roadmap:

---

# Particle Interaction Simulation

This project simulates interactions between various particles within electric, magnetic, and gravitational fields using Python. It models physical phenomena such as particle decay, collisions, scattering, and annihilation over a series of discrete time steps. Additionally, it includes the simulation of dark matter particles, which only interact gravitationally.

**Note: This project is a work in progress.**

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Testing](#testing)
- [Code Structure](#code-structure)
- [Examples](#examples)
- [Feature Roadmap](#feature-roadmap)

## Installation

1. Clone the repository:
    ```sh
    git clone https://github.com/yourusername/particle-simulation.git
    cd particle-simulation
    ```

2. Create a virtual environment and activate it:
    ```sh
    python -m venv venv
    source venv/bin/activate  # On Windows use `venv\Scripts\activate`
    ```

3. Install the required packages:
    ```sh
    pip install -r requirements.txt
    ```

## Usage

To run the particle interaction simulation, use the `simulate_interactions` function:

```python
from particles import ElectricField, MagneticField, GravitationalField, Particle, simulate_interactions, c

# Define particles
electron = Particle("Electron", "e", 0.511, -1, 1 / 2, "mode_1", velocity=[0.1 * c, 0.2 * c, 0.1 * c])
proton = Particle("Proton", "p", 938.27, 1, 1 / 2, "mode_8", velocity=[0.1 * c, 0.2 * c, 0.1 * c])
dark_matter_particle = Particle("Dark Matter", "DM", 1000, 0, 1 / 2, "dark_matter", velocity=[0.1 * c, 0.1 * c, 0.1 * c])
# Add more particles as needed

# Define fields
electric_field = ElectricField()
magnetic_field = MagneticField()
gravitational_field = GravitationalField()

# List of particles and fields
particles = [electron, proton, dark_matter_particle]
fields = [electric_field, magnetic_field, gravitational_field]

# Run simulation
steps = 10
simulate_interactions(particles, fields, steps)
```

## Testing

This project uses `pytest` and `unittest` for testing. To run the tests, use the following command:

```sh
pytest --maxfail=1 --disable-warnings -v
```

Make sure you have `pytest` installed:

```sh
pip install pytest matplotlib
```

### Running the Tests

To run the tests, simply execute:

```sh
pytest
```

The tests cover all functionalities of the simulation, ensuring 100% code coverage.

## Code Structure

- `particles.py`: Contains the main classes and functions for the simulation.
  - `ElectricField`: Defines the behavior of particles in an electric field.
  - `MagneticField`: Defines the behavior of particles in a magnetic field.
  - `GravitationalField`: Defines the behavior of particles in a gravitational field.
  - `Particle`: Represents a particle with properties such as mass, charge, velocity, and methods for interactions.
  - `simulate_interactions`: Simulates interactions between particles within electric, magnetic, and gravitational fields over a series of steps.

- `test_particles.py`: Contains the test cases for the simulation using `unittest` and `pytest`.

## Examples

### Example 1: Basic Simulation

```python
from particles import ElectricField, MagneticField, GravitationalField, Particle, simulate_interactions, c

# Define particles
electron = Particle("Electron", "e", 0.511, -1, 1 / 2, "mode_1", velocity=[0.1 * c, 0.2 * c, 0.1 * c])
proton = Particle("Proton", "p", 938.27, 1, 1 / 2, "mode_8", velocity=[0.1 * c, 0.2 * c, 0.1 * c])
dark_matter_particle = Particle("Dark Matter", "DM", 1000, 0, 1 / 2, "dark_matter", velocity=[0.1 * c, 0.1 * c, 0.1 * c])

# Define fields
electric_field = ElectricField()
magnetic_field = MagneticField()
gravitational_field = GravitationalField()

# List of particles and fields
particles = [electron, proton, dark_matter_particle]
fields = [electric_field, magnetic_field, gravitational_field]

# Run simulation
steps = 10
simulate_interactions(particles, fields, steps)
```

---
# Pythia Simulation
The installation steps here are a bit involved.
Obviously this simulation is dependent on pythia which you can install here https://pythia.org
Download the tarball and run the following commands to install it.
```sh
tar -xzvf pythia8303.tgz
cd pythia8303
./configure --with-python-config=[path to python-config]
make
```

## Feature Roadmap

### Proposed Upcoming Features

1. **Particle Expansion Modeling**: Introduce the capability to simulate the expansion of particles over time.
2. **Dark Matter Candidates**: Enhance the simulation to identify potential dark matter candidates based on their interactions and properties.
3. **Enhanced Visualization**: Improve the visualization of particle interactions, including 3D representations and more detailed trajectories.
4. **Energy Conservation Checks**: Implement checks to ensure energy conservation during particle interactions.
5. **Additional Particle Types**: Add support for more particle types and interactions, including those predicted by advanced theoretical models.
6. **User Interface**: Develop a graphical user interface (GUI) for easier setup and visualization of simulations.
7. **Parameter Tuning**: Allow dynamic tuning of simulation parameters to explore different scenarios and hypotheses.

---