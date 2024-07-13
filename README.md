# Particle Interaction Simulation

This project simulates interactions between various particles within electric and magnetic fields using Python. It models physical phenomena such as particle decay, collisions, scattering, and annihilation over a series of discrete time steps.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Testing](#testing)
- [Code Structure](#code-structure)
- [Examples](#examples)


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
from particles import ElectricField, MagneticField, Particle, simulate_interactions, c

# Define particles
electron = Particle("Electron", "e", 0.511, -1, 1 / 2, "mode_1", velocity=(0.1 * c, 0.2 * c))
proton = Particle("Proton", "p", 938.27, 1, 1 / 2, "mode_8", velocity=(0.1 * c, 0.2 * c))
# Add more particles as needed

# Define fields
electric_field = ElectricField((1e5, 0))
magnetic_field = MagneticField((0, 0, 1))

# List of particles and fields
particles = [electron, proton]
fields = [electric_field, magnetic_field]

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
  - `Particle`: Represents a particle with properties such as mass, charge, velocity, and methods for interactions.
  - `simulate_interactions`: Simulates interactions between particles within electric and magnetic fields over a series of steps.

- `test_particles.py`: Contains the test cases for the simulation using `unittest` and `pytest`.

## Examples

### Example 1: Basic Simulation

```python
from particles import ElectricField, MagneticField, Particle, simulate_interactions, c

# Define particles
electron = Particle("Electron", "e", 0.511, -1, 1 / 2, "mode_1", velocity=(0.1 * c, 0.2 * c))
proton = Particle("Proton", "p", 938.27, 1, 1 / 2, "mode_8", velocity=(0.1 * c, 0.2 * c))

# Define fields
electric_field = ElectricField((1e5, 0))
magnetic_field = MagneticField((0, 0, 1))

# List of particles and fields
particles = [electron, proton]
fields = [electric_field, magnetic_field]

# Run simulation
steps = 10
simulate_interactions(particles, fields, steps)
```
