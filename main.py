import math
import random
import matplotlib.pyplot as plt

c = 299792458  # speed of light in m/s
G = 6.67430e-11  # gravitational constant in m^3 kg^-1 s^-2
k_e = 8.9875517873681764e9  # Coulomb constant in N m^2 C^-2
k_s = 1e38  # Strong force constant (arbitrary large value for simplicity)
strong_force_range = 1e-15  # Strong force range in meters (1 femtometer)


class ElectricField:
    def __init__(self):
        self.field_contributions = []

    def update(self, particles):
        self.field_contributions = []
        for particle in particles:
            if particle.charge != 0 and particle.string_mode != 'dark_matter':
                self.field_contributions.append(particle)

    def field_at(self, position):
        Ex, Ey, Ez = 0, 0, 0
        for particle in self.field_contributions:
            dx, dy, dz = [position[i] - particle.positions[-1][i] for i in range(3)]
            r = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
            if r == 0:
                continue
            force_magnitude = k_e * particle.charge / r ** 2
            Ex += force_magnitude * dx / r
            Ey += force_magnitude * dy / r
            Ez += force_magnitude * dz / r
        return [Ex, Ey, Ez]

    def apply(self, particle, time_step, field_strength):
        if particle.string_mode == 'dark_matter':
            return
        force = [particle.charge * f for f in field_strength]
        if particle.mass > 0:
            acceleration = [f / particle.mass for f in force]
            new_velocity = [particle.velocity[i] + acceleration[i] * time_step for i in range(len(particle.velocity))]
            particle.update_velocity(new_velocity)


class MagneticField:
    def __init__(self):
        self.field_contributions = []

    def update(self, particles):
        self.field_contributions = []
        for particle in particles:
            if particle.charge != 0 and particle.mass != 0 and particle.string_mode != 'dark_matter':
                self.field_contributions.append(particle)

    def field_at(self, position):
        Bx, By, Bz = 0, 0, 0
        for particle in self.field_contributions:
            vx, vy, vz = particle.velocity
            dx, dy, dz = [position[i] - particle.positions[-1][i] for i in range(3)]
            r = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
            if r == 0:
                continue
            force_magnitude = k_e * particle.charge / r ** 2
            Bx += force_magnitude * (vy * dz - vz * dy) / r
            By += force_magnitude * (vz * dx - vx * dz) / r
            Bz += force_magnitude * (vx * dy - vy * dx) / r
        return [Bx, By, Bz]

    def apply(self, particle, time_step, field_strength):
        if particle.string_mode == 'dark_matter':
            return
        if len(particle.velocity) != 3:
            raise ValueError("MagneticField currently only supports 3D velocity")
        vx, vy, vz = particle.velocity
        Bx, By, Bz = field_strength
        force_x = particle.charge * (vy * Bz - vz * By)
        force_y = particle.charge * (vz * Bx - vx * Bz)
        force_z = particle.charge * (vx * By - vy * Bx)

        if particle.mass > 0:
            ax = force_x / particle.mass
            ay = force_y / particle.mass
            az = force_z / particle.mass

            new_vx = particle.velocity[0] + ax * time_step
            new_vy = particle.velocity[1] + ay * time_step
            new_vz = particle.velocity[2] + az * time_step
            particle.update_velocity([new_vx, new_vy, new_vz])


class GravitationalField:
    def __init__(self):
        self.field_contributions = []

    def update(self, particles):
        self.field_contributions = []
        for particle in particles:
            if particle.mass != 0:
                self.field_contributions.append(particle)

    def field_at(self, position):
        Fx, Fy, Fz = 0, 0, 0
        for particle in self.field_contributions:
            dx, dy, dz = [position[i] - particle.positions[-1][i] for i in range(3)]
            r = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
            if r == 0:
                continue
            force_magnitude = G * particle.mass / r ** 2
            Fx += force_magnitude * dx / r
            Fy += force_magnitude * dy / r
            Fz += force_magnitude * dz / r
        return [Fx, Fy, Fz]

    def apply(self, particle, time_step, field_strength):
        force = [particle.mass * f for f in field_strength]
        if particle.mass > 0:
            acceleration = [f / particle.mass for f in force]
            new_velocity = [particle.velocity[i] + acceleration[i] * time_step for i in range(len(particle.velocity))]
            particle.update_velocity(new_velocity)


class StrongForceField:
    def __init__(self):
        self.field_contributions = []

    def update(self, particles):
        self.field_contributions = []
        for particle in particles:
            if particle.color_charge != 0:
                self.field_contributions.append(particle)

    def field_at(self, position):
        Sx, Sy, Sz = 0, 0, 0
        for particle in self.field_contributions:
            dx, dy, dz = [position[i] - particle.positions[-1][i] for i in range(3)]
            r = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
            if r == 0 or r > strong_force_range:
                continue
            force_magnitude = k_s * particle.color_charge / r ** 2
            Sx += force_magnitude * dx / r
            Sy += force_magnitude * dy / r
            Sz += force_magnitude * dz / r
        return [Sx, Sy, Sz]

    def apply(self, particle, time_step, field_strength):
        if particle.color_charge == 0:
            return
        force = [particle.color_charge * f for f in field_strength]
        if particle.mass > 0:
            acceleration = [f / particle.mass for f in force]
            new_velocity = [particle.velocity[i] + acceleration[i] * time_step for i in range(len(particle.velocity))]
            particle.update_velocity(new_velocity)


class Particle:
    def __init__(self, name, symbol, mass, charge, spin, string_mode, color_charge=0, decay_constant=0.0, decays_to=None, velocity=None,
                 dimensions=3):
        self.name = name
        self.symbol = symbol
        self.mass = mass
        self.charge = charge
        self.spin = spin
        self.string_mode = string_mode
        self.color_charge = color_charge
        self.decay_constant = decay_constant
        self.decays_to = decays_to if decays_to else []
        self.dimensions = dimensions
        self.velocity = velocity if velocity else [random.uniform(-0.5 * c, 0.5 * c) for _ in range(dimensions)]
        self.energy = self.calculate_energy()
        self.momentum = self.calculate_momentum()
        self.positions = [[0] * dimensions]

    def __repr__(self):
        return f"{self.name} ({self.symbol}): Mass={self.mass} MeV/c^2, Charge={self.charge}e, Spin={self.spin}, String Mode={self.string_mode}, Velocity={self.velocity}, Energy={self.energy} MeV, Momentum={self.momentum} MeV/c"

    def lorentz_factor(self):
        speed_squared = sum(v ** 2 for v in self.velocity)
        if self.mass == 0:
            return 1
        if speed_squared >= c ** 2:
            speed_squared = c ** 2 * 0.99
        return 1 / math.sqrt(1 - speed_squared / c ** 2)

    def calculate_energy(self):
        speed_squared = sum(v ** 2 for v in self.velocity)
        if self.mass == 0:
            return math.sqrt(speed_squared) * c
        gamma = self.lorentz_factor()
        return gamma * self.mass * c ** 2

    def calculate_momentum(self):
        gamma = self.lorentz_factor()
        return [gamma * self.mass * v for v in self.velocity]

    def update_velocity(self, new_velocity):
        speed_squared = sum(v ** 2 for v in new_velocity)
        if speed_squared >= c ** 2:
            scale = math.sqrt(c ** 2 * 0.99 / speed_squared)
            new_velocity = [v * scale for v in new_velocity]
        self.velocity = new_velocity
        self.energy = self.calculate_energy()
        self.momentum = self.calculate_momentum()

    def update_position(self, time_step):
        new_position = [self.positions[-1][i] + self.velocity[i] * time_step for i in range(self.dimensions)]
        self.positions.append(new_position)

    def decay(self, time_step):
        if self.decay_constant == 0:
            return [self]
        decay_probability = 1 - math.exp(-self.decay_constant * time_step)
        if random.random() < decay_probability:
            decay_products = random.choice(self.decays_to)
            total_energy = self.energy
            total_momentum = self.momentum
            num_products = len(decay_products)
            for product in decay_products:
                product.energy = total_energy / num_products
                product.momentum = [m / num_products for m in total_momentum]
                if product.mass != 0:
                    p_velocity_magnitude = math.sqrt(sum(m ** 2 for m in product.momentum) / product.mass ** 2)
                    if p_velocity_magnitude >= c:
                        scale = (c * 0.99) / p_velocity_magnitude
                        p_velocity_magnitude *= scale
                    product.update_velocity([m / product.mass for m in product.momentum])
                else:
                    product.velocity = [m / product.energy * c for m in product.momentum]
            print(f"{self.name} decays into {', '.join([p.symbol for p in decay_products])}")
            return decay_products
        else:
            return [self]

    def collide(self, other):
        if self.string_mode == other.string_mode:
            new_mass = self.mass + other.mass
            new_charge = self.charge + other.charge
            new_velocity = [(self.velocity[i] + other.velocity[i]) / 2 for i in range(self.dimensions)]
            new_particle = Particle(f"NewParticle({self.symbol}+{other.symbol})", f"{self.symbol}{other.symbol}",
                                    new_mass, new_charge, 0, self.string_mode, velocity=new_velocity,
                                    dimensions=self.dimensions)
            print(f"{self.name} and {other.name} collide to form {new_particle.name}")
            return new_particle
        else:
            print(f"{self.name} and {other.name} have different string modes and cannot collide")
            return None

    def scatter(self, other):
        if self.string_mode == other.string_mode:
            self_velocity = self.velocity
            other_velocity = other.velocity
            self.update_velocity(other_velocity)
            other.update_velocity(self_velocity)
            print(f"{self.name} and {other.name} scatter, exchanging velocities")
        else:
            print(f"{self.name} and {other.name} have different string modes and cannot scatter")

    def annihilate(self, other):
        if self.charge + other.charge == 0:
            if self.string_mode == other.string_mode:
                total_energy = self.energy + other.energy
                photon1 = Particle("Photon", "γ", 0, 0, 1, self.string_mode,
                                   velocity=[c if i == 0 else 0 for i in range(self.dimensions)],
                                   dimensions=self.dimensions)
                photon2 = Particle("Photon", "γ", 0, 0, 1, self.string_mode,
                                   velocity=[-c if i == 0 else 0 for i in range(self.dimensions)],
                                   dimensions=self.dimensions)
                photon1.energy = total_energy / 2
                photon2.energy = total_energy / 2
                print(f"{self.name} and {other.name} annihilate to form two photons")
                return [photon1, photon2]
            else:
                print(f"{self.name} and {other.name} have different string modes and cannot annihilate")
                return [self, other]
        else:
            return [self, other]


def simulate_interactions(particles, fields, steps):
    for step in range(steps):
        print(f"\nStep {step + 1}")

        particles_copy = particles[:]
        interactions_done = set()
        for i in range(len(particles_copy)):
            for j in range(i + 1, len(particles_copy)):
                p1 = particles_copy[i]
                p2 = particles_copy[j]
                if (p1, p2) in interactions_done or (p2, p1) in interactions_done:
                    continue
                interactions_done.add((p1, p2))
                if p1.charge + p2.charge == 0 and random.random() < 0.1:
                    new_particles = p1.annihilate(p2)
                    particles.extend(new_particles)
                    if p1 in particles:
                        particles.remove(p1)
                    if p2 in particles:
                        particles.remove(p2)
                elif p1.string_mode == p2.string_mode and random.random() < 0.1:
                    new_particle = p1.collide(p2)
                    if new_particle:
                        particles.append(new_particle)
                        if p1 in particles:
                            particles.remove(p1)
                        if p2 in particles:
                            particles.remove(p2)
                elif p1.string_mode == p2.string_mode and random.random() < 0.1:
                    p1.scatter(p2)

        for field in fields:
            field.update(particles)  # Update field contributions based on particle positions

        for field in fields:
            for particle in particles:
                if isinstance(field, ElectricField) or isinstance(field, MagneticField) or isinstance(field,
                                                                                                      GravitationalField):
                    field_strength = field.field_at(particle.positions[-1])
                    field.apply(particle, 1 / steps, field_strength)

        for particle in particles[:]:
            new_velocity = [v + random.uniform(-0.01 * c, 0.01 * c) for v in particle.velocity]
            particle.update_velocity(new_velocity)
            particle.update_position(1 / steps)
            decay_products = particle.decay(1 / steps)
            if decay_products != [particle]:
                particles.remove(particle)
                particles.extend(decay_products)

        print("Current particles:")
        for particle in particles:
            print(particle)


def visualize_particles(particles):
    fig, ax = plt.subplots()
    for particle in particles:
        positions = list(zip(*particle.positions))
        ax.plot(*positions[:2], label=particle.symbol)  # Only plot in 2D
    ax.set_xlabel("X Position (m)")
    ax.set_ylabel("Y Position (m)")
    ax.set_title("Particle Trajectories")
    ax.legend()
    plt.show()


def electric_field_function(position):
    # Example of a position-dependent electric field
    x, y, z = position
    Ex = 1e5 * math.sin(2 * math.pi * x / 1e6)
    Ey = 1e5 * math.sin(2 * math.pi * y / 1e6)
    Ez = 1e5 * math.sin(2 * math.pi * z / 1e6)
    return [Ex, Ey, Ez]


def magnetic_field_function(position):
    # Example of a position-dependent magnetic field
    x, y, z = position
    Bx = 1 * math.cos(2 * math.pi * x / 1e6)
    By = 1 * math.cos(2 * math.pi * y / 1e6)
    Bz = 1 * math.cos(2 * math.pi * z / 1e6)
    return [Bx, By, Bz]


def simulate():
    electron = Particle("Electron", "e", 0.511, -1, 1 / 2, "mode_known", velocity=[0.1 * c, 0.2 * c, 0.1 * c],
                        dimensions=3)
    electron_neutrino = Particle("Electron Neutrino", "νe", 0.0000022, 0, 1 / 2, "mode_known",
                                 velocity=[0.1 * c, 0.2 * c, 0.1 * c], dimensions=3)
    muon_neutrino = Particle("Muon Neutrino", "νμ", 0.17, 0, 1 / 2, "mode_known", velocity=[0.1 * c, 0.2 * c, 0.1 * c],
                             dimensions=3)
    tau_neutrino = Particle("Tau Neutrino", "ντ", 15.5, 0, 1 / 2, "mode_known", velocity=[0.1 * c, 0.2 * c, 0.1 * c],
                            dimensions=3)
    muon = Particle("Muon", "μ", 105.66, -1, 1 / 2, "mode_known", 0.1,
                    decays_to=[[electron, electron_neutrino, muon_neutrino]], velocity=[0.05 * c, 0.1 * c, 0.05 * c],
                    dimensions=3)
    tau = Particle("Tau", "τ", 1776.86, -1, 1 / 2, "mode_known", 0.2, decays_to=[[muon, muon_neutrino, tau_neutrino]],
                   velocity=[0.01 * c, 0.02 * c, 0.01 * c], dimensions=3)
    photon = Particle("Photon", "γ", 0, 0, 1, "mode_known", velocity=[c, c, c], dimensions=3)
    proton = Particle("Proton", "p", 938.27, 1, 1 / 2, "mode_known", velocity=[0.1 * c, 0.2 * c, 0.1 * c], dimensions=3)
    graviton = Particle("Graviton", "G", 0, 0, 2, "mode_unknown", velocity=[0.5 * c, 0.5 * c, 0.5 * c], dimensions=3)
    axion = Particle("Axion", "a", 1e-5, 0, 0, "mode_unknown", velocity=[0.05 * c, 0.1 * c, 0.05 * c], dimensions=3)
    dark_matter_particle = Particle("Dark Matter", "DM", 1000, 0, 1 / 2, "dark_matter",
                                    velocity=[0.1 * c, 0.1 * c, 0.1 * c], dimensions=3)

    particles = [electron, proton, tau, muon, photon, graviton, axion, dark_matter_particle]
    fields = [ElectricField(), MagneticField(), GravitationalField(), StrongForceField()]
    steps = 10

    print("Simulating interactions between particles in fields:")
    simulate_interactions(particles, fields, steps)
    visualize_particles(particles)


if __name__ == "__main__":
    simulate()
