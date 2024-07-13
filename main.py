import math
import random

import matplotlib.pyplot as plt

c = 299792458  


class ElectricField:
    def __init__(self, field_strength):
        self.field_strength = field_strength

    def apply(self, particle, time_step):
        force_x = particle.charge * self.field_strength[0]
        force_y = particle.charge * self.field_strength[1]

        if particle.mass > 0:
            ax = force_x / particle.mass
            ay = force_y / particle.mass

            new_vx = particle.velocity[0] + ax * time_step
            new_vy = particle.velocity[1] + ay * time_step
            particle.update_velocity(new_vx, new_vy)


class MagneticField:
    def __init__(self, field_strength):
        self.field_strength = field_strength

    def apply(self, particle, time_step):
        vx, vy = particle.velocity
        Bz = self.field_strength[2]
        force_x = particle.charge * vy * Bz
        force_y = -particle.charge * vx * Bz

        if particle.mass > 0:
            ax = force_x / particle.mass
            ay = force_y / particle.mass

            new_vx = particle.velocity[0] + ax * time_step
            new_vy = particle.velocity[1] + ay * time_step
            particle.update_velocity(new_vx, new_vy)


class Particle:
    def __init__(
        self,
        name,
        symbol,
        mass,
        charge,
        spin,
        string_mode,
        decay_constant=0.0,
        decays_to=None,
        velocity=None,
    ):
        self.name = name
        self.symbol = symbol
        self.mass = mass
        self.charge = charge
        self.spin = spin
        self.string_mode = string_mode
        self.decay_constant = decay_constant
        self.decays_to = decays_to if decays_to else []
        self.velocity = (
            velocity
            if velocity
            else (
                random.uniform(-0.5 * c, 0.5 * c),
                random.uniform(-0.5 * c, 0.5 * c),
            )
        )
        self.energy = self.calculate_energy()
        self.momentum = self.calculate_momentum()
        self.positions = [(0, 0)]  

    def __repr__(self):
        return f"{self.name} ({self.symbol}): Mass={self.mass} MeV/c^2, Charge={self.charge}e, Spin={self.spin}, String Mode={self.string_mode}, Velocity={self.velocity}, Energy={self.energy} MeV, Momentum={self.momentum} MeV/c"

    def lorentz_factor(self):
        vx, vy = self.velocity
        speed_squared = vx**2 + vy**2
        if self.mass == 0:
            return 1
        if speed_squared >= c**2:
            speed_squared = c**2 * 0.99
        return 1 / math.sqrt(1 - speed_squared / c**2)

    def calculate_energy(self):
        vx, vy = self.velocity
        speed_squared = vx**2 + vy**2
        if self.mass == 0:
            return math.sqrt(speed_squared) * c
        gamma = self.lorentz_factor()
        return gamma * self.mass * c**2

    def calculate_momentum(self):
        vx, vy = self.velocity
        speed_squared = vx**2 + vy**2
        if self.mass == 0:
            energy = math.sqrt(speed_squared) * c
            return (
                vx / math.sqrt(speed_squared) * energy / c,
                vy / math.sqrt(speed_squared) * energy / c,
            )
        gamma = self.lorentz_factor()
        return gamma * self.mass * vx, gamma * self.mass * vy

    def update_velocity(self, vx, vy):
        speed_squared = vx**2 + vy**2
        if speed_squared >= c**2:
            scale = math.sqrt(c**2 * 0.99 / speed_squared)
            vx *= scale
            vy *= scale
        self.velocity = (vx, vy)
        self.energy = self.calculate_energy()
        self.momentum = self.calculate_momentum()

    def update_position(self, time_step):
        x, y = self.positions[-1]
        vx, vy = self.velocity
        new_x = x + vx * time_step
        new_y = y + vy * time_step
        self.positions.append((new_x, new_y))

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
                product.momentum = (
                    total_momentum[0] / num_products,
                    total_momentum[1] / num_products,
                )
                if product.mass != 0:
                    p_velocity_magnitude = math.sqrt(
                        (product.momentum[0] ** 2 + product.momentum[1] ** 2)
                        / product.mass**2
                    )
                    if p_velocity_magnitude >= c:
                        scale = (c * 0.99) / p_velocity_magnitude
                        p_velocity_magnitude *= scale
                    product.update_velocity(
                        product.momentum[0] / product.mass,
                        product.momentum[1] / product.mass,
                    )
                else:
                    product.velocity = (
                        product.momentum[0] / product.energy * c,
                        product.momentum[1] / product.energy * c,
                    )
            print(
                f"{self.name} decays into {', '.join([p.symbol for p in decay_products])}"
            )
            return decay_products
        else:
            return [self]

    def collide(self, other):
        new_mass = self.mass + other.mass
        new_charge = self.charge + other.charge
        new_velocity = (
            (self.velocity[0] + other.velocity[0]) / 2,
            (self.velocity[1] + other.velocity[1]) / 2,
        )
        new_particle = Particle(
            f"NewParticle({self.symbol}+{other.symbol})",
            f"{self.symbol}{other.symbol}",
            new_mass,
            new_charge,
            0,
            "mode_new",
            velocity=new_velocity,
        )
        print(
            f"{self.name} and {other.name} collide to form {new_particle.name}"
        )
        return new_particle

    def scatter(self, other):
        self_velocity = self.velocity
        other_velocity = other.velocity
        self.update_velocity(*other_velocity)
        other.update_velocity(*self_velocity)
        print(f"{self.name} and {other.name} scatter, exchanging velocities")

    def annihilate(self, other):
        if self.charge + other.charge == 0:
            total_energy = self.energy + other.energy
            photon1 = Particle(
                "Photon", "γ", 0, 0, 1, "mode_γ", velocity=(c, 0)
            )
            photon2 = Particle(
                "Photon", "γ", 0, 0, 1, "mode_γ", velocity=(-c, 0)
            )
            photon1.energy = total_energy / 2
            photon2.energy = total_energy / 2
            print(
                f"{self.name} and {other.name} annihilate to form two photons"
            )
            return [photon1, photon2]
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
                if (p1, p2) in interactions_done or (
                    p2,
                    p1,
                ) in interactions_done:
                    continue
                interactions_done.add((p1, p2))
                if p1.charge + p2.charge == 0 and random.random() < 0.1:
                    new_particles = p1.annihilate(p2)
                    particles.extend(new_particles)
                    if p1 in particles:
                        particles.remove(p1)
                    if p2 in particles:
                        particles.remove(p2)
                elif random.random() < 0.1:
                    new_particle = p1.collide(p2)
                    particles.append(new_particle)
                    if p1 in particles:
                        particles.remove(p1)
                    if p2 in particles:
                        particles.remove(p2)
                elif random.random() < 0.1:
                    p1.scatter(p2)

        for field in fields:
            for particle in particles:
                field.apply(particle, 1 / steps)

        for particle in particles[:]:
            new_vx = particle.velocity[0] + random.uniform(-0.01 * c, 0.01 * c)
            new_vy = particle.velocity[1] + random.uniform(-0.01 * c, 0.01 * c)
            particle.update_velocity(new_vx, new_vy)
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
        positions = zip(*particle.positions)
        ax.plot(*positions, label=particle.symbol)
    ax.set_xlabel("X Position (m)")
    ax.set_ylabel("Y Position (m)")
    ax.set_title("Particle Trajectories")
    ax.legend()
    plt.show()


def simulate():
    electron = Particle(
        "Electron",
        "e",
        0.511,
        -1,
        1 / 2,
        "mode_1",
        velocity=(0.1 * c, 0.2 * c),
    )
    electron_neutrino = Particle(
        "Electron Neutrino",
        "νe",
        0.0000022,
        0,
        1 / 2,
        "mode_2",
        velocity=(0.1 * c, 0.2 * c),
    )
    muon_neutrino = Particle(
        "Muon Neutrino",
        "νμ",
        0.17,
        0,
        1 / 2,
        "mode_3",
        velocity=(0.1 * c, 0.2 * c),
    )
    tau_neutrino = Particle(
        "Tau Neutrino",
        "ντ",
        15.5,
        0,
        1 / 2,
        "mode_4",
        velocity=(0.1 * c, 0.2 * c),
    )
    muon = Particle(
        "Muon",
        "μ",
        105.66,
        -1,
        1 / 2,
        "mode_5",
        0.1,
        decays_to=[[electron, electron_neutrino, muon_neutrino]],
        velocity=(0.05 * c, 0.1 * c),
    )
    tau = Particle(
        "Tau",
        "τ",
        1776.86,
        -1,
        1 / 2,
        "mode_6",
        0.2,
        decays_to=[[muon, muon_neutrino, tau_neutrino]],
        velocity=(0.01 * c, 0.02 * c),
    )
    photon = Particle("Photon", "γ", 0, 0, 1, "mode_7", velocity=(c, c))
    proton = Particle(
        "Proton", "p", 938.27, 1, 1 / 2, "mode_8", velocity=(0.1 * c, 0.2 * c)
    )
    graviton = Particle(
        "Graviton", "G", 0, 0, 2, "mode_9", velocity=(0.5 * c, 0.5 * c)
    )
    axion = Particle(
        "Axion", "a", 1e-5, 0, 0, "mode_10", velocity=(0.05 * c, 0.1 * c)
    )
    particles = [electron, proton, tau, muon, photon, graviton, axion]
    fields = [ElectricField((1e5, 0)), MagneticField((0, 0, 1))]
    steps = 10

    print("Simulating interactions between particles in fields:")
    simulate_interactions(particles, fields, steps)
    visualize_particles(particles)


if __name__ == "__main__":  
    simulate()
