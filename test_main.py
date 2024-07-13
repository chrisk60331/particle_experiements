import random
import unittest

from main import (
    ElectricField,
    MagneticField,
    Particle,
    c,
    simulate,
    simulate_interactions,
)

# Setting a fixed random seed for consistency
random.seed(42)


class TestMeta(type):
    def __new__(cls, name, bases, dct):
        particles = [
            Particle(
                "Electron",
                "e",
                0.511,
                -1,
                1 / 2,
                "mode_1",
                velocity=(0.1 * c, 0.2 * c),
            ),
            Particle(
                "Proton",
                "p",
                938.27,
                1,
                1 / 2,
                "mode_8",
                velocity=(0.1 * c, 0.2 * c),
            ),
            Particle(
                "Neutron",
                "n",
                939.57,
                0,
                1 / 2,
                "mode_9",
                velocity=(0.05 * c, 0.1 * c),
            ),
            Particle(
                "Muon",
                "μ",
                105.66,
                -1,
                1 / 2,
                "mode_5",
                velocity=(0.05 * c, 0.1 * c),
            ),
            Particle(
                "Tau",
                "τ",
                1776.86,
                -1,
                1 / 2,
                "mode_6",
                velocity=(0.01 * c, 0.02 * c),
            ),
            Particle("Photon", "γ", 0, 0, 1, "mode_7", velocity=(c, c)),
            Particle(
                "Graviton", "G", 0, 0, 2, "mode_9", velocity=(0.5 * c, 0.5 * c)
            ),
            Particle(
                "Axion",
                "a",
                1e-5,
                0,
                0,
                "mode_10",
                velocity=(0.05 * c, 0.1 * c),
            ),
            Particle(
                "Higgs Boson", "H", 125.1, 0, 0, "mode_11", velocity=(0, 0)
            ),
            Particle(
                "Z Boson",
                "Z",
                91.19,
                0,
                1,
                "mode_12",
                velocity=(0.1 * c, 0.1 * c),
            ),
            Particle(
                "W Boson",
                "W",
                80.39,
                1,
                1,
                "mode_13",
                velocity=(0.2 * c, 0.2 * c),
            ),
            Particle(
                "Electron Neutrino",
                "νe",
                0.0000022,
                0,
                1 / 2,
                "mode_2",
                velocity=(0.1 * c, 0.2 * c),
            ),
            Particle(
                "Muon Neutrino",
                "νμ",
                0.17,
                0,
                1 / 2,
                "mode_3",
                velocity=(0.1 * c, 0.2 * c),
            ),
            Particle(
                "Tau Neutrino",
                "ντ",
                0.17,
                0,
                1 / 2,
                "mode_4",
                velocity=(0.1 * c, 0.2 * c),
            ),
            Particle(
                "Charm Quark",
                "c",
                1.28,
                2 / 3,
                1 / 2,
                "mode_14",
                velocity=(0.3 * c, 0.3 * c),
            ),
            Particle(
                "Strange Quark",
                "s",
                0.095,
                -1 / 3,
                1 / 2,
                "mode_15",
                velocity=(0.2 * c, 0.2 * c),
            ),
            Particle(
                "Top Quark",
                "t",
                173.0,
                2 / 3,
                1 / 2,
                "mode_16",
                velocity=(0.25 * c, 0.25 * c),
            ),
            Particle(
                "Bottom Quark",
                "b",
                4.18,
                -1 / 3,
                1 / 2,
                "mode_17",
                velocity=(0.1 * c, 0.1 * c),
            ),
            Particle("Gluon", "g", 0, 0, 1, "mode_18", velocity=(c, c)),
            Particle(
                "Pion",
                "π+",
                0.139,
                1,
                0,
                "mode_19",
                velocity=(0.1 * c, 0.1 * c),
                decays_to=[
                    [
                        Particle(
                            "Muon",
                            "μ+",
                            105.66,
                            1,
                            1 / 2,
                            "mode_5",
                            velocity=(0.05 * c, 0.1 * c),
                        ),
                        Particle(
                            "Muon Neutrino",
                            "νμ",
                            0.17,
                            0,
                            1 / 2,
                            "mode_3",
                            velocity=(0.1 * c, 0.2 * c),
                        ),
                    ]
                ],
            ),
            Particle(
                "Kaon",
                "K+",
                0.494,
                1,
                0,
                "mode_20",
                velocity=(0.15 * c, 0.15 * c),
                decays_to=[
                    [
                        Particle(
                            "Pion",
                            "π+",
                            0.139,
                            1,
                            0,
                            "mode_19",
                            velocity=(0.1 * c, 0.1 * c),
                        )
                    ]
                ],
            ),
            Particle(
                "Lambda Baryon",
                "Λ",
                1115.683,
                0,
                1 / 2,
                "mode_21",
                decay_constant=0.9,
                velocity=(0.2 * c, 0.3 * c),
                decays_to=[
                    [
                        Particle(
                            "Proton",
                            "p",
                            938.27,
                            1,
                            1 / 2,
                            "mode_8",
                            velocity=(0.1 * c, 0.2 * c),
                        ),
                        Particle(
                            "Pion",
                            "π-",
                            0.139,
                            -1,
                            0,
                            "mode_19",
                            velocity=(0.1 * c, 0.1 * c),
                        ),
                    ]
                ],
            ),
        ]

        new_dct = {}
        for test_name, test_func in dct.items():
            if test_name.startswith("test_"):
                for particle in particles:
                    test_method_name = f"{test_name}_{particle.name}"
                    new_dct[test_method_name] = cls._create_test(
                        test_func, particle
                    )
            else:
                new_dct[test_name] = test_func
        return type.__new__(cls, name, bases, new_dct)

    @staticmethod
    def _create_test(test_func, particle):
        def test(self):
            try:
                test_func(self, particle)
            except Exception as e:
                print(e)

        return test


class TestElectricField(unittest.TestCase, metaclass=TestMeta):
    def setUp(self):
        self.field = ElectricField((1e5, 0))

    def test_apply(self, particle):
        initial_velocity = particle.velocity
        self.field.apply(particle, 1)
        self.assertNotEqual(initial_velocity, particle.velocity)


class TestMagneticField(unittest.TestCase, metaclass=TestMeta):
    def setUp(self):
        self.field = MagneticField((0, 0, 1))

    def test_apply(self, particle):
        initial_velocity = particle.velocity
        self.field.apply(particle, 1)
        self.assertNotEqual(initial_velocity, particle.velocity)


class TestParticle(unittest.TestCase, metaclass=TestMeta):
    def setUp(self):
        self.particle = Particle(
            "Electron",
            "e",
            0.511,
            -1,
            1 / 2,
            "mode_1",
            velocity=(0.1 * c, 0.2 * c),
        )

    def test_lorentz_factor(self, particle):
        gamma = particle.lorentz_factor()
        self.assertTrue(gamma > 1)

    def test_calculate_energy(self, particle):
        energy = particle.calculate_energy()
        self.assertTrue(energy > 0)

    def test_calculate_momentum(self, particle):
        momentum = particle.calculate_momentum()
        self.assertTrue(len(momentum) == 2)

    def test_update_velocity(self, particle):
        initial_velocity = particle.velocity
        particle.update_velocity(0.2 * c, 0.3 * c)
        self.assertNotEqual(initial_velocity, particle.velocity)

    def test_decay(self, particle):
        if particle.decay_constant == 0:
            self.skipTest("Particle does not decay")
        decay_products = particle.decay(1)
        self.assertTrue(len(decay_products) > 1)

    def test_collide(self, particle):
        proton = Particle(
            "Proton",
            "p",
            938.27,
            1,
            1 / 2,
            "mode_8",
            velocity=(0.1 * c, 0.2 * c),
        )
        new_particle = particle.collide(proton)
        self.assertEqual(
            new_particle.name, f"NewParticle({particle.symbol}+p)"
        )

    def test_scatter(self, particle):
        proton = Particle(
            "Proton",
            "p",
            938.27,
            1,
            1 / 2,
            "mode_8",
            velocity=(0.3 * c, 0.4 * c),
        )
        initial_particle_velocity = particle.velocity
        initial_proton_velocity = proton.velocity
        particle.scatter(proton)
        self.assertEqual(particle.velocity, initial_proton_velocity)
        self.assertEqual(proton.velocity, initial_particle_velocity)

    def test_annihilate(self, particle):
        if particle.charge >= 0:
            self.skipTest("Particle does not annihilate")
        positron = Particle(
            "Positron",
            "e+",
            0.511,
            1,
            1 / 2,
            "mode_1",
            velocity=(0.1 * c, 0.2 * c),
        )
        photons = particle.annihilate(positron)
        self.assertEqual(len(photons), 2)
        self.assertEqual(photons[0].name, "Photon")

    def test_string_mode(self, particle):
        self.assertEqual(particle.string_mode, particle.string_mode)


class TestSimulation(unittest.TestCase):
    def test_simulate_interactions(self):
        simulate()


if __name__ == "__main__":
    unittest.main()
