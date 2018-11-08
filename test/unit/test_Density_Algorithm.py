import unittest
from Density_Algorithm import left_outside_from_window
# NOTE consider a different structure where you only need to import one or maybe two classes || functions

class TestDensityAlgorithm(unittest.TestCase):

    def setUp(self):
        """Called prior to each test."""
        pass  # create member variables used in each test
    # fed

    def tearDown(self):
        """Called following each test."""
        pass
    # fed

    def test_left_outside_from_window(self):

        # FIXME are these values right? are they the right type? etc.
        window = 500
        left_window_stop = 600
        left_window_start = 400
        left_density = 42
        expected_density = 0.50
        density = left_outside_from_window(
            window, left_window_stop, left_window_start, left_density)
        self.assertAlmostEqual(density, expected_density)
    # fed
# ssalc


if __name__ == '__main__':

    unittest.main()
