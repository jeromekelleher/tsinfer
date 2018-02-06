"""
Test cases for the evaluation code.
"""

import unittest

import msprime
import numpy as np

import tsinfer


class TestKCMetric(unittest.TestCase):
    """
    Tests on the KC metric distances.
    """
    def test_same_tree_zero_distance(self):
        for n in range(2, 10):
            for seed in range(1, 10):
                ts = msprime.simulate(n, random_seed=seed)
                tree = ts.first()
                self.assertEqual(tsinfer.kc_distance(tree, tree), 0)
                ts = msprime.simulate(n, random_seed=seed)
                tree2 = ts.first()
                self.assertEqual(tsinfer.kc_distance(tree, tree2), 0)

    def test_sample_2_zero_distance(self):
        # All trees with 2 leaves must be equal distance from each other.
        for seed in range(1, 10):
            tree1 = msprime.simulate(2, random_seed=seed).first()
            tree2 = msprime.simulate(2, random_seed=seed + 1).first()
            self.assertEqual(tsinfer.kc_distance(tree1, tree2), 0)

    def test_different_samples_error(self):
        tree1 = msprime.simulate(10, random_seed=1).first()
        tree2 = msprime.simulate(2, random_seed=1).first()
        self.assertRaises(ValueError, tsinfer.kc_distance, tree1, tree2)

    # TODO add more tests checking actual examples.


class TestTreeSequenceCompare(unittest.TestCase):
    """
    Tests of the method to compare to tree sequences.
    """
    def test_same_ts(self):
        n = 15
        for seed in range(1, 10):
            ts = msprime.simulate(n, recombination_rate=10, random_seed=seed)
            self.assertGreater(ts.num_trees, 1)
            bp, distance = tsinfer.compare(ts, ts)
            self.assertEqual(list(bp), list(ts.breakpoints()))
            self.assertEqual(distance.shape, (bp.shape[0] - 1,))
            self.assertTrue(np.all(distance == 0))

    def test_single_tree(self):
        n = 15
        for seed in range(1, 10):
            ts1 = msprime.simulate(n, random_seed=seed)
            ts2 = msprime.simulate(n, random_seed=seed + 1)
            bp, distance = tsinfer.compare(ts1, ts2)
            self.assertEqual(list(bp), [0, 1])
            self.assertEqual(distance.shape, (1,))

    def test_single_tree_many_trees(self):
        n = 5
        for seed in range(1, 10):
            ts1 = msprime.simulate(n, recombination_rate=5, random_seed=seed)
            ts2 = msprime.simulate(n, random_seed=seed + 1)
            self.assertGreater(ts1.num_trees, 1)
            bp, distance = tsinfer.compare(ts1, ts2)
            self.assertEqual(list(bp), list(ts1.breakpoints()))
            self.assertEqual(distance.shape, (ts1.num_trees,))

    def test_single_many_trees(self):
        n = 5
        for seed in range(1, 10):
            ts1 = msprime.simulate(n, recombination_rate=5, random_seed=seed)
            ts2 = msprime.simulate(n, recombination_rate=5, random_seed=seed + 1)
            self.assertGreater(ts1.num_trees, 1)
            self.assertGreater(ts2.num_trees, 1)
            bp, distance = tsinfer.compare(ts1, ts2)
            breakpoints = set(ts1.breakpoints()) | set(ts2.breakpoints())
            self.assertEqual(list(bp), sorted(breakpoints))
            self.assertEqual(distance.shape, (len(breakpoints) - 1,))
    # TODO add some examples testing for specific instances.
