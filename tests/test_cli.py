#
# Copyright (C) 2018 University of Oxford
#
# This file is part of tsinfer.
#
# tsinfer is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tsinfer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with tsinfer.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Tests for the tsinfer CLI.
"""

import io
import os.path
import pathlib
import sys
import tempfile
import unittest
import unittest.mock as mock

import numpy as np
import msprime

import tsinfer
import tsinfer.cli as cli
import tsinfer.__main__ as main
import tsinfer.exceptions as exceptions


def capture_output(func, *args, **kwargs):
    """
    Runs the specified function and arguments, and returns the
    tuple (stdout, stderr) as strings.
    """
    buffer_class = io.BytesIO
    if sys.version_info[0] == 3:
        buffer_class = io.StringIO
    stdout = sys.stdout
    sys.stdout = buffer_class()
    stderr = sys.stderr
    sys.stderr = buffer_class()

    try:
        func(*args, **kwargs)
        stdout_output = sys.stdout.getvalue()
        stderr_output = sys.stderr.getvalue()
    finally:
        sys.stdout.close()
        sys.stdout = stdout
        sys.stderr.close()
        sys.stderr = stderr
    return stdout_output, stderr_output


class TestMain(unittest.TestCase):
    """
    Simple tests for the main function.
    """
    def test_cli_main(self):
        with mock.patch("argparse.ArgumentParser.parse_args") as mocked_parse:
            cli.tsinfer_main()
            mocked_parse.assert_called_once_with(None)

    def test_main(self):
        with mock.patch("argparse.ArgumentParser.parse_args") as mocked_parse:
            main.main()
            mocked_parse.assert_called_once_with(None)


class TestDefaultPaths(unittest.TestCase):
    """
    Tests for the default path creation routines.
    """
    def test_get_default_path(self):
        # The second argument is ignored if the input path is specified.
        for path in ["a", "a/b/c", "a.stuff"]:
            self.assertEqual(cli.get_default_path(path, "a", "b"), path)
        self.assertEqual(cli.get_default_path(None, "a", ".x"), "a.x")
        self.assertEqual(cli.get_default_path(None, "a.y", ".z"), "a.z")
        self.assertEqual(cli.get_default_path(None, "a/b/c/a.y", ".z"), "a/b/c/a.z")

    def test_get_ancestors_path(self):
        self.assertEqual(cli.get_ancestors_path(None, "a"), "a.ancestors")

    def test_get_ancestors_trees_path(self):
        self.assertEqual(cli.get_ancestors_trees_path(None, "a"), "a.ancestors.trees")

    def test_get_output_trees_path(self):
        self.assertEqual(cli.get_output_trees_path(None, "a"), "a.trees")


class TestCli(unittest.TestCase):
    """
    Parent of all CLI test cases.
    """
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory(prefix="tsinfer_cli_test")
        self.sample_file = str(pathlib.Path(self.tempdir.name, "input-data.samples"))
        self.ancestor_file = str(pathlib.Path(self.tempdir.name, "input-data.ancestors"))
        self.input_ts = msprime.simulate(
            10, mutation_rate=10, recombination_rate=10, random_seed=10)
        sample_data = tsinfer.SampleData(
            sequence_length=self.input_ts.sequence_length, path=self.sample_file)
        for var in self.input_ts.variants():
            sample_data.add_site(var.site.position, var.genotypes, var.alleles)
        sample_data.finalise()
        tsinfer.generate_ancestors(sample_data, path=self.ancestor_file, chunk_size=10)
        sample_data.close()


class TestCommandsDefaults(TestCli):
    """
    Tests that the basic commands work if we provide the default arguments.
    """
    def verify_output(self, output_path):
        output_trees = msprime.load(output_path)
        self.assertEqual(output_trees.num_samples, self.input_ts.num_samples)
        self.assertEqual(output_trees.sequence_length, self.input_ts.sequence_length)
        self.assertEqual(output_trees.num_sites, self.input_ts.num_sites)
        self.assertGreater(output_trees.num_sites, 1)
        self.assertTrue(np.array_equal(
            output_trees.genotype_matrix(), self.input_ts.genotype_matrix()))

    # Need to mock out setup_logging here or we spew logging to the console
    # in later tests.
    @mock.patch("tsinfer.cli.setup_logging")
    def run_command(self, command, mock_setup_logging):
        stdout, stderr = capture_output(cli.tsinfer_main, command)
        self.assertEqual(stderr, "")
        self.assertEqual(stdout, "")
        self.assertTrue(mock_setup_logging.called)

    def test_infer(self):
        output_trees = os.path.join(self.tempdir.name, "output.trees")
        self.run_command(["infer", self.sample_file, "-O", output_trees])
        self.verify_output(output_trees)

    def test_nominal_chain(self):
        output_trees = os.path.join(self.tempdir.name, "output.trees")
        self.run_command(["generate-ancestors", self.sample_file])
        self.run_command(["match-ancestors", self.sample_file])
        self.run_command(["match-samples", self.sample_file, "-O", output_trees])
        self.verify_output(output_trees)


class TestList(TestCli):
    """
    Tests cases for the list command.
    """
    # Need to mock out setup_logging here or we spew logging to the console
    # in later tests.
    @mock.patch("tsinfer.cli.setup_logging")
    def run_command(self, command, mock_setup_logging):
        stdout, stderr = capture_output(cli.tsinfer_main, command)
        self.assertEqual(stderr, "")
        self.assertTrue(mock_setup_logging.called)
        return stdout

    def test_list_samples(self):
        output1 = self.run_command(["list", self.sample_file])
        self.assertGreater(len(output1), 0)
        output2 = self.run_command(["ls", self.sample_file])
        self.assertEqual(output1, output2)

    def test_list_samples_storage(self):
        output1 = self.run_command(["list", "-s", self.sample_file])
        self.assertGreater(len(output1), 0)
        output2 = self.run_command(["list", "--storage", self.sample_file])
        self.assertEqual(output1, output2)

    def test_list_ancestors(self):
        output1 = self.run_command(["list", self.ancestor_file])
        self.assertGreater(len(output1), 0)
        output2 = self.run_command(["ls", self.ancestor_file])
        self.assertEqual(output1, output2)

    def test_list_ancestors_storage(self):
        output1 = self.run_command(["list", "-s", self.ancestor_file])
        self.assertGreater(len(output1), 0)
        output2 = self.run_command(["list", "--storage", self.ancestor_file])
        self.assertEqual(output1, output2)

    def test_list_unknown_files(self):
        zero_file = os.path.join(self.tempdir.name, "zeros")
        with open(zero_file, "wb") as f:
            f.write(bytearray(100))
        for bad_file in ["/", zero_file]:
            self.assertRaises(
                exceptions.FileFormatError, self.run_command, ["list", bad_file])
