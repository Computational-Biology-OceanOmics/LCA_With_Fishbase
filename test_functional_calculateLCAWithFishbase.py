#!/usr/bin/env python3
"""
Functional tests for calculateLCAWithFishbase.py

This module contains end-to-end functional tests that run the complete
calculateLCAWithFishbase.py script with real input data and verify
the output against expected results.
"""

import unittest
import subprocess
import tempfile
import shutil
import sys
from pathlib import Path
import pandas as pd


class TestCalculateLCAFunctional(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(__file__).parent / "test_data/"
        self.script_path = Path(__file__).parent / "calculateLCAWithFishbase.py"
        self.test_input = self.test_dir / "test_data.tsv"
        self.expected_output = self.test_dir / "test_data.asv.tsv"

        self.temp_dir = Path(tempfile.mkdtemp())
        self.test_output = self.temp_dir / "test_output.tsv"
        self.test_missing = self.temp_dir / "test_missing.csv"

        self.assertTrue(
            self.script_path.exists(), f"Script not found: {self.script_path}"
        )
        self.assertTrue(
            self.test_input.exists(), f"Test input not found: {self.test_input}"
        )
        self.assertTrue(
            self.expected_output.exists(),
            f"Expected output not found: {self.expected_output}",
        )

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir)

    def test_end_to_end_lca_calculation(self):
        """Test the complete LCA calculation pipeline end-to-end."""
        cmd = [
            sys.executable,
            str(self.script_path),
            "-f",
            str(self.test_input),
            "-o",
            str(self.test_output),
            "--missing_out",
            str(self.test_missing),
            "--log_level",
            "WARNING",  # Reduce log noise in tests
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.test_dir)

        # did we exit successfully? would be weird if not
        self.assertEqual(
            result.returncode,
            0,
            f"Script failed with return code {result.returncode}. "
            f"STDOUT: {result.stdout}. STDERR: {result.stderr}",
        )

        self.assertTrue(self.test_output.exists(), "Output file was not created")

        expected_df = pd.read_csv(self.expected_output, sep="\t")
        actual_df = pd.read_csv(self.test_output, sep="\t")

        # check that we have the same number of rows
        self.assertEqual(
            len(actual_df),
            len(expected_df),
            f"Different number of rows: expected {len(expected_df)}, got {len(actual_df)}",
        )

        # check that all expected ASV names are present
        expected_asvs = set(expected_df["ASV_name"])
        actual_asvs = set(actual_df["ASV_name"])
        self.assertEqual(
            expected_asvs,
            actual_asvs,
            f"Missing ASVs: {expected_asvs - actual_asvs}. "
            f"Extra ASVs: {actual_asvs - expected_asvs}",
        )

        # check column structure
        expected_cols = set(expected_df.columns)
        actual_cols = set(actual_df.columns)
        self.assertEqual(
            expected_cols,
            actual_cols,
            f"Column mismatch: expected {expected_cols}, got {actual_cols}",
        )

        # sort both dataframes by ASV_name for consistent comparison
        expected_df = expected_df.sort_values("ASV_name").reset_index(drop=True)
        actual_df = actual_df.sort_values("ASV_name").reset_index(drop=True)

        # compare each row
        for i, (expected_row, actual_row) in enumerate(
            zip(expected_df.itertuples(), actual_df.itertuples())
        ):
            with self.subTest(row=i, asv=expected_row.ASV_name):
                self.assertEqual(actual_row.ASV_name, expected_row.ASV_name)
                self.assertEqual(actual_row.Class, expected_row.Class)
                self.assertEqual(actual_row.Order, expected_row.Order)
                self.assertEqual(actual_row.Family, expected_row.Family)
                self.assertEqual(actual_row.Genus, expected_row.Genus)
                self.assertEqual(actual_row.Species, expected_row.Species)

                self.assertAlmostEqual(
                    float(actual_row.PercentageID),
                    float(expected_row.PercentageID),
                    places=2,
                )
                self.assertAlmostEqual(
                    float(actual_row.Coverage), float(expected_row.Coverage), places=2
                )

                # species_In_LCA might have different ordering but same content
                expected_species = (
                    set(expected_row.Species_In_LCA.split(", "))
                    if expected_row.Species_In_LCA
                    else set()
                )
                actual_species = (
                    set(actual_row.Species_In_LCA.split(", "))
                    if actual_row.Species_In_LCA
                    else set()
                )
                self.assertEqual(actual_species, expected_species)

                # sources might have different ordering but same content
                expected_sources = (
                    set(expected_row.Sources.split(", "))
                    if expected_row.Sources
                    else set()
                )
                actual_sources = (
                    set(actual_row.Sources.split(", ")) if actual_row.Sources else set()
                )
                self.assertEqual(actual_sources, expected_sources)

    def test_script_with_different_cutoffs(self):
        """Test script with different percentage and coverage cutoffs."""
        cmd = [
            sys.executable,
            str(self.script_path),
            "-f",
            str(self.test_input),
            "-o",
            str(self.test_output),
            "--cutoff",
            "2.0",
            "--pident",
            "85.0",
            "--min_coverage",
            "85.0",
            "--missing_out",
            str(self.test_missing),
            "--log_level",
            "WARNING",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.test_dir)

        self.assertEqual(
            result.returncode,
            0,
            f"Script failed with return code {result.returncode}. "
            f"STDERR: {result.stderr}",
        )

        self.assertTrue(self.test_output.exists())
        actual_df = pd.read_csv(self.test_output, sep="\t")
        self.assertGreater(len(actual_df), 0, "No results with different cutoffs")

    def test_missing_output_creation(self):
        """Test that missing species file is created properly."""
        cmd = [
            sys.executable,
            str(self.script_path),
            "-f",
            str(self.test_input),
            "-o",
            str(self.test_output),
            "--missing_out",
            str(self.test_missing),
            "--log_level",
            "WARNING",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.test_dir)

        self.assertEqual(result.returncode, 0)

        # missing file should exist (even if empty)
        self.assertTrue(self.test_missing.exists(), "Missing file was not created")

    def test_invalid_input_file(self):
        """Test script behavior with non-existent input file."""
        non_existent_file = self.temp_dir / "non_existent.tsv"

        cmd = [
            sys.executable,
            str(self.script_path),
            "-f",
            str(non_existent_file),
            "-o",
            str(self.test_output),
            "--log_level",
            "WARNING",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        self.assertNotEqual(
            result.returncode, 0, "Script should fail with non-existent input file"
        )

    def test_invalid_parameters(self):
        """Test script behavior with invalid parameters."""
        cmd = [
            sys.executable,
            str(self.script_path),
            "-f",
            str(self.test_input),
            "-o",
            str(self.test_output),
            "--pident",
            "150.0",  # Invalid: > 100
            "--log_level",
            "WARNING",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        self.assertNotEqual(
            result.returncode, 0, "Script should fail with invalid pident > 100"
        )

        cmd = [
            sys.executable,
            str(self.script_path),
            "-f",
            str(self.test_input),
            "-o",
            str(self.test_output),
            "--cutoff",
            "-1.0",  # Invalid: negative
            "--log_level",
            "WARNING",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        self.assertNotEqual(
            result.returncode, 0, "Script should fail with negative cutoff"
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
