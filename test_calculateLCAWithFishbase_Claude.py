#!/usr/bin/env python3
"""
Automated tests for calculateLCAWithFishbase_Claude.py

This module contains comprehensive tests for the BLAST LCA Analysis Tool,
covering all major components and edge cases.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import pandas as pd
from collections import OrderedDict

from calculateLCAWithFishbase_Claude import (
    Config, TaxonomicLineage, LCAResult, NCBITaxdumpParser,
    DatabaseManager, SpeciesNameCorrector, TaxonomicAssigner,
    LCACalculator, BLASTLCAAnalyzer
)


class TestConfig(unittest.TestCase):
    """Test Configuration class."""
    
    def test_config_defaults(self):
        """Test that config has expected default values."""
        config = Config()
        self.assertEqual(config.DEFAULT_CUTOFF, 1.0)
        self.assertEqual(config.DEFAULT_PIDENT_CUTOFF, 90.0)
        self.assertEqual(len(config.BLAST_COLUMNS), 22)
        self.assertEqual(config.PIDENT_COLUMN_INDEX, 6)
        self.assertIn("ftp.ncbi.nlm.nih.gov", config.NCBI_TAXDUMP_URL)


class TestTaxonomicLineage(unittest.TestCase):
    """Test TaxonomicLineage dataclass."""
    
    def test_taxonomic_lineage_creation(self):
        """Test creating and converting taxonomic lineage."""
        lineage = TaxonomicLineage(
            class_name="Actinopterygii",
            order="Perciformes", 
            family="Gobiidae",
            genus="Gobius",
            species="Gobius niger"
        )
        
        lineage_list = lineage.to_list()
        expected = [
            ("C", "Actinopterygii"),
            ("O", "Perciformes"),
            ("F", "Gobiidae"),
            ("G", "Gobius"),
            ("S", "Gobius niger")
        ]
        self.assertEqual(lineage_list, expected)


class TestLCAResult(unittest.TestCase):
    """Test LCAResult dataclass."""
    
    def test_lca_result_creation(self):
        """Test LCA result creation with all fields."""
        taxa = {"Gobius niger", "Gobius cobitis"}
        result = LCAResult(
            percentage=95.5,
            assignment="Gobius",
            included_taxa=taxa
        )
        
        self.assertEqual(result.percentage, 95.5)
        self.assertEqual(result.assignment, "Gobius")
        self.assertEqual(result.included_taxa, taxa)


class TestNCBITaxdumpParser(unittest.TestCase):
    """Test NCBI Taxdump Parser."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = Path(tempfile.mkdtemp())
        self.parser = NCBITaxdumpParser(self.temp_dir)
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir)
    
    def test_parser_initialization(self):
        """Test parser initialization."""
        self.assertEqual(self.parser.cache_dir, self.temp_dir)
        self.assertTrue(self.temp_dir.exists())
        self.assertEqual(len(self.parser.taxid_to_lineage), 0)
    
    @patch('calculateLCAWithFishbase_Claude.urlretrieve')
    @patch('calculateLCAWithFishbase_Claude.tarfile.open')
    def test_download_and_extract_taxdump(self, mock_tarfile, mock_urlretrieve):
        """Test downloading and extracting taxdump."""
        mock_tar = MagicMock()
        mock_tarfile.return_value.__enter__.return_value = mock_tar
        
        result = self.parser.download_and_extract_taxdump()
        
        expected_dir = self.temp_dir / "taxdump"
        self.assertEqual(result, expected_dir)
        mock_urlretrieve.assert_called_once()
        mock_tar.extractall.assert_called_once_with(expected_dir)
    
    def test_parse_nodes_file(self):
        """Test parsing nodes.dmp file."""
        nodes_content = "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tsuperkingdom\t|\n"
        nodes_file = self.temp_dir / "nodes.dmp"
        nodes_file.write_text(nodes_content)
        
        self.parser.parse_nodes_file(nodes_file)
        
        self.assertEqual(self.parser.taxid_to_parent["1"], "1")
        self.assertEqual(self.parser.taxid_to_parent["2"], "1")
        self.assertEqual(self.parser.taxid_to_rank["1"], "no rank\t|")
        self.assertEqual(self.parser.taxid_to_rank["2"], "superkingdom\t|")


class TestSpeciesNameCorrector(unittest.TestCase):
    """Test Species Name Corrector."""
    
    def test_default_corrections(self):
        """Test default corrections are applied."""
        corrector = SpeciesNameCorrector()
        line = "ASV1\tseq123\t12345\tPetroschmidtia albonotatus\tcommon\tEukaryota\t95.0"
        corrected = corrector.correct_line(line)
        self.assertIn("Petroschmidtia albonotata", corrected)
        self.assertNotIn("Petroschmidtia albonotatus", corrected)
    
    def test_no_corrections_needed(self):
        """Test line with no corrections needed."""
        corrector = SpeciesNameCorrector()
        line = "ASV1\tseq123\t12345\tGobius niger\tcommon\tEukaryota\t95.0"
        corrected = corrector.correct_line(line)
        self.assertEqual(line, corrected)


class TestLCACalculator(unittest.TestCase):
    """Test LCA Calculator."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.calculator = LCACalculator(cutoff=1.0)
    
    def test_single_taxon_lca(self):
        """Test LCA calculation with single taxon."""
        entries = [(95.0, "Gobius niger")]
        result = self.calculator.calculate_lca(entries)
        
        self.assertEqual(result.percentage, 95.0)
        self.assertEqual(result.assignment, "Gobius niger")
        self.assertEqual(result.included_taxa, {"Gobius niger"})
    
    def test_multiple_taxa_within_cutoff(self):
        """Test LCA calculation with multiple taxa within cutoff."""
        entries = [(95.0, "Gobius niger"), (94.5, "Gobius cobitis")]
        result = self.calculator.calculate_lca(entries)
        
        self.assertEqual(result.percentage, 95.0)  # CHANGE: should be the max now. 
        self.assertEqual(result.assignment, "dropped")
        self.assertEqual(result.included_taxa, {"Gobius niger", "Gobius cobitis"})
    
    def test_empty_entries(self):
        """Test LCA calculation with empty entries."""
        entries = []
        result = self.calculator.calculate_lca(entries)
        
        self.assertEqual(result.percentage, 0.0)
        self.assertEqual(result.assignment, "no_hits")
        self.assertEqual(result.included_taxa, set())


class TestTaxonomicAssigner(unittest.TestCase):
    """Test Taxonomic Assigner."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Mock Fishbase data
        self.fishbase_genera = {
            "Gobius": ["Gobiidae", "Perciformes", "Actinopterygii"]
        }
        self.fishbase_speccode = {
            123: "Gobius niger"
        }
        self.fishbase_synonyms = {
            "Gobius old_name": 123
        }
        
        # Mock WoRMS data
        self.worms_genera = {
            "Pagrus": ["Sparidae", "Perciformes", "Actinopterygii"]
        }
        self.worms_species = {"Pagrus pagrus"}
        
        # Mock database manager
        self.db_manager = Mock()
        
        self.assigner = TaxonomicAssigner(
            self.fishbase_genera, self.fishbase_speccode, self.fishbase_synonyms,
            self.worms_genera, self.worms_species, self.db_manager
        )
    
    def test_fishbase_genus_match(self):
        """Test finding species in Fishbase by genus."""
        elements = ["some", "prefix", "Gobius", "niger", "suffix"]
        result = self.assigner.find_species_info(elements)
        
        self.assertIsNotNone(result)
        genus, species, source, lineage = result
        self.assertEqual(genus, "Gobius")
        self.assertEqual(species, "niger")
        self.assertEqual(source, "fishbase")
        self.assertEqual(lineage.genus, "Gobius")
    
    def test_worms_genus_match(self):
        """Test finding species in WoRMS when not in Fishbase."""
        elements = ["some", "prefix", "Pagrus", "pagrus", "suffix"]
        result = self.assigner.find_species_info(elements)
        
        self.assertIsNotNone(result)
        genus, species, source, lineage = result
        self.assertEqual(genus, "Pagrus")
        self.assertEqual(species, "pagrus")
        self.assertEqual(source, "worms")
        self.assertEqual(lineage.genus, "Pagrus")
    
    def test_ncbi_fallback(self):
        """Test falling back to NCBI when not found in other databases."""
        elements = ["unknown", "genus", "species"]
        mock_lineage = TaxonomicLineage("TestClass", "TestOrder", "TestFamily", "TestGenus", "TestSpecies")
        self.db_manager.query_ncbi_taxonomy.return_value = mock_lineage
        
        result = self.assigner.find_species_info(elements, taxid="12345")
        
        self.assertIsNotNone(result)
        genus, species, source, lineage = result
        self.assertEqual(source, "ncbi")
        self.assertEqual(lineage, mock_lineage)


class TestBLASTLCAAnalyzer(unittest.TestCase):
    """Test main BLAST LCA Analyzer."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.config = Config()
        self.analyzer = BLASTLCAAnalyzer(self.config)
        self.temp_dir = Path(tempfile.mkdtemp())
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir)
    
    def test_analyzer_initialization(self):
        """Test analyzer initialization."""
        self.assertEqual(self.analyzer.config, self.config)
        self.assertIsNotNone(self.analyzer.db_manager)
        self.assertIsNotNone(self.analyzer.corrector)
        self.assertIsNotNone(self.analyzer.lca_calculator)
    
    @patch.object(BLASTLCAAnalyzer, 'load_databases')
    def test_process_blast_file_valid_input(self, mock_load_db):
        """Test processing valid BLAST file."""
        # Create test BLAST file
        blast_content = "ASV1\tseq1\t12345\tGobius niger\tgoby\tEukaryota\t95.0\t100\t150\t200\t2\t0\t0\t1\t100\t1\t200\ttitle\t1e-50\t180\t90\t85\n"
        blast_file = self.temp_dir / "test_blast.txt"
        blast_file.write_text(blast_content)
        
        missing_file = self.temp_dir / "missing.txt"
        
        # Mock the assigner
        mock_lineage = TaxonomicLineage("TestClass", "TestOrder", "TestFamily", "Gobius", "Gobius niger")
        self.analyzer.assigner = Mock()
        self.analyzer.assigner.find_species_info.return_value = ("Gobius", "niger", "fishbase", mock_lineage)
        
        result = self.analyzer.process_blast_file(blast_file, 90.0, missing_file)
        
        self.assertIn("ASV1", result)
        self.assertEqual(len(result["ASV1"]), 1)
        source, pident, lineage = result["ASV1"][0]
        self.assertEqual(source, "fishbase")
        self.assertEqual(pident, 95.0)
    
    def test_calculate_lca_assignments(self):
        """Test LCA assignment calculation."""
        # Create mock hits data
        mock_lineage = TaxonomicLineage("Actinopterygii", "Perciformes", "Gobiidae", "Gobius", "Gobius niger")
        asv_hits = {
            "ASV1": [("fishbase", 95.0, mock_lineage)]
        }
        
        results = self.analyzer.calculate_lca_assignments(asv_hits)
        
        self.assertEqual(len(results), 1)
        result = results[0]
        self.assertEqual(result["ASV_name"], "ASV1")
        self.assertEqual(result["Class"], "Actinopterygii")
        self.assertEqual(result["Order"], "Perciformes")
        self.assertEqual(result["Family"], "Gobiidae")
        self.assertEqual(result["Genus"], "Gobius")
        self.assertEqual(result["Species"], "Gobius niger")
    
    def test_write_results(self):
        """Test writing results to file."""
        results = [{
            'ASV_name': 'ASV1',
            'Class': 'Actinopterygii',
            'Order': 'Perciformes',
            'Family': 'Gobiidae',
            'Genus': 'Gobius',
            'Species': 'Gobius niger',
            'PercentageID': '95.00',
            'Species_In_LCA': 'Gobius niger',
            'Sources': 'fishbase'
        }]
        
        output_file = self.temp_dir / "output.txt"
        self.analyzer.write_results(results, output_file)
        
        self.assertTrue(output_file.exists())
        content = output_file.read_text()
        lines = content.strip().split('\n')
        self.assertEqual(len(lines), 2)  # header + 1 data line
        self.assertIn("ASV1", content)
        self.assertIn("Gobius niger", content)


if __name__ == '__main__':
    unittest.main()
