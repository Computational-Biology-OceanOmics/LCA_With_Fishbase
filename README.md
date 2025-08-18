[![Tests](https://img.shields.io/badge/tests-43%20passed-brightgreen.svg)](.)
[![Coverage](https://img.shields.io/badge/coverage-73%25-yellow.svg)](.)
[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

```
Fishbase  ──>  WoRMS  ──>  NCBI  ──>  LCA
   🐟           🌊         🧬
Multi-Database Taxonomic Classification
```

This tool calculates Lowest Common Ancestor (LCA) assignments for BLAST results using multiple taxonomy databases in order of preference: Fishbase → WoRMS → NCBI Taxonomy.

Most LCA scripts rely solely on NCBI's Taxonomy database. However, fish taxonomy changes frequently, and NCBI often contains outdated classifications. Fishbase and WoRMS are updated more regularly and provide more accurate lineage information for marine organisms. Plus, Fishbase follows mostly Betancur-R, and what more could one ask for?

## Quick Start

### Installation

```bash
pip install pandas pyarrow fastparquet
```

### Basic Usage

```bash
python calculateLCAWithFishbase.py -f blast_results.tsv -o lca_results.tsv --pident 97
```

### FAIRe Compatible Usage

```bash
python calculateLCAWithFishbase_FAIReCompatible.py -f blast_results.tsv -o lca_results.tsv \
    --pident 90 --worms_file worms_species.txt.gz --asv_table asv_count_table.tsv \
    --raw_output taxaRaw.tsv --final_output taxaFinal.tsv
```

## How It Works

The tool processes BLAST results through a three-step approach:

1. **Fishbase first**: Queries Fishbase taxonomy for each BLAST hit
2. **WoRMS fallback**: For hits not found in Fishbase, queries the World Register of Marine Species (WoRMS) database
3. **NCBI backup**: Only uses NCBI Taxonomy when species aren't found in either Fishbase or WoRMS

When a query hits several species in different databases, a mix of the above may be used.

### Species Identification Method

The script:
- Goes through tabular BLAST results and checks if every word is a valid genus in: 1) Fishbase, 2) Fishbase synonyms, 3) WoRMS
- Uses the genus name to query Fishbase for taxonomy. If not in Fishbase, queries WoRMS
- If there's a genus name, assumes the next element is the species name
- If neither Fishbase nor WoRMS are found, uses the NCBI taxonomy ID from the third column
- Writes unmatched species to a missing species CSV file

### LCA Calculation Method

The LCA calculation works similarly to [eDNAFlow](https://github.com/mahsa-mousavi/eDNAFlow)'s approach:

1. Given a group of potential species for an ASV, take the species with the highest percentage base-pair identity, subtract 1 from the identity, and include all species above that cutoff in the LCA
2. **Coverage adjustment**: BP identity is adjusted by query coverage (e.g., 99% coverage × 100% identity = 99% adjusted identity)
3. **LCA grouping**: If multiple species fall within the cutoff, the species is set to 'dropped' and the algorithm moves up one taxonomic level, repeating for genus, family, order, and class

## Input Requirements

### BLAST Output Format

The input must be BLAST tabular output using this format:

```bash
-outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp"
```

### Standard Script Inputs

- **Input file**: Table of BLAST hits
- **WoRMS file** (optional): Path to WoRMS species file (default: `worms_species.txt.gz`)

### FAIRe Compatible Script Inputs

- **Input file**: Table of BLAST hits  
- **WoRMS file** (optional): Path to WoRMS species file (default: `worms_species.txt.gz`)
- **ASV table**: First column 'ASV', middle columns are samples, last column 'ASV_sequence'
- **taxaRaw/taxaFinal**: FAIRe standard format with detailed taxonomic and identification metadata

## Output

### Standard Output

Tab-delimited table with columns:
- ASV_name, Class, Order, Family, Genus, Species
- PercentageID, Coverage, Species_In_LCA, Source

Example:
```
ASV_name        Class      Order           Family         Genus        Species                    PercentageID  Coverage  Species_In_LCA                          Source
ASV_17067       Teleostei  Ophidiiformes   Ophidiidae     dropped      dropped                    89.60         100.00    Ventichthys biospeedoi, Bassozetus...   worms
ASV_17079       Teleostei  Ovalentaria...  Pomacentridae  Acanthochromis  Acanthochromis polyacanthus  79.03      90.00     Acanthochromis polyacanthus            worms
```

### FAIRe Compatible Output

- **Main output**: 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'OTU', 'numberOfUnq_BlastHits', '%ID', 'species_in_LCA', 'sources', plus sample columns

## Command Line Options

### Standard Script

```
usage: calculateLCAWithFishbase.py [-h] -f FILE -o OUTPUT [--cutoff CUTOFF] [--pident PIDENT] [--min_coverage MIN_COVERAGE] [--missing_out MISSING_OUT] [--worms_file WORMS_FILE]
                                   [--log_level {ERROR,WARNING,INFO,DEBUG}] [--normalise_identity]

Parses BLAST-tabular output and produces LCAs using Fishbase and WoRMS APIs

options:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Input file of BLAST results
  -o OUTPUT, --output OUTPUT
                        Output file of LCAs (tab-delimited)
  --cutoff CUTOFF       Basepair identity percentage cutoff for LCA calculation (default: 1.0)
  --pident PIDENT       Minimum percentage identity for BLAST hits (default: 90.0)
  --min_coverage MIN_COVERAGE
                        Minimum query coverage identity for BLAST hits (default: 90.0)
  --missing_out MISSING_OUT
                        File to write missing species to (default: missing.csv)
  --worms_file WORMS_FILE
                        Path to WoRMS species file (optional). Default is worms_species.txt.gz, included in the Github repository.
  --log_level {ERROR,WARNING,INFO,DEBUG}
                        Logging level (default: INFO)
  --normalise_identity
                        Disable identity normalisation by coverage (default: normalisation enabled). Otherwise bp identity is multiplied by coverage.

```

Example command:

    python calculateLCAWithFishbase.py -f blast_results.tsv -o lca_results.tsv --pident 97


### Key Parameters Explained

- **--cutoff**: Controls LCA stringency - larger values include more species in LCA calculations
- **--pident**: Minimum identity threshold for including BLAST hits (default: 90%)  
- **--min_coverage**: Minimum query coverage threshold (default: 90%)
- **--no_normalise_identity**: Use only BP identity without coverage adjustment (matches eDNAFlow behavior)

# Turning on coverage normalisation

    python calculateLCAWithFishbase.py -f input.txt -o output.txt --normalise_identity

There's an optional flag that lets you include the query coverage in the LCA calculation. In this case, it multiplies bp identity by coverage before calculating the LCA.

## Data sources

Fishbase: RopenSci hosts Parquet files of Fishbase species, families, and synonyms. The Python script accesses those directly.

Worms: I downloaded a relatively recent dump of WoRMS from GBIF at https://www.gbif.org/dataset/2d59e5db-57ad-41ff-97d6-11f5fb264527 and extracted the species names from the file `taxon.txt`. That file is included here (data/worms_species.txt). `grep -P '\tSpecies\t' taxon.txt | grep -w 'accepted' | cut -d'    ' -f 7,8,11,12,13,14,15,16,17,18 | gzip > worms_species.txt.gz`

NCBI: The script will download the most recent taxdump from NCBI.

When the script runs it will cache these files in a folder named `cache/`. Delete that folder to start with a clean slate.


# Tests

You can use pytest to run all tests:

    pytest

or, alternatively, you can run unittests and functional tests separately: 

    python -m unittest test_calculateLCAWithFishbase.py

    python -m unittest test_calculateLCAWithFishbase_FAIReCompatible.py

    pyton test_functional_calculateLCAWithFishbase.py

with appropriate test data in `test_data/`. These are initially written by Claude Code, then modified and expanded by me.

There's a Github CI integration to run these tests on `push` in the `.github/` folder.

# *IMPORTANT*

- The BLAST results need the taxonomy ID. Make sure that the taxonomy ID is included in the BLAST output and not just N/A. 
- Make sure that the BLAST results are formatted correctly; see below in the Input section.
- Adjust the percent identity (--pident), by default this script includes everything with >= 90% identity. That may be too lenient. Same for query coverage, where the default is also 90%.
- Fishbase, WoRMS, and NCBI Taxonomy change often. Write down the date you ran this tool.
- Some sequences on NCBI or other databases do not have specific taxonomic labels, such as 'Carangidae sp.'. These lead to very high-level LCAs, obviously. Consider removing them before running this script.

## FAQ

- How long does this run for?

In my tests for a whole eDNA dataset with a few thousand ASVs, about 2 minutes. Subsequent runs will be faster as it won't download Fishbase and NCBI data a second time.

- I have weird taxonomic sub-levels that I am interested in, why are they dropped?

The script only works with domain, phylum, class, order, family, genus, species. That way I can make the LCA calculation fairly lazy: instead of having to build a graph of all taxonomic levels and doing weird graph-based magic, I can just calculate the sets of unique species, unique genera, unique families, unique orders, unique classes, and check set size after filtering. If the set size is > 1, set this taxonomic level's taxonomic label to 'dropped'. Easier than breaking my head trying to come up with recursive tree-walking algorithms for the sake of methodological complexity you can publish in a paper, I'd rather have results.

- I may have unacccepted taxonomic names in my results.

Fishbase is nice in that it knows about many, but not all unaccepted names. If an unaccepted name made it into Fishbase the API returns the accepted version of this name along with the taxonomic lineage of the new, accepted name. We then use that name for the LCA. Neat, isn't it!  
The script then checks whether the species is in WoRMS - lastly, your misspelled species name may be in the NCBI taxonomy. Check your missing names csv to see if it got picked up.

- I get an error ValueError: Wrong number of dimensions. values.ndim > ndim [2 > 1]

There seems to be a strange bug with older versions of pyarrow and newer versions of pandas. Running `pip install pyarrow pandas -U` fixed it for me.

- I have checked my missing.csv results and realised that I have a species that is neither on Fishbase, WoRMS, or NCBI Taxonomy. What can I do?

The easiest solution is to manually add this species to the worms_species.txt file following the WoRMS taxonomy. Here's an example for a new fish species:

    Melanophorichthys penicillus    Melanophorichthys       Animalia        Chordata        Teleostei       Perciformes     Gobiesocidae    Melanophorichthys       penicillus

- I have more questions!

Please contact me at pbayer AT minderoo.org

## Acknowledgments

- Thanks to Dr Shannon Corrigan for pushing me to use the *correct* taxonomy.

# AI statement

I've written an initial prototype (minus the NCBI Taxonomy!!!) of this code and then gave that to Claude for hardening and making production-ready. I then fixed a bunch of Claude-introduced bugs, removed most of the unnecessary comments, and added some minor bugfixes that my original implementation had already. The original non-AI code used to live in calculateLCAWithFishbase.py, and the AI code lived in calculateLCAWithFishbase_Claude.py but that got confusing, so since version 0.04 it's only one script.
Most of the tests were written by Claude Code, as well.

## CHANGELOG

- v0.04: moved Claude code to main script, now filters by 90% query coverage and normalises bp identity by coverage before LCA calculation
- v0.03: added script for better compatibility with the FAIRe metadata standards
- v0.02: adding NCBI taxonomy
- v0.01: initial release a.k.a. 'works on my system'


## TODO

- parsing out the genus name from the lines is very iffy. I'm sure there will be cases where it's just 'H. sapiens', or '[bla]Homo sapiens'. I'll probably find many cases where I have to add more rules to cleaning the strings. I could use taxonerd but that feels a bit like overkill. Could also store all fishbase/worms/NCBI species names in a lookup table and then align with every BLAST result line.
- the tests rely on unittest. Could go 'pure' pytest, too.... not too hard?
