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

## 📋 Table of Contents
- [Quick Start](#usage)
- [Dependencies](#dependencies)
- [Input/Output](#inputoutput)
- [Method](#method-species-id)
- [Data Sources](#data-sources)
- [FAQ](#faq)
- [Tests](#tests)
- [Changelog](#changelog)

The tool processes BLAST results through a three-step approach:

**Fishbase first**: Queries Fishbase taxonomy for each BLAST hit. If all hits for a query sequence are found, calculates LCA using Fishbase lineages.  
**WoRMS fallback**: For hits not found in Fishbase, queries the World Register of Marine Species (WoRMS) database.  
**NCBI backup**: Only uses NCBI Taxonomy when species aren't found in either Fishbase or WoRMS.  

When a query hits several species in different databases, a mix of the above may be used.

# Input/Output

**Input**: a table of BLAST hits (see below for format)

**Worms file**: Path to WoRMS species file (optional). Default is worms_species.txt.gz, included in the Github repository.

**Output**: a table of LCAs for every query in the BLAST table.

# Input/Output when using calculateLCAWithFishbase_FAIReCompatible.py

This requires more input data (the asv_count_table.tsv and the final and raw taxa tables) but is more useful if you want output in the FAIRe standard.

**Input**: a table of BLAST hits (see below for format)

**Worms file**: Path to WoRMS species file (optional). Default is worms_species.txt.gz, included in the Github repository.

**ASV_table**: First column should be 'ASV', the next columns should be the samples, the last column should be 'ASV_sequence'

**Output**: a table of LCAs for every query in the BLAST table.

**Raw Output**: A taxaRaw table following the FAIRe standards.

**Final Output**: A taxaFinal table following the FAIRe standards.

# Usage

    python calculateLCAWithFishbase.py -f blast_results.tsv -o lca_results.tsv --pident 97

# Usage for FAIRe compatibility

    python calculateLCAWithFishbase_FAIReCompatible.py -f blast_results.tsv -o lca_results.tsv --pident 90 --worms_file worms_species.txt.gz --asv_table asv_count_table.tsv --raw_output taxaRaw.tsv --final_output taxaFinal.tsv

# Dependencies

    pip install pandas pyarrow fastparquet

Any fairly recent version (2023-2025) should be fine, I believe.

# Turning off coverage normalisation

    python calculateLCAWithFishbase.py -f input.txt -o output.txt --no_normalise_identity

There's an optional flag that lets you ignore the query coverage in the LCA calculation. In this case, it uses only bp identity to calculate the LCA.

# AI statement

I've written an initial prototype (minus the NCBI Taxonomy!!!) of this code and then gave that to Claude for hardening and making production-ready. I then fixed a bunch of Claude-introduced bugs, removed most of the unnecessary comments, and added some minor bugfixes that my original implementation had already. The original non-AI code used to live in calculateLCAWithFishbase.py, and the AI code lived in calculateLCAWithFishbase_Claude.py but that got confusing, so since version 0.04 it's only one script.

# Tests

You can use pytest to run all tests:

    pytest

or, alternatively, you can run unittests and functional tests separately: 

    python -m unittest test_calculateLCAWithFishbase.py

    python -m unittest test_calculateLCAWithFishbase_FAIReCompatible.py

    pyton test_functional_calculateLCAWithFishbase.py

with appropriate test data in `test_data/`. These are initially written by Claude Code, then modified and expanded by me.

# *IMPORTANT*

- The BLAST results need the taxonomy ID. Make sure that the taxonomy ID is included in the BLAST output and not just N/A. 
- Make sure that the BLAST results are formatted correctly; see below in the Input section.
- Adjust the percent identity (--pident), by default this script includes everything with >= 90% identity. That may be too lenient. Same for query coverage, where the default is also 90%.
- Fishbase, WoRMS, and NCBI Taxonomy change often. Write down the date you ran this tool.
- Some sequences on NCBI or other databases do not have specific taxonomic labels, such as 'Carangidae sp.'. These lead to very high-level LCAs, obviously. Consider removing them before running this script.

## Method: Species ID:

What the script does:
- go through the tabular blast results, check if every word is a valid genus in 1. Fishbase, 2. Fishbase synonyms, 3. Worms. We want to trust the Fishbase taxonomy the most but not every species we hit is in Fishbase. Mammals etc. will instead hit into Worms.
- using the genus name, ask Fishbase what the taxonomy for that genus is. If the genus is not in Fishbase, ask Worms.
- if both Fishbase and WoRMS were not found in the row, assume that the third column is the NCBI taxonomy ID. Use that to look up the lineage instead.
- If neither Fishbase nor WoRMS nor NCBI have the species or genus, write the entire line of BLAST results to the missing species CSV.

## Method: LCA calculation

The LCA calculation works almost the same way as [eDNAFlow](https://github.com/mahsa-mousavi/eDNAFlow)'s LCA calculation, except that we also normalise by query coverage.

1) Given a group of potential species for an ASV, take the species with the highest percentage base-pair identity, subtract 1 from the identity, and then include all species above that cutoff in the LCA.
2) There is one change in the way the eDNAFlow script works: we also include query coverage. In the past we used only a 100% query coverage so this was moot, but we found many ASVs that overlapped 98%, 99% etc., not 100%. So now we adjust the bp identity by the query coverage, and then use that adjusted bp identity in the LCA calculation. The adjustment is multiplication: a query coverage of 99% and a bp identity of 100% means that the adjusted bp identity becomes 0.99 * 1 = 0.99 = 99%.
3) The LCA itself is just a grouping: if there are several species within the cutoff, then the species is set to 'dropped' and we go up one taxonomic level, repeat for the genus, repeat for the family, repeat for the class, repeat for the order. There's no LCA voting or similar, though that's not hard to add.


## Data sources

Fishbase: RopenSci hosts Parquet files of Fishbase species, families, and synonyms. The Python script accesses those directly.

Worms: I downloaded a relatively recent dump of WoRMS from GBIF at https://www.gbif.org/dataset/2d59e5db-57ad-41ff-97d6-11f5fb264527 and extracted the species names from the file `taxon.txt`. That file is included here (data/worms_species.txt). `grep -P '\tSpecies\t' taxon.txt | grep -w 'accepted' | cut -d'    ' -f 7,8,11,12,13,14,15,16,17,18 | gzip > worms_species.txt.gz`

NCBI: The script will download the most recent taxdump from NCBI.

When the script runs it will cache these files in a folder named `cache/`. Delete that folder to start with a clean slate.

## Input

The input is blast-output, tabular, using this output format:

     -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp"

--cutoff changes how lenient the LCA calculation it is - the larger the cutoff, the more species are included in an LCA. By default this is 1 - meaning that a species with 98% identity and another species with 98.5% identity are both included in the LCA, as they are within 1% of each other's identities.

--pident changes how the BLAST results are parsed, hits below that cutoff will never make it into the LCA. Default is 90.

--min_coverage changes how the BLAST results are parsed, hits below that query coverage will be ignored. Default is 90.

--missing_out changes the filename of the file missing species are written to, by default 'missing.csv'. Missing species are BLAST result lines where we couldn't find anything in Fishbase, nor in WoRMS, nor in the NCBI Taxonomy. Ideally this file should be empty - if there are many rows in this file, something may have gone wrong (missing NCBI taxonomy IDs in the BLAST output?).

--no_normalise_identity changes the behaviour of the LCA calculation and turns off the normalisation using the query coverage. Instead, this flag use only the bp identity. Using this flag means that the script will behave identically to the eDNAFlow LCA script.

## Output

Looks like this:

```
ASV_name        Class   Order   Family  Genus   Species PercentageID    Coverage	Species_In_LCA	Source
ASV_17067       Teleostei       Ophidiiformes   Ophidiidae      dropped dropped 89.60   100.00	Ventichthys biospeedoi, Bassozetus zenkevitchi	worms
ASV_17079       Teleostei       Ovalentaria incertae sedis      Pomacentridae   Acanthochromis  Acanthochromis polyacanthus     79.03   90.00	Acanthochromis polyacanthus	worms
ASV_17100       Teleostei       Centrarchiformes        Aplodactylidae  Crinodus        Aplodactylus lophodon   100.00  100.00	Aplodactylus lophodon	worms
ASV_17102       Teleostei       Anguilliformes  Muraenidae      Gymnothorax     Gymnothorax prasinus    99.02   99.00	Gymnothorax prasinus	fishbase
ASV_17176       Teleostei       Myctophiformes  Myctophidae     Symbolophorus   Symbolophorus evermanni 89.11   100.00	Symbolophorus evermanni	fishbase
ASV_17291       Teleostei       Stomiiformes    Sternoptychidae Valenciennellus Valenciennellus tripunctulatus  83.87   90.00	Valenciennellus tripunctulatus	ncbi
ASV_17546       Teleostei       Ophidiiformes   Ophidiidae      dropped dropped 76.22   80.00	Lepophidium profundorum, Genypterus chilensis, Genypterus tigerinus, Genypterus blacodes, Genypterus capensis, Apagesoma australe	ncbi
```

A tab-delimited table, one row per unique query in the BLAST results, showing which Fishbase taxonomic levels were included, and which were dropped. It also shows the highest BLAST identity of the species-hits included in the LCA, the highest query coverage, and the species that were included in the LCA. *IMPORTANT*: By default BLAST does not report queries with no hits. That means the output table of this script will not contain all queries.

## Output when using calculateLCAWithFishbase_FAIReCompatible.py

The output file contains the columns:
'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'OTU', 'numberOfUnq_BlastHits', '%ID', 'species_in_LCA', 'sources', and one column for each sample.

The taxaRaw and taxaFinal files both contain the columns:
'seq_id', 'dna_sequence', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'specificEpithet', 'scientificName', 'scientificNameAuthorship', 'taxonRank', 'taxonID', 'taxonID_db', 'verbatimIdentification', 'accession_id', 'accession_id_ref_db', 'percent_match', 'percent_query_cover', 'confidence_score', and 'identificationRemarks'


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

## CHANGELOG

- v0.04: moved Claude code to main script, now filters by 90% query coverage and normalises bp identity by coverage before LCA calculation
- v0.03: added script for better compatibility with the FAIRe metadata standards
- v0.02: adding NCBI taxonomy
- v0.01: initial release a.k.a. 'works on my system'


## TODO

- parsing out the genus name from the lines is very iffy. I'm sure there will be cases where it's just 'H. sapiens', or '[bla]Homo sapiens'. I'll probably find many cases where I have to add more rules to cleaning the strings. I could use taxonerd but that feels a bit like overkill. Could also store all fishbase/worms/NCBI species names in a lookup table and then align with every BLAST result line.
- the tests rely on unittest. Could go 'pure' pytest, too.... not too hard?
