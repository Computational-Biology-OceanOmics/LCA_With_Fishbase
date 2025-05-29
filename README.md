# Calculate LCA using Fishbase, then WoRMS and then NCBI

This tool calculates Lowest Common Ancestor (LCA) assignments for BLAST results using multiple taxonomy databases in order of preference: Fishbase → WoRMS → NCBI Taxonomy.

Most LCA scripts rely solely on NCBI's Taxonomy database. However, fish taxonomy changes frequently, and NCBI often contains outdated classifications. Fishbase and WoRMS are updated more regularly and provide more accurate lineage information for marine organisms. Plus, Fishbase follows mostly Betancur-R, and what more could one ask for?

The tool processes BLAST results through a three-step approach:

**Fishbase first**: Queries Fishbase taxonomy for each BLAST hit. If all hits for a query sequence are found, calculates LCA using Fishbase lineages.  
**WoRMS fallback**: For hits not found in Fishbase, queries the World Register of Marine Species (WoRMS) database.  
**NCBI backup**: Only uses NCBI Taxonomy when species aren't found in either Fishbase or WoRMS.  

When a query hits several species in different databases, a mix of the above may be used.

# Input/Output

**Input**: a table of BLAST hits (see below for format)

**Output**: a table of LCAs for every query in the BLAST table.

# Usage

    python calculateLCAWithFishbase_Claude.py -f blast_results.tsv -o lca_results.tsv --pident 97

# Dependencies

    pip install pandas pyarrow fastparquet

Any fairly recent version (2023-2025) should be fine, I believe.

# Why that Claude in the filename?

I've written an initial prototype (minus the NCBI Taxonomy!!!) of this code and then gave that to Claude for hardening and making production-ready. I then fixed a bunch of Claude-introduced bugs, removed most of the unnecessary comments, and added some minor bugfixes that my original implementation had already. The orignal, non-AI code lives at calculateLCAWithFishbase.py. Same usage if you don't want to trust our AI overlords.

     python calculateLCAWithFishbase.py -f blast_results.tsv -o lca_results.tsv --pident 97

# *IMPORTANT*

- The BLAST results need the taxonomy ID. Make sure that the taxonomy ID is included in the BLAST output and not just N/A. 
- Make sure that the BLAST results are formatted correctly; see below in the Input section.
- Adjust the percent identity (--pident), by default this script includes everything with >= 90% identity. That may be too lenient.
- Fishbase, WoRMS, and NCBI Taxonomy change often. Write down the date you ran this tool.
- Some sequences on NCBI or other databases do not have specific taxonomic labels, such as 'Carangidae sp.'. These lead to very high-level LCAs, obviously. Consider removing them before running this script.
- The script does no filtering for alignment length or anything else, just percent identity of the alignment. Consider filtering by query coverage % at least.

## Method: Species ID:

What the script does:
- go through the tabular blast results, check if every word is a valid genus in 1. Fishbase, 2. Fishbase synonyms, 3. Worms. We want to trust the Fishbase taxonomy the most but not every species we hit is in Fishbase. Mammals etc. will instead hit into Worms.
- using the genus name, ask Fishbase what the taxonomy for that genus is. If the genus is not in Fishbase, ask Worms. If the genus is not in Worms, write to the MISSING_OUT file and ignore from LCA calculation.
- if both Fishbase and WoRMS were not found in the row, assume that the third column is the NCBI taxonomy ID. Use that to look up the lineage instead.

## Method: LCA calculation

The LCA calculation works the same way as [eDNAFlow](https://github.com/mahsa-mousavi/eDNAFlow)'s LCA calculation. Given a group of potential species for an ASV, take the species with the highest percentage base-pair identity, subtract 1 from the identity, and then include all species above that cutoff in the LCA. The LCA itself is just a grouping: if there are several species within the cutoff, then the species is set to 'dropped' and we go up one taxonomic level, repeat for the genus, repeat for the family, repeat for the class, repeat for the order. There's no LCA voting or similar, though that's not hard to add.


## Data sources

Fishbase: RopenSci hosts Parquet files of Fishbase species, families, and synonyms. The Python script accesses those directly.

Worms: I downloaded a relatively recent dump of WoRMS from GBIF at https://www.gbif.org/dataset/2d59e5db-57ad-41ff-97d6-11f5fb264527 and extracted the species names from the file `taxon.txt`. That file is included here (worms_species.txt). `grep -P '\tSpecies\t' taxon.txt | grep -w 'accepted' | cut -d'    ' -f 7,8,11,12,13,14,15,16,17,18 | gzip > worms_species.txt.gz`

NCBI: The script will download the most recent taxdump from NCBI.

## Input

The input is blast-output, tabular, using this output format:

     -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp"

--cutoff changes how lenient the LCA calculation it is - the larger the cutoff, the more species are included in an LCA. By default this is 1 - meaning that a species with 98% identity and another species with 98.5% identity are both included in the LCA, as they are within 1% of each other's identities.

--pident changes how the BLAST results are parsed, hits below that cutoff will never make it into the LCA.

--missing_out changes the filename of the file missing species are written to, by default 'missing.csv'. Missing species are BLAST result lines where we couldn't find anything in Fishbase, nor in WoRMS, nor in the NCBI Taxonomy. Ideally this file should be empty - if there are many rows in this file, something may have gone wrong (missing NCBI taxonomy IDs in the BLAST output?).

## Output

Looks like this:

```
ASV_name        Class   Order   Family  Genus   Species PercentageID    Species_In_LCA	Source
ASV_17067       Teleostei       Ophidiiformes   Ophidiidae      dropped dropped 89.60   Ventichthys biospeedoi, Bassozetus zenkevitchi	worms
ASV_17079       Teleostei       Ovalentaria incertae sedis      Pomacentridae   Acanthochromis  Acanthochromis polyacanthus     79.03   Acanthochromis polyacanthus	worms
ASV_17100       Teleostei       Centrarchiformes        Aplodactylidae  Crinodus        Aplodactylus lophodon   100.00  Aplodactylus lophodon	worms
ASV_17102       Teleostei       Anguilliformes  Muraenidae      Gymnothorax     Gymnothorax prasinus    99.02   Gymnothorax prasinus	fishbase
ASV_17176       Teleostei       Myctophiformes  Myctophidae     Symbolophorus   Symbolophorus evermanni 89.11   Symbolophorus evermanni	fishbase
ASV_17291       Teleostei       Stomiiformes    Sternoptychidae Valenciennellus Valenciennellus tripunctulatus  83.87   Valenciennellus tripunctulatus	ncbi
ASV_17546       Teleostei       Ophidiiformes   Ophidiidae      dropped dropped 76.22   Lepophidium profundorum, Genypterus chilensis, Genypterus tigerinus, Genypterus blacodes, Genypterus capensis, Apagesoma australe	ncbi
```

A tab-delimited table, one row per unique query in the BLAST results, showing which Fishbase taxonomic levels were included, and which were dropped. It also shows the average BLAST identity of the species-hits included in the LCA, and the species that were included in the LCA. *IMPORTANT*: By default BLAST does not report queries with no hits. That means the output table of this script will not contain all queries.

## FAQ

- How long does this run for?

In my tests for a whole eDNA dataset with a few thousand ASVs, about 2 minutes.

- I have weird taxonomic sub-levels that I am interested in, why are they dropped?

The script only works with class, order, family, genus, species. That way I can make the LCA calculation fairly lazy: instead of having to build a graph of all taxonomic levels and doing weird graph-based magic, I can just calculate the sets of unique species, unique genera, unique families, unique orders, unique classes, and check set size after filtering. If the set size is > 1, set this taxonomic level's taxonomic label to 'dropped'. Easier than breaking my head trying to come up with recursive tree-walking algorithms for the sake of methodological complexity you can publish in a paper, I'd rather have results.

- I may have unacccepted taxonomic names in my results.

Fishbase is nice in that it knows about many, but not all unaccepted names. If an unaccepted name made it into Fishbase the API returns the accepted version of this name along with the taxonomic lineage of the new, accepted name. We then use that name for the LCA. Neat, isn't it!  
The script then checks whether the species is in WoRMS - lastly, your misspelled species name may be in the NCBI taxonomy.

- I have more questions!

Please contact me at pbayer AT minderoo.org

## CHANGELOG

- v0.02: adding NCBI taxonomy
- v0.01: initial release a.k.a. 'works on my system'


## TODO

- parsing out the genus name from the lines is very iffy. I'm sure there will be cases where it's just 'H. sapiens', or '[bla]Homo sapiens'. I'll probably find many cases where I have to add more rules to cleaning the strings. I could use taxonerd but that feels a bit like overkill.
