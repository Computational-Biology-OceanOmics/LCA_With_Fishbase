# Calculate LCA using Fishbase, then WoRMS and then NCBI

Calculates LCA using the Fishbase taxonomy - if that fails, the WoRMS taxonomy - if that fails, the NCBI Taxonomy.

Most LCA pipelines out there rely on NCBI's Taxonomy database. However, that database is large and fish are notoriously prone to change. Other databases like Fishbase are updated more often, so here is a tool that takes a table of BLAST results, takes the hits for each ASV, and queries the Fishbase database to ask for the 'updated'/'current' lineage of the hit, then uses the Fishbase lineages to calculate LCAs. Sometimes several similar species have wrong, outdated families on NCBI Taxonomy that are correct on Fishbase, and when you use NCBI Taxonomy, you will get an order-level LCA, while Fishbase-based LCA will have a family or genus-level LCA. Fishbase also tracks outdated synonyms, so we might get more accurate species-level names.

What this script does:

1. For every BLAST hit, check the lineage using Fishbase. If all hits for a query are in Fishbase, calculate a simple LCA using those hits.
2. If the hit is not in Fishbase, check WoRMS, and calculate the LCA.
3. If the hit is in neither Fishbase nor WoRMS, use NCBI Taxonomy.

Since these databases often changes you better write down the date you ran this tool with your data.

# Usage

    python calculateLCAWithFishbase_Claude.py -f blast_results.tsv -o lca_results.tsv --pident 97

# Dependencies

    pip install pandas pyarrow fastparquet

Any version should be fine?

# Why that Claude in the filename?

Most of this code is AI-written. I've written an initial prototype (minus the NCBI Taxonomy!!!) and gave that to Claude for hardening and making production-ready. The orignal, non-AI code lives at calculateLCAWithFishbase.py. Same usage if you don't want to trust our AI overlords (I did fix a few bugs Claude had introduced).

     python calculateLCAWithFishbase.py -f blast_results.tsv -o lca_results.tsv --pident 97


# *IMPORTANT*

- The BLAST results need the taxonomy ID. Make sure that the taxonomy ID is included in the BLAST output and not just N/A. 
- Make sure that the BLAST results are formatted correctly; see below in the Input section.
- Adjust the percent identity (--pident), by default this script includes everything with >= 90% identity. That may be too lenient.

## Method: Species ID:

What the script does:
- go through the tabular blast results, check if every word is a valid genus in 1. Fishbase, 2. Fishbase synonyms, 3. Worms. We want to trust the Fishbase taxonomy the most but not every species we hit is in Fishbase. Mammals etc. will instead hit into Worms.
- if the word is a genus, it follows that the next word is a species.
- using the genus name, ask Fishbase what the taxonomy for that genus is. If the genus is not in Fishbase, ask Worms. If the genus is not in Worms, write to the MISSING_OUT file and ignore from LCA calculation.
- if both Fishbase and WoRMS were not found in the row, assume that the third column is the NCBI taxonomy ID. Use that to look up the lineage instead.

## Method: LCA calculation

The LCA calculation works the same way as [eDNAFlow](https://github.com/mahsa-mousavi/eDNAFlow)'s LCA calculation. Given a group of potential species for an ASV, take the species with the highest percentage base-pair identity, subtract 1 from the identity, and then include all species above that cutoff in the LCA. The LCA itself is just a grouping: if there are several species within the cutoff, then the species is set to 'dropped' and we go up one taxonomic level, repeat for the genus, repeat for the family, repeat for the class, repeat for the order. There's no LCA voting or similar, though that's not hard to add.


## Data sources

Fishbase: RopenSci hosts Parquet files of Fishbase species, families, and synonyms. The Python script accesses those directly.

Worms: I downloaded a relatively recent dump of WoRMS from GBIF at https://www.gbif.org/dataset/2d59e5db-57ad-41ff-97d6-11f5fb264527 and extracted the species names from the file `taxon.txt`. That file is included here (worms_species.txt). `grep -P '\tSpecies\t' taxon.txt | cut -d'    ' -f 7,8,11,12,13,14,15,16,17,18 | gzip > worms_species.txt.gz`

NCBI: The script will download the most recent taxdump from NCBI.

## Input

The input is blast-output, tabular, using this output format:

     -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp"


--cutoff changes how lenient the LCA calculation it is - the larger the cutoff, the more species are included in an LCA.

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

## Installation

Any fairly recent Python should be OK - we only need `pandas` installed.

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
