# Calculate LCA using Fishbase

Most LCA pipelines out there rely on NCBI's Taxonomy database. However, that database is large and fish are notoriously prone to change. Other databases like Fishbase are updated more often, so here is a tool that takes a table of BLAST results, takes the hits for each ASV, and queries the Fishbase API to ask for the 'updated'/'current' lineage of the hit, then uses the Fishbase lineages to calculate LCAs. Sometimes several similar species have wrong, outdated families on NCBI Taxonomy that are correct on Fishbase, and when you use NCBI Taxonomy, you will get an order-level LCA, while Fishbase-based LCA will have a family or genus-level LCA. Fishbase also tracks outdated synonyms, so 

Be careful, though: Fishbase includes *fish* will return nothing for non-fish species. In these cases we go back to World Register of Marine Species (WoRMS) for the whales and the sealions. If you have a mix of marine and non-marine species in your results better use a different database. We work with marine eDNA so filtering out non-marine hits is very useful for us, lots of similar looking freshwater fish out there!

Since Fishbase often changes you better write down the date you ran this tool with your data.

## Method: Species ID

Normally, you can rely on NCBI taxonomy IDs and it's simple enough to pull them from BLAST results. However, we don't have NCBI taxonomy IDs for many species we work with. We work with external partners' internal databases that have not been published yet, and those databases include many firsts.

Instead this script does this:

- go through the tabular blast results, check if every word is a valid genus in 1. Fishbase, 2. Fishbase synonyms, 3. Worms. We want to trust the Fishbase taxonomy the most but not every species we hit is in Fishbase. Mammals etc. will instead hit into Worms.
- if the word is a genus, it follows that the next word is a species.
- using the genus name, ask Fishbase what the taxonomy for that genus is. If the genus is not in Fishbase, ask Worms. If the genus is not in Worms, write to the MISSING_OUT file and ignore from LCA calculation.

## Method: LCA calculation

The LCA calculation works the same way as [eDNAFlow](https://github.com/mahsa-mousavi/eDNAFlow)'s LCA calculation. Given a group of potential species for an ASV, take the species with the highest percentage base-pair identity, subtract 1 from the identity, and then include all species above that cutoff in the LCA. The LCA itself is just a grouping: if there are several species within the cutoff, then the species is set to 'dropped' and we go up one taxonomic level, repeat for the genus, repeat for the family, repeat for the class, repeat for the order. There's no LCA voting or similar, though that's not hard to add.


## Data sources

Fishbase: RopenSci hosts Parquet files of Fishbase species, families, and synonyms. The Python script accesses those directly.

Worms: I downloaded a relatively recent dump of WoRMS from GBIF at https://www.gbif.org/dataset/2d59e5db-57ad-41ff-97d6-11f5fb264527 and extracted the species names from the file `taxon.txt`. That file is included here (worms_species.txt). `grep -P '\tSpecies\t' taxon.txt | cut -d'    ' -f 7,8,11,12,13,14,15,16,17,18 | gzip > worms_species.txt.gz`

## Input

The input is blast-output, tabular, using this output format:

     -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp"

I think it's OK as long as the third column is the NCBI Taxonomy ID and the seventh column is the percent identity.

## Usage


python calculateLCAWithFishbase.py -f asv_blastn_results.tsv -o LCA_results.tsv


```
usage: calculateLCAWithFishbase.py [-h] -f FILE -o OUTPUT [--cutoff CUTOFF] [--pident PIDENT] [--missing_out MISSING_OUT]

Parses a BLAST-tabular output file and produces LCAs by asking the Fishbase API for each hit's lineage. BLAST formatting assumed is: -outfmt "6
qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue
bitscore qcovs qcovhsp"

options:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  input file of BLAST results
  -o OUTPUT, --output OUTPUT
                        Output file of LCAs. Tab delimited.
  --cutoff CUTOFF       OPTIONAL: Percentage cutoff between best BLAST hit and followup to be considered in LCA. Only species within this
                        percentage identity cutoff will be included in LCA calculation. Default: 1 in line with eDNAFlow's LCA script.
  --pident PIDENT       OPTIONAL: Percentage cutoff for BLAST hits. Hits below this cutoff will be ignored for LCA calculation. Default:
                        Initially consider all BLAST hits, but then filter in line with --cutoff.
  --missing_out MISSING_OUT
                        OPTIONAL: Filename to write missing species (not in Fishbase) to. Default: missing.csv.
```

--cutoff changes how lenient the LCA calculation it is - the larger the cutoff, the more species are included in an LCA.

--pident changes how the BLAST results are parsed, hits below that cutoff will never make it into the LCA.

--missing_out changes the filename of the file missing species are written to, by default 'missing.csv'. Missing species are species present in the BLAST results and on NCBI Taxonomy, but are not in Fishbase. These species are a mix of non-marine species and species with sp./cf. names.

## Output

Looks like this:

```
ASV_name        Class   Order   Family  Genus   Species PercentageID    Species_In_LCA
ASV_17067       Teleostei       Ophidiiformes   Ophidiidae      dropped dropped 89.60   Ventichthys biospeedoi, Bassozetus zenkevitchi
ASV_17079       Teleostei       Ovalentaria incertae sedis      Pomacentridae   Acanthochromis  Acanthochromis polyacanthus     79.03   Acanthochromis polyacanthus
ASV_17100       Teleostei       Centrarchiformes        Aplodactylidae  Crinodus        Aplodactylus lophodon   100.00  Aplodactylus lophodon
ASV_17102       Teleostei       Anguilliformes  Muraenidae      Gymnothorax     Gymnothorax prasinus    99.02   Gymnothorax prasinus
ASV_17176       Teleostei       Myctophiformes  Myctophidae     Symbolophorus   Symbolophorus evermanni 89.11   Symbolophorus evermanni
ASV_17291       Teleostei       Stomiiformes    Sternoptychidae Valenciennellus Valenciennellus tripunctulatus  83.87   Valenciennellus tripunctulatus
ASV_17546       Teleostei       Ophidiiformes   Ophidiidae      dropped dropped 76.22   Lepophidium profundorum, Genypterus chilensis, Genypterus tigerinus, Genypterus blacodes, Genypterus capensis, Apagesoma australe
```

A tab-delimited table, one row per unique query in the BLAST results, showing which Fishbase taxonomic levels were included, and which were dropped. It also shows the average BLAST identity of the species-hits included in the LCA, and the species that were included in the LCA. *IMPORTANT*: By default BLAST does not report queries with no hits. That means the output table of this script will not contain all queries.

## Installation

    conda install -c bioconda pytaxonkit
    pip install pyworms

## FAQ

- Why does this depend on pytaxonkit? I thought it ignores NCBI Taxonomy?

Yes, but I have not found a nice way to pull out the species name from a row of BLAST hits. I have played with [taxonerd](https://github.com/nleguillarme/taxonerd) but it feels like overkill for this problem. Relying on the NCBI Taxonomy ID for the species name, then pulling out the Fishbase lineage for the NCBI species name is the best middle-ground in my opinion. I *could* base this on the sscinames column produced by BLAST, but when BLAST does not have the NCBI taxonomy database present when running that column will be all NA.

- How long does this run for?

In my tests for a whole eDNA dataset with a few thousand ASVs, about 15 minutes. The Fishbase API accepts only 500 species at a time so it takes most of the running time to pull down all lineages.

- You could make this faster using a dump of the Fishbase database!?

Yes but then I have to distribute that database dump with the code, and that dump will be outdated within a week, but others will find the database dump on Github and base their own scripts on the outdated dump and I'd rather not have that happen.

- Can this take other taxonomic databases?

I have thought about adding AFD (Australian Faunal Directory) as another source. Maybe later.

- I have weird taxonomic sub-levels that I am interested in, why are they dropped?

The script only works with class, order, family, genus, species. That way I can make the LCA calculation fairly lazy: instead of having to build a graph of all taxonomic levels and doing weird graph-based magic, I can just calculate the sets of unique species, unique genera, unique families, unique orders, unique classes, and check set size after filtering. If the set size is > 1, set this taxonomic level's taxonomic label to 'dropped'. Easier than breaking my head trying to come up with recursive tree-walking algorithms for the sake of methodological complexity you can publish in a paper, I'd rather have results.

- My BLAST results have NAs for taxonomy IDs so everything breaks.

Your BLAST database was not made using taxonomy IDs. Make a table of sequence IDs and their associated taxonomy IDs, then run this:

    makeblastdb -dbtype nucl -in your_database.fasta -parse_seqids -taxid_map your_database.taxids.txt

- I may have unacccepted taxonomic names in my results.

Fishbase is nice in that it knows about many, but not all unaccepted names. If an unaccepted name made it into Fishbase the API returns the accepted version of this name along with the taxonomic lineage of the new, accepted name. We then use that name for the LCA. Neat, isn't it!  
On the other hand, if your species is not in Fishbase it will be dropped from the LCA calculation and it will get written to missing.csv.

- I have more questions!

Please contact me at pbayer AT minderoo.org


## CHANGELOG

- v0.01: initial release a.k.a. 'works on my system'
