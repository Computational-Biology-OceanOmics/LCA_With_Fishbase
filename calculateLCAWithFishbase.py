from collections import OrderedDict, Counter
import pandas as pd
from argparse import ArgumentParser
import statistics

CUTOFF = 1

def get_lca(entries):
    # given a list of entries -could be several species, or several genera,
    # and their percentage identity,
    # pull out the highest percentage and remove all hits with pident - 1 (CUTOFF)
    # an entry can be a species or a genus, it comes with its percent identity
    sorted_spec = sorted(entries)
    top_spec = sorted_spec[-1]
    bottom_spec = sorted_spec[0]
    top_perc = top_spec[0]
    perc_range = top_perc - cutoff
    new_spec = set()
    all_percentages = []
    for s in sorted_spec:
        perc, entry = s
        #print(perc, entry)
        if perc >= perc_range:
            new_spec.add(entry)
            all_percentages.append(perc)

    lca_perc = statistics.mean(all_percentages)
    # do we have only one species/genus/family/order?
    if len(new_spec) == 1:
        # easy!
        lca = list(new_spec)[0]
    else:
        # we have several - go up one level
        lca = 'dropped'
    return lca_perc, lca, new_spec


parser = ArgumentParser(description = 'Parses a BLAST-tabular output file and produces LCAs by asking the Fishbase API for each hit\'s lineage. BLAST formatting assumed is: -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp"')
parser.add_argument('-f', '--file', help = 'input file of BLAST results', required = True)
parser.add_argument('-o', '--output', help = 'Output file of LCAs. Tab delimited.', required = True)
parser.add_argument('--cutoff', help = 'OPTIONAL: Percentage cutoff between best BLAST hit and followup to be considered in LCA. Only species within this percentage identity cutoff will be included in LCA calculation.\nDefault: %(default)s in line with eDNAFlow\'s LCA script.', default = CUTOFF, type = float)
parser.add_argument('--pident', help = 'OPTIONAL: Percentage cutoff for BLAST hits. Hits below this cutoff will be ignored for LCA calculation.\nDefault: Initially consider all BLAST hits, but then filter in line with --cutoff.', default = 0, type = float)
parser.add_argument('--missing_out', help = 'OPTIONAL: Filename to write missing species (not in Fishbase) to.\nDefault: %(default)s.', default = 'missing.csv')
args = parser.parse_args()

cutoff = args.cutoff
pident_cutoff = args.pident

assert pident_cutoff >= 0 and pident_cutoff <= 100, 'ERROR: Parameter --pident_cutoff must be between 0 and 100'

asv_hits = OrderedDict()
missing_c = Counter()

# download all species from fishbase
all_fishbase = pd.read_parquet("https://fishbase.ropensci.org/fishbase/species.parquet") 
# and get their taxonomic ranks, too
fams = pd.read_parquet("https://fishbase.ropensci.org/fishbase/families.parquet")

synonyms = pd.read_parquet('https://fishbase.ropensci.org/fishbase/synonyms.parquet')

merged = pd.merge(all_fishbase, fams, on = 'FamCode')[ ['SpecCode', 'Genus', 'Species_x', 'Family', 'Order', 'Class']]
merged = merged.rename(columns = {"Species_x": 'Species'})
merged['Species'] = merged['Genus'] + ' ' + merged['Species']
merged = merged.reset_index() # just doubly making sure

fishbase_genera_to_lineage = {}
for index, row in merged.iterrows():
    fishbase_genera_to_lineage[row['Genus']] = [row['Family'], row['Order'], row['Class']]

speccode_to_row = dict(zip(merged['SpecCode'], merged['Species']))

all_synonyms = dict(zip(synonyms['SynGenus'] + ' ' + synonyms['SynSpecies'], synonyms['SpecCode']))


worms = pd.read_csv('worms_species.txt.gz', sep = '\t', header = None, names = ['Species', 'Genus', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'g', 'trash', 'sp'])
worms_genera_to_lineage = {}
for index, row in worms.iterrows():
    worms_genera_to_lineage[row['Genus']] = [row['Family'], row['Order'], row['Class']]

worms_species = set(worms['Species'])
# Aaadonta angaurana      Aaadonta        Animalia        Mollusca        Gastropoda      Stylommatophora Endodontidae    Aaadonta                angaurana

look_up = {}

with open(args.missing_out, 'w') as out:
    for line in open(args.file):
        # LOADS of typos in our species names. fixing them here.
        line = line.replace('Petroschmidtia albonotatus', 'Petroschmidtia albonotata')

        ll = line.rstrip().split('\t')

        pident = float(ll[6])
        if pident < pident_cutoff:
            continue
        ll = line.rstrip().split()

        # We need to find the scientific name of this species. The problem is that I've seen every possible variant of species
        # names in this text.
        # we therefore go through the line word by word, and check if it is in our big genus list. Then the next item has to be a species - in some rare cases, there are two species names.
        # in other cases the species is also not in fishbase.
        found_genus = False
        found_species = False
        in_fishbase = False
        in_worms = False
        for index, element in enumerate(ll):
            if element in fishbase_genera_to_lineage:
                genus = element
                thisspecies = ll[index+1]
                found_genus = True
                in_fishbase = True
                break
            if index+1 == len(ll):
                # are we at the end?
                break

            # different case - we have synonyms! those are not in fishbase' species table
            # but in fishbase synonyms table. so we have to take those hits
            # and link them back to the fishbase species table
            if (element + ' ' + ll[index+1]) in all_synonyms:
                correct_spec_code = all_synonyms[ element + ' ' + ll[index+1]]
                correct_species = speccode_to_row[correct_spec_code]
                genus, thisspecies = correct_species.split(' ')
                found_genus = True
                in_fishbase = True
                break

            # found in nothing - is it in Worms?
            if element in worms_genera_to_lineage:
                genus = element
                thisspecies = ll[index+1]
                found_genus = True
                in_worms = True

                break

        if not found_genus:
            out.write(line)
            continue

        species = f'{genus} {thisspecies}'

        source = 'none'
        if in_fishbase:
            #lineage = merged[merged['Genus'] == genus].head(1)
            lineage = fishbase_genera_to_lineage[genus]
            # the lineage's species may be wrong but that's ok - we only look by genus
            source = 'fishbase'
        if in_worms:
            lineage = worms_genera_to_lineage[genus]
            #lineage = worms[worms['Genus'] == genus].head(1)

            source = 'worms'
        

        family, order, thisclass = lineage

        lineage = [ ("C", thisclass),
                            ("O", order),
                            ("F", family),
                            ("G", genus),
                            ("S", species) ]

        if ll[0] not in asv_hits:
            asv_hits[ll[0]] = [ (source, pident, lineage) ]
        else:
            asv_hits[ll[0]].append( (source, pident, lineage) )


with open(args.output, 'w') as out:
    out.write('ASV_name\tClass\tOrder\tFamily\tGenus\tSpecies\tPercentageID\tSpecies_In_LCA\n')
    for asv_name in asv_hits:
        orf_hits = asv_hits[asv_name]

        # get all the species
        species = set()
        genera = set()
        families = set()
        orders = set()
        classes = set()
        sources = set()
        for a in orf_hits:
            # this is one HIT it has all the levels
            source, pident, lineage = a
            thisclass, thisorder, thisfamily, thisgenus, thisspecies = [x[1] for x in lineage]
            classes.add( (pident, thisclass) )
            orders.add( (pident, thisorder) )
            families.add( (pident, thisfamily) )
            genera.add( (pident, thisgenus) )
            species.add( (pident, thisspecies) )
            sources.add(source)


        #[('C', 'Actinopterygii'), ('O', 'Ophidiiformes'), ('F', 'Ophidiidae'), ('G', 'Ventichthys'), ('S', 'Ventichthys biospeedoi')]
        # oK now we have all the species and genera

        # we need to find the specs that has >1% difference in pident
      
        lca_spec_perc, lca_spec, included_spec = get_lca(species)
        lca_genus_perc, lca_genus, included_genera = get_lca(genera)
        lca_fam_perc, lca_fam, _ = get_lca(families)
        lca_order_perc, lca_order, _ = get_lca(orders)
        lca_class_perc, lca_class, _ = get_lca(classes)
        out.write(f'{asv_name}\t{lca_class}\t{lca_order}\t{lca_fam}\t{lca_genus}\t{lca_spec}\t{lca_spec_perc:.2f}\t{", ".join(included_spec)}\t{", ".join(sources)}\n')
