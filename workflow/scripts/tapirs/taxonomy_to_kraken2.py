# adapted from taxopy core.py functions (https://github.com/apcamargo/taxopy)
# includes merged taxid dictionary for taxid corrections adapted from Simple-LCA
# updated to robustly parse NCBI taxdump files using "|" delimiters

from itertools import dropwhile

nodes_dmp  = snakemake.params.taxdump + "/nodes.dmp"
names_dmp  = snakemake.params.taxdump + "/names.dmp"
merged_dmp = snakemake.params.taxdump + "/merged.dmp"

kraken     = snakemake.input.kraken2
kraken_tax = snakemake.output.kraken2_tax


def import_nodes():
    taxid2parent = {}
    taxid2rank   = {}
    with open(nodes_dmp, "r") as file:
        for line in file:
            parts = [x.strip() for x in line.split("|")]
            if len(parts) < 3:
                continue
            try:
                taxid  = int(parts[0])
                parent = int(parts[1])
                rank   = parts[2]
            except ValueError:
                continue
            taxid2parent[taxid] = parent
            taxid2rank[taxid]   = rank
    return taxid2parent, taxid2rank


def import_names():
    taxid2name = {}
    with open(names_dmp, "r") as file:
        for line in file:
            parts = [x.strip() for x in line.split("|")]
            if len(parts) < 4:
                continue
            try:
                taxid = int(parts[0])
            except ValueError:
                continue
            name       = parts[1]
            name_class = parts[3]
            if name_class == "scientific name":
                taxid2name[taxid] = name
    return taxid2name


def merged_taxonomy():
    merged_dict = {}
    with open(merged_dmp, "r") as merged:
        for line in merged:
            parts     = [x.strip() for x in line.split("|")]
            if len(parts) < 2:
                continue
            old_taxid = parts[0]
            new_taxid = parts[1]
            if old_taxid and new_taxid:
                merged_dict[old_taxid] = new_taxid
    return merged_dict


def get_lineage(taxid):
    lineage       = []
    current_taxid = taxid
    while True:
        lineage.append(current_taxid)
        if current_taxid not in taxid2parent:
            break
        parent = taxid2parent[current_taxid]
        if parent == current_taxid:
            break
        current_taxid = parent
    return lineage


def rank_name_dictionary(taxid):
    lineage        = get_lineage(taxid)
    rank_name_dict = {}
    for lineage_taxid in lineage:
        rank = taxid2rank.get(lineage_taxid, "no rank")
        name = taxid2name.get(lineage_taxid)
        if rank != "no rank" and name is not None:
            rank_name_dict[rank] = name
    return rank_name_dict


def taxonomy_string(taxid):
    ranks   = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    taxdict = rank_name_dictionary(taxid)

    if "species" in taxdict:
        taxdict["species"] = "_".join(taxdict["species"].split(" ")[0:2])
        if "_sp." in taxdict["species"]:
            del taxdict["species"]

    if "species" in taxdict and "genus" not in taxdict:
        del taxdict["species"]

    tax_ranks = []
    for rank in ranks:
        if rank in taxdict:
            tax_ranks.append(taxdict[rank])
        else:
            tax_ranks.append("unknown")
    return tax_ranks


def strip_unknown(l):
    l = list(dropwhile(lambda x: x == "unknown", l[::-1]))
    return l[::-1]


taxonomy = ("domain", "phylum", "class", "order", "family", "genus", "species")
header   = "query\ttax_rank\totu_id\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\n"

if len([1 for line in open(kraken)]) > 0:
    taxid2parent, taxid2rank = import_nodes()
    taxid2name               = import_names()
    merged_dict              = merged_taxonomy()

    with open(kraken, "r") as kraken_in, open(kraken_tax, "w") as output:
        output.write(header)
        for hit in kraken_in:
            fields = hit.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            if fields[0] == "C":
                taxid = fields[2]
                if taxid == "1":
                    continue
                try:
                    tax_str = strip_unknown(taxonomy_string(int(taxid)))
                except (KeyError, ValueError):
                    if taxid in merged_dict:
                        try:
                            tax_str = strip_unknown(taxonomy_string(int(merged_dict[taxid])))
                        except Exception:
                            continue
                    else:
                        continue
                if not tax_str:
                    continue
                tax_rank = taxonomy[len(tax_str) - 1]
                otu_id   = tax_str[-1]
                tax_str  = tax_str + (["unidentified"] * (7 - len(tax_str)))
                tax_str  = "\t".join([str(x) for x in tax_str])
                output.write("%s\t%s\t%s\t%s\n" % (fields[1].strip(), tax_rank, otu_id, tax_str))
else:
    with open(kraken_tax, "w") as output:
        output.write(header)
