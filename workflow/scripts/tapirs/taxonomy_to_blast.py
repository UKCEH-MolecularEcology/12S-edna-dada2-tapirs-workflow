# adapted from taxopy core.py functions (https://github.com/apcamargo/taxopy)
# includes merged taxid dictionary for taxid corrections adapted from Simple-LCA
# updated to robustly parse NCBI taxdump files using "|" delimiters

nodes_dmp  = snakemake.params.taxdump + "/nodes.dmp"
names_dmp  = snakemake.params.taxdump + "/names.dmp"
merged_dmp = snakemake.params.taxdump + "/merged.dmp"

blast     = snakemake.input.blast
blast_tax = snakemake.output.blast_tax


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
    lineage      = []
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
        taxdict["species"] = taxdict["species"].split("/")[0]
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


if len([1 for line in open(blast)]) > 0:
    taxid2parent, taxid2rank = import_nodes()
    taxid2name               = import_names()
    merged_dict              = merged_taxonomy()

    with open(blast, "r") as blasthits, open(blast_tax, "w") as output:
        for hit in blasthits:
            fields = hit.rstrip("\n").split("\t")
            if len(fields) < 4:
                output.write(hit.strip() + "\tunknown/unknown/unknown/unknown/unknown/unknown/unknown\n")
                continue
            taxid = fields[3]
            if taxid == "N/A":
                output.write(hit.strip() + "\tunknown/unknown/unknown/unknown/unknown/unknown/unknown\n")
            else:
                if ";" not in str(taxid):
                    try:
                        output.write(hit.strip() + "\t" + "/".join(taxonomy_string(int(taxid))) + "\n")
                    except (KeyError, ValueError):
                        try:
                            taxid = merged_dict[taxid]
                            output.write(hit.strip() + "\t" + "/".join(taxonomy_string(int(taxid))) + "\n")
                        except (KeyError, ValueError):
                            output.write(hit.strip() + "\tunknown/unknown/unknown/unknown/unknown/unknown/unknown\n")
                else:
                    output.write(hit.strip() + "\tunknown/unknown/unknown/unknown/unknown/unknown/unknown\n")
else:
    open(blast_tax, "w").close()
