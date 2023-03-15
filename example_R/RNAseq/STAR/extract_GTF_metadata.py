

"""
When GenomicFeatures::makeTxDbFromGFF is run to create an exons by gene file for summarizeOverlaps, it only 
keeps a small subset of the metadata associated with each feature.
It also uses the "gene_id" as the primary identifier.
This script pulls all unique gene_id's from a gtf, then finds all of the associated keys and greats a lookup.
https://uswest.ensembl.org/info/website/upload/gff.html
seqname, source, feature, start, end, score, strand, frame, attribute
"""


fields = []
gene_ids = {}


i = 0

for line in open("GCF_002892825.2_ISE6_asm2.2_deduplicated_genomic.gtf", 'r'):
    i += 1
    if i % 100000 == 0: print(i)
    if len(line) == 0: continue  # skip empty
    if line[0] == '#': continue  # skip header / comment lines
    l2 = line.strip().split('\t')
    attributes = l2[8].split(';')
    line_values = {}
    for attribute in attributes:
        a = attribute.split(' "')
        if len(a) == 2:
            k = a[0].strip()
            v = a[1].strip().strip('"')
            line_values[k] = v
            if k not in fields: fields.append(k)
    if 'gene_id' in line_values:
        gene_id = line_values['gene_id']
        if gene_id not in gene_ids: gene_ids[gene_id] = {}
        for k,v in line_values.items():
            if k not in gene_ids[gene_id]:
                gene_ids[gene_id][k] = v
        
print(i)

with open("gene_id-lookup.tsv", 'w') as outf:
    outf.write('\t'.join(fields) + '\n')
    for k in gene_ids.keys():
        outf.write('\t'.join([gene_ids[k].get(x, '') for x in fields]) + '\n')

