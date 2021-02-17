import os
from Bio import Entrez
import multiprocessing
import argparse
import sys

#### Argparse
parser = argparse.ArgumentParser(description="Assign taxonomy to fasta files by from the consensus of top blast hits.\n \
    EXAMPLE USAGE: python bca2.py example.fasta -b blastx -evalue 1e-20 -n 5")
parser.add_argument("query", help="input fasta protein or nucleotide file")
parser.add_argument("-b","-blastcmd", default="blastn", help="blast command; blastn, blastp, blastx, etc.\
    The script chooses the correct nt/nr database based on your choice.\
    Default is blastn.")
parser.add_argument("-e", "-evalue", default="1", help="e value cutoff, default is 1")
parser.add_argument("-m", "-max_target_seqs", type=int, default=10, help="The number of top hits that will be compared, default is 10")
parser.add_argument("-o", "-out", help="custom output file location")
parser.add_argument("-n", "-num_threads", type=int, default=1, help="number of threads, default is 1")
args = parser.parse_args()

#### Additional command line argument parsing
dbdict = {"blastn": "nt", "blastp": "nr", "blastx": "nr", "tblastx": "nt", "tblastn": "nt"}
try:
    db = dbdict[args.b.lower()]
except KeyError:
     sys.exit("Error parsing blast command \"{}\"; you might have a typo?".format(args.b))
if args.o is None:
    outfile = args.query.rsplit(".", 1)[0] + ".tax"
else:
    outfile = args.o

blastout = outfile.rsplit(".", 1)[0] + ".blastout"
cpu_count = multiprocessing.cpu_count()

#### Constructs blast command
print("Running the following command: \n \
{blastcmd} -query {query} -db {db} -outfmt '6 qseqid sacc staxids evalue' \
-max_target_seqs {m} -out {out} -num_threads {threads} -evalue {e}\
".format(blastcmd=args.b.lower(), query=args.query, db=db, m=args.m, out=blastout, threads=args.n, e=args.e))
sys.stdout.flush()

os.system("{blastcmd} -query {query} -db {db} -outfmt '6 qseqid sacc staxids evalue' \
-max_target_seqs {m} -out {out} -num_threads {threads} -evalue {e}\
".format(blastcmd=args.b.lower(), query=args.query, db=db, m=args.m, out=blastout, threads=args.n, e=args.e))


# parses the blastout file into two dictionaries:
# uniq_subject_taxids: all taxids, returns lineage list or list of corresponding ids
# uniq_query_ids: all query sequence names, returns list of subject taxids that were hit by that query

## parsing blastout
uniq_subject_taxids = {}
uniq_query_ids = {}

with open(blastout, "r") as b:
    for line in b:
        l = line.rstrip("\n").split("\t")
        l[2] = l[2].split(";")[0]   # sometimes blast returns multiple basically identical taxids separated by a semicolon, e.g. 1750886;1750895
        if l[2] not in uniq_subject_taxids.keys():
            uniq_subject_taxids[l[2]] = {'taxid_list': [], 'name_list': []}
        if l[0] not in uniq_query_ids.keys():
            uniq_query_ids[l[0]] = set([l[2]])
        else:
            uniq_query_ids[l[0]].add(l[2])
Entrez.email = "jason.vailionis@uconn.edu"
handle = Entrez.efetch(db="Taxonomy", id=list(uniq_subject_taxids.keys()), retmode="xml")
record = Entrez.parse(handle)
for i in record:
    taxid_list = []
    name_list = []
    current_taxid = i['TaxId']
    for j in i['LineageEx']:
        taxid_list.append(j['TaxId'])
        name_list.append(j['ScientificName'])
    taxid_list.append(current_taxid)
    name_list.append(i['ScientificName'])

    uniq_subject_taxids[current_taxid]['taxid_list'] = taxid_list
    uniq_subject_taxids[current_taxid]['name_list'] = name_list

# for each uniq query, find consensus lineage using set intersection and write output
# Current bug: since queries are based on the blastout file, seqs that didn't blast to anything won't show up in the output
with open(outfile, "w") as o:
    for i in uniq_query_ids.keys():
        if len(uniq_query_ids[i]) == 1:
            out = uniq_subject_taxids[uniq_query_ids[i].pop()]
            o.write("\t".join([i, out['name_list'][-1], out['taxid_list'][-1]]) + "\n")
        else:
            #currently this gets ruined for genes like 16S by "uncultured bacterium" stuff..
            #maybe just add an option to filter those and/or add some e_val filtering too
            setlist = [set(uniq_subject_taxids[j]['name_list']) for j in uniq_query_ids[i]]
            shared = setlist[0].intersection(*setlist[1:])
            reference_list = uniq_subject_taxids[uniq_query_ids[i].pop()]
            indices = [reference_list['name_list'].index(z) for z in shared]

            if len(shared) == 0:
                o.write("\t".join([i, "N/A", "N/A"]) + "\n")
            else:
                o.write("\t".join([i, reference_list['name_list'][max(indices)], reference_list['taxid_list'][max(indices)]]) + "\n")
