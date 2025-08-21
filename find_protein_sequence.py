import re
import sys

reference_rna=""
codonmap={}
domains={}

del_mut_pattern=re.compile("(Del):(\\d+)(?:-(\\d+))?")
ins_mut_pattern=re.compile("(Ins):(\\d+):([ACGT]+)")
nuc_mut_pattern=re.compile("Nuc:(\\d+)([ACGT])")

def load_codon_map():
  f=open("data/dna_codon.csv","r")
  for line in f.readlines():
    data=line[:-1].split(",")
    aa=data[1]
    codon=data[0]
    codonmap[codon]=aa
  f.close()

def load_reference_strain_rna():
  global reference_rna
  f=open("data/reference_sequence.dat","r")
  rnadata=""
  for l in f.readlines():
    rnadata="%s%s"%(rnadata,l[:-1])
  f.close()
  reference_rna=rnadata

def load_domain_file(fname):
  global domains
  f=open(f"data/{fname}","r")
  lin=f.readline()
  while (len(lin)>0):
    parts=lin[:-1].split(",")
    domains[parts[0]]=(int(parts[1]),int(parts[2]))
    lin=f.readline()
  f.close()

def load_domains():
  load_domain_file("domains.csv")
  load_domain_file("gisaid_domains.csv")
  load_domain_file("non_nextclade_domains.csv")
  domains['ORF1ab']=(domains['ORF1a'][0],domains['ORF1b'][1])

def get_mutations_path(lineage):
  path0=re.sub(r"([A-Z])",r"/\1",lineage)
  path=re.sub(r"[.]",r"/",path0)
  return f"../sars-cov-2-lineage-dominant-mutations{path}/{lineage}-muts.txt"

def convert_to_protein(rna):
  result=""
  for p in range(0,len(rna),3):
    result="%s%s"%(result,codonmap[rna[p:p+3]])
  return result

def load_lineage_nuc_mutations(lineage,positions,frame):
  mutations_path=get_mutations_path(lineage)
  f=open(mutations_path,"r")
  indels={}
  lineage_rna=reference_rna
  output_start=positions[0]-1
  output_end=positions[1]-1
  if (frame=='NSP12') or (frame=='ORF1ab'):
    # For frames that span the ribosome slip from ORF1a to ORF1b,
    # add a deletion for the slip effect
    indels[13466]=('Del',13466,13467)
  for l in f.readlines():
    del_mut_match=del_mut_pattern.match(l[:-1])
    if del_mut_match:
      start_pos=int(del_mut_match.group(2))
      indels[start_pos]=(del_mut_match.group(1),del_mut_match.group(2),del_mut_match.group(3))
    ins_mut_match=ins_mut_pattern.match(l[:-1])
    if ins_mut_match:
      start_pos=int(ins_mut_match.group(2))
      indels[start_pos]=(ins_mut_match.group(1),ins_mut_match.group(2),ins_mut_match.group(3))
    nuc_mut_match=nuc_mut_pattern.match(l[:-1])
    if nuc_mut_match:
      mutPos=int(nuc_mut_match.group(1))
      lineage_rna=lineage_rna[:mutPos-1]+nuc_mut_match.group(2)+lineage_rna[mutPos:]

  f.close()

  indel_pos_list=list(indels.keys())
  indel_pos_list.sort()
  indel_pos_list.reverse()

  for p in indel_pos_list:
     delta=0
     indel=indels[p]
     if (indel[0]=='Del'):
       del_start=int(indel[1])
       if indel[2] is not None:
         del_end=int(indel[2])
       else:
         del_end=del_start
       delta=(del_start-del_end)-1
       lineage_rna=lineage_rna[:(del_start-1)]+lineage_rna[del_end:]
     if (indel[0]=='Ins'):
       ins_pos=int(indel[1])
       ins_seq=indel[2]
       delta=len(ins_seq)
       lineage_rna=lineage_rna[:ins_pos]+ins_seq+lineage_rna[ins_pos:]
     if (p<=output_end):
       output_end=output_end+delta
     if (p<=output_start):
       output_start=output_start+delta

  print(convert_to_protein(lineage_rna[output_start:output_end+1]))

load_reference_strain_rna()
load_codon_map()
load_domains()


if len(sys.argv)>2:
   load_lineage_nuc_mutations(sys.argv[2],domains[sys.argv[1]],sys.argv[1])
else:
   usage="""
Usage:

   python3 find_protein_sequence.py [protein] [lineage]
   
Available protein specifications are:
   """
   proteins=list(domains.keys())
   proteins.sort()
   print(usage)
   for p in proteins:
      print(f"   {p}")
   warranty="""

Use this program at your own risk as it has no warranty or guarantee of correctness.
   """
   print(warranty)


