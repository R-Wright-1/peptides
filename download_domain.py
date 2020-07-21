#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 10:55:58 2020

@author: robynwright

Modified from the scripts at https://github.com/fischuu/Kraken_db_install_scripts/

This requires installation of Biopython: conda install -c conda-forge biopython
and also pandas: conda install pandas

Note that this seems to struggle with finding paths if they have any spaces. If this is the case, you will get something like the following error:
    FileNotFoundError: [Errno 2] No such file or directory: '/path/to/directory with spaces/archaea'
To fix, just make the directory already. In this case, I would need to make the archaea directory inside the one the script is running from.

The script will also create a log file (.txt) that you can search to check whether there were any sequences that didn't download/unzip/get the taxids changes to kraken ones.

"""

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='This script is to download all sequences from a certain domain that are present in NCBI refseq (ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/).\nThis requires the additional packages pandas and Biopython\nAs long as the script is able to download the assembly summary file, then it will create a log_file that tells you about whether each sequence was downloaded or not')
parser.add_argument('--domain', dest='domain', default='bacteria',
                    help='pick which domain to download, deafult is bacteria. Other options include: archaea, fungi, protozoa, viral, vertebrate_mammalian (or any others in NCBI refseq)')
parser.add_argument('--human', dest='human', default=False, 
                    help="add true here to download only human samples, default is False. Note that this only works if the domain is set to vertebrate_mammalian")
parser.add_argument('--ext', dest='ext', default='DNA', 
                    help="choose whether to download DNA or protein sequences, default is DNA. Valid options are DNA, dna, protein or Protein")
parser.add_argument('--complete', dest='complete', default=False, 
                    help="choose whether to only download complete genomes, or all genomes. Default is False, meaning all genomes are downloaded")


args = parser.parse_args()

wd = os.getcwd()+"/"
domain = args.domain
human = args.human
if args.ext == 'DNA' or args.ext == 'dna':
    ext1, ext2 = '_genomic', '.fna'
elif args.ext == 'Protein' or args.ext == 'protein':
    ext1, ext2 = '_protein', '.faa'
else:
    print('Unknown extension given. Valid options are DNA, dna, Protein or protein')
    quit()
complete = args.complete

#Create a directory
if not os.path.isdir(wd+domain):
  os.system("mkdir "+wd+domain)
os.chdir(wd+domain)

#Get the assembly file
try:
  if not os.path.exists("assembly_summary.txt"):
    os.system("wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/"+domain+"/assembly_summary.txt")
except:
  print("Unable to download assembly_summary.txt for "+domain)
  quit()

#Parse the file
try:
  f = pd.read_csv("assembly_summary.txt", header=1, index_col=0, sep="\t")
except:
  print("Unable to open assembly_summary.txt for "+domain)
  quit()

#Parse each of the lines in the file, getting the sequences and changing the IDs of all sequences to match kraken tax ids
#Output what happens for each line of the file - whether it works, and also which point it fails at
log_file = []
for l in list(f.index.values):
  line = f.loc[l, :]
  if human:
      if line["organism_name"] != "Homo sapiens":
          continue
  if complete:
      if not human:
          if not line['assembly_level'] == 'Complete Genome':
              continue
  ftppath = line['ftp_path']
  aname = ftppath.split('/')[-1]
  if os.path.exists(aname+ext1+".tax."+ext2):
    log_file.append("Already got this sequence: "+aname+ext1+".tax."+ext2+"\n")
    print("Already had this sequence: "+aname+ext1+".tax."+ext2)
    continue
  else:
      fullpath = ftppath+"/"+aname+ext1+ext2+".gz"
      if not os.path.exists(aname+ext1+ext2):
          try:
            os.system("wget -q "+fullpath)
          except:
            log_file.append("Didn't download "+aname+ext1+ext2+"\n")
            print("Didn't download "+aname)
            continue
          try:
            os.system("gunzip "+aname+ext1+ext2+".gz")
          except:
            log_file.append("Didn't unzip "+aname+ext1+ext2+".gz\n")
            print("Didn't unzip "+aname+ext1+ext2+".gz")
            continue
      else:
          log_file.append("Already had the unzipped file: "+aname+ext1+ext2+".gz\n")
          print("Already had the unzipped file: "+aname+ext1+ext2+".gz")
      try:
          taxid = line['taxid']  
          new_records = []
          records = SeqIO.parse(aname+ext1+ext2, "fasta")
          for record in records:
              newid = record.id+"|kraken:taxid|"+str(taxid)
              newseq = SeqRecord(record.seq, id=newid, description=record.description)
              new_records.append(newseq)
          SeqIO.write(new_records, aname+ext1+".tax"+ext2, "fasta")
          os.system("rm "+aname+ext1+ext2)
          log_file.append("Got this sequence with kraken taxid: "+aname+ext1+".tax"+ext2+"\n")
          print("Got this sequence with kraken taxid: "+aname+ext1+".tax"+ext2)
      except:
          log_file.append("Didn't manage to change the taxids for this file: "+aname+ext1+ext2+"\n")
          print("Didn't manage to change the taxids for this file: "+aname+ext1+ext2)
os.chdir(wd)
with open("logfile_"+domain+".txt", 'w') as f:
    for row in log_file:
        f.write(row)