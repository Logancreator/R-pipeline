#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: cjy-
#Date: 2022-07-12
#Description:Filter reads from fastq.gz 

import gzip
import argparse
parser = argparse.ArgumentParser(description='filter reads from fastq.gz')
parser.add_argument('--fastq', '-q', dest='fastq', help='input a fastq.gz file')
parser.add_argument('--idlist', '-i', dest='idlist', help='input idlist.gz file')
parser.add_argument('--outfile', '-o', dest='outfile', help='input outfile name,end by gz')
args = parser.parse_args()

fastqdict = {}
fastqid = None
first = False

with gzip.open(args.fastq, 'rb') as fastq:
    for line in fastq:
        if line.startswith(b'@FP'):
            fastqid = line.split(b"/1")[0]+b"\n"
            # print(fastqid)
            # break
            fastqdict[fastqid] = b''
        else:
            fastqdict[fastqid] += line
        
outfile = gzip.open(args.outfile, 'wb')
    
with gzip.open(args.idlist, 'rb') as idfile:
    for line in idfile:
        readsid = line.split(b"\r")[0]+b"\n"
        # print(readsid)
        # break
        if not fastqdict.get(readsid, False):
            continue
        res =  readsid + fastqdict[readsid]
        outfile.write(res)

outfile.close()