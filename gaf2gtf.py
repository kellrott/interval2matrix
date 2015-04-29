#!/usr/bin/env python

import sys
import os

COL_ENTRY_NUMBER = 0
COL_FEATURE_ID = 1
COL_FEATURE_TYPE = 2
COL_FEATURE_DB_SOURCE = 3
COL_FEATURE_DB_VERSION = 4
COL_FEATURE_DB_DATE = 5
COL_FEATURE_SEQ_FILENAME = 6
COL_COMPOSITE_ID = 7
COL_COMPOSITE_TYPE = 8
COL_COMPOSITE_DB_SOURCE = 9
COL_COMPOSITE_DB_VERSION = 10
COL_COMPOSITE_DB_DATE = 11
COL_ALIGNMENT_TYPE = 12
COL_FEATURE_COORDINATES = 13
COL_COMPOSITE_COORDINATES = 14
COL_GENE = 15
COL_GENE_LOCUS = 16
COL_FEATURE_ALIASES = 17
COL_FEATURE_INFO = 18


if __name__ == "__main__":
	if sys.argv[1] == "-":
		src_name = "stdin" 
		handle = sys.stdin
	else:		
		src_name = os.path.basename(sys.argv[1])
		handle = open(sys.argv[1])
	ohandle = sys.stdout
	
	for line in handle:
		if not line.startswith("#"):
			row = line.rstrip().split("\t")
			#only working on genes right now
			if row[COL_FEATURE_TYPE] == "gene":		
				
				#use the composite coordinates to create an array of exons, and 
				#calculate the frame windows for the different exons
				seq_str = row[COL_COMPOSITE_COORDINATES]
				chrome, locs, strand = seq_str.split(":")
				pos_list = []
				if strand == '+':
					cds_pos = 0
					for span in locs.split(","):
						start, stop = span.split("-")
						pos_list.append( (start, stop, cds_pos) )
						cds_pos += int(stop) + 1 - int(start)
				elif strand == '-':
					#for the reverse strand, we'll have to do this backwards, then 
					#flip the results
					cds_pos = 0
					tmp = locs.split(",")
					tmp.reverse()
					for span in tmp:
						start, stop = span.split("-")
						pos_list.append( (start, stop, cds_pos) )
						cds_pos += int(stop) + 1 - int(start)
					pos_list.reverse()
				
				#cycle through each exon and output a CDS and exon entry
				for i, pos in enumerate(pos_list):
					start, stop, cds_pos = pos
					
					kv = [
						["gene_id", row[COL_FEATURE_ID]],
						["transcript_id", row[COL_FEATURE_ID]],
						["exon_number", str(i+1)],
						["exon_id", "%s.%d" % (row[COL_FEATURE_ID], i+1)]						
					]
					
					#generate the key-value array to be put into attributes field
					kv_prep = ( "%s \"%s\"" % (a[0], a[1]) for a in kv )
					
					out = [
						chrome,
						src_name,
						"exon",
						start,
						stop,
						".",
						strand,
						".",
						"; ".join(kv_prep)
					]					
					ohandle.write("%s\n" % ("\t".join(out)))
					out[2] = "CDS"
					out[7] = "%d" % ( (3 - (cds_pos % 3)) % 3 )
					ohandle.write("%s\n" % ("\t".join(out)))
					
				#output the cds_start and cds_end fields using the gene locations
				seq_str = row[COL_GENE_LOCUS]
				if seq_str.count(";") == 0:
					chrome, span, strand = seq_str.split(":")
					gene_start, gene_stop = span.split("-")
					
					kv = [
						["gene_id", row[COL_FEATURE_ID]],
						["transcript_id", row[COL_FEATURE_ID]],
						["exon_number", str(i+1)],
						["exon_id", "%s.%d" % (row[COL_FEATURE_ID], i+1)]						
					]
						
					kv_prep = ( "%s \"%s\"" % (a[0], a[1]) for a in kv )
					kv_str = "; ".join(kv_prep)
					out = [
						chrome,
						src_name,
						"start_codon",
						"0",
						"0",
						".",
						strand,
						"0",
						kv_str
					]				
					if strand == '+':
						out[3] = gene_start
						out[4] = str(int(gene_start) + 2)
					else:
						out[3] = str(int(gene_stop) - 2)
						out[4] = gene_stop
					ohandle.write("%s\n" % ("\t".join(out)))
	
					out = [
						chrome,
						src_name,
						"stop_codon",
						"0",
						"0",
						".",
						strand,
						"0",
						kv_str
					]
					if strand == '+':
						out[3] = str(int(gene_stop) + 1 - 2)
						out[4] = str(int(gene_stop) + 1)
					else:
						out[3] = str(int(gene_start) - 2)
						out[4] = str(int(gene_start))
					ohandle.write("%s\n" % ("\t".join(out)))
	
					
