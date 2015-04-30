#!/usr/bin/env python

import re
import sys
import vcf
import argparse
import math

from bx.intervals.intersection import Intersecter, Interval

reAttr = re.compile(r'([^ ]+) \"(.*)\"')


class GTFLine:
    def __init__(self, seqname, source, feature, start, end, score, strand, frame, attr):
        self.seqname = seqname
        self.source = source
        self.feature = feature
        self.start = long(start)
        self.end = long(end)
        try:
            self.score = float(score)
        except ValueError:
            self.score = None
        self.strand = strand
        try:
            self.frame = int(frame)
        except ValueError:
            self.frame = None
        self.attr = attr

    def __str__(self):
        return "\t".join( [
            self.seqname,
            self.source,
            self.feature,
            str(self.start),
            str(self.end),
            '.' if self.score is None else str(self.score) ,
            self.strand,
            '.' if self.frame is None else str(self.frame),
            str(self.attr)
            ])

    def __repr__(self):
        return "%s:%s:%s-%s" % (self.feature, self.seqname, self.start, self.end)

    def gene_id(self):
        return self.attr['gene_id']

    def get_type(self):
        return self.feature

    def __getitem__(self, name):
        return self.attr[name]



def gtfParse(handle):
    for line in handle:
        row = line.split("\t")
        attr_str = row[8]
        attr = {}
        for s in attr_str.split("; "):
            res = reAttr.search(s)
            if res:
                attr[res.group(1)] = res.group(2)
        yield GTFLine(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], attr)

class GTFGene:
    def __init__(self, gene_id):
        self.elements = {}
        self.gene_id = gene_id
        self.start = None
        self.end = None
        self.seqname = None

    def append(self, elem):
        if self.start is None or elem.start < self.start:
            self.start = elem.start
        if self.end is None or elem.end > self.end:
            self.end = elem.end
        self.seqname = elem.seqname
        t = elem.get_type()
        if t not in self.elements:
            self.elements[t] = []
        self.elements[t].append(elem)

    def gene_total_length(self):
        if 'start_codon' in self.elements and 'stop_codon' in self.elements:
            end = self.elements['stop_codon'][0].end
            start = self.elements['start_codon'][0].start
            return max(start,end) - min(start,end)

    def exon_total_length(self):
        total = 0
        if 'exon' in self.elements:
            for ex in self.elements['exon']:
                total += ex.end - ex.start
        return total



class GTFMap:
    def __init__(self):
        self.gene_map = {}
        self.chrom_map = {}

    def read(self, handle):
        for line in gtfParse(handle):
            g = line.gene_id()

            if line.seqname not in self.chrom_map:
                self.chrom_map[line.seqname] = []

            if g not in self.gene_map:
                gene = GTFGene(g)
                self.gene_map[g] = gene
                self.chrom_map[line.seqname].append(gene)

            self.gene_map[g].append(line)

    def __iter__(self):
        for k in self.gene_map:
            yield k

    def __getitem__(self, n):
        return self.gene_map[n]

    def keys(self):
        return self.gene_map.keys()

    def has_item(self, gene, item):
        if gene not in self.gene_map:
            return False
        if item not in self.gene_map[gene]:
            return False
        return True

    def has_chrom(self, chrome):
        return chrome in self.chrom_map

    def get_chrom(self, chrome):
        return self.chrom_map[chrome]

    def get_item_len(self, gene, item):
        total = 0
        for i in self.gene_map[gene][item]:
            total += i.end - i.start
        return total

    def get_item_span(self, gene, item):
        start = None
        end = None
        for i in self.gene_map[gene][item]:
            if start is None:
                start = i.start
            else:
                start = min(start, i.start)
            if end is None:
                end = i.end
            else:
                end = max(end, i.end)
        return start, end

BED_SEQ    = 0
BED_START  = 1
BED_STOP   = 2
BED_SAMPLE = 3
BED_VALUE  = 4

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-gtf', dest='gft', help="GTF definition of genes", required=True)
    parser.add_argument('-exon', dest='exon', action="store_true", default=False)
    parser.add_argument('--fix-chrome', action="store_true", default=False)
    parser.add_argument('--log2-ave', action="store_true", help="The values are in log2 space, average 2^value", default=False)
    parser.add_argument('-vcf', dest='vcf', default=None)
    parser.add_argument('-bed', dest='bed', default=None)
    args = parser.parse_args()

    #load the gene or exon intervals
    handle = open(args.gft)
    g = GTFMap()
    g.read(handle)

    i_map = {}
    if args.exon:
        for gene_name in g:
            gene = g[gene_name]
            if args.fix_chrome and gene.seqname.startswith("chr"):
                seqname = gene.seqname[3:]
            else:
                seqname = gene.seqname
            if seqname not in i_map:
                i_map[seqname] = Intersecter()
            for e in gene.elements['exon']:
                i_map[seqname].add_interval( Interval(e.start, e.end, value=gene) )
    else:
        for gene_name in g:
            gene = g[gene_name]
            if args.fix_chrome and gene.seqname.startswith("chr"):
                seqname = gene.seqname[3:]
            else:
                seqname = gene.seqname
            if seqname not in i_map:
                i_map[seqname] = Intersecter()
            i_map[seqname].add_interval( Interval(gene.start, gene.end, value=gene) )


    if args.vcf is not None:
        handle = open(args.vcf)

        missing = {}
        samples = {}
        v_read = vcf.Reader(handle)
        for rec in v_read:
            chrom = rec.CHROM
            if args.fix_chrome and rec.CHROM.startswith("chr"):
                chrom = rec.CHROM[3:]
            else:
                chrom = rec.CHROM
            if chrom in i_map:
                for hit in i_map[chrom].find(rec.POS, rec.POS+1):
                    gene = hit.value
                    gene_name = gene.gene_id
                    if gene_name not in mut_count:
                        mut_count[gene_name] = {}
                    for sam in rec.get_hets():
                        samples[sam.sample] = True
                        mut_count[gene_name][sam.sample] = mut_count[gene_name].get(sam.sample, 0 ) + 1
                    for sam in rec.get_hom_alts():
                        samples[sam.sample] = True
                        mut_count[gene_name][sam.sample] = mut_count[gene_name].get(sam.sample, 0 ) + 2
            else:
                if chrom not in missing:
                    sys.stderr.write("Missing Chrome %s\n" %(chrom))
                    missing[chrom] = True

        out = sys.stdout
        head = sorted(samples.keys())
        out.write("probe\t%s\n" % ("\t".join(head) ))
        for symbol in mut_count:
            cur = mut_count[symbol]
            row = []
            for c in head:
                row.append("%d" % (cur.get(c, 0)))
            out.write("%s\t%s\n" % (symbol, "\t".join(row)))

    if args.bed is not None:
        missing = {}
        samples = {}
        values = {}
        spans = {}
        with open(args.bed) as handle:
            for line in handle:
                row = line.split("\t")
                chrom = row[BED_SEQ]
                if args.fix_chrome and chrom.startswith('chr'):
                    chrom = chrom[3:]
                if chrom in i_map:
                    start = long(row[BED_START])
                    stop = long(row[BED_STOP])
                    for hit in i_map[chrom].find(start, stop+1):
                        gene = hit.value
                        gene_name = gene.gene_id
                        sample = row[BED_SAMPLE]
                        if gene_name not in values:
                            values[gene_name] = {}
                            spans[gene_name] = {}
                        samples[sample] = True
                        new_value = float(row[BED_VALUE])
                        span = float(min(stop, hit.end) - max(start, hit.start)) / float(hit.end - hit.start)
                        if span > 0.0:
                            if args.log2_ave:
                                new_value = math.pow(2,new_value)
                            values[gene_name][sample] = values[gene_name].get(sample, []) + [ new_value ]
                            spans[gene_name][sample] = spans[gene_name].get(sample, []) + [span]

                else:
                    if chrom not in missing:
                        sys.stderr.write("Missing Chrome %s\n" %(chrom))
                        missing[chrom] = True

        out = sys.stdout
        head = sorted(samples.keys())
        out.write("probe\t%s\n" % ("\t".join(head) ))
        for symbol in sorted(values.keys()):
            cur_values = values[symbol]
            cur_spans = spans[symbol]
            row = []
            for c in head:
                if c in cur_values:
                    cv = cur_values.get(c)
                    cs = cur_spans.get(c)
                    assert sum(cs) <= 1.0
                    value = 0.0
                    for v, s in zip( cv, cs ):
                        value += v * s
                    value /= sum(cs)
                    if args.log2_ave:
                        value = math.log(value,2)
                    row.append("%f" % (value))
                else:
                    row.append('NA')
            out.write("%s\t%s\n" % (symbol, "\t".join(row)))
