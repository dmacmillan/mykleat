import time
import random,string
import copy
import argparse
import logging
import os
import sys
import re
import subprocess
# External modules below
import pysam

parser = argparse.ArgumentParser(description='this program tries to find polya cleavage sites through short-read assembly.it is expected that contigs are aligned to contigs, and reads aligned to contigs. these 2 alignment steps can be performed by trans-abyss. the aligners used are gmap for contig-genome and bwa-sw for read-contig alignments. annotations files for ensembl, knowngenes, refseq, and aceview are downloaded from ucsc. est data(optional) are also downloaded from ucsc. the analysis can be composed of 2 phases: 1. contig-centric phase - cleavage sites per contig are captured 2. coordinate-centric phase - contigs capturing the same cleavage site are consolidated into 1 report where expression/evidence-related data are summed. customized filtering based on evidence data can be performed.')
parser.add_argument('c2g', metavar='<contig-to-genome>', help='The contig-to-genome alignment file in bam format.')
parser.add_argument('contigs', metavar='<contigs>', help='The file containing all contigs in fasta format.')
parser.add_argument('ref_genome', metavar='<reference_genome>', help='The path to the reference genome to use.')
parser.add_argument('annot', metavar='<annotations>', help='The annotations file to use with the reference in gtf format.')
parser.add_argument('r2c', metavar='<reads-to-contigs>', help='The contigs-to-genome alignment file.')
parser.add_argument('out', metavar='<output-file>', help='The file to output results.')
parser.add_argument('-t', '--trim', dest='trim_reads', help='Trim bases of quality <= this value. default is 3.', type=int, default=3)
parser.add_argument('-ss', '--strand_specific', action='store_true', help='Enable if reads are strand specific.')
parser.add_argument('--use_tmp', help='Use tmp space', action='store_true', default=False)
parser.add_argument('-e', '--overlap_est', help='Overlap expressed sequence tags', action='store_true', default=False)
parser.add_argument('--min_at', help='Minimum number of a|t bases in tail. default is 4.', type=int, default=4)
parser.add_argument('--max_diff', help='Maximum rate of x to y bases allowed in tail. where x are non a|t bases and y are a|t bases. default: 1 5.', nargs=2, default=[1,5])
parser.add_argument('--max_diff_link', help='Maximum number of non A|T bases in entire link read. Default is 2.', type=int, default=2)
parser.add_argument('--min_bridge_size', help='Minimum size of bridge. Default is 1.', type=int, default=1)
parser.add_argument('--max_dist', help='The maximum distance that a putative cleavage site may be from an annotated end. Default is 5000', default=5000)
parser.add_argument('-k', '--track', metavar=('[name]','[description]'), help='Name and description of BED graph track to output.', nargs=2)
parser.add_argument('--rgb', help='RGB value of BED graph. Default is 0,0,255', default='0,0,255')
parser.add_argument('-c', help='Specify a contig/s to look at.', nargs='+')
parser.add_argument('--link', action='store_true', help='Enable searching for cleavage site link evidence. This will substantially increase runtime.')

args = parser.parse_args()
#logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('polyA_logger')

#fh = logging.FileHandler('info.log')
#fh.setLevel(logging.INFO)
#logger.addHandler(fh)

# Reference genome sequence (refseq)
logger.debug("Loading reference genome via pysam 0.8.1")
refseq = pysam.FastaFile(args.ref_genome)
logger.debug("Reference genome successfully loaded!")

# Contigs to genome alignments (aligns)
logger.debug("Checking contig to genome alignment is in BAM format. Splitting file by '.': %s", args.c2g.split('.'))
if (args.c2g.split('.')[-1].lower() == 'bam'):
    logger.debug("Contig to genome alignment file is in BAM format, loading BAM file...")
    aligns = pysam.AlignmentFile(args.c2g, "rb")
    logger.debug("Contig to genome alignment file loaded successfully")
elif (args.c2g.split('.')[-1].lower() == 'sam'):
    aligns = pysam.AlignmentFile(args.c2g, "r")
else:
    sys.exit("Unrecognized format for <aln>, must be BAM format file. Exiting.")

# Contig sequences
contigs = pysam.FastaFile(args.contigs)

# Features (features)
features = pysam.TabixFile(args.annot, parser=pysam.asGTF())

# Feature dictionary
feature_dict = {}
f = pysam.tabix_iterator(open(args.annot),parser=pysam.asGTF())
for c in f:
    chrom = c.contig
    tid = c.asDict()['transcript_id']
    if chrom not in feature_dict:
        feature_dict[chrom] = {}
    if tid not in feature_dict[chrom]:
        feature_dict[chrom][tid] = {'feats':[c], 'cstart':None, 'cend':None, 
                                    'i':0, 'start_codon': None, 'stop_codon': None,
                                    'maxfeat': 0, 'utr3': [], 'utr5': [], 'strand': None,
                                    'cleavage_sites': []}
    else:
        feature_dict[chrom][tid]['feats'].append(c)
    if (not feature_dict[chrom][tid]['strand']):
        feature_dict[chrom][tid]['strand'] = c.strand
    if (feature_dict[chrom][tid]['feats'][feature_dict[chrom][tid]['i']].end > feature_dict[chrom][tid]['feats'][feature_dict[chrom][tid]['maxfeat']].end):
        feature_dict[chrom][tid]['maxfeat'] = feature_dict[chrom][tid]['i']
    if (c.feature == 'CDS') and (not feature_dict[chrom][tid]['cstart']):
        feature_dict[chrom][tid]['cstart'] = feature_dict[chrom][tid]['i']
    elif (c.feature == 'CDS'):
        feature_dict[chrom][tid]['cend'] = feature_dict[chrom][tid]['i']
    if (c.feature == 'start_codon'):
        feature_dict[chrom][tid]['start_codon'] = feature_dict[chrom][tid]['i']
    elif (c.feature == 'stop_codon'):
        feature_dict[chrom][tid]['stop_codon'] = feature_dict[chrom][tid]['i']
    feature_dict[chrom][tid]['i'] += 1

# Go through the feature_dict and try to determine 3utrs and 5utrs
for chrom in feature_dict:
    for tid in feature_dict[chrom]:
        current = feature_dict[chrom][tid]
        if (current['feats'][0].strand == '+'):
            if (current['cend']) and (current['feats'][current['maxfeat']].end > (current['feats'][current['cend']].end+3)):
                # Plus 3 to first coordinate to account for stop codon
                current['utr3'] = [current['feats'][current['cend']].end+3,current['feats'][current['maxfeat']].end]
        else:
            if (current['cstart']) and (current['feats'][0].start < (current['feats'][current['cstart']].start-3)):
                # Minus 3 to last coordinate to account for stop codon
                current['utr3'] = [current['feats'][0].start,current['feats'][current['cstart']].start-3]

# Reads to contigs alignment (r2c)
r2c = pysam.AlignmentFile(args.r2c, "rb")

def int_to_base_qual(qual_int, offset):
    """Converts integer to base quality base"""
    if offset == 64 or offset == 33:
        return chr(qual_int + offset)

# Poor quals (poor_quals)
poor_quals = ''
if args.trim_reads:
    for i in range(args.trim_reads + 1):
        poor_qual = int_to_base_qual(i, 33)
        if poor_qual is not None:
            poor_quals += poor_qual

output_fields=['gene','transcript','transcript_strand','coding','contig','chromosome','cleavage_site','within_UTR','distance_from_annotated_site','ESTs','length_of_tail_in_contig','number_of_tail_reads','number_of_bridge_reads','max_bridge_read_tail_length','bridge_read_identities','tail+bridge_reads','number_of_link_pairs','max_link_pair_length','link_pair_identities','hexamer_loc+id','3UTR_start_end']
if args.c:
    args.out += '.'+('_').join(args.c)
prefix = os.path.splitext(args.out)[0]
basedir = os.path.dirname(args.out)
if not os.path.exists(basedir):
    os.makedirs(basedir)
#out_link_pairs = open(os.path.join(basedir,prefix+'-link.fa'), 'a')
reads_to_check = open(os.path.join(basedir,'.reads_to_check'), 'w')

# Filters (filters)
global_filters = {}
global_filters['min_at'] = int(args.min_at)
global_filters['max_diff'] = [int(x) for x in args.max_diff]
global_filters['max_diff_link'] = args.max_diff_link
global_filters['min_bridge_size'] = args.min_bridge_size
global_filters['feature_dict'] = feature_dict
global_filters['hexamer_colours'] = ["255,0,0", "255,100,100", "255,150,150", "255,200,200",
                                     "0,255,0", "100,255,100", "150,255,150", "200,255,200",
                                     "0,0,255", "100,100,255", "150,150,255", "200,200,255",
                                     "255,0,255", "255,100,255", "255,150,255", "255,200,255"]
global_filters['binding_sites'] = ['AATAAA','ATTAAA','AGTAAA','TATAAA',
                                   'CATAAA','GATAAA','AATATA','AATACA',
                                   'AATAGA','AAAAAG','ACTAAA','AAGAAA',
                                   'AATGAA','TTTAAA','AAAACA','GGGGCT']
global_filters['all_results'] = {}

def ucsc_chroms(genome):
    """Extracts conversion of UCSC chromosome names
    eg. hg19"""
    package_dir = "/".join(os.path.abspath(__file__).split("/")[:-2])
    conversion_file = package_dir + '/annotations/' + genome + "/ucsc_chr.txt"

    conversions = {}
    if os.path.exists(conversion_file):
        for line in open(conversion_file, 'r'):
            chr_from, chr_to = line.rstrip('\n').split()
            conversions[chr_from] = chr_to
            
    return conversions

def compare_chr(chr1, chr2):
    """For sorting chromosome names ignoring 'chr'"""
    if chr1[:3].lower() == 'chr':
        chr1 = chr1[3:]
    if chr2[:3].lower() == 'chr':
        chr2 = chr2[3:]
    
    if re.match('^\d+$', chr1) and not re.match('^\d+$', chr2):
        return -1
    
    elif not re.match('^\d+$', chr1) and re.match('^\d+$', chr2):
        return 1
    
    else:
        if re.match('^\d+$', chr1) and re.match('^\d+$', chr2):
            chr1 = int(chr1)
            chr2 = int(chr2)
            
        if chr1 < chr2:
            return -1
        elif chr1 > chr2:
            return 1
        else:
            return 0

def get_coding_type(transcript):
    """Returns transcript type: CODING/NONCODING/NA
    CODING when cdsStart != cdsEnd
    """
    #if transcript['cstart'] == None or transcript['cend'] == None:
    #    return 'unknown'
    #elif transcript['cstart'] and transcript['cend'] and (transcript['cstart'] != transcript['cend']):
    if transcript['cstart'] and transcript['cend']:
        return 'yes'
    else:
        return 'no'

def cantorPairing(a,b):
    return (0.5*(a+b)*(a+b+1))+b

def group_and_filter(lines_result, out_file, filters=None, make_track=None, rgb='0,0,0'):
    #print 'grouping and filtering'
    #print 'lines_result: {}'.format(lines_result)
    global output_fields
    hexamer_colours = ["255,0,0", "255,100,100", "255,150,150", "255,200,200",
                       "0,255,0", "100,255,100", "150,255,150", "200,255,200",
                       "0,0,255", "100,100,255", "150,150,255", "200,200,255",
                       "255,0,255", "255,100,255", "255,150,255", "255,200,255"]
    binding_sites = ['AATAAA','ATTAAA','AGTAAA','TATAAA',
                     'CATAAA','GATAAA','AATATA','AATACA',
                     'AATAGA','AAAAAG','ACTAAA','AAGAAA',
                     'AATGAA','TTTAAA','AAAACA','GGGGCT']
    """Consolidates contig-centric results into coordinate-centric results
    
    path = directory of tab-delimited results files
    """
    groups = {}
    
    lines_result = [x.split('\t') for x in lines_result.splitlines()]
    for cols in lines_result:
        chrom, cleavage_site = cols[5], cols[6]
        if not groups.has_key(chrom):
            groups[chrom] = {}
        if not groups[chrom].has_key(cleavage_site):
            groups[chrom][cleavage_site] = []
        groups[chrom][cleavage_site].append(cols)
            
    stats = {'gene': {},
             'transcript': {'coding': {}, 'noncoding':{}, 'unknown': {}},
             'screened_out': 0,
             'cleavage_site': 0,
             'with_tail': 0,
             'tail_and_bridge_and_link': 0,
             'tail_and_bridge': 0,
             'tail_and_link': 0,
             'bridge_and_link': 0,
             'just_tail': 0,
             'just_bridge': 0,
             'just_link': 0,
             }
            
    out = open(out_file, 'w')
    out.write('%s\n' % '\t'.join(output_fields))
    
    track_plus = []
    track_minus = []
    hexamers = []
    utrs = []
    uniqueutrs = set()
    uniquehexamers = {}
    for chrom in sorted(groups.keys(), cmp=compare_chr):
        for cleavage in sorted(groups[chrom].keys()):
            results = groups[chrom][cleavage]
            # if more than one contig reports same cleavage site, add up the support numbers
            if len(results) > 1:
                result = merge_results(results)
            else:
                result = results[0]
            if (result[10] != '-') and (result[12] != '-') and (result[16] != '-'):
                if (int(result[10]) == 0) and (int(result[12]) == 0) and (int(result[16]) == 0):
                    continue
                                    
            if result[12] != '-':
                if (filters) and ('min_bridge_size' in filters) and ((int(result[12]) > 0) and (int(result[13]) < filters['min_bridge_size'])):
                    continue
                
            out.write('%s\n' % '\t'.join(result))
            if make_track is not None:
                if result[2] == '+':
                    track_plus.append(show_expression(result))
                    for thing in show_hexamer(result,binding_sites,hexamer_colours):
                        
                        hexamers.append(thing)
                else:
                    track_minus.append(show_expression(result))
                    for thing in show_hexamer(result,binding_sites,hexamer_colours):
                        hexamers.append(thing)

                utr = show_utr(result)
                if utr:
                    #print 'utr: {}'.format(utr)
                    a,b = [int(x) for x in utr.split('\t')[1:3]]
                    pairing = cantorPairing(a,b)
                    if pairing not in uniqueutrs:
                        utrs.append(utr)
                        uniqueutrs.add(pairing)
            
            # stats
            update_stats(stats, result)
                
    out.close()
        
    # prefix used for track and stats files
    prefix = os.path.splitext(out_file)[0]
    # output track
    #randstr = ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(10))
    if make_track is not None:
        track_file = prefix + '.bg'
        output_track2(prefix+'.+.bg', make_track[0]+'.+', make_track[1], rgb, track_plus)
        output_track2(prefix+'.-.bg', make_track[0]+'.-', make_track[1], rgb, track_minus)
        output_hexamer(prefix + '.HEXAMERS.bed', hexamers)
        output_track_utr(prefix+'.3UTR.bed', make_track[0]+'.3UTRs', make_track[1], rgb, utrs)
        
    # output stats file
    stats_file = prefix + '.stats'
    output_stats(stats, stats_file)

def prepare_track_header(name, desc, rgb):
    """Creates header for track"""
    return 'track type=bedGraph name="%s" description="%s" visibility=full color=%s' % (name, desc, rgb)

def output_hexamer(out_file, track=[]):
    out = open(out_file, 'w')
    out.write('track name="hexamer_track" description="Track containing all CPSF hexamer binding sites" visibility=2 itemRgb="On"\n')
    for line in track:
        out.write('%s\n' % line)
    out.close()

def output_track_utr(out_file, name, desc, rgb, track=[]):
    out = open(out_file, 'w')
    out.write('track name="{}" description="3\'UTR" visibility=full itemRgb="On"\n'.format(name, rgb))
    for line in track:
        out.write('{}\n'.format(line))

def output_track2(out_file, name, desc, rgb, track=[]):
    out = open(out_file, 'w')
    out.write('track type=bedGraph name="{}" description="{}" visibility=full color={}\n'.format(name, desc, rgb))
    for line in track:
        out.write('{}\n'.format(line))

def merge_results(results):        
    """Merges results from different contigs of same cleavage site into single result"""
    # join fields: contig, bridge_name, link_name
    join = [4, 14, 18]              
    # add fields: num_tail_reads, num_bridge_reads, tail+bridge, num_link_pairs
    add = [11, 12, 15, 16]
    # max fields: tail_len, bridge_len, link_len
    biggest = [10, 13, 17]
    
    merged = []
    for i in range(len(results[0])):
        if i in add:
            data = [int(r[i]) for r in results if r[i].isdigit()]
            if data:
                merged.append(str(sum(data)))
            else:
                merged.append('-')              
        elif i in biggest:
            data = [int(r[i]) for r in results if r[i].isdigit()]
            if data:
                merged.append(str(max(data)))
            else:
                merged.append('-')              
        elif i in join:
            data = [r[i] for r in results if r[i] != '-']
            
            # if there is data other than '-', then list all items
            if data:
                merged.append(','.join([r[i] for r in results if r[i] != '-']))
            else:
                merged.append('-')                  
        else:
            merged.append(results[0][i])

    return merged

def update_stats(stats, result):
    """Updates summary stats with result
    
    stats = dictionary of final stats results
    result = list of values of each output line
    """
    has_tail = has_bridge = has_link = False
    if result[10].isdigit() and int(result[10]) > 0:
        has_tail = True
    if result[12].isdigit() and int(result[12]) > 0:
        has_bridge = True
    if result[16].isdigit() and int(result[16]) > 0:
        has_link = True
        
    if not (has_tail or has_bridge or has_link):
        return
        
    stats['cleavage_site'] += 1
    
    if not stats['gene'].has_key(result[0]):
        stats['gene'][result[0]] = 0
    stats['gene'][result[0]] += 1
    
    if result[3] == 'yes':
        stats['transcript']['coding'][result[1]] = True
    elif result[3] == 'no':
        stats['transcript']['noncoding'][result[1]] = True
    else:
        stats['transcript']['unknown'][result[1]] = True
                        
    if has_tail and has_bridge and has_link:
        stats['tail_and_bridge_and_link'] += 1
    elif has_tail and has_bridge:
        stats['tail_and_bridge'] += 1
    elif has_bridge and has_link:
        stats['bridge_and_link'] += 1
    elif has_tail and has_link:
        stats['tail_and_link'] += 1
    elif has_tail:
        stats['just_tail'] += 1
    elif has_bridge:
        stats['just_bridge'] += 1
    elif has_link:
        stats['just_link'] += 1

def output_stats(stats, out_file):
    """Outputs stats data into output file"""
    out = open(out_file, 'w')
    out.write('total cleavage sites: %d\n' % stats['cleavage_site'])
    out.write('cleavage sites with tail, bridge, link support: %d\n' % stats['tail_and_bridge_and_link'])
    out.write('cleavage sites with tail, bridge support: %d\n' % stats['tail_and_bridge'])
    out.write('cleavage sites with tail, link support: %d\n' % stats['tail_and_link'])
    out.write('cleavage sites with bridge, link support: %d\n' % stats['bridge_and_link'])
    out.write('cleavage sites with only tail support: %d\n' % stats['just_tail'])
    out.write('cleavage sites with only bridge support: %d\n' % stats['just_bridge'])
    out.write('cleavage sites with only link support: %d\n' % stats['just_link'])
    out.write('total genes: %d\n' % len(stats['gene'].keys()))
    try:
        out.write('average cleavage sites per gene: %.1f\n' % (float(stats['cleavage_site'])/len(stats['gene'].keys())))
    except ZeroDivisionError:
        out.write('average cleavage sites per gene: error')
    out.write('total transcripts: %d\n' % (len(stats['transcript']['coding'].keys()) + len(stats['transcript']['noncoding'].keys())))
    out.write('total coding transcripts: %d\n' % len(stats['transcript']['coding'].keys()))
    out.write('total noncoding transcripts: %d\n' % len(stats['transcript']['noncoding'].keys()))
    out.close()


def show_expression(result):
    """Creates bed-graph line depicting expression of cleavage site"""
    return '%s\t%s\t%s\t%s' % (result[5], int(result[6]) - 1, result[6], result[15])

def show_hexamer(result, binding_sites, rbgs):
    r = []
    sites = result[19].split(';')
    sites = [x.split(':') for x in sites]
    for site in reversed(sites):
        if len(site) < 2:
            break
        if result[2] == '+':
            r.append('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(result[5], site[0], int(site[0])+6, binding_sites[int(site[1])-1], int(site[1])*62.5, result[2], site[0], int(site[0])+6, rbgs[int(site[1])-1]))
        else:
            r.append('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(result[5], int(site[0])-6, site[0], binding_sites[int(site[1])-1], int(site[1])*62.5, result[2], int(site[0])-6, site[0], rbgs[int(site[1])-1]))
    return r

def show_utr(result):
    if result[20] == '-':
        return None
    start, end = result[20].split('-')
    if (start == 'None') or (end == 'None'):
        return None
    rgb='255,0,0'
    if result[20][0] == 'N':
        start = start[1:]
        rgb='0,255,0'
    #if result[2] == '-':
    #    return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(result[5], end, start, result[4], 0, result[2], end, start, rgb)
    #else:
    #    return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(result[5], start, end, result[4], 0, result[2], start, end, rgb)
    return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(result[5], start, end, result[4], 0, result[2], start, end, rgb)

# chrom_proper
chrom_proper = ucsc_chroms(args.ref_genome)

def tpos_to_qpos(a,tpos):
    blocks = a['qblocks']
    tblocks = a['align'].blocks
    for i in xrange(len(tblocks)):
        tb0 = tblocks[i][0]
        tb1 = tblocks[i][1]
        if ((tpos >= tb0) and (tpos <= tb1)) or ((tpos <= tb0) and (tpos >= tb1)):
            block = i
            break
    qpos = None
    try:
        if (a['strand'] == '+'):
            qpos = blocks[block][0]-1 + tpos - tblocks[block][0] # added the +1
        else:
            qpos = blocks[block][1] - (tpos - tblocks[block][1])
    except NameError:
        pass
    return qpos

def qpos_to_tpos(a, qpos):
    blocks = align.blocks
    if align.is_reverse:
        strand = '-'
    else:
        strand = '+'
    for i in xrange(len(a['qblocks'])):
        qb0 = a['qblocks'][i][0]
        qb1 = a['qblocks'][i][1]
        if ((qpos >= qb0) and (qpos <= qb1)) or ((qpos <= qb0) and (qpos >= qb1)):
            block = i
            break
    tpos = None
    try:
        if (a['strand'] == '+'):
            tpos = blocks[block][0]+1 + qpos - a['qblocks'][block][0] # added the +1
        else:
            tpos = blocks[block][1] - (qpos - a['qblocks'][block][1])
    except NameError:
        pass
    return tpos

def check_freq(seq):
    """Returns frequency of each base in given sequence"""
    freq = {}
    for nt in seq:
        if not freq.has_key(nt.upper()):
            freq[nt.upper()] = 0
        freq[nt.upper()] += 1
    return freq

def is_polyA_tail(seq, expected_base, min_len, max_nonAT_allowed):
    """Determines if sequence can be possible tail
    
    min_len = minimum number of expected base in seq
    max_nonAT_allowed(N,M) = N base(s) other than expected base allowed
                             per stretch of M bases
    """
    #print '{}\t{}\t{}\t{}'.format(seq,expected_base,min_len,max_nonAT_allowed)
    result = True

    if seq is None or seq == '' or expected_base is None or not expected_base in seq:
        return False
        
    # minimum number of As or Ts
    freq = check_freq(seq)
    #print freq
    if freq[expected_base] < min_len:
        return False

    if (max_nonAT_allowed[0] == 0):
        if (freq[expected_base] < len(seq)):
            return False
        
    non_base = sum([freq[b] for b in freq if b != expected_base])
    if (non_base == 0):
        return True
    ratio = float(freq[expected_base])/non_base
    #print 'ratio: {}'.format(ratio)
    #print 'max_nonAT: {}'.format(max_nonAT_allowed[1]/max_nonAT_allowed[0])
    if (ratio < max_nonAT_allowed[1]/max_nonAT_allowed[0]):
        return False

#    for i in range(0, len(seq), max_nonAT_allowed[1]):
#        subseq = seq[i:i+10]
#    
#        freq = check_freq(subseq)
#        non_expected_freq = 0
#        for base, f in freq.iteritems():
#            if base.upper() != expected_base.upper():
#                non_expected_freq += f
#        
#        if non_expected_freq > max_nonAT_allowed[0]:
#            print 'returning false because non_expected_feq > max_nonAT_allowed[0]'
#            result = False
            
    return result

def getQueryLenFromCigar(cigar):
    lens = []
    # if op is not 'D'(deletion), 'N'(skipped region), 'P'(padding)
    query_lens = [int(op[1]) for op in cigar if op[0] != 2 and op[0] != 3 and op[0] != 6]
    return sum(query_lens)

def cigarToBlocks(cigar, tstart, strand):
    query_len = getQueryLenFromCigar(cigar)
    qstart = 1 if strand == '+' else query_len
    tblocks = []
    qblocks = []          
    for i in range(len(cigar)):
        op, length = cigar[i]
        # 'D' (deletion)
        if op == 2 and i == 0:
            return None, None
        # 'S' or 'H' (clips)
        if op == 4 or op == 5:
            if i == 0:
                qstart = qstart + length if strand == '+' else qstart - length
                continue
        tblock = None
        qblock = None
        if not tblocks and op != 0:
            return None, None
        # match
        if op == 0:
            tend = tstart + length - 1
            qend = qstart + length - 1 if strand == '+' else qstart - length + 1
            tblock = [tstart, tend]
            qblock = [qstart, qend]
        # intron ('N'), skipped reference or deletion in reference ('D')
        elif op == 2 or op == 3:
            #intron
            if op == 3 and length < 3:
                continue
            qend = qend
            tend = tstart + length - 1 
        # insertion ('I') to reference
        elif op == 1:
            tend = tend
            qend = qstart + length - 1 if strand == '+' else qstart - length + 1
        if tblock:
            tblocks.append(tblock)
        if qblock:
            qblocks.append(qblock)
        tstart = tend + 1
        qstart = qend + 1 if strand == '+' else qend - 1
    return tblocks, qblocks

def revComp(seq):
    seq = seq.upper()
    ndic = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    revcomp = ''
    for nuc in seq:
        revcomp += ndic[nuc]
    return revcomp[::-1]

def getQueryBlocks(align):
    query_blocks = []
    start = align.qstart
    for block in align.blocks:
        diff = block[1]-block[0]
        query_blocks.append([start,start+diff])
        start = start+diff
    return query_blocks

def inferStrand(align, feature):
    if feature:
        return feature.strand

def find_link_pairs(a, utr3, homo_len=20, max_mismatch=0):
    """Finds reads pairs where one mate is mapped to contig and the other mate (mapped elsewhere or unmapped)
    is a potential polyA tail
    
    The mate that is mapped to the contig (anchor) should be pointing outwards to the edge of the contig.
    The mate that is potential polyA tail can be unmapped or mapped to another contig.  If the mate is unmapped,
    it's expected to be stored under the same contig.  Unfortunately this is not the behaviour of BWA-SW so all the 
    unmapped mates are ignored for BWA-SW alignments.
    
    A matching transcript is provided.  The only place this is used is in the check of whether the vicinity of the 
    cleavage site has a homopolyer run.  The transcript strand is used for deciding whether the upstream or downstream
    region of the cleavage site should be checked.  If a polyT or polyA is in the neighborhood (200bp), then the case
    won't be further explored.
    """
    # determine if 3'UTR is first or last query block
    anchor_read_strand = None
    if utr3['clipped_pos'] == 'start':
        anchor_read_strand = '-'
    elif utr3['clipped_pos'] == 'end':
        anchor_read_strand = '+'
    
    if anchor_read_strand is None:
        return []
    
    # check if genomic region has polyA - if so, no good
    genome_buffer = 200
    if utr3['feature'].strand == '-':
        span = (int(utr3['cleavage_site']) - genome_buffer, int(utr3['cleavage_site']))
    else:
        span = (int(utr3['cleavage_site']), int(utr3['cleavage_site']) + genome_buffer)
    genome_seq = refseq.fetch(a['target'], span[0]-1, span[1])
    #genome_seq = self.refseq.GetSequence(align.target, span[0], span[1])
    if re.search('A{%s,}' % (homo_len), genome_seq, re.IGNORECASE) or re.search('T{%s,}' % (homo_len), genome_seq, re.IGNORECASE):
        sys.stdout.write('genome sequence has polyAT tract - no reliable link pairs can be retrieved %s %s %s:%s-%s\n' % 
                         (a['contig_seq'], utr3['cleavage_site'], a['target'], span[0], span[1]))
        return []

    mate_loc = {}
    #for read in self.bam.bam.fetch(align.query):
    for read in r2c.fetch(a['align'].qname, multiple_iterators=True):
        # skip when both mates mapped to same contig
        if not read.mate_is_unmapped and read.tid == read.rnext:
            continue
        
        # anchor read must be mapped entirely within contig
        if len(read.cigar) > 1:
            continue
        
        if anchor_read_strand == '+' and (read.is_reverse or read.pos + 1 > utr3['last_matched']):
            continue
        if anchor_read_strand == '-' and (not read.is_reverse or read.pos + read.rlen < utr3['last_matched']):
            continue
        
        if read.rnext >= 0:
            #mate_contig = self.bam.bam.getrname(read.rnext)             
            mate_contig = r2c.getrname(read.next_reference_id)
            if not mate_loc.has_key(mate_contig):
                mate_loc[mate_contig] = {}
            mate_loc[mate_contig][read.qname] = read
        else:
            pass
            #print 'cannot find unmapped mate %s' % read.qname
            
    link_pairs = []
    for contig in mate_loc.keys():
        #for read in self.bam.bam.fetch(contig):
        for mate in r2c.fetch(mate_contig, read.next_reference_start, read.next_reference_start+1):#,multiple_iterators=True):
            if not mate_loc[contig].has_key(mate.qname):
                continue
            
            trimmed_seq = mate.seq
            if args.trim_reads:
                trimmed_seq = trim_bases(mate.seq, mate.qual)

            if trimmed_seq:
                for base in ('A', 'T'):
                    if is_bridge_read_good(trimmed_seq, base, len(trimmed_seq) - max_mismatch, mismatch=[max_mismatch, len(trimmed_seq)]):
                        #print 'Read {} from contig {} has mate:'.format(read.query_name, r2c.getrname(read.reference_id))
                        #print 'Read {} from contig {}'.format(mate.query_name, mate_contig)
                        link_pairs.append([read, mate_loc[contig][read.qname], trimmed_seq])
                        break
    
    return link_pairs

def is_bridge_read_good(clipped_seq, base, min_len, mismatch):
    """Determines if clipped sequence is possible polyA tail
    
    If clipped_seq is composed of base only, then it is automatically 
    considered a potential bridge read regardless of length
    Otherwise will check the frequecy of 'the other bases' using is_polyA_tail()
    to determine whether it's acceptable
    """
    good = False
    if clipped_seq[0].upper() == base and len(re.sub(r'(.)\1+', r'\1', clipped_seq)) == 1:
        good = True
    elif is_polyA_tail(clipped_seq, base, min_len=min_len, max_nonAT_allowed=mismatch):
        good = True
    #if good:
        #print '{}-{}'.format(clipped_seq, good)
    return good

#def trim_bases(seq, qual, qual_dict, end=None):
#    tstart, tend = None, None

def trim_bases(seq, qual, end=None):
    """Trim poor quality bases from read sequence"""
    if poor_quals is None or poor_quals == '':
        return seq
    
    match_end = match_start = None
    match_end = re.search(r'[%s]+$' % poor_quals, qual)            
    match_start = re.search(r'^[%s]+' % poor_quals, qual)
    
    if match_start and not match_end:
        if end is None or end == 'start':
            return seq[match_start.end():]
        
    elif match_end and not match_start:
        if end is None or end == 'end':
            return seq[:match_end.start()]
        
    elif match_start and match_end:
        if end == 'start':
            return seq[match_start.end():]
        
        elif end == 'end':
            return seq[:match_end.start()]
        
        else:
            if len(match_start.group()) > len(match_end.group()):
                return seq[match_start.end():]
            elif len(match_end.group()) > len(match_start.group()):
                return seq[:match_end.start()]
    
    return seq

def show_trimmed_bases(seq, trimmed_bases):
    """Shows (link) read sequence with trimmed bases"""
    match_start = re.search('^' + trimmed_bases, seq)
    match_end = re.search(trimmed_bases + '$', seq)
    
    if match_start or match_end:
        if match_start:
            match = match_start
        else:
            match = match_end
        return seq[:match.start()].lower() + seq[match.start():match.end()].upper() + seq[match.end():].lower()
    
    else:
        return seq.lower()

def output_link_pairs(align, reads):
    """Outputs link pairs in FASTA format"""
    out = ''
    for r in reads:
        trimmed_seq = show_trimmed_bases(r[0].seq, r[2])
        num_trimmed_bases = r[0].rlen - len(r[2])
        if r[1].is_reverse:
            anchor_direction = 'L'
        else:
            anchor_direction = 'R'
        if r[0].is_unmapped:
            mate_contig = 'unmapped'
        else:
            mate_contig = r2c.getrname(r[0].tid)
        out += '>%s %s %s %s %s trimmed:%s\n%s\n' % (r[0].qname, align.qname, r[1].pos, anchor_direction, 
                                                     mate_contig, num_trimmed_bases, trimmed_seq) 
    return out

def annotate_cleavage_site(a, feature_list, cleavage_site, clipped_pos, base, fd, min_txt_match_percent=0.6):
    """Finds transcript where proposed cleavage site makes most sense, and also fetches matching ESTs

    This method assesses whether the cleavage site makes sense with
    a) the clipped position in the contig (start or end)
    b) the alignment strand
    c) the clipped base
    It determines what the strand the transcript should be on based on a) and b).
    Then it checks to see if the transcript strand agrees with c)
    (As an alternative, it can use splice sites, if the contig is large enough, to determine the 
    transcript strand, and then see whether both the clipped position and clipped base make sense)

    If the above makes sense, it moves to find transcript(s) that overlaps the alignment, given the
    'expected' strand and that cleavage site must exist 3' to the transcript's 3'UTR start, 
    using get_txts_with_min_match().
    If an annotated transcript whose end matches the cleavage site, that transcript is chosen.
    If not and more than one transcripts overlap, the one that is closest to the annotated 'end' is chosen.
    
    If a candidate transcript is chosen, it will try to find ESTs which end exactly at the 
    given cleavage site. (This is only done if the code is asked to at the command prompt)
    
    If no trancript can be found given the clipped postion, clipped base, and alignment, it will
    return None
    """     
    result = None
    
    chrom = proper_chrom(a['target'], chrom_proper=chrom_proper)
    
    # determine which transcript strand can the cleavage site come from
    if clipped_pos == 'start':
        if a['strand'] == '+':
            txt_strand = '-'
        else:
            txt_strand = '+'
    else:
        if a['strand'] == '+':
            txt_strand = '+'
        else:
            txt_strand = '-'
            
    #print 'txt_strand: {}'.format(txt_strand)
    #print 'base: {}'.format(base)
    if (txt_strand == '+' and base == 'A') or\
       (txt_strand == '-' and base == 'T'):     
        ests = []
        
        txts_screened = feature_list
        flength = len(txts_screened)
        for feature in txts_screened:
            if feature.strand != txt_strand:
                flength -= 1
#        print '{}\t{}'.format(cleavage_site,flength)

        if flength == 0:
            return result
        
        within_utr, identical = False,False
        closest_tid = None
        min_dist = 9000000
        closest = None

        if a['strand'] == '+':
            closest = sorted(a['close'],key=lambda(x):abs(cleavage_site - x[2].end))
            #for i in xrange(len(closest)):
                #closest[i][1] = abs(cleavage_site - closest[i][2].end)
        else:
            closest = sorted(a['close'],key=lambda(x):abs(cleavage_site - x[2].start))
            #for i in xrange(len(closest)):
                # Plus one because pysam grabs 1 before the actual start
                #closest[i][1] = abs(cleavage_site - closest[i][2].start + 1)
        #print 'closest: {}'.format(closest)
        #raw_input('*')
        if not closest:
            return result
        closest_within_utr3 = [x for x in closest if x[3]]
        if closest_within_utr3 and (closest_within_utr3[0][1] <= 20):
            closest = closest_within_utr3[0]
        else:
            closest = closest[0]
        if a['strand'] == '+':
            min_dist = abs(cleavage_site - closest[2].end)
        else:
            min_dist = abs(cleavage_site - (closest[2].start + 1))
        closest_tid = closest[0]
        if closest[1] == 0:
            identical = True
        if fd[a['target']][closest_tid]['utr3']:
            within_utr = True

        fd[a['target']][closest_tid]['cleavage_sites'].append(cleavage_site)

        if (min_dist <= thresh_dist) or (any([abs(cleavage_site - x) <= thresh_dist for x in fd[a['target']][closest_tid]['cleavage_sites']])):
            a['report_closest'] = False

        result = {
                'ests': ests,
                'txt': closest_tid,
                'novel': not identical, 
                'within_utr': within_utr,
                'coord': '%s:%d' % (a['target'], cleavage_site),
                'cleavage_site': cleavage_site,
                'from_end': min_dist,
                'which_end': clipped_pos,
                'base': base
            }

        #print '\tcs: {}\ttxt: {}'.format(cleavage_site,closest_tid)
    #print 'annotate_result: {}'.format(result)
    return result

def find_polyA_cleavage(a,gf,fd):
    gf = global_filters
    """Finds PolyA cleavage sites of a given aligned contig
    
    This method first checks if the given contig captures a polyA tail (find_tail_contig),
    and then tries to capture bridge reads (find_bridge_reads) given the above result.
    These methods may return multiple polyA tails, which will then be assessed one-by-one
    against the annotated gene models to determine if each tail is plausible or not 
    (annotate_cleavage_site).  The final result is kept in a dictionary where the key is
    either the 'start' or the 'end' of the contig.
    If plausible tails are observed from both the 'start' and 'end' of the contig, the
    contig is dismissed.        
    """

    tail = find_bridge_reads(a, gf['min_at'], gf['max_diff'], gf, tail=find_tail_contig(a, gf['min_at'], gf['max_diff']))
    results = []
    #print 'tail: {}'.format(tail)
    for clipped_pos in tail.keys():
        for event in tail[clipped_pos]:
            # contig coordinate of cleavage site                
            last_matched, cleavage_site, base, tail_seq, num_tail_reads, bridge_reads, bridge_clipped_seq = event
            result = annotate_cleavage_site(a, a['feature_list'], cleavage_site, clipped_pos, base, fd)
            if result:
                if bridge_reads:
                    result['num_bridge_reads'] = len(bridge_reads)
                    result['bridge_reads'] = bridge_reads
                    result['bridge_clipped_seq'] = bridge_clipped_seq
                result['base'] = base
                result['last_matched'] = last_matched
                result['clipped_pos'] = clipped_pos
                result['tail_seq'] = tail_seq
                result['num_tail_reads'] = num_tail_reads
                results.append(result)
    return results

def find_extended_bridge_reads(a, reads_to_screen, min_len, mismatch, genome_buffer=1000):
    global global_filters
    """Finds bridge reads where only ending portion represent polyA tail"""
    query_seqs = {}
    for reads in reads_to_screen.values():
        for read in reads:
            query_seqs[read.qname] = read.seq
    
    #entirely_mapped = self.align_transcript_seq(align, query_seqs, 'extended', self.get_full_blat_aln)
    
    clipped_reads = {'start':{}, 'end':{}}
    for clipped_pos, reads in reads_to_screen.iteritems():
        if reads:
            if (clipped_pos == 'start' and a['strand'] == '+') or\
               (clipped_pos == 'end' and a['strand'] == '-'):
                target_coord = [int(a['align'].reference_start)+1 - genome_buffer, int(a['align'].reference_end)]
            else:
                target_coord = [int(a['align'].reference_start)+1, int(a['align'].reference_end) + genome_buffer]
            
            query_seqs = dict((read.qname, read.seq) for read in reads)# if not entirely_mapped.has_key(read.qname))
            partial_aligns = align_genome_seq(a, query_seqs, target_coord, 'extended-bridge-genome', get_partial_blat_aln)
            
            read_objs = dict((read.qname, read) for read in reads)
            #print 'read_objs:\n{}'.format(read_objs)

            for read_name, mapped_coord in partial_aligns.iteritems():
                if mapped_coord[0] == 0:
                    clipped_seq = read_objs[read_name].seq[mapped_coord[1]:]
                else:
                    clipped_seq = read_objs[read_name].seq[:mapped_coord[0]]
                    
                if global_filters is not None and global_filters.has_key('min_bridge_size') and len(clipped_seq) < global_filters['min_bridge_size']:
                    continue
                                        
                # reverse complement to be in agreement with reference instead of contig
                clipped_seq_genome = clipped_seq
                if a['strand'] == '-':
                    clipped_seq_genome = revComp(clipped_seq)
                    
                if mapped_coord[0] == 0:
                    last_matched = read_objs[read_name].pos + mapped_coord[1]
                    pos_genome = target_coord[0] + mapped_coord[3] - 1
                else:
                    last_matched = read_objs[read_name].pos - mapped_coord[0]
                    pos_genome = target_coord[0] + mapped_coord[2]

                for base in ('A', 'T'):
                    if is_bridge_read_good(clipped_seq, base, min_len, mismatch):
                        #print 'possible extended bridge reads', align.query, read_objs[read_name].qname, read_objs[read_name].seq, is_seed, clipped_seq_genome
                        if not clipped_reads[clipped_pos].has_key(last_matched):
                            clipped_reads[clipped_pos][last_matched] = {}
                        if not clipped_reads[clipped_pos][last_matched].has_key(base):
                            clipped_reads[clipped_pos][last_matched][base] = []
                        clipped_reads[clipped_pos][last_matched][base].append([read_objs[read_name], clipped_seq_genome, pos_genome])
                            
    return clipped_reads

def merge_clipped_reads(clipped_reads, extended_clipped_reads):
    """Merges clipped_reads and extended_clipped_reads"""
    for clipped_pos in extended_clipped_reads.keys():
        for pos in extended_clipped_reads[clipped_pos].keys():
            if not clipped_reads[clipped_pos].has_key(pos):
                clipped_reads[clipped_pos][pos] = extended_clipped_reads[clipped_pos][pos]
            else:
                for base in extended_clipped_reads[clipped_pos][pos]:
                    if not clipped_reads[clipped_pos][pos].has_key(base):
                        clipped_reads[clipped_pos][pos][base] = extended_clipped_reads[clipped_pos][pos][base]
                    else:
                        clipped_reads[clipped_pos][pos][base].extend(extended_clipped_reads[clipped_pos][pos][base])

def find_bridge_reads(a, min_len, mismatch, gf, genome_buffer=1000, tail=None):
    #print "Running find_bridge_reads"
    """Finds bridge reads
    
    It first checks to see clipped reads if the entire clipped portion is A's or T's.
    Then it will check, through find_exteneded_bridge_reads() if the ending portion of 
    the clipped sequence is A's or T's.
    The 2 results will be merged together.
    """
    # used for check if read is mapped to the aligned portion of the contig
    query_bounds = sorted([int(a['qstart']), int(a['qend'])])
    
    # identify clipped reads that are potential pA/pT
    clipped_reads = {'start':{}, 'end':{}}
    second_round = {'start':[], 'end':[]}
    #for read in self.bam.bam.fetch(align.query):
    for read in r2c.fetch(a['align'].query_name):
        if not read.cigar or len(read.cigar) != 2:
            continue
        if (read.cigar[0][0] == 4 or read.cigar[0][0] == 5) or\
           (read.cigar[-1][0] == 4 or read.cigar[-1][0] == 5):
            # clipped at start
            if read.cigar[0][0] == 4 or read.cigar[0][0] == 5:
                clipped_pos = 'start'
                #print read
                #test = raw_input('Press any key to continue')
                last_matched = read.pos + 1
                # to see whether the clipped seq is a pA tail
                clipped_seq = read.seq[:read.cigar[0][1]]
                # for trimming of poor quality bases if so descired
                clipped_qual = read.qual[:read.cigar[0][1]]
                
            # clipped at end
            else:
                clipped_pos = 'end'
                last_matched = read.pos + read.alen
                # to see whether the clipped seq is a pA tail
                clipped_seq = read.seq[-1 * read.cigar[-1][1]:]
                # for trimming of poor quality bases if so descired
                clipped_qual = read.qual[-1 * read.cigar[-1][1]:]
                
            # Experimental feature
            offset = 33
            if args.trim_reads:
                if any([(ord(x)-offset) <= args.trim_reads for x in read.qual]):
                    continue
            # if last_match is beyond the limit of the alignment, adjust last_matched
            if last_matched < query_bounds[0] or last_matched > query_bounds[1]:
                if last_matched < query_bounds[0]:
                    diff = query_bounds[0] - last_matched
                else:
                    diff = last_matched - query_bounds[1]
                if clipped_pos == 'start':
                    last_matched = last_matched + diff
                    clipped_seq = read.seq[:read.cigar[0][1] + diff]
                    clipped_qual = read.qual[:read.cigar[0][1] + diff]
                else:
                    last_matched = last_matched - diff
                    clipped_seq = read.seq[-1 * (read.cigar[-1][1] + diff):]
                    clipped_qual = read.qual[-1 * (read.cigar[-1][1] + diff):]
                                    
            # trim poor quality base if desired
#            print read.qname
            #print 'clipped_seq pre-trim: {}'.format(clipped_seq)
            if args.trim_reads:
                clipped_seq = trim_bases(clipped_seq, clipped_qual, clipped_pos)
            #print 'clipped_seq post-trim: {}'.format(clipped_seq)
                
            #print '{}\t{}\t{}\t{}'.format(read.qname,read.cigar,read.seq,read.is_reverse)
            if len(clipped_seq) < 1:
                continue
            if (len(clipped_seq) < global_filters['min_bridge_size']):
                continue
            
            # reverse complement to be in agreement with reference instead of contig
            clipped_seq_genome = clipped_seq
            if a['strand'] == '-':
                clipped_seq_genome = revComp(clipped_seq)
            #print clipped_seq_genome
            
            # check for possible tail (stretch of A's or T's)
            pos_genome = qpos_to_tpos(a, last_matched)
            #print 'pos_genome:\n{}'.format(pos_genome)
            picked = False
            #print 'clipped_seq_genome:\n{}'.format(clipped_seq_genome)
            for base in ('A', 'T'):
                #print 'base: {}'.format(base)
                #print 'min_len: {}'.format(min_len)
                #print 'mismatch: {}'.format(mismatch)
                #raw_input('*'*20)
                #print '{}\t{}\t{}\t{}'.format(read.qname,base,clipped_pos,clipped_seq_genome)
                if is_bridge_read_good(clipped_seq_genome, base, min_len, mismatch):
                    #print 'yes'
#                    print 'below is good bridge_read'
#                    print '{}\t{}\t{}\t{}'.format(read.qname,base,clipped_pos,clipped_seq_genome)
                    if not clipped_reads[clipped_pos].has_key(last_matched):
                        clipped_reads[clipped_pos][last_matched] = {}
                    if not clipped_reads[clipped_pos][last_matched].has_key(base):
                        clipped_reads[clipped_pos][last_matched][base] = []
                    clipped_reads[clipped_pos][last_matched][base].append([read, clipped_seq_genome, pos_genome])
                    picked = True
                    reads_to_check.write('>{}\n{}\n'.format(read.qname,read.seq))#read_name,read_objs[read_name].seq))
                    
            if not picked:
                second_round[clipped_pos].append(read)
                #extended.write('>{}\t{}\n{}\n'.format(a['align'].qname, read.qname, read.seq))
    #print 'clipped_reads:\n{}'.format(clipped_reads)

    #extended_clipped_reads = find_extended_bridge_reads(a, second_round, min_len, mismatch)
    #merge_clipped_reads(clipped_reads, extended_clipped_reads)
    #print 'clipped_reads:\n{}'.format(clipped_reads)

    # filter events against reference sequence
    #filter_vs_reference(align, target, clipped_reads)
    
    # translate into results
    if tail is None:
        results = {}
        tail_search = None
    else:
        results = tail
        tail_search = {'start':{}, 'end':{}}
        for clipped_pos in tail.keys():
            for i in range(len(tail[clipped_pos])):
                cleavage_site = tail[clipped_pos][i][1]
                tail_search[clipped_pos][cleavage_site] = i
                            
    #filtered = copy.deepcopy(clipped_reads)
    for clipped_pos in clipped_reads:
        if clipped_pos not in results:
            results[clipped_pos] = []
        for pos in clipped_reads[clipped_pos]:
            skip = False
            for base in clipped_reads[clipped_pos][pos]:
                # check if clipped sequence is genomic/transcriptomic, otherwise skip
                is_genomic = is_transcriptome = bad_neighbor = False
    
                for read in clipped_reads[clipped_pos][pos][base]:
                    # Homopolymer neighbour check
                    if read[2] is not None and in_homopolymer_neighbor(a['target'], read[2], read[1], base):
                        skip = True
                        break
                        #del filtered[clipped_pos][pos][base]

                if (skip == True):
                    continue

                pos_genome = clipped_reads[clipped_pos][pos][base][0][2]
                #print 'pos_genome: {}'.format(pos_genome)
                
                if pos_genome is not None:
                    if tail_search is not None and tail_search[clipped_pos].has_key(pos_genome):
                        results_idx = tail_search[clipped_pos][pos_genome]
                        results[clipped_pos][results_idx][5] = [r[0] for r in clipped_reads[clipped_pos][pos][base]]
                        results[clipped_pos][results_idx][6] = [r[1] for r in clipped_reads[clipped_pos][pos][base]]
                    else:
                        results[clipped_pos].append([pos,
                                                     pos_genome, 
                                                     base,
                                                     None, 
                                                     None,
                                                     [r[0] for r in clipped_reads[clipped_pos][pos][base]], 
                                                     [r[1] for r in clipped_reads[clipped_pos][pos][base]]
                                                     ])      
    #print 'find_bridge_reads results:\n{}'.format(results)
    return results

def proper_chrom(chrom, genome=None, chrom_proper=None):
    """Returns proper chromosome name
    UCSC-format if available
    """
    if not chrom_proper and genome:
        chrom_proper = ucsc_chroms(genome)
    
    if chrom_proper:
        if chrom_proper.has_key(chrom):
            chrom = chrom_proper[chrom]
        elif chrom[:3].lower() == 'chr' and chrom_proper.has_key(chrom[3:]):
            chrom = chrom_proper[chrom[3:]]
    
    if not re.match('^(chr|scaffold)', chrom, re.IGNORECASE):
        chrom = 'chr' + chrom
            
    return chrom

def in_homopolymer_neighbor(chrom, pos, tail, base):
    """Checks if tail is juxtaposed against a homopolymer of the given base in genome
    
    The window size for checking is 5 or 2*length of the tail sequence (tail), whichever larger.
    The homopolyer run is constructed by making a string of the given base (base) with length 
    equal to the frequency of that base in the tail sequence.
    If the homopolyer is found immedicately before or after the given position (pos),
    then the result is True.
    """
    min_len = 5
    length = max(min_len, 2 * len(tail))
    try:
        neighbor_seq = refseq.fetch(chrom, pos-length, pos+length)
    except IndexError:
        return False
    #neighbor_seq = self.refseq.GetSequence(chrom, pos-length, pos+length)
    homo = ''
    freq = check_freq(tail)
    
    if freq.has_key(base):
        for i in range(length):
            homo += base.upper()
                
        m = re.search(homo, neighbor_seq, re.IGNORECASE)
        if m:               
            # if homopolymer is touching the middle base of neighbor_seq
            if m.start() <= len(neighbor_seq)/2 + 1 and m.end() - 1 >= len(neighbor_seq)/2 - 1:
                return True
                    
    return False

def get_num_tail_reads(align, last_matched):
    """Reports number of reads spanning cleavage site in contig"""
    num = 0
    for read in r2c.fetch(align.qname):
        if not read.cigar or len(read.cigar) != 1:
            continue
        
        if read.pos + 1 <= last_matched and read.pos + read.alen > last_matched:
            num += 1
    
    return num

def find_tail_contig(a, min_len, mismatch):
    """Finds contigs that have polyA tail reconstructed"""
    # size of stretch of contig sequence to be added to polyA tail
    # against trasncript sequences to see it polyA tail is genomic
    junction_buffer=50
    results = {}
    
    clipped = {'start':False, 'end':False}
    if int(a['qstart']) > 1:
        clipped['start'] = True
    
    if int(a['qend']) < int(a['align'].infer_query_length(True)):
        clipped['end'] = True

    for clipped_pos in ('start', 'end'):
        if clipped[clipped_pos]:
            if clipped_pos == 'start':
                last_matched = int(a['qstart'])
                clipped_seq = a['contig_seq'][:int(a['qstart'])-1] # used to have -1
                junction_seq = a['contig_seq'][:int(a['qstart']) -1 + junction_buffer] # used to have -1
            else:
                last_matched = int(a['qend'])
                clipped_seq = a['contig_seq'][int(a['qend']):]
                junction_seq = a['contig_seq'][int(a['qend']) - junction_buffer:]
                                    
            cleavage_site = qpos_to_tpos(a, last_matched)
            clipped_seq_genome = clipped_seq
            if a['strand'] == '-':
                clipped_seq_genome = revComp(clipped_seq)
            
            matched_transcript = in_homopolymer = False
            for base in ('A', 'T'):
                if matched_transcript or in_homopolymer:
                    continue
                                        
                perfect = is_polyA_tail(clipped_seq_genome, base, min_len=1, max_nonAT_allowed=[0, 1])
                imperfect = is_polyA_tail(clipped_seq_genome, base, min_len=min_len, max_nonAT_allowed=mismatch)
                #print 'perfect: {}'.format(perfect)
                #print 'imperfect: {}'.format(imperfect)
                
                if perfect or imperfect:
                    # don't need to do the following 2 checks if it's not a potential tail
                    if len(clipped_seq) == 1 and in_homopolymer_neighbor(a['target'], cleavage_site, clipped_seq_genome, clipped_seq[0]):
#                        print '%s : clipped seq in middle of homopolyer run (%s) %s' % (align.qname, clipped_pos, clipped_seq)
                        in_homopolymer = True
                        continue
                    
                    #entirely_mapped = align_transcript_seq(align, target, {align.qname:junction_seq}, 'junction', get_full_blat_aln)
                    #if entirely_mapped.has_key(align.qname):
                    #    print '%s : clipped seq is just transcript seq (%s) %s %d' % (align.qname, clipped_pos, clipped_seq, len(clipped_seq))
                    #    matched_transcript = True
                    #    continue
                    
                    # find reads corresponding to tail
                    num_tail_reads = get_num_tail_reads(a['align'], last_matched)
                                                                    
                    if not results.has_key(clipped_pos):
                        results[clipped_pos] = []
                    results[clipped_pos].append([last_matched,
                                                 cleavage_site,
                                                 base,
                                                 clipped_seq_genome,
                                                 num_tail_reads,
                                                 None,
                                                 None,
                                                 ]
                                                )
            
                                        
    if not clipped['start'] and not clipped['end']:
        results = None
        
#    print 'find_tail_contig results:\n{}'.format(results)
    return results

def align_transcript_seq(align, target, query_seqs, label, parse_fn):
    print 'aligning transcript seq'
    start = time.time()
    """Aligns(BLAT) query sequences to transcripts overlapping alignment
    
    Query sequences is given in a hash (query_seq) where
    key=query name, value=query sequence
    Targets are all transcript sequences overlapping 'tstart' and 'tend'
    of the alignment.
    'label' will be added in addition to query name to all the 
    temporary files
    The BLAT alignments will be processed by parse_fn() to return the results
    """
    result = None
    
    #for cleanup
    tmp_files = []
    
    if args.use_tmp:
        path = '/tmp'
    else:
        path = os.path.dirname(os.path.abspath(args.out))
    
    # create transcript sequence file for Blatting
    target_file = '%s/%s-target-%s.fa' % (path, align.qname, label)
    if os.path.exists(target_file):
        os.remove(target_file)
    target_empty = extract_transcript_seq(align, target, target_file)
    tmp_files.append(target_file)

    # create query for Blatting
    query_file = '%s/%s-query-%s.fa' % (path, align.qname, label)
    if os.path.exists(query_file):
        os.remove(query_file)
    out = open(query_file, 'w')
    for query, seq in query_seqs.iteritems():
        out.write('>%s\n%s\n' % (query, seq))
    out.close()
    tmp_files.append(query_file)

    # align query against target
    aln_file = '%s/%s-%s.psl' % (path, align.qname, label)
    if os.path.exists(aln_file):
        os.remove(aln_file)
    tmp_files.append(aln_file)
    try:
        FNULL = open(os.devnull, 'w')
        task = subprocess.Popen(['blat', target_file, query_file, aln_file], stdout=FNULL)
        task.communicate()
    except CalledProcessError as err:
        sys.stderr.write('error running blat:%s' % err.cmd)
    else:
        if os.path.exists(aln_file):
            result = parse_fn(aln_file)
            
    # clean up temporary alignment files
    for ff in tmp_files:
        if os.path.exists(ff):
            #print 'cleanup', ff
            os.remove(ff)
            
    return result

def extract_transcript_seq(align, target, out_file):
    """Extracts transcripts overlapping an alignment and outputs their
    sequences to the given file
    
    Returns empty=True if it fails to find overlapping transcripts
    """
    
    chrom = proper_chrom(target, chrom_proper=chrom_proper)
    feats = features.fetch(chrom, align.reference_start, align.reference_end)
    out = open(out_file, 'w')
    empty = True
    transcripts = {}
    for feature in feats:
        if (feature.transcript_id not in transcripts):
            transcripts[feature.transcript_id] = ''
        if feature.feature == 'exon':
            empty = False
            transcripts[feature.transcript_id] += refseq.fetch(chrom, feature.start, feature.end).upper()
    for tid in transcripts:
        #print "transcript_seq: {}".format(transcripts[tid])
        out.write('>%s\n%s\n' % (tid, transcripts[tid]))
    out.close()
    
    return empty

def get_full_blat_aln(aln_file):
    """Extracts full hits from BLAT aligments
    
    This is used for removing false-positive bridge reads where their entirety in
    sequence can be mapped to a single transcript
    """
    fully_aligned = {}
    for line in open(aln_file, 'r'):
        if not re.search('^\d', line):
            continue
        cols = line.rstrip('\n').split('\t')
        query, qsize, qstart, qend, target = cols[9:14]
        block_count = cols[17]
        #print 'ver2 aln', query, qsize, qstart, qend, block_count, target
        
        if int(qstart) == 0 and int(qsize) == int(qend) and int(block_count) == 1:
            fully_aligned[query] = True
        
    return fully_aligned

def get_partial_blat_aln(aln_file):
    """Extracts single-block, partial hits from BLAT aligments
    
    This is for capturing the polyA tails of extended bridge reads/
    The clipped portion of extended bridge reads contains both genomic sequence
    and the polyA tail. By alignment these reads to the genmomic region, the 
    polyA tail should be unaligned whereas the genomic portion would align.
    The return variable is a dictionary, where:
    key = query(read) name
    value = [qstart, qend, tstart, tend]
    """

    partially_aligned = {}
    for line in open(aln_file, 'r'):
        if not re.search('^\d', line):
            continue
        cols = line.rstrip('\n').split('\t')
        query, qsize, qstart, qend, target, tsize, tstart, tend = cols[9:17]
        block_count = cols[17]
        
        if int(block_count) == 1 and \
            ((int(qstart) == 0 and int(qend) < int(qsize)) or\
             (int(qstart) > 0 and int(qsize) == int(qend))):
            partially_aligned[query] = [int(qstart), int(qend), int(tstart), int(tend)]
        
    return partially_aligned

def align_genome_seq(a, query_seqs, coord, label, parse_fn):
    start = time.time()
    """Aligns(BLAT) query sequences to genomic sequence of given coordinate
    
    Query sequences is given in a hash (query_seq) where
    key=query name, value=query sequence
    Target is genomic sequence between coord[0] and coord[1]
    'label' will be added in addition to query name to all the 
    temporary files
    The BLAT alignments will be processed by parse_fn() to return the results
    """
    result = None
    
    #target_seq = self.refseq.GetSequence(align.target, coord[0], coord[1])
    target_seq = refseq.fetch(a['target'], max(0, coord[0]-1), coord[1]).upper()
    
    #for cleanup
    tmp_files = []
    
    if args.use_tmp:
        path = '/tmp'
    else:
        path = os.path.dirname(os.path.abspath(args.out))
    
    # create genome reference file for Blatting
    target_file = '%s/%s-genome-%s.fa' % (path, a['align'].qname, label)
    if os.path.exists(target_file):
        os.remove(target_file)
    out = open(target_file, 'w')
    out.write('>%s:%d-%d\n%s\n' % (a['target'], coord[0], coord[1], target_seq))
    out.close()
    tmp_files.append(target_file)

    # create query for Blatting
    query_file = '%s/%s-query-%s.fa' % (path, a['align'].qname, label)
    if os.path.exists(query_file):
        os.remove(query_file)
    out = open(query_file, 'w')
    for query, seq in query_seqs.iteritems():
        out.write('>%s\n%s\n' % (query, seq))
    out.close()
    tmp_files.append(query_file)

    # align query against target
    aln_file = '%s/%s-%s.psl' % (path, a['align'].qname, label)
    if os.path.exists(aln_file):
        os.remove(aln_file)
    tmp_files.append(aln_file)
    try:
        FNULL = open(os.devnull, 'w')
        task = subprocess.Popen(['blat', target_file, query_file, aln_file], stdout=FNULL)
        task.communicate()
    except CalledProcessError as err:
        sys.stderr.write('error running blat:%s' % err.cmd)
    else:
        if os.path.exists(aln_file):
            result = parse_fn(aln_file)
            
    # clean up temporary alignment files
    for ff in tmp_files:
        if os.path.exists(ff):
            #print 'cleanup', ff
            os.remove(ff)
            
    return result

def output_result(a, result, output_fields, fd, link_pairs=[]):
    """Outputs main results in tab-delimited format
    
    The output fields are specified in class variable 'output_fields'
    """
    data = {}
    data['contig'] = a['align'].query_name
    
#    if ('feature' in result['txt']):
#        data['transcript'] = result['txt']['feature'].transcript_id#['transcript_id']
#        data['transcript_strand'] = result['txt']['feature'].strand
#        if result['txt']['feature'].gene_id:
#            data['gene'] = result['txt']['feature'].gene_id
#    else:
#        data['transcript'] = result['txt']['transcript_id']
#        data['transcript_strand'] = result['txt']['strand']
    data['transcript'] = result['txt']
    data['transcript_strand'] = fd[a['target']][result['txt']]['feats'][0].strand
    data['gene'] = fd[a['target']][result['txt']]['feats'][0].asDict()['gene_id']
    data['coding'] = get_coding_type(fd[a['target']][result['txt']])
#    if coding_type == 'CODING':
#        data['coding'] = 'yes'
#    elif coding_type == 'NONCODING':
#        data['coding'] = 'no'
#    else:
#        data['coding'] = 'unknown'
        
    if result['within_utr']:
        data['within_UTR'] = 'yes'
    else:
        data['within_UTR'] = 'no'
        
    data['chromosome'] = a['target']
    data['cleavage_site'] = result['cleavage_site']
    data['distance_from_annotated_site'] = result['from_end']
    
    if 'tail_seq' in result:
        if result['tail_seq'] is None:
            data['length_of_tail_in_contig'] = 0
            data['number_of_tail_reads'] = 0
        else:
            data['length_of_tail_in_contig'] = len(result['tail_seq'])
            data['number_of_tail_reads'] = result['num_tail_reads']
    else:
        data['length_of_tail_in_contig'] = data['number_of_tail_reads'] = '-'
        
    if ('bridge_reads' in result) and (result['bridge_reads']):
        data['number_of_bridge_reads'] = len(result['bridge_reads'])
        data['max_bridge_read_tail_length'] = max([len(s) for s in result['bridge_clipped_seq']])
        data['bridge_read_identities'] = ','.join([read.qname for read in result['bridge_reads']])
    else:
        data['number_of_bridge_reads'] = data['max_bridge_read_tail_length'] = 0
        data['bridge_read_identities'] = '-'
        
    if link_pairs is not None:
        if link_pairs:
            data['number_of_link_pairs'] = len(link_pairs)
            data['link_pair_identities'] = ','.join(r[0].qname for r in link_pairs)
            data['max_link_pair_length'] = max([len(r[-1]) for r in link_pairs])
        else:
            data['number_of_link_pairs'] = data['max_link_pair_length'] = 0
            
    elif not args.link:
        data['number_of_link_pairs'] = data['max_link_pair_length'] = 0
            
    if result['ests'] is not None:
        data['ESTs'] = '%s' % (len(result['ests']))
    else:
        data['ESTs'] = '-'
            
    data['tail+bridge_reads'] = 0
    if data['number_of_tail_reads'] != '-':
        data['tail+bridge_reads'] += data['number_of_tail_reads']
    if data['number_of_bridge_reads'] != '-':
        data['tail+bridge_reads'] += data['number_of_bridge_reads']

    if a['binding_sites']:
        data['hexamer_loc+id'] = (';').join([str(x[0])+':'+str(x[1]) for x in a['binding_sites']])
    
    # A random utr3 is chosen out of the list
    if (result['txt']) and (result['txt'] in a['utr3s']):
        data['3UTR_start_end'] = str(a['utr3s'][result['txt']][0])+'-'+str(a['utr3s'][result['txt']][1])
#        N = ''
#        if ('novel' in utr3) and (utr3['novel'] == True):
#            N = 'N'
#        if ('feature' in utr3):
#            data['3UTR_start_end'] = N+str(utr3['feature'].start)+'-'+str(utr3['feature'].end)
#        else:
#            data['3UTR_start_end'] = N+str(utr3['start'])+'-'+str(utr3['end'])
    
    cols = []
    for field in output_fields:
        if field in data:
            cols.append(str(data[field]))
        else:
            cols.append('-')
            
    result = '%s\n' % '\t'.join(cols)
    return result

def output(report_lines, bridge_lines=None, link_lines=None):
    """Writes output lines to files"""      
    if out_result is not None and report_lines:
        for line in report_lines:
            out_result.write(line)
    out_result.close()
            
#    if bridge_lines and out_bridge_reads is not None:
#        for line in bridge_lines:
#            out_bridge_reads.write(line)
#    if out_bridge_reads is not None:
#        out_bridge_reads.close()
#        
#    if link_lines and out_link_pairs is not None:
#        for line in link_lines:
#            out_link_pairs.write(line)
#    if out_link_pairs is not None:
#        out_link_pairs.close()

def fetchUtr3(alignd, feature_list):
    utr3 = {'start': None, 'end': None}
    if (alignd['strand'] == '+'):
        utr3['strand'] = '+'
        CDS = False
        for f in reversed(feature_list):
            if (utr3['end'] == None) and (f['feature'].feature == 'exon'):
                utr3['end'] = f['feature'].end
            if (f['feature'].feature == 'stop_codon'):
                utr3['start'] = f['feature'].end + 1
            elif (f['feature'].feature == 'CDS') and (CDS == False):
                utr3['start'] = f['feature'].end + 4
                CDS = True
            if (utr3['start'] and utr3['end']) and (utr3['start'] <= utr3['end']):
                return utr3
    if (alignd['strand'] == '-'):
        utr3['strand'] = '-'
        CDS = False
        for f in feature_list:
            if (utr3['start'] == None) and (f['feature'].feature == 'exon'):
                utr3['start'] = f['feature'].start
            if (f['feature'].feature == 'stop_codon'):
                utr3['end'] = f['feature'].start - 1
            elif (f['feature'].feature == 'CDS') and (CDS == False):
                utr3['end'] = f['feature'].start - 4
                CDS = True
            if (utr3['start'] and utr3['end']) and (utr3['start'] <= utr3['end']):
                return utr3
    if (utr3['start'] > utr3['end']):
        utr3['start'] = utr3['end'] = None
    return utr3

def findBindingSites(a, cleavage_site):
    binding_sites = {'AATAAA':1,'ATTAAA':2,'AGTAAA':3,'TATAAA':4,
                     'CATAAA':5,'GATAAA':6,'AATATA':7,'AATACA':8,
                     'AATAGA':9,'AAAAAG':10,'ACTAAA':11,'AAGAAA':12,
                     'AATGAA':13,'TTTAAA':14,'AAAACA':15,'GGGGCT':16}
    results = []
    if not args.strand_specific:
        seq = refseq.fetch(a['target'],cleavage_site-50,cleavage_site+50).upper()
        for i in xrange(len(seq)):
            hexamer = seq[i:i+6]
            rev = revComp(seq[i:i+6])
            if (hexamer in binding_sites):
                results.append([i+cleavage_site-50,binding_sites[hexamer]])
            if (rev in binding_sites):
                results.append([cleavage_site+i+6,binding_sites[rev]])
        return results
    else:
        if a['strand'] == '+':
            seq = refseq.fetch(a['target'],cleavage_site-50,cleavage_site)
            for i in xrange(len(seq)):
                hexamer = seq[i:i+6]
                if (hexamer in binding_sites):
                    results.append([i+cleavage_site-50,binding_sites[hexamer]])
        else:
            seq = refseq.fetch(a['target'],cleavage_site,cleavage_site+50)
            for i in xrange(len(seq)):
                hexamer = revComp(seq[i:i+6])
                if (hexamer in binding_sites):
                    results.append([cleavage_site+i+6,binding_sites[hexamer]])
    # Enable this to fix BTL-513
    if results:
        results = [sorted(results, key=lambda(x):x[1])[0]]
    return results

def findNovel3UTR(a):
    #print 'Finding novel 3\'UTR'
    stops = ['TGA','TAA','TAG']
    if a['strand'] == '+':
        utr3 = {'start': None, 'end': a['align'].reference_end, 'novel': True}
    else:
        utr3 = {'start': None, 'end': a['align'].reference_start, 'novel': True}
    # 60 is sort of the minimum length for a human 3'UTR
    seqlen = len(a['contig_seq'])
    for i in xrange(seqlen-60,max(seqlen-2000,0),-1):
        codon = a['contig_seq'][i:i+3]
        if (codon in stops):
            for j in xrange(max(0,i-(seqlen/4)),0,-3):
                start = a['contig_seq'][j:j+3]
                if (start == 'ATG'):
                    utr3['start'] = qpos_to_tpos(a,i+5)
                    if (utr3['start'] > utr3['end']):
                        temp = utr3['start']
                        utr3['start'] = utr3['end']
                        utr3['end'] = temp
                    if (utr3['start']) and (utr3['end']):
                        return utr3
    return None

def loadResults(outfile):
    with open(outfile, 'r') as o:
        line = o.readline()
        return ('').join(o.readlines())

def fetchUtrs(a, feature_list, first_cds, last_cds):
    utr3s,utr3e,utr5s,utr5e = [None]*4
    fl = feature_list
    maxval = 0
    for i in xrange(len(fl)-1,-1,-1):
        if (fl[i]['feature'].end > maxval):
            maxval = fl[i]['feature'].end
    if (a['strand'] == '+'):
        utr3s = fl[last_cds]['feature'].end
        utr3e = maxval
        utr5s = fl[0]['feature'].start
        utr5e = fl[first_cds]['feature'].start
    else:
        utr5s = fl[last_cds]['feature'].end
        utr5e = maxval
        utr3s = fl[0]['feature'].start
        utr3e = fl[first_cds]['feature'].start
    return [[utr3s,utr3e],[utr5s,utr5e]]

def filter_contig_sites(contig_sites,fd):
    res = []
    for event in contig_sites:
        if (all([(abs(event['cleavage_site'] - cs ) >= 20) for cs in fd[event['a']['target']][event['txt']]['cleavage_sites']])):
            res.append(event)
    return res

# Short variable for variable containing all result info
#ar = global_filters['all_results']
# If the distance between any transcript end and the contig end is
# less than this value, the transcript end should be reported as a cleavage event
contig_sites = []
thresh_dist = 20
lines_result = lines_bridge = lines_link = ''
for align in aligns:
    # If contigs are specified only look at those
    if args.c:
        if align.query_name not in args.c:
            continue
    #sys.stdout.write('{}-{}-{}{}\n'.format(align.qname,align.reference_start, align.reference_end,'*'*10))
    print '{}\t{}\t{}'.format(align.qname,align.reference_start, align.reference_end)
    # If the contig has no start or no end coordinate, we can't
    # do any analysis on it, so we must skip it
    if (align.reference_start == None) or (align.reference_end == None):
        continue
    # tids              = Set of transcript ids that overlap contig
    # feature_list      = List of features overlapping contig
    # closest_tid       = Set the closest transcript to the end of the contig
    # closest_feat      = Set the closest feature to the end of the contig
    # report_closest    = Whether to report the closest transcript end as a cs
    # min_dist          = The minimum distance between any transcript and the contig
    a = {'align': align,'closest_tid': None,'closest_feat': None,
         'report_closest': False, 'feature_list': [], 'tids': set(),
         'min_dist': 1000000, 'utr3s': {}, 'utr5s': {}, 'close':[],'base': None}
    # Get target/chromosome
    a['target'] = aligns.getrname(align.tid)
    # Get the sequence of the contig
    a['contig_seq'] = contigs.fetch(align.query_name)
    # Get the overlapping features
    try:
        feats = features.fetch(a['target'], align.reference_start, align.reference_end)
    # If fails, skip this contig
    except ValueError:
        continue
    # If the library is strand specific, assume the contig strand is correct
    if args.strand_specific:
        if (align.is_reverse):
            a['strand'] = '-'
        elif not (align.is_reverse):
            a['strand'] = '+'
    # Store all transcript id's
    for f in feats:
        tid = f.asDict()['transcript_id']
        if (args.strand_specific):
            if (f.strand != a['strand']):
                continue
        a['feature_list'].append(f)
        a['tids'].add(tid)
    # If not a strand specific library, infer the strand by looking at the overlapping features
    # First, count how many overlapping transcripts are + and -
    likely_strand = {'-': 0, '+': 0}
    for tid in a['tids']:
        likely_strand[feature_dict[a['target']][tid]['strand']] += 1
    if not a['strand']:
        a['strand'] = max(likely_strand, key=lambda x: likely_strand[x])
    # If the feature_list is empty, skip this contig
    if not a['feature_list']:
        continue
    # Go through the list of transcripts and find the one closest
    # to the end of the contig
    for t in a['tids']:
        has_utr3 = False
        # If the transcript is of a different strand than the contig, skip it
        if (feature_dict[a['target']][t]['feats'][0].strand != a['strand']):
            continue
        if (a['strand'] == '+'):
            subclosest_feat = feature_dict[a['target']][t]['feats'][feature_dict[a['target']][t]['maxfeat']]
            dist = abs(subclosest_feat.end - align.reference_end)
        else:
            subclosest_feat = feature_dict[a['target']][t]['feats'][0]
            dist = abs(subclosest_feat.start - align.reference_start)
        if feature_dict[a['target']][t]['utr3']:
            has_utr3 = True
        a['close'].append([t,dist,subclosest_feat,has_utr3])
    # Get 3utrs for all overlapping transcripts
    for t in a['tids']:
        utr3 = feature_dict[a['target']][t]['utr3']
        if (utr3):
            a['utr3s'][t] = utr3
    #a['close'] = [x for x in a['close'] if x[1] < thresh_dist]
    if (a['close']):
        a['close'] = sorted(a['close'], key=lambda(x):x[1])
        #print 'a[close]: {}'.format([[x[0],x[2].start,x[2].end] for x in a['close']])
        within_utr3 = [x for x in a['close'] if x[3] and (x[1] <= thresh_dist)]
        if within_utr3:
            within_utr3 = within_utr3[0]
            a['closest_tid'],a['min_dist'],a['closest_feat'] = within_utr3[:-1]
        else:
            a['closest_tid'],a['min_dist'],a['closest_feat'] = a['close'][0][:-1]
    if (a['min_dist'] <= thresh_dist):
        a['report_closest'] = True
    # Skip contig if there is no feature close to it
    if not a['closest_tid']:
        continue
    # Get query blocks
    a['qblocks'] = cigarToBlocks(align.cigar, align.reference_start, a['strand'])[1]
    if not a['qblocks']:
        continue
    a['qstart'] = min(a['qblocks'][0][0], a['qblocks'][0][1], a['qblocks'][-1][0], a['qblocks'][-1][1])
    a['qend'] = max(a['qblocks'][0][0], a['qblocks'][0][1], a['qblocks'][-1][0], a['qblocks'][-1][1])
    result_link = link_pairs = None
    #for k in a:
    #    logger.debug(k)
    #    logger.debug(a[k])
    results = find_polyA_cleavage(a,global_filters,feature_dict)
    #print 'results: {}'.format(results)
    if (a['report_closest']):
        if (a['strand'] == '+'):
            cs = align.reference_end
        else:
            cs = align.reference_start+1
        res = {'txt': a['closest_tid'], 'cleavage_site': cs, 'within_utr': True,
               'from_end': a['min_dist'], 'ests': None, 'a': {'target': a['target'],
               'align': align, 'utr3s': a['utr3s']}}
        try:
            res['a']['binding_sites'] = findBindingSites(a, res['cleavage_site'])
        except TypeError:
            res['a']['binding_sites'] = None
        contig_sites.append(res)
    if results:
        for result in results:
            # If there is already a cs close to the end, we don't need the implied one
            try:
                a['binding_sites'] = findBindingSites(a, result['cleavage_site'])
            except TypeError:
                a['binding_sites'] = None
            lines_result += output_result(a, result, output_fields, feature_dict, link_pairs=link_pairs)
            # check if chrom is in all_results
                
contig_sites = filter_contig_sites(contig_sites,feature_dict)
for result in contig_sites:
    #print result
    #print output_result(result['a'], result, output_fields, feature_dict, link_pairs=link_pairs)
    #raw_input('^'*20)
    lines_result += output_result(result['a'], result, output_fields, feature_dict, link_pairs=link_pairs)

print lines_result
# close output streams
reads_to_check.close()
FNULL = open(os.devnull, 'w')
bstart = [time.time(),time.strftime("%c")]
print "Blatting extended bridge reads..."
#if not os.path.isfile('./.blat_alignment'):
reads_to_check = os.path.join(os.path.dirname(args.out),'.reads_to_check')
blat_alignment = os.path.join(os.path.dirname(args.out),'.blat_alignment')
task = subprocess.Popen(['blat', args.ref_genome, reads_to_check, blat_alignment], stdout=FNULL)
task.communicate()
print "Blat complete"
#bend = [time.time(), time.strftime("%c")]
blat_results = get_full_blat_aln(blat_alignment)
#print 'blat_results: {}'.format(blat_results)
lines_result = lines_result.splitlines()
keep = []
print 'lines_result: {}'.format(lines_result)
for result in lines_result:
    result = result.split('\t')
    bridge_reads = result[14].split(',')
    for read in bridge_reads:
        if read in blat_results:
            result[12] = str(int(result[12]) - 1)
            bridge_reads.remove(read)
            if (int(result[12]) == 0) and (result[10] == '-'):
                continue
    if bridge_reads:
        result[14] = bridge_reads
    else:
        result[14] = '-'
    keep.append(('\t').join(result))
lines_result = ('\n').join(keep)
group_and_filter(lines_result, args.out+'.KLEAT', filters=global_filters, make_track=args.track, rgb=args.rgb)
