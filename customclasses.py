class Transcript:
    
    def __init__(self, name, chrom, tid, strand, cstart, cend, tstart, tend):
        self.name = name
        self.chrom = chrom
        self.tid = tid
        self.strand = strand
        self.cstart = cstart
        self.cend = cend
        self.tstart = tstart
        self.tend = tend
        self.seq = None

    def get_utr3(self):
        if self.cend and (self.tend > self.cend+3):
            return [self.cend+3,self.tend]
        else:
            return None

class Contig:
    
    def __init__(self,name,target=None,qstart=None,qend=None,tstart=None,tend=None,closest_tid=None,tblocks=None,tids=None,strand=None,cigar=None,seq=None,utr3s=None,polya_signals=None):
        self.name = name
        self.target = target
        self.tstart = tstart
        self.tend = tend
        self.qstart = qstart
        self.qend = qend
        self.closest_tid = closest_tid
        self.tblocks = tblocks
        self.tids = tids
        self.strand = strand
        self.cigar = cigar
        self.seq = seq
        self.utr3s = utr3s
        self.polya_signals = polya_signals

    def get_qblocks(self):
        qblocks = []
        start = self.qstart
        for block in self.tblocks:
            diff = block[1]-block[0]
            qblocks.append([start, start + diff])
            start = start + diff
        return qblocks

class Read:
    
    def __init__(self,name,clipped_pos,seq,qblocks,cigar,strand,target,qstart,qend,tstart,tend,clipped_seq,is_bridge,potential_bridge):
        self.name = name
        self.clipped_pos = clipped_pos
        self.seq = seq
        self.qblocks = qblocks
        self.cigar = cigar
        self.strand = strand
        self.target = target
        self.qstart = qstart
        self.qend = qend
        self.tstart = tstart
        self.tend = tend
        self.clipped_seq = clipped_seq
        self.is_bridge = is_bridge
        self.potential_bridge = potential_bridge

class Cleavage_event:

    def __init__(self,gene,transcript,transcript_strand,coding,contig,chromosome,coordinate,within_utr3,distance_from_annot,len_contig_tail,tail_ids,bridge_ids,num_link_pairs,max_link_len,link_ids,polya_signals,utr3_coords):
        self.gene = gene
        self.transcript = transcript
        self.transcript_strand = transcript_strand
        self.coding = coding
        self.contig = contig
        self.chromosome = chromosome
        self.coordinate = coordinate
        self.within_utr3 = within_utr3
        self.distance_from_annot = distance_from_annot
        self.len_contig_tail = len_contig_tail
        self.tail_ids = tail_ids
        self.bridge_ids = bridge_ids
        self.num_link_pairs = num_link_pairs
        self.max_link_len = max_link_len
        self.link_ids = link_ids
        self.polya_signals = polya_signals
        self.utr3_coords = utr3_coords
