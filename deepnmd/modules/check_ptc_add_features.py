#!/usr/bin/env python3
"""
PTC (Premature Termination Codon) Analyzer for VEP output - WITH CAI, m6A, DOWNSTREAM AUG, GENE EXPRESSION, AND TRANSLATIONAI FASTA OUTPUT
Analyzes variants in protein-coding regions to identify PTCs and calculate various features
Supports multiple VEP output formats: VCF with CSQ, tab-delimited with Extra column, and standard tab-delimited

"""

import argparse
import sys
import re
import gzip
import math
import pandas as pd
import pysam
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def find_downstream_augs(modified_sequence, ptc_cds_pos):
    """
    Find the first in-frame and out-of-frame AUG codons downstream of a PTC.
    
    Args:
        modified_sequence: CDS sequence with the PTC already present
        ptc_cds_pos: 1-based position of the PTC in the modified CDS
    
    Returns:
        (inframe_distance, outframe_distance)
        where distances are from PTC start position to AUG start position (in nucleotides)
        Returns None for each distance if no corresponding AUG is found
        Note: Calling code will fill None values with 100000
    """
    # Convert to 0-based position
    ptc_position = ptc_cds_pos - 1
    
    # Find which codon the PTC is in
    ptc_codon_start = (ptc_position // 3) * 3
    
    # Make sure it is a complete codon
    if ptc_codon_start + 2 >= len(modified_sequence):
        return None, None
    
    # Search for AUGs starting after the PTC codon
    search_start = ptc_codon_start + 3  

    if search_start >= len(modified_sequence):
        return None, None
    
    # Look for all ATG occurrences downstream
    aug_pattern = re.compile(r'ATG', re.IGNORECASE)
    
    inframe_distance = None
    outframe_distance = None
    
    # Search through the sequence for AUGs
    for match in aug_pattern.finditer(modified_sequence, search_start):
        aug_position = match.start()
        
        # Check if this AUG is in-frame relative to the original CDS
        # In-frame: aug_position % 3 == 0 
        is_inframe = (aug_position % 3) == 0
        
        if is_inframe and inframe_distance is None:
            # First in-frame AUG found
            inframe_distance = aug_position - ptc_codon_start
        elif not is_inframe and outframe_distance is None:
            # First out-of-frame AUG found
            outframe_distance = aug_position - ptc_codon_start
        
        # Stop searching if found both
        if inframe_distance is not None and outframe_distance is not None:
            break
    
    return inframe_distance, outframe_distance


class StandaloneCAICalculator:
    """Standalone CAI calculator implementing Sharp & Li (1987) algorithm."""
    
    def __init__(self):
        """Initialize with human codon usage reference table."""
        self.stop_codons = {'TAA', 'TAG', 'TGA'}
        self.genetic_code = self._get_genetic_code()
        self.codon_weights = self._build_codon_weights()
    
    def _get_kazusa_human_frequencies(self):
        """
        Get exact human codon usage frequencies from HIVE-CUT database.
        Source: 123,938 human CDS sequences with 80,831,647 codons.
        Frequencies are per thousand codons (converted from RNA to DNA).
        """
        kazusa_frequencies = {
            # UUN → TTN codons (Phe, Leu)
            'TTT': 17.14, 'TTC': 17.48, 'TTA': 8.71, 'TTG': 13.44,
            # UCN → TCN codons (Ser)  
            'TCT': 16.93, 'TCC': 17.32, 'TCA': 14.14, 'TCG': 4.03,
            # UAN → TAN codons (Tyr, stop)
            'TAT': 12.11, 'TAC': 13.49, 'TAA': 0.44, 'TAG': 0.35,
            # UGN → TGN codons (Cys, stop, Trp)
            'TGT': 10.40, 'TGC': 10.81, 'TGA': 0.79, 'TGG': 11.60,
            
            # CUN → CTN codons (Leu)
            'CTT': 14.08, 'CTC': 17.81, 'CTA': 7.44, 'CTG': 36.10,
            # CCN → CCN codons (Pro)
            'CCT': 19.31, 'CCC': 19.11, 'CCA': 18.92, 'CCG': 6.22,
            # CAN → CAN codons (His, Gln)
            'CAT': 11.83, 'CAC': 14.65, 'CAA': 14.06, 'CAG': 35.53,
            # CGN → CGN codons (Arg)
            'CGT': 4.55, 'CGC': 8.71, 'CGA': 6.42, 'CGG': 10.79,
            
            # AUN → ATN codons (Ile, Met)
            'ATT': 16.48, 'ATC': 18.67, 'ATA': 8.08, 'ATG': 21.53,
            # ACN → ACN codons (Thr)
            'ACT': 14.26, 'ACC': 17.85, 'ACA': 16.52, 'ACG': 5.59,
            # AAN → AAN codons (Asn, Lys)
            'AAT': 18.43, 'AAC': 18.30, 'AAA': 27.48, 'AAG': 31.77,
            # AGN → AGN codons (Ser, Arg)
            'AGT': 14.05, 'AGC': 19.69, 'AGA': 13.28, 'AGG': 12.13,
            
            # GUN → GTN codons (Val)
            'GTT': 11.74, 'GTC': 13.44, 'GTA': 7.66, 'GTG': 25.87,
            # GCN → GCN codons (Ala)
            'GCT': 18.99, 'GCC': 25.84, 'GCA': 17.04, 'GCG': 5.91,
            # GAN → GAN codons (Asp, Glu)
            'GAT': 24.03, 'GAC': 24.27, 'GAA': 33.65, 'GAG': 39.67,
            # GGN → GGN codons (Gly)
            'GGT': 10.83, 'GGC': 19.79, 'GGA': 17.12, 'GGG': 15.35
        }
        
        # Remove stop codons for CAI calculation
        return {codon: freq for codon, freq in kazusa_frequencies.items() 
                if codon not in self.stop_codons}
    
    def _get_genetic_code(self):
        """Standard genetic code mapping codons to amino acids."""
        return {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
    
    def _build_codon_weights(self):
        """
        Build codon weights (relative adaptiveness) from reference frequencies.
        Implements Sharp & Li (1987) method: w_ij = f_ij / f_max for each amino acid.
        """
        frequencies = self._get_kazusa_human_frequencies()
        genetic_code = self._get_genetic_code()
        
        # Group codons by amino acid
        aa_codons = defaultdict(list)
        for codon, aa in genetic_code.items():
            if aa != '*':  # Skip stop codons
                aa_codons[aa].append(codon)
        
        # Calculate relative adaptiveness for each codon
        codon_weights = {}
        
        for aa, codons in aa_codons.items():
            # Get frequencies for all codons of this amino acid
            aa_frequencies = {codon: frequencies.get(codon, 0) for codon in codons}
            
            # Find maximum frequency (optimal codon)
            max_freq = max(aa_frequencies.values())
            
            if max_freq > 0:
                # Calculate relative adaptiveness: w_ij = f_ij / f_max
                for codon in codons:
                    codon_weights[codon] = aa_frequencies[codon] / max_freq
            else:
                # If no frequency data, assign equal weights
                for codon in codons:
                    codon_weights[codon] = 1.0
        
        return codon_weights
    
    def validate_and_clean_sequence(self, sequence):
        """
        Validate and clean DNA sequence for CAI calculation.
        Returns (cleaned_sequence, codon_count) or (None, 0) if invalid.
        """
        if pd.isna(sequence) or sequence == '' or sequence is None:
            return None, 0
        
        # Handle 'NA' values
        if str(sequence).upper() == 'NA':
            return None, 0
        
        # Clean sequence: uppercase, remove whitespace and non-ATCG characters
        clean_seq = re.sub(r'[^ATCG]', '', str(sequence).upper().strip())
        
        if len(clean_seq) == 0:
            return None, 0
        
        # Ensure sequence is multiple of 3 (complete codons only)
        remainder = len(clean_seq) % 3
        if remainder != 0:
            print(f"ERROR: Sequence length {len(clean_seq)} is not a multiple of 3.")
            return None, 0
        codon_count = len(clean_seq) // 3
        return clean_seq if codon_count > 0 else None, codon_count
    
    def extract_upstream_codons(self, full_sequence, num_codons):
        """
        Extract the last N codons from a full CDS sequence (closest to PTC/NTC).
        
        Parameters:
        - full_sequence: Full CDS sequence upstream of PTC/NTC
        - num_codons: Number of codons to extract from the end
        
        Returns: (extracted_sequence, original_codon_count, extracted_codon_count)
        """
        clean_seq, total_codons = self.validate_and_clean_sequence(full_sequence)
        
        if clean_seq is None or total_codons == 0:
            return None, 0, 0
        
        if num_codons >= total_codons:
            # Return entire sequence if requesting more codons than available
            return clean_seq, total_codons, total_codons
        
        # Extract last N codons (closest to PTC/NTC)
        start_position = (total_codons - num_codons) * 3
        extracted_seq = clean_seq[start_position:]
        
        return extracted_seq, total_codons, num_codons
    
    def calculate_cai(self, sequence):
        """
        Calculate Codon Adaptation Index for a DNA sequence.
        
        Parameters:
        - sequence: DNA sequence string
        
        Returns: CAI value (float) or np.nan if invalid
        """
        clean_seq, codon_count = self.validate_and_clean_sequence(sequence)
        
        if clean_seq is None or codon_count == 0:
            return float('nan')
        
        # Extract codons
        codons = [clean_seq[i:i+3] for i in range(0, len(clean_seq), 3)]
        
        # Calculate geometric mean of codon weights (excluding stop codons)
        log_sum = 0.0
        valid_codons = 0
        
        for codon in codons:
            # Skip stop codons for CAI calculation
            if codon in self.stop_codons:
                continue
            
            # Get codon weight
            weight = self.codon_weights.get(codon, None)
            
            if weight is None or weight <= 0:
                # Unknown codon or zero weight
                continue
            
            log_sum += math.log(weight)
            valid_codons += 1
        
        if valid_codons == 0:
            return float('nan')
        
        # Calculate CAI = geometric mean = exp(mean of log values)
        cai = math.exp(log_sum / valid_codons)
        
        return cai


class TranscriptAnnotation:
    """Class to store transcript information from GTF"""
    def __init__(self, transcript_id, chrom, strand, biotype=None, gene_id=None):
        self.transcript_id = transcript_id
        self.chrom = chrom
        self.strand = strand
        self.biotype = biotype
        self.gene_id = gene_id
        self.exons = []
        self.cds = []
        self.stop_codon = []
        self.utr5 = []
        self.utr3 = []
    
    def get_cds_length(self):
        """Calculate total CDS length accounting for frame offset"""
        if not self.cds:
            return 0
        
        cds_sorted = sorted(self.cds, key=lambda x: x[0])
        
        if self.strand == '+':
            biological_5prime_idx = 0
        else:
            biological_5prime_idx = len(cds_sorted) - 1
        
        total_length = 0
        for i, (start, end, frame) in enumerate(cds_sorted):
            segment_length = end - start + 1
            
            if i == biological_5prime_idx and frame > 0:
                segment_length -= frame
            
            total_length += segment_length
        
        return total_length
        
    def get_cds_sequence(self, genome_fasta):
        """Extract CDS sequence from genome accounting for frame"""
        if not self.cds:
            return ""
            
        cds_seq = ""
        cds_sorted = sorted(self.cds, key=lambda x: x[0])
        
        if self.strand == '+':
            biological_5prime_idx = 0
        else:
            biological_5prime_idx = len(cds_sorted) - 1
        
        for i, (start, end, frame) in enumerate(cds_sorted):
            seg_seq = genome_fasta.fetch(self.chrom, start-1, end)
            
            if i == biological_5prime_idx and frame > 0:
                if self.strand == '+':
                    seg_seq = seg_seq[frame:]
                else:
                    seg_seq = seg_seq[:-frame] if frame < len(seg_seq) else ""
            
            cds_seq += seg_seq
        
        if self.strand == '-':
            cds_seq = str(Seq(cds_seq).reverse_complement())
        
        return cds_seq
    
    def get_utr3_length(self):
        """Calculate total 3' UTR length"""
        if self.utr3:
            return sum(end - start + 1 for start, end in self.utr3)
        return None
    
    def get_exon_junctions(self):
        """Get positions of exon-exon junctions in CDS coordinates"""
        junctions = []
        cds_sorted = sorted(self.cds, key=lambda x: x[0])
        
        if self.strand == '-':
            cds_sorted = cds_sorted[::-1]
        
        if self.strand == '+':
            biological_5prime_idx = 0
        else:
            biological_5prime_idx = 0
        
        cumulative_length = 0
        for i, (start, end, frame) in enumerate(cds_sorted[:-1]):
            if i == biological_5prime_idx and frame > 0:
                segment_length = (end - start + 1) - frame
            else:
                segment_length = end - start + 1
            
            cumulative_length += segment_length
            junctions.append(cumulative_length)
        
        return junctions

    def get_cds_segments_with_frame(self):
        """Get CDS segments properly ordered for coordinate mapping"""
        cds_sorted = sorted(self.cds, key=lambda x: x[0])
        if self.strand == '-':
            cds_sorted = cds_sorted[::-1]
        return cds_sorted
    
    def build_genomic_to_transcript_map(self):
        """Build mapping from genomic positions to transcript positions (0-based)"""
        genomic_to_transcript = {}
        transcript_pos = 0
        
        # Build transcript regions in 5' to 3' order
        regions = []
        
        # Add 5' UTR
        regions.extend([('5utr', start, end) for start, end in sorted(self.utr5)])
        
        # Add CDS
        regions.extend([('cds', start, end, frame) for start, end, frame in self.cds])
        
        # Add 3' UTR
        regions.extend([('3utr', start, end) for start, end in sorted(self.utr3)])
        
        # Sort by genomic coordinates
        if self.strand == '+':
            regions.sort(key=lambda x: x[1])
        else:
            regions.sort(key=lambda x: x[1], reverse=True)
        
        # Build mapping
        for region_info in regions:
            region_type = region_info[0]
            start = region_info[1]
            end = region_info[2]
            
            if self.strand == '+':
                for genomic_pos in range(start, end + 1):
                    genomic_to_transcript[genomic_pos] = transcript_pos
                    transcript_pos += 1
            else:
                for genomic_pos in range(end, start - 1, -1):
                    genomic_to_transcript[genomic_pos] = transcript_pos
                    transcript_pos += 1
        
        return genomic_to_transcript


def parse_gtf(gtf_file):
    """Parse GTF file to extract transcript annotations with frame information"""
    transcripts = {}
    
    with gzip.open(gtf_file, 'rt') if gtf_file.endswith('.gz') else open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            start, end = int(start), int(end)
            
            if frame == '.':
                frame = 0
            else:
                frame = int(frame)
            
            attr_dict = {}
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr:
                    match = re.search(r'(\w+)\s+"([^"]+)"', attr)
                    if match:
                        attr_dict[match.group(1)] = match.group(2)
            
            if 'transcript_id' not in attr_dict:
                continue
            
            transcript_id = attr_dict['transcript_id']
            gene_id = attr_dict.get('gene_id', None)
            
            if transcript_id not in transcripts:
                biotype = attr_dict.get('transcript_biotype') or attr_dict.get('transcript_type')
                transcripts[transcript_id] = TranscriptAnnotation(transcript_id, chrom, strand, biotype, gene_id)
            
            if feature == 'exon':
                transcripts[transcript_id].exons.append((start, end))
            elif feature == 'CDS':
                transcripts[transcript_id].cds.append((start, end, frame))
            elif feature == 'stop_codon':
                transcripts[transcript_id].stop_codon.append((start, end))
            elif feature == 'five_prime_utr':
                transcripts[transcript_id].utr5.append((start, end))
            elif feature == 'three_prime_utr':
                transcripts[transcript_id].utr3.append((start, end))
            elif feature == 'UTR':
                if strand == '+':
                    if transcripts[transcript_id].cds and start < min(c[0] for c in transcripts[transcript_id].cds):
                        transcripts[transcript_id].utr5.append((start, end))
                    else:
                        transcripts[transcript_id].utr3.append((start, end))
                else:
                    if transcripts[transcript_id].cds and end > max(c[1] for c in transcripts[transcript_id].cds):
                        transcripts[transcript_id].utr5.append((start, end))
                    else:
                        transcripts[transcript_id].utr3.append((start, end))
    
    return transcripts


def load_m6a_sites_by_gene(m6a_file):
    """Load m6A sites and group them by gene ID."""
    print("Loading m6A sites and grouping by gene...")
    m6a_df = pd.read_csv(m6a_file, sep='\t')
    
    # Dictionary to store m6A sites by gene ID
    gene_m6a = defaultdict(list)
    
    processed_sites = 0
    
    for _, row in m6a_df.iterrows():
        processed_sites += 1
        
        # Extract gene ID and genomic position
        gene_id = row.get('Ensembl_ID', '')
        chrom = row['seqnames']
        genomic_pos = int(row['start'])
        
        # Skip if no gene ID
        if not gene_id:
            continue
        
        # Store genomic position for this gene
        gene_m6a[gene_id].append({
            'chromosome': chrom,
            'genomic_pos': genomic_pos
        })
    
    # Sort positions for each gene
    for gene_id in gene_m6a:
        gene_m6a[gene_id].sort(key=lambda x: x['genomic_pos'])
    
    print(f"Processed {processed_sites} m6A sites")
    print(f"Found m6A sites for {len(gene_m6a)} genes")
    
    return gene_m6a


def count_m6a_in_range(m6a_positions, start_pos, end_pos):
    """Count m6A sites within a range (inclusive)."""
    if not m6a_positions:
        return 0
    
    count = 0
    for pos in m6a_positions:
        if start_pos <= pos <= end_pos:
            count += 1
        elif pos > end_pos:
            break  # Since positions are sorted
    
    return count


def load_cds_sequences(cds_fasta_file):
    """
    Load CDS sequences from FASTA file and remove stop codon if present.
    Also extracts just the 3-nt stop codon for each transcript (memory efficient).
    
    Returns:
        cds_sequences: dict of transcript_id -> CDS sequence (without stop)
        stop_codons_dict: dict of transcript_id -> 3-nt stop codon string
    """
    cds_sequences = {}
    stop_codons_dict = {}
    stop_codons = ['TAA', 'TAG', 'TGA']
    
    opener = gzip.open if cds_fasta_file.endswith('.gz') else open
    mode = 'rt' if cds_fasta_file.endswith('.gz') else 'r'
    
    with opener(cds_fasta_file, mode) as f:
        for record in SeqIO.parse(f, "fasta"):
            transcript_id = record.id.split()[0]
            seq = str(record.seq).upper()
            
            # Extract stop codon 
            natural_stop = ''  
            if len(seq) >= 3 and seq[-3:] in stop_codons:
                natural_stop = seq[-3:]
                seq = seq[:-3]
            
            cds_sequences[transcript_id] = seq
            stop_codons_dict[transcript_id] = natural_stop
            
            base_id = transcript_id.split('.')[0]
            if base_id not in cds_sequences:
                cds_sequences[base_id] = seq
                stop_codons_dict[base_id] = natural_stop
    
    print(f"Loaded {len(cds_sequences)} CDS sequences (stop codons removed if present)")
    return cds_sequences, stop_codons_dict


def load_gene_expression(expression_file):
    """
    Load gene expression data from GTEx CSV file.
    
    Args:
        expression_file: Path to CSV with 'Clean_Gene_ID' and 'Mean_Expression' columns
    
    Returns:
        Dictionary mapping gene IDs to mean expression values
    """
    print("Loading gene expression data...")
    expression_df = pd.read_csv(expression_file)
    
    # Check required columns exist
    if 'Clean_Gene_ID' not in expression_df.columns or 'Mean_Expression' not in expression_df.columns:
        raise ValueError("Expression file must have 'Clean_Gene_ID' and 'Mean_Expression' columns")
    
    # Create dictionary mapping gene ID to expression
    gene_expression = dict(zip(expression_df['Clean_Gene_ID'], expression_df['Mean_Expression']))
    
    print(f"Loaded expression data for {len(gene_expression)} genes")
    return gene_expression


def detect_vep_format(vep_file):
    """Detect the format of VEP output file"""
    with open(vep_file, 'r') as f:
        for line in f:
            if line.startswith('##fileformat=VCF'):
                return 'vcf'
            elif line.startswith('#') and not line.startswith('##'):
                if '\tExtra' in line or line.endswith('Extra\n'):
                    return 'extra'
                else:
                    return 'standard'
    return 'standard'


def passes_consequence_filter(consequence_str):
    """Check if consequence contains stop_gained or frameshift_variant"""
    if not consequence_str or consequence_str == '-':
        return False
    return 'stop_gained' in consequence_str or 'frameshift_variant' in consequence_str


def passes_biotype_filter(biotype_str):
    """Check if biotype is protein_coding"""
    if not biotype_str or biotype_str == '-':
        return None  # Unknown, will check GTF later
    return biotype_str == 'protein_coding'


def parse_vcf_csq_format(vep_file):
    """Parse VCF format with CSQ annotations - WITH PRE-FILTERING"""
    records = []
    csq_format = None
    
    print("  Pre-filtering VCF for stop_gained/frameshift_variant consequences...")
    
    with open(vep_file, 'r') as f:
        for line in f:
            if line.startswith('##INFO=<ID=CSQ'):
                match = re.search(r'Format: ([^"]+)', line)
                if match:
                    csq_format = match.group(1).split('|')
            elif line.startswith('#CHROM'):
                continue
            elif line.startswith('#'):
                continue
            else:
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue
                
                chrom, pos, var_id, ref, alt, qual, filt, info = fields[:8]
                
                csq_match = re.search(r'CSQ=([^;\t]+)', info)
                if not csq_match:
                    continue
                
                csq_data = csq_match.group(1)
                
                for annotation in csq_data.split(','):
                    values = annotation.split('|')
                    if csq_format and len(values) == len(csq_format):
                        # PRE-FILTER: Check consequence first
                        consequence_idx = csq_format.index('Consequence') if 'Consequence' in csq_format else None
                        if consequence_idx is not None:
                            consequence = values[consequence_idx]
                            if not passes_consequence_filter(consequence):
                                continue
                        
                        # PRE-FILTER: Check biotype if available
                        biotype_idx = csq_format.index('BIOTYPE') if 'BIOTYPE' in csq_format else None
                        if biotype_idx is not None:
                            biotype = values[biotype_idx]
                            biotype_ok = passes_biotype_filter(biotype)
                            if biotype_ok is False:  # Explicitly non-protein_coding
                                continue
                        
                        record = {
                            'CHR': chrom,
                            'POS': pos,
                            'ID': var_id,
                            'REF_ALLELE': ref,
                            'ALT_ALLELE': alt,
                            'Uploaded_variation': f"{chrom}_{pos}_{ref}/{alt}"
                        }
                        
                        for i, field_name in enumerate(csq_format):
                            record[field_name] = values[i] if values[i] else '-'
                        
                        record['Location'] = f"{chrom}:{pos}"
                        
                        if 'Feature' in record:
                            record['transcript_id'] = record['Feature']
                        
                        records.append(record)
    
    return pd.DataFrame(records) if records else pd.DataFrame()


def parse_extra_column_format(vep_file):
    """Parse tab-delimited format with Extra column - WITH PRE-FILTERING"""
    records = []
    header = None
    
    print("  Pre-filtering for stop_gained/frameshift_variant consequences...")
    
    with open(vep_file, 'r') as f:
        for line in f:
            if line.startswith('#') and not line.startswith('##'):
                header = line.lstrip('#').strip().split('\t')
            elif not line.startswith('#'):
                if header is None:
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) != len(header):
                    continue
                
                record = {}
                for i, field_name in enumerate(header):
                    record[field_name] = fields[i]
                
                # PRE-FILTER: Check consequence
                if 'Consequence' in record:
                    if not passes_consequence_filter(record['Consequence']):
                        continue
                
                # PRE-FILTER: Check biotype if available
                if 'BIOTYPE' in record:
                    biotype_ok = passes_biotype_filter(record['BIOTYPE'])
                    if biotype_ok is False:
                        continue
                
                # Parse Extra column
                if 'Extra' in record and record['Extra'] != '-':
                    extra_pairs = record['Extra'].split(';')
                    for pair in extra_pairs:
                        if '=' in pair:
                            key, value = pair.split('=', 1)
                            record[key.strip()] = value.strip()
                
                records.append(record)
    
    return pd.DataFrame(records) if records else pd.DataFrame()


def parse_standard_format(vep_file):
    """Parse standard tab-delimited format - WITH PRE-FILTERING"""
    records = []
    header = None
    
    print("  Pre-filtering for stop_gained/frameshift_variant consequences...")
    
    with open(vep_file, 'r') as f:
        for line in f:
            if line.startswith('#') and not line.startswith('##'):
                header = line.lstrip('#').strip().split('\t')
            elif not line.startswith('#'):
                if header is None:
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) != len(header):
                    continue
                
                record = {}
                for i, field_name in enumerate(header):
                    record[field_name] = fields[i]
                
                # PRE-FILTER: Check consequence
                if 'Consequence' in record:
                    if not passes_consequence_filter(record['Consequence']):
                        continue
                
                # PRE-FILTER: Check biotype if available
                if 'BIOTYPE' in record:
                    biotype_ok = passes_biotype_filter(record['BIOTYPE'])
                    if biotype_ok is False:
                        continue
                
                records.append(record)
    
    return pd.DataFrame(records) if records else pd.DataFrame()


def find_variant_induced_ptc(original_seq, modified_seq, variant_cds_pos):
    """Find variant-induced PTC in modified sequence"""
    stop_codons = ['TAA', 'TAG', 'TGA']
    
    variant_codon_start = ((variant_cds_pos - 1) // 3) * 3
    for i in range(0, min(variant_codon_start, len(original_seq)-2), 3):
        codon = original_seq[i:i+3]
        if codon in stop_codons:
            print(f"Warning! Stop codon found upstream of variant at position {i+1}!")
            return None
    
    for i in range(variant_codon_start, len(modified_seq) - 2, 3):
        codon = modified_seq[i:i+3]
        
        if codon in stop_codons:
            return i + 1
    
    return None


def apply_variant_to_sequence(ref_seq, variant_pos, ref_allele, alt_allele, suppress_warnings=False):
    """Apply variant to sequence and return modified sequence"""
    pos_0based = variant_pos - 1
    
    if ref_seq[pos_0based:pos_0based+len(ref_allele)] != ref_allele:
        if not suppress_warnings:
            print(f"Warning: Reference allele does not match at position {variant_pos}. Expected {ref_allele}, found {ref_seq[pos_0based:pos_0based+len(ref_allele)]}")
        return None
    
    new_seq = ref_seq[:pos_0based] + alt_allele + ref_seq[pos_0based+len(ref_allele):]
    return new_seq


def calculate_variant_length_change(ref_allele, alt_allele):
    """Calculate the length change introduced by the variant"""
    return len(alt_allele) - len(ref_allele)


def calculate_gc_content(sequence):
    """Calculate GC content percentage"""
    if not sequence:
        return 0
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence))


def map_cds_position_to_genomic(transcript, cds_pos):
    """Map CDS position back to genomic coordinates (accounting for frame)"""
    if not transcript.cds:
        return None
    
    cds_sorted = sorted(transcript.cds, key=lambda x: x[0])
    
    if transcript.strand == '-':
        cds_sorted = cds_sorted[::-1]
    
    if transcript.strand == '+':
        biological_5prime_idx = 0
    else:
        biological_5prime_idx = 0
    
    cumulative_pos = 0
    for i, (start, end, frame) in enumerate(cds_sorted):
        if i == biological_5prime_idx and frame > 0:
            effective_length = (end - start + 1) - frame
            effective_start = start + frame if transcript.strand == '+' else start
        else:
            effective_length = end - start + 1
            effective_start = start
        
        if cumulative_pos < cds_pos <= cumulative_pos + effective_length:
            pos_in_segment = cds_pos - cumulative_pos
            
            if transcript.strand == '+':
                if i == biological_5prime_idx and frame > 0:
                    return start + frame + pos_in_segment - 1
                else:
                    return start + pos_in_segment - 1
            else:
                if i == biological_5prime_idx and frame > 0:
                    return end - frame - pos_in_segment + 1
                else:
                    return end - pos_in_segment + 1
        
        cumulative_pos += effective_length
    
    return None


def map_genomic_position_to_cds(transcript, genomic_pos):
    """Map genomic position to CDS coordinates (accounting for frame)"""
    if not transcript.cds:
        return None
    
    cds_sorted = sorted(transcript.cds, key=lambda x: x[0])
    
    if transcript.strand == '-':
        cds_sorted = cds_sorted[::-1]
    
    if transcript.strand == '+':
        biological_5prime_idx = 0
    else:
        biological_5prime_idx = 0
    
    cumulative_pos = 0
    for i, (start, end, frame) in enumerate(cds_sorted):
        if start <= genomic_pos <= end:
            if i == biological_5prime_idx and frame > 0:
                if transcript.strand == '+':
                    effective_start = start + frame
                    if genomic_pos < effective_start:
                        return None
                    pos_in_segment = genomic_pos - effective_start + 1
                else:
                    effective_end = end - frame
                    if genomic_pos > effective_end:
                        return None
                    pos_in_segment = effective_end - genomic_pos + 1
            else:
                if transcript.strand == '+':
                    pos_in_segment = genomic_pos - start + 1
                else:
                    pos_in_segment = end - genomic_pos + 1
            
            return cumulative_pos + pos_in_segment
        
        if i == biological_5prime_idx and frame > 0:
            segment_length = (end - start + 1) - frame
        else:
            segment_length = end - start + 1
        cumulative_pos += segment_length
    
    return None


def analyze_ptc(transcript, ptc_cds_pos, variant_cds_pos, variant_length_change, 
               original_cds_length, modified_cds_seq, original_cds_seq, 
               cai_calculator=None, num_codons=50, min_codons=10,
               gene_m6a=None, ptc_genomic_pos=None):
    """
    Analyze PTC and calculate all features including CAI, m6A (with length-normalized), and downstream AUG.
    """
    features = {}
    
    modified_cds_length = original_cds_length + variant_length_change
    
    junctions = transcript.get_exon_junctions()
    if junctions:
        last_junction = junctions[-1]
        if variant_cds_pos <= last_junction:
            adjusted_last_junction = last_junction + variant_length_change
        else:
            adjusted_last_junction = last_junction
        
        features['50nt_rule'] = (adjusted_last_junction - ptc_cds_pos) <= 50
    else:
        features['50nt_rule'] = True
    
    features['CDS_position'] = ptc_cds_pos
    features['variant_CDS_POS'] = variant_cds_pos
    
    features['distance_to_stop'] = modified_cds_length - ptc_cds_pos + 1
    
    features['gc_content'] = calculate_gc_content(modified_cds_seq)
    
    utr3_length = transcript.get_utr3_length()
    if utr3_length is not None:
        features['dis_to_3utr_end'] = features['distance_to_stop'] + utr3_length
    else:
        features['dis_to_3utr_end'] = 'NA'

    if ptc_cds_pos > variant_cds_pos:
        original_ptc_cds_pos = ptc_cds_pos - variant_length_change
    else:
        original_ptc_cds_pos = ptc_cds_pos
    
    ptc_genomic_pos_calc = map_cds_position_to_genomic(transcript, original_ptc_cds_pos)
    variant_genomic_pos = map_cds_position_to_genomic(transcript, variant_cds_pos)

    if ptc_genomic_pos_calc:
        exons_sorted = sorted(transcript.exons) if transcript.exons else [(s, e) for s, e, f in sorted(transcript.cds, key=lambda x: x[0])]
        if transcript.strand == '-':
            exons_sorted = exons_sorted[::-1]
        
        exon_idx = None
        actual_exon = None
        variant_in_same_exon = False
        for i, (start, end) in enumerate(exons_sorted):
            if start <= ptc_genomic_pos_calc <= end:
                exon_idx = i
                actual_exon = (start, end)
                if variant_genomic_pos and start <= variant_genomic_pos <= end:
                    variant_in_same_exon = True
                break
        
        if exon_idx is not None and actual_exon:
            features['downstream_exons'] = len(exons_sorted) - exon_idx - 1
            features['upstream_exons'] = exon_idx
            base_exon_length = actual_exon[1] - actual_exon[0] + 1
            if variant_in_same_exon:
                features['exon_length'] = base_exon_length + variant_length_change
            else:
                features['exon_length'] = base_exon_length
            
            if transcript.strand == '+':
                features['dis_to_exon_end'] = actual_exon[1] - ptc_genomic_pos_calc + 1
            else:
                features['dis_to_exon_end'] = ptc_genomic_pos_calc - actual_exon[0] + 1
        else:
            features['downstream_exons'] = 0
            features['upstream_exons'] = 0
            features['exon_length'] = 0
            features['dis_to_exon_end'] = 0
    else:
        features['downstream_exons'] = 0
        features['upstream_exons'] = 0
        features['exon_length'] = 0
        features['dis_to_exon_end'] = 0

    features['PTC_POS'] = ptc_genomic_pos_calc if ptc_genomic_pos_calc is not None else 'NA'
    
    # ===== CAI CALCULATION =====
    if cai_calculator is not None:
        cai_col_suffix = f'{num_codons}codon'
        features[f'CAI_{cai_col_suffix}_upstream_PTC'] = 'NA'
        features[f'CAI_{cai_col_suffix}_upstream_NTC'] = 'NA'
        features[f'CAI_{cai_col_suffix}_upstream_diff'] = 'NA'
        features['upstream_PTC_length_codons'] = 0
        features['upstream_NTC_length_codons'] = 0
        
        ptc_pos_0based = ptc_cds_pos - 1
        upstream_ptc_seq = modified_cds_seq[:ptc_pos_0based]
        
        if len(upstream_ptc_seq) >= min_codons * 3:
            extracted_ptc, total_ptc_codons, extracted_ptc_codons = cai_calculator.extract_upstream_codons(
                upstream_ptc_seq, num_codons
            )
            
            features['upstream_PTC_length_codons'] = total_ptc_codons
            
            if extracted_ptc is not None and extracted_ptc_codons >= min_codons:
                cai_ptc = cai_calculator.calculate_cai(extracted_ptc)
                if not math.isnan(cai_ptc):
                    features[f'CAI_{cai_col_suffix}_upstream_PTC'] = cai_ptc
        
        upstream_ntc_seq = original_cds_seq
        
        if len(upstream_ntc_seq) >= min_codons * 3:
            extracted_ntc, total_ntc_codons, extracted_ntc_codons = cai_calculator.extract_upstream_codons(
                upstream_ntc_seq, num_codons
            )
            
            features['upstream_NTC_length_codons'] = total_ntc_codons
            
            if extracted_ntc is not None and extracted_ntc_codons >= min_codons:
                cai_ntc = cai_calculator.calculate_cai(extracted_ntc)
                if not math.isnan(cai_ntc):
                    features[f'CAI_{cai_col_suffix}_upstream_NTC'] = cai_ntc
        
        cai_ptc_val = features[f'CAI_{cai_col_suffix}_upstream_PTC']
        cai_ntc_val = features[f'CAI_{cai_col_suffix}_upstream_NTC']
        
        if cai_ptc_val != 'NA' and cai_ntc_val != 'NA':
            features[f'CAI_{cai_col_suffix}_upstream_diff'] = cai_ntc_val - cai_ptc_val
    
    # ===== m6A CALCULATION =====
    features['m6A_CDS'] = 0
    features['m6A_all_transcript'] = 0
    features['m6A_CDS_length_normalized_unconstrained'] = 0
    features['m6A_all_length_normalized_unconstrained'] = 0
    
    if gene_m6a is not None and transcript.gene_id and transcript.gene_id in gene_m6a:
        # Build genomic to transcript coordinate mapping
        genomic_to_transcript = transcript.build_genomic_to_transcript_map()
        
        # Calculate total transcript length
        transcript_length = len(genomic_to_transcript) if genomic_to_transcript else 0
        
        # Map m6A sites to transcript coordinates
        m6a_transcript_positions = []
        chrom_with_chr = f"chr{transcript.chrom}" if not transcript.chrom.startswith('chr') else transcript.chrom
        chrom_without_chr = transcript.chrom[3:] if transcript.chrom.startswith('chr') else transcript.chrom
        
        for site in gene_m6a[transcript.gene_id]:
            site_chrom = site['chromosome']
            # Handle chromosome naming variations
            if site_chrom == chrom_with_chr or site_chrom == chrom_without_chr or \
               site_chrom == transcript.chrom:
                transcript_pos = genomic_to_transcript.get(site['genomic_pos'])
                if transcript_pos is not None:
                    m6a_transcript_positions.append(transcript_pos)
        
        m6a_transcript_positions.sort()
        
        # Calculate m6A in CDS (from start codon to PTC)
        if ptc_genomic_pos is not None:
            ptc_transcript_pos = genomic_to_transcript.get(ptc_genomic_pos)
            
            if ptc_transcript_pos is not None:
                # Find start of CDS in transcript coordinates
                cds_start_transcript = None
                if transcript.cds:
                    cds_sorted = sorted(transcript.cds, key=lambda x: x[0])
                    first_cds_start = cds_sorted[0][0]
                    if transcript.strand == '+':
                        cds_start_genomic = first_cds_start
                    else:
                        cds_start_genomic = cds_sorted[-1][1]
                    
                    cds_start_transcript = genomic_to_transcript.get(cds_start_genomic)
                
                if cds_start_transcript is not None:
                    # Count m6A from CDS start to PTC 
                    features['m6A_CDS'] = count_m6a_in_range(
                        m6a_transcript_positions, 
                        cds_start_transcript, 
                        ptc_transcript_pos - 1
                    )
                    
                    # Calculate length-normalized m6A_CDS
                    # Normalize by CDS position of PTC (length from CDS start to PTC)
                    if features['m6A_CDS'] != 'NA' and ptc_cds_pos > 0:
                        features['m6A_CDS_length_normalized_unconstrained'] = features['m6A_CDS'] / ptc_cds_pos
        
        # Calculate m6A in entire transcript
        if m6a_transcript_positions:
            features['m6A_all_transcript'] = len(m6a_transcript_positions)
            
            # Calculate length-normalized m6A_all_transcript
            # Normalize by total transcript length
            if transcript_length > 0:
                features['m6A_all_length_normalized_unconstrained'] = features['m6A_all_transcript'] / transcript_length
    
    # ===== DOWNSTREAM AUG CALCULATION =====
    # Use 100000 as fill value for cases where no downstream AUG is found
    features['dis_to_first_inframeAUG'] = 100000
    features['dis_to_first_outframeAUG'] = 100000
    features['has_downstream_inframeAUG'] = False
    features['has_downstream_outframeAUG'] = False
    
    inframe_dist, outframe_dist = find_downstream_augs(modified_cds_seq, ptc_cds_pos)
    
    if inframe_dist is not None:
        features['dis_to_first_inframeAUG'] = inframe_dist
        features['has_downstream_inframeAUG'] = True
    
    if outframe_dist is not None:
        features['dis_to_first_outframeAUG'] = outframe_dist
        features['has_downstream_outframeAUG'] = True
    
    return features


def process_vep_line(row, transcripts, genome_fasta, cds_sequences=None, stop_codons_dict=None,
                    force_alt=False, suppress_warnings=False, cai_calculator=None, num_codons=50, 
                    min_codons=10, gene_m6a=None, af_column=None, store_sequence=False, skip_check=False,
                    warnings_list=None):
    """
    Process a single VEP output line and return PTC analysis results.
    If store_sequence=True, stores modified_sequence WITH natural stop codon for TranslationAI FASTA output.
    If skip_check=True, skips PTC creation verification and directly calculates features.
    If warnings_list is provided, all validation failures are logged to it.
    If af_column is provided, extracts AF from that column.
    """
    results = []
    
    def log_warning(msg, variant_info=""):
        """Helper to log warnings to the warnings list"""
        if warnings_list is not None:
            if variant_info:
                warnings_list.append(f"[{variant_info}] {msg}")
            else:
                warnings_list.append(msg)
    
    def safe_extract_value(row, column_name, default=''):
        """Safely extract a scalar value from a pandas Series row or namedtuple"""
        try:
            if hasattr(row, column_name):
                val = getattr(row, column_name)
            elif hasattr(row, 'get'):
                val = row.get(column_name, default)
            else:
                return default
            
            if val is None or (hasattr(val, '__iter__') and not isinstance(val, str)):
                return default
            str_val = str(val)
            if str_val in ['nan', 'NaN', 'None', '<NA>']:
                return default
            return str_val
        except:
            return default
    
    location = safe_extract_value(row, 'Location')
    allele = safe_extract_value(row, 'Allele') or safe_extract_value(row, 'ALT_ALLELE')
    
    transcript_id = safe_extract_value(row, 'Feature') or safe_extract_value(row, 'transcript_id')
    
    uploaded_variation = safe_extract_value(row, 'Uploaded_variation')
    cds_position_str = safe_extract_value(row, 'CDS_position')
    consequence = safe_extract_value(row, 'Consequence')
    biotype = safe_extract_value(row, 'BIOTYPE')
    canonical = safe_extract_value(row, 'CANONICAL')
    
    # Extract AF from specified column if available
    af_value = None
    if af_column:
        af_value = safe_extract_value(row, af_column)
    
    # Create variant identifier for logging
    variant_id = f"{location}_{transcript_id}_{allele}" if location and transcript_id and allele else "unknown_variant"
    
    if not location or not allele or not transcript_id or not cds_position_str:
        missing_fields = []
        if not location: missing_fields.append("Location")
        if not allele: missing_fields.append("Allele")
        if not transcript_id: missing_fields.append("transcript_id")
        if not cds_position_str: missing_fields.append("CDS_position")
        log_warning(f"Missing required fields: {', '.join(missing_fields)}", variant_id)
        return results

    var_cds_pos = None
    if cds_position_str and cds_position_str != '-' and cds_position_str != 'NA':
        try:
            if '-' in cds_position_str and cds_position_str != '-':
                var_cds_pos = int(cds_position_str.split('-')[0])
            elif cds_position_str.isdigit():
                var_cds_pos = int(cds_position_str)
            else:
                match = re.match(r'(\d+)', cds_position_str)
                if match:
                    var_cds_pos = int(match.group(1))
        except (ValueError, IndexError) as e:
            log_warning(f"Failed to parse CDS position '{cds_position_str}': {e}", variant_id)
            return results
    
    if var_cds_pos is None:
        log_warning(f"Invalid CDS position: '{cds_position_str}'", variant_id)
        return results
    
    try:
        if ':' in location:
            chrom, pos_str = location.split(':')
            if '-' in pos_str:
                pos = int(pos_str.split('-')[0])
            else:
                pos = int(pos_str)
        else:
            log_warning(f"Invalid location format: '{location}'", variant_id)
            return results
    except (ValueError, IndexError) as e:
        log_warning(f"Failed to parse location '{location}': {e}", variant_id)
        return results
    
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    
    ref_allele = safe_extract_value(row, 'REF_ALLELE') or safe_extract_value(row, 'REF')
    alt_allele = allele
    
    if not ref_allele and uploaded_variation and '/' in uploaded_variation:
        try:
            parts = uploaded_variation.split('_')
            if len(parts) >= 3:
                alleles_part = parts[-1]
                if '/' in alleles_part:
                    alleles = alleles_part.split('/')
                    if len(alleles) == 2:
                        ref_allele = alleles[0]
        except:
            pass
    
    if ref_allele == '-':
        ref_allele = ''
    if alt_allele == '-':
        alt_allele = ''
    
    transcript_with_version = None
    transcript_id_base = transcript_id.split('.')[0]
    
    if transcript_id in transcripts:
        transcript_with_version = transcript_id
    elif transcript_id_base in transcripts:
        transcript_with_version = transcript_id_base
    else:
        for gtf_transcript_id in transcripts:
            if gtf_transcript_id.startswith(transcript_id_base + '.'):
                transcript_with_version = gtf_transcript_id
                break
    
    if not transcript_with_version:
        log_warning(f"Transcript '{transcript_id}' not found in GTF", variant_id)
        return results
    
    transcript = transcripts[transcript_with_version]
    
    effective_biotype = None
    if biotype and biotype != '-':
        effective_biotype = biotype
    elif transcript.biotype:
        effective_biotype = transcript.biotype
    
    if effective_biotype and effective_biotype != 'protein_coding':
        if not skip_check:
            log_warning(f"Biotype '{effective_biotype}' is not protein_coding", variant_id)
            return results
        else:
            log_warning(f"Biotype '{effective_biotype}' is not protein_coding (processing anyway due to --skip-check)", variant_id)
    
    transcript_chrom = transcript.chrom
    if transcript_chrom.startswith('chr'):
        transcript_chrom = transcript_chrom[3:]
    
    if transcript_chrom != chrom:
        if not skip_check:
            log_warning(f"Chromosome mismatch: variant on '{chrom}', transcript on '{transcript_chrom}'", variant_id)
            return results
        else:
            log_warning(f"Chromosome mismatch: variant on '{chrom}', transcript on '{transcript_chrom}' (processing anyway due to --skip-check)", variant_id)
    
    original_cds_seq = None
    natural_stop_codon = None  # Will store the 3-nt natural stop codon for TranslationAI
    
    if cds_sequences:
        cds_lookup_id = None
        if transcript_id in cds_sequences:
            cds_lookup_id = transcript_id
        elif transcript_id_base in cds_sequences:
            cds_lookup_id = transcript_id_base
        else:
            for fasta_id in cds_sequences:
                if fasta_id.startswith(transcript_id_base):
                    cds_lookup_id = fasta_id
                    break
        
        if cds_lookup_id:
            original_cds_seq = cds_sequences[cds_lookup_id]
            
            # Get natural stop codon from dictionary 
            if store_sequence and stop_codons_dict and cds_lookup_id in stop_codons_dict:
                natural_stop_codon = stop_codons_dict[cds_lookup_id]
    
    if original_cds_seq is None:
        seq_from_genome = transcript.get_cds_sequence(genome_fasta)
        
        # Always remove stop codon from sequence for analysis
        if len(seq_from_genome) >= 3:
            stop_codons = ['TAA', 'TAG', 'TGA']
            potential_stop = seq_from_genome[-3:]
            if potential_stop in stop_codons:
                if store_sequence:
                    natural_stop_codon = potential_stop
                original_cds_seq = seq_from_genome[:-3]
            else:
                original_cds_seq = seq_from_genome
                # No stop codon found, use default TAA for TranslationAI if needed
                if store_sequence:
                    natural_stop_codon = 'TAA'
        else:
            original_cds_seq = seq_from_genome
    
    original_cds_length = len(original_cds_seq)
    
    if var_cds_pos > original_cds_length:
        if not skip_check:
            log_warning(f"CDS position {var_cds_pos} exceeds CDS length {original_cds_length}", variant_id)
            return results
        else:
            log_warning(f"CDS position {var_cds_pos} exceeds CDS length {original_cds_length} (processing anyway due to --skip-check)", variant_id)
    
    variant_length_change = calculate_variant_length_change(ref_allele, alt_allele)
    
    ref_mismatch = False
    ref_mismatch_msg = None
    try:
        if transcript.strand == '-':
            ref_rc = str(Seq(ref_allele).reverse_complement()) if ref_allele else ''
            alt_rc = str(Seq(alt_allele).reverse_complement()) if alt_allele else ''
            modified_seq = apply_variant_to_sequence(original_cds_seq, var_cds_pos, ref_rc, alt_rc, suppress_warnings=True)
        else:
            modified_seq = apply_variant_to_sequence(original_cds_seq, var_cds_pos, ref_allele, alt_allele, suppress_warnings=True)
        
        if modified_seq is None:
            ref_mismatch = True
            pos_0based = var_cds_pos - 1
            if transcript.strand == '-':
                ref_rc = str(Seq(ref_allele).reverse_complement()) if ref_allele else ''
                expected = ref_rc
                found = original_cds_seq[pos_0based:pos_0based+len(ref_rc)]
            else:
                expected = ref_allele
                found = original_cds_seq[pos_0based:pos_0based+len(ref_allele)]
            
            ref_mismatch_msg = f"Reference allele mismatch at CDS position {var_cds_pos} in {transcript_with_version}. Expected '{expected}', found '{found}'"
            
            if force_alt or skip_check:
                if skip_check and not force_alt:
                    log_warning(f"{ref_mismatch_msg} (forcing ALT due to --skip-check)", variant_id)
                else:
                    log_warning(ref_mismatch_msg, variant_id)
                    
                if transcript.strand == '-':
                    ref_rc = str(Seq(ref_allele).reverse_complement()) if ref_allele else ''
                    alt_rc = str(Seq(alt_allele).reverse_complement()) if alt_allele else ''
                    pos_0based = var_cds_pos - 1
                    modified_seq = original_cds_seq[:pos_0based] + alt_rc + original_cds_seq[pos_0based+len(ref_rc):]
                else:
                    pos_0based = var_cds_pos - 1
                    modified_seq = original_cds_seq[:pos_0based] + alt_allele + original_cds_seq[pos_0based+len(ref_allele):]
            else:
                log_warning(ref_mismatch_msg, variant_id)
                return results
    except Exception as e:
        log_warning(f"Exception during sequence modification: {e}", variant_id)
        if not skip_check:
            return results
        else:
            log_warning(f"Exception during sequence modification: {e} (continuing due to --skip-check, may fail later)", variant_id)
    
    # If we still don't have a modified sequence, we cannot proceed
    if modified_seq is None:
        log_warning(f"Failed to create modified sequence (modified_seq is None)", variant_id)
        return results
    
    # Check for PTC creation or skip check if requested
    if skip_check:
        variant_codon_start = ((var_cds_pos - 1) // 3) * 3
        ptc_pos = variant_codon_start + 1  # 1-based position of the PTC
        log_warning(f"PTC check skipped (--skip-check enabled), assuming PTC at position {ptc_pos}", variant_id)
    else:
        ptc_pos = find_variant_induced_ptc(original_cds_seq, modified_seq, var_cds_pos)
        if not ptc_pos:
            log_warning(f"No PTC found in modified sequence despite stop_gained/frameshift annotation", variant_id)
    
    if ptc_pos:
        try:
            features = analyze_ptc(transcript, ptc_pos, var_cds_pos, variant_length_change, 
                                 original_cds_length, modified_seq, original_cds_seq,
                                 cai_calculator=cai_calculator, num_codons=num_codons, 
                                 min_codons=min_codons, gene_m6a=gene_m6a, ptc_genomic_pos=pos)
            
            result = {
                'CHR': chrom,
                'POS': pos,
                'REF_ALLELE': ref_allele,
                'ALT_ALLELE': alt_allele,
                'Strand': transcript.strand,
                'variant_id': f"{chrom}:{pos}:{ref_allele}>{alt_allele}",
                'transcript_id': transcript_with_version,
                'gene_id': transcript.gene_id if transcript.gene_id else 'NA',
                'is_canonical': canonical == 'YES',
                '_ref_mismatch': ref_mismatch,
                '_ref_mismatch_msg': ref_mismatch_msg,
                **features
            }
            
            # Only store modified_sequence if needed for FASTA output
            # For TranslationAI, append the natural stop codon to the modified sequence
            if store_sequence:
                seq_for_fasta = modified_seq
                if natural_stop_codon:
                    seq_for_fasta = modified_seq + natural_stop_codon
                result['modified_sequence'] = seq_for_fasta
            
            # Only add AF if af_column was specified/detected
            if af_column:
                result['AF'] = af_value if af_value and af_value != '-' else 'NA'
            
            results.append(result)
        except Exception as e:
            log_warning(f"Failed to calculate features: {e}", variant_id)
            pass
    
    return results


def filter_transcripts_by_priority(results, transcripts):
    """Filter transcripts based on priority AFTER checking for PTCs"""
    variant_groups = defaultdict(list)
    for result in results:
        variant_groups[result['variant_id']].append(result)
    
    filtered_results = []
    
    for variant_id, variant_results in variant_groups.items():
        if len(variant_results) == 1:
            filtered_results.extend(variant_results)
            continue
        
        canonical_results = [r for r in variant_results if r['is_canonical']]
        
        if canonical_results:
            filtered_results.extend(canonical_results)
        else:
            max_length = 0
            longest_result = None
            
            for result in variant_results:
                transcript_id = result['transcript_id']
                if transcript_id in transcripts:
                    cds_length = transcripts[transcript_id].get_cds_length()
                    if cds_length > max_length:
                        max_length = cds_length
                        longest_result = result
            
            if longest_result:
                filtered_results.append(longest_result)
            else:
                filtered_results.append(variant_results[0])
    
    return filtered_results


def write_fasta_output(results_df, fasta_file):
    """
    Write results to FASTA format for TranslationAI.
    
    Header format: >chr1:2304-2304(+)(ENST1000001)(567, 456,)
    Where: (downstream_start_pos, ptc_pos,) - both 0-based
    Note: downstream_start_pos = 0 if no start codon found
    
    The sequences include:
    - Modified CDS with the PTC-causing variant applied
    - Original natural stop codon at the end 
    
    Requires 'modified_sequence' column in results_df.
    """
    print(f"\nWriting TranslationAI FASTA output to {fasta_file}...")
    
    # Check if modified_sequence column exists
    if 'modified_sequence' not in results_df.columns:
        print("Error: modified_sequence column not found. This is a programming error.")
        print("The --translationAI-fasta option should ensure sequences are stored.")
        return
    
    written_count = 0
    skipped_count = 0
    
    with open(fasta_file, 'w') as f:
        for _, row in results_df.iterrows():
            # Get sequence - skip if missing
            sequence = row.get('modified_sequence', '')
            if not sequence or sequence == 'NA':
                skipped_count += 1
                continue
            
            # Extract required fields
            chr_name = f"chr{row['CHR']}" if not str(row['CHR']).startswith('chr') else str(row['CHR'])
            pos = int(row['POS'])
            strand = row['Strand']
            transcript_id = row['transcript_id']
            
            # Calculate 0-based PTC position
            ptc_pos = int(row['CDS_position'])
            ptc_pos_0based = ptc_pos - 1
            
            # Calculate downstream start codon position (0-based)
            # If inframe AUG exists and distance is not 100000, calculate absolute position
            if row['has_downstream_inframeAUG'] and row['dis_to_first_inframeAUG'] != 100000:
                downstream_start_pos = ptc_pos_0based + int(row['dis_to_first_inframeAUG'])
            else:
                downstream_start_pos = 0
            
            # Write FASTA header
            header = f">{chr_name}:{pos}-{pos}({strand})({transcript_id})({downstream_start_pos}, {ptc_pos_0based},)"
            f.write(header + '\n')
            f.write(sequence + '\n')
            
            written_count += 1
    
    print(f"TranslationAI FASTA output complete!")
    print(f"  Written: {written_count} sequences to {fasta_file}")
    if skipped_count > 0:
        print(f"  Skipped: {skipped_count} sequences (missing sequence)")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze VEP output for PTCs with CAI, m6A, downstream AUG, gene expression, AF, and optional TranslationAI FASTA output (OPTIMIZED for large files)',
        epilog="""
Note: dis_to_first_inframeAUG and dis_to_first_outframeAUG are set to 100000 when no downstream AUG is found.
        """
    )
    parser.add_argument('vep_file', help='VEP output file (VCF or tab-delimited)')
    parser.add_argument('gtf', help='Reference GTF file')
    parser.add_argument('fasta', help='Reference genome FASTA file')
    parser.add_argument('--cds-fasta', help='CDS sequences FASTA file (optional)')
    parser.add_argument('--m6a-file', help='m6A sites file (optional, TSV format with Ensembl_ID, seqnames, start columns)')
    parser.add_argument('--expression-file', help='GTEx gene expression file (optional, CSV with Clean_Gene_ID and Mean_Expression columns)')
    parser.add_argument('-o', '--output', default='ptc_analysis.tsv', help='Output TSV file')
    parser.add_argument('--warning-log', help='Warning log file (default: <output>.warnings.log)')
    parser.add_argument('--translationAI-fasta', help='Output FASTA file for TranslationAI (optional)')
    parser.add_argument('--canonical', action='store_true',
                    help='Apply transcript prioritization: canonical > longest')
    parser.add_argument('--force-alt', action='store_true',
                    help='Force apply alternative allele even when reference allele mismatches')
    parser.add_argument('--num-codons', type=int, default=25,
                    help='Number of codons to extract for CAI calculation (default: 25)')
    parser.add_argument('--min-codons', type=int, default=10,
                    help='Minimum codons required for CAI calculation (default: 10)')
    parser.add_argument('--skip-cai', action='store_true',
                    help='Skip CAI calculation (faster processing)')
    parser.add_argument('--skip-check', action='store_true',
                    help='Skip validation checks: biotype, chromosome match, CDS bounds, ref allele match, and PTC sequence verification. Forces processing of all variants that pass basic parsing.')
    parser.add_argument('--af-col', help='Column name to use for allele frequency (default: auto-detect gnomAD_AF or gnomADg_AF)')
    
    args = parser.parse_args()
    
    # Set up warning log file
    if args.warning_log:
        warning_log_path = args.warning_log
    else:
        # Default: output_file.warnings.log
        import os
        base_name = os.path.splitext(args.output)[0]
        warning_log_path = f"{base_name}.warnings.log"
    
    if args.num_codons <= 0:
        print("Error: --num-codons must be positive", file=sys.stderr)
        sys.exit(1)
    
    print("Detecting VEP output format...")
    vep_format = detect_vep_format(args.vep_file)
    print(f"Detected format: {vep_format}")
    
    print("Loading GTF annotations...")
    transcripts = parse_gtf(args.gtf)
    print(f"Loaded {len(transcripts)} transcripts")
    
    print("Loading genome FASTA...")
    genome_fasta = pysam.FastaFile(args.fasta)
    
    cds_sequences = None
    stop_codons_dict = None
    if args.cds_fasta:
        print("Loading CDS sequences...")
        cds_sequences, stop_codons_dict = load_cds_sequences(args.cds_fasta)
        print("CDS FASTA loaded - sequences and natural stop codons stored separately")
    else:
        print("No CDS FASTA provided - will extract CDS from genomic coordinates")
    
    # Load m6A sites if provided
    gene_m6a = None
    if args.m6a_file:
        gene_m6a = load_m6a_sites_by_gene(args.m6a_file)
    else:
        print("No m6A file provided - m6A columns (including normalized features) will be 'NA'")
    
    # Load gene expression data if provided
    gene_expression = None
    if args.expression_file:
        gene_expression = load_gene_expression(args.expression_file)
    else:
        print("No expression file provided - Mean_Expression column will be 'NA'")
    
    # Initialize CAI calculator
    cai_calculator = None
    if not args.skip_cai:
        print("Initializing CAI calculator...")
        cai_calculator = StandaloneCAICalculator()
        print(f"CAI will be calculated for {args.num_codons} codons upstream of PTC/NTC")
        print(f"Minimum {args.min_codons} codons required for CAI calculation")
    else:
        print("CAI calculation disabled (--skip-cai)")
    
    # Check if PTC verification should be skipped
    if args.skip_check:
        print("WARNING: Validation checks are RELAXED (--skip-check)")
        print("         Skipping: biotype check, chromosome match, CDS bounds, ref allele match, PTC verification")
        print("         All variants passing basic parsing will be processed")
    
    print("Reading and pre-filtering VEP output file...")
    
    # Parse with pre-filtering
    if vep_format == 'vcf':
        vep_df = parse_vcf_csq_format(args.vep_file)
    elif vep_format == 'extra':
        vep_df = parse_extra_column_format(args.vep_file)
    else:
        vep_df = parse_standard_format(args.vep_file)
    
    if vep_df.empty:
        print("No relevant data found after pre-filtering")
        sys.exit(1)
    
    print(f"After pre-filtering: {len(vep_df)} annotations to process")
    
    # Determine which AF column to use
    af_column = None
    if args.af_col:
        # User specified a column
        if args.af_col in vep_df.columns:
            af_column = args.af_col
            print(f"Using user-specified AF column: '{af_column}'")
        else:
            print(f"WARNING: Specified AF column '{args.af_col}' not found in VEP output")
            print(f"Available columns: {list(vep_df.columns)}")
            print("AF column will not be included")
    else:
        # Auto-detect: try gnomAD_AF first, then gnomADg_AF
        if 'gnomAD_AF' in vep_df.columns:
            af_column = 'gnomAD_AF'
            print(f"Auto-detected AF column: '{af_column}'")
        elif 'gnomADg_AF' in vep_df.columns:
            af_column = 'gnomADg_AF'
            print(f"Auto-detected AF column: '{af_column}'")
        else:
            print("No AF column found (tried 'gnomAD_AF' and 'gnomADg_AF')")
            print("AF column will not be included in output")
    
    has_af_col = af_column is not None
    
    # Check required columns
    required_cols = ['Location', 'CDS_position', 'Consequence']
    if args.canonical:
        required_cols.append('CANONICAL')
    
    has_transcript_col = any(col in vep_df.columns for col in ['Feature', 'transcript_id'])
    has_allele_col = any(col in vep_df.columns for col in ['Allele', 'ALT_ALLELE'])
    
    missing_cols = [col for col in required_cols if col not in vep_df.columns]
    
    if missing_cols or not has_transcript_col or not has_allele_col:
        print(f"Error: Required columns missing or not found")
        print(f"Available columns: {list(vep_df.columns)}")
        sys.exit(1)
    
    has_biotype = 'BIOTYPE' in vep_df.columns
    if not has_biotype:
        print("Warning: BIOTYPE column not found in VEP output.")
        print("         Will use transcript_biotype from GTF file for filtering.")
    
    print(f"Processing {len(vep_df)} pre-filtered annotations...")
    
    all_results = []
    all_warnings = []  # Collect all warnings for logging
    
    # Determine if need to store sequences for FASTA output
    store_sequences = args.translationAI_fasta is not None
    
    if store_sequences:
        print(f"TranslationAI FASTA output enabled: sequences will include natural stop codons")
    
    # Use itertuples for faster iteration
    for idx, row in enumerate(vep_df.itertuples(index=False), 1):
        if idx % 10000 == 0:
            print(f"  Processed {idx}/{len(vep_df)} annotations...")
        
        results = process_vep_line(row, transcripts, genome_fasta, cds_sequences, stop_codons_dict,
                                  args.force_alt, suppress_warnings=True, cai_calculator=cai_calculator, 
                                  num_codons=args.num_codons, min_codons=args.min_codons,
                                  gene_m6a=gene_m6a, af_column=af_column, 
                                  store_sequence=store_sequences, skip_check=args.skip_check,
                                  warnings_list=all_warnings)
        all_results.extend(results)
    
    print(f"Processed {len(vep_df)} annotations total")
    print(f"Found {len(all_results)} PTC-introducing variant-transcript pairs")
    
    # Write all warnings to log file
    if all_warnings:
        print(f"\nWriting {len(all_warnings)} warnings to {warning_log_path}...")
        with open(warning_log_path, 'w') as log_f:
            log_f.write(f"# Warning log for PTC analysis\n")
            log_f.write(f"# Total warnings: {len(all_warnings)}\n")
            log_f.write(f"# Command: {' '.join(sys.argv)}\n\n")
            for warning in all_warnings:
                log_f.write(f"{warning}\n")
        print(f"Warnings written to {warning_log_path}")
    else:
        print(f"\nNo warnings generated - all variants passed validation")
    
    if args.canonical and len(all_results) > 0:
        print("Applying transcript prioritization (canonical > longest)...")
        all_results = filter_transcripts_by_priority(all_results, transcripts)
        print(f"After filtering: {len(all_results)} variant-transcript pairs")
    
    warnings_found = False
    for result in all_results:
        if result.get('_ref_mismatch', False) and result.get('_ref_mismatch_msg'):
            if not warnings_found:
                print("\nWarnings for kept variant-transcript pairs:")
                warnings_found = True
            print(f"  {result['_ref_mismatch_msg']}")
    
    if warnings_found and args.force_alt:
        print("Note: --force-alt was used, so these variants were processed despite reference mismatch\n")
    
    if all_results:
        df = pd.DataFrame(all_results)
        
        # Add gene expression data if available
        if gene_expression is not None:
            print("\nMerging gene expression data...")
            # Map expression to gene_id
            df['Mean_Expression'] = df['gene_id'].map(gene_expression)
            
            # Report merge statistics
            matched_genes = df['Mean_Expression'].notna().sum()
            total_genes = len(df)
            print(f"Matched expression data: {matched_genes}/{total_genes} variants ({100*matched_genes/total_genes:.1f}%)")
        else:
            df['Mean_Expression'] = 'NA'
        
        # Remove internal columns before saving TSV
        # Always drop modified_sequence from TSV - it's only needed for FASTA
        columns_to_drop = ['is_canonical', '_ref_mismatch', '_ref_mismatch_msg', 'modified_sequence']
        columns_to_drop = [col for col in columns_to_drop if col in df.columns]
        df_tsv = df.drop(columns=columns_to_drop)
        
        df_tsv.to_csv(args.output, sep='\t', index=False)
        print(f"Results saved to {args.output}")
        
        # Write FASTA output if requested
        if args.translationAI_fasta:
            write_fasta_output(df, args.translationAI_fasta)

    else:
        print("No PTC-introducing variants found")
    
    genome_fasta.close()


if __name__ == '__main__':
    main()