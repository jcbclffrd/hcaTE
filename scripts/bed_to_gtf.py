#!/usr/bin/env python3
"""
Convert RepeatMasker BED file to GTF format for STAR indexing
This allows TE annotations to be included alongside gene annotations
"""

import sys

def bed_to_gtf(bed_file, gtf_file):
    """Convert BED to GTF format"""
    with open(bed_file, 'r') as infile, open(gtf_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#') or line.strip() == '':
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue
            
            chrom = fields[0]
            start = int(fields[1]) + 1  # BED is 0-based, GTF is 1-based
            end = fields[2]
            name = fields[3]
            score = fields[4] if len(fields) > 4 else '.'
            strand = fields[5] if len(fields) > 5 else '.'
            
            # Parse the name field (format: repName_repClass_repFamily)
            name_parts = name.split('_')
            rep_name = name_parts[0] if len(name_parts) > 0 else name
            rep_class = name_parts[1] if len(name_parts) > 1 else 'Unknown'
            rep_family = name_parts[2] if len(name_parts) > 2 else 'Unknown'
            
            # Create GTF attributes
            attributes = f'gene_id "{rep_name}"; transcript_id "{rep_name}"; gene_name "{rep_name}"; gene_type "TE"; repClass "{rep_class}"; repFamily "{rep_family}";'
            
            # Write GTF line (format: chr, source, feature, start, end, score, strand, frame, attributes)
            gtf_line = f'{chrom}\tRepeatMasker\texon\t{start}\t{end}\t{score}\t{strand}\t.\t{attributes}\n'
            outfile.write(gtf_line)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: bed_to_gtf.py input.bed output.gtf")
        sys.exit(1)
    
    bed_to_gtf(sys.argv[1], sys.argv[2])
    print(f"Converted {sys.argv[1]} to {sys.argv[2]}")
