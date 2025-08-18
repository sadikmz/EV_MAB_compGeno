#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import pybedtools
import re

class CDSPresenceAbsence:
    def __init__(self, bam_file, output_dir, gene_bed, reads, genome_prefix, ref_cds_fasta, query_genome, genome_file, cpus=28, bedtools_path=None):
        """Initialize with BAM file, output directory, gene BED file, reads prefix, genome prefix, reference CDS FASTA, query genome, genome file, CPUs, and bedtools path."""
        self.bam_file = bam_file
        self.output_dir = output_dir
        self.gene_bed = gene_bed
        self.reads = reads
        self.genome_prefix = genome_prefix
        self.ref_cds_fasta = ref_cds_fasta
        self.query_genome = query_genome
        self.genome_file = genome_file
        self.cpus = cpus
        self.bedtools_path = shutil.which('bedtools') or bedtools_path or 'bedtools'

    def generate_gene_coverage_dry_run(self, memory_size):
        """Print command to generate coverage for all mapped reads in genic regions."""
        cov_file = os.path.join(self.output_dir, f"{self.reads}.allMapped.reads_gene.cov.bed")
        cmd = f"{self.bedtools_path} bamtobed -i {self.bam_file} | {self.bedtools_path} coverage -a {self.gene_bed} -iobuf {memory_size} -b - > {cov_file}"
        return [{'command': cmd, 'type': 'gene_coverage'}]

    def generate_fraction_coverage_dry_run(self):
        """Print commands to generate fraction coverage for genic regions < 0.25 and extended PAV output."""
        cov_file = os.path.join(self.output_dir, f"{self.reads}.allMapped.reads_gene.cov.bed")
        frac_file = os.path.join(self.output_dir, f"{self.reads}.allMapped.gene.l25p.reads_cov.bed")
        commands = [
            {
                'command': f"{self.bedtools_path} intersect -a {cov_file} -b {self.gene_bed} -wo | awk '{{if($5<0.24) print $1,$2,$3,$9,$5}}' OFS='\\t' > {frac_file}",
                'type': 'fraction_coverage'
            }
        ]
        return commands

    def extract_cds_sequences_dry_run(self):
        """Print command to extract CDS sequences using gene IDs from fraction coverage."""
        frac_file = os.path.join(self.output_dir, f"{self.reads}.allMapped.gene.l25p.reads_cov.bed")
        cds_fasta = os.path.join(self.output_dir, f"{self.reads}.allMapped.gene.l25p.cov.cds.fasta")
        cmd = f"pybedtools extract CDS sequences from {self.query_genome} using gene IDs from {frac_file} to {cds_fasta}"
        return [{'command': cmd, 'type': 'extract_cds'}]

class NucmerAlignment:
    def __init__(self, output_dir, reads, genome_prefix, ref_cds_fasta, query_genome, reads_mapping_cds_dir, cpus=28, nucmer_path=None, bedtools_path=None):
        """Initialize with output directory, reads prefix, genome prefix, reference CDS FASTA, query genome, reads mapping CDS directory, CPUs, and tool paths."""
        self.output_dir = output_dir
        self.reads = reads
        self.genome_prefix = genome_prefix
        self.ref_cds_fasta = ref_cds_fasta
        self.query_genome = query_genome
        self.reads_mapping_cds_dir = reads_mapping_cds_dir
        self.cpus = cpus
        self.nucmer_path = shutil.which('nucmer') or nucmer_path or 'nucmer'
        self.bedtools_path = shutil.which('bedtools') or bedtools_path or 'bedtools'

    def align_nucmer_dry_run(self, ref_genome_prefix):
        """Print commands for Nucmer alignment and processing."""
        cds_fasta = os.path.join(self.output_dir, f"{self.reads}.allMapped.gene.l25p.cov.cds.fasta")
        prefix = f"{self.genome_prefix}_{ref_genome_prefix}"
        delta_file = f"{prefix}.delta"
        coords_file = f"{prefix}.coords"
        nucmer_bed = f"{ref_genome_prefix}_{self.genome_prefix}.nucmer.bed"
        nucmer_sorted_bed = f"{ref_genome_prefix}_{self.genome_prefix}.nucmer.sorted.merged.bed"
        ref_specific_bed = f"{ref_genome_prefix}_specific_CDS_vs_{self.genome_prefix}.nucmer_CDS.bed"
        query_nucmer_bed = f"{self.genome_prefix}_{ref_genome_prefix}.nucmer.bed"
        query_nucmer_sorted_bed = f"{self.genome_prefix}_{ref_genome_prefix}.nucmer.sorted.merged.bed"
        query_specific_bed = f"{self.genome_prefix}_specific_CDS_vs_{ref_genome_prefix}.nucmer_cds.bed"
        commands = [
            {
                'command': f"{self.nucmer_path} --prefix {prefix} {self.ref_cds_fasta} {cds_fasta} --batch 1 --threads {self.cpus}",
                'type': 'nucmer'
            },
            {
                'command': f"show-coords -c -d -l -r -o -T {delta_file} > {coords_file}",
                'type': 'show_coords'
            },
            {
                'command': f"cat {coords_file} | grep -v '=\|/\|NUCMER' | sed 's/|//g' | grep -v \"\\[\" | awk '{{print $14,$1,$2,$5,$7,$10,$8}}' | grep \"[a-zA-Z]\" | awk '{{ if ($2>$3) print $1,$3,$2,\"-\", $4,$5,$6,$7; else print $1,$2,$3,\"+\", $4,$5,$6,$7}}' OFS='\\t' > {nucmer_bed}",
                'type': 'coords_to_bed_ref'
            },
            {
                'command': f"{self.bedtools_path} sort -i {nucmer_bed} | {self.bedtools_path} merge > {nucmer_sorted_bed}",
                'type': 'sort_merge_ref'
            },
            {
                'command': f"{self.bedtools_path} intersect -a {self.reads_mapping_cds_dir}/EV_mazia_panma_cov.sorted.merged.bed -b {nucmer_sorted_bed} -wao | awk '{{if ($7/($3-$2+1)*100 < 24.4) print $1,$2,$3,$4,$5,$6,$7,$7/($3-$2)*100}}' OFS='\\t' > {ref_specific_bed}",
                'type': 'intersect_ref'
            },
            {
                'command': f"cat {coords_file} | grep -v '=\|/\|NUCMER' | sed 's/|//g' | grep -v \"\\[\" | awk '{{print $15,$3,$4,$6,$7,$11,$9}}' | grep \"[a-zA-Z]\" | awk '{{ if ($2>$3) print $1,$3,$2,\"-\", $4,$5,$6,$7; else print $1,$2,$3,\"+\", $4,$5,$6,$7}}' OFS='\\t' > {query_nucmer_bed}",
                'type': 'coords_to_bed_query'
            },
            {
                'command': f"{self.bedtools_path} sort -i {query_nucmer_bed} | {self.bedtools_path} merge | sed 's/:/\\t/g' | sed 's/-/\t/g' | awk '{{print $1, $2+$4-1,$2+$5-1}}' OFS='\\t' | {self.bedtools_path} sort | {self.bedtools_path} merge > {query_nucmer_sorted_bed}",
                'type': 'sort_merge_query'
            },
            {
                'command': f"{self.bedtools_path} intersect -a {self.reads_mapping_cds_dir}/MA_panev_cov.sorted.merged.bed -b {query_nucmer_sorted_bed} -wao | awk '{{if ($7/($3-$2+1)*100 < 25) print $1,$2,$3,$4,$5,$6,$7,$7/($3-$2)*100}}' OFS='\\t' > {query_specific_bed}",
                'type': 'intersect_query'
            }
        ]
        return commands

def main():
    parser = argparse.ArgumentParser(description="Dry-run: Print bash commands for CDS fraction coverage and Nucmer alignment (Sections 2.13-2.15)")
    parser.add_argument("--query_genome", required=True, help="Path to query genome FASTA")
    parser.add_argument("--ref_cds_fasta", required=True, help="Path to reference CDS FASTA")
    parser.add_argument("--genome_prefix", required=True, help="Query genome prefix (max 4 alphabetic letters, capitalized)", type=str)
    parser.add_argument("--ref_genome_prefix", required=True, help="Reference genome prefix (max 4 alphabetic letters, capitalized)", type=str)
    parser.add_argument("--memory_size", required=True, help="Memory size for bedtools (e.g., 50G)")
    parser.add_argument("--sra_csv", required=True, help="Path to SRA CSV file (e.g., EV_sra.txt.csv)")
    parser.add_argument("--gene_bed", required=True, help="Path to gene BED file")
    parser.add_argument("--genome_file", required=True, help="Path to genome file for bedtools coverage")
    parser.add_argument("--reads_mapping_cds_dir", required=True, help="Directory containing CDS mapping files (e.g., EV_mazia_panma_cov.sorted.merged.bed)")
    parser.add_argument("--workdir", default=os.getcwd(), help="Working directory")
    parser.add_argument("--cpus", default=28, type=int, help="Number of CPUs")
    parser.add_argument("--bedtools", help="Path to bedtools executable")
    parser.add_argument("--nucmer", help="Path to nucmer executable")
    parser.add_argument("--extended_PAV_out", help="Path to extended PAV output file (tab-separated)")
    args = parser.parse_args()

    # Validate genome prefixes: max 4 alphabetic letters, capitalize
    if not re.match(r'^[A-Za-z]{1,4}$', args.genome_prefix):
        raise ValueError("genome_prefix must be 1-4 alphabetic letters")
    args.genome_prefix = args.genome_prefix.upper()
    if not re.match(r'^[A-Za-z]{1,4}$', args.ref_genome_prefix):
        raise ValueError("ref_genome_prefix must be 1-4 alphabetic letters")
    args.ref_genome_prefix = args.ref_genome_prefix.upper()

    # Append 'G' to memory_size if not present
    memory_size = args.memory_size if args.memory_size.endswith('G') else f"{args.memory_size}G"

    # Set output directory
    output_dir = os.path.join(args.workdir, f"pan{args.genome_prefix}_out") if args.workdir else os.path.join(os.getcwd(), f"pan{args.genome_prefix}_out")

    # Load SRA CSV and filter to rows 2 and 3
    sra_df = pd.read_csv(args.sra_csv, dtype=str).dropna(how='all')
    sra_df = sra_df[sra_df['sra'].notna() & (sra_df['sra'] != '#REF!')]
    sra_df = sra_df.iloc[1:3]  # Rows 2 and 3: Jungleseed, Bedadeti

    # Collect commands
    commands = []

    # Process each SRA
    for _, row in sra_df.iterrows():
        sra_ids = row['sra'].split(',') if ',' in str(row['sra']) else [row['sra']]
        for sra_id in sra_ids:
            sra_id = sra_id.strip()
            bam_file = os.path.join(output_dir, f"{args.genome_prefix}_{sra_id}.allMapped.sorted.markdup.bam")
            cds_pav = CDSPresenceAbsence(bam_file, output_dir, args.gene_bed, sra_id, args.genome_prefix, args.ref_cds_fasta, args.query_genome, args.genome_file, args.cpus, args.bedtools)
            commands.extend(cds_pav.generate_gene_coverage_dry_run(memory_size))
            commands.extend(cds_pav.generate_fraction_coverage_dry_run())
            commands.extend(cds_pav.extract_cds_sequences_dry_run())
            nucmer = NucmerAlignment(output_dir, sra_id, args.genome_prefix, args.ref_cds_fasta, args.query_genome, args.reads_mapping_cds_dir, args.cpus, args.nucmer, args.bedtools)
            commands.extend(nucmer.align_nucmer_dry_run(args.ref_genome_prefix))

    # Extended PAV output (placeholder command, requires execution for actual data)
    if args.extended_PAV_out:
        commands.append({
            'command': f"Generate extended PAV table with gene_id, start, end, coverage per SRA, and cov_fraction to {args.extended_PAV_out}",
            'type': 'extended_pav'
        })

    # Save commands to CSV
    commands_df = pd.DataFrame(commands)
    commands_df.to_csv(os.path.join(output_dir, 'cds_nucmer_commands_dry_run.csv'), index=False)

    # Print commands
    print("Dry-run: Bash Commands for CDS Fraction Coverage and Nucmer Alignment")
    for _, row in commands_df.iterrows():
        print(f"Type: {row['type']}")
        print(f"Command: {row['command']}\n")

if __name__ == "__main__":
    main()