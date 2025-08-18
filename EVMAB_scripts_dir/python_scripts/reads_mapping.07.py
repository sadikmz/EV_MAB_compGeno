#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import shutil
import pybedtools
import re
import subprocess

class ToolChecker:
    def __init__(self, tool_paths, picard_jar):
        """Initialize with dictionary of tool paths and Picard JAR path."""
        self.tool_paths = tool_paths
        self.picard_jar = picard_jar
        self.tools = ['bwa-mem2', 'samtools', 'bedtools', 'qualimap', 'bamCoverage']
        self.report = []

    def check_tools(self):
        """Check tool availability, prioritizing PATH, and generate stdout report."""
        missing_tools = []
        # Check Picard JAR
        if self.picard_jar and os.path.exists(self.picard_jar):
            self.report.append({'tool': 'picard_jar', 'status': f'Found: {self.picard_jar}'})
        else:
            picard_path = shutil.which('picard')
            if picard_path:
                self.report.append({'tool': 'picard_jar', 'status': f'Found picard in PATH: {picard_path}'})
            else:
                self.report.append({'tool': 'picard_jar', 'status': 'Not found'})
                missing_tools.append('picard_jar')
        # Check other tools
        for tool in self.tools:
            path = shutil.which(tool) or self.tool_paths.get(tool)
            if path:
                self.report.append({'tool': tool, 'status': f'Found: {path}'})
            else:
                self.report.append({'tool': tool, 'status': 'Not found'})
                missing_tools.append(tool)
        # Print stdout report for missing tools
        if missing_tools:
            print("Missing tools:", ", ".join(missing_tools))
        else:
            print("All tools are installed.")
        return pd.DataFrame(self.report)

class GFFtoBED:
    def __init__(self, gff_file, gene_bed):
        """Initialize with GFF file and output gene BED file."""
        self.gff_file = gff_file
        self.gene_bed = gene_bed

    def convert(self):
        """Convert GFF to BED for gene features using pybedtools."""
        try:
            bed = pybedtools.BedTool(self.gff_file).filter(lambda x: x[2] == 'gene').to_dataframe(names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])[['seqid', 'start', 'end']]
            bed.to_csv(self.gene_bed, sep='\t', header=False, index=False)
            return {'status': f'Converted GFF to BED: {self.gene_bed}'}
        except Exception as e:
            return {'status': f'Failed to convert GFF to BED: {e}'}

class ReadMapping:
    def __init__(self, ref_genome, fread, rread, output_dir, genome_prefix, cpus=28, bwa_mem2_path=None, samtools_path=None):
        """Initialize with reference genome, FASTQ files, output directory, prefix, CPUs, and tool paths."""
        self.ref_genome = ref_genome
        self.fread = fread
        self.rread = rread
        self.output_dir = output_dir
        self.genome_prefix = genome_prefix
        self.cpus = cpus
        self.bwa_mem2_path = shutil.which('bwa-mem2') or bwa_mem2_path or 'bwa-mem2'
        self.samtools_path = shutil.which('samtools') or samtools_path or 'samtools'
        self.out = f"{self.genome_prefix}_{os.path.basename(self.fread).split('_')[0]}"

    def index_genome(self):
        """Execute commands for indexing genome."""
        try:
            subprocess.run(f"{self.samtools_path} faidx {self.ref_genome}", shell=True, check=True)
            subprocess.run(f"{self.bwa_mem2_path} index {self.ref_genome}", shell=True, check=True)
            return {'status': f'Indexed genome: {self.ref_genome}'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to index genome: {e}'}

    def map_reads(self):
        """Execute command for mapping paired-end reads."""
        try:
            cmd = f"{self.bwa_mem2_path} mem -t {self.cpus} {self.ref_genome} {self.fread} {self.rread} | {self.samtools_path} view - -Sb -@{self.cpus} | {self.samtools_path} view -b -@{self.cpus} -F 4 | {self.samtools_path} sort - -@{self.cpus} -o {self.output_dir}/{self.out}.allMapped.sorted.bam"
            subprocess.run(cmd, shell=True, check=True)
            return {'status': f'Mapped reads: {self.output_dir}/{self.out}.allMapped.sorted.bam'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to map reads: {e}'}

    def remove_duplicates(self, memory_size, picard_jar):
        """Execute commands for removing duplicates with Picard or samtools."""
        bam_file = f"{self.output_dir}/{self.out}.allMapped.sorted.bam"
        try:
            if picard_jar and os.path.exists(picard_jar):
                cmd = f"java -Xmx{memory_size} -jar {picard_jar} MarkDuplicates INPUT={bam_file} O={self.output_dir}/{self.out}.allMapped.sorted.markdup.bam M={self.output_dir}/{self.out}.allMapped.sorted.markdup.bam.metrics.txt REMOVE_DUPLICATES=True"
                subprocess.run(cmd, shell=True, check=True)
                return {'status': 'Duplicates removed with Picard'}
            else:
                print("Picard JAR not provided or invalid, using samtools for duplicates removal")
                subprocess.run(f"{self.samtools_path} collate -o {self.output_dir}/{self.out}.allMapped.sorted.collate.bam {bam_file} -@{self.cpus}", shell=True, check=True)
                os.remove(bam_file)
                subprocess.run(f"{self.samtools_path} fixmate -m {self.output_dir}/{self.out}.allMapped.sorted.collate.bam {self.output_dir}/{self.out}.allMapped.sorted.fixmate.bam -@{self.cpus}", shell=True, check=True)
                os.remove(f"{self.output_dir}/{self.out}.allMapped.sorted.collate.bam")
                subprocess.run(f"{self.samtools_path} sort -o {self.output_dir}/{self.out}.allMapped.sorted.02.bam {self.output_dir}/{self.out}.allMapped.sorted.fixmate.bam -@{self.cpus}", shell=True, check=True)
                os.remove(f"{self.output_dir}/{self.out}.allMapped.sorted.fixmate.bam")
                subprocess.run(f"{self.samtools_path} markdup {self.output_dir}/{self.out}.allMapped.sorted.02.bam {self.output_dir}/{self.out}.allMapped.sorted.markdup.bam -@{self.cpus}", shell=True, check=True)
                os.remove(f"{self.output_dir}/{self.out}.allMapped.sorted.02.bam")
                return {'status': 'Duplicates removed with samtools'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to remove duplicates: {e}'}

    def index_bam(self):
        """Execute command for indexing BAM file."""
        bam_file = f"{self.output_dir}/{self.out}.allMapped.sorted.markdup.bam"
        try:
            subprocess.run(f"{self.samtools_path} index {bam_file} -@ {self.cpus}", shell=True, check=True)
            return {'status': f'Indexed BAM: {bam_file}'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to index BAM: {e}'}

class CoverageAnalysis:
    def __init__(self, bam_file, output_dir, gene_bed, reads, genome_prefix, cpus=28, qualimap_path=None, bamcoverage_path=None):
        """Initialize with BAM file, output directory, gene BED file, reads prefix, genome prefix, CPUs, and tool paths."""
        self.bam_file = bam_file
        self.output_dir = output_dir
        self.gene_bed = gene_bed
        self.reads = reads
        self.genome_prefix = genome_prefix
        self.cpus = cpus
        self.qualimap_path = shutil.which('qualimap') or qualimap_path or 'qualimap'
        self.bamcoverage_path = shutil.which('bamCoverage') or bamcoverage_path or 'bamCoverage'

    def run_qualimap(self, memory_size):
        """Execute commands for Qualimap analysis using user-provided memory size."""
        qualimap_dir = os.path.join(self.output_dir, 'qualimap')
        qualimap_by_region_dir = os.path.join(self.output_dir, 'qualimap_by_region')
        try:
            cmd = f"{self.qualimap_path} bamqc -bam {self.bam_file} -outdir {qualimap_dir} -outfile {self.reads}_{self.genome_prefix}.qualimap -sd -c -nt {self.cpus} -outformat PDF:HTML -ip --java-mem-size={memory_size}"
            subprocess.run(cmd, shell=True, check=True)
            os.rename(os.path.join(qualimap_dir, 'genome_results.txt'), os.path.join(qualimap_dir, f"{self.reads}_genome_results.txt"))
            cmd2 = f"{self.qualimap_path} bamqc -bam {self.bam_file} -outdir {qualimap_by_region_dir} -outfile {self.reads}_{self.genome_prefix}.qualimap_genes -sd -c -nt {self.cpus} -gff {self.gene_bed} -oc {self.reads}_{self.genome_prefix}.qualimap_genes_cov.txt -os {self.reads}_{self.genome_prefix}.qualimap_repeats -outformat PDF:HTML -ip --java-mem-size={memory_size}"
            subprocess.run(cmd2, shell=True, check=True)
            os.rename(os.path.join(qualimap_by_region_dir, 'genome_results.txt'), os.path.join(qualimap_by_region_dir, f"{self.reads}_genome_results.txt"))
            return {'status': 'Qualimap analysis completed'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed Qualimap analysis: {e}'}

    def run_bamcoverage(self):
        """Execute command for bamCoverage."""
        bw_file = os.path.join(self.output_dir, f"{os.path.basename(self.bam_file).replace('.bam', '')}.coverage.bw")
        try:
            cmd = f"{self.bamcoverage_path} -b {self.bam_file} --numberOfProcessors {self.cpus} -o {bw_file}"
            subprocess.run(cmd, shell=True, check=True)
            return {'status': f'Coverage bigwig generated: {bw_file}'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed bamCoverage: {e}'}

def main():
    parser = argparse.ArgumentParser(description="Execute genomic analysis (Sections 2.12-2.15)")
    parser.add_argument("--genome", required=True, help="Path to reference genome FASTA")
    parser.add_argument("--genome_prefix", required=True, help="Reference genome prefix (max 4 alphabetic letters, capitalized)", type=str)
    parser.add_argument("--memory_size", required=True, help="Memory size for Picard and Qualimap (e.g., 50G)")
    parser.add_argument("--bwa_mem2", help="Path to bwa-mem2 executable")
    parser.add_argument("--samtools", help="Path to samtools executable")
    parser.add_argument("--bedtools", help="Path to bedtools executable")
    parser.add_argument("--qualimap", help="Path to qualimap executable")
    parser.add_argument("--picard_jar", help="Path to picard.jar file")
    parser.add_argument("--bamcoverage", help="Path to bamCoverage executable")
    parser.add_argument("--sra_csv", required=True, help="Path to SRA CSV file (e.g., EV_sra.txt.csv)")
    parser.add_argument("--gff", required=True, help="Path to reference genome GFF")
    parser.add_argument("--workdir", default=os.getcwd(), help="Working directory")
    parser.add_argument("--cpus", default=28, type=int, help="Number of CPUs")
    args = parser.parse_args()

    # Validate genome_prefix: max 4 alphabetic letters, capitalize
    if not re.match(r'^[A-Za-z]{1,4}$', args.genome_prefix):
        raise ValueError("genome_prefix must be 1-4 alphabetic letters")
    args.genome_prefix = args.genome_prefix.upper()

    # Append 'G' to memory_size if not present
    memory_size = args.memory_size if args.memory_size.endswith('G') else f"{args.memory_size}G"

    # Set output directory and create it
    output_dir = os.path.join(args.workdir, f"pan{args.genome_prefix}_out") if args.workdir else os.path.join(os.getcwd(), f"pan{args.genome_prefix}_out")
    os.makedirs(output_dir, exist_ok=True)
    gene_bed = f"{args.genome_prefix}.gene.bed"

    # Set tool paths dictionary
    tool_paths = {
        'bwa-mem2': args.bwa_mem2,
        'samtools': args.samtools,
        'bedtools': args.bedtools,
        'qualimap': args.qualimap,
        'bamCoverage': args.bamcoverage
    }

    # Load SRA CSV and filter to rows 2 and 3
    sra_df = pd.read_csv(args.sra_csv, dtype=str).dropna(how='all')
    sra_df = sra_df[sra_df['sra'].notna() & (sra_df['sra'] != '#REF!')]
    sra_df = sra_df.iloc[1:3]  # Rows 2 and 3: Jungleseed, Bedadeti

    # Collect status reports
    status = []

    # 1. Tool checking
    checker = ToolChecker(tool_paths, args.picard_jar)
    status.append({'step': 'tool_check', 'status': checker.check_tools().to_dict('records')})

    # 2. GFF to BED
    gff_to_bed = GFFtoBED(args.gff, gene_bed)
    status.append({'step': 'gff2bed', 'status': gff_to_bed.convert()})

    # 3. Read mapping for each SRA
    for _, row in sra_df.iterrows():
        sra_ids = row['sra'].split(',') if ',' in str(row['sra']) else [row['sra']]
        for sra_id in sra_ids:
            sra_id = sra_id.strip()
            fread = os.path.join(args.workdir, f"{args.genome_prefix}_trim_out", f"{sra_id}_1_val_1.fq.gz")
            rread = os.path.join(args.workdir, f"{args.genome_prefix}_trim_out", f"{sra_id}_2_val_2.fq.gz")
            read_mapping = ReadMapping(args.genome, fread, rread, output_dir, args.genome_prefix, args.cpus, args.bwa_mem2, args.samtools)
            status.append({'step': f'index_genome_{sra_id}', 'status': read_mapping.index_genome()})
            status.append({'step': f'map_reads_{sra_id}', 'status': read_mapping.map_reads()})
            status.append({'step': f'remove_duplicates_{sra_id}', 'status': read_mapping.remove_duplicates(memory_size, args.picard_jar)})
            status.append({'step': f'index_bam_{sra_id}', 'status': read_mapping.index_bam()})

            # 4. Coverage analysis
            bam_file = os.path.join(output_dir, f"{args.genome_prefix}_{sra_id}.allMapped.sorted.markdup.bam")
            coverage_analysis = CoverageAnalysis(bam_file, output_dir, gene_bed, sra_id, args.genome_prefix, args.cpus, args.qualimap, args.bamcoverage)
            status.append({'step': f'qualimap_{sra_id}', 'status': coverage_analysis.run_qualimap(memory_size)})
            status.append({'step': f'bamcoverage_{sra_id}', 'status': coverage_analysis.run_bamcoverage()})

    # Save status to CSV
    status_df = pd.DataFrame(status)
    status_df.to_csv(os.path.join(output_dir, 'tool_and_analysis_status.csv'), index=False)

    # Print status
    print("Genomic Analysis Status")
    print(status_df)

if __name__ == "__main__":
    main()