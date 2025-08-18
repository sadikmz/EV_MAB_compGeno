#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import shutil
import pybedtools
import re

class ToolChecker:
    def __init__(self, tools, tools_path=None):
        """Initialize with list of tools and optional tools path."""
        self.tools = tools
        self.tools_path = tools_path
        self.report = []

    def check_tools_dry_run(self):
        """Print commands to check tool availability and generate stdout report for missing tools."""
        commands = []
        missing_tools = []
        for tool in self.tools:
            if tool == 'picard':
                # Pythonic check for picard.jar
                picard_path = shutil.which('picard', path=self.tools_path)
                if picard_path:
                    picard_dir = os.path.dirname(picard_path)
                    picard_jar = os.path.join(picard_dir, '..', 'share', 'picard', 'picard.jar')
                    if os.path.exists(picard_jar):
                        commands.append({'tool': tool, 'command': f"echo 'Picard JAR found: {picard_jar}'", 'type': 'check'})
                        os.environ['PICARD'] = picard_jar
                    else:
                        commands.append({'tool': tool, 'command': f"echo 'Picard JAR not found'", 'type': 'check'})
                        missing_tools.append(tool)
                else:
                    commands.append({'tool': tool, 'command': f"echo 'Picard not installed'", 'type': 'check'})
                    missing_tools.append(tool)
            else:
                path = shutil.which(tool, path=self.tools_path)
                check_cmd = f"which {tool}"
                if self.tools_path:
                    check_cmd = f"PATH={self.tools_path}:$PATH which {tool}"
                commands.append({'tool': tool, 'command': check_cmd, 'type': 'check'})
                if not path:
                    missing_tools.append(tool)
        # Print stdout report for missing tools
        if missing_tools:
            print("Missing tools:", ", ".join(missing_tools))
        else:
            print("All tools are installed.")
        return commands

class GFFtoBED:
    def __init__(self, gff_file, output_bed):
        """Initialize with GFF file and output BED file."""
        self.gff_file = gff_file
        self.output_bed = output_bed

    def convert_dry_run(self):
        """Print command for converting GFF to BED using pybedtools."""
        cmd = f"pybedtools filter gene features from {self.gff_file} to {self.output_bed} with columns seqid,start,end"
        return [{'command': cmd, 'type': 'gff2bed'}]

class ReadMapping:
    def __init__(self, ref_genome, fread, rread, output_dir, genome_prefix, cpus=28):
        """Initialize with reference genome, FASTQ files, output directory, prefix, and CPUs."""
        self.ref_genome = ref_genome
        self.fread = fread
        self.rread = rread
        self.output_dir = output_dir
        self.genome_prefix = genome_prefix
        self.cpus = cpus
        self.out = f"{self.genome_prefix}_{os.path.basename(self.fread).split('_')[0]}"

    def index_genome_dry_run(self):
        """Print commands for indexing genome."""
        commands = [
            {'command': f"samtools faidx {self.ref_genome}", 'type': 'index_faidx'},
            {'command': f"bwa-mem2 index {self.ref_genome}", 'type': 'index_bwa'}
        ]
        return commands

    def map_reads_dry_run(self):
        """Print command for mapping paired-end reads."""
        cmd = f"bwa-mem2 mem -t {self.cpus} {self.ref_genome} {self.fread} {self.rread} | samtools view - -Sb -@{self.cpus} | samtools view -b -@{self.cpus} -F 4 | samtools sort - -@{self.cpus} -o {self.output_dir}/{self.out}.allMapped.sorted.bam"
        return [{'command': cmd, 'type': 'map_reads'}]

    def remove_duplicates_dry_run(self, memory_size):
        """Print commands for removing duplicates with Picard or samtools."""
        bam_file = f"{self.output_dir}/{self.out}.allMapped.sorted.bam"
        commands = []
        picard_cmd = f"java -Xmx{memory_size} -jar $PICARD MarkDuplicates INPUT={bam_file} O={self.output_dir}/{self.out}.allMapped.sorted.markdup.bam M={self.output_dir}/{self.out}.allMapped.sorted.markdup.bam.metrics.txt REMOVE_DUPLICATES=True"
        commands.append({'command': picard_cmd, 'type': 'picard_markdup'})

        # Samtools fallback commands
        samtools_cmds = [
            f"samtools collate -o {self.output_dir}/{self.out}.allMapped.sorted.collate.bam {bam_file} -@{self.cpus}",
            f"rm {bam_file}",
            f"samtools fixmate -m {self.output_dir}/{self.out}.allMapped.sorted.collate.bam {self.output_dir}/{self.out}.allMapped.sorted.fixmate.bam -@{self.cpus}",
            f"rm {self.output_dir}/{self.out}.allMapped.sorted.collate.bam",
            f"samtools sort -o {self.output_dir}/{self.out}.allMapped.sorted.02.bam {self.output_dir}/{self.out}.allMapped.sorted.fixmate.bam -@{self.cpus}",
            f"rm {self.output_dir}/{self.out}.allMapped.sorted.fixmate.bam",
            f"samtools markdup {self.output_dir}/{self.out}.allMapped.sorted.02.bam {self.output_dir}/{self.out}.allMapped.sorted.markdup.bam -@{self.cpus}",
            f"rm {self.output_dir}/{self.out}.allMapped.sorted.02.bam"
        ]
        for cmd in samtools_cmds:
            commands.append({'command': cmd, 'type': 'samtools_markdup'})
        return commands

    def index_bam_dry_run(self):
        """Print command for indexing BAM file."""
        bam_file = f"{self.output_dir}/{self.out}.allMapped.sorted.markdup.bam"
        cmd = f"samtools index {bam_file} -@ {self.cpus}"
        return [{'command': cmd, 'type': 'index_bam'}]

class CoverageAnalysis:
    def __init__(self, bam_file, output_dir, gff_bed, reads, genome_prefix, cpus=28):
        """Initialize with BAM file, output directory, BED file, reads prefix, genome prefix, and CPUs."""
        self.bam_file = bam_file
        self.output_dir = output_dir
        self.gff_bed = gff_bed
        self.reads = reads
        self.genome_prefix = genome_prefix
        self.cpus = cpus

    def run_qualimap_dry_run(self, memory_size):
        """Print commands for Qualimap analysis using user-provided memory size."""
        qualimap_dir = os.path.join(self.output_dir, 'qualimap')
        qualimap_by_region_dir = os.path.join(self.output_dir, 'qualimap_by_region')
        commands = [
            {
                'command': f"qualimap bamqc -bam {self.bam_file} -outdir {qualimap_dir} -outfile {self.reads}_{self.genome_prefix}.qualimap -sd -c -nt {self.cpus} -outformat PDF:HTML -ip --java-mem-size={memory_size}",
                'type': 'qualimap_genome'
            },
            {
                'command': f"mv {qualimap_dir}/genome_results.txt {qualimap_dir}/{self.reads}_genome_results.txt",
                'type': 'qualimap_rename'
            },
            {
                'command': f"qualimap bamqc -bam {self.bam_file} -outdir {qualimap_by_region_dir} -outfile {self.reads}_{self.genome_prefix}.qualimap_genes -sd -c -nt {self.cpus} -gff {self.gff_bed} -oc {self.reads}_{self.genome_prefix}.qualimap_genes_cov.txt -os {self.reads}_{self.genome_prefix}.qualimap_repeats -outformat PDF:HTML -ip --java-mem-size={memory_size}",
                'type': 'qualimap_genes'
            },
            {
                'command': f"mv {qualimap_by_region_dir}/genome_results.txt {qualimap_by_region_dir}/{self.reads}_genome_results.txt",
                'type': 'qualimap_rename'
            }
        ]
        return commands

    def run_bamcoverage_dry_run(self):
        """Print command for bamCoverage."""
        bw_file = os.path.join(self.output_dir, f"{os.path.basename(self.bam_file).replace('.bam', '')}.coverage.bw")
        cmd = f"bamCoverage -b {self.bam_file} --numberOfProcessors {self.cpus} -o {bw_file}"
        return [{'command': cmd, 'type': 'bamcoverage'}]

def main():
    parser = argparse.ArgumentParser(description="Dry-run: Print bash commands for genomic analysis (Sections 2.12-2.15)")
    parser.add_argument("--genome", required=True, help="Path to reference genome FASTA")
    parser.add_argument("--genome_prefix", required=True, help="Reference genome prefix (max 4 alphabetic letters, capitalized)", type=str)
    parser.add_argument("--memory_size", required=True, help="Memory size for Picard and Qualimap (e.g., 50G)")
    parser.add_argument("--tools_path", help="Optional path to tools directory")
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

    # Set output directory
    output_dir = os.path.join(args.workdir, f"pan{args.genome_prefix}_out") if args.workdir else os.path.join(os.getcwd(), f"pan{args.genome_prefix}_out")
    bed_file = f"{args.genome_prefix}.gene.bed"

    # Load SRA CSV and filter to rows 2 and 3
    sra_df = pd.read_csv(args.sra_csv, dtype=str).dropna(how='all')
    sra_df = sra_df[sra_df['sra'].notna() & (sra_df['sra'] != '#REF!')]
    sra_df = sra_df.iloc[1:3]  # Rows 2 and 3: Jungleseed, Bedadeti

    # Collect all commands
    commands = []

    # 1. Tool checking
    tools = ['bwa-mem2', 'samtools', 'bedtools', 'qualimap', 'picard', 'bamCoverage']
    checker = ToolChecker(tools, args.tools_path)
    commands.extend(checker.check_tools_dry_run())

    # 2. GFF to BED
    gff_to_bed = GFFtoBED(args.gff, bed_file)
    commands.extend(gff_to_bed.convert_dry_run())

    # 3. Read mapping for each SRA
    for _, row in sra_df.iterrows():
        sra_ids = row['sra'].split(',') if ',' in str(row['sra']) else [row['sra']]
        for sra_id in sra_ids:
            sra_id = sra_id.strip()
            fread = os.path.join(args.workdir, f"{args.genome_prefix}_trim_out", f"{sra_id}_1_val_1.fq.gz")
            rread = os.path.join(args.workdir, f"{args.genome_prefix}_trim_out", f"{sra_id}_2_val_2.fq.gz")
            read_mapping = ReadMapping(args.genome, fread, rread, output_dir, args.genome_prefix, args.cpus)
            commands.extend(read_mapping.index_genome_dry_run())
            commands.extend(read_mapping.map_reads_dry_run())
            commands.extend(read_mapping.remove_duplicates_dry_run(memory_size))
            commands.extend(read_mapping.index_bam_dry_run())

            # 4. Coverage analysis
            bam_file = os.path.join(output_dir, f"{args.genome_prefix}_{sra_id}.allMapped.sorted.markdup.bam")
            coverage_analysis = CoverageAnalysis(bam_file, output_dir, bed_file, sra_id, args.genome_prefix, args.cpus)
            commands.extend(coverage_analysis.run_qualimap_dry_run(memory_size))
            commands.extend(coverage_analysis.run_bamcoverage_dry_run())

    # Save commands to CSV
    commands_df = pd.DataFrame(commands)
    commands_df.to_csv(os.path.join(output_dir, 'tool_and_analysis_commands_dry_run.csv'), index=False)

    # Print commands
    print("Dry-run: Bash Commands for Genomic Analysis")
    for _, row in commands_df.iterrows():
        print(f"Type: {row['type']}")
        print(f"Command: {row['command']}\n")

if __name__ == "__main__":
    main()