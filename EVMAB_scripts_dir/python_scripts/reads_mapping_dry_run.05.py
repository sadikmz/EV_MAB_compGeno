#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import shutil
import pybedtools
import re

class ToolChecker:
    def __init__(self, tool_paths, picard_jar):
        """Initialize with dictionary of tool paths and Picard JAR path."""
        self.tool_paths = tool_paths
        self.picard_jar = picard_jar
        self.tools = ['bwa-mem2', 'samtools', 'bedtools', 'qualimap', 'bamCoverage']
        self.report = []

    def check_tools_dry_run(self):
        """Print commands to check tool availability and generate stdout report for missing tools."""
        commands = []
        missing_tools = []
        # Check Picard JAR
        if self.picard_jar:
            if os.path.exists(self.picard_jar):
                commands.append({'tool': 'picard_jar', 'command': f"echo 'Picard JAR found: {self.picard_jar}'", 'type': 'check'})
                os.environ['PICARD'] = self.picard_jar
            else:
                commands.append({'tool': 'picard_jar', 'command': f"echo 'Picard JAR not found: {self.picard_jar}'", 'type': 'check'})
                missing_tools.append('picard_jar')
        else:
            commands.append({'tool': 'picard_jar', 'command': f"echo 'Picard JAR not provided'", 'type': 'check'})
            missing_tools.append('picard_jar')
        # Check other tools
        for tool in self.tools:
            path = self.tool_paths.get(tool, shutil.which(tool))
            check_cmd = f"which {tool}"
            if self.tool_paths.get(tool):
                check_cmd = f"PATH={self.tool_paths[tool]}:$PATH which {tool}"
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
    def __init__(self, ref_genome, fread, rread, output_dir, genome_prefix, cpus=28, bwa_mem2_path=None, samtools_path=None):
        """Initialize with reference genome, FASTQ files, output directory, prefix, CPUs, and tool paths."""
        self.ref_genome = ref_genome
        self.fread = fread
        self.rread = rread
        self.output_dir = output_dir
        self.genome_prefix = genome_prefix
        self.cpus = cpus
        self.bwa_mem2_path = bwa_mem2_path or 'bwa-mem2'
        self.samtools_path = samtools_path or 'samtools'
        self.out = f"{self.genome_prefix}_{os.path.basename(self.fread).split('_')[0]}"

    def index_genome_dry_run(self):
        """Print commands for indexing genome."""
        commands = [
            {'command': f"{self.samtools_path} faidx {self.ref_genome}", 'type': 'index_faidx'},
            {'command': f"{self.bwa_mem2_path} index {self.ref_genome}", 'type': 'index_bwa'}
        ]
        return commands

    def map_reads_dry_run(self):
        """Print command for mapping paired-end reads."""
        cmd = f"{self.bwa_mem2_path} mem -t {self.cpus} {self.ref_genome} {self.fread} {self.rread} | {self.samtools_path} view - -Sb -@{self.cpus} | {self.samtools_path} view -b -@{self.cpus} -F 4 | {self.samtools_path} sort - -@{self.cpus} -o {self.output_dir}/{self.out}.allMapped.sorted.bam"
        return [{'command': cmd, 'type': 'map_reads'}]

    def remove_duplicates_dry_run(self, memory_size):
        """Print commands for removing duplicates with Picard or samtools."""
        bam_file = f"{self.output_dir}/{self.out}.allMapped.sorted.bam"
        commands = []
        picard_cmd = f"java -Xmx{memory_size} -jar $PICARD MarkDuplicates INPUT={bam_file} O={self.output_dir}/{self.out}.allMapped.sorted.markdup.bam M={self.output_dir}/{self.out}.allMapped.sorted.markdup.bam.metrics.txt REMOVE_DUPLICATES=True"
        commands.append({'command': picard_cmd, 'type': 'picard_markdup'})

        # Samtools fallback commands
        samtools_cmds = [
            f"{self.samtools_path} collate -o {self.output_dir}/{self.out}.allMapped.sorted.collate.bam {bam_file} -@{self.cpus}",
            f"rm {bam_file}",
            f"{self.samtools_path} fixmate -m {self.output_dir}/{self.out}.allMapped.sorted.collate.bam {self.output_dir}/{self.out}.allMapped.sorted.fixmate.bam -@{self.cpus}",
            f"rm {self.output_dir}/{self.out}.allMapped.sorted.collate.bam",
            f"{self.samtools_path} sort -o {self.output_dir}/{self.out}.allMapped.sorted.02.bam {self.output_dir}/{self.out}.allMapped.sorted.fixmate.bam -@{self.cpus}",
            f"rm {self.output_dir}/{self.out}.allMapped.sorted.fixmate.bam",
            f"{self.samtools_path} markdup {self.output_dir}/{self.out}.allMapped.sorted.02.bam {self.output_dir}/{self.out}.allMapped.sorted.markdup.bam -@{self.cpus}",
            f"rm {self.output_dir}/{self.out}.allMapped.sorted.02.bam"
        ]
        for cmd in samtools_cmds:
            commands.append({'command': cmd, 'type': 'samtools_markdup'})
        return commands

    def index_bam_dry_run(self):
        """Print command for indexing BAM file."""
        bam_file = f"{self.output_dir}/{self.out}.allMapped.sorted.markdup.bam"
        cmd = f"{self.samtools_path} index {bam_file} -@ {self.cpus}"
        return [{'command': cmd, 'type': 'index_bam'}]

class CoverageAnalysis:
    def __init__(self, bam_file, output_dir, gff_bed, reads, genome_prefix, cpus=28, qualimap_path=None, bamcoverage_path=None):
        """Initialize with BAM file, output directory, BED file, reads prefix, genome prefix, CPUs, and tool paths."""
        self.bam_file = bam_file
        self.output_dir = output_dir
        self.gff_bed = gff_bed
        self.reads = reads
        self.genome_prefix = genome_prefix
        self.cpus = cpus
        self.qualimap_path = qualimap_path or 'qualimap'
        self.bamcoverage_path = bamcoverage_path or 'bamCoverage'

    def run_qualimap_dry_run(self, memory_size):
        """Print commands for Qualimap analysis using user-provided memory size."""
        qualimap_dir = os.path.join(self.output_dir, 'qualimap')
        qualimap_by_region_dir = os.path.join(self.output_dir, 'qualimap_by_region')
        commands = [
            {
                'command': f"{self.qualimap_path} bamqc -bam {self.bam_file} -outdir {qualimap_dir} -outfile {self.reads}_{self.genome_prefix}.qualimap -sd -c -nt {self.cpus} -outformat PDF:HTML -ip --java-mem-size={memory_size}",
                'type': 'qualimap_genome'
            },
            {
                'command': f"mv {qualimap_dir}/genome_results.txt {qualimap_dir}/{self.reads}_genome_results.txt",
                'type': 'qualimap_rename'
            },
            {
                'command': f"{self.qualimap_path} bamqc -bam {self.bam_file} -outdir {qualimap_by_region_dir} -outfile {self.reads}_{self.genome_prefix}.qualimap_genes -sd -c -nt {self.cpus} -gff {self.gff_bed} -oc {self.reads}_{self.genome_prefix}.qualimap_genes_cov.txt -os {self.reads}_{self.genome_prefix}.qualimap_repeats -outformat PDF:HTML -ip --java-mem-size={memory_size}",
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
        bw_file = os.path.join(self.output_dir, f"{os.path.basename(self.bam_file).replace('.bam', '')}.