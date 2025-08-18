#!/usr/bin/env python3
import argparse
import os
import subprocess
import shutil
import pandas as pd
import numpy as np

class ToolChecker:
    def __init__(self, tools):
        self.tools = tools
        self.report = []

    def check_tools(self):
        for tool in self.tools:
            if tool == 'picard':
                try:
                    picard_path = shutil.which('picard')
                    if picard_path:
                        picard_dir = os.path.dirname(picard_path)
                        picard_jar = os.path.join(picard_dir, '..', 'share', 'picard', 'picard.jar')  # Adjust path as needed
                        if os.path.exists(picard_jar):
                            self.report.append({'tool': tool, 'status': 'Installed', 'path': picard_jar})
                            os.environ['PICARD'] = picard_jar
                        else:
                            self.report.append({'tool': tool, 'status': 'Not Installed', 'path': 'N/A'})
                except:
                    self.report.append({'tool': tool, 'status': 'Not Installed', 'path': 'N/A'})
            else:
                path = shutil.which(tool)
                if path:
                    self.report.append({'tool': tool, 'status': 'Installed', 'path': path})
                else:
                    self.report.append({'tool': tool, 'status': 'Not Installed', 'path': 'N/A'})
        report_df = pd.DataFrame(self.report)
        report_df.to_csv('tool_report.csv', index=False)
        print("Tool Check Report:")
        print(report_df)
        return report_df

class GFFtoBED:
    def __init__(self, gff_file, output_bed):
        self.gff_file = gff_file
        self.output_bed = output_bed

    def convert(self):
        cmd = f"gff2bed < {self.gff_file} | awk '{{if($8==\"gene\") print $1,$2,$3}}' OFS='\\t' > {self.output_bed}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"GFF converted to BED: {self.output_bed}")

class ReadMapping:
    def __init__(self, ref_genome, fread, rread, output_dir, cpus=28):
        self.ref_genome = ref_genome
        self.fread = fread
        self.rread = rread
        self.output_dir = output_dir
        self.cpus = cpus
        self.out = f"{os.path.basename(self.ref_genome)}_{os.path.basename(self.fread).split('_')[0]}"

    def index_genome(self):
        subprocess.run(f"samtools faidx {self.ref_genome}", shell=True, check=True)
        subprocess.run(f"bwa-mem2 index {self.ref_genome}", shell=True, check=True)
        print(f"Indexed genome: {self.ref_genome}")

    def map_reads(self):
        cmd = f"bwa-mem2 mem -t {self.cpus} {self.ref_genome} {self.fread} {self.rread} | samtools view - -Sb -@{self.cpus} | samtools view -b -@{self.cpus} -F 4 | samtools sort - -@{self.cpus} -o {self.output_dir}/{self.out}.allMapped.sorted.bam"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Mapped reads: {self.output_dir}/{self.out}.allMapped.sorted.bam")

    def remove_duplicates(self):
        bam_file = f"{self.output_dir}/{self.out}.allMapped.sorted.bam"
        if 'PICARD' in os.environ:
            picard = os.environ['PICARD']
            cmd = f"java -Xmx100G -jar {picard} MarkDuplicates INPUT={bam_file} O={self.output_dir}/{self.out}.allMapped.sorted.markdup.bam M={self.output_dir}/{self.out}.allMapped.sorted.markdup.bam.metrics.txt REMOVE_DUPLICATES=True"
            subprocess.run(cmd, shell=True, check=True)
            print("Duplicates removed with Picard")
        else:
            print("Picard not installed, using samtools for removing duplicates")
            coll_cmd = f"samtools collate -o {self.output_dir}/{self.out}.allMapped.sorted.collate.bam {bam_file} -@{self.cpus}"
            subprocess.run(coll_cmd, shell=True, check=True)
            os.remove(bam_file)
            fixmate_cmd = f"samtools fixmate -m {self.output_dir}/{self.out}.allMapped.sorted.collate.bam {self.output_dir}/{self.out}.allMapped.sorted.fixmate.bam -@{self.cpus}"
            subprocess.run(fixmate_cmd, shell=True, check=True)
            os.remove(f"{self.output_dir}/{self.out}.allMapped.sorted.collate.bam")
            sort_cmd = f"samtools sort -o {self.output_dir}/{self.out}.allMapped.sorted.02.bam {self.output_dir}/{self.out}.allMapped.sorted.fixmate.bam -@{self.cpus}"
            subprocess.run(sort_cmd, shell=True, check=True)
            os.remove(f"{self.output_dir}/{self.out}.allMapped.sorted.fixmate.bam")
            markdup_cmd = f"samtools markdup {self.output_dir}/{self.out}.allMapped.sorted.02.bam {self.output_dir}/{self.out}.allMapped.sorted.markdup.bam -@{self.cpus}"
            subprocess.run(markdup_cmd, shell=True, check=True)
            os.remove(f"{self.output_dir}/{self.out}.allMapped.sorted.02.bam")
            print("Duplicates removed with samtools")

    def index_bam(self):
        bam_file = f"{self.output_dir}/{self.out}.allMapped.sorted.markdup.bam"
        subprocess.run(f"samtools index {bam_file} -@ {self.cpus}", shell=True, check=True)
        print(f"Indexed BAM: {bam_file}")

class CoverageAnalysis:
    def __init__(self, bam_file, output_dir, gff_bed, reads, ref_genome_prefix, cpus=28):
        self.bam_file = bam_file
        self.output_dir = output_dir
        self.gff_bed = gff_bed
        self.reads = reads
        self.ref_genome_prefix = ref_genome_prefix
        self.cpus = cpus

    def run_qualimap(self):
        qualimap_dir = os.path.join(self.output_dir, 'qualimap')
        cmd = f"qualimap bamqc -bam {self.bam_file} -outdir {qualimap_dir} -outfile {self.reads}_{self.ref_genome_prefix}.qualimap -sd -c -nt {self.cpus} -outformat PDF:HTML -ip --java-mem-size=200G"
        subprocess.run(cmd, shell=True, check=True)
        os.rename(os.path.join(qualimap_dir, 'genome_results.txt'), os.path.join(qualimap_dir, f"{self.reads}_genome_results.txt"))

        qualimap_by_region_dir = os.path.join(self.output_dir, 'qualimap_by_region')
        cmd2 = f"qualimap bamqc -bam {self.bam_file} -outdir {qualimap_by_region_dir} -outfile {self.reads}_{self.ref_genome_prefix}.qualimap_genes -sd -c -nt {self.cpus} -gff {self.gff_bed} -oc {self.reads}_{self.ref_genome_prefix}.qualimap_genes_cov.txt -os {self.reads}_{self.ref_genome_prefix}.qualimap_repeats -outformat PDF:HTML -ip --java-mem-size=200G"
        subprocess.run(cmd2, shell=True, check=True)
        os.rename(os.path.join(qualimap_by_region_dir, 'genome_results.txt'), os.path.join(qualimap_by_region_dir, f"{self.reads}_genome_results.txt"))
        print("Qualimap analysis completed")

    def run_bamcoverage(self):
        bw_file = os.path.join(self.output_dir, f"{os.path.basename(self.bam_file).replace('.bam', '')}.coverage.bw")
        cmd = f"bamCoverage -b {self.bam_file} --numberOfProcessors {self.cpus} -o {bw_file}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Coverage bigwig generated: {bw_file}")

def main():
    parser = argparse.ArgumentParser(description="Genomic Analysis Pipeline")
    parser.add_argument("--query_genome_prefix", required=True, help="Query genome prefix")
    parser.add_argument("--ref_genome_prefix", required=True, help="Reference genome prefix")
    parser.add_argument("--ref_genome", required=True, help="Reference genome FASTA")
    parser.add_argument("--ref_genome_gff", required=True, help="Reference genome GFF")
    parser.add_argument("--reads", required=True, help="Reads prefix")
    parser.add_argument("--workdir", default=os.getcwd(), help="Working directory")
    parser.add_argument("--cpus", default=28, type=int, help="Number of CPUs")
    args = parser.parse_args()

    # Set paths
    workdir = args.workdir
    path_wgs = os.path.join(workdir, f"{args.query_genome_prefix}_trim_out")
    path_genomes = "/home/data/salixomics/sadik/evmab_dir/genome_dir"  # Hard-coded from script
    appsdir = "/home/data/salixomics/sadik/apps/"  # Hard-coded from script
    output_dir = os.path.join(workdir, f"panEV_{args.ref_genome_prefix}_out")
    os.makedirs(output_dir, exist_ok=True)
    fread = os.path.join(path_wgs, f"{args.reads}_1_val_1.fq.gz")
    rread = os.path.join(path_wgs, f"{args.reads}_2_val_2.fq.gz")
    out = f"{args.ref_genome}_{args.reads}"
    bed_file = f"{args.ref_genome_prefix}.gene.bed"

    # 1. Checking tools
    tools = ['bwa-mem2', 'samtools', 'bedtools', 'qualimap', 'picard', 'bamCoverage']
    checker = ToolChecker(tools)
    checker.check_tools()

    # 2. Convert GFF to BED
    gff_to_bed = GFFtoBED(args.ref_genome_gff, bed_file)
    gff_to_bed.convert()

    # 3. Read Mapping
    read_mapping = ReadMapping(args.ref_genome, fread, rread, output_dir, args.cpus)
    read_mapping.index_genome()
    read_mapping.map_reads()
    read_mapping.remove_duplicates()
    read_mapping.index_bam()

    # 4. Coverage Analysis
    bam_file = os.path.join(output_dir, f"{out}.allMapped.sorted.markdup.bam")
    coverage_analysis = CoverageAnalysis(bam_file, output_dir, bed_file, args.reads, args.ref_genome_prefix, args.cpus)
    coverage_analysis.run_qualimap()
    coverage_analysis.run_bamcoverage()

if __name__ == "__main__":
    main()