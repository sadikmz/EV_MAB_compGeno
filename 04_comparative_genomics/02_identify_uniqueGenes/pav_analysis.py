#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
import tempfile

import pandas as pd
import pybedtools

# Coverage threshold for PAV detection (fractional; genes below this are called absent)
LOW_COVERAGE_THRESHOLD = 0.245

def _sh(cmd):
    subprocess.run(cmd, shell=True, check=True)

def parse_attributes(attr_str):
    """Parse GFF attributes into a dictionary."""
    return {k.strip(): v.strip() for a in attr_str.split(";") if (s := a.strip()) and '=' in s for k, v in [s.split("=", 1)]}

def get_gene(gff_path):
    """Return set of genic regions (chrom, start, end, strand, gene_id) from a GFF file."""
    try:
        return {'status': 'success', 'cds_gene': {
            (f.chrom, f.start, f.end, f.strand, gid)
            for f in pybedtools.BedTool(gff_path).filter(lambda x: x[2] == "gene")
            if (gid := f.attrs.get("ID", ''))
        }}
    except Exception as e:
        return {'status': f'Failed to process {gff_path}: {e}', 'cds_gene': None}

def get_genome_file(genome, output_dir, samtools_path='samtools', genome_prefix=''):
    """Generate a tab-separated file of sequence ID and length using samtools faidx."""
    try:
        genome_file = os.path.join(output_dir, f"{genome_prefix}_genome_file.txt")
        subprocess.run([samtools_path, "faidx", genome], check=True)
        with open(f"{genome}.fai", 'r') as fai, open(genome_file, 'w') as out:
            out.writelines(f"{f[0]}\t{f[1]}\n" for line in fai if len(f := line.strip().split("\t")) >= 2)
        return {'status': f'Successfully generated {genome_file}', 'genome_file': genome_file}
    except subprocess.CalledProcessError as e:
        return {'status': f'Failed to generate genome file: {e}', 'genome_file': None}

class ToolChecker:
    def __init__(self, tool_paths, picard_jar):
        self.tool_paths = tool_paths
        self.picard_jar = picard_jar
        self.tool = ['bwa-mem2','samtools','bedtools','qualimap','bamCoverage','nucmer']
        self.report = []

    def check_tool(self):
        """Check tools availability, prioritizing PATH, generate stdout report."""
        missing_tools = []
        picard_path = shutil.which('picard')
        self.report.append({'tool': 'picard_jar', 'status': f'Found: {self.picard_jar}'} if (self.picard_jar and os.path.exists(self.picard_jar))
                           else {'tool': 'picard', 'status': f'Found picard in PATH: {picard_path}'} if picard_path
                           else {'tool': 'picard_jar', 'status': 'Failed to find picard in the PATH'})
        for tool in self.tool:
            path = shutil.which(tool) or self.tool_paths.get(tool)
            self.report.append({'tool': tool, 'status': f'Found {tool} in PATH: {path}' if path else f'Failed to find {tool} in PATH'})
            if not path: missing_tools.append(tool)
        print("Missing tools:", ",".join(missing_tools) if missing_tools else "All tools are installed")
        return pd.DataFrame(self.report)

class GFF2BED:
    def __init__(self, gff_file):
        self.gff_file = gff_file

    def convert(self):
        return get_gene(self.gff_file)

class ReadMapping:
    def __init__(self, query_genome, query_genome_prefix, fread, rread, output_dir, cpus=28, bwa_mem2_path=None, samtools_path=None):
        self.query_genome = query_genome
        self.query_genome_prefix = query_genome_prefix
        self.fread = fread
        self.rread = rread
        self.output_dir = output_dir
        self.cpus = cpus
        self.bwa_mem2_path = shutil.which(bwa_mem2_path) or bwa_mem2_path or 'bwa-mem2'
        self.samtools_path = shutil.which(samtools_path) or samtools_path or 'samtools'
        self.out = f"{self.query_genome_prefix}_{os.path.basename(self.fread).split('_')[0]}"

    def index_genome(self):
        """Index genome with samtools faidx and bwa-mem2."""
        try:
            subprocess.run([self.samtools_path, "faidx", self.query_genome], check=True)
            subprocess.run([self.bwa_mem2_path, "index", self.query_genome], check=True)
            return {'status': f'{self.query_genome} indexed'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to index genome: {e}'}

    def map_reads(self):
        """Map reads with bwa-mem2, filter, sort, and output BAM."""
        try:
            _sh(f"{self.bwa_mem2_path} mem -t {self.cpus} {self.query_genome} {self.fread} {self.rread} "
                f"| {self.samtools_path} view - -Sb -@{self.cpus} "
                f"| {self.samtools_path} view -b -@{self.cpus} -F 4 "
                f"| {self.samtools_path} sort - -@{self.cpus} -o {self.output_dir}/{self.out}.allMapped.sorted.bam")
            return {'status': f'Mapped reads: {self.output_dir}/{self.out}.allMapped.sorted.bam'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to map reads: {e}'}

    def remove_duplicates(self, memory_size, picard_jar):
        """Remove duplicate reads using Picard (if available) or SAMtools."""
        bam = f"{self.output_dir}/{self.out}.allMapped.sorted.bam"
        o, st, s, c = self.output_dir, self.samtools_path, self.out, self.cpus
        try:
            if picard_jar and os.path.exists(picard_jar):
                _sh(f"java -Xmx{memory_size} -jar {picard_jar} MarkDuplicates "
                    f"INPUT={bam} O={o}/{s}.allMapped.sorted.markdup.bam "
                    f"M={o}/{s}.allMapped.sorted.markdup.bam.metrics.txt REMOVE_DUPLICATES=True")
                return {'status': 'PCR duplicates removed with Picard Duplicates'}
            print('Picard jar is not provided or missing, using SAMtools for duplicates removal')
            _sh(f"{st} collate -o {o}/{s}.allMapped.sorted.collate.bam {bam} -@{c}")
            os.remove(bam)
            _sh(f"{st} fixmate -m {o}/{s}.allMapped.sorted.collate.bam {o}/{s}.allMapped.sorted.fixmate.bam -@{c}")
            os.remove(f'{o}/{s}.allMapped.sorted.collate.bam')
            _sh(f"{st} sort -o {o}/{s}.allMapped.sorted.02.bam {o}/{s}.allMapped.sorted.fixmate.bam -@{c}")
            os.remove(f'{o}/{s}.allMapped.sorted.fixmate.bam')
            _sh(f"{st} markdup {o}/{s}.allMapped.sorted.02.bam {o}/{s}.allMapped.sorted.markdup.bam -@{c}")
            os.remove(f"{o}/{s}.allMapped.sorted.02.bam")
            return {'status': 'Duplicates removed with SAMtools'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to remove duplicates: {e}'}

    def index_bam(self):
        """Index the markdup BAM file."""
        try:
            _sh(f"{self.samtools_path} index {self.output_dir}/{self.out}.allMapped.sorted.markdup.bam -@{self.cpus}")
            return {'status': f'Indexed bam file: {self.output_dir}/{self.out}.allMapped.sorted.markdup.bam'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to index bam file: {e}'}

class CoverageAnalysis:
    def __init__(self, bam_file, output_dir, cds_gene, sra_ids, query_genome_prefix, cpus=28, bedtools_path=None, qualimap_path=None, bamcoverage_path=None):
        self.bam_file = bam_file
        self.output_dir = output_dir
        self.sra_ids = sra_ids
        self.cds_gene = cds_gene
        self.query_genome_prefix = query_genome_prefix
        self.cpus = cpus
        self.bamcoverage_path = shutil.which(bamcoverage_path) or bamcoverage_path or 'bamcoverage'
        self.qualimap_path, self.bedtools_path = shutil.which(qualimap_path) or qualimap_path or 'qualimap', shutil.which(bedtools_path) or bedtools_path or 'bedtools'

    def _write_temp_bed(self):
        """Write cds_gene to a temporary BED file, return its path."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
            tmp.writelines(f"{c}\t{s}\t{e}\t{st}\t{g}\n" for c, s, e, st, g in self.cds_gene)
            return tmp.name

    def run_qualimap(self, memory_size):
        """Run qualimap bamqc on the whole genome and by gene regions."""
        qd, qdr, sid, qgp = os.path.join(self.output_dir, 'qualimap'), os.path.join(self.output_dir, 'qualimap_by_regions'), self.sra_ids, self.query_genome_prefix
        try:
            tmp = self._write_temp_bed()
            _sh(f"{self.qualimap_path} bamqc -bam {self.bam_file} -output {qd} "
                f"-outfile {sid}_{qgp}.qualimap -sd -c -nt {self.cpus} "
                f"-outformat PDF:HTML ip --java-mem-size={memory_size}")
            os.rename(os.path.join(qd, 'genome_results.txt'), os.path.join(qd, f"{sid}_{qgp}.genome_results.txt"))
            _sh(f"{self.qualimap_path} bamqc -bam {self.bam_file} "
                f"-output {qdr}.qualimap_genes "
                f"-outfile {sid}_{qgp}.qualimap -sd -c -nt {self.cpus} "
                f"-gff {tmp} -oc {sid}_{qgp}.qualimap_genes_cov.txt "
                f"-os {sid}_{qgp}.qualimap_repeats "
                f"-outformat PDF:HTML ip --java-mem-size={memory_size}")
            os.rename(os.path.join(qdr, 'genome_results.txt'), os.path.join(qdr, f"{sid}_{qgp}.genome_results.txt"))
            os.remove(tmp)
            return {'status': 'Qualimap analysis completed'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to run qualimap command: {e}'}

    def generate_bigwig_coverage(self):
        """Generate bigWig coverage file with bamCoverage."""
        try:
            _sh(f"{self.bamcoverage_path} -b {self.bam_file} --numberOfProcessors {self.cpus} -o {os.path.join(self.output_dir, os.path.basename(self.bam_file).replace('.bam', '') + '.coverage.bw')}")
            return {'status': 'BigWig coverage generated'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to generate bigWig coverage: {e}'}

    def generate_gene_coverage(self, memory_size):
        """Compute per-gene read coverage using bedtools coverage."""
        try:
            tmp = self._write_temp_bed()
            _sh(f"{self.bedtools_path} bamtobed -i {self.bam_file} "
                f"| {self.bedtools_path} coverage -a {tmp} -iobuf {memory_size} -b - "
                f"> {os.path.join(self.output_dir, f'{self.sra_ids}_{self.query_genome_prefix}.allMapped.reads_gene.cov.bed')}")
            return {'status': 'Gene coverage generated'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to generate gene coverage: {e}'}

class CDSPresenceAbsence:
    def __init__(self, output_dir, query_genome_prefix, query_genome, query_genome_gff,
                 sra_ids, cds_gene, coverage_threshold=LOW_COVERAGE_THRESHOLD, cpus=28, bedtools_path=None):
        self.output_dir = output_dir
        self.query_genome_prefix = query_genome_prefix
        self.query_genome = query_genome
        self.gff_file = query_genome_gff
        self.sra_ids = sra_ids
        self.cpus = cpus
        self.cds_gene = cds_gene
        self.coverage_threshold = coverage_threshold
        self.bedtools_path = shutil.which("bedtools") or bedtools_path or "bedtools"

    def extract_reads_pav_cds_sequences(self):
        """Extract CDS coords and FASTA for genes with reads coverage below the coverage threshold."""
        reads_mapping_CDS_frac_file = os.path.join(self.output_dir, f"{self.sra_ids}_{self.query_genome_prefix}.l25p_reads_cov_gene.bed")
        try:
            low_cov_gene = set()
            cov_file = os.path.join(self.output_dir, f"{self.sra_ids}_{self.query_genome_prefix}.allMapped.reads_gene.cov.bed")
            if os.path.exists(cov_file):
                cov_data = pd.read_csv(cov_file, sep="\t", header=None)
                low_cov_gene.update(cov_data[cov_data.iloc[:, -1] < self.coverage_threshold].iloc[:, 4].unique())
            if not low_cov_gene:
                return {'status': f'No genes with coverage < {self.coverage_threshold} found'}
            gff = pybedtools.BedTool(self.gff_file).filter(lambda x: x[2] == "gene" and x.attrs.get("ID", '') in low_cov_gene)
            cds_bed = gff.to_dataframe(names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])[['seqid', 'start', 'end', 'attributes']]
            cds_bed['gene_id'] = cds_bed['attributes'].apply(lambda x: parse_attributes(x).get("ID", ''))
            cds_bed = cds_bed[cds_bed['gene_id'] != ''][['seqid', 'start', 'end', 'gene_id']]
            cds_bed.to_csv(reads_mapping_CDS_frac_file, sep='\t', index=False)
            _sh(f"{self.bedtools_path} getfasta -fi {self.query_genome} -bed {reads_mapping_CDS_frac_file} -fo {os.path.join(self.output_dir, self.query_genome_prefix + '.l25p_reads_cov_gene.fasta')}")
            return {'status': 'CDS sequences extracted'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to extract CDS sequences from {e}'}

class CDSAlignment:
    def __init__(self, output_dir, sra_ids, query_genome_prefix, query_genome, ref_genome_prefix, ref_genome,
                 ref_genome_gff, reads_mapping_cds_dir, query_cov_bed=None, ref_cov_bed=None,
                 coverage_threshold=LOW_COVERAGE_THRESHOLD,
                 cpus=28, nucmer_path=None, lastz_path=None, bedtools_path=None, samtools_path=None):
        self.output_dir = output_dir
        self.sra_ids = sra_ids
        self.query_genome_prefix = query_genome_prefix
        self.query_genome = query_genome
        self.ref_genome_prefix = ref_genome_prefix
        self.ref_genome = ref_genome
        self.ref_genome_gff = ref_genome_gff
        self.reads_mapping_cds_dir = reads_mapping_cds_dir
        # Pre-computed merged coverage BED filenames (looked up in reads_mapping_cds_dir)
        self.query_cov_bed = query_cov_bed or f"{query_genome_prefix}_panref_cov.sorted.merged.bed"
        self.ref_cov_bed = ref_cov_bed or f"{ref_genome_prefix}_panquery_cov.sorted.merged.bed"
        self.coverage_threshold = coverage_threshold
        self.cpus = cpus
        self.nucmer_path, self.lastz_path = shutil.which("nucmer") or nucmer_path or "nucmer", shutil.which("lastz") or lastz_path or "lastz"
        self.bedtools_path, self.samtools_path = shutil.which("bedtools") or bedtools_path or "bedtools", shutil.which("samtools") or samtools_path or "samtools"

    def mask_repeats(self, cds_bed, genome):
        """Mask non-CDS regions of a genome FASTA using bedtools complement + maskfasta."""
        try:
            genome_file_result = get_genome_file(genome, self.output_dir, self.samtools_path)
            genome_file = genome_file_result.get('genome_file')
            if not genome_file:
                return {'status': genome_file_result['status']}
            cds_complement = os.path.join(self.output_dir, f"{os.path.basename(genome)}.cds_complement.bed")
            masked_fasta = os.path.join(self.output_dir, os.path.basename(genome) + '.repeats_masked.fasta')
            if cds_bed and genome_file:
                _sh(f"{self.bedtools_path} complement -i {cds_bed} -g {genome_file} > {cds_complement}")
                _sh(f"{self.bedtools_path} maskfasta -fi {genome} -bed {cds_complement} -fo {masked_fasta}")
                return {'status': 'Repeats masked genome generated'}
            return {'status': 'No CDS BED or genome file provided'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to generate genome mask: {e}'}

    def generate_ref_cds_fasta(self):
        """Extract reference CDS regions to BED and FASTA."""
        try:
            ref_cds_bed = os.path.join(self.output_dir, self.ref_genome_prefix + "_ref_cds.bed")
            gff_result = get_gene(self.ref_genome_gff)
            if not gff_result['cds_gene']:
                return {'status': f'No gene features found in: {self.ref_genome_prefix}'}
            with open(ref_cds_bed, 'w') as f:
                f.writelines(f'{c}\t{s}\t{e}\t{st}\t{g}\n' for c, s, e, st, g in gff_result['cds_gene'])
            _sh(f"{self.bedtools_path} getfasta -fi {self.ref_genome} -bed {ref_cds_bed} -fo {os.path.join(self.output_dir, self.ref_genome_prefix + '_ref_cds.fasta')}")
            return {'status': f'CDS sequences extracted for: {self.ref_genome_prefix}'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to generate reference CDS fasta: {e}'}

    def align_nucmer(self, ref_genome_prefix):
        prefix = os.path.join(self.output_dir, f"{self.query_genome_prefix}_{self.ref_genome_prefix}")
        p = lambda s: f"{prefix}.{s}"
        delta_file, coords_file = p("delta"), p("coords")
        ref_nucmer_bed, nucmer_sorted_bed, ref_specific_bed = p("nucmer.bed"), p("nucmer.sorted.bed"), p("ref_specific.cds.bed")
        query_nucmer_bed, query_nucmer_sorted_bed, query_specific_bed = p("query_nucmer.bed"), p("query_nucmer.sorted.bed"), p("query_specific.cds.bed")
        bt, pct = self.bedtools_path, self.coverage_threshold * 100
        qcov = os.path.join(self.reads_mapping_cds_dir, self.query_cov_bed)
        rcov = os.path.join(self.reads_mapping_cds_dir, self.ref_cov_bed)
        try:
            _sh(f"{self.nucmer_path} --prefix {prefix} {os.path.join(self.output_dir, self.ref_genome_prefix + '_ref_cds.fasta')} {os.path.join(self.output_dir, self.query_genome_prefix + '.l25p_reads_cov_gene.fasta')} --batch 1 --threads {self.cpus}")
            _sh(f"show-coords -c -d -l -r -o -T {delta_file} > {coords_file}")
            _sh(f"cat {coords_file} | grep -v '=\\|/\\|NUCMER' | sed 's/|//g' | grep -v \"\\[\" | awk '{{print $14,$1,$2,$5,$7,$10,$8}}' | grep \"[a-zA-Z]\" | awk '{{ if ($2>$3) print $1,$3,$2,\"-\", $4,$5,$6,$7; else print $1,$2,$3,\"+\", $4,$5,$6,$7}}' OFS='\\t' > {ref_nucmer_bed}")
            _sh(f"{bt} sort -i {ref_nucmer_bed} | {bt} merge > {nucmer_sorted_bed}")
            _sh(f"{bt} intersect -a {qcov} -b {nucmer_sorted_bed} -wao | awk '{{if ($3>$2 && $7/($3-$2)*100 < {pct}) print $1,$2,$3,$4,$5,$6,$7,$7/($3-$2)*100}}' OFS='\\t' > {ref_specific_bed}")
            _sh(f"cat {coords_file} | grep -v '=\\|/\\|NUCMER' | sed 's/|//g' | grep -v \"\\[\" | awk '{{print $15,$3,$4,$6,$7,$11,$9}}' | grep \"[a-zA-Z]\" | awk '{{ if ($2>$3) print $1,$3,$2,\"-\", $4,$5,$6,$7; else print $1,$2,$3,\"+\", $4,$5,$6,$7}}' OFS='\\t' > {query_nucmer_bed}")
            _sh(f"{bt} sort -i {query_nucmer_bed} | {bt} merge | sed 's/:/\\t/g' | sed 's/-/\t/g' | awk '{{print $1, $2+$4-1,$2+$5-1}}' OFS='\\t' | {bt} sort | {bt} merge > {query_nucmer_sorted_bed}")
            _sh(f"{bt} intersect -a {rcov} -b {query_nucmer_sorted_bed} -wao | awk '{{if ($3>$2 && $7/($3-$2)*100 < {pct}) print $1,$2,$3,$4,$5,$6,$7,$7/($3-$2)*100}}' OFS='\\t' > {query_specific_bed}")
            return {'status': f'Nucmer alignment completed: {ref_specific_bed}, {query_specific_bed}'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed Nucmer alignment: {e}'}

    def align_lastz(self):
        """Align low-coverage CDS sequences using LASTZ (--strand=both --gapped --format=PAF)."""
        prefix = f"{self.query_genome_prefix}_{self.ref_genome_prefix}"
        p = lambda s: os.path.join(self.output_dir, f"{prefix}.{s}")
        lastz_paf, lastz_bed = p("lastz.paf"), p("lastz.bed")
        lastz_sorted_bed, lastz_specific_bed = p("lastz.sorted.merged.bed"), p("lastz_specific.cds.bed")
        bt, pct = self.bedtools_path, self.coverage_threshold * 100
        qcov = os.path.join(self.reads_mapping_cds_dir, self.query_cov_bed)
        try:
            _sh(f"{self.lastz_path} {os.path.join(self.output_dir, self.ref_genome_prefix + '_ref_cds.fasta')}[multiple] {os.path.join(self.output_dir, self.query_genome_prefix + '.l25p_reads_cov_gene.fasta')} --strand=both --gapped --format=PAF > {lastz_paf}")
            _sh(f"awk '{{print $6,$8,$9,$1,$12,$5}}' OFS='\\t' {lastz_paf} | awk '{{if ($2>$3) print $1,$3,$2,$4,$5,$6; else print}}' OFS='\\t' > {lastz_bed}")
            _sh(f"{bt} sort -i {lastz_bed} | {bt} merge > {lastz_sorted_bed}")
            _sh(f"{bt} intersect -a {qcov} -b {lastz_sorted_bed} -wao | awk '{{if ($3>$2 && $7/($3-$2)*100 < {pct}) print $1,$2,$3,$4,$5,$6,$7,$7/($3-$2)*100}}' OFS='\\t' > {lastz_specific_bed}")
            return {'status': f'LASTZ alignment completed: {lastz_specific_bed}'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed LASTZ alignment: {e}'}

class LongReadMapping:
    """Map long reads (PacBio HiFi, PacBio CLR, ONT) using minimap2."""
    PRESET_MAP = {'hifi': 'map-hifi', 'pacbio': 'map-pb', 'ont': 'map-ont'}

    def __init__(self, query_genome, query_genome_prefix, reads, read_type, output_dir, cpus=28, minimap2_path=None, samtools_path=None):
        self.query_genome = query_genome
        self.query_genome_prefix = query_genome_prefix
        self.reads = reads
        self.read_type = read_type
        self.output_dir = output_dir
        self.cpus = cpus
        self.minimap2_path, self.samtools_path = shutil.which(minimap2_path) or minimap2_path or 'minimap2', shutil.which(samtools_path) or samtools_path or 'samtools'
        self.preset = self.PRESET_MAP.get(read_type, 'map-hifi')
        self.out = f"{query_genome_prefix}_{os.path.basename(reads).split('.')[0]}_{read_type}"

    def map_reads(self):
        """Map long reads with minimap2 using appropriate preset."""
        bam_out = os.path.join(self.output_dir, f"{self.out}.allMapped.sorted.bam")
        try:
            _sh(f"{self.minimap2_path} -ax {self.preset} -t {self.cpus} "
                f"{self.query_genome} {self.reads} "
                f"| {self.samtools_path} view -b -F 4 -@ {self.cpus} "
                f"| {self.samtools_path} sort -@ {self.cpus} -o {bam_out}")
            return {'status': f'Long reads mapped: {bam_out}'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to map long reads: {e}'}

    def index_bam(self):
        """Index the mapped BAM file."""
        bam_file = os.path.join(self.output_dir, f"{self.out}.allMapped.sorted.bam")
        try:
            subprocess.run([self.samtools_path, 'index', bam_file], check=True)
            return {'status': f'Indexed: {bam_file}'}
        except subprocess.CalledProcessError as e:
            return {'status': f'Failed to index BAM: {e}'}


def main():
    parser = argparse.ArgumentParser(description='PAV (Presence/Absence Variation) analysis pipeline for comparative genomics.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--query-genome', required=True, help='Query genome FASTA (e.g. Ensete Mazia)')
    parser.add_argument('--query-prefix', required=True, help='Short prefix for query genome (e.g. EVMZ)')
    parser.add_argument('--query-gff', required=True, help='Query genome GFF3 annotation')
    parser.add_argument('--ref-genome', required=True, help='Reference genome FASTA (e.g. Musa acuminata)')
    parser.add_argument('--ref-prefix', required=True, help='Short prefix for reference genome (e.g. MA)')
    parser.add_argument('--ref-gff', required=True, help='Reference genome GFF3 annotation')
    parser.add_argument('--sra-ids', required=True, nargs='+', help='SRA accession IDs for short reads')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--cpus', type=int, default=28, help='Number of CPU threads')
    parser.add_argument('--coverage-threshold', type=float, default=0.245, help='Fractional coverage threshold below which gene is called absent')
    parser.add_argument('--picard-jar', default=None, help='Path to Picard JAR (optional; samtools fallback used if absent)')
    parser.add_argument('--memory-size', default='50G', help='Java/tool memory allocation (e.g. 50G)')
    parser.add_argument('--long-reads', default=None, help='Long reads FASTA/FASTQ for minimap2 mapping')
    parser.add_argument('--long-read-type', choices=['hifi', 'pacbio', 'ont'], default='hifi', help='Long read technology type')
    parser.add_argument('--skip-mapping', action='store_true', help='Skip read mapping (use pre-existing BAMs)')
    parser.add_argument('--skip-lastz', action='store_true', help='Skip LASTZ step (use NUCmer only)')
    parser.add_argument('--reads-mapping-cds-dir', default=None, help='Directory with pre-computed merged coverage BED files')
    parser.add_argument('--query-cov-bed', default=None, help='Filename (in --reads-mapping-cds-dir) of query merged coverage BED')
    parser.add_argument('--ref-cov-bed', default=None, help='Filename (in --reads-mapping-cds-dir) of reference merged coverage BED')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(ToolChecker({}, args.picard_jar).check_tool().to_string(index=False))

    print(f"\n[1/6] Parsing GFF: {args.query_gff}")
    if (gff_result := get_gene(args.query_gff))['status'] != 'success':
        raise RuntimeError(f"GFF parsing failed: {gff_result['status']}")
    cds_gene = gff_result['cds_gene']
    print(f"      Found {len(cds_gene)} gene features")

    bam_file = os.path.join(args.output_dir, f"{args.query_prefix}_{args.sra_ids[0]}.allMapped.sorted.markdup.bam")

    if not args.skip_mapping:
        print(f"\n[2/6] Mapping short reads with bwa-mem2")
        for sra_id in args.sra_ids:
            mapper = ReadMapping(args.query_genome, args.query_prefix, f"{sra_id}_1.fastq.gz", f"{sra_id}_2.fastq.gz", args.output_dir, cpus=args.cpus)
            print(mapper.index_genome()['status'])
            print(mapper.map_reads()['status'])
            print(mapper.remove_duplicates(args.memory_size, args.picard_jar)['status'])
            print(mapper.index_bam()['status'])

        if args.long_reads:
            print(f"\n[2b/6] Mapping long reads with minimap2 ({args.long_read_type})")
            lr_mapper = LongReadMapping(args.query_genome, args.query_prefix, args.long_reads, args.long_read_type, args.output_dir, cpus=args.cpus)
            print(lr_mapper.map_reads()['status'])
            print(lr_mapper.index_bam()['status'])

    print(f"\n[3/6] Coverage analysis")
    cov = CoverageAnalysis(bam_file, args.output_dir, cds_gene, args.sra_ids[0], args.query_prefix, cpus=args.cpus)
    print(cov.run_qualimap(args.memory_size)['status'])
    print(cov.generate_bigwig_coverage()['status'])
    print(cov.generate_gene_coverage(args.memory_size)['status'])

    print(f"\n[4/6] Identifying PAV genes (threshold: {args.coverage_threshold})")
    pav = CDSPresenceAbsence(args.output_dir, args.query_prefix, args.query_genome, args.query_gff, args.sra_ids[0], cds_gene, coverage_threshold=args.coverage_threshold, cpus=args.cpus)
    print(pav.extract_reads_pav_cds_sequences()['status'])

    reads_mapping_cds_dir = args.reads_mapping_cds_dir or args.output_dir
    print(f"\n[5/6] CDS alignment with NUCmer")
    aligner = CDSAlignment(args.output_dir, args.sra_ids[0], args.query_prefix, args.query_genome, args.ref_prefix, args.ref_genome, args.ref_gff, reads_mapping_cds_dir, query_cov_bed=args.query_cov_bed, ref_cov_bed=args.ref_cov_bed, coverage_threshold=args.coverage_threshold, cpus=args.cpus)
    print(aligner.generate_ref_cds_fasta()['status'])
    print(aligner.align_nucmer(args.ref_prefix)['status'])

    if not args.skip_lastz:
        print(f"\n[6/6] CDS alignment with LASTZ (genes missed by NUCmer)")
        print(aligner.align_lastz()['status'])

    print("\nPAV analysis complete.")

if __name__ == '__main__':
    main()
