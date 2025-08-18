#!/usr/bin/env python3
import argparse
import os
import pandas as pd

class GenomicAnalyzer:
    def __init__(self, sra_csv, output_dir, cpus=48):
        """Initialize with SRA CSV file, output directory, and CPU count."""
        # Load CSV and filter to rows 2 and 3 (0-based index: 1 and 2)
        self.sra_df = pd.read_csv(sra_csv, dtype=str).dropna(how='all')
        self.sra_df = self.sra_df[self.sra_df['sra'].notna() & (self.sra_df['sra'] != '#REF!')]
        self.sra_df = self.sra_df.iloc[1:3]  # Select rows 2 and 3
        self.output_dir = output_dir
        self.cpus = cpus

    def section_2_11_quality_control_dry_run(self):
        """2.11: Print bash commands for SRA download and trimming for rows 2 and 3."""
        commands = []
        for _, row in self.sra_df.iterrows():
            genotype = row['genotypes']
            sra_ids = row['sra'].split(',') if ',' in str(row['sra']) else [row['sra']]
            bioproject = row['bioproject']
            sra_pattern = '|'.join([sra_id.strip() for sra_id in sra_ids])
            sra_list = row['sra']  # Full comma-separated SRA list for download row
            # Construct esearch command for download
            download_cmd = f"esearch -db sra -query {bioproject} | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | grep -E '{sra_pattern}' | xargs -n 1 -P {self.cpus} fastq-dump --split-files --gzip --skip-technical"
            commands.append({'genotype': genotype, 'sra': sra_list, 'bioproject': bioproject, 'command': download_cmd, 'type': 'download'})

            # Construct trim_galore command for each SRA
            for sra_id in sra_ids:
                sra_id = sra_id.strip()
                fastq_1 = os.path.join(self.output_dir, f"{sra_id}_1.fastq.gz")
                fastq_2 = os.path.join(self.output_dir, f"{sra_id}_2.fastq.gz")
                trim_cmd = f"trim_galore --cores {self.cpus} --quality 30 --fastqc --output_dir {os.path.join(self.output_dir, 'trim_out')} --paired {fastq_1} {fastq_2}"
                commands.append({'genotype': genotype, 'sra': sra_id, 'bioproject': bioproject, 'command': trim_cmd, 'type': 'trim'})

        # Save commands to CSV for reference
        commands_df = pd.DataFrame(commands)
        commands_df.to_csv(os.path.join(self.output_dir, 'sra_commands_dry_run.csv'), index=False)
        return commands_df

def main():
    parser = argparse.ArgumentParser(description="Dry-run: Print bash commands for SRA download and trim (Section 2.11)")
    parser.add_argument("--sra_csv", required=True, help="Path to SRA CSV file (e.g., EV_sra.txt.csv)")
    parser.add_argument("--output_dir", default="results", help="Output directory")
    parser.add_argument("--cpus", default=48, type=int, help="Number of CPUs for parallel processing")
    args = parser.parse_args()

    analyzer = GenomicAnalyzer(args.sra_csv, args.output_dir, args.cpus)

    # Run dry-run for section 2.11
    commands_df = analyzer.section_2_11_quality_control_dry_run()
    print("Section 2.11: Bash Commands for SRA Download and Trim")
    for _, row in commands_df.iterrows():
        print(f"Genotype: {row['genotype']}, SRA: {row['sra']}, Type: {row['type']}")
        print(f"Command: {row['command']}\n")

if __name__ == "__main__":
    main()