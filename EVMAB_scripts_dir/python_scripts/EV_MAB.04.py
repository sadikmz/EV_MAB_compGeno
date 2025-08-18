#!/usr/bin/env python3
import argparse
import subprocess
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
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, 'trim_out'), exist_ok=True)

    def section_2_11_quality_control(self):
        """2.11: Download SRA reads and perform quality control/trimming for rows 2 and 3."""
        status_data = []
        for _, row in self.sra_df.iterrows():
            genotype = row['genotypes']
            sra_ids = row['sra'].split(',') if ',' in str(row['sra']) else [row['sra']]
            bioproject = row['bioproject']
            sra_pattern = '|'.join([sra_id.strip() for sra_id in sra_ids])
            sra_list = row['sra']  # Full comma-separated SRA list for download row
            # Construct esearch command for download
            download_cmd = f"esearch -db sra -query {bioproject} | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | grep -E '{sra_pattern}' | xargs -n 1 -P {self.cpus} fastq-dump --split-files --gzip --skip-technical"
            print(f"Executing download command: {download_cmd}")
            try:
                subprocess.run(download_cmd, shell=True, check=True)
                status = "Downloaded"
            except subprocess.CalledProcessError as e:
                status = f"Download Failed: {e}"
                status_data.append({'genotype': genotype, 'sra': sra_list, 'bioproject': bioproject, 'status': status})
                continue

            # Construct trim_galore command for each SRA
            for sra_id in sra_ids:
                sra_id = sra_id.strip()
                fastq_1 = os.path.join(self.output_dir, f"{sra_id}_1.fastq.gz")
                fastq_2 = os.path.join(self.output_dir, f"{sra_id}_2.fastq.gz")
                trim_cmd = f"trim_galore --cores {self.cpus} --quality 30 --fastqc --output_dir {os.path.join(self.output_dir, 'trim_out')} --paired {fastq_1} {fastq_2}"
                print(f"Executing trim command: {trim_cmd}")
                try:
                    subprocess.run(trim_cmd, shell=True, check=True)
                    status = "Trimmed"
                except subprocess.CalledProcessError as e:
                    status = f"Trim Failed: {e}"
                status_data.append({'genotype': genotype, 'sra': sra_id, 'bioproject': bioproject, 'status': status})

        # Save status to CSV for reference
        status_df = pd.DataFrame(status_data)
        status_df.to_csv(os.path.join(self.output_dir, 'sra_commands_status.csv'), index=False)
        return status_df

def main():
    parser = argparse.ArgumentParser(description="Execute: Download and trim SRA reads (Section 2.11)")
    parser.add_argument("--sra_metadata", required=True, help="Path to SRA CSV file (e.g., EV_sra.txt.csv)")
    parser.add_argument("--output_dir", default="results", help="Output directory")
    parser.add_argument("--cpus", default=48, type=int, help="Number of CPUs for parallel processing")
    args = parser.parse_args()

    analyzer = GenomicAnalyzer(args.sra_csv, args.output_dir, args.cpus)

    # Run execution for section 2.11
    status_df = analyzer.section_2_11_quality_control()
    print("Section 2.11: SRA Download and Trim Status")
    print(status_df)

if __name__ == "__main__":
    main()