import argparse
import sys

def parse_info(info_str):
    """Parses the INFO field of a VCF record and returns a dictionary."""
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True # For flags
    return info_dict

def main():
    parser = argparse.ArgumentParser(description="Analyze and filter VCF files.")
    parser.add_argument("-i", "--vcf_file", required=True, help="Path to the input VCF file.")
    parser.add_argument("-o", "--output_vcf", help="Path to the output filtered VCF file.")
    parser.add_argument("-r", "--summary_report", help="Path to the output summary report file (tab-delimited).")
    parser.add_argument("-q", "--min_qual", type=float, default=30.0, help="Minimum variant quality score to keep a record.")
    parser.add_argument("-d", "--min_dp", type=int, default=10, help="Minimum total read depth (DP field in INFO) to keep a record.")

    args = parser.parse_args()

    if not args.output_vcf and not args.summary_report:
        parser.error("At least one of --output_vcf or --summary_report must be specified.")

    print(f"Processing VCF file: {args.vcf_file}")

    header_lines = []
    filtered_data_lines = []
    report_data = []

    try:
        with open(args.vcf_file, 'r') as vcf_reader:
            for line in vcf_reader:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("##"):
                    header_lines.append(line)
                    continue
                if line.startswith("#CHROM"):
                    header_lines.append(line) # VCF header line
                    report_header = line[1:].split("\t")[:5] + ["QUAL", "DP"] # CHROM, POS, ID, REF, ALT, QUAL, DP
                    continue

                fields = line.split('\t')
                if len(fields) < 8:
                    print(f"Warning: Skipping malformed VCF line: {line}", file=sys.stderr)
                    continue

                chrom, pos, record_id, ref, alt, qual_str, _, info_str = fields[:8]
                
                try:
                    qual = float(qual_str) if qual_str != '.' else -1.0
                except ValueError:
                    print(f"Warning: Skipping record with invalid QUAL value {qual_str} at {chrom}:{pos}", file=sys.stderr)
                    continue

                info_dict = parse_info(info_str)
                dp_value_str = info_dict.get('DP')
                dp_value = -1
                if dp_value_str is not None:
                    try:
                        dp_value = int(dp_value_str)
                    except ValueError:
                        print(f"Warning: Skipping record with invalid DP value {dp_value_str} at {chrom}:{pos}", file=sys.stderr)
                        # Continue to QUAL check, maybe it passes QUAL and we still want to report it if only DP is bad?
                        # For now, let's treat invalid DP as not passing the DP filter.
                        dp_value = -1


                passes_qual = qual >= args.min_qual
                passes_dp = dp_value >= args.min_dp
                
                if passes_qual and passes_dp:
                    filtered_data_lines.append(line)
                    report_data.append({
                        "CHROM": chrom,
                        "POS": pos,
                        "ID": record_id,
                        "REF": ref,
                        "ALT": alt, # ALT is already a string, potentially comma-separated
                        "QUAL": f"{qual:.1f}" if qual != -1.0 else ".", # Format QUAL to one decimal place for consistency
                        "DP": str(dp_value) if dp_value != -1 else "."
                    })

    except FileNotFoundError:
        print(f"Error: Input VCF file not found at {args.vcf_file}", file=sys.stderr)
        exit(1)
    except Exception as e:
        print(f"Error processing VCF file: {e}", file=sys.stderr)
        exit(1)

    if args.output_vcf:
        try:
            with open(args.output_vcf, 'w') as vcf_writer:
                for header_line in header_lines:
                    vcf_writer.write(header_line + "\n")
                for data_line in filtered_data_lines:
                    vcf_writer.write(data_line + "\n")
        except Exception as e:
            print(f"Error writing filtered VCF to {args.output_vcf}: {e}", file=sys.stderr)

    if args.summary_report:
        try:
            with open(args.summary_report, 'w') as report_file:
                # Write header, ensure it matches the expected format
                report_file.write("CHROM\tPOS\tID\tREF\tALT\tQUAL\tDP\n")
                # Write data
                for row in report_data:
                    report_file.write(f"{row['CHROM']}\t{row['POS']}\t{row['ID']}\t{row['REF']}\t{row['ALT']}\t{row['QUAL']}\t{row['DP']}\n")
        except Exception as e:
            print(f"Error writing summary report to {args.summary_report}: {e}", file=sys.stderr)

    print("Filtering complete.")
    if args.output_vcf:
        print(f"Filtered VCF written to {args.output_vcf}")
    if args.summary_report:
        print(f"Summary report written to {args.summary_report}")

if __name__ == '__main__':
    main()
