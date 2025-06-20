#!/bin/bash

# Get the directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Define paths relative to SCRIPT_DIR
SCRIPT_PATH="$SCRIPT_DIR/../analyze_vcf.py"
INPUT_VCF="$SCRIPT_DIR/test_data/input.vcf"

# Test Case 1: Filtered VCF
OUTPUT_FILTERED_VCF_ACTUAL="$SCRIPT_DIR/test_data/actual_filtered_q45_d20.vcf"
EXPECTED_FILTERED_VCF="$SCRIPT_DIR/test_data/expected_filtered_q45_d20.vcf"
MIN_QUAL_1=45
MIN_DP_1=20

# Test Case 2: Summary Report
OUTPUT_SUMMARY_REPORT_ACTUAL="$SCRIPT_DIR/test_data/actual_report_q30_d35.tsv"
EXPECTED_SUMMARY_REPORT="$SCRIPT_DIR/test_data/expected_report_q30_d35.tsv"
MIN_QUAL_2=30
MIN_DP_2=35

# Ensure the python script is available
if [ ! -f "$SCRIPT_PATH" ]; then
    echo "FAIL: Analysis script not found at $SCRIPT_PATH"
    exit 1
fi

# Run Test Case 1
echo "Running Test Case 1: Filtered VCF (QUAL >= $MIN_QUAL_1, DP >= $MIN_DP_1)"
python3 "$SCRIPT_PATH" -i "$INPUT_VCF" -q "$MIN_QUAL_1" -d "$MIN_DP_1" -o "$OUTPUT_FILTERED_VCF_ACTUAL"
if [ $? -ne 0 ]; then
    echo "FAIL: Script execution failed for Test Case 1"
    exit 1
fi

# Compare filtered VCF files
diff_output_vcf=$(diff "$OUTPUT_FILTERED_VCF_ACTUAL" "$EXPECTED_FILTERED_VCF")
if [ -z "$diff_output_vcf" ]; then
    echo "PASS: Filtered VCF matches expected output."
    TEST_1_PASSED=true
else
    echo "FAIL: Filtered VCF does not match expected output."
    echo "Differences:"
    echo "$diff_output_vcf"
    TEST_1_PASSED=false
fi

# Run Test Case 2
echo -e "\nRunning Test Case 2: Summary Report (QUAL >= $MIN_QUAL_2, DP >= $MIN_DP_2)"
python3 "$SCRIPT_PATH" -i "$INPUT_VCF" -q "$MIN_QUAL_2" -d "$MIN_DP_2" -r "$OUTPUT_SUMMARY_REPORT_ACTUAL"
if [ $? -ne 0 ]; then
    echo "FAIL: Script execution failed for Test Case 2"
    # Clean up Test Case 1 actual file before exiting
    rm -f "$OUTPUT_FILTERED_VCF_ACTUAL"
    exit 1
fi

# Compare summary report files
diff_output_report=$(diff "$OUTPUT_SUMMARY_REPORT_ACTUAL" "$EXPECTED_SUMMARY_REPORT")
if [ -z "$diff_output_report" ]; then
    echo "PASS: Summary report matches expected output."
    TEST_2_PASSED=true
else
    echo "FAIL: Summary report does not match expected output."
    echo "Differences:"
    echo "$diff_output_report"
    TEST_2_PASSED=false
fi

# Clean up actual files
echo -e "\nCleaning up actual output files..."
rm -f "$OUTPUT_FILTERED_VCF_ACTUAL"
rm -f "$OUTPUT_SUMMARY_REPORT_ACTUAL"
echo "Cleanup complete."

# Final result
echo -e "\n----- Summary -----"
if $TEST_1_PASSED; then
    echo "Test Case 1 (Filtered VCF): PASS"
else
    echo "Test Case 1 (Filtered VCF): FAIL"
fi

if $TEST_2_PASSED; then
    echo "Test Case 2 (Summary Report): PASS"
else
    echo "Test Case 2 (Summary Report): FAIL"
fi

if $TEST_1_PASSED && $TEST_2_PASSED; then
    echo -e "\nAll tests passed!"
    exit 0
else
    echo -e "\nSome tests failed."
    exit 1
fi
