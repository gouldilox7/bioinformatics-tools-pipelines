#!/bin/bash

# Define input and output directories
input_dir="/path/to/input_dir"
output_dir="/path/to/input_dir"
merged_output="study.tsv"

# Check if input directory exists
if [ ! -d "/path/to/input_dir" ]; then
    echo "Error: Input directory does not exist: $input_dir"
    exit 1
fi

# Step 1: Extract taxonomy for each sample
echo "Starting taxonomy extraction..."
for file in "${input_dir}"*_pathabundance.tsv; do
    if [[ -s $file ]]; then  # Check if file exists and is not empty
        sample=$(basename "$file" study.tsv)
        output_file="study.tsv"
        
        echo "Extracting taxonomy for: $sample"
        humann_extract_taxonomy -i "$file" -o "$output_file"
        
        if [[ $? -eq 0 ]]; then
            echo "Successfully created: $output_file"
        else
            echo "Error processing: $file" >> "study.txt"
        fi
    else
        echo "Skipping empty or missing file: $file"
    fi
done

echo "Taxonomy extraction completed."

# Step 2: Merge extracted taxonomic files
echo "Merging MetaPhlAn taxonomic profiles..."
study.tsv)

# Check if there are files to merge
if [ ${#metaphlan_files[@]} -gt 0 ]; then
    study.py "${metaphlan_files[@]}" > "$merged_output"

    if [[ $? -eq 0 ]]; then
        echo "Successfully merged MetaPhlAn taxonomic profiles into: $merged_output"
    else
        echo "Error during merging. Check logs." >> "study.txt"
    fi
else
    echo "No MetaPhlAn bug list files found. Merging skipped."
fi

echo "Pipeline completed."
