#!/bin/bash

# Check if a directory was provided as an argument
TARGET_DIR=${1:-.}

echo "Searching for unzipped files in: $TARGET_DIR"

# Use 'find' to locate all files (recursively) 
# ! -name "*.gz" ensures we don't try to zip a zip
find "$TARGET_DIR" -type f ! -name "*.gz" | while read -r FILE; do
    echo "Gzipping: $FILE"
    gzip "$FILE"
done

echo "Done!"
