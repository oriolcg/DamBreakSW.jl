#!/bin/bash

# Usage: ./replace_text.sh "old_text" "new_text" /path/to/directory

if [ $# -ne 3 ]; then
    echo "Usage: $0 old_text new_text directory"
    exit 1
fi

OLD_TEXT=$1
NEW_TEXT=$2
DIRECTORY=$3

# Find all files in the specified directory
for file in $(find "$DIRECTORY" -type f); do
    echo "Processing $file ..."
    # Use sed with | as a delimiter to avoid issues with special characters
    sed -i "s|$OLD_TEXT|$NEW_TEXT|g" "$file"
done

echo "Text replacement completed."
