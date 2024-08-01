#!/bin/bash

# Check if the input file exists
if [ ! -f "related_pairs.txt" ]; then
  echo "File not found!"
  exit 1
fi

# Create a working copy of the related pairs file
cp related_pairs.txt related_pairs_working.txt

# Extract all IDs from the file and count their occurrences
awk '{print $2; print $4}' related_pairs_working.txt | sort | uniq -c | sort -nr > id_count.txt

# Initialize an empty list for removal
> to_remove.txt

# Keep removing the individual with the highest count until no pairs remain
while [ -s related_pairs_working.txt ]; do
  # Get the individual with the highest count
  to_remove=$(head -n 1 id_count.txt | awk '{print $2}')
  
  # Add this individual to the removal list with FID and IID columns
  echo "$to_remove $to_remove" >> to_remove.txt
  
  # Remove all pairs involving this individual from related_pairs_working.txt
  awk -v remove="$to_remove" '$2 != remove && $4 != remove' related_pairs_working.txt > temp.txt
  mv temp.txt related_pairs_working.txt
  
  # Recount the occurrences of each individual
  awk '{print $2; print $4}' related_pairs_working.txt | sort | uniq -c | sort -nr > id_count.txt
done

echo "Individuals to be removed (FID IID):"
cat to_remove.txt

