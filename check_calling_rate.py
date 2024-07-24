#!/usr/bin/env python3

# Read the related_pairs.txt file
related_pairs_file = "related_pairs.txt"
related_pairs = []
with open(related_pairs_file, 'r') as f:
    for line in f:
        cols = line.strip().split()
        related_pairs.append((cols[0], cols[1], cols[2], cols[3]))

# Read the missing_report.imiss file and store calling rates
missing_report_file = "missing_report.imiss"
calling_rates = {}
with open(missing_report_file, 'r') as f:
    next(f)  # Skip header
    for line in f:
        cols = line.strip().split()
        calling_rates[(cols[0], cols[1])] = float(cols[4])

# Find pairs with at least one individual having a calling rate below 0.2
low_call_rate_pairs = set()
for pair in related_pairs:
    if calling_rates.get((pair[0], pair[1]), 1.0) < 0.2 or calling_rates.get((pair[2], pair[3]), 1.0) < 0.2:
        low_call_rate_pairs.add(pair)

# Write the results to 0.2_low_call_rate_pihat.txt
output_file = "0.2_low_call_rate_pihat.txt"
with open(output_file, 'w') as f:
    for pair in low_call_rate_pairs:
        f.write(f"{pair[0]} {pair[1]} {pair[2]} {pair[3]}\n")

