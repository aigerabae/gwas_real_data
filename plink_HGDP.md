1) make HGDP into plink binary file with 38 build

Ped file
a) getting p1 from phenotypes metadata, encoding male and female as 1 and 2, sorting
```bash
awk '{print $1"\t" $1"\t" "0\t" "0\t" $2"\t" "-9\t"}' metadata.txt | tail -n +2 > HGDP_p1_1.txt
awk -F'\t' '{ if ($5 == "male") $5 = 1; else if ($5 == "female") $5 = 2; print }' OFS='\t' HGDP_p1_1.txt > HGDP_p1_2.txt
sort -k1,1 HGDP_p1_2.txt > HGDP_p1_3.txt
```

b) Transposing for p2 (genotypes)
```bash
split -l 90000 HGDP1.txt file_
```

I renamed them to gtReport1 to gtReport8:
```bash
mv file_aa ./gtReport1
mv file_ab ./gtReport2
mv file_ac ./gtReport3
mv file_ad ./gtReport4
mv file_ae ./gtReport5
mv file_af ./gtReport6
mv file_ag ./gtReport7
mv file_ah ./gtReport8
```

then using the code below removed spaces and replaced them with tabs, then transposed it, sorted the file, removed empty lines from genotypes file, removed all intermediate files
(takes about 1 hour to run)

```bash
process_file() {
  file=$1
  awk '{printf $1"\t"; for(i=2;i<=NF;i++){printf $i"\t"}; print ""}' "$file" > "${file}_converted.tsv"
  awk -F'\t' 'NF' "${file}_converted.tsv" > "${file}_converted_cleaned.tsv"
  datamash transpose < "${file}_converted_cleaned.tsv" > "${file}_transposed.tsv"
}
export -f process_file
parallel -j 8 process_file ::: gtReport{1..8}
for i in {1..8}; do sed -i '' 's/-/0/g' gtReport${i}_transposed.tsv; done
paste gtReport{1..8}_transposed.tsv > HGDP_transposed.txt
sort -k1,1 HGDP_transposed.txt | tail -n +2 > HGDP_p2_1.txt
sed -E 's/([^\t]+)\t/\1\t/;s/\t(.)/\t\1\t/g' HGDP_p2_1.txt > HGDP_p2_2.txt
rm *gtRe*
```

c) joining p1 and p2
```bash
join -t$'\t' HGDP_p1_3.txt HGDP_p2_2.txt | sed 's/\t\t*/\t/g' > HGDP.ped
```

d) Map file
```bash
awk '{print $2"\t" $1 "\t0\t" $3}' HGDP_Map.txt > HGDP.map
```

e) Ped/map to plink binary
```bash
plink --file HGDP --missing-code -9,0,NA,na,- --make-bed --out HGDP
```

f) remap 36.1 to 38 build (for some reason cut wasnt working properly with chromosomes; i dont know why)
```bash
comm -12 <(awk '{print $2}' kaz12_autosomal.bim | sort) <(awk '{print $2}' HGDP.bim | sort) > common_snps.txt
cat kaz12_autosomal.bim | awk '{print $2"\t" $1}' > dictionary_chr
cat kaz12_autosomal.bim | cut -f 2,4 > dictionary_pos
plink --bfile HGDP --extract common_snps.txt --make-bed --out HGDP1
plink2 --bfile HGDP1 --update-chr dictionary_chr --update-map dictionary_pos --sort-vars --make-pgen --out HGDP2
plink2 --pfile HGDP2 --make-bed --out HGDP3
```

2) QC of HGDP data
```bash
plink --bfile HGDP3 --geno 0.02 --make-bed --out HGDP4
plink --bfile HGDP4 --mind 0.02 --make-bed --out HGDP5
plink --bfile HGDP5 --maf 0.001 --make-bed --out HGDP6
plink --bfile HGDP6 --genome --min 0.2 --out pihat_min0.2
plink --bfile HGDP6 --missing --out missing_report
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
```

```bash
nano remove_relatives.sh
```

This script will remove relatives
```bash
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
```

```bash
plink --bfile HGDP6 --remove to_remove.txt --make-bed --out HGDP7
plink --bfile HGDP7 --snps-only 'just-acgt' --make-bed --out HGDP8
```

4) PCA

Only keeping selected ethicities
```bash
cut -f 1,5 metadata.txt > ethicities.txt
awk 'NR==FNR {ids[$1]; next} $1 in ids {print $1"\t" $2}' HGDP8.fam ethicities.txt > matching_ethnicities.txt
cat matching_ethnicities.txt | grep -e "Uygur" -e "Hazara" -e "Russian" -e "French" -e "Basque" -e "Bergamo" -e "Pathan" -e "Sindhi" -e "Kalash" -e "Adygei" -e "Bedouin" -e "Mozabite" -e "Japanese" -e "Northern" -e "Mongolian" -e "Yakut" -e "Han" | cut -f 1 | awk '{print $1"\t" $1}' > selected_ethnicities.txt
plink --bfile HGDP8 --keep selected_ethnicities.txt --biallelic-only strict --make-bed --out HGDP9
```

Merging kazakh and HGDP data (first - deal with multiallelic variants) and doing PCA;

Problem - need to merge it in a way that avoids the multiallelic problem
```bash

plink --bfile kaz12_autosomal --bmerge HGDP9.bed HGDP9.bim HGDP9.fam --make-bed --out merged1
plink --bfile HGDP9 --exclude merged_dataset-merge.missnp --biallelic-only strict --make-bed --out HGDP10
plink --bfile kaz12_autosomal --exclude merged_dataset-merge.missnp --biallelic-only strict --make-bed --out kaz13_autosomal
plink --bfile kaz13_autosomal --bmerge HGDP10.bed HGDP10.bim HGDP10.fam --make-bed --out merged2
plink --bfile merged2 --geno 0.02 --make-bed --out merged3
plink --bfile merged3 --mind 0.02 --make-bed --out merged4
echo -e "rs9267522\nrs11229\nrs75412302\nrs12660719" > duplicates.txt
plink --bfile merged4 --exclude duplicates.txt --make-bed --out merged5

cat kaz12_autosomal.fam | cut -f 1 -d ' ' >> temp_ethicities.txt 
cat temp_ethicities.txt | awk '{$2 = "\tKazakh"; print }' >> ethicities.txt
rm temp_ethicities.txt
```

```bash
nano outliers.txt
```

```bash
1210510 1210510
HGDP00621       HGDP00621
HGDP01270       HGDP01270
HGDP01271       HGDP01271
HGDP00175       HGDP00175
HGDP00953       HGDP00953
HGDP00949       HGDP00949
HGDP00969       HGDP00969
HGDP00959       HGDP00959
HGDP00621       HGDP00621
```

```bash
plink --bfile merged5 --remove outliers.txt --make-bed --out merged6
plink2 --bfile merged6 --pca 10 --out all_pca 
```

```bash
nano plot_eigenvec.py 
```

```python
import sys
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.spatial import ConvexHull
import numpy as np

# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python plot_eigenvec.py <eigenvec_file> <ethnicities_file>")
    sys.exit(1)

# Read the eigenvec file, skipping the first line
eigenvec_file = sys.argv[1]
with open(eigenvec_file, 'r') as f:
    lines = f.readlines()[1:]

# Read the ethnicities file
ethnicities_file = sys.argv[2]
with open(ethnicities_file, 'r') as f:
    ethnicities = {line.split()[0]: line.split()[1] for line in f.readlines()}

# Extract data from the file
data = {
    'PC1': [float(line.split()[2]) for line in lines],
    'PC2': [float(line.split()[3]) for line in lines],
    'sample_id': [line.split()[0] for line in lines],
    'ethnicity': [ethnicities[line.split()[0]] for line in lines]
}

df = pd.DataFrame(data)

# Define marker shapes and colors based on ethnicity
shape_map = {
    'Uygur': 'circle-open', 'Kazakh': 'circle-open', 'Hazara': 'circle-open',
    'Russian': 'square', 'French': 'square', 'Basque': 'square', 'Bergamo': 'square',
    'Pathan': 'triangle-up', 'Sindhi': 'triangle-up', 'Kalash': 'triangle-up',
    'Adygei': 'diamond',
    'Bedouin': 'square-open', 'Mozabite': 'square-open',
    'Japanese': 'triangle-down', 'Northern': 'triangle-down', 'Mongolian': 'triangle-down', 'Yakut': 'triangle-down', 'Han': 'triangle-down'
}

color_map = {
    'circle-open': 'rgba(150,0,100,0.6)',   # Darker Pink
    'square': 'rgba(100,100,0,0.6)',        # Darker Blue
    'triangle-up': 'rgba(100,0,150,0.6)',    # Darker Green
    'diamond': 'rgba(0,150,50,0.6)',       # Darker Yellow
    'square-open': 'rgba(50,200,100,0.6)',   # Darker Beige
    'triangle-down': 'rgba(200,200,100,0.6)'  # Darker Khaki
}

# Create interactive scatter plot
fig = px.scatter(df, x='PC1', y='PC2', color='ethnicity', symbol='ethnicity', hover_data=['sample_id', 'ethnicity'])

# Update marker shapes and colors based on ethnicity
for ethnicity, shape in shape_map.items():
    fig.update_traces(marker=dict(symbol=shape), selector=dict(name=ethnicity))

# Function to add convex hulls with custom names
def add_convex_hull(df, shape, color, legend_name):
    shape_df = df[df['ethnicity'].map(lambda x: shape_map[x] == shape)]
    if shape_df.shape[0] < 3:  # Need at least 3 points to form a convex hull
        return
    
    points = shape_df[['PC1', 'PC2']].values
    try:
        hull = ConvexHull(points)
        # Add convex hull trace
        hull_points = np.append(points[hull.vertices], [points[hull.vertices][0]], axis=0)  # close the hull
        fig.add_trace(go.Scatter(
            x=hull_points[:, 0],
            y=hull_points[:, 1],
            mode='lines',
            line=dict(color=color, width=1),  # Thinner lines
            name=legend_name  # Set custom legend name
        ))
    except:
        print(f"Warning: Convex hull computation failed for shape: {shape}")

# Add convex hulls for each shape type with custom legend names
shape_legends = {
    'circle-open': 'Central Asia',
    'square': 'Europe',
    'triangle-up': 'South Asia',
    'triangle-down': 'East Asia',
    'square-open': 'Middle East'
}

for shape, legend_name in shape_legends.items():
    add_convex_hull(df, shape, color_map[shape], legend_name)

# Update layout for larger font sizes
fig.update_layout(
    title='PCA Plot',
    title_font_size=24,
    xaxis_title='PC1',
    xaxis_title_font_size=20,
    yaxis_title='PC2',
    yaxis_title_font_size=20,
    legend_title_font_size=18,
    font=dict(size=16)
)

# Save the plot as an HTML file
fig.write_html('interactive_plot.html')

# Show the plot in a browser
fig.show()
```

```bash
python plot_eigenvec.py all_pca.eigenvec ethicities.txt
```

6) runs of homozygosity (ROH)
```bash
plink --bfile kaz12_autosomal --homozyg-density 60 --homozyg-gap 500 --homozyg-window-snp 100 --homozyg-window-het 0
plink --bfile HGDP7 --homozyg-density 60 --homozyg-gap 500 --homozyg-window-snp 100 --homozyg-window-het 0
```

tasks for thursday:
- understand whats Fst
- check how to visualize ROH
- check how to combine HGDP and kaz data
- redo PCA: include only chosen populations into the final dataset, see how many of them will be present, rewrite the python script for making the PCA graph

use all populations from HDGP and SGDP for ROH and Fst and only selected eurasian populations for PCA

Find ALDH2 gene in kazakh and other populations and see whether we absorb alcohol better or rose than other central asians or europeans