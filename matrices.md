This script contains work I did to process matrices for GEO database to remove unnecwsary samples (300->224)

Filtering matrices
columns_to_remove.tsv has column numbers that I want to remove

```bash
columns_to_remove=$(awk '{print $1-1}' columns_to_remove.tsv | tr '\n' ' ')
awk_command=$(printf 'BEGIN{FS=OFS="\t"; %s} {for(i=1; i<=NF; i++) if(!(i in remove)) printf $i (i<NF?OFS:ORS)}' "$(printf 'remove[%d];' ${columns_to_remove})")
awk "$awk_command" GSE199322_Matrix_processed.txt > GSE199322_Matrix_processed_filtered.txt
```

```bash
columns_to_remove_signals=$(awk '{print $1-1}' columns_to_remove_signals.tsv | tr '\n' ' ')
awk_command=$(printf 'BEGIN{FS=OFS="\t"; %s} {for(i=1; i<=NF; i++) if(!(i in remove)) printf $i (i<NF?OFS:ORS)}' "$(printf 'remove[%d];' ${columns_to_remove_signals})")
awk "$awk_command" GSE199322_Matrix_signals.txt > GSE199322_Matrix_signals_filtered.txt
```
