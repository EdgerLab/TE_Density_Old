awk -F, '{a[$2]++} END {for (b in a) {print b}}' GO_Groupings.csv > filtered.txt
