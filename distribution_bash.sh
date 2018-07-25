# To gather the unique values from a column:
awk '{a[$3]++} END { for (b in a) {printb}}' file

awk -F , '{ a[$2]++ } END { for (b in a) {print b} }' file > out.txt

# To gather the counts of those unique values
awk '{count[$3]++} END {for (word in count) print word, count[word]}' file
