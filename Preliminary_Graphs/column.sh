

for f in 500*.csv; do sed -i "s/}/,'windowSize': 500}/" ${f}; done
python cleanup.py

# filetype:500
