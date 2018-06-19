for f in 500*.csv; do sed -i "s/}/,'windowSize': 500}/" ${f}; done

for f in 1000*.csv; do sed -i "s/}/,'windowSize': 1000}/" ${f}; done

for f in 1500*.csv; do sed -i "s/}/,'windowSize': 1500}/" ${f}; done

for f in 2000*.csv; do sed -i "s/}/,'windowSize': 2000}/" ${f}; done

for f in 2500*.csv; do sed -i "s/}/,'windowSize': 2500}/" ${f}; done

for f in 3000*.csv; do sed -i "s/}/,'windowSize': 3000}/" ${f}; done

for f in 3500*.csv; do sed -i "s/}/,'windowSize': 3500}/" ${f}; done

for f in 4000*.csv; do sed -i "s/}/,'windowSize': 4000}/" ${f}; done

for f in 4500*.csv; do sed -i "s/}/,'windowSize': 4500}/" ${f}; done

for f in 5000*.csv; do sed -i "s/}/,'windowSize': 5000}/" ${f}; done





python cleanup.py
