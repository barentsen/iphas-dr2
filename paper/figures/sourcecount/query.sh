PSQL="ionice -c3 nice -n19 psql -p 5433 -t -A -F',' -c"

$PSQL "SELECT ROUND(l) AS lonbin, COUNT(*), AVG(r-i), AVG(r-ha) FROM iphas2 WHERE b BETWEEN -5 and +5 GROUP BY lonbin ORDER BY lonbin" > count-all.csv
$PSQL "SELECT ROUND(l) AS lonbin, COUNT(*), AVG(r-i), AVG(r-ha) FROM iphas2 WHERE b BETWEEN -5 and +5 AND NOT deblend GROUP BY lonbin ORDER BY lonbin" > count-not-deblend.csv

