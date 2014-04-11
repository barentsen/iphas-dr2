PSQL="ionice -c3 nice -n19 psql -p 5433 -t -A -F',' -c"

for band in r i ha
do
    echo "Now querying ${band}"
    $PSQL "SELECT ROUND(${band}*2)/2 AS magbin, AVG(ABS(${band}-${band}2)), STDDEV(ABS(${band}-${band}2)) FROM iphas2 WHERE errBits < 16 AND errBits2 < 16 GROUP BY magbin ORDER BY magbin" > ${band}m${band}2_errbits16.csv
    $PSQL "SELECT ROUND(${band}*2)/2 AS magbin, AVG(ABS(${band}-${band}2)), STDDEV(ABS(${band}-${band}2)) FROM iphas2 WHERE veryReliable AND errBits2 = 0 GROUP BY magbin ORDER BY magbin" > ${band}m${band}2_errbits0_veryreliable.csv
done
