# Get magnitude distributions from the database
PSQL="ionice -c3 nice -n19 psql -p 5433 -t -A -F',' -c"

for band in r i ha
do
    echo "Now querying "$band
    $PSQL "SELECT ROUND(${band}*2)/2 AS magbin, COUNT(*) FROM iphas2 GROUP BY magbin ORDER BY magbin;" > ${band}.csv
    $PSQL "SELECT ROUND(${band}*2)/2 AS magbin, COUNT(*) FROM iphas2 WHERE a10 GROUP BY magbin ORDER BY magbin;" > ${band}_a10.csv
    $PSQL "SELECT ROUND(${band}*2)/2 AS magbin, COUNT(*) FROM iphas2 WHERE a10point GROUP BY magbin ORDER BY magbin;" > ${band}_a10point.csv
done

