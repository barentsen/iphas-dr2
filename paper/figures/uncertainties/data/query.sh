# Computes the mean Poissonian uncertainties as a function of magnitude
PSQL="ionice -c3 nice -n19 psql -p 5433 -t -A -F',' -c"
for band in r i ha
do
    echo "Now querying ${band}"
    $PSQL "SELECT ROUND(${band}*2)/2 AS magbin, AVG(${band}Err), STDDEV(${band}Err) FROM iphas2 GROUP BY magbin ORDER BY magbin" > ${band}_uncertainties.csv
done
