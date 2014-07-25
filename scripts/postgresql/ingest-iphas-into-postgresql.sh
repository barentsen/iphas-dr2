#!/bin/bash
#
# Ingest the IPHAS DR2 catalogue into a PostgreSQL database table
# Author: Geert Barentsen

# Configuration
DATAPATH="/home/gb/tmp/iphas-dr2-rc6/concatenated/full-augmented"  # Path to catalogue FITS tables
TABLENAME="iphas2"
DBNAME="gb"
DBPORT=5433
DBUSER="gb"
PSQL="psql $DBNAME --port $DBPORT --user $DBUSER "
STILTS="java -jar ../../dr2/lib/stilts.jar"

# Delete the previous version of the table
$PSQL -c "DROP TABLE IF EXISTS $TABLENAME;"

# Define the columns
$PSQL -c "CREATE TABLE $TABLENAME ( \
name text, \
ra double precision, \
dec double precision, \
sourceID text, \
posErr real, \
l double precision, \
b double precision, \
mergedClass smallint, \
mergedClassStat real, \
pStar real, \
pGalaxy real, \
pNoise real, \
rmi real, \
rmha real, \
r real, \
rErr real, \
rPeakMag real, \
rPeakMagErr real, \
rAperMag1 real, \
rAperMag1Err real, \
rAperMag3 real, \
rAperMag3Err real, \
rGauSig real, \
rEll real, \
rPA real, \
rClass smallint, \
rClassStat real, \
rDeblend boolean, \
rSaturated boolean, \
rMJD double precision, \
rSeeing real, \
rDetectionID text, \
rX real, \
rY real, \
i real, \
iErr real, \
iPeakMag real, \
iPeakMagErr real, \
iAperMag1 real, \
iAperMag1Err real, \
iAperMag3 real, \
iAperMag3Err real, \
iGauSig real, \
iEll real, \
iPA real, \
iClass smallint, \
iClassStat real, \
iDeblend boolean, \
iSaturated boolean, \
iMJD double precision, \
iSeeing real, \
iDetectionID text, \
iX real, \
iY real, \
iXi real, \
iEta real, \
ha real, \
haErr real, \
haPeakMag real, \
haPeakMagErr real, \
haAperMag1 real, \
haAperMag1Err real, \
haAperMag3 real, \
haAperMag3Err real, \
haGauSig real, \
haEll real, \
haPA real, \
haClass smallint, \
haClassStat real, \
haDeblend boolean, \
haSaturated boolean, \
haMJD double precision, \
haSeeing real, \
haDetectionID text, \
haX real, \
haY real, \
haXi real, \
haEta real, \
brightNeighb boolean, \
deblend boolean, \
saturated boolean, \
nBands smallint, \
a10 boolean, \
a10point boolean, \
fieldID text, \
fieldGrade text, \
night integer, \
seeing real, \
ccd smallint, \
nObs smallint, \
sourceID2 text, \
fieldID2 text, \
r2 real, \
rErr2 real, \
i2 real, \
iErr2 real, \
ha2 real, \
haErr2 real, \
errBits2 integer \
);"

# Insert all catalogue tables, one by one
for FILE in `find $DATAPATH -name "*.fits*" | sort`; do
    echo "Now ingesting $FILE"
    # Pipe the FITS tables as CSV into a pg COPY instruction
    ionice -c3 nice -n19 $STILTS tcat in=$FILE ofmt=csv \
        | $PSQL -c "COPY $TABLENAME FROM STDIN WITH CSV HEADER"
done

# Create indices
echo "Creating iphas2_ra_idx"
$PSQL -c "CREATE INDEX iphas2_ra_idx ON iphas2(ra);"
echo "Creating iphas2_dec_idx"
$PSQL -c "CREATE INDEX iphas2_dec_idx ON iphas2(dec);"
echo "Creating iphas2_l_idx"
$PSQL -c "CREATE INDEX iphas2_l_idx ON iphas2(l);"
echo "Creating iphas2_b_idx"
$PSQL -c "CREATE INDEX iphas2_b_idx ON iphas2(b);"
echo "Creating iphas2_a10_idx"
$PSQL -c "CREATE INDEX iphas2_a10_idx ON iphas2(a10);"

# Analyze
echo "Analyzing"
$PSQL -c "ANALYZE iphas2;"

# Clean up
$PSQL -c "VACUUM;"

