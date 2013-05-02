IPHAS Data Release Processing Scripts
=====================================

Scope
-----
The scripts carry out four processing stages:
 1) detections: converts pipeline output to less obfuscated catalogues.
 2) bandmerging: merges the r/i/H-alpha catalogues at same-epoch pointings.
 3) seaming: identifies duplicate detections amongst the bandmerged catalogues.
 4) concatenating: combines the seamed catalogues into the final products.

Output
------
The processing scripts produce binary FITS tables in the following directories:
 
    iphas-dr2/
              detected/    [258 GB; ~96 CPU hours]
              bandmerged/  [182 GB; ~144 CPU hours]
              seamed/      [TBC]
              concatenated/
                           iphas-dr2-psc-glon025.fits
                           iphas-dr2-psc-glon030.fits
                           ...
                           iphas-dr2-psc-glon215.fits





IPHAS Seaming Algorithm
=======================

Goal
----
The vast majority of objects in IPHAS are detected more than once. This is a result of the survey pointing strategy, which shows significant overlaps to aid the photometric calibration and to fill in the gaps between the detectors. 

The 'seaming' algorithm serves to compare the different detections of each unique source, and decides which is the best one (called the 'primary detection'.)

Literature
----------
Warren et al. (2007), UKIDSS: "The seaming algorithm deals with sources that are detected in separate (two or more) multiframes that overlap. To decide which detection will be designated primary a sequence of discriminants is considered, until a preference is established: first the number of filters in which the source is detected, second the count from ppErrBits, and finally the distance of the source from the detector edge."

Hambly et al. (2008), UKIDSS: "A source is considered to be duplicated when an adjacent frame set contains a source within 0.8 arcsec using the same pairing/handshaking procedure described earlier. Briefly, the decision logic behind the choice of the best source examines each set of duplicates (there may be two or more to choose between) on a source-by-source basis. Source records that have the most complete passband coverage are favoured primarily; when two or more source records all have the same number of passband measures, the choice of primary source is based on position relative to the edges of the corresponding image (detections farthest from the edges are favoured) amongst the set of duplicates that have the fewest quality error bit flags set."

Algorithm outline
-----------------
For row in table(matched sourceIDs):
   choose best source => assign priOrSec
   check if partner was obtained in same hour => assign partnerID

For row in catalogue:
    add priOrSec + partnerID

Choosing the best detection
---------------------------
Sources are ranked in four steps, which are executed consecutively until only one winner remains. Only the best-ranked detections are carried forward after each step.

The ranking criteria are
 1. filter coverage (3 > 2 > 1);
 2. error flags;
 3. worst seeing amongst the filters (keep all detections within 20% of the best value);
 4. distance from the detector+ccd edge.

The morphological classification and the ellipticity are not taking into account.
