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
                           strip30.fits   [2.1 GB]
                           strip130.fits  [0.6 GB]
                           strip140.fits  [0.4 GB]
                           strip160.fits  [0.5 GB]
                           strip200.fits  [0.5 GB]
                           strip210.fits  [0.3 GB]

Open issues
-----------

Critical:
 * Fix WCS of ~50 outstanding CCDs. [1 day]
 * Fix bad column issue in dec2003 data. [1 day]
 * Relax vignetting constraint. [0.1 day]
 * Add IPHAS object identifier strings. [0.1 day]

Necessary:
 * Automated tests (validation) of data products. [2 days?]
 * Discard fields with poor calibration. [2 days?]
 * Document the files and columns on a webpage. [1 day]

Optional:
 * Add (r2, i2, ha2) columns. [1 day]
 * Optimize bright neighbour flag. [1 day]

 Challenges
 ----------
  * Computing power with high I/O throughput.
  * Disk space on desktop.
  * Hosting 20+ GB of catalogues locally.