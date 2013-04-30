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


