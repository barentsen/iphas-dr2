IPHAS Data Release Scripts
==========================

Contents:
 * detections: converts pipeline output to user-friendly single-band catalogues.
 * bandmerging: merges the r/i/H-alpha catalogues.
 * seaming: identifies duplicate detections amongst the bandmerged catalogues.
 * concatenating: combined the seamed catalogues into the final catalogue products.


 Running these scripts will produce the following data output
 
 iphas-dr2/detected

          /bandmerged

          /seamed/strip30
                 /strip40
                 ...
                 /strip210

          /concatenated
