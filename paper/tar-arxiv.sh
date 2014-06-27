# Warning: make sure iphas-dr2-arxiv.bbl exists and is an up-to-date bbl file
tar -cvf arxiv.tar iphas-dr2-arxiv.tex iphas-dr2-arxiv.bbl mn2e.bst mn2e.cls aas_macros.sty tables/columns.tex figures/footprint/footprint_small.png figures/depth/depth_r.pdf figures/depth/depth_i.pdf figures/depth/depth_h.pdf figures/seeing/seeing_r.pdf figures/seeing/seeing_i.pdf figures/seeing/seeing_ha.pdf figures/caldiagram/ccd-uncalibrated.pdf figures/caldiagram/ccd-calibrated.pdf figures/calibration/APASS-IPHAS-DR2_ishift.pdf figures/calibration/APASS-IPHAS-DR2_rshift.pdf figures/calibration/colourbar_apass_r.pdf figures/calibration/colourbar_apass_i.pdf figures/calibration/SDSS-IPHAS_rshift.pdf figures/calibration/colourbar_sdss_r.pdf figures/calibration/SDSS-IPHAS_ishift.pdf figures/calibration/colourbar_sdss_i.pdf figures/magdist/magdist-r.pdf figures/uncertainties/uncertainties.pdf figures/repeatability/repeatability.pdf figures/repeatability/repeatability-reliable.pdf figures/sourcecount/sourcecount.pdf figures/diagrams/ccd-* figures/diagrams/cmd-* figures/sh2-82/sh2-82-*
gzip arxiv.tar