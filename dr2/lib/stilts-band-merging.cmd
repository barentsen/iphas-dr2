######
# The commands below are executed as the 'stilts ocmd' command after 
# band-merging the r/i/Halpha detection catalogues.
#
# It is assumed that the order of the bands given to stilts is:
# in1=r, in2=i, in3=H-alpha
######

# sourceID (use the first non-NULL value of detectionID_1/2/3)
addcol sourceID "detectionID_1"
replacecol sourceID "NULL_sourceID?detectionID_2:sourceID"
replacecol sourceID -desc "Unique identification of the band-merged detection. Identical to rDetectionID if the source was detected in the r-band; identical to iDetectionID or haDetectionID otherwise." "NULL_sourceID?detectionID_3:sourceID"

# Rename RA/Dec
colmeta -name ra ra_1
colmeta -name iRA ra_2
colmeta -name haRA ra_3
colmeta -name dec dec_1
colmeta -name iDec dec_2
colmeta -name haDec dec_3
# Fill in RA/Dec if not provided by the r-band
replacecol ra "NULL_ra?iRA:ra"
replacecol ra -desc "J2000 Right Ascension with respect to the 2MASS PSC reference frame (which is consistent with ICRS to within 0.1 arcsec). The coordinate given is obtained from the astrometric measurement in the r'-band exposure. If the source is undetected in r', then the i' or H-alpha coordinate is given."  "NULL_ra?haRA:ra"
replacecol dec "NULL_dec?iDec:dec"
replacecol dec -desc "J2000 Declination. See comments for ra."  "NULL_dec?haDec:dec"
colmeta -name posErr posErr_1
replacecol posErr "NULL_posErr?posErr_2:posErr"
replacecol posErr -desc "Astrometric fit error (RMS). Be aware that the error might be significantly larger than this near CCD edges." "NULL_posErr?posErr_3:posErr"
# Add galactic coordinates
addskycoords -inunit deg -outunit deg icrs galactic ra dec l b;

# Position offsets
addcol -desc "Position offset of the i-band detection in RA" -units "Arcsec" iXi "toFloat(3600.0*(iRA-ra))"
addcol -desc "Position offset of the i-band detection in DEC" -units "Arcsec" iEta "toFloat(3600.0*(iDec-dec))"
addcol -desc "Position offset of the Ha-band detection in RA" -units "Arcsec" haXi "toFloat(3600.0*(haRA-ra))"
addcol -desc "Position offset of the Ha-band detection in DEC" -units "Arcsec" haEta "toFloat(3600.0*(haDec-dec))"

# Rename r-band columns
colmeta -name r -desc "Default r-band magnitude using a 2.3 arcsec aperture diameter. Calibrated in the Vega system." aperMag2_1
colmeta -name rErr -desc "Uncertainty for r. Does not include systematic errors." aperMag2Err_1
colmeta -name rPeakMag -desc "Alternative r-band magnitude derived from the peak pixel height (i.e. a 0.3x0.3 arcsec square aperture). Calibrated in the Vega system." peakMag_1
colmeta -name rPeakMagErr -desc "Uncertainty in rPeakMag. Does not include systematics." peakMagErr_1
colmeta -name rAperMag1 -desc "Alternative r-band magnitude using a 1.2 arcsec aperture diameter. Calibrated in the Vega system." aperMag1_1
colmeta -name rAperMag1Err -desc "Uncertainty in rAperMag1. Does not include systematics." aperMag1Err_1
colmeta -name rAperMag3 -desc "Alternative r-band magnitude using a 3.3 arcsec aperture diameter. Calibrated in the Vega system." aperMag3_1
colmeta -name rAperMag3Err -desc "Uncertainty in rAperMag3. Does not include systematics." aperMag3Err_1
colmeta -name rGauSig -desc "RMS of axes of ellipse fit in r." gauSig_1
colmeta -name rEll -desc "Ellipticity in the r-band." ell_1
colmeta -name rPA -desc "Position angle in the r-band." pa_1
colmeta -name rClass -desc "Discrete image classification flag (1=galaxy, 0=noise, -1=star, -2=probableStar, -3=probableGalaxy, -9=saturated)." class_1
colmeta -name rClassStat -desc "N(0,1) stellarness-of-profile statistic." classStat_1
colmeta -name rBrightNeighb brightNeighb_1
colmeta -name rDeblend deblend_1
colmeta -name rSaturated saturated_1
colmeta -name rVignetted vignetted_1
colmeta -name rTruncated truncated_1
colmeta -name rBadPix badPix_1
colmeta -name rErrBits -desc "Bitmask used to flag a bright neighbour (1), source blending (2), saturation (8), vignetting (64), truncation (128) and bad pixels (32768)." errBits_1
colmeta -name rMJD -desc "Modified Julian Date at the start of the r-band exposure." mjd_1
colmeta -name rSeeing -desc "Average FWHM in the r-band exposure." seeing_1
colmeta -name rDetectionID -desc "Unique identifier of the r-band detection. Composed of the INT telescope run number, CCD number and a sequential source detection number, separated by dashes." detectionID_1
colmeta -name rX -desc "Pixel coordinate of the source in the r-band exposure, in the coordinate system of the CCD." x_1
colmeta -name rY -desc "Pixel coordinate of the source in the r-band exposure, in the coordinate system of the CCD." y_1
colmeta -name rPlaneX planeX_1
colmeta -name rPlaneY planeY_1

# Rename i-band columns
colmeta -name i -desc "Default i-band magnitude using a 2.3 arcsec aperture diameter. Calibrated in the Vega system." aperMag2_2
colmeta -name iErr -desc "Uncertainty for i. Does not include systematic errors." aperMag2Err_2
colmeta -name iPeakMag -desc "Alternative i-band magnitude derived from the peak pixel height (i.e. a 0.3x0.3 arcsec square aperture). Calibrated in the Vega system." peakMag_2
colmeta -name iPeakMagErr -desc "Uncertainty in iPeakMag. Does not include systematics." peakMagErr_2
colmeta -name iAperMag1 -desc "Alternative i-band magnitude using a 1.2 arcsec aperture diameter. Calibrated in the Vega system." aperMag1_2
colmeta -name iAperMag1Err -desc "Uncertainty in iAperMag1. Does not include systematics." aperMag1Err_2
colmeta -name iAperMag3 -desc "Alternative i-band magnitude using a 3.3 arcsec aperture diameter. Calibrated in the Vega system." aperMag3_2
colmeta -name iAperMag3Err -desc "Uncertainty in iAperMag3. Does not include systematics." aperMag3Err_2
colmeta -name iGauSig -desc "RMS of axes of ellipse fit in i." gauSig_2
colmeta -name iEll -desc "Ellipticity in the i-band." ell_2
colmeta -name iPA -desc "Position angle in the i-band." pa_2
colmeta -name iClass -desc "Discrete image classification flag (1=galaxy, 0=noise, -1=star, -2=probableStar, -3=probableGalaxy, -9=saturated)." class_2
colmeta -name iClassStat -desc "N(0,1) stellarness-of-profile statistic." classStat_2
colmeta -name iBrightNeighb brightNeighb_2
colmeta -name iDeblend deblend_2
colmeta -name iSaturated saturated_2
colmeta -name iVignetted vignetted_2
colmeta -name iTruncated truncated_2
colmeta -name iBadPix badPix_2
colmeta -name iErrBits -desc "Bitmask used to flag a bright neighbour (1), source blending (2), saturation (8), vignetting (64), truncation (128) and bad pixels (32768)." errBits_2
colmeta -name iMJD -desc "Modified Julian Date at the start of the i-band exposure." mjd_2
colmeta -name iSeeing -desc "Average FWHM in the i-band exposure." seeing_2
colmeta -name iDetectionID -desc "Unique identifier of the i-band detection. Composed of the INT telescope run number, CCD number and a sequential source detection number, separated by dashes." detectionID_2
colmeta -name iX -desc "Pixel coordinate of the source in the i-band exposure, in the coordinate system of the CCD." x_2
colmeta -name iY -desc "Pixel coordinate of the source in the i-band exposure, in the coordinate system of the CCD." y_2
colmeta -name iPlaneX planeX_2
colmeta -name iPlaneY planeY_2

# Rename H-alpha columns
colmeta -name ha -desc "Default H-alpha magnitude using a 2.3 arcsec aperture diameter. Calibrated in the Vega system." aperMag2_3
colmeta -name haErr -desc "Uncertainty for ha. Does not include systematic errors." aperMag2Err_3
colmeta -name haPeakMag -desc "Alternative H-alpha magnitude derived from the peak pixel height (i.e. a 0.3x0.3 arcsec square aperture). Calibrated in the Vega system." peakMag_3
colmeta -name haPeakMagErr -desc "Uncertainty in haPeakMag. Does not include systematics." peakMagErr_3
colmeta -name haAperMag1 -desc "Alternative H-alpha magnitude using a 1.2 arcsec aperture diameter. Calibrated in the Vega system." aperMag1_3
colmeta -name haAperMag1Err -desc "Uncertainty in haAperMag1. Does not include systematics." aperMag1Err_3
colmeta -name haAperMag3 -desc "Alternative H-alpha magnitude using a 3.3 arcsec aperture diameter. Calibrated in the Vega system." aperMag3_3
colmeta -name haAperMag3Err -desc "Uncertainty in haAperMag3. Does not include systematics." aperMag3Err_3
colmeta -name haGauSig -desc "RMS of axes of ellipse fit in H-alpha." gauSig_3
colmeta -name haEll -desc "Ellipticity in H-alpha." ell_3
colmeta -name haPA -desc "Position angle in H-alpha." pa_3
colmeta -name haClass -desc "Discrete image classification flag (1=galaxy, 0=noise, -1=star, -2=probableStar, -3=probableGalaxy, -9=saturated)." class_3
colmeta -name haClassStat -desc "N(0,1) stellarness-of-profile statistic." classStat_3
colmeta -name haBrightNeighb brightNeighb_3
colmeta -name haDeblend deblend_3
colmeta -name haSaturated saturated_3
colmeta -name haVignetted vignetted_3
colmeta -name haTruncated truncated_3
colmeta -name haBadPix badPix_3
colmeta -name haErrBits -desc "Bitmask used to flag a bright neighbour (1), source blending (2), saturation (8), vignetting (64), truncation (128) and bad pixels (32768)." errBits_3
colmeta -name haMJD -desc "Modified Julian Date at the start of the H-alpha exposure." mjd_3
colmeta -name haSeeing -desc "Average FWHM in the H-alpha exposure."  seeing_3
colmeta -name haDetectionID -desc "Unique identifier of the H-alpha detection. Composed of the INT telescope run number, CCD number and a sequential source detection number, separated by dashes." detectionID_3
colmeta -name haX -desc "Pixel coordinate of the source in the H-alpha exposure, in the coordinate system of the CCD." x_3
colmeta -name haY -desc "Pixel coordinate of the source in the H-alpha exposure, in the coordinate system of the CCD." y_3
colmeta -name haPlaneX planeX_3
colmeta -name haPlaneY planeY_3


# Single-band pStar
addcol rPStar "toFloat(NULL_rClass?1:rClass==-9?0.00:rClass==-3?0.25:rClass==-2?0.70:rClass==-1?0.90:rClass==0?0.05:rClass==1?0.05:NULL)"
addcol iPStar "toFloat(NULL_iClass?1:iClass==-9?0.00:iClass==-3?0.25:iClass==-2?0.70:iClass==-1?0.90:iClass==0?0.05:iClass==1?0.05:NULL)"
addcol haPStar "toFloat(NULL_haClass?1:haClass==-9?0.00:haClass==-3?0.25:haClass==-2?0.70:haClass==-1?0.90:haClass==0?0.05:haClass==1?0.05:NULL)"

# Single-band pGalaxy
addcol rPGalaxy "toFloat(NULL_rClass?1:rClass==-9?0.00:rClass==-3?0.70:rClass==-2?0.25:rClass==-1?0.05:rClass==0?0.05:rClass==1?0.90:NULL)"
addcol iPGalaxy "toFloat(NULL_iClass?1:iClass==-9?0.00:iClass==-3?0.70:iClass==-2?0.25:iClass==-1?0.05:iClass==0?0.05:iClass==1?0.90:NULL)"
addcol haPGalaxy "toFloat(NULL_haClass?1:haClass==-9?0.00:haClass==-3?0.70:haClass==-2?0.25:haClass==-1?0.05:haClass==0?0.05:haClass==1?0.90:NULL)"

# Single-band pNoise
addcol rPNoise "toFloat(NULL_rClass?1:rClass==0?0.90:0.05)"
addcol iPNoise "toFloat(NULL_iClass?1:iClass==0?0.90:0.05)"
addcol haPNoise "toFloat(NULL_haClass?1:haClass==0?0.90:0.05)"

# Single-band pSaturate
addcol rPSaturated "toFloat(NULL_rClass?1:rClass==-9?0.95:0.00)"
addcol iPSaturated "toFloat(NULL_iClass?1:iClass==-9?0.95:0.00)"
addcol haPSaturated "toFloat(NULL_haClass?1:haClass==-9?0.95:0.00)"

# Merged class flag
addcol pStarProd   "(rPStar * iPStar * haPStar)"
addcol pGalaxyProd "(rPGalaxy * iPGalaxy * haPGalaxy)"
addcol pNoiseProd  "(rPNoise * iPNoise * haPNoise)"
addcol pSaturatedProd  "(rPSaturated * iPSaturated * haPSaturated)"

addcol -desc "Probability the source is stellar." pStar   "pStarProd /(pStarProd+pGalaxyProd+pNoiseProd+pSaturatedProd)"
addcol -desc "Probability the source is a galaxy." pGalaxy "pGalaxyProd /(pStarProd+pGalaxyProd+pNoiseProd+pSaturatedProd)"
addcol -desc "Probability the source is noise." pNoise "pNoiseProd /(pStarProd+pGalaxyProd+pNoiseProd+pSaturatedProd)"
addcol -desc "Probability the source is saturated." pSaturated "pSaturatedProd / (pStarProd+pGalaxyProd+pNoiseProd+pSaturatedProd)"

# Assign the actual mergedClass flag, following the procedure described at
# http://surveys.roe.ac.uk/wsa/www/gloss_m.html#gpssource_mergedclass
addcol mergedClass -desc "See http://surveys.roe.ac.uk/wsa/www/gloss_m.html#gpssource_mergedclass" "toShort( pStar>0.899?-1:pStar>=0.699?-2:pGalaxy>=0.899?1:pGalaxy>=0.699?-3:pNoise>=0.899?0:pSaturated>=0.899?-9:NULL )"

# mergedClassStat
addcol mergedClassStat -desc "See http://surveys.roe.ac.uk/wsa/www/gloss_m.html#gpssource_mergedclassstat" "toFloat( mean(array(rClassStat, iClassStat, haClassStat)) * sqrt(count(array(rClassStat, iClassStat, haClassStat))) )"


# merged quality flags
addcol brightNeighb  -desc "True if a very bright star is nearby."  "(NULL_rBrightNeighb ? false : rBrightNeighb) || (NULL_iBrightNeighb ? false : iBrightNeighb) || (NULL_haBrightNeighb ? false : haBrightNeighb)"

addcol deblend  -desc "True if the source was blended with a nearby neighbour."  "(NULL_rDeblend ? false : rDeblend) || (NULL_iDeblend ? false : iDeblend) || (NULL_haDeblend ? false : haDeblend)"

addcol saturated  -desc "True if saturated in one or more bands."  "(NULL_rSaturated ? false : rSaturated) || (NULL_iSaturated ? false : iSaturated) || (NULL_haSaturated ? false : haSaturated) || r<13 || i<12 || ha<12.5"

addcol vignetted  -desc "True if in part of focal plane where image quality is poor."  "(NULL_rVignetted ? false : rVignetted) || (NULL_iVignetted ? false : iVignetted) || (NULL_haVignetted ? false : haVignetted)"

addcol truncated  -desc "True if close to the CCD boundary."  "(NULL_rTruncated ? false : rTruncated) || (NULL_iTruncated ? false : iTruncated) || (NULL_haTruncated ? false : haTruncated)"

addcol badPix  -desc "True if one or more bad pixel(s) in the aperture."  "(NULL_rBadPix ? false : rBadPix>=1) || (NULL_iBadPix ? false : iBadPix>=1) || (NULL_haBadPix ? false : haBadPix>=1)"

addcol errBits -desc "Bitmask indicating bright neighbour (1), source blending (2), saturation (8), vignetting (64), truncation (128) and bad pixels (32768)." "toInteger( maximum(array(NULL_rErrBits ? 0 : rErrBits, NULL_iErrBits ? 0 : iErrBits, NULL_haErrBits ? 0 : haErrBits)) )"

# Number of bands with a detection
addcol nBands -desc "Number of bands in which the source is detected." "toShort(sum(array(NULL_r?0:1,NULL_i?0:1,NULL_ha?0:1)))"

# Reliable
addcol reliable -desc "True if the source is detected in all three bands with SNR better than 10-sigma, errBits <= 2, and the AperMag1 magnitudes are consistent with the default magnitudes." "errBits <= 2 & nBands == 3 & r>13 & i>12 & ha>12.5 & rErr<0.1 & iErr<0.1 & haErr<0.1 & (abs(r-rAperMag1) < 3*hypot(rErr,rAperMag1Err)+0.03) & (abs(i-iAperMag1) < 3*hypot(iErr,iAperMag1Err)+0.03) & (abs(ha-haAperMag1) < 3*hypot(haErr,haAperMag1Err)+0.03)"

# night
colmeta -name night night_1
replacecol night "NULL_night?night_2:night"
replacecol -desc "Night of the observation (YYYYMMDD)." night "NULL_night?night_3:night"

# seeing
addcol -desc "Worst seeing amongst the three bands." seeing "toFloat( maximum(array(NULL_rSeeing ? 0 : rSeeing, NULL_iSeeing ? 0 : iSeeing, NULL_haSeeing ? 0 : haSeeing)) )"

# fieldID
addcol fieldID "param$fieldID"
addcol fieldGrade "param$fieldGrade"

# CCD column
colmeta -name ccd ccd_1
replacecol ccd "NULL_ccd?ccd_2:ccd"
replacecol ccd -desc "CCD-chip number on the INT Wide Field Camera of the detection." "NULL_ccd?ccd_3:ccd"

# Colours
addcol -desc "(r' - i') colour" -units "Magnitude" rmi "r - i"
addcol -desc "(r' - Ha) colour" -units "Magnitude" rmha "r - ha"

# Distance from optical axis
addcol rAxisDist "sqrt(pow(rPlaneX,2)+pow(rPlaneY,2))"
addcol iAxisDist "sqrt(pow(iPlaneX,2)+pow(iPlaneY,2))"
addcol haAxisDist "sqrt(pow(haPlaneX,2)+pow(haPlaneY,2))"
addcol -desc "Distance from the optical axis" rAxis "toFloat( maximum(array(NULL_rAxisDist?0:rAxisDist, NULL_iAxisDist?0:iAxisDist, NULL_haAxisDist?0:haAxisDist)) )"

# Remove obsolete columns
keepcols 'ra dec sourceID posErr l b mergedClass mergedClassStat pStar pGalaxy pNoise pSaturated rmi rmha r rErr rPeakMag rPeakMagErr rAperMag1 rAperMag1Err rAperMag3 rAperMag3Err rGauSig rEll rPA rClass rClassStat rErrBits rMJD rSeeing rDetectionID rX rY rPlaneX rPlaneY i iErr iPeakMag iPeakMagErr iAperMag1 iAperMag1Err iAperMag3 iAperMag3Err iGauSig iEll iPA iClass iClassStat iErrBits iMJD iSeeing iDetectionID iX iY iPlaneX iPlaneY iXi iEta ha haErr haPeakMag haPeakMagErr haAperMag1 haAperMag1Err haAperMag3 haAperMag3Err haGauSig haEll haPA haClass haClassStat haErrBits haMJD haSeeing haDetectionID haX haY haPlaneX haPlaneY haXi haEta brightNeighb deblend saturated vignetted truncated badPix errBits nBands reliable fieldID fieldGrade night seeing rAxis ccd'

