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
replacecol sourceID "NULL_sourceID?detectionID_3:sourceID"

# Rename RA/Dec
colmeta -name ra ra_1
colmeta -name iRA ra_2
colmeta -name haRA ra_3
colmeta -name dec dec_1
colmeta -name iDec dec_2
colmeta -name haDec dec_3
# Fill in RA/Dec if not provided by the r-band
replacecol ra "NULL_ra?iRA:ra"
replacecol ra "NULL_ra?haRA:ra"
replacecol dec "NULL_dec?iDec:dec"
replacecol dec "NULL_dec?haDec:dec"
colmeta -name posErr posErr_1
replacecol posErr "NULL_posErr?posErr_2:posErr"
replacecol posErr "NULL_posErr?posErr_3:posErr"
# Add galactic coordinates
addskycoords -inunit deg -outunit deg icrs galactic ra dec l b;

# Position offsets
addcol iXi "toFloat(3600.0*(iRA-ra))"
addcol iEta "toFloat(3600.0*(iDec-dec))"
addcol haXi "toFloat(3600.0*(haRA-ra))"
addcol haEta "toFloat(3600.0*(haDec-dec))"

# Rename r-band columns
colmeta -name r aperMag2_1
colmeta -name rErr aperMag2Err_1
colmeta -name rPeakMag peakMag_1
colmeta -name rPeakMagErr peakMagErr_1
colmeta -name rAperMag3 aperMag3_1
colmeta -name rAperMag3Err aperMag3Err_1
colmeta -name rGauSig gauSig_1
colmeta -name rEll ell_1
colmeta -name rPA pa_1
colmeta -name rClass class_1
colmeta -name rClassStat classStat_1
colmeta -name rBadPix badPix_1
colmeta -name rDeblend deblend_1
colmeta -name rSaturated saturated_1
colmeta -name rTruncated truncated_1
colmeta -name rBrightNeighb brightNeighb_1
colmeta -name rMJD mjd_1
colmeta -name rSeeing seeing_1
colmeta -name rDetectionID detectionID_1
colmeta -name rX x_1
colmeta -name rY y_1
colmeta -name rPlaneX Xn_1
colmeta -name rPlaneY Xi_1

# Rename i-band columns
colmeta -name i aperMag2_2
colmeta -name iErr aperMag2Err_2
colmeta -name iPeakMag peakMag_2
colmeta -name iPeakMagErr peakMagErr_2
colmeta -name iAperMag3 aperMag3_2
colmeta -name iAperMag3Err aperMag3Err_2
colmeta -name iGauSig gauSig_2
colmeta -name iEll ell_2
colmeta -name iPA pa_2
colmeta -name iClass class_2
colmeta -name iClassStat classStat_2
colmeta -name iBadPix badPix_2
colmeta -name iDeblend deblend_2
colmeta -name iSaturated saturated_2
colmeta -name iTruncated truncated_2
colmeta -name iBrightNeighb brightNeighb_2
colmeta -name iMJD mjd_2
colmeta -name iSeeing seeing_2
colmeta -name iDetectionID detectionID_2
colmeta -name iX x_2
colmeta -name iY y_2
colmeta -name iPlaneX Xn_2
colmeta -name iPlaneY Xi_2

# Rename H-alpha columns
colmeta -name ha aperMag2_3
colmeta -name haErr aperMag2Err_3
colmeta -name haPeakMag peakMag_3
colmeta -name haPeakMagErr peakMagErr_3
colmeta -name haAperMag3 aperMag3_3
colmeta -name haAperMag3Err aperMag3Err_3
colmeta -name haGauSig gauSig_3
colmeta -name haEll ell_3
colmeta -name haPA pa_3
colmeta -name haClass class_3
colmeta -name haClassStat classStat_3
colmeta -name haBadPix badPix_3
colmeta -name haDeblend deblend_3
colmeta -name haSaturated saturated_3
colmeta -name haTruncated truncated_3
colmeta -name haBrightNeighb brightNeighb_3
colmeta -name haMJD mjd_3
colmeta -name haSeeing seeing_3
colmeta -name haDetectionID detectionID_3
colmeta -name haX x_3
colmeta -name haY y_3
colmeta -name haPlaneX Xn_3
colmeta -name haPlaneY Xi_3


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

addcol pStar   "pStarProd /(pStarProd+pGalaxyProd+pNoiseProd+pSaturatedProd)"
addcol pGalaxy "pGalaxyProd /(pStarProd+pGalaxyProd+pNoiseProd+pSaturatedProd)"
addcol pNoise  "pNoiseProd /(pStarProd+pGalaxyProd+pNoiseProd+pSaturatedProd)"
addcol pSaturated "pSaturatedProd / (pStarProd+pGalaxyProd+pNoiseProd+pSaturatedProd)"

# Assign the actual mergedClass flag, following the procedure described at
# http://surveys.roe.ac.uk/wsa/www/gloss_m.html#gpssource_mergedclass
addcol mergedClass "toShort( pStar>0.899?-1:pStar>=0.699?-2:pGalaxy>=0.899?1:pGalaxy>=0.699?-3:pNoise>=0.899?0:pSaturated>=0.899?-9:NULL )"

# mergedClassStat
addcol mergedClassStat "toFloat( mean(array(rClassStat, iClassStat, haClassStat)) * sqrt(count(array(rClassStat, iClassStat, haClassStat))) )"


# Vignetted if > 3700 pixels from rotation axis
addcol rVignetted "sqrt(pow(rPlaneX,2)+pow(rPlaneY,2)) > 3700"
addcol iVignetted "sqrt(pow(iPlaneX,2)+pow(iPlaneY,2)) > 3700"
addcol haVignetted "sqrt(pow(haPlaneX,2)+pow(haPlaneY,2)) > 3700"
addcol vignetted "(NULL_rVignetted ? false : rVignetted) || (NULL_iVignetted ? false : iVignetted) || (NULL_haVignetted ? false : haVignetted)"

# Quality flags
addcol badPix "(NULL_rBadPix ? false : rBadPix>=1) || (NULL_iBadPix ? false : iBadPix>=1) || (NULL_haBadPix ? false : haBadPix>=1)"

addcol deblend "(NULL_rDeblend ? false : rDeblend) || (NULL_iDeblend ? false : iDeblend) || (NULL_haDeblend ? false : haDeblend)"

addcol saturated "(NULL_rSaturated ? false : rSaturated) || (NULL_iSaturated ? false : iSaturated) || (NULL_haSaturated ? false : haSaturated)"

addcol truncated "(NULL_rTruncated ? false : rTruncated) || (NULL_iTruncated ? false : iTruncated) || (NULL_haTruncated ? false : haTruncated)"

addcol brightNeighb "(NULL_rBrightNeighb ? false : rBrightNeighb) || (NULL_iBrightNeighb ? false : iBrightNeighb) || (NULL_haBrightNeighb ? false : haBrightNeighb)"

addcol reliable "(NULL_rErr?false:rErr<0.198) & (NULL_iErr?false:iErr<0.198) & (NULL_haErr?false:haErr<0.198) & ! saturated & ! truncated & ! brightNeighb & ! badPix"

addcol reliableStar "reliable & (NULL_mergedClass ? false: ((mergedClass < -0.5) & (mergedClass > -2.5)))"

# night
colmeta -name night night_1
replacecol night "NULL_night?night_2:night"
replacecol night "NULL_night?night_3:night"

# fieldID
addcol fieldID "param$fieldID"


# Remove obsolete columns
keepcols 'sourceID ra dec posErr l b mergedClass mergedClassStat pStar pGalaxy pNoise pSaturated r rErr rPeakMag rPeakMagErr rAperMag3 rAperMag3Err rGauSig rEll rPA rClass rClassStat rBadPix rDeblend rSaturated rTruncated rBrightNeighb rMJD rSeeing rDetectionID rX rY rPlaneX rPlaneY i iErr iPeakMag iPeakMagErr iAperMag3 iAperMag3Err iGauSig iEll iPA iClass iClassStat iBadPix iDeblend iSaturated iTruncated iBrightNeighb iMJD iSeeing iDetectionID iX iY iPlaneX iPlaneY iXi iEta ha haErr haPeakMag haPeakMagErr haAperMag3 haAperMag3Err haGauSig haEll haPA haClass haClassStat haBadPix haDeblend haSaturated haTruncated haBrightNeighb haMJD haSeeing haDetectionID haX haY haPlaneX haPlaneY haXi haEta badPix deblend saturated truncated vignetted brightNeighb reliable reliableStar night fieldID'