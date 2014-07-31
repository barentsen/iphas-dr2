"""Plots the IPHAS survey footprint."""
import abc
import numpy as np
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.io import fits

def equ2gal(equ):
    """Convert single Equatorial J2000d to Galactic coordinates."""
    import math as m
    from math import sin, cos, atan, asin, floor
    
    ra, dec = equ
    OB = m.radians(23.4333334);
    dec = m.radians(dec)
    ra = m.radians(ra)
    a = 27.128251 # The RA of the North Galactic Pole
    d = 192.859481 # The declination of the North Galactic Pole
    l = 32.931918 # The ascending node of the Galactic plane on the equator
    sdec = sin(dec)
    cdec = cos(dec)
    sa = sin(m.radians(a))
    ca = cos(m.radians(a))

    GT = asin(cdec * ca * cos(ra - m.radians(d)) + sdec * sa)
    GL = m.degrees(atan((sdec - sin(GT) * sa) / (cdec * sin(ra - m.radians(d)) * ca)))
    TP = sdec - sin(GT) * sa
    BT = cdec * sin(ra - m.radians(d)) * ca
    if (BT < 0):
        GL += 180
    else:
        if (TP < 0):
            GL += 360
           
    GL += l
    if (GL > 360):
        GL -= 360
   
    LG = floor(GL)
    LM = floor((GL - floor(GL)) * 60)
    LS = ((GL - floor(GL)) * 60 - LM) * 60
    GT = m.degrees(GT)

    D = abs(GT)
    if (GT > 0):
        BG = floor(D)
    else:
        BG = -1*floor(D)
        
    BM = floor((D - floor(D)) * 60)
    BS = ((D - floor(D)) * 60 - BM) * 60
    if (GT < 0):
        BM = -BM
        BS = -BS
    
    #if GL > 180:
    #    GL -= 360
    return [GL,GT]
    
def equ2gal_array(equ):
    """Convert Equatorial J2000d to Galactic coordinates."""
    out = []
    for e in equ:
        out.append( equ2gal(e) )
    return np.array(out)


class Survey(object):
    """Common base class for photometric surveys."""
    
    def __init__(self, name):
        self.name = name
    
    
    @abc.abstractmethod
    def get_pointings(self):
        """Retrieve ndarray of survey pointings giving name/ra/dec/etc"""
        return
    
    def _footprint(self, pixel_polygon, target, coordsys="J2000d"):
        """Convert a pixel-space polygon into a sky footprint.""" 
        if coordsys != "J2000d" and coordsys != "Gal":
            raise NotImplementedError("Unknown coordinate system "+coordsys)
        
        # Converts pixel space to WCS-defined sky system
        wcs = self.get_footprint_wcs(target)
        sky = wcs.wcs_pix2world(pixel_polygon, 1)
        # We assume WCS is in J2000d, apply conversion if Galactic is requested
        if coordsys == "Gal":
            for i, coord in enumerate(sky):
                sky[i] = equ2gal(coord)
        return sky
    
    
    def footprint(self, target, coordsys="J2000d"):
        """Lists camera footprint corners for a given target pointing."""
        return self._footprint(self._CAMERA_FOOTPRINT, target, coordsys)     
        
    
    def ccd_footprint(self, target, coordsys="J2000d"):
        """Lists CCD chip footprint corners; returns dictionary with CCDs.""" 
        out = {}
        
        for name in self._CAMERA_CCD_FOOTPRINTS.keys():
            polygon = self._CAMERA_CCD_FOOTPRINTS[name]
            out[name] = self._footprint(polygon, target, coordsys)     
        return out
        
    def footprint_patch(self, target, coordsys="J2000d", **args):
        """Create a matplotlib patch for a pointing."""
        from matplotlib import patches
        points = self.footprint(target, coordsys)
        return patches.Polygon(points, **args)
    
    def ccd_footprint_patches(self, target, coordsys="J2000d", **args):
        """Create a list of matplotlib patches for each CCD in a footprint."""
        from matplotlib import patches
        polygons = self.ccd_footprint(target, coordsys)
        out = []
        for name in polygons.keys():
            out.append( patches.Polygon(polygons[name], **args) )
        return out
  

class FootprintPlot():  
    
    def __init__(self, xlim, ylim):
        self.fig = plt.figure()
        self.xlim = xlim
        self.ylim = ylim
        self.prepare_axes()
    
    def savefig(self, filename, dpi=200):
        self.fig.savefig(filename, dpi=dpi)

    def prepare_axes(self):
        """Prepare a matplotlib figure() object for a survey footprint"""
        for i, lim in enumerate(self.xlim): 
            # Axes for subplot
            ax = self.fig.add_subplot(len(self.xlim), 1, i+1)
            # Limits
            ax.set_xlim(lim)
            ax.set_ylim(self.ylim)
            # Show degree symbol on tick labels
            
            if abs(lim[1] - lim[0]) < 10:
                xdigits = 1
            else:
                xdigits = 0
            x_formatter = plt.FormatStrFormatter(u"%%.%df\N{DEGREE SIGN}" % xdigits) 

            from matplotlib.ticker import FuncFormatter
            def glatformatter(x, pos):
                if x <= 0:
                    return u"%.1f\N{DEGREE SIGN}" % x
                else:
                    return u"%+.1f\N{DEGREE SIGN}" % x
            y_formatter = FuncFormatter(glatformatter) 

            ax.xaxis.set_major_formatter(x_formatter)
            ax.yaxis.set_major_formatter(y_formatter)


class GalacticFootprintPlot(FootprintPlot):
    
    def prepare_axes(self):
        FootprintPlot.prepare_axes(self)
        # Y label for middle plot            
        i = int(np.floor(len(self.xlim) / 2))
        #self.fig.axes[i].set_ylabel("Galactic latitude (b)")
        # X label for final subplot
        self.fig.axes[-1].set_xlabel("Galactic longitude ($l$)")
        
        # Ticks
        for ax in self.fig.axes:
            ax.xaxis.set_minor_locator(plt.MultipleLocator(base=1.0))
            ax.yaxis.set_minor_locator(plt.MultipleLocator(base=1.0))
            ax.yaxis.set_major_locator(plt.MultipleLocator(base=5.0))
                

class EquatorialFootprintPlot(FootprintPlot):
    
    def prepare_axes(self):
        FootprintPlot.prepare_axes(self)
        # Y label for middle plot
        for ax in self.fig.axes:
            ax.set_ylabel("Dec.")
        # X label for final subplot
        self.fig.axes[-1].set_xlabel("R.A.")    

        # Ticks
        for ax in self.fig.axes:
            ax.xaxis.set_major_locator(plt.MultipleLocator(base=0.25))
            ax.xaxis.set_minor_locator(plt.MultipleLocator(base=0.05))
            ax.yaxis.set_major_locator(plt.MultipleLocator(base=0.25))
            ax.yaxis.set_minor_locator(plt.MultipleLocator(base=0.05))
            


class IPHAS(Survey):
    """INT Photometric H-Alpha Survey"""
        
    def __init__(self):
        Survey.__init__(self, "IPHAS")
        
    
    def get_pointings(self):
        """Retrieve ndarray of survey pointings giving name/ra/dec/etc"""
        if not hasattr(self, "_pointings"):
            self._pointings = np.array(fits.getdata("../footprint/iphas_planner.fits", 1))
        return self._pointings
    
    # Corners of the four INT/WFC CCD detectors in pixel space.
    # These coordinates refer to values obtained after running CASU "mosaic" tool.
    _CAMERA_CCD_FOOTPRINTS = {
                    1: [ [4173.5, 4079.5], [6213.5, 4076.5], 
                         [6202.5, 9.5], [4171.5, 0.5] ],
                    2: [ [2163.5, 6151.5], [6216.5, 6148.5], 
                         [6213.5, 4123.5], [2132.5, 4171.5] ],
                    3: [ [0.5, 4114.5], [2018.5, 4126.5], 
                         [2032.5, 53.5], [21.5, 66.5] ],
                    4: [ [2065.5, 4095.5], [4109.5, 4096.5], 
                         [4104.5, 17.5], [2069.5, 22.5] ]
                    }
    
    _CAMERA_FOOTPRINT = [ _CAMERA_CCD_FOOTPRINTS[2][0], 
                          _CAMERA_CCD_FOOTPRINTS[2][1], 
                          _CAMERA_CCD_FOOTPRINTS[2][2], 
                          _CAMERA_CCD_FOOTPRINTS[1][1], 
                          _CAMERA_CCD_FOOTPRINTS[1][2], 
                          _CAMERA_CCD_FOOTPRINTS[1][3], 
                          _CAMERA_CCD_FOOTPRINTS[4][2], 
                          _CAMERA_CCD_FOOTPRINTS[4][3], 
                          _CAMERA_CCD_FOOTPRINTS[3][2], 
                          _CAMERA_CCD_FOOTPRINTS[3][3], 
                          _CAMERA_CCD_FOOTPRINTS[3][0], 
                          _CAMERA_CCD_FOOTPRINTS[3][1], 
                          _CAMERA_CCD_FOOTPRINTS[4][0],
                          _CAMERA_CCD_FOOTPRINTS[2][3] ]    

    def get_footprint_wcs(self, target):
        # Camera mosaic WCS to enable footprints
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [3108., 3075.5]
        w.wcs.cd = np.array([[ -1.24777420e-06, -9.24371710e-05],
                             [ -9.24378040e-05, 1.29159580e-06]])
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        # Offset between target coords & WCS center
        # pointing = np.array([324.341004166667, 57.6998555555556]) # Field 6686
        # crval = np.array([324.172579083079, 57.6976344549348])
        # offset = crval - pointing        
        offset = np.array([-0.16842508, -0.0022211 ])
        w.wcs.crval = target + offset
        return w

    def plot_galactic_footprint(self, filename="/tmp/iphas_footprint.png"):
        xlimits = [[230,160],[160,90], [90,20]]
        ylimit = [-7.5,+7.5]
        plot = GalacticFootprintPlot(xlimits, ylimit)
        plot.fig.suptitle("IPHAS survey footprint", fontsize=28)
        
        pointings = self.get_pointings()
        
        for i, xbin in enumerate(xlimits): 
            select = (pointings["l"] < xbin[0]+10) & (pointings["l"] > xbin[1]-10)
            for p in pointings[select]:
                target = np.array([p["ra"], p["dec"]])
                offset = np.array([-5 / 60.0, -5 / 60.0]) # offset pointing
                for mytarget in [target, target+offset]:
                    patches = self.ccd_footprint_patches(target, "Gal", 
                                                          edgecolor="None", 
                                                          facecolor="red", 
                                                          alpha=0.3)
                    for patch in patches: 
                        plot.fig.axes[i].add_patch(patch)
    
        plot.savefig(filename)


    def plot_dr2_footprint(self, filename="/tmp/iphas_footprint.png"):
        #xlimits = [[230,160],[160,90], [90,20]]
        #ylimit = [-7.,+7.]
        xlimits = [[220,120],[120,20]]
        ylimit = [-7.5, +7.5]
        
        plot = GalacticFootprintPlot(xlimits, ylimit)
        
        #pointings = self.get_pointings()
        d = fits.getdata('/home/gb/dev/iphas-qc/qcdata/iphas-qc.fits', 1)
        pointings = d[d['is_dr2']]
        
        for i, xbin in enumerate(xlimits): 
            select = (pointings["l"] < xbin[0]+10) & (pointings["l"] > xbin[1]-10)
            for p in pointings[select]:
                target = np.array([p["ra"], p["dec"]])
                offset = np.array([-5 / 60.0, -5 / 60.0]) # offset pointing
                for mytarget in [target, target+offset]:
                    patches = self.ccd_footprint_patches(target, "Gal", 
                                                          edgecolor="None",
                                                          lw=0, 
                                                          facecolor="black", 
                                                          alpha=0.2)
                    for patch in patches:
                        plot.fig.axes[i].add_patch(patch)

        ax = plot.fig.add_subplot(1, 1, 1, frameon=False)
        ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
        ax.set_ylabel('Galactic latitude ($b$)')

        plt.tight_layout(pad=0, w_pad=0.5, h_pad=1)
        plot.savefig(filename, dpi=1000)

    def plot_field_footprint(self, fieldname, filename="/tmp/iphas_field.png"):
        # Get pointing details
        p = self.get_pointings()
        select = (p["name"] == fieldname)
        p = p[select][0]
        # Plot the field
        #xlimits = [[p["ra"]+0.7, p["ra"]-1.1]]
        ylimit = [p["dec"]-0.45, p["dec"]+0.35]
        xlimits = [[104.55, 103.55]]
        plot = EquatorialFootprintPlot(xlimits, ylimit)
        target = np.array([p["ra"], p["dec"]])
        offset = np.array([-5 / 60.0, -5 / 60.0])
        for mytarget in [target, target+offset]:
            patches = self.ccd_footprint_patches(mytarget, "J2000d", 
                                                              edgecolor="None",
                                                              lw=2,
                                                              facecolor="black", 
                                                              alpha=0.5)
             
            [plot.fig.axes[0].add_patch(patch) for patch in patches]
        plot.fig.tight_layout(pad=0.4)
        plot.savefig(filename)
        

    def plot_local_footprint(self, xlim, ylim, filename="/tmp/iphas_local.png"):
        # Get pointing details
        f = self.get_pointings() 
        select = ( (f["ra"] < xlim[0]+1) & (f["ra"] > xlim[1]-1) 
              & (f["dec"] > ylim[0]-1) & (f["dec"] < ylim[1]+1) )
        
        plot = EquatorialFootprintPlot([xlim], ylim)
        for p in f[select]:
            target = np.array([p["ra"], p["dec"]])
            offset = np.array([-5 / 60.0, -5 / 60.0])
            for mytarget in [target, target+offset]:
                
                patches = self.ccd_footprint_patches(mytarget, "J2000d", 
                                                          edgecolor="None", 
                                                          facecolor="red", 
                                                          alpha=0.3)
        
                [plot.fig.axes[0].add_patch(patch) for patch in patches]
        plot.fig.suptitle("Combined IPHAS exposures", fontsize=28)
        plot.savefig(filename)   


if __name__ == "__main__":
 
    iphas = IPHAS()
    iphas.plot_field_footprint("intphas_4000", filename='fieldpair.pdf')   
    #iphas.plot_dr2_footprint(filename='fieldpair.pdf')
