import qupath.lib.gui.scripting.QPEx
import qupath.lib.roi.ROIs
import qupath.lib.roi.GeometryTools
import qupath.lib.objects.PathObjects
import qupath.lib.regions.RegionRequest
import qupath.lib.regions.ImagePlane

import ij.process.ByteProcessor
import ij.process.BinaryProcessor
import ij.plugin.filter.RankFilters

import java.awt.geom.Point2D

/*
Script for creating concentric layers in an automatically detected dark rectangular implant. 

Automatic detection is based on thresholding and PCA-based rectangle fitting.
It downsamples the image to speed up calculations, converts it to grayscale (using max(R,G,B) to be able to detect fully dark areas), 
and applies a threshold followed by inversion and morphological operations to clean up the mask.

The script then fits a rectangle with a fixed long side to the detected points, 
computes layer thicknesses based on the specified radial thickness and the width of the implant at the section plane, 
and creates annotation objects for each layer with measurements of layer thickness in mm.

We assume that the implant area was created by cutting layered cylindrical implant perpendicular to the long axis.
The thickness of layers is calculated based on the offset of the cut from the centre. 

Input: 
- An image with a dark rectangular area. 
    The implant area should be darker than the surrounding tissue and should be the only dark object in the image.

Output: 
- A series of concentric rectangular layers as annotation objects, 
    with a measurement of layer thickness in mm.
Layers are named "Layer 1", "Layer 2", etc. starting from the outer layer.

Written by: Iva Svecova (isvecova47@gmail.com)
Date: 2026-04-05
*/

// ------------------------------------------------------------
// 0) Specify parameters
// ------------------------------------------------------------
// Implant parameters
int nLayers = 5 // number of layers on each side from the centre
double radialLayerThicknessMM = 0.4  // thickness of one layer in the radial direction in mm
double implantHeightMM = 7.5  // height of the implant (longer side of the ROI) in mm

// Parameters for automatic rectangle detection
double downsample = 20.0
int thresh = 5

// Outer layer thickness in mm 
outerLayerThicknessMM = 0.5 // thickness of the outer layer in mm

// ------------------------------------------------------------
// 1) load image
// ------------------------------------------------------------
def imageData = QPEx.getCurrentImageData()
def server = imageData.getServer()

// Downsampling added to speed up calculations, as we only need an approximate rectangle fit
def request = RegionRequest.createInstance(server.getPath(), downsample, 0, 0, server.getWidth(), server.getHeight())
def img = server.readRegion(request)

int w = img.getWidth()
int h = img.getHeight()

def bp = new ByteProcessor(w, h)

// Convert to grayscale using max(R,G,B) to be able to detect fully dark areas 
for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
        int rgb = img.getRGB(x, y)
        int r = (rgb >> 16) & 0xff
        int g = (rgb >> 8) & 0xff
        int b = (rgb) & 0xff
        int gray = Math.max(r, Math.max(g, b))
        bp.set(x, y, gray)
    }
}

// ------------------------------------------------------------
// 2) segmentation
// ------------------------------------------------------------

// We are interested in the dark areas - therefore we use invert
def mask = bp.duplicate()
mask.threshold(thresh)
mask.invert()

// Clean up the mask and remove border-touching objects 
def rf = new RankFilters()
rf.rank(mask, 10.0, RankFilters.MIN)

def bin = new BinaryProcessor(mask)
bin.erode()

clearBorder(bin)

// ------------------------------------------------------------
// 3) extract points
// ------------------------------------------------------------
def pts = new ArrayList<Point2D.Double>()
for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
        if (bin.get(x, y) != 0)
            pts.add(new Point2D.Double(x, y))
    }
}

if (pts.isEmpty()) {
    print "No foreground.\n"
    return
}

// ------------------------------------------------------------
// 4) calibration
// ------------------------------------------------------------
def cal = server.getPixelCalibration()
if (cal == null || !cal.hasPixelSizeMicrons()) {
    print "No calibration.\n"
    return
}

double pixelSizeMicrons = cal.getAveragedPixelSizeMicrons()
print "Pixel size: " + pixelSizeMicrons + " microns\n"

// ------------------------------------------------------------
// 5) rectangle fit
// ------------------------------------------------------------
double fixedLongPx = (implantHeightMM * 1000 / pixelSizeMicrons) / downsample
def rectInfo = fitRectWithFixedLongSide(pts, fixedLongPx)

// ------------------------------------------------------------
// 6) compute layer thicknesses 
//      We assume that the implant was cut perpendicular to the long axis,
//       but with a possible offset relative to the centre.
//
//      The layers are concentric and have constant thickness in the radial direction.
//      As the cut might be offset from the centre, the layer thickness in the plane of the section will be different from the radial thickness.
// ------------------------------------------------------------

double widthMM = (rectInfo.shortHalf * 2) * pixelSizeMicrons * downsample / 1000.0

double k = Math.sqrt(Math.pow(nLayers*radialLayerThicknessMM, 2) - Math.pow(widthMM / 2, 2)) / radialLayerThicknessMM

def layerWidthsPx = []
double inside = 0

for (int l = 1; l <= nLayers; l++) {

    if (l <= k) {
        // layerWidthsPx.add(0.0)   // No need to add zero widths, just skip
        continue
    }

    double x = Math.sqrt(Math.pow((l*radialLayerThicknessMM), 2) - Math.pow((k*radialLayerThicknessMM), 2))
    double layerWidthMM = (x - inside)

    double layerWidthPx = (layerWidthMM * 1000.0 / pixelSizeMicrons) / downsample

    layerWidthsPx.add(layerWidthPx)
    inside = x
}
layerWidthsPx = layerWidthsPx.reverse()  // reverse to start from outer layers

// ------------------------------------------------------------
// 7) build rectangles
// ------------------------------------------------------------
def rectangles = []

double currentHalf = rectInfo.shortHalf

// Omit the last layer because it will be added as the final inner rectangle rather than a difference between two rectangles
for (int i = 0; i < layerWidthsPx.size() - 1; i++) {

    double layerWidth = layerWidthsPx[i]

    if (layerWidth <= 0)
        continue

    rectangles.add([
        center   : rectInfo.center,
        ux       : rectInfo.ux,
        uy       : rectInfo.uy,
        longHalf : rectInfo.longHalf,
        shortHalf: currentHalf
    ])

    double nextHalf = currentHalf - layerWidth

    if (nextHalf <= 0)
        break

    currentHalf = nextHalf
}

// final inner rectangle
rectangles.add([
    center   : rectInfo.center,
    ux       : rectInfo.ux,
    uy       : rectInfo.uy,
    longHalf : rectInfo.longHalf,
    shortHalf: currentHalf
])

// Build also the outer rectangle based on the specified outer layer thickness
outerRectangle = [
    center   : rectInfo.center,
    ux       : rectInfo.ux,
    uy       : rectInfo.uy,
    longHalf : rectInfo.longHalf,
    shortHalf: (rectInfo.shortHalf + (outerLayerThicknessMM * 1000.0 / pixelSizeMicrons) / downsample)
]

// ------------------------------------------------------------
// 8) scale to full resolution
// ------------------------------------------------------------
def rectanglesFull = rectangles.collect { r -> scaleRectInfo(r, downsample) }
def outerRectangleFull = scaleRectInfo(outerRectangle, downsample)

// ------------------------------------------------------------
// 9) build ROIs
// ------------------------------------------------------------
int z = (server.nZSlices() > 1) ? 1 : 0     // If this is a z-stack, switch to the second plane
def plane = ImagePlane.getPlane(z, 0)

def makePolyROI = { polyPts ->
    double[] xs = new double[polyPts.size()]
    double[] ys = new double[polyPts.size()]
    for (int i = 0; i < polyPts.size(); i++) {
        xs[i] = polyPts[i].x
        ys[i] = polyPts[i].y
    }
    ROIs.createPolygonROI(xs, ys, plane)
}

def makeLayerROI = { outerROI, innerROI ->
    def gOuter = GeometryTools.roiToGeometry(outerROI)
    def gInner = GeometryTools.roiToGeometry(innerROI)
    def gRing = gOuter.difference(gInner)

    PathObjects.createAnnotationObject(
        GeometryTools.geometryToROI(gRing, plane)
    )
}

// create ROIs
def rois = rectanglesFull.collect { r ->
    def corners = rectCorners(r.center, r.ux, r.uy, r.longHalf, r.shortHalf)
    makePolyROI(corners)
}

// Create outer layer ROI based on the specified outer layer thickness
def outerCorners = rectCorners(outerRectangleFull.center, outerRectangleFull.ux, outerRectangleFull.uy, outerRectangleFull.longHalf, outerRectangleFull.shortHalf)
def outerROI = makePolyROI(outerCorners)

def outer = makeLayerROI(outerROI, rois[0])
outer.setName("Outer layer")

outer.getMeasurementList().put("Layer thickness mm", outerLayerThicknessMM)

addObject(outer)


// create layers with calculated thickness
for (int i = 0; i < rois.size() - 1; i++) {

    def obj = makeLayerROI(rois[i], rois[i + 1])
    obj.setName("Layer " + (i + 1))

    double thicknessPx = layerWidthsPx[i]
    double thicknessMM = thicknessPx * downsample * pixelSizeMicrons / 1000.0

    obj.getMeasurementList().put("Layer thickness mm", thicknessMM)

    addObject(obj)
}

// inner layer
def inner = PathObjects.createAnnotationObject(rois[-1])
inner.setName("Layer " + rois.size())

double remainingPx = currentHalf
double remainingMM = remainingPx * downsample * pixelSizeMicrons / 1000.0

inner.getMeasurementList().put("Layer thickness mm", remainingMM)

addObject(inner)

selectObjects(inner)

print "Done.\n"


// ============================================================
// Helper functions
// ============================================================

static void clearBorder(BinaryProcessor bp) {
    int w = bp.getWidth(), h = bp.getHeight()
    boolean[] visited = new boolean[w*h]
    int[] qx = new int[w*h]
    int[] qy = new int[w*h]
    int qs=0, qe=0

    def push = {x,y->
        int i=y*w+x
        if (!visited[i] && bp.get(x,y)!=0){
            visited[i]=true
            qx[qe]=x; qy[qe]=y; qe++
        }
    }

    for (int x=0;x<w;x++){ push(x,0); push(x,h-1) }
    for (int y=0;y<h;y++){ push(0,y); push(w-1,y) }

    while (qs<qe){
        int x=qx[qs], y=qy[qs]; qs++
        if(x>0)push(x-1,y)
        if(x<w-1)push(x+1,y)
        if(y>0)push(x,y-1)
        if(y<h-1)push(x,y+1)
    }

    for (int i=0;i<visited.length;i++){
        if(visited[i]){
            int x=i%w, y=(int)(i/w)
            bp.set(x,y,0)
        }
    }
}

// Fit a rectangle with a fixed long side to a set of 2D points using PCA
static def fitRectWithFixedLongSide(List<Point2D.Double> pts, double fixedLongSidePx) {
    // 1. Compute centroid (mean x, y)
    double mx = pts.collect{it.x}.sum() / pts.size()
    double my = pts.collect{it.y}.sum() / pts.size()

    // 2. Compute covariance components (for PCA)
    double sxx=0, syy=0, sxy=0
    for (p in pts) {
        double dx = p.x - mx
        double dy = p.y - my
        sxx += dx*dx      // variance in x
        syy += dy*dy      // variance in y
        sxy += dx*dy      // covariance
    }

    // 3. Calculate principal axis orientation (theta)
    double theta = 0.5 * Math.atan2(2*sxy, (sxx - syy))

    // 4. Get unit vectors for principal (u) and orthogonal (v) axes
    double ux = Math.cos(theta)
    double uy = Math.sin(theta)

    double vx = -uy
    double vy = ux

    // 5. Project points onto principal axes to find extents
    double minU=1e12, maxU=-1e12, minV=1e12, maxV=-1e12

    for (p in pts) {
        double dx = p.x - mx
        double dy = p.y - my

        double u = dx*ux + dy*uy // projection onto principal axis
        double v = dx*vx + dy*vy // projection onto orthogonal axis

        minU = Math.min(minU, u)
        maxU = Math.max(maxU, u)
        minV = Math.min(minV, v)
        maxV = Math.max(maxV, v)
    }

    double extentU = maxU - minU // length along principal axis
    double extentV = maxV - minV // length along orthogonal axis

    // 6. Ensure U is the long axis (swap if needed)
    if (extentV > extentU) {
        double tmp
        tmp = ux; ux = vx; vx = tmp
        tmp = uy; uy = vy; vy = tmp
        tmp = extentU; extentU = extentV; extentV = tmp
    }

    // 7. Compute center in (u, v) space and map back to (x, y)
    double centerU = (minU + maxU)/2
    double centerV = (minV + maxV)/2

    double cx = mx + centerU*ux + centerV*vx
    double cy = my + centerU*uy + centerV*vy

    // 8. Return rectangle parameters: center, orientation, half-lengths
    [
        center   : new Point2D.Double(cx, cy),
        ux       : ux, // unit vector x-component (long axis)
        uy       : uy, // unit vector y-component (long axis)
        longHalf : fixedLongSidePx / 2.0, // half of fixed long side
        shortHalf: extentV / 2.0          // half of fitted short side
    ]
}

static def scaleRectInfo(def r,double s){
    [
        center:new Point2D.Double(r.center.x*s,r.center.y*s),
        ux:r.ux, uy:r.uy,
        longHalf:r.longHalf*s,
        shortHalf:r.shortHalf*s
    ]
}

// Given rectangle parameters, compute the four corners in order
static List<Point2D.Double> rectCorners(center, ux, uy, longHalf, shortHalf) {
    // Compute orthogonal unit vector (v) from (u)
    double vx = -uy
    double vy = ux

    // Return corners in order (clockwise or counterclockwise)
    [
        // Bottom-left
        new Point2D.Double(center.x - longHalf*ux - shortHalf*vx, center.y - longHalf*uy - shortHalf*vy),
        // Bottom-right
        new Point2D.Double(center.x + longHalf*ux - shortHalf*vx, center.y + longHalf*uy - shortHalf*vy),
        // Top-right
        new Point2D.Double(center.x + longHalf*ux + shortHalf*vx, center.y + longHalf*uy + shortHalf*vy),
        // Top-left
        new Point2D.Double(center.x - longHalf*ux + shortHalf*vx, center.y - longHalf*uy + shortHalf*vy)
    ]
}