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
Script for creating concentric layers in a rectangular section. 
Takes in a polygonal ROI and calculates the width of each layer based on a cylindrical model. 

We assume that the rectangle was created by cutting layered cylinder perpendicular to the long axis.
The thickness of layers is calculated based on the offset of the cut from the centre. 

Input: 
- A selected polygonal annotation object representing the rectangular section.  

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

// ------------------------------------------------------------
// 1) get selected ROI
// ------------------------------------------------------------
def selected = getSelectedObject()

if (selected == null || selected.getROI() == null) {
    print "No ROI selected.\n"
    return
}

def roi = selected.getROI()

// ------------------------------------------------------------
// 2) extract points from ROI
// ------------------------------------------------------------
def points = roi.getAllPoints()
if (points == null || points.isEmpty()) {
    print "ROI has no points.\n"
    return
}

def pts = points.collect { p -> new Point2D.Double(p.getX(), p.getY()) }


// ------------------------------------------------------------
// 3) calibration
// ------------------------------------------------------------
def server = getCurrentImageData().getServer()
def cal = server.getPixelCalibration()

if (cal == null || !cal.hasPixelSizeMicrons()) {
    print "No calibration.\n"
    return
}

double pixelSizeMicrons = cal.getAveragedPixelSizeMicrons()
print "Pixel size: " + pixelSizeMicrons + " microns\n"


// ------------------------------------------------------------
// 4) get rectangle parameters (PCA)
// ------------------------------------------------------------
def rectInfo = getRectFromPoints(pts)


// ------------------------------------------------------------
// 6) compute layer thicknesses 
//      We assume that the implant was cut perpendicular to the long axis,
//       but with a possible offset relative to the centre.
//
//      The layers are concentric and have constant thickness in the radial direction.
//      As the cut might be offset from the centre, the layer thickness in the plane of the section will be different from the radial thickness.
// ------------------------------------------------------------

double widthMM = (rectInfo.shortHalf * 2) * pixelSizeMicrons / 1000.0
print "Width at section plane [mm]: " + widthMM

double k = 0
if (widthMM < 2 * nLayers * radialLayerThicknessMM) {
    k = Math.sqrt(Math.pow(nLayers*radialLayerThicknessMM, 2) - Math.pow(widthMM / 2, 2)) / radialLayerThicknessMM
} else {
    print "Warning: The specified number of layers and thickness exceeds the width of the implant at the section plane."
}

def layerWidthsPx = []
double inside = 0

for (int l = 1; l <= nLayers; l++) {

    if (l <= k) {
        // layerWidthsPx.add(0.0)   // No need to add zero widths, just skip
        continue
    }

    double x = Math.sqrt(Math.pow((l*radialLayerThicknessMM), 2) - Math.pow((k*radialLayerThicknessMM), 2))
    double layerWidthMM = (x - inside)

    double layerWidthPx = (layerWidthMM * 1000.0 / pixelSizeMicrons)

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
def rois = rectangles.collect { r ->
    def corners = rectCorners(r.center, r.ux, r.uy, r.longHalf, r.shortHalf)
    makePolyROI(corners)
}

// create layers with calculated thickness
for (int i = 0; i < rois.size() - 1; i++) {

    def obj = makeLayerROI(rois[i], rois[i + 1])
    obj.setName("Layer " + (i + 1))

    double thicknessPx = layerWidthsPx[i]
    double thicknessMM = thicknessPx * pixelSizeMicrons / 1000.0

    obj.getMeasurementList().put("Layer thickness mm", thicknessMM)

    addObject(obj)
}

// inner layer
def inner = PathObjects.createAnnotationObject(rois[-1])
inner.setName("Layer " + rois.size())

double remainingPx = currentHalf
double remainingMM = remainingPx * pixelSizeMicrons / 1000.0

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
static def getRectFromPoints(List<Point2D.Double> pts) {
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
        longHalf : extentU / 2.0, // half of fitted long side
        shortHalf: extentV / 2.0  // half of fitted short side
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
