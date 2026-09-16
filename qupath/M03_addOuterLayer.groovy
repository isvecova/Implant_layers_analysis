import qupath.lib.gui.scripting.QPEx
import qupath.lib.roi.ROIs
import qupath.lib.roi.GeometryTools
import qupath.lib.objects.PathObjects

import java.awt.geom.Point2D

/*
Script for adding an outer layer annotation around an implant rectangle.

Input:
- Select the existing Layer 1 annotation before running the script.
- The selected annotation is expected to be a layer ROI made from two parallel
  rotated rectangles, representing the current outermost implant layer.

Output:
- An intermediate rectangle annotation named "Convex hull of implant", reconstructed
  from the extreme points of the selected annotation.
- A new ring annotation named "Outer layer", created by expanding the implant
  rectangle by outerLayerThicknessMM along dimensionToSegment and subtracting the
  intermediate implant rectangle.

Note:
- The intermediate rectangle annotation is intentionally added to the image so it
  can be inspected after running the script.
*/

// IMPORTANT: Select the layer 1 annotation before running the script

// Implant parameters
String dimensionToSegment = "shorter" // "shorter" or "longer" dimension of the rectangle to use for layer thickness calculation and layer creation
double outerLayerThicknessMM = 0.5 // thickness of the outer layer in mm

// Get the selected object, exit if no object is selected
def selected = getSelectedObject()

if (selected == null || selected.getROI() == null) {
    print "No ROI selected.\n"
    return
}

// The ROI most likely consists of two parallel rotated rectangles, corresponding to the outermost layer of the implant

// Extract points from the ROI and check if they are valid
def roi = selected.getROI()
def points = roi.getAllPoints()

if (points == null || points.isEmpty()) {
    print "ROI has no points.\n"
    return
}

// Extract the minX, maxX, minY, maxY values - these correspond to the extreme points of the ROI in the X and Y directions
// We also get indices of the extreme values and from these indices, we get the corresponding coordinate (yLow, yHigh, xLow, xHigh)
def ptsX = points.collect { p -> p.getX() }

def minX = Collections.min(ptsX)
def maxX = Collections.max(ptsX)

def minIndexX = ptsX.findIndexOf { it == minX }
def maxIndexX = ptsX.findIndexOf { it == maxX }

def ptsY = points.collect { p -> p.getY() }

def minY = Collections.min(ptsY)
def maxY = Collections.max(ptsY)

def minIndexY = ptsY.findIndexOf { it == minY }
def maxIndexY = ptsY.findIndexOf { it == maxY }

def yLow = ptsY[minIndexX]
def yHigh = ptsY[maxIndexX]

def xLow = ptsX[minIndexY]
def xHigh = ptsX[maxIndexY]

// Get the image data
def imageData = QPEx.getCurrentImageData()
def server = imageData.getServer()
def plane = roi.getImagePlane()

dimensionToSegment = dimensionToSegment?.trim()?.toLowerCase()
if (!(dimensionToSegment in ["shorter", "longer"])) {
    print "Invalid dimensionToSegment: " + dimensionToSegment + ". Use 'shorter' or 'longer'.\n"
    return
}
boolean segmentLongDimension = dimensionToSegment == "longer"

def cal = server.getPixelCalibration()
if (cal == null || !cal.hasPixelSizeMicrons()) {
    print "Pixel calibration missing.\n"
    return
}
double pixelSizeMicrons = cal.getAveragedPixelSizeMicrons()

// Given extreme X and Y values and their corresponding other coordinates, get the four corners

// =============== Define functions needed to create the rectangles =============== 

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

// Function to define the wrapping rotated rectangle of the selected outermost implant layer
static List<Point2D.Double> rectCornersFromCoordinates(minX, maxX, minY, maxY, xLow, xHigh, yLow, yHigh) {

    // Return corners in order (clockwise or counterclockwise)
    [
        // Bottom-left
        new Point2D.Double(minX, yLow),
        // Bottom-right
        new Point2D.Double(xHigh, maxY),
        // Top-right
        new Point2D.Double(maxX, yHigh),
        // Top-left
        new Point2D.Double(xLow, minY)
    ]
}


// Get corners of a rectangle given its center, orientation (ux, uy), and half-lengths along the long and short axes
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

// Create a polygon ROI from a list of points
def makePolyROI = { polyPts ->
    double[] xs = new double[polyPts.size()]
    double[] ys = new double[polyPts.size()]
    for (int i = 0; i < polyPts.size(); i++) {
        xs[i] = polyPts[i].x
        ys[i] = polyPts[i].y
    }
    ROIs.createPolygonROI(xs, ys, plane)
}

// Create a new ROI that represents the outer layer by subtracting the inner ROI from the outer ROI
def makeLayerROI = { outerROI, innerROI ->
    def gOuter = GeometryTools.roiToGeometry(outerROI)
    def gInner = GeometryTools.roiToGeometry(innerROI)
    def gRing = gOuter.difference(gInner)

    PathObjects.createAnnotationObject(
        GeometryTools.geometryToROI(gRing, plane)
    )
}

// ==============================================================


// Reconstruct the rectangle outlining the implant from the extreme points of the user-selected outermost layer
def corners = rectCornersFromCoordinates(minX, maxX, minY, maxY, xLow, xHigh, yLow, yHigh)
def implantROI = makePolyROI(corners)
def inner = PathObjects.createAnnotationObject(implantROI)
inner.setName("Convex hull of implant")
addObject(inner)

// Extract information about the rectangle
def implantPoints = implantROI.getAllPoints()
def pts = implantPoints.collect { p -> new Point2D.Double(p.getX(), p.getY()) }
def rectInfo = getRectFromPoints(pts)


// Create outer layer ROI of a given thickness
outerRectangle = [
    center   : rectInfo.center,
    ux       : rectInfo.ux,
    uy       : rectInfo.uy,
    longHalf : segmentLongDimension ? (rectInfo.longHalf + outerLayerThicknessMM * 1000.0 / pixelSizeMicrons) : rectInfo.longHalf,
    shortHalf: segmentLongDimension ? rectInfo.shortHalf : (rectInfo.shortHalf + outerLayerThicknessMM * 1000.0 / pixelSizeMicrons)
]
def outerCorners = rectCorners(outerRectangle.center, outerRectangle.ux, outerRectangle.uy, outerRectangle.longHalf, outerRectangle.shortHalf)
def outerROI = makePolyROI(outerCorners)

// Subtract the implant ROI from the extended ROI to create the outer layer ROI
def outer = makeLayerROI(outerROI, implantROI)
outer.setName("Outer layer")

outer.getMeasurementList().put("Layer thickness mm", outerLayerThicknessMM)

addObject(outer)
