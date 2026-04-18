import qupath.lib.gui.scripting.QPEx
import qupath.lib.roi.ROIs
import qupath.lib.objects.PathObjects
import qupath.lib.regions.ImagePlane

import java.awt.geom.Point2D

/*
Script for fitting a rotated rectangle to a polygonal annotation 
and creating a new annotation with the specified long side length.

Input: 
- A selected polygonal annotation object representing the rectangular section.  

Output: 
- An annotation object: 
    - The rotation is determined by PCA of the points in the original polygon.
    - The long side is fixed to the specified height value (default 7.5 mm).
    - The shord side is set to bind the original polygon.

Written by: Iva Svecova (isvecova47@gmail.com)
Date: 2026-04-05
*/

// ------------------------------------------------------------
// 0) Specify parameters
// ------------------------------------------------------------
// Implant parameter
double implantHeightMM = 7.5  // mm

// ------------------------------------------------------------
// 1) get selected ROI
// ------------------------------------------------------------
def selected = getSelectedObject()

if (selected == null || selected.getROI() == null) {
    print "No ROI selected.\n"
    return
}

def roi = selected.getROI()

// Ensure polygon-like ROI
def points = roi.getAllPoints()

if (points == null || points.isEmpty()) {
    print "ROI has no points.\n"
    return
}

def pts = points.collect { p -> new Point2D.Double(p.getX(), p.getY()) }

// ------------------------------------------------------------
// 2) calibration → convert given height to pixels
// ------------------------------------------------------------
def server = getCurrentImageData().getServer()
def cal = server.getPixelCalibration()

if (cal == null || !cal.hasPixelSizeMicrons()) {
    print "Pixel calibration missing.\n"
    return
}

double pixelSizeMicrons = cal.getAveragedPixelSizeMicrons()
double fixedLongPx = implantHeightMM * 1000 / pixelSizeMicrons   // full resolution

println "Fixed long side [px]: ${fixedLongPx}"

// ------------------------------------------------------------
// 3) fit rotated rectangle (PCA-based)
// ------------------------------------------------------------
def rectInfo = fitRectWithFixedLongSide(pts, fixedLongPx)

// ------------------------------------------------------------
// 4) build rectangle ROI
// ------------------------------------------------------------
def corners = rectCorners(
    rectInfo.center,
    rectInfo.ux,
    rectInfo.uy,
    rectInfo.longHalf,
    rectInfo.shortHalf
)

// Use same plane as original ROI
def plane = roi.getImagePlane()

double[] xs = new double[corners.size()]
double[] ys = new double[corners.size()]

for (int i = 0; i < corners.size(); i++) {
    xs[i] = corners[i].x
    ys[i] = corners[i].y
}

def rectROI = ROIs.createPolygonROI(xs, ys, plane)

// ------------------------------------------------------------
// 5) create annotation
// ------------------------------------------------------------
def rectObj = PathObjects.createAnnotationObject(rectROI)
rectObj.setName("Rotated rectangle (" + String.format("%.1f mm", implantHeightMM) + " long)")

addObject(rectObj)
selectObjects(rectObj)

print "Done.\n"


// ============================================================
// Helper functions
// ============================================================

// Fit a rectangle with a fixed long side to a set of 2D points using PCA
static def fitRectWithFixedLongSide(List<Point2D.Double> pts, double fixedLongSidePx) {
    // 1. Compute centroid (mean x, y)
    double mx = pts.collect{it.x}.sum() / pts.size()
    double my = pts.collect{it.y}.sum() / pts.size()

    // 2. Compute covariance components (for PCA)
    double sxx = 0, syy = 0, sxy = 0
    for (p in pts) {
        double dx = p.x - mx
        double dy = p.y - my
        sxx += dx * dx      // variance in x
        syy += dy * dy      // variance in y
        sxy += dx * dy      // covariance
    }

    // 3. Calculate principal axis orientation (theta)
    double theta = 0.5 * Math.atan2(2 * sxy, (sxx - syy))

    // 4. Get unit vectors for principal (u) and orthogonal (v) axes
    double ux = Math.cos(theta)
    double uy = Math.sin(theta)

    double vx = -uy
    double vy = ux

    // 5. Project points onto principal axes to find extents
    double minU = Double.POSITIVE_INFINITY
    double maxU = Double.NEGATIVE_INFINITY
    double minV = Double.POSITIVE_INFINITY
    double maxV = Double.NEGATIVE_INFINITY

    for (p in pts) {
        double dx = p.x - mx
        double dy = p.y - my

        double u = dx * ux + dy * uy // projection onto principal axis
        double v = dx * vx + dy * vy // projection onto orthogonal axis

        if (u < minU) minU = u
        if (u > maxU) maxU = u
        if (v < minV) minV = v
        if (v > maxV) maxV = v
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
    double centerU = (minU + maxU) / 2.0
    double centerV = (minV + maxV) / 2.0

    double cx = mx + centerU * ux + centerV * vx
    double cy = my + centerU * uy + centerV * vy

    // 8. Return rectangle parameters: center, orientation, half-lengths
    [
        center   : new Point2D.Double(cx, cy),
        ux       : ux, // unit vector x-component (long axis)
        uy       : uy, // unit vector y-component (long axis)
        longHalf : fixedLongSidePx / 2.0, // half of fixed long side
        shortHalf: extentV / 2.0          // half of fitted short side
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
        new Point2D.Double(center.x - longHalf * ux - shortHalf * vx, center.y - longHalf * uy - shortHalf * vy),
        // Bottom-right
        new Point2D.Double(center.x + longHalf * ux - shortHalf * vx, center.y + longHalf * uy - shortHalf * vy),
        // Top-right
        new Point2D.Double(center.x + longHalf * ux + shortHalf * vx, center.y + longHalf * uy + shortHalf * vy),
        // Top-left
        new Point2D.Double(center.x - longHalf * ux + shortHalf * vx, center.y - longHalf * uy + shortHalf * vy)
    ]
}