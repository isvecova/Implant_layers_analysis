# coding: utf-8

"""
Build rectangular layer segments for implant ROIs.

This script takes a tagged implant ROI in arivis and creates nested rectangular
segments that represent implant layers.

Processing steps:
1. Read ROI contour points from the source segment.
2. Fit an oriented rectangle to the ROI using PCA.
3. Compute visible layer widths from the cylindrical implant model.
4. Generate outer-to-inner rectangular layer segments.
5. Save output segments with tags and custom features.

Main parts of the script:
- Logger: Handles console and optional file logging.
- ProcessingConfig: Stores and validates user parameters.
- Geometry helpers: Rectangle fitting and layer-width calculations.
- Main pipeline: arivis I/O, segment creation, and feature assignment.

Author: Iva Svecova (svecovaiva01@gmail.com)
Last edited: 2026-04-18
Created: 2026-04-06

This script was adapted from the Matryoshka script by Maurizio Abbate: 
# NAME: Matrioshka Doll 
# FILE: MATRYOSHKA_DOLL_revE1(4_4)_OP
# REVISION : 1.20  - 2023-12-15
# AUTHOR : Maurizio Abbate
# Copyright(c) 2025 Carl Zeiss Microscopy Software Center Rostock GmbH,Germany. 
#              All Rights Reserved.
"""

import math
import time
import os
import json
from datetime import datetime
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional

import arivis
import arivis_core as core
import arivis_objects as objects


# ============================================================================
# Configuration & Logging
# ============================================================================

@dataclass
class ProcessingConfig:
    """Configuration for implant layer processing.
    
    Attributes:
        input_tag: Tag to find source segments in the Arivis store.
        n_layers: Nominal number of radial layers in the implant model.
            Used for layer width calculation, not the actual count of generated
            layers (which may be fewer if the cross-section is smaller).
        radial_layer_thickness_mm: Radial thickness of each model layer in mm.
            Used in layer width calculation based on cylindrical geometry.
        implant_height_mm: Nominal implant height in mm.
            Used to set the long axis of output rectangles; not measured from
            the input ROI.
        verbose_path: Directory path for verbose log output.
            Empty string disables file logging. Use os.path.expanduser() for
            paths with ~ on Unix-like systems.
    """
    input_tag: str
    n_layers: int
    radial_layer_thickness_mm: float
    implant_height_mm: float
    verbose_path: str

    def __post_init__(self) -> None:
        """Validate configuration parameters after initialization."""
        if self.n_layers <= 0:
            raise ValueError(f"n_layers must be positive, got {self.n_layers}")
        if self.radial_layer_thickness_mm <= 0:
            raise ValueError(
                f"radial_layer_thickness_mm must be positive, "
                f"got {self.radial_layer_thickness_mm}"
            )
        if self.implant_height_mm <= 0:
            raise ValueError(
                f"implant_height_mm must be positive, "
                f"got {self.implant_height_mm}"
            )


class Logger:
    """Simple logger that writes to console and optional file.
    
    Encapsulates logging behavior to avoid threading parameters through
    function calls. Supports both console output and file appending.
    
    Attributes:
        module_name: Prefix for log messages (typically the calling module name).
        verbose_path: Directory for log file output. If empty/None, file logging
            is disabled.
        verbose_filename: Name of log file. If not provided, a filename is
            generated from script name and run timestamp.
    """

    def __init__(
        self,
        module_name: str,
        verbose_path: Optional[str] = None,
        verbose_filename: Optional[str] = None
    ) -> None:
        """Initialize logger.
        
        Args:
            module_name: Prefix for all log messages.
            verbose_path: Directory for log file. If empty or None, file logging
                is disabled.
            verbose_filename: Name of log file in verbose_path. If None, uses
                '<script_name>_<YYYYMMDD_HHMMSS>.log'.
        """
        self.module_name = module_name
        self.verbose_path = verbose_path
        self.run_datetime = datetime.today()
        self.script_name = os.path.basename(__file__)

        if verbose_filename:
            self.verbose_filename = verbose_filename
        else:
            script_base, _ = os.path.splitext(self.script_name)
            timestamp = self.run_datetime.strftime("%Y%m%d_%H%M%S")
            self.verbose_filename = f"{script_base}_{timestamp}.log"
        self.enabled = bool(verbose_path)

        if self.enabled:
            self._initialize_log_file()

    def _initialize_log_file(self) -> None:
        """Create or reset the verbose log file with header."""
        if not os.path.isdir(self.verbose_path):
            self.log(f"Warning: verbose path does not exist: {self.verbose_path}")
            return

        fullpath = os.path.join(self.verbose_path, self.verbose_filename)
        try:
            with open(fullpath, "w") as f:
                f.write(f"Run datetime: {self.run_datetime.isoformat()}\n")
                f.write(f"Script: {self.script_name}\n")
                f.write(f"Module: {self.module_name}\n")
        except IOError as e:
            self.log(f"Warning: could not initialize log file: {e}")

    def log(self, message: str) -> None:
        """Write message to console and optionally to file.
        
        Args:
            message: The message to log.
        """
        formatted = f"[{self.module_name}] {message}"
        print(formatted)

        if not self.enabled or not os.path.isdir(self.verbose_path):
            return

        fullpath = os.path.join(self.verbose_path, self.verbose_filename)
        try:
            with open(fullpath, "a") as f:
                f.write(f"{formatted}\n")
        except IOError as e:
            print(f"Warning: could not write to log file: {e}")


# ============================================================================
# Geometry Helpers
# ============================================================================

def get_polygon_points_from_object(obj: objects.Segment) -> Tuple[Optional[object], Optional[List[Tuple[float, float]]]]:
    """Extract 2D contour points from the first polygon plane of an object.
    
    Assumes the input object is effectively a 2D polygon in a single plane.
    
    Args:
        obj: An Arivis Segment object.
    
    Returns:
        Tuple of (plane, points_list) where:
        - plane: The Arivis plane object containing the polygon, or None on error.
        - points_list: List of (x, y) tuples representing the polygon contour,
          or None if no valid contour found.
    """
    planes = obj.get_all_planes()
    if planes is None or len(planes) == 0:
        return None, None

    plane = planes[0]
    polygons = obj.get_polygons(plane)
    if polygons is None or len(polygons) == 0:
        return plane, None

    poly = polygons[0]
    contour = poly.get_contour()
    if contour is None or len(contour) == 0:
        return plane, None

    pts = [(float(p.x), float(p.y)) for p in contour]
    return plane, pts


def get_rect_from_points(pts: List[Tuple[float, float]]) -> Dict[str, float]:
    """Fit an oriented rectangle to 2D points using PCA.
    
    Uses Principal Component Analysis to find the principal axes of the
    points, then computes the bounding rectangle aligned with these axes.
    
    Algorithm:
    1. Compute covariance matrix from points (centered at centroid).
    2. Find principal axes via eigendecomposition (using atan2 for angle).
    3. Project points onto principal axes.
    4. Compute axis-aligned bounding box in the rotated space.
    5. Ensure longer axis is 'u' (longHalf >= shortHalf).
    
    Args:
        pts: List of (x, y) tuples representing polygon vertices.
    
    Returns:
        Dictionary with keys:
        - 'center': Tuple (cx, cy) of rectangle center.
        - 'ux', 'uy': Unit vector for long axis (length 1).
        - 'longHalf': Half-length of long axis (pixels).
        - 'shortHalf': Half-length of short axis (pixels).
    """
    # Compute centroid
    mx = sum(p[0] for p in pts) / float(len(pts))
    my = sum(p[1] for p in pts) / float(len(pts))

    # Compute covariance matrix elements
    sxx = sxy = syy = 0.0
    for p in pts:
        dx = p[0] - mx
        dy = p[1] - my
        sxx += dx * dx
        syy += dy * dy
        sxy += dx * dy

    # PCA: compute principal angles using eigenvalue decomposition
    # For 2x2 covariance, eigenvalues/eigenvectors can be found analytically.
    # Angle from atan2 gives the principal direction.
    theta = 0.5 * math.atan2(2.0 * sxy, (sxx - syy))

    ux = math.cos(theta)
    uy = math.sin(theta)
    vx = -uy  # Perpendicular axis
    vy = ux

    # Project points onto u,v axes and find extents
    min_u = min_v = 1e12
    max_u = max_v = -1e12

    for p in pts:
        dx = p[0] - mx
        dy = p[1] - my
        u = dx * ux + dy * uy
        v = dx * vx + dy * vy

        min_u, max_u = min(min_u, u), max(max_u, u)
        min_v, max_v = min(min_v, v), max(max_v, v)

    extent_u = max_u - min_u
    extent_v = max_v - min_v

    # Ensure u is the long axis
    if extent_v > extent_u:
        ux, vx, uy, vy = vx, ux, vy, uy
        extent_u, extent_v = extent_v, extent_u

    # Rectangle center in original space
    center_u = 0.5 * (min_u + max_u)
    center_v = 0.5 * (min_v + max_v)
    cx = mx + center_u * ux + center_v * vx
    cy = my + center_u * uy + center_v * vy

    return {
        "center": (cx, cy),
        "ux": ux,
        "uy": uy,
        "longHalf": 0.5 * extent_u,
        "shortHalf": 0.5 * extent_v
    }


def rect_corners(
    center: Tuple[float, float],
    ux: float,
    uy: float,
    long_half: float,
    short_half: float
) -> List[Tuple[float, float]]:
    """Compute the four corner coordinates of an oriented rectangle.
    
    Args:
        center: (cx, cy) center of rectangle.
        ux, uy: Unit vector for long axis.
        long_half: Half-length of long axis.
        short_half: Half-length of short axis.
    
    Returns:
        List of four (x, y) tuples in counterclockwise order starting from
        the bottom-left (-long, -short) corner.
    """
    cx, cy = center
    vx = -uy  # Perpendicular to (ux, uy)
    vy = ux

    return [
        (cx - long_half * ux - short_half * vx, cy - long_half * uy - short_half * vy),
        (cx + long_half * ux - short_half * vx, cy + long_half * uy - short_half * vy),
        (cx + long_half * ux + short_half * vx, cy + long_half * uy + short_half * vy),
        (cx - long_half * ux + short_half * vx, cy - long_half * uy + short_half * vy)
    ]


def make_polygon_from_points(points_xy: List[Tuple[float, float]]) -> objects.Polygon:
    """Create an Arivis Polygon from a list of (x, y) coordinates.
    
    Args:
        points_xy: List of (x, y) tuples.
    
    Returns:
        An Arivis Polygon object with the contour set.
    """
    poly = objects.Polygon()
    contour = [core.Point2D(float(x), float(y)) for x, y in points_xy]
    poly.set_contour(contour)
    return poly


def create_segment_from_outer_inner_rectangles(
    template_object: objects.Segment,
    plane: object,
    outer_pts: List[Tuple[float, float]],
    inner_pts: Optional[List[Tuple[float, float]]] = None
) -> objects.Segment:
    """Create an Arivis Segment representing a rectangular layer.
    
    If inner_pts is provided, creates a hollow rectangle (annulus) with
    the inner rectangle as a hole. Otherwise creates a solid rectangle.
    
    Args:
        template_object: Source segment to copy timepoint and color from.
        plane: Arivis plane object for the segment.
        outer_pts: List of (x, y) tuples for outer rectangle corners.
        inner_pts: List of (x, y) tuples for inner rectangle (hole), or None
            for a solid rectangle.
    
    Returns:
        An Arivis Segment object with the polygon(s) added.
    """
    seg = objects.Segment()
    seg.set_timepoint(template_object.get_timepoint())
    seg.set_color(template_object.get_color())

    poly = make_polygon_from_points(outer_pts)
    if inner_pts is not None:
        hole = [core.Point2D(float(x), float(y)) for x, y in inner_pts]
        poly.set_holes([hole])

    seg.add_polygon(poly, plane=plane)
    return seg


def compute_layer_widths_px(
    short_half_px: float,
    pixel_size_xy_um: float,
    n_layers: int,
    radial_layer_thickness_mm: float,
    logger: Logger
) -> Optional[List[float]]:
    """Compute visible layer widths for a cylindrical implant model.

    Given a rectangular cross-section with width w (from shortHalf), this
    computes how much of each theoretical radial layer is visible in that
    cross-section.
    
    Geometry Model:
    - Implant is a cylinder with nominal radius = n_layers * radial_layer_thickness_mm.
    - Each concentric layer has thickness = radial_layer_thickness_mm.
    - The cross-section is a rectangle with width w = 2 * shortHalf.
    - The offset k = sqrt(R^2 - (w/2)^2) / t determines which layers are visible.
    - Only layers l where l > k contribute to the visible section.
    
    Args:
        short_half_px: Half-width of rectangular cross-section (pixels).
        pixel_size_xy_um: Pixel size in micrometers.
        n_layers: Nominal number of radial layers.
        radial_layer_thickness_mm: Thickness of each layer in mm.
        logger: Logger instance for messages.
    
    Returns:
        List of visible layer widths in pixels, ordered from outermost to
        innermost, or None if no layers are visible.
    """
    width_mm = (short_half_px * 2.0) * pixel_size_xy_um / 1000.0
    logger.log(f"[compute_layer_widths_px] Section width [mm]: {width_mm}")

    k = 0.0
    max_radius_mm = n_layers * radial_layer_thickness_mm

    if width_mm < 2.0 * max_radius_mm:
        inside_sqrt = (max_radius_mm ** 2) - ((width_mm / 2.0) ** 2)
        if inside_sqrt < 0:
            inside_sqrt = 0
        k = math.sqrt(inside_sqrt) / radial_layer_thickness_mm
    else:
        logger.log(
            "[compute_layer_widths_px] Warning: width exceeds nominal implant "
            "diameter for supplied n_layers and thickness."
        )

    layer_widths_px: List[float] = []
    inside = 0.0

    for l in range(1, n_layers + 1):
        if l <= k:
            continue

        inside_sqrt = ((l * radial_layer_thickness_mm) ** 2) - \
                      ((k * radial_layer_thickness_mm) ** 2)
        if inside_sqrt < 0:
            inside_sqrt = 0

        x = math.sqrt(inside_sqrt)
        layer_width_mm = (x - inside)
        layer_width_px = layer_width_mm * 1000.0 / pixel_size_xy_um

        layer_widths_px.append(layer_width_px)
        inside = x

    layer_widths_px.reverse()
    return layer_widths_px if layer_widths_px else None


# ============================================================================
# Arivis Integration
# ============================================================================

def get_object_type(
    store: objects.Store,
    id_obj: int,
    outline: bool = True,
    str_tag: str = "",
    logger: Logger = None
) -> Tuple[str, Optional[objects.Segment]]:
    """Retrieve an object from the store and determine its type.
    
    Args:
        store: Arivis Store to query.
        id_obj: Object ID to retrieve.
        outline: Whether to get outline (True) or filled object (False).
        str_tag: If non-empty, verify the object has this tag.
        logger: Logger instance for error messages.
    
    Returns:
        Tuple of (object_type, object) where:
        - object_type: String type (e.g., "Segment"), or empty string on error.
        - object: The Arivis object, or None on error.
    """
    seg_type = ""
    if store is None:
        if logger:
            logger.log("[get_object_type] No store available")
        return seg_type, None

    obj1 = store.get_object(id_obj, outline)
    if obj1 is None:
        if logger:
            logger.log("[get_object_type] Error: invalid object instance")
        return seg_type, None

    tags = obj1.get_tags()
    if str_tag != "" and str_tag not in tags:
        if logger:
            logger.log(f"[get_object_type] Tag '{str_tag}' not on object")
        return seg_type, None

    feature_desc = store.get_feature_descriptor("Type")
    if feature_desc is None:
        if logger:
            logger.log("[get_object_type] Feature descriptor 'Type' not found")
        return seg_type, obj1

    feature_value = feature_desc.get_values()
    if feature_value is None or len(feature_value) == 0:
        if logger:
            logger.log("[get_object_type] No values in feature descriptor 'Type'")
        return seg_type, obj1

    feature_value_1 = feature_value[0].get_name()
    feature_1 = store.get_feature_for_object(feature_desc, obj1)
    if feature_1 is None:
        if logger:
            logger.log("[get_object_type] Feature 'Type' missing for object")
        return seg_type, obj1

    seg_type = feature_1.get_value(feature_value_1)
    return seg_type, obj1


def create_custom_feature(
    store: objects.Store,
    description: str,
    name: str,
    feature_type: int = objects.ValueDescriptor.TYPE_FLOAT,
    logger: Logger = None
) -> bool:
    """Create a stored feature descriptor if it does not already exist.
    
    Args:
        store: Arivis Store to add feature to.
        description: Long name/description of the feature.
        name: Short name for the feature value.
        feature_type: Arivis ValueDescriptor type constant (TYPE_FLOAT, TYPE_STRING, etc.).
        logger: Logger instance for messages.
    
    Returns:
        True if feature was created or already exists, False on error.
    """
    if store is None:
        if logger:
            logger.log("[create_custom_feature] No store available")
        return False

    list_descriptors = store.get_feature_descriptors()
    for descriptor in list_descriptors:
        if descriptor.get_name() == description:
            if logger:
                logger.log(f"[create_custom_feature] Feature already exists: {description}")
            return True

    feature_descriptor = objects.FeatureDescriptor()
    feature_descriptor.set_name(description)

    val_desc = objects.ValueDescriptor()
    val_desc.set_name(name)
    val_desc.set_type(feature_type)
    feature_descriptor.add_value(val_desc)

    if not store.create_stored_feature(feature_descriptor):
        if logger:
            logger.log(f"[create_custom_feature] Failed to create: {description}")
    else:
        if logger:
            logger.log(f"[create_custom_feature] Created: {description}")

    return True


def insert_custom_feature(
    store: objects.Store,
    value: object,
    obj1: objects.Segment,
    description: str,
    name: str,
    logger: Logger = None
) -> bool:
    """Insert a feature value for an object.
    
    Args:
        store: Arivis Store.
        value: Value to set (will be converted as needed).
        obj1: Arivis object to set feature on.
        description: Long feature name.
        name: Short feature name.
        logger: Logger instance for error messages.
    
    Returns:
        True on success, False on error.
    """
    if store is None or obj1 is None:
        if logger:
            logger.log("[insert_custom_feature] Error: store/object not available")
        return False

    feature_desc = store.get_feature_descriptor(description)
    if feature_desc is None:
        if logger:
            logger.log(f"[insert_custom_feature] Feature not found: {description}")
        return False

    feature = objects.Feature()
    feature.set_value(name, value)
    store.set_feature_for_object(feature, feature_desc, obj1)
    return True


# ============================================================================
# Main Processing
# ============================================================================

def build_layers_for_object(
    store_out: objects.Store,
    obj1: objects.Segment,
    source_id: int,
    config: ProcessingConfig,
    pixel_size_xy_um: float,
    logger: Logger
) -> bool:
    """Generate concentric layer segments for a single ROI object.
    
    Fits an oriented rectangle to the input object's polygon, then generates
    concentric rectangular layer segments based on the cylindrical implant model.
    
    Args:
        store_out: Arivis Store to add output segments to.
        obj1: Source Arivis Segment containing the ROI polygon.
        source_id: Object ID to record in feature metadata.
        config: Processing configuration.
        pixel_size_xy_um: Pixel size in micrometers.
        logger: Logger instance.
    
    Returns:
        True on success, False on error.
    """
    plane, pts = get_polygon_points_from_object(obj1)
    if pts is None or len(pts) < 4:
        logger.log("[build_layers_for_object] Input has insufficient polygon points")
        return False

    # Fit oriented rectangle to contour points
    rect_info = get_rect_from_points(pts)
    logger.log(f"Fitted rectangle shortHalf [px]: {rect_info['shortHalf']:.2f}")
    logger.log(f"Fitted rectangle longHalf [px]: {rect_info['longHalf']:.2f}")

    # Use user-specified height instead of fitted rectangle's height
    long_half_px = config.implant_height_mm * 1000.0 / (2.0 * pixel_size_xy_um)
    measured_height_mm = (2.0 * rect_info["longHalf"]) * pixel_size_xy_um / 1000.0
    logger.log(
        f"Measured ROI height: {measured_height_mm:.2f} mm; "
        f"using specified height: {config.implant_height_mm:.2f} mm"
    )

    # Compute visible layer widths
    layer_widths_px = compute_layer_widths_px(
        short_half_px=rect_info["shortHalf"],
        pixel_size_xy_um=pixel_size_xy_um,
        n_layers=config.n_layers,
        radial_layer_thickness_mm=config.radial_layer_thickness_mm,
        logger=logger
    )

    if layer_widths_px is None or len(layer_widths_px) == 0:
        logger.log("[build_layers_for_object] No visible layers computed")
        return False

    # Build successive rectangles from outer to inner
    rectangles: List[Dict[str, float]] = []
    current_half = rect_info["shortHalf"]

    # If the initial half-width exceeds the maximum based on the model, cap it
    max_start_half = (config.n_layers * config.radial_layer_thickness_mm * 1000.0 / pixel_size_xy_um)
    if current_half > max_start_half:
        logger.log(
            f"Capping initial half-width from {current_half:.2f} px "
            f"to {max_start_half:.2f} px (n_layers * user layer thickness)"
        )
        current_half = max_start_half

    for i in range(0, len(layer_widths_px) - 1):
        layer_width = layer_widths_px[i]
        if layer_width <= 0:
            continue

        rectangles.append({
            "center": rect_info["center"],
            "ux": rect_info["ux"],
            "uy": rect_info["uy"],
            "longHalf": long_half_px,
            "shortHalf": current_half
        })

        next_half = current_half - layer_width
        if next_half <= 0:
            current_half = 0
            break

        current_half = next_half

    # Add innermost rectangle
    rectangles.append({
        "center": rect_info["center"],
        "ux": rect_info["ux"],
        "uy": rect_info["uy"],
        "longHalf": long_half_px,
        "shortHalf": current_half
    })

    # Convert rectangles to corner points
    rois_pts = [
        rect_corners(r["center"], r["ux"], r["uy"], r["longHalf"], r["shortHalf"])
        for r in rectangles
    ]

    # Create ring layers (annuli)
    created_count = 0
    for i in range(0, len(rois_pts) - 1):
        outer_pts = rois_pts[i]
        inner_pts = rois_pts[i + 1]
        seg = create_segment_from_outer_inner_rectangles(obj1, plane, outer_pts, inner_pts)
        seg.add_tag(f"Layer_{i + 1}")
        seg.add_tag("_LAYER_MODEL_")
        store_out.add_object(seg)

        thickness_px = layer_widths_px[i]
        thickness_mm = thickness_px * pixel_size_xy_um / 1000.0

        insert_custom_feature(store_out, str(i + 1), seg, "Layer_Index", "L_IDX", logger)
        insert_custom_feature(store_out, float(thickness_mm), seg, "Layer_Thickness_mm", "L_TH_MM", logger)
        insert_custom_feature(store_out, str(source_id), seg, "Source_Object_ID", "SRC_ID", logger)

        created_count += 1

    # Add innermost solid layer
    inner_pts = rois_pts[-1]
    inner_seg = create_segment_from_outer_inner_rectangles(obj1, plane, inner_pts, None)
    inner_seg.add_tag(f"Layer_{len(rois_pts)}")
    inner_seg.add_tag("_LAYER_MODEL_")
    store_out.add_object(inner_seg)

    remaining_mm = current_half * pixel_size_xy_um / 1000.0
    insert_custom_feature(store_out, str(len(rois_pts)), inner_seg, "Layer_Index", "L_IDX", logger)
    insert_custom_feature(store_out, float(remaining_mm), inner_seg, "Layer_Thickness_mm", "L_TH_MM", logger)
    insert_custom_feature(store_out, str(source_id), inner_seg, "Source_Object_ID", "SRC_ID", logger)

    created_count += 1

    logger.log(f"Created {created_count} layer segments for object {source_id}")
    return True


def main(config: ProcessingConfig) -> None:
    """Main processing entry point for generating implant layer segments.
    
    Orchestrates the full pipeline:
    1. Initialize Arivis context (viewer, document, image set, stores).
    2. Extract pixel calibration.
    3. Create custom feature descriptors.
    4. Process each tagged object to generate layer segments.
    
    Args:
        config: ProcessingConfig object with all parameters.
    """
    start_time = time.time()
    logger = Logger("RECT_LAYER_MODEL", config.verbose_path)

    logger.log("Starting implant layer generation...")

    # Validate inputs
    try:
        config.__post_init__()  # Runs validation
    except ValueError as e:
        logger.log(f"Configuration error: {e}")
        return

    # Get active Arivis context
    viewer = arivis.App.get_active_viewer()
    if not viewer:
        logger.log("Error: No active viewer found")
        return

    doc = viewer.get_document()
    if not doc:
        logger.log("Error: No document in viewer")
        return

    image_set = doc.get_default_imageset()
    if not image_set:
        logger.log("Error: No imageset in document")
        return

    # Find object store
    store_in = None
    store_out = None
    id_list: List[int] = []
    for store_kind in (objects.Store.DOCUMENT_STORE, objects.Store.ANALYSIS_STORE, None):
        candidate = doc.get_store(image_set, store_kind)
        if not candidate:
            continue
        ids = candidate.get_object_ids(config.input_tag)
        if ids:
            store_in = candidate
            store_out = store_out or store_in
            id_list = ids
            logger.log(f"Found {len(id_list)} objects with tag '{config.input_tag}'")
            break

    if not store_in or not store_out or not id_list:
        logger.log("Error: No store or objects found for tag")
        return

    # Get pixel calibration
    try:
        pixel_size = image_set.get_pixel_size()
        pixel_size_x_um = float(pixel_size[0])
        pixel_size_y_um = float(pixel_size[1])
        # Average pixel size (assumes approximately isotropic imaging)
        pixel_size_xy_um = 0.5 * (pixel_size_x_um + pixel_size_y_um)
    except Exception as e:
        logger.log(f"Error reading pixel size: {e}")
        return

    logger.log(f"Pixel size: {pixel_size_xy_um:.4f} um")

    # Create custom features
    create_custom_feature(store_out, "Layer_Index", "L_IDX", objects.ValueDescriptor.TYPE_STRING, logger)
    create_custom_feature(store_out, "Layer_Thickness_mm", "L_TH_MM", objects.ValueDescriptor.TYPE_FLOAT, logger)
    create_custom_feature(store_out, "Source_Object_ID", "SRC_ID", objects.ValueDescriptor.TYPE_STRING, logger)

    # Process each object
    for obj_id in id_list:
        logger.log(f"Processing object ID: {obj_id}")

        seg_type, obj = get_object_type(
            store_in, obj_id, outline=True,
            str_tag=config.input_tag, logger=logger
        )
        if obj is None:
            logger.log(f"Skipped: Invalid object {obj_id}")
            continue

        if seg_type != "Segment":
            logger.log(f"Skipped: Object {obj_id} is type '{seg_type}', not Segment")
            continue

        try:
            build_layers_for_object(
                store_out, obj, obj_id,
                config, pixel_size_xy_um, logger
            )
        except Exception as e:
            logger.log(f"Error processing object {obj_id}: {e}")

    elapsed = time.time() - start_time
    logger.log(f"Completed in {elapsed:.2f} seconds")


# ============================================================================
# Script Execution
# ============================================================================

def load_config_from_file(filepath: str) -> Optional[ProcessingConfig]:
    """Load ProcessingConfig from a JSON file.
    
    Args:
        filepath: Path to the JSON configuration file.
    
    Returns:
        ProcessingConfig on success, None on error.
    """
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
        
        config = ProcessingConfig(
            input_tag=data.get("input_tag", "implant"),
            n_layers=data.get("n_layers", 5),
            radial_layer_thickness_mm=data.get("radial_layer_thickness_mm", 0.4),
            implant_height_mm=data.get("implant_height_mm", 7.5),
            verbose_path=data.get("verbose_path", "")
        )
        print(f"[CONFIG] Loaded configuration from: {filepath}")
        return config
    except FileNotFoundError:
        return None
    except json.JSONDecodeError as e:
        print(f"[CONFIG] Error parsing JSON: {e}")
        return None
    except Exception as e:
        print(f"[CONFIG] Error loading config file: {e}")
        return None


def load_config() -> ProcessingConfig:
    """Load configuration with two-tier fallback:
    1. Default config file in script directory (config_default.json)
    2. Hardcoded defaults (if file doesn't exist)
    
    Returns:
        ProcessingConfig object.
    """
    # Tier 1: Default config file in script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_config_path = os.path.join(script_dir, "config_default.json")
    
    if os.path.exists(default_config_path):
        config = load_config_from_file(default_config_path)
        if config:
            return config
    
    # Tier 2: Hardcoded defaults
    print("[CONFIG] Using hardcoded default configuration")
    config = ProcessingConfig(
        input_tag="implant",
        n_layers=5,
        radial_layer_thickness_mm=0.4,
        implant_height_mm=7.5,
        verbose_path="D:/logs"
    )
    return config


if __name__ == "__main__":
    config = load_config()
    main(config)
