# Implant offset layer creation

## Scripts
This repository provides scripts for creating concentric layers in a rectangular implant section for image analysis in QuPath. There are two alrernative workflows:

### 1. Automatic workflow

- [**01_detectAndMakeLayers.groovy**](01_detectAndMakeLayers.groovy)  
  Automatically detects a dark rectangular implant in the image, fits a rectangle using thresholding and PCA, and generates concentric rectangular layers as annotation objects.  
  - **Input:** Image with a single, dark rectangular implant area.
  - **Output:** Concentric rectangular layers with thickness measurements in mm.

### 2. Manual workflow

If automatic detection is not suitable, use the manual path:

- [**M01_polygonToRectangle.groovy**](M01_polygonToRectangle.groovy)  
  Converts a user-drawn polygonal ROI (around the area of interest) into a fitted rotated rectangle with a fixed long side (implant height).
  - **Input:** Selected polygonal annotation.
  - **Output:** Rotated rectangle annotation matching the implant section.

- [**M02_makeLayers.groovy**](M02_makeLayers.groovy)  
  Takes the fitted rectangle (from M01) and creates concentric rectangular layers, calculating thickness based on a cylindrical model.
  - **Input:** Selected rectangle annotation.
  - **Output:** Concentric rectangular layers with thickness measurements in mm.

#### Manual path steps:
1. Draw a polygon ROI around the implant area.
2. Run **M01_polygonToRectangle.groovy** to fit a rectangle.
3. Run **M02_makeLayers.groovy** to generate layers.

## Rationale behind layer thickness calculation
The script assumes that the region of interest comes from cutting a cylinder perpendicular and offset to the main axis. 


[**Offset_layer_thickness_computation.ipynb**](Offset_layer_thickness_computation.ipynb) contains calculation of the thickness based on the trigonometrical model shown on the attached figure: 

![Implant layers calculation](images/Offset_cut_model_dark.png)

## How to use

The groovy scripts are meant to be [ran from QuPath](https://qupath.readthedocs.io/en/stable/docs/scripting/workflows_to_scripts.html). 
