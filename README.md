## gridutils ##

**gridutils** provides C library functions and command line utilities for working with curvilinear grids.
**gridutils** has been developed and used mainly for grids generated by [gridgen](https://github.com/sakov/gridgen-c), but can be used to handle arbitrary 2D structured quadrilateral simply connected multi-corner grids.
The recent addition (from v. 1.02) of kd-tree based mappings makes it now possible to map (xy <-> ij) non-simply-connected grids, such as tri-polar ORCA grids.

Contains the following utilities:

  - **getbound** - extracts the bounding polygon
  - **getnodes** - validates and converts grids between double density, center and corner formats
  - **gridbathy** - interpolates scattered data onto the grid
  - **insertgrid** - inserts nodes from one grid into another at a specified location in index space
  - **setbathy** - given a three-column data file (X Y Z), sets Z values in specified index range to the specified value
  - **subgrid** - extracts nodes for a subgrid specified by index range
  - **xy2ij** - converts point coordinates between physical and index spaces

To build the interpolation utility **gridbathy** one needs to install [nn](https://github.com/sakov/nn-c), a Natural 
Neighbours interpolation library and [csa](https://github.com/sakov/nn-c), a cubic spline approximation library. Apart 
from that **gridutils** does not depend on other codes.

Checkout **gridutils** by running `git clone https://github.com/sakov/gridutils-c`.
