#!/bin/bash
set -eu

if [ ! -x ../getnodes ]
then
    echo "error: no ../getnodes found"
    echo 'run "./configure" and "make" in the upper level directory'
    exit 1
fi

echo
echo "input double density grid: gridnodes.txt"
echo

echo "1. Validating and printing some stats:"
echo
../getnodes gridpoints_DD-raw.txt -i DD -o DD -v > gridpoints_DD.txt
echo
echo "     (gridpoints_DD-raw.txt -> gridpoints_DD.txt)"
echo

echo "2. Testing conversion between the grid types:"
echo -n "   Extracting cell corner nodes from a double-density grid..."
../getnodes gridpoints_DD.txt -i DD -o CO > gridpoints_CO.txt
echo "done"
echo "     (gridpoints_DD.txt -> gridpoints_CO.txt)"
echo -n "   Extracting centre nodes from corner nodes..."
../getnodes gridpoints_CO.txt -i CO -o CE > gridpoints_CE.txt
echo "done"
echo "     (gridpoints_CO.txt -> gridpoints_CE.txt)"
echo -n "   Recovering corner nodes from centre nodes..."
../getnodes gridpoints_CE.txt -i CE -o CO > gridpoints_CO2.txt
echo "done"
echo "     (gridpoints_CE.txt -> gridpoints_CO2.txt)"
echo -n "   Recovering double-density nodes from centre nodes..."
../getnodes gridpoints_CE.txt -i CE -o DD > gridpoints_DD2.txt
echo "done"
echo "     (gridpoints_CE.txt -> gridpoints_DD2.txt)"
echo

echo -n "3. Getting the boundary..."
../getbound gridpoints_DD.txt > bound.txt
echo "done"
echo "     (gridpoints_DD.txt -> bound.txt)"
echo

echo -n "4. Getting the boundary in index space..."
../getbound gridpoints_DD.txt -r > bound-r.txt
echo "done"
echo "     (gridpoints_DD.txt -> bound-r.txt)"
echo

echo -n "5. As above, with compacting..."
../getbound gridpoints_DD.txt -c -r > bound-c-r.txt
echo "done"
echo "     (gridpoints_DD.txt -> bound-c-r.txt)"
echo

echo "6. Converting a few points to and from index space:"
echo "   point 1:"
echo -n '     513252.3881 5186890.274 -> '
echo "513252.3881 5186890.274" | ../xy2ij -g gridpoints_DD.txt -o stdin
echo "     and back:"
echo -n "     (index) "
echo "513252.3881 5186890.274" | ../xy2ij -g gridpoints_DD.txt -o stdin |tr -d "\n"
echo -n '-> '
echo "513252.3881 5186890.274" | ../xy2ij -g gridpoints_DD.txt -o stdin | ../xy2ij -g gridpoints_DD.txt -o stdin -r
echo "   point 2:"
echo -n '     (index) 20.5 10.5 -> '
echo "20.5 10.5" | ../xy2ij -g gridpoints_DD.txt -o stdin -r
echo "     and back:"
echo -n "     "`echo "20.5 10.5" | ../xy2ij -g gridpoints_DD.txt -o stdin -r |tr -d "\n"`
echo -n '-> '
echo "20.5 10.5" | ../xy2ij -g gridpoints_DD.txt -o stdin -r | ../xy2ij -g gridpoints_DD.txt -o stdin
echo

echo "7. As p.6, using mapping via kd-tree:"
echo "   point 1:"
echo -n '     513252.3881 5186890.274 -> '
echo "513252.3881 5186890.274" | ../xy2ij -g gridpoints_DD.txt -o stdin -k
echo "     and back:"
echo -n "     (index) "
echo "513252.3881 5186890.274" | ../xy2ij -g gridpoints_DD.txt -o stdin -k |tr -d "\n"
echo -n '-> '
echo "513252.3881 5186890.274" | ../xy2ij -g gridpoints_DD.txt -o stdin -k | ../xy2ij -g gridpoints_DD.txt -o stdin -r
echo "   point 2:"
echo -n '     (index) 20.5 10.5 -> '
echo "20.5 10.5" | ../xy2ij -g gridpoints_DD.txt -o stdin -r -k
echo "     and back:"
echo -n "     "`echo "20.5 10.5" | ../xy2ij -g gridpoints_DD.txt -o stdin -r -k |tr -d "\n"`
echo -n '-> '
echo "20.5 10.5" | ../xy2ij -g gridpoints_DD.txt -o stdin -r -k | ../xy2ij -g gridpoints_DD.txt -o stdin
echo

if [ -x ../gridbathy ]
then
    echo -n "8. Interpolating bathymetry with bivariate cubic spline..."
    ../gridbathy -b bathy.txt -g gridpoints_DD.txt > bathy-cs.txt
    echo "done"
    echo "     (bathy.txt -> bathy-cs.txt)"
    echo

    echo -n "9. Interpolating bathymetry with linear interpolation..."
    ../gridbathy -b bathy.txt -g gridpoints_DD.txt -a 3 > bathy-l.txt
    echo "done"
    echo "     (bathy.txt -> bathy-l.txt)"
    echo

    echo -n "10. Interpolating bathymetry with Natural Neighbours interpolation..."
    ../gridbathy -b bathy.txt -g gridpoints_DD.txt -a 2 > bathy-nn.txt
    echo "done"
    echo "     (bathy.txt -> bathy-nn.txt)"
    echo

    echo -n "11. Interpolating bathymetry with Non-Sibsonian NN interpolation..."
    ../gridbathy -b bathy.txt -g gridpoints_DD.txt -a 1 > bathy-ns.txt
    echo "done"
    echo "     (bathy.txt -> bathy-ns.txt)"
    echo
else
    echo "no ../gridbathy found"
    echo "omitting tests for gridbathy"
    echo
fi

echo "  - to visualise a grid, in matlab run \"viewgrid('<grid fname>')\", e.g.:"
echo "      >>viewgrid('gridpoints_CO.txt')"
echo "  - to draw one grid over another, first run \"viewgrid('<grid1 fname>),\""
echo "    and then \"viewgrid('<grid2 fname>', <figure #>)\", e.g.:"
echo "      >>viewgrid('gridpoints_CE.txt')"
echo "      >>viewgrid('gridpoints_CO.txt', 1)"
echo "  - to visualise a grid boundary, in matlab run \"viewbound('bound.txt')\""
echo "    (also \"viewbound('bound-r.txt')\", \"viewbound('bound-c-r.txt')\")"
echo "  - to visualise an interpolated bathymetry, in matlab run \"viewbathy\""
echo "  - to visualise input bathymetry data, in matlab run \"viewbathydata('bathy.txt', 60)\""
