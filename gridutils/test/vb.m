function vb(fname_z, fname_x, fname_y, nce1, nce2, zmax);
    if ~exist(fname_z, 'file')
        fprintf(1, '"%s" not found\n', fname_z);
        return
    end
    if ~exist(fname_x, 'file')
        fprintf(1, '"%s" not found\n', fname_x);
        return
    end
    if ~exist(fname_y, 'file')
        fprintf(1, '"%s" not found\n', fname_y);
        return
    end
  
    z = load(fname_z);
    x = load(fname_x);
    y = load(fname_y);
  
    nx = nce1 * 2 + 1;
    ny = nce2 * 2 + 1;
  
    if size(z) ~= nce1 * nce2
        error('size(z) ~= nce1 * nce2');
    end
    if size(x) ~= nx * ny
        error('size(x) ~= (nce1 * 2 + 1) * (nce2 * 2 + 1)');
    end
    if size(y) ~= nx * ny
        error('size(y) ~= (nce1 * 2 + 1) * (nce2 * 2 + 1)');
    end
  
    if ~exist('zmax', 'var')
        zmax = max(z);
    end
  
    hold on;
    map = thecolormap;
  
    iz = 1;
    ix = 1;
    for j = 1 : nce2
        for i = 1 : nce1
            i1 = ix;
            i2 = ix + 2;
            i3 = ix + nx * 2 + 2;
            i4 = i3 - 2;
            xx = [x(i1) x(i2) x(i3) x(i4)];
            yy = [y(i1) y(i2) y(i3) y(i4)];
            sum = x(i1) + x(i2) + x(i3) + x(i4) + z(iz);
            if ~isnan(sum)
                col = map(cindex(z(iz), zmax), :);
                fill(xx, yy, col);
            end
            ix = ix + 2;
            iz = iz + 1;
        end
        ix = ix + nx + 1;
    end
    axis equal;
    set(gca, 'XTick', [], 'YTick', []);
    xrange = [min(x) max(x)];
    xlim([xrange(1) - (xrange(2) - xrange(1)) * 0.05, xrange(2) + (xrange(2) - xrange(1)) * 0.05]);
    yrange = [min(y) max(y)];
    ylim([yrange(1) - (yrange(2) - yrange(1)) * 0.05, yrange(2) + (yrange(2) - yrange(1)) * 0.05]);
    box on;
    hold off;

    return;
  
function index = cindex(z, zmax)
  
    if z <= 0 | isnan(z)
        index = 1;
    elseif z >= zmax
        index = 88;
    else
        index = floor(z / zmax * 88 + 1);
    end

    return

function map = thecolormap();  
    map = [ 
        1.0000    1.0000    1.0000
        0.0000    0.0000    0.5625
        0.0000    0.0000    0.6250
        0.0000    0.0000    0.6875
        0.0000    0.0000    0.7500
        0.0000    0.0000    0.8125
        0.0000    0.0000    0.8750
        0.0000    0.0000    0.9375
        0.0000    0.0000    1.0000
        0.0000    0.0625    1.0000
        0.0000    0.1250    1.0000
        0.0000    0.1875    1.0000
        0.0000    0.2500    1.0000
        0.0000    0.3125    1.0000
        0.0000    0.3750    1.0000
        0.0000    0.4375    1.0000
        0.0000    0.5000    1.0000
        0.0000    0.5625    1.0000
        0.0000    0.6250    1.0000
        0.0000    0.6875    1.0000
        0.0000    0.7500    1.0000
        0.0000    0.8125    1.0000
        0.0000    0.8750    1.0000
        0.0000    0.9375    1.0000
        0.0000    1.0000    1.0000
        0.0625    1.0000    1.0000
        0.1250    1.0000    0.9375
        0.1875    1.0000    0.8750
        0.2500    1.0000    0.8125
        0.3125    1.0000    0.7500
        0.3750    1.0000    0.6875
        0.4375    1.0000    0.6250
        0.5000    1.0000    0.5625
        0.5625    1.0000    0.5000
        0.6250    1.0000    0.4375
        0.6875    1.0000    0.3750
        0.7500    1.0000    0.3125
        0.8125    1.0000    0.2500
        0.8750    1.0000    0.1875
        0.9375    1.0000    0.1250
        1.0000    1.0000    0.0625
        1.0000    1.0000    0.0000
        1.0000    0.9375    0.0000
        1.0000    0.8750    0.0000
        1.0000    0.8125    0.0000
        1.0000    0.7500    0.0000
        1.0000    0.6875    0.0000
        1.0000    0.6250    0.0000
        1.0000    0.5625    0.0000
        1.0000    0.5000    0.0000
        1.0000    0.4375    0.0000
        1.0000    0.3750    0.0000
        1.0000    0.3125    0.0000
        1.0000    0.2500    0.0000
        1.0000    0.1875    0.0000
        1.0000    0.1250    0.0000
        1.0000    0.0625    0.0000
        1.0000    0.0000    0.0000
        1.0000    0.0000    0.0000
        0.9375    0.0000    0.0000
        0.9375    0.0000    0.0000
        0.9375    0.0000    0.0000
        0.8750    0.0000    0.0000
        0.8750    0.0000    0.0000
        0.8750    0.0000    0.0000
        0.8125    0.0000    0.0000
        0.8125    0.0000    0.0000
        0.7500    0.0000    0.0000
        0.7500    0.0000    0.0000
        0.7500    0.0000    0.0000
        0.6875    0.0000    0.0000
        0.6875    0.0000    0.0000
        0.6875    0.0000    0.0000
        0.6250    0.0000    0.0000     
        0.6250    0.0000    0.0000     
        0.6250    0.0000    0.0000 
        0.5625    0.0000    0.0000
        0.5625    0.0000    0.0000
        0.5625    0.0000    0.0000
        0.5000    0.0000    0.0000
        0.5000    0.0000    0.0000
        0.5000    0.0000    0.0000
        0.4375    0.0000    0.0000
        0.4375    0.0000    0.0000
        0.4375    0.0000    0.0000
        0.3750    0.0000    0.0000
        0.3750    0.0000    0.0000
        0.3750    0.0000    0.0000
          ];
    
    return;
    