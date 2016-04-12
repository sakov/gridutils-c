function viewgrid(fname, f, i0, j0)
    [nx, ny, x, y] = fload(fname);

    xmin = min(x);
    xmax = max(x);
    ymin = min(y);
    ymax = max(y);
    xdiff = xmax - xmin;
    ydiff = ymax - ymin;
    diff = max([xdiff ydiff]) / 20.0;
    xmin = xmin - diff;
    xmax = xmax + diff;
    ymin = ymin - diff;
    ymax = ymax + diff;
  
    if ~exist('f', 'var') | isempty(f)
        figure;
        ht = title(fname);
        set(ht, 'interpreter', 'none');
        hold on;
        axis([xmin xmax ymin ymax]);
        axis image;
      
        fill([xmin xmin xmax xmax], [ymin ymax ymax ymin], [0.7 0.7 0.7]);
        col = [0 0 0];
        alpha = 1;
    else
        figure(f);
        col = [1 0 0];
        alpha = 0.2;
    end
  
    if ~exist('i0', 'var')
        i0 = -1;
        j0 = -1;
    end
  
    xx0 = [];
    yy0 = [];
    ix = 1;
    for j = 1:ny-1
        for i = 1:nx-1
            i1 = ix;
            i2 = ix + 1;
            i3 = ix + nx + 1;
            i4 = i3 - 1;
            xx = [x(i1) x(i2) x(i3) x(i4)];
            yy = [y(i1) y(i2) y(i3) y(i4)];
            sum = x(i1) + x(i2) + x(i3) + x(i4);
            if ~isnan(sum)
                fill(xx, yy, [1 1 1], 'edgecolor', col, 'facealpha', alpha);
                if i == i0 & j == j0
                    xx0 = xx;
                    yy0 = yy;
                end
            end
            ix = ix + 1;
        end
        ix = ix + 1;
    end
    if ~isempty(xx0)
        fill(xx0, yy0, [1 1 1], 'edgecolor', [0 1 0], 'facealpha', alpha);
    end

    return

function [nx, ny, x, y] = fload(fname);
    f = fopen(fname);
    if (f < 0)
        return;
    end
  
    nx = 0;
    ny = 0;
    x = [];
    y = [];
  
    str = fgetl(f);
    [n, count] = sscanf(str, '## %d x %d');

    if (count == 2 & n(1) > 0 & n(2) > 0)
        nx = n(1);
        ny = n(2);
        [xy] = fscanf(f, '%f', [2,inf]);
        x = xy(1,:);
        y = xy(2,:);
    end
  
    fclose(f);
    return
