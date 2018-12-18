subplot(2, 2, 1)
vb('bathy-cs.txt', 'x.txt', 'y.txt', 102, 140, 60);
set(gca, 'position', [0.10 0.55 0.40 0.40]);
title('Approximation by C_3^1 spline');

subplot(2, 2, 2)
vb('bathy-nn.txt', 'x.txt', 'y.txt', 102, 140, 60);
set(gca, 'position', [0.50 0.55 0.40 0.40]);
title('Natural Neighbours interpolation');

subplot(2, 2, 3)
vb('bathy-ns.txt', 'x.txt', 'y.txt', 102, 140, 60);
set(gca, 'position', [0.10 0.05 0.40 0.40]);
title('Non-Sibsoninan NN interpolation');

subplot(2, 2, 4)
vb('bathy-l.txt', 'x.txt', 'y.txt', 102, 140, 60);
set(gca, 'position', [0.50 0.05 0.40 0.40]);
title('Linear interpolation');
