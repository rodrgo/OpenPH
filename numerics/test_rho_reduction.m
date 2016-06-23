
init;

figure_dir = './figures/';
figure_tag = 'morozov';

% Plot params
LW = 'LineWidth';
MS = 'MarkerSize';
markers = '+o*.xsd^v><ph';

% Algorithms to test
order = 4;

stream = examples.MorozovCubicTimeExample.getMorozovCubicTimeExample(order);
D = BoundaryMatrix(stream, 'plain');

[lows_std, t] = D.reduce('std_sparse');
[lows_rho, t] = D.reduce('rho_reduction');

norm(sum(lows_std ~= lows_rho))

