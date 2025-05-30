% % Howe Truss
% Model B Assignment – Structures
clear all; close all; clc;
%% MATERIAL AND GEOMETRY DATA
E = 11e9;            % Approximate Young's modulus for wood (Pa)
rho = 500;           % Wood density (kg/m³)
g = 9.81;            % Gravity acceleration (m/s²)
A = 0.01;            % Cross-sectional area (m²)
%% NODES AND ELEMENTS
nodes = [0 0; 2 0; 4 0; 6 0; 8 0; 2 1; 6 1; 4 2];
elements = [1 2; 2 3; 3 4; 4 5;
           1 6; 6 8; 8 7; 7 5;
           2 6; 3 8; 4 7; 6 3; 7 3];
nn = size(nodes,1);
ne = size(elements,1);
ndof = 2;                  % DOF per node (u_x and u_y)
ngdl = nn * ndof;          % total DOF in the structure
K = zeros(ngdl);           % global stiffness matrix (ngdl×ngdl)
f = zeros(ngdl,1);         % global force vector (ngdl×1)
% Gauss points and weights (2-point integration)
xi_gauss = [-1/sqrt(3),  1/sqrt(3)];
w_gauss = [1, 1];
%% ASSEMBLE GLOBAL STIFFNESS MATRIX
for e = 1:ne
  ni = elements(e,1); nj = elements(e,2);
  xi = nodes(ni,1); yi = nodes(ni,2);
  xj = nodes(nj,1); yj = nodes(nj,2);
  L = norm([xj - xi, yj - yi]);
  c = (xj - xi)/L;
  s = (yj - yi)/L;
  J = L / 2;
  % Local stiffness matrix integration
  k_local = zeros(2);
  for gp = 1:length(xi_gauss)
      wi = w_gauss(gp);
      dN_dxi = [-0.5; 0.5];
      B = dN_dxi / J;
      k_local = k_local + (B * B') * E * A * J * wi;
  end
  % Rotate to global coordinates and assemble
  T_rot = [ c,  s, 0, 0;
            0,  0, c, s];
  k_global = T_rot' * k_local * T_rot;
  gdls = [2*ni-1, 2*ni, 2*nj-1, 2*nj];
  K(gdls, gdls) = K(gdls, gdls) + k_global;
  % Self-weight load as vertical nodal forces
  w_line = rho * A * g;
  fy = w_line * L / 2;
  f(2*ni) = f(2*ni) - fy;
  f(2*nj) = f(2*nj) - fy;
end
% Additional vertical loads (snow and slate, 2kN/m)
w_extra = 1000 + 1000;  % N/m
roof_bars = [5 6 7 8];
for e = roof_bars
  ni = elements(e,1); nj = elements(e,2);
  xi = nodes(ni,1); yi = nodes(ni,2);
  xj = nodes(nj,1); yj = nodes(nj,2);
  L = norm([xj - xi, yj - yi]);
  fy = w_extra * L / 2;
  f(2*ni) = f(2*ni) - fy;
  f(2*nj) = f(2*nj) - fy;
end
%% BOUNDARY CONDITIONS (node 1 fixed, node 5 on roller)
fixed = [1 2 10];  % Node 1 fixed in X and Y (1,2) and node 5 fixed only in Y (10)
free = setdiff(1:ngdl, fixed);
u = zeros(ngdl,1);
u(free) = K(free, free) \ f(free);
%% CALCULATION OF AXIAL FORCES IN MEMBERS
axial_force = zeros(ne,1);
for e = 1:ne
  ni = elements(e,1); nj = elements(e,2);
  xi = nodes(ni,1); yi = nodes(ni,2);
  xj = nodes(nj,1); yj = nodes(nj,2);
  L = norm([xj - xi, yj - yi]);
  c = (xj - xi)/L; s = (yj - yi)/L;
  gdls = [2*ni-1 2*ni 2*nj-1 2*nj];
  u_elem = u(gdls);
  delta = [-c -s c s] * u_elem;
  axial_force(e) = E * A / L * delta;
end
%% GLOBAL RESULTS
u_disp = u(1:2:end);
v_disp = u(2:2:end);
total_disp = sqrt(u_disp.^2 + v_disp.^2);
[max_disp, idx_max] = max(total_disp);
axial_max = max(axial_force);
axial_min = min(axial_force);
fprintf('Axial max: %+8.2f N\n', axial_max);
fprintf('Axial min: %+8.2f N\n', axial_min);
H = ne + length(fixed) - 2*nn;
fprintf('Hyperstatic degree H = %d\n', H);
if H == 0
  fprintf('The structure is STATCIALLY DETERMINATE (H = 0).\n');
elseif H > 0
  fprintf('The structure is HYPERSTATIC with %d redundancies.\n', H);
else
  fprintf('The structure is a MECHANISM (missing %d constraints).\n', -H);
end
%% FIGURE 1: ORIGINAL VS. DEFORMED STRUCTURE (WITHOUT COLORS)
scale = 100;
u_disp = u(1:2:end);
v_disp = u(2:2:end);
nodes_def = nodes + scale * [u_disp v_disp];
figure(1); clf; hold on; axis equal;
title(['Figure 1: Original and deformed structure (x' num2str(scale) ')']);
xlabel('X [m]'); ylabel('Y [m]');
set(gca, 'FontSize', 10);
for e = 1:ne
  ni = elements(e,1);
  nj = elements(e,2);
  % Coordinates
  x_orig = [nodes(ni,1), nodes(nj,1)];
  y_orig = [nodes(ni,2), nodes(nj,2)];
  x_def = [nodes_def(ni,1), nodes_def(nj,1)];
  y_def = [nodes_def(ni,2), nodes_def(nj,2)];
  % Plot original and deformed
  plot(x_orig, y_orig, 'k--', 'LineWidth', 1);
  plot(x_def, y_def, 'r-', 'LineWidth', 2);
end
legend('Original', 'Deformed', 'Location', 'best');
figure(11); clf; hold on; axis equal;
title('Figure 1.1: Axial force map on deformed structure');
xlabel('X [m]'); ylabel('Y [m]');
set(gca,'FontSize',10);
colormap jet
cmap     = colormap;
n_colors = size(cmap,1);
axial_abs_max = max(abs(axial_force));  % absolute maximum
% Min and max line widths
lw_min = 0.5;
lw_max = 5;
for e = 1:ne
  ni = elements(e,1); nj = elements(e,2);
  x_def = [nodes_def(ni,1), nodes_def(nj,1)];
  y_def = [nodes_def(ni,2), nodes_def(nj,2)];
   % Normalize force between –1 and +1
  t = axial_force(e) / axial_abs_max;
  % Color index
  idx = round(((t + 1)/2)*(n_colors-1)) + 1;
  col = cmap(min(max(idx,1),n_colors),:);
   % Line width proportional to |force|
  lw = lw_min + (lw_max - lw_min)*(abs(axial_force(e)) / axial_abs_max);
   % Plot
  plot(x_def, y_def, '-', 'LineWidth', lw, 'Color', col);
   % Label value
  xm = mean(x_def); ym = mean(y_def);
  text(xm, ym, sprintf('%+.0f N', axial_force(e)), ...
       'FontSize', 8, 'Color', 'k', 'HorizontalAlignment','center');
end
cb = colorbar;
caxis([-axial_abs_max, axial_abs_max]);
ylabel(cb, 'Axial force [N]');
%% FIGURE 2: NODAL DISPLACEMENT MAP (ARROWS)
figure(2); clf; hold on; axis equal;
title('Figure 2: Nodal displacement map');
xlabel('X [m]'); ylabel('Y [m]');
set(gca, 'FontSize', 10);
% Plot original structure
for e = 1:ne
  ni = elements(e,1); nj = elements(e,2);
  plot([nodes(ni,1), nodes(nj,1)], [nodes(ni,2), nodes(nj,2)], 'k-', 'LineWidth', 1);
end
% Arrow amplification factor
scale_arrow = 500;
% Draw displacements as arrows
for i = 1:nn
  x = nodes(i,1);
  y = nodes(i,2);
  ux = u(2*i - 1);
  uy = u(2*i);
  quiver(x, y, scale_arrow*ux, scale_arrow*uy, 0, ...
      'Color','r', 'MaxHeadSize', 0.5, 'LineWidth', 1.5);
end
% Optional node labels
for i = 1:nn
  text(nodes(i,1), nodes(i,2)+0.15, sprintf('N%d', i), ...
      'FontSize', 9, 'HorizontalAlignment','center');
end
legend('Original structure', 'Displacements (amplified)');
%% FIGURE 3: IMPROVED NODAL DISPLACEMENT TABLE
figure(3); clf; axis off;
title('Figure 3: Nodal displacements table', 'FontSize', 12);
% Headers
table_txt = cell(nn + 3, 1);
table_txt{1} = 'NODAL DISPLACEMENTS';
table_txt{2} = '----------------------------';
table_txt{3} = sprintf('%-6s %-15s %-15s %-18s', 'Node', 'u_x [m]', 'u_y [m]', 'Total disp. [m]');
% Table content
for i = 1:nn
  ux = u(2*i - 1);
  uy = u(2*i);
  utot = sqrt(ux^2 + uy^2);
  table_txt{i + 3} = sprintf('%-6d %-15.4e %-15.4e %-18.4e', i, ux, uy, utot);
end
% Display as text
text(0.05, 1, table_txt, 'FontName','Courier', 'FontSize', 10, 'VerticalAlignment','top');
%% FIGURE 4: ENHANCED AXIAL FORCE TABLE WITH LENGTH AND TYPE
tol_axial = 1e-3;  % Threshold to consider zero force
figure(4); clf; axis off;
title('Figure 4: Axial forces in members table', 'FontSize', 12);
% Headers
table_axial = cell(ne + 3, 1);
table_axial{1} = 'AXIAL FORCES IN MEMBERS';
table_axial{2} = '--------------------------------------------------------------';
table_axial{3} = sprintf('%-6s %-8s %-8s %-15s %-12s %-12s', ...
                       'Elem', 'Node_i', 'Node_j', 'Axial [N]', 'Length [m]', 'Type');
% Table rows
for e = 1:ne
  ni = elements(e,1);
  nj = elements(e,2);
  xi = nodes(ni,1); yi = nodes(ni,2);
  xj = nodes(nj,1); yj = nodes(nj,2);
  L = norm([xj - xi, yj - yi]);
  N = axial_force(e);
  type = 'U';  % Undetermined by default
  if abs(N) < tol_axial
      N = 0;
      type = 'U';  % U = Unloaded
  elseif N > 0
      type = 'T';  % Tension
  elseif N < 0
      type = 'C';  % Compression
  end
  table_axial{e + 3} = sprintf('%-6d %-8d %-8d %+15.4e %-12.4f %-12s', ...
                               e, ni, nj, N, L, type);
end
% Display as text
text(0.05, 1, table_axial, 'FontName','Courier', 'FontSize', 10, 'VerticalAlignment','top');
%% FIGURE 5: LOAD SCHEME
figure(5); clf; hold on; axis equal;
title('Figure 5: Applied load scheme');
xlabel('X [m]'); ylabel('Y [m]');
set(gca, 'FontSize', 10);
% Draw structure
for e = 1:ne
  ni = elements(e,1); nj = elements(e,2);
  plot([nodes(ni,1), nodes(nj,1)], [nodes(ni,2), nodes(nj,2)], 'k-', 'LineWidth', 1.2);
end
% Members with distributed vertical load
roof_bars = [5 6 7 8];
w_total = 2049.05;  % N/m (self-weight + snow + slate)
for e = roof_bars
  ni = elements(e,1);
  nj = elements(e,2);
  xi = nodes(ni,1); yi = nodes(ni,2);
  xj = nodes(nj,1); yj = nodes(nj,2);
  xm = mean([xi, xj]);  % midpoint
  ym = mean([yi, yj]);
  % Draw downward arrow
  quiver(xm, ym + 0.3, 0, -0.3, 0, 'Color', 'blue', 'LineWidth', 1.5, 'MaxHeadSize', 2);
   % Label with value
  text(xm, ym + 0.35, sprintf('%.0f N/m', w_total), ...
       'HorizontalAlignment','center', 'FontSize', 9, 'Color','blue');
end
% Label the load type
text(0.5, max(nodes(:,2))+0.8, 'Vertical load = Self-weight + Snow + Slate', ...
   'FontSize', 10, 'Color', 'black');
% Optional: node numbering
for i = 1:nn
  text(nodes(i,1), nodes(i,2) - 0.2, sprintf('N%d', i), ...
      'FontSize', 8, 'HorizontalAlignment','center');
end
%% FIGURE 6: SUPPORTS AND BOUNDARY CONDITIONS SCHEME
figure(6); clf; hold on; axis equal;
title('Figure 6: Supports and boundary conditions');
xlabel('X [m]'); ylabel('Y [m]');
set(gca, 'FontSize', 10);
% Draw members (structure)
h_bars = gobjects(1,1);
for e = 1:ne
  ni = elements(e,1);
  nj = elements(e,2);
  h_bars = plot([nodes(ni,1), nodes(nj,1)], [nodes(ni,2), nodes(nj,2)], 'k-', 'LineWidth', 1.2);
end
% Draw nodes
h_nodes = plot(nodes(:,1), nodes(:,2), 'ko', 'MarkerFaceColor','k');
% Fixed support – Node 1
n1 = 1; x1 = nodes(n1,1); y1 = nodes(n1,2);
h_fixed = plot(x1, y1, 's', 'MarkerSize', 10, 'MarkerEdgeColor','b', 'MarkerFaceColor','b');
text(x1, y1 - 0.25, 'Fixed support (u_x, u_y)', 'HorizontalAlignment','center', 'FontSize', 9, 'Color','b');
% Roller support – Node 5
n5 = 5; x5 = nodes(n5,1); y5 = nodes(n5,2);
h_roller = plot(x5, y5, 'o', 'MarkerSize', 10, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
text(x5, y5 - 0.25, 'Roller support (u_y)', 'HorizontalAlignment','center', 'FontSize', 9, 'Color','r');
% Node labels
for i = 1:nn
  text(nodes(i,1), nodes(i,2) + 0.2, sprintf('N%d', i), 'HorizontalAlignment','center', 'FontSize', 9);
end
% Clear legend
legend([h_bars, h_nodes, h_fixed, h_roller], ...
     {'Members', 'Nodes', 'Fixed support', 'Roller support'}, ...
     'Location', 'northeastoutside');
%% FIGURE 7: GLOBAL RESULTS SUMMARY
figure(7); clf; axis off;
title('Figure 7: Global results summary', 'FontSize', 12);
% Summary content
summary = {
  'GLOBAL RESULTS SUMMARY';
  '--------------------------------------------';
  sprintf('Node with max displacement: %d', idx_max);
  sprintf('Max total displacement: %.4e m', max_disp);
  sprintf('Max axial force: %+10.2f N', axial_max);
  sprintf('Min axial force: %+10.2f N', axial_min);
};
% Display as text
text(0.05, 1, summary, 'FontName','Courier', 'FontSize', 10, 'VerticalAlignment','top');
%% CALCULATION OF HYPERSTATIC DEGREE
b = ne;               % number of members
r = length(fixed);    % number of reaction unknowns (fixed components)
j = nn;               % number of nodes
H = b + r - 2*j;      % hyperstaticity formula
fprintf('\nHyperstatic degree H = %d\n', H);
if H == 0
  fprintf('The structure is STATCIALLY DETERMINATE (H = 0).\n');
elseif H > 0
  fprintf('The structure is HYPERSTATIC with %d redundancies.\n', H);
else
  fprintf('The structure is a MECHANISM (missing %d constraints).\n', -H);
end
%% Tim Zurbiggen & Roger Marqués – 20 May 2025
% Structural Engineering – Master’s in Civil Engineering

