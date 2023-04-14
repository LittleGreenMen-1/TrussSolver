clear; clc;

% Truss data
global nodes;
global graph;
global tension;

% Material properties
global A;
global E;

% FEM vectors
global K;
global u;
global F;

% Diagram data
global height;
global width;

A = 2; % mm^2
E = 1; % MPa

filename = input("Truss file = ", "s");

loadTruss(filename);

% Calculate the global stiffness matrix3
for i = 1:size(nodes, 1)
    for j = 1:i
        if graph(i, j) == 1
            Kij = getK(i, j);

            ui = 2 * i - 1;
            vi = 2 * i;

            uj = 2 * j - 1;
            vj = 2 * j;

            indexes = [ui, vi, uj, vj];

            for ind1 = 1:4
                for ind2 = 1:4
                    K(indexes(ind1), indexes(ind2)) = K(indexes(ind1), indexes(ind2)) + Kij(ind1, ind2);
                end
            end
        end
    end
end

% Formulation of FE equation

% Remove the rows and columns of the global stiffness matrix that correspond to the fixed supports (0 in the displacement vector)
% Thus, we get a system of equations that gives the displacements of the nodes that are not fixed
sys_A = K;
sys_B = F;

for i = size(u, 1):-1:1
    if u(i) == 0
        sys_A(i, :) = [];
        sys_A(:, i) = [];

        sys_B(i) = [];
    end
end

% Solve the system of equations
X = linsolve(sys_A, sys_B);

% Add the displacements to the displacement vector
ind = 1;

for i = 1:size(u, 1)
    if u(i) == Inf
        u(i) = X(ind);

        ind = ind + 1;
    end
end

% Calculate the nodal forces
% K * u = F

F = K * u;

% Calculate the tensions in the members
for i = 1:size(nodes, 1)
    for j = 1:i
        if graph(i, j) == 1
            stress = getStress(i, j);

            tension(i, j) = stress * A;

            % If tension is neglijable, set it to 0
            if abs(tension(i, j)) < 0.001
                tension(i, j) = 0;
            end

            tension(j, i) = tension(i, j);
        end
    end
end

% Display the calculated data
fprintf("---------------------\n\n");
fprintf("Global stiffness matrix (N/mm):\n\n");
disp(K);

fprintf("---------------------\n\n");
fprintf("Nodal forces (N):\n\n");
for i = 1:2:size(F, 1)
    ind = floor(i / 2) + 1;

    fprintf("F%dx = %6.3fN\n", ind, F(i));
    fprintf("F%dy = %6.3fN\n", ind, F(i + 1));
end

fprintf("---------------------\n\n");
fprintf("Nodal displacements (m):\n\n");
for i = 1:2:size(F, 1)
    ind = floor(i / 2) + 1;

    fprintf("d%dx = %6.3f\n", ind, u(i));
    fprintf("d%dy = %6.3f\n", ind, u(i + 1));
end

fprintf("---------------------\n\n");
fprintf("Member tensions (N):\n\n");
for i = 1:size(nodes, 1)
    for j = 1:i
        if graph(i, j) == 1
            fprintf("T(%d, %d) = %6.3fN\n", j, i, tension(i, j));
        end
    end
end

% Show the solved truss
showTruss();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load a truss from a data file
% Parameters:
%   filename - the name of the data file
%
function loadTruss(filename)
    global nodes;
    global graph;
    global tension;

    global K;
    global u;
    global F;

    global height;
    global width;

    file = fopen(filename, "r");
    
    n = fscanf(file, "%d\n", 1);    % No. of nodes (joints)
    m = 0;                          % No. of members
    r = 0;                          % No. of 'degrees of freedom'

    % Initialize truss arrays
    nodes = zeros(n, 2);
    graph = zeros(n, n);
    tension = Inf(n, n);
 
    % Initialize FEM vectors
    K = zeros(size(nodes, 1) * 2, size(nodes, 1) * 2); % Global stiffness matrix
    u = Inf(n * 2, 1);                                 % Displacement vector
    F = zeros(n * 2, 1);                               % Force vector

    % These will help us compute the height and width of the diagram
    max_x = -Inf;
    max_y = -Inf;

    min_x = Inf;
    min_y = Inf;

    % Read the nodes
    for i = 1:n
        node = fscanf(file, "%d: %f %f %f %f %c\n", [1, 6]);

        % Add the node to the list of nodes
        nodes(node(1), :) = [node(2), node(3)];

        %
        if node(2) > max_x max_x = node(2); end
        if node(3) > max_y max_y = node(3); end

        if node(2) < min_x min_x = node(2); end
        if node(3) < min_y min_y = node(3); end        

        % Add the forces to the force vector
        if node(4)
            F(2 * node(1) - 1) = node(4); % x-component of force
        end

        if node(5)
            F(2 * node(1)) = node(5);    % y-component of force
        end

        % Add the supports to the support vector
        if (node(6) == '*') % Fixed support
            u(2 * node(1) - 1) = 0;
            u(2 * node(1)) = 0;
        end

        if (node(6) == 'x') % Horizontal support
            u(2 * node(1)) = 0;
        end
    end

    height = max_y - min_y;
    width = max_x - min_x;

    fscanf(file, "\n");

    % Read the connections
    for i = 1:n
        node = fscanf(file, "%d: ", 1);

        line = fgetl(file);
        line = split(line);

        for j = 1:length(line)
            connection = str2num(line{j});
            m = m + 1;

            % Add the connection to the list of connections
            graph(node, connection) = 1;
            graph(connection, node) = 1;
        end
    end

    m = m / 2;
    r = sum(~u(:)); % Number of 0 elements in the displacement vector

    fprintf("m = %f\nr = %f\n2 * n = %f\n", m, r, 2 * n);
    
    % % Check if the truss is statically determinate
    % if m + r == 2 * n
    %     fprintf("The truss is statically determinate.\n");
    % else
    %     error("The truss is statically indeterminate.");
    % end

    fclose(file);
end


% Get the stiffness matrix of a member between two nodes
% Parameters:
%   nodeA - the index of the first node
%   nodeB - the index of the second node
%
% Returns:
%   K - the stiffness matrix of the member
%
function K = getK(nodeA, nodeB)
    global nodes;
    global graph;

    global A;
    global E;

    if ~graph(nodeA, nodeB) % If the nodes are not connected
        K = zeros(4, 4);
        return ;
    end

    length = mag(nodes(nodeA, :) - nodes(nodeB, :));  % |AB|

    l = (nodes(nodeB, 1) - nodes(nodeA, 1)) / length; % cos(theta)
    m = (nodes(nodeB, 2) - nodes(nodeA, 2)) / length; % sin(theta)

    M = [l^2, l * m, -l^2, -l * m;
         l * m, m^2, -l * m, -m^2;
         -l^2, -l * m, l^2, l * m;
         -l * m, -m^2, l * m, m^2];

    K = A * E / length .* M;
end


% Get the stress in a member between two nodes
% Parameters:
%   nodeA - the index of the first node
%   nodeB - the index of the second node
%
% Returns:
%   sig - the stress in the member
%
function sig = getStress(nodeA, nodeB)
    global nodes;
    global graph;

    global A;
    global E;

    global u;

    if ~graph(nodeA, nodeB) % If the nodes are not connected
        sig = 0;
        return ;
    end

    length = mag(nodes(nodeA, :) - nodes(nodeB, :));  % |AB|

    l = (nodes(nodeB, 1) - nodes(nodeA, 1)) / length; % cos(theta)
    m = (nodes(nodeB, 2) - nodes(nodeA, 2)) / length; % sin(theta)

    v = [-l -m l m];
    c = [u(2 * nodeA - 1); u(2 * nodeA); u(2 * nodeB - 1); u(2 * nodeB)];

    sig = (E / length) .* (v * c);
end


% Get the magnitude of a vector
% Parameters:
%    v - the vector
%
% Returns:
%    ret - the magnitude of the vector
%
function ret = mag(v)
    ret = sqrt(sum(v .* v));
end


% Display the truss structure
%
function showTruss()
    global nodes;
    global graph;
    global tension;

    global u;
    global F;

    global height;
    global width;

    colors = ['b', 'g', 'r']; % Colors for the members :
                              % b - member is in compression
                              % g - there is no tension in the member
                              % r - member is in tension

    max_force = max(abs(F));

    % Plot truss structure
    for i = 1:size(nodes, 1)
        axis off;

        nodeA = nodes(i, :);

        % Joint
        plot(nodeA(1), nodeA(2), '.k', 'MarkerSize', 20);

        % Members
        for j = 1:i
            if graph(i, j) == 1
                nodeB = nodes(j, :);

                % Plot member
                color_index = sign(tension(i, j)) + 2;
                plot([nodeA(1), nodeB(1)], [nodeA(2), nodeB(2)], colors(color_index), 'LineWidth', 3);
                
                % Plot text to show tension in member
                mid = [(nodeA(1) + nodeB(1)), (nodeA(2) + nodeB(2))] ./ 2;

                % Get the angle of the member -> rotate text by this angle
                theta = atan2(nodeA(2) - nodeB(2), nodeA(1) - nodeB(1)) * 180 / pi;
                if theta > 90 theta = theta - 180; end % if the text would be upside down, rotate it 180 degrees
                
                text(mid(1), mid(2), sprintf("%.2fN", tension(i, j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Rotation', theta);
            end
        end

        % Forces
        f = F(2 * i - 1:2 * i);

        h_dir = ['right', 'left'];
        v_dir = ['top', 'bottom'];

        % Plot Ox force
        % if abs(f(1)) > 0.01
        %     if f(1) < 0 dir = -width; text_dir = 'right'; else dir = width; text_dir = 'left'; end

        %     dir = dir * abs(f(1)) / max_force;
            
        %     quiver(nodes(i, 1), nodes(i, 2), .4 * dir, 0, 0);
        %     text(nodes(i, 1) + .4 * dir, nodes(i, 2), sprintf("%.1fN", f(1)), 'HorizontalAlignment', text_dir, 'VerticalAlignment', 'middle', 'FontSize', 10);
        % end

        % Plot Oy force
        if abs(f(2)) > 0.01
            if f(2) < 0 dir = -height; text_dir = 'top'; else dir = height; text_dir = 'bottom'; end
            
            dir = dir * abs(f(2)) / max_force;

            quiver(nodeA(1), nodeA(2), 0, .4 * dir, 0);
            text(nodeA(1), nodeA(2) + .4 * dir, sprintf("%.1fN", f(2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', text_dir, 'FontSize', 10);
        end

        hold on;
    end

    % Plot supports
    for i = 1:2:size(u, 1)
        if u(i) == 0 || u(i + 1) == 0            
            node = (i + 1) / 2;
    
            p = nsidedpoly(3, 'Center', [nodes(node, 1), nodes(node, 2) - 0.1 * sin(pi / 3) * 2 / 3], 'SideLength', 0.1);
            plot(p);
        end
    end

    axis equal;
    hold off;
end
