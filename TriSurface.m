classdef TriSurface
    % TriSurface is a class that defines attributes and methods for obtaining the
    % triangulation of a 2D manifold via its level sets embedded in R^3.
    
    properties
        % Triangulation parameters
        max_length = 2;         % Maximum allowed edge length
        max_ratio = 0.00001;  % Maximum allowed edge ratio
        dist_threshold = 1e-4;  % Point-to-line distance threshold between first and last points (eliminates bumps).
        
        % Smoothing parameters
        lambda = 0.12;          % Positive step size for Laplacian smoothing
        mu = -0.125;            % Negative step size for inverse Laplacian smoothing
        smoothing_iter = 15;    % Number of smoothing iterations

        % Lighting parameters 
        DiffuseStrength = 0.3;
        AmbientStrength = 0.7;
        
        % Mesh data
        interp_points;          % Nx3 array storing arclength mesh
        num_points;             % Number of points in the interp_points set corresponding to each level set
        triangles;              % Connectivity array. 

        % Triangulation object
        triangulation
    end

    methods
        function obj = TriSurface(interp_points, num_points)
            % Constructor

            if nargin > 0
                obj.interp_points = interp_points;
                obj.num_points = num_points;
            end
        end

        % Getting and setter methods
        function obj = set_interp_points(obj, new)
            % Setter method to change current mesh and exchange it with 'new'
            obj.interp_points = new;
        end
        
        function obj = set_triangulation(obj, new)
            % Setter method to set triganulation

            obj.triangulation = new;
        end
        
        function obj = set_smoothing_par(obj, new_lambda, new_mu, new_iter)
            if nargin < 3
                new_iter = obj.smoothing_iter;
            end
            obj.smoothing_iter = new_iter;
            obj.lambda = new_lambda;
            obj.mu = new_mu;
        end

        function obj = set_tri_parameters(obj, new_max_length, new_max_ratio, new_dist_tresh)
            obj.max_length = new_max_length;
            obj.max_ratio = new_max_ratio;
            obj.dist_threshold = new_dist_tresh;
        end

        function s = get_triangulation(obj)
            s = obj.triangulation;
        end

        % Main methods
        function obj = apply_mesh_smoothing(obj)
            % Applies mesh smoothing to the attribute interp_points
        
            new_mesh = obj.taubin_smoothing();
            obj = obj.set_interp_points(new_mesh);
        end

        function plot_triangulation(obj, color)
            % Plots the triangulation of the 2D manifold. Also acts as a
            % setter to store the triangulation in a attribute 
            % 

            if nargin < 2
                color = '0.2902 0.7098 0.7412';
            end

            % Builds the triangles via the connectivity set method
            obj.triangles = obj.build_connectivity_set(obj.interp_points, obj.num_points); 
            % Apply mesh smoothing on interp_points 
            obj = obj.apply_mesh_smoothing();  
            % Plots the triangulated surface 
            surf_handle = trisurf(obj.triangles, obj.interp_points(:,1), obj.interp_points(:,2), obj.interp_points(:,3)/(2*pi), 'FaceColor', color, 'EdgeColor', 'black');
            obj.set_triangulation(surf_handle);  % Store surface object.
        end

        function triangles = build_connectivity_set(obj, points, num_points_per_curve)
            % Build connectivity by taking the Delaunay triangulation then
            % eliminating triangles via threshold and ratio criteria

            % Project the 3D points onto the 2D x-y plane by setting z-coordinates to 0
            projected_points = points;
            projected_points(:, 3) = 0;

            % Perform Delaunay triangulation in 2D (x-y plane)
            dt = delaunayTriangulation(projected_points(:, 1:2));

            % Get the connectivity list from the triangulation
            triangles = dt.ConnectivityList;

            % Find the nearest neighbors for each point in the x-y plane
            nearest_neighbors_all = knnsearch(projected_points, projected_points, 'K', num_points_per_curve);
            nearest_neighbors_all = nearest_neighbors_all(:, 2:end); % Exclude the first column, as it is the point itself.
            
            % Initialize a boolean array to mark valid triangles
            valid_triangle_indices = false(size(triangles, 1), 1);
            
            for i = 1:size(triangles, 1)
                triangle_points = triangles(i, :);
                curve_indices = ceil(triangle_points / num_points_per_curve);
            
                % Check if at least two points are on the same curve. 
                % Mesh is well-connected along the curves, and that the triangles connecting the curves are meaningful.
                if max(histcounts(curve_indices)) >= 2
                    valid_triangle_indices(i) = true;
                end
%                 % Check if the triangle meets the threshold and ratio criteria
                  edge_lengths = sqrt(sum((points(triangle_points, :) - points(triangle_points([2, 3, 1]))).^2, 2));
                max_edge_length = max(edge_lengths);
                min_edge_length = min(edge_lengths);
        
                % Set a threshold for the maximum allowed edge length
                max_allowed_edge_length = obj.max_length * mean(edge_lengths(:));
        
                % Set a threshold for the minimum allowed edge length ratio
                % between shortest and longest edge.
                min_allowed_edge_length_ratio = obj.max_ratio;
        
                % Check if the triangle meets the threshold and ratio criteria
                if max_edge_length <= max_allowed_edge_length && min_edge_length / max_edge_length >= min_allowed_edge_length_ratio
                    valid_triangle_indices(i) = true;
                end
            end
            
            % Filter the triangles based on the valid_triangle_indices array
            triangles = triangles(valid_triangle_indices, :);
        end

        function smoothed_points = taubin_smoothing(obj)
            % Input:
            % points: original point coordinates (Nx3 matrix)
            % triangles: connectivity list (Mx3 matrix)
            % num_iterations: number of times to apply the smoothing
            % lambda: positive step size for Laplacian smoothing
            % mu: negative step size for inverse Laplacian smoothing
            
            % Output:
            % points: smoothed point coordinates (Nx3 matrix)
            smoothed_points = obj.interp_points;

            % Loop through the desired number of iterations
            for iter = 1:obj.smoothing_iter
                % Apply Laplacian smoothing with the given lambda
                smoothed_points = laplacian_step(smoothed_points, obj.triangles, obj.lambda);

                % Apply inverse Laplacian smoothing with the given mu
                smoothed_points = laplacian_step(smoothed_points, obj.triangles, obj.mu);
            end
            

            function new_positions = laplacian_step(points, triangles, step_size)
                % Calculate new positions for each point in points array by taking a
                % weighted average of the positions of neighboring vertices and
                % applying the given step size
                
                new_positions = zeros(size(points));
                
                for i = 1:size(points, 1)
                       % Find the neighboring vertices for the current point
                    neighbor_indices = unique(triangles(any(triangles == i, 2), :));
                    neighbor_indices(neighbor_indices == i) = [];

                    % Compute the weights for each neighbor based on the inverse of
                    % the Euclidean distance
                    weights = 1 ./ vecnorm(points(neighbor_indices, :) - points(i, :), 2, 2);
                    weights = weights / sum(weights);

                    % Calculate the new position for the current point by taking a
                    % weighted average of its neighbors' positions
                    new_positions(i, :) = points(i, :) + step_size * (weights' * (points(neighbor_indices, :) - points(i, :)));
                end
            end
        end
        
    end
end
        
