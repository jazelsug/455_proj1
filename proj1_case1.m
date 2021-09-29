% Name: proj1_case1.m
% Author: Jazel A. Suguitan
% Last Modified: Sept. 29, 2021

clc,clear
close all

% CASE 1: Algorithm 1 (MSN Fragmentation)

% --- SET PARAMETERS ---

d = 15; % Set desired distance among sensor nodes
k_scale = 1.2;  % Set the scale of MSN: 
r = k_scale * d;  % Set the active range
r_prime = .22 * k_scale * r;    % Set the active range of beta agent
epsilon = 0.1;  % Set a constant for sigma norm
num_nodes = 100;    % Set number of sensor nodes
n = 2;  % Set number of dimensions
%nodes = load('node_distribution2.dat'); % distributed in 2D
nodes = 150.*rand(num_nodes,n)+150.*repmat([0 1],num_nodes,1);  % Randomly generate initial positions of MSN
p_nodes = zeros(num_nodes,n);   % Set initial velocties of MSN
delta_t_update = 0.008;  % Set time step
t = 0:delta_t_update:7; % Set simulation time

nodes_old = nodes; %KEEP privious positions of MSN
q_mean = zeros(size(t,2),n);%Save positions of COM (Center of Mass)
p_mean = zeros(size(t,2),n);%Save velocities of COM (Center of Mass)
Connectivity = []; %save connectivity of MSN
q_nodes_all = cell(size(t,2),num_nodes);
p_nodes_all = cell(size(t,2),num_nodes);
nFrames = 20; %set number of frames for the movie
mov(1:nFrames) = struct('cdata', [],'colormap', []); % Preallocate movie structure.

% %TEST - DELETE LATER
% [Nei_agent, A] = findNeighbors(nodes, r);
% celldisp(Nei_agent)
% A(1:10)
% nodes(1,:), nodes(2,:)
% aij_test = aij(nodes(1,:), nodes(2,:), epsilon, r)
% nij_test = nij(nodes(1,:), nodes(2,:), epsilon)
% [Ui] = inputcontrol_Algorithm1(nodes, Nei_agent, num_nodes, epsilon, r, d, p_nodes, n); % CHECK


%[Nei_agent, Nei_beta_agent, p_ik, q_ik, A] = findNeighbors(nodes_old,nodes,r, r_prime,obstacles, Rk, n, p_nodes,delta_t_update)
% r_prime, obstacles, Rk doesn't matter for case 1

for iteration =1:length(t)
    [Nei_agent, A] = findNeighbors(nodes, r);
    [Ui] = inputcontrol_Algorithm1(nodes, Nei_agent, num_nodes, epsilon, r, d, p_nodes, n); % CHECK
    p_nodes = (nodes - nodes_old)/delta_t_update; %COMPUTE velocities of sensor nodes
    p_nodes_all{iteration} = p_nodes; %SAVE VELOCITY OF ALL NODES
    nodes_old = nodes;
    nodes = nodes_old + p_nodes*delta_t_update  + Ui*delta_t_update*delta_t_update /2;
    q_mean(iteration,:) = mean(nodes); %Compute position of COM of MSN
    plot(q_mean(:,1),q_mean(:,2),'ro','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','k','MarkerSize',4.2)
    hold on      
    %p_mean(iteration,:) = mean(p_nodes); %Compute velocity of COM of MSN
    q_nodes_all{iteration} = nodes;
    Connectivity(iteration)= (1/(num_nodes))*rank(A);
    
    % Plot
    plot(nodes(:,1),nodes(:,2), '.')
    hold on
    plot(nodes(:,1),nodes(:,2), 'k>','LineWidth',.2,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)
    hold off
    for node_i = 1:num_nodes
        tmp=nodes(Nei_agent{node_i},:);
        for j = 1:size(nodes(Nei_agent{node_i},1))
            line([nodes(node_i,1),tmp(j,1)],[nodes(node_i,2),tmp(j,2)]) 
        end
    end
mov(iteration) = getframe;
hold off
end  
VideoWriter(mov, 'flocking.avi', 'Compression', 'None');

% --- FUNCTIONS ---

function [Nei_agent, A] = findNeighbors(nodes, range)
    num_nodes = size(nodes,1);
    Nei_agent = cell(num_nodes,1);
    
    for i = 1:num_nodes
        for j = 1:num_nodes
           q1 = [nodes(i,1) nodes(i,2)];
           q2 = [nodes(j,1) nodes(j,2)];
           dist = norm(q1-q2);
           if i~= j && dist <= range
              Nei_agent{i} = [Nei_agent{i} j];  %Add j to list of i's neighbors
           end
        end
    end
    
    A = adjMatrix(nodes, Nei_agent);
end

function [A] = adjMatrix(nodes, Nei_agent)
%     Returns a matrix with values 0 & 1 corresponding with the adjacency of the nodes from the nodes input.
% 
%     Inputs:
%     nodes : double matrix
%     Nei_agent : cell array
%     
%     Output:
%     [A] : adjacency matrix for nodes input
    num_nodes = size(nodes, 1);
    A = zeros(num_nodes);

    for i = 1:num_nodes
       for j = 1:num_nodes
           if ismember(j, Nei_agent{i})
                A(i,j) = 1;
           end
       end
    end
end

%[Ui] = inputcontrol_Algorithm1(nodes_old,nodes,Nei_agent,n,epsilon,r,r_prime,d,k_scale,Nei_beta_agent,p_ik,q_ik,obstacles,qt1(iteration,:),pt1(iteration,:), p_nodes);
function [Ui] = inputcontrol_Algorithm1(nodes, Nei_agent, num_nodes, epsilon, r, d, p_nodes, dimensions)
    c1_alpha = 30;
    c2_alpha = 2*sqrt(c1_alpha);
    Ui = zeros(num_nodes, dimensions);
    gradient = 0.;
    consensus = 0.;
    
    % need c1*gradient + c2*consensus
    
    for i = 1:num_nodes
        for j = 1:size(Nei_agent{i},1)
            % i refers to node i
            % j refers to the jth neighbor of node i
            phi_alpha_in = sigmaNorm(nodes(Nei_agent{i}(j),:) - nodes(i,:), epsilon);
            gradient = phi_alpha(phi_alpha_in, r, d, epsilon) * nij(nodes(i,:), nodes(Nei_agent{i}(j),:), epsilon);
            consensus = aij(nodes(i,:), nodes(Nei_agent{i}(j),:), epsilon, r) * (p_nodes(j,:) - p_nodes(i,:));
        end
        Ui(i,:) = (c1_alpha * gradient) + (c2_alpha * consensus);
    end
end

function result = sigmaNorm(z, epsilon)
%     Returns the sigma norm of a given value/vector z and an epsilon contant value.
%     
%     Inputs:
%     z : double
%     epsilon : double, constant
%     
%     Output:
%     result : double
    result = (1/epsilon) * (sqrt(1 + epsilon*(norm(z))^2)-1);
end

function result = nij(i, j, epsilon)
%     Used for the gradient term in Algorithm 1 Ui.
%     
%     Inputs:
%     i : 1x2 double, position of node i
%     j : 1x2 double, position of node j
%     epsilon : double, constant
%     
%     Output:
%     result : 1x2 double
    
    numerator = j - i;
    denominator = sqrt(1 + epsilon * (norm(j-i))^2);
    result = numerator/denominator;
end

function result = bump(z)
%     A scalar function varying between 0 and 1. Used for construction of smooth potential functions with finite cut-offs and smooth adj. matrices.
%     
%     Inputs:
%     z : double
%     
%     Output:
%     result : double
    
    h = 0.2;    % Set constant h
    
    if z >= 0 && z < h
        result = 1;
    elseif z >= h && z <= 1
        result = 0.5 * (1 + cos(pi*(z-h)/(1-h)));  % CHECK - cosine in degrees or radians?
    else
        result = 0;
    end
end

function result = phi_alpha(z, r, d, epsilon)
%     Used in gradient term of Alg. 1 Ui.
%     
%     Inputs:
%     z : double
%     r : double, interaction range of nodes in MSN
%     d : double, desired distance of nodes in MSN
%     espilon : double, constant
%     
%     Output:
%     result : double
    r_alpha = sigmaNorm(r, epsilon);
    d_alpha = sigmaNorm(d, epsilon);    %CHECK - is this what d alpha is?
    result = bump(z/r_alpha) * phi(z-d_alpha);
end

function result = phi(z)
%     Used for phi_alpha function.
%     
%     Input:
%     z : double
%     
%     Output:
%     result : double

    %Set constants
    a = 5;
    b = 5;
    c = abs(a-b) / sqrt(4*a*b);
    
    sigma1 = z / sqrt(1+z^2);
    result = 0.5*((a+b)*sigma1*(z+c) + (a-b));
end

function result = aij(i, j, epsilon, r)
%     Returns the spatial adjacency matrix given the positions of two nodes, i and j.
%     
%     Inputs:
%     i : 1x2 double, position of node i
%     j : 1x2 double, position of node j
%     epsilon : constant for sigma norm
%     r : double, interaction range for nodes in MSN
%     
%     Output:
%     result : 1x2 double
    
    r_alpha = sigmaNorm(r, epsilon);
    input_to_bump = sigmaNorm(j-i, epsilon) / r_alpha;
    result = zeros(size(i));    % result is a 1x2 matrix
    result = bump(input_to_bump);
end

% function result = aij(nodes, epsilon, r, dimensions)
% %     Returns the spatial adjacency matrix given the positions of two nodes, i and j.
% %     
% %     Inputs:
% %     i : 1x2 double, position of node i
% %     j : 1x2 double, position of node j
% %     epsilon : constant for sigma norm
% %     r : double, interaction range for nodes in MSN
% %     
% %     Output:
% %     result : 1x2 double
%     
%     r_alpha = sigmaNorm(r, epsilon);
%     num_nodes = size(nodes, 1);
%     
%     result = zeros(num_nodes, num_nodes, dimensions);    % result is a 1x2 matrix
%     
%     for i = 1:num_nodes
%        for j = 1:num_nodes
%            input_to_bump = sigmaNorm(nodes(j,:)-nodes(i,:), epsilon) / r_alpha;
%            bump_val = bump(input_to_bump);
%            
%            if i ~= j
%                result(i,j,:) = bump_val;
%            end
%        end
%     end
% end
