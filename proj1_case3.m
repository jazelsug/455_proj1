% Name: proj1_case3.m
% Author: Jazel A. Suguitan
% Last Modified: Oct. 1, 2021

clc,clear
close all

% CASE 3: Algorithm 2 (MSN Quasi-Lattice Formation w/ Dynamic Target)

%================= SET PARAMETERS ===============

d = 5; % Set desired distance among sensor nodes - ORIGINALLY 15
k_scale = 1.2;  % Set the scale of MSN - ORIGINALLY 1.2
r = k_scale * d;  % Set the active range
r_prime = .22 * k_scale * r;    % Set the active range of beta agent
epsilon = 0.1;  % Set a constant for sigma norm
num_nodes = 100;    % Set number of sensor nodes
n = 2;  % Set number of dimensions
%nodes = load('node_distribution2.dat'); % distributed in 2D
nodes = 150.*rand(num_nodes,n)+150.*repmat([0 1],num_nodes,1);  % Randomly generate initial positions of MSN
p_nodes = zeros(num_nodes,n);   % Set initial velocties of MSN
delta_t_update = 0.0108;  % Set time step - ORIGINALLY 0.008, THEN 0.04
t = 0:delta_t_update:7; % Set simulation time

%================= SET A STATIC TARGET ===============
qt1 = [150,150]; %Set position of the static target (gamma agent)
pt1= [0,0]; % %Set initial velocity of the target

%================= SET NODES AND CHECKS DURING ITERATIONS ===============

nodes_old = nodes; %KEEP privious positions of MSN
q_mean = zeros(size(t,2),n);%Save positions of COM (Center of Mass)
p_mean = zeros(size(t,2),n);%Save velocities of COM (Center of Mass)
Connectivity = []; %save connectivity of MSN
q_nodes_all = cell(size(t,2),num_nodes);
p_nodes_all = cell(size(t,2),num_nodes);
nFrames = 20; %set number of frames for the movie
mov(1:nFrames) = struct('cdata', [],'colormap', []); % Preallocate movie structure.

%================= START ITERATION ===============

for iteration =1:length(t)
%   Sinewave Trajectory of a moving target 
    qt_x1 = 50 + 50*t(iteration);
    qt_y1 = 295 - 50*sin(t(iteration));
    
%   Circle Trajectory of a moving target 
%    qt_x1 = 310 - 160*cos(t(iteration));
%    qt_y1 = 255 + 160*sin(t(iteration));
%    Line Trajectory of a moving target 
%     qt_x1 = 200 + 130*t(iteration);
%     qt_y1 = 200+1*t(iteration); 

    %compute position of target 
    qt1(iteration,:) = [qt_x1, qt_y1];
    %compute velocities of target
    if iteration >1
    pt1(iteration,:)=(qt1(iteration,:)-qt1(iteration-1,:))/delta_t_update;
    else
        continue
    end  
    plot(qt1(:,1),qt1(:,2),'ro','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r', 'MarkerSize',4.2)
    hold on
    
%     [Nei_agent, Nei_beta_agent, p_ik, q_ik, A] = findneighbors(nodes_old,nodes,r, r_prime,obstacles, Rk,n, p_nodes,delta_t_update);
%     [Ui] = inputcontrol_Algorithm2(nodes_old,nodes,Nei_agent,n,epsilon,r,r_prime,d,k_scale,Nei_beta_agent,p_ik,q_ik,obstacles,qt1(iteration,:),pt1(iteration,:), p_nodes);
    
    [Nei_agent, A] = findNeighbors(nodes, r);
    [Ui] = inputcontrol_Algorithm2(nodes, Nei_agent, num_nodes, epsilon, r, d, p_nodes, n, qt1(iteration,:), pt1(iteration,:));
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
    
    %================= PLOT and LINK SENSOR TOGETHER ===============
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

%======================== MAKE AVI ===========================
% % used to be movie2avi
% v = VideoWriter('flocking.avi');
% open(v)
% writeVideo(v, mov)
% close(v)

%======================== PLOT VELOCITY OF MSN ===========================
p_each_nodes = [];
for i = 2:size(t,2)                    
    tmp7 = p_nodes_all{i};
    for j = 1:num_nodes
     if j ==1 %Plot velociy of sensor node 1; you can change this number to plot for other nodes
       p_each_nodes(i) =  norm(tmp7(j,:));
    end
    end
end
figure(3), plot(p_each_nodes, 'b')
hold on
figure(4),plot(Connectivity)
grid on
%======================= PLOT TRAJECTORY OF SENSOR NODES ===============
for i = 2:length(q_nodes_all)
    tmp8 = q_nodes_all{i};
    figure(5), plot(tmp8(:,1), tmp8(:,2), 'k.')
    hold on
end
hold on
plot(nodes(:,1),nodes(:,2), 'm>','LineWidth',.2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5)


%================= FUNCTIONS ===============

function [Ui] = inputcontrol_Algorithm2(nodes, Nei_agent, num_nodes, epsilon, r, d, p_nodes, dimensions, q_mt, p_mt)
%     Function for generating the Ui controller of the MSN.
%     
%     Parameters
%     -------------
%     nodes : double matrix (100x2)
%         Matrix of node positions in x-y coordinates
%     Nei_agent : cell array (100x1)
%         A container holding the neighbor indices for each node
%     num_nodes : double
%         The number of nodes in the MSN
%     epsilon : double
%         A constant for sigma norm
%     r : double
%         The interaction range of the nodes in the MSN
%     d : double
%         The desired distance among nodes in the MSN
%     p_nodes : double matrix (100x2)
%         The velocities of nodes, given in x and y directions
%     dimensions : double
%         The number of dimensions in which the MSN is operating
%     q_mt : double array (1x2)
%         The position of the moving target.
%     p_mt : double array (1x2)
%         The velocity of the moving target.
%         
%     Returns
%     -------------
%     [Ui] : double matrix (100x2)
%         Controls the positions of the nodes in the MSN as time progresses

    % Set constants
    c1_alpha = 30;
    c2_alpha = 2*sqrt(c1_alpha);
    c1_mt = 1.1;    % ORIGINALLY 1.1
    c2_mt = 2*sqrt(c1_mt);
    Ui = zeros(num_nodes, dimensions);  % initialize Ui matrix to all 0's
    gradient = 0.;  % Initialize gradient part of Ui equation
    consensus = 0.; % Initialize consensus part of Ui equation
    feedback = 0.;  % Initialize navigational feedback of Ui equation
    
    % Sum gradient and consensus values for each node i
    for i = 1:num_nodes
        for j = 1:size(Nei_agent{i},1)
            % i refers to node i
            % j refers to the jth neighbor of node i
            phi_alpha_in = sigmaNorm(nodes(Nei_agent{i}(j),:) - nodes(i,:), epsilon);
            gradient = phi_alpha(phi_alpha_in, r, d, epsilon) * nij(nodes(i,:), nodes(Nei_agent{i}(j),:), epsilon);
            consensus = aij(nodes(i,:), nodes(Nei_agent{i}(j),:), epsilon, r) * (p_nodes(j,:) - p_nodes(i,:));
        end
        feedback = -(c1_mt * (nodes(i,:) - q_mt)) - (c2_mt * (p_nodes(i,:) - p_mt));
        Ui(i,:) = (c1_alpha * gradient) + (c2_alpha * consensus) + feedback;   % Set Ui for node i using gradient, consensus, and feedback
    end
end

function [Nei_agent, A] = findNeighbors(nodes, range)
%     Function for determining the neighbors of each node in a collection of nodes.
%     
%     Parameters
%     ------------
%     nodes : double matrix (100x2)
%         Matrix of node positions in x-y coordinates
%     range : double
%         The interaction range of the nodes in the MSN
%         
%     Returns
%     ------------
%     Nei_agent : cell array (100x1)
%         A container holding the neighbor indices for each node
%     A : double matrix (100x100)
%         The adjacency matrix of nodes

    num_nodes = size(nodes, 1);
    Nei_agent = cell(num_nodes, 1);  % Initialize cell array to hold indices of neighbors
    
    % Iterate through each node i
    for i = 1:num_nodes
        for j = 1:num_nodes
           % Check each node j if it's a neighbor of node i
           q1 = [nodes(i,1) nodes(i,2)];    % Set q1 with node i values
           q2 = [nodes(j,1) nodes(j,2)];    % Set q2 with node j values
           dist = norm(q1-q2);  % Euclidean norm of q1 and q2
           if i~= j && dist <= range && dist ~= 0
              Nei_agent{i} = [Nei_agent{i} j];  %Add j to list of i's neighbors
           end
        end
    end
    
    A = adjMatrix(nodes, Nei_agent); % Use adjMatrix function to obtain adjacency matrix
end

function [A] = adjMatrix(nodes, Nei_agent)
%     Function for obtaining the adjacency matrix for a set of nodes. Used
%     for calculating Connectivity in the MSN.
% 
%     Parameters
%     -------------
%     nodes : double matrix (100x2)
%           Matrix of node positions in x-y coordinates
%     Nei_agent : cell array
%           A container holding the neighbor indices for each node
%     
%     Returns
%     -------------
%     [A] : double matrix (100x100)
%           The adjacency matrix of nodes

    num_nodes = size(nodes, 1);
    A = zeros(num_nodes);   % Initialize matrix with 0s

    for i = 1:num_nodes
       for j = 1:num_nodes
           if ismember(j, Nei_agent{i})
               % Node i and node j are neighbors
                A(i,j) = 1; % Set value in matrix to 1
           end
       end
    end
end

function result = sigmaNorm(z, epsilon)
%     Returns the sigma norm of a given value/vector z and an epsilon contant value.
%     
%     Parameters
%     -------------
%     z : double
%           Vector of which to take the sigma norm
%     epsilon : double
%           A constant used to calculate the sigma norm
%     
%     Returns
%     -----------
%     result : double
%           The sigma norm value of vector z

    result = (1/epsilon) * (sqrt(1 + epsilon*(norm(z))^2)-1);
end

function result = nij(i, j, epsilon)
%     Function for obstaining the vector along the line connecting two nodes.
%     Used in calculating the gradient-based term in Algorithm 1 Ui.
%     
%     Parameters
%     -------------
%     i : double array (1x2)
%           Position of node i
%     j : double array (1x2)
%           Position of node j
%     epsilon : double
%           A constant
%     
%     Returns
%     -----------
%     result : double array (1x2)
%           The vector along the line connecting node i and node j
    
    numerator = j - i;
    denominator = sqrt(1 + epsilon * (norm(j-i))^2);
    result = numerator/denominator;
end

function result = bump(z)
%     A scalar function varying between 0 and 1. Used for construction of smooth potential functions with finite cut-offs and smooth adj. matrices.
%     
%     Parameters
%     -------------
%     z : double
%           The input to be smoothened
%     
%     Returns
%     -----------
%     result : double
%           The 0 or 1 value
    
    h = 0.2;    % Set constant h
    
    if z >= 0 && z < h
        result = 1;
    elseif z >= h && z <= 1
        result = 0.5 * (1 + cos(pi*(z-h)/(1-h)));
    else
        result = 0;
    end
end

function result = phi_alpha(z, r, d, epsilon)
%     The action function used to construct a smooth pairwise potential
%     with finite cut-off in the gradient-based term of the Alg.1 Ui.
%     
%     Parameters
%     -------------
%     z : double
%           Sigma norm value of two nodes
%     r : double
%           Interaction range of nodes in MSN
%     d : double
%           Desired distance of nodes in MSN
%     espilon : double
%           A constant for the sigma norm
%     
%     Returns
%     -----------
%     result : double
%           Value to be used in Ui gradient-based term

    r_alpha = sigmaNorm(r, epsilon);
    d_alpha = sigmaNorm(d, epsilon);    %CHECK - is this what d alpha is?
    result = bump(z/r_alpha) * phi(z-d_alpha);
end

function result = phi(z)
%     An uneven sigmoidal function, used in the phi_alpha function.
%     
%     Parameters
%     -------------
%     z : double
%     
%     Returns
%     -----------
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
%     Parameters
%     -------------
%     i : double array (1x2)
%           Position of node i
%     j : double array (1x2)
%           Position of node j
%     epsilon : double
%           Constant for sigma norm
%     r : double
%           Interaction range for nodes in MSN
%     
%     Returns
%     -----------
%     result : double array (1x2)
%           The spatial adjacency matrix
    
    r_alpha = sigmaNorm(r, epsilon);
    input_to_bump = sigmaNorm(j-i, epsilon) / r_alpha;
    result = zeros(size(i));    % result is a 1x2 matrix
    result = bump(input_to_bump);
end
