% Name: proj1_case5.m
% Author: Jazel A. Suguitan
% Last Modified: Oct. 8, 2021

clc,clear
close all

% CASE 5: Algorithm 3 (MSN Quasi-Lattice Formation w/ Obstacle Avoidance)

%================= SET PARAMETERS ===============

d = 15; % Set desired distance among sensor nodes - ORIGINALLY 15
k_scale = 1.2;  % Set the scale of MSN - ORIGINALLY 1.2
r = k_scale * d;  % Set the active range
r_prime = .22 * k_scale * r;    % Set the active range of beta agent
epsilon = 0.1;  % Set a constant for sigma norm
num_nodes = 100;    % Set number of sensor nodes
n = 2;  % Set number of dimensions
%nodes = load('node_distribution2.dat'); % distributed in 2D
nodes = 50.*rand(num_nodes,n)+50.*repmat([0 1],num_nodes,1);  % Randomly generate initial positions of MSN
p_nodes = zeros(num_nodes,n);   % Set initial velocties of MSN
delta_t_update = 0.008;  % Set time step - ORIGINALLY 0.008, THEN 0.04, THEN 0.0108
t = 0:delta_t_update:8; % Set simulation time

%================= SET A STATIC TARGET ===============
qt1 = [200,25]; %Set position of the static target (gamma agent)
pt1= [0,0]; % %Set initial velocity of the target

%================= SET OBSTACLES ===============
%obstacles = [100, 25];
%Rk = 15;
obstacles = [100, 150; 150 80; 200, 230; 280 150 ]; %set positions of obstacles - first ORIGINALLY 50,100
Rk = [20; 10; 15; 8]; %Radii of obstacles
num_obstacles = size(obstacles,1); %Find number of obstacles

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
%     qt_x1 = 50 + 50*t(iteration);
%     qt_y1 = 295 - 50*sin(t(iteration));
    
%   Circle Trajectory of a moving target 
%    qt_x1 = 310 - 160*cos(t(iteration));
%    qt_y1 = 255 + 160*sin(t(iteration));

%    Line Trajectory of a moving target 
    qt_x1 = 200 + 15*t(iteration); %ORIGINALLY + 130t
    qt_y1 = 200+1*t(iteration); 

%     %compute position of target 
    qt1(iteration,:) = [qt_x1, qt_y1];
%     %compute velocities of target
    if iteration >1
    pt1(iteration,:)=(qt1(iteration,:)-qt1(iteration-1,:))/delta_t_update;
    else
        continue
    end  
    plot(qt1(:,1),qt1(:,2),'ro','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r', 'MarkerSize',4.2)
    hold on
 
    [Nei_agent, Nei_beta_agent, A] = findNeighbors(nodes, r, obstacles);
    [Ui] = inputcontrol_Algorithm3(nodes, Nei_agent, num_nodes, epsilon, r, d, p_nodes, n, qt1(iteration,:), pt1(iteration,:), obstacles, Rk, Nei_beta_agent);
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
    
    %================= DRAW OBSTACLES ===============
    phiVal = 0:.1:2*pi;
    for k = 1:num_obstacles
        X = Rk(k)*cos(phiVal);
        Y = Rk(k)*sin(phiVal);
        plot(X+obstacles(k,1),Y+obstacles(k,2),'r',nodes(:,1),nodes(:,2), 'g>')
        fill(X+obstacles(k,1),Y+obstacles(k,2),'r')
        %axis([0 250 0 80]);
        hold on
    end

    
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
%      if j ==1 %Plot velociy of sensor node 1; you can change this number to plot for other nodes
%        p_each_nodes(i) =  norm(tmp7(j,:));
%     end
        p_each_nodes(i,j) =  norm(tmp7(j,:));
    end
end
figure(3), plot(p_each_nodes)%, 'b')
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

%========================PLOT TRAJECTORY OF COM AND TARGET===============                
figure(6), plot(q_mean(:,1),q_mean(:,2),'k.')
hold on
plot(qt1(:,1), qt1(:,2),'r.')


%================= FUNCTIONS ===============

function [Ui] = inputcontrol_Algorithm3(nodes, Nei_agent, num_nodes, epsilon, r, d, p_nodes, dimensions, q_mt, p_mt, obstacles, Rk, Nei_beta_agent)
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
    c1_alpha = 40;  %ORIGINALLY 30
    c2_alpha = 2*sqrt(c1_alpha);
    c1_mt = 1.7;    % ORIGINALLY 1.1
    c2_mt = 2*sqrt(c1_mt);
    c1_beta = 1500;
    c2_beta = 2*sqrt(c1_beta);
    Ui = zeros(num_nodes, dimensions);  % initialize Ui matrix to all 0's
    gradient = 0.;  % Initialize gradient part of Ui equation
    consensus = 0.; % Initialize consensus part of Ui equation
    feedback = 0.;  % Initialize navigational feedback of Ui equation
    beta_term1 = 0.;
    beta_term2 = 0.;
    
    % Sum gradient and consensus values for each node i
    for i = 1:num_nodes
        for j = 1:length(Nei_agent{i})
            % i refers to node i
            % j refers to the jth neighbor of node i
            phi_alpha_in = sigmaNorm(nodes(Nei_agent{i}(j),:) - nodes(i,:), epsilon);
            gradient = gradient + phi_alpha(phi_alpha_in, r, d, epsilon) * nij(nodes(i,:), nodes(Nei_agent{i}(j),:), epsilon);
            consensus = consensus + aij(nodes(i,:), nodes(Nei_agent{i}(j),:), epsilon, r) * (p_nodes(Nei_agent{i}(j),:) - p_nodes(i,:));
        end
        
        for k = 1:length(Nei_beta_agent{i})
            % k refers to the kth beta agent of node i
            phi_beta_in = sigmaNorm(qik(nodes(i,:), obstacles(Nei_beta_agent{i}(k),:), Rk(Nei_beta_agent{i}(k)))-nodes(i,:), epsilon);
            beta_term1 = beta_term1 + phi_beta(phi_beta_in, d, epsilon) * nik(nodes(i,:), obstacles(Nei_beta_agent{i}(k),:), epsilon, Rk(Nei_beta_agent{i}(k)));
            beta_term2 = beta_term2 + bik(nodes(i,:), obstacles(Nei_beta_agent{i}(k),:), d, Rk(Nei_beta_agent{i}(k)), epsilon) * (pik(nodes(i,:), obstacles(Nei_beta_agent{i}(k),:), Rk(Nei_beta_agent{i}(k)), p_nodes(i,:)) - p_nodes(i,:));
        end
        
        feedback = -(c1_mt * (nodes(i,:) - q_mt)) - (c2_mt * (p_nodes(i,:) - p_mt));
        beta_total = (c1_beta * beta_term1) + (c2_beta * beta_term2);
        Ui(i,:) = (c1_alpha * gradient) + (c2_alpha * consensus) + feedback + beta_total;   % Set Ui for node i using gradient, consensus, and feedback
        gradient = 0;
        consensus = 0;
        beta_term1 = 0;
        beta_term2 = 0;
        beta_total = 0;
        feedback = 0;
    end
end

function [Nei_agent, Nei_beta_agent, A] = findNeighbors(nodes, range, obstacles)
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
    %INITIAL BETA AGENT INCORPORATION
    Nei_beta_agent = cell(num_nodes, 1);
    num_obstacles = size(obstacles, 1); %EDIT - need to pass in obstacles to findNeighbors
    
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
    
    % Iterate through each node i for Nei_beta_agent
    for i = 1:num_nodes
        for j = 1:num_obstacles
           % Check if obstacle j is in the interaction range of node i
           beta_pos = qik(nodes(i,:), obstacles(j,:), range);
           if norm(beta_pos - nodes(i,:)) <= range
              Nei_beta_agent{i} = [Nei_beta_agent{i} j];
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
    
    result = sigmaE(j-i, epsilon);
end

function result = sigmaE(z, epsilon)
%     Function to be used in nij function.
%     
%     Parameters
%     -------------
%     z : double (1x2)
%           A vector
%     epsilon : double
%           A constant
%     
%     Returns
%     -----------
%     result : double (1x2)
%           The resulting vector

    result = z / (1 + epsilon * sigmaNorm(z, epsilon));
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
    d_alpha = sigmaNorm(d, epsilon);
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
    
    sigmaZ = sigma1(z+c);
    result = 0.5*((a+b)*sigmaZ + (a-b));
end

function result = sigma1(z)
%     A function to be used in the phi function.
%     
%     Parameters
%     -------------
%     z : double
%     
%     Returns
%     -----------
%     result : double

    result = z / sqrt(1+z^2);
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
    
    result = zeros(size(i));    % result is a 1x2 matrix - CHECK
    
    if ~isequal(i,j)
        r_alpha = sigmaNorm(r, epsilon);
        input_to_bump = sigmaNorm(j-i, epsilon) / r_alpha;
        result = bump(input_to_bump);
    end
end

%=======================START OF CASE 5 FUNCTIONS==============================

function result = nik(i, k, epsilon, radius)
%     Function to be used in the beta term of Ui.
%     
%     Parameters
%     ------------
%     i : double (1x2)
%         The position of node i
%     k : double (1x2)
%         The position of obstacle k
%     epsilon : double
%         A constant for the sigma norm
%     radius : double
%         The radius of obstacle k
%         
%     Returns
%     ------------
%     result : double (1x2)
%         n_i,k

    q_ik = qik(i, k, radius);
    
    numerator = q_ik - i;
    denominator = sqrt(1 + epsilon * (norm(q_ik-i))^2);
    
    result = numerator/denominator;
end

function result = qik(i, k, radius)
%     Function for the position of a beta agent for node i to obstacle k.
%     
%     Parameters
%     ------------
%     i : double (1x2)
%         The position of node i
%     k : double (1x2)
%         The position of obstacle k
%     radius : double
%         The radius of obstacle k
%         
%     Returns
%     ------------
%     result : double (1x)
%         The position of the beta agent.

    mu1 = mu(i, k, radius);
    result = mu1*i + (1-mu1)*k;
end

function result = mu(i, k, radius)
%     Function for obtaining mu, to be used in qik and pik.
%     
%     Parameters
%     ------------
%     i : double (1x2)
%         The position of node i
%     k : double (1x2)
%         The position of obstacle k
%     radius : double
%         The radius of obstacle k
%         
%     Returns
%     ------------
%     result : double
%         mu

    result = radius/norm(i-k);
end

function result = bik(i, k, d, radius, epsilon)
%     Function for obtaining the smoothening of the beta Ui term.
%     
%     Parameters
%     ------------
%     i : double (1x2)
%         The position of node i
%     k : double (1x2)
%         The position of obstacle k
%     d : double
%         The desired distance among nodes in the MSN
%     radius : double
%         The radius of obstacle k
%     epsilon : double
%         A constant for the sigma norm
%         
%     Returns
%     ------------
%     result : double
%         b_i,k(q)

    q_ik = qik(i, k, radius);
    d_beta = sigmaNorm(d, epsilon);
    sig = sigmaNorm(q_ik - i, epsilon);
    
    result = bump(sig/d_beta);
end

function result = pik(i, k, radius, p_i)
%     Function for obtaining the velocitiy of the beta agent for node i to obstacle k.
%     
%     Parameters
%     ------------
%     i : double (1x2)
%         The position of node i
%     k : double (1x2)
%         The position of obstacle k
%     radius : double
%         The radius of obstacle k
%     p_i : double (1x2)
%         The velocity of node i
%         
%     Returns
%     -------------
%     result : double (1x2)
%         The velocity of the beta agent
    
    result = mu(i, k, radius) * p_mat(i, k) * p_i;
end

function result = p_mat(i, k)
%     Function for obtaining P, to be used in the pik function.
%     
%     Parameters
%     ------------
%     i : double (1x2)
%         The position of node i
%     k : double (1x2)
%         The position of obstacle k
%         
%     Returns
%     ------------
%     result : double matrix
%         The P matrix to be used in finding a beta agent's velocity

    a_k = ak(i,k);
    ak_akT = a_k * transpose(a_k);
    I = eye(size(ak_akT,1), size(ak_akT,2));
    result = I - ak_akT;
end

function result = ak(i, k)
%     Function for getting the unit normal.
%     
%     Parameters
%     ------------
%     i : double (1x2)
%         The position of node i
%     k : double (1x2)
%         The position of obstacle k
%         
%     Returns
%     ------------
%     result : double matrix
%         The unit normal
    
    result = (i-k)/norm(i-k);
end

function result = phi_beta(z, d, epsilon)
%     The action function used to construct a smooth pairwise potential
%     with finite cut-off in the beta term of the Alg.3 Ui.
%     
%     Parameters
%     -------------
%     z : double
%           Sigma norm value of two nodes
%     d : double
%           Desired distance of nodes in MSN
%     espilon : double
%           A constant for the sigma norm
%     
%     Returns
%     -----------
%     result : double
%           Value to be used in Ui beta term

    d_beta = sigmaNorm(d, epsilon);
    
    result = bump(z/d_beta) * (sigma1(z-d_beta) - 1);
end