
clear all; clc;

rng("default")


% field = 1;
for field = 1:11
    field
for n = 1:300
    n
    clear x0 y0 f
    for i = 1:10
        x0(i) = -10 + 20 * rand();
        y0(i) = -10 + 20 * rand();
    end
    x0 = x0';
    y0 = y0';
    % Define the time span for the simulation
    % tspan = [0 5+20*rand()];
    tspan = linspace(0,30*rand(),10);
    
    % Define the initial positions of the particles
    % For example, let's start 10 particles on a line at y = 1
    % y0 = ones(10, 1);
    % x0 = linspace(-10, 10, 10)';
    initial_conditions = [x0 y0];
    
    
    f = figure('visible', 'off');
    % f = figure;
    % Run the ODE solver for each particle
    for i = 1:size(initial_conditions, 1)
        % Solve the ODE using ode45
        switch field
            case 1
            [t, Y] = ode45(@(t, y) spiral_sink_field(t, y), tspan, initial_conditions(i, :));
            case 2
            [t, Y] = ode45(@(t, y) spiral_source_field(t, y), tspan, initial_conditions(i, :));   
            case 3
            [t, Y] = ode45(@(t, y) center_velocity_field(t, y), tspan, initial_conditions(i, :));
            case 4
            [t, Y] = ode45(@(t, y) saddle_field(t, y), tspan, initial_conditions(i, :));  
            case 5
            [t, Y] = ode45(@(t, y) shear_field(t, y), tspan, initial_conditions(i, :));   
            case 6
            [t, Y] = ode45(@(t, y) focus_converging_field(t, y), tspan, initial_conditions(i, :));  
            case 7
            [t, Y] = ode45(@(t, y) focus_diverging_field(t, y), tspan, initial_conditions(i, :)); 
            case 8
            [t, Y] = ode45(@(t, y) node_converging_field(t, y), tspan, initial_conditions(i, :)); 
            case 9
            [t, Y] = ode45(@(t, y) node_diverging_field(t, y), tspan, initial_conditions(i, :));  
            case 10
            [t, Y] = ode45(@(t, y) improper_node_convergent(t, y), tspan, initial_conditions(i, :)); 
            case 11
            [t, Y] = ode45(@(t, y) improper_node_divergent(t, y), tspan, initial_conditions(i, :)); 
            otherwise
            % costmize
            [t, Y] = ode45(@(t, y) node_converging_field_plus_shear(t, y), tspan, initial_conditions(i, :)); 
        end
        % Plot the trajectory of each particle
        plot(Y(:, 1), Y(:, 2),'k','LineWidth',2);
        hold on;
        plot(Y(end,1), Y(end,2), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'k','MarkerEdgeColor','k');
    end
    
    % Set the axes to be equal and labels
    axis equal;
    axis off ;
    % xlabel('X');
    % ylabel('Y');
    % title('Particle Trajectories in a Center Velocity Field');
    % legend('Trajectories');
    switch field
        case 1
        folderName = 'spiral_sink_field'; % Specify the name of the folder you want to check/create
        % Check if the folder exists
        if ~exist(folderName, 'dir')
        % Folder does not exist so create it
            mkdir(folderName);
        end
        % filename = sprintf('./spiral_sink_field/sample_%d.png', n);
        filename = sprintf('./spiral_sink_sample_%d.png', n);
        case 2
        folderName = 'spiral_source_field'; % Specify the name of the folder you want to check/create
        % Check if the folder exists
        if ~exist(folderName, 'dir')
        % Folder does not exist so create it
            mkdir(folderName);
        end
        % filename = sprintf('./spiral_source_field/sample_%d.png', n);
        filename = sprintf('./spiral_source_sample_%d.png', n);
        case 3
        folderName = 'center_velocity_field'; % Specify the name of the folder you want to check/create
        % Check if the folder exists
        if ~exist(folderName, 'dir')
        % Folder does not exist so create it
            mkdir(folderName);
        end
        % filename = sprintf('./center_velocity_field/sample_%d.png', n);
        filename = sprintf('./center_velocity_sample_%d.png', n);
        case 4
        folderName = 'saddle_field'; % Specify the name of the folder you want to check/create
        % Check if the folder exists
        if ~exist(folderName, 'dir')
        % Folder does not exist so create it
            mkdir(folderName);
        end
        % filename = sprintf('./saddle_field/sample_%d.png', n);
        filename = sprintf('./saddle_sample_%d.png', n);
        case 5
        folderName = 'shear_field'; % Specify the name of the folder you want to check/create
        % Check if the folder exists
        if ~exist(folderName, 'dir')
        % Folder does not exist so create it
            mkdir(folderName);
        end
        % filename = sprintf('./shear_field/sample_%d.png', n);
        filename = sprintf('./shear_sample_%d.png', n);
        case 6
        folderName = 'focus_converging_field'; % Specify the name of the folder you want to check/create
        % Check if the folder exists
        if ~exist(folderName, 'dir')
        % Folder does not exist so create it
            mkdir(folderName);
        end
        % filename = sprintf('./focus_converging_field/sample_%d.png', n);
        filename = sprintf('./focus_converging_sample_%d.png', n);
        case 7
        folderName = 'focus_diverging_field'; % Specify the name of the folder you want to check/create
        % Check if the folder exists
        if ~exist(folderName, 'dir')
        % Folder does not exist so create it
            mkdir(folderName);
        end
        % filename = sprintf('./focus_diverging_field/sample_%d.png', n);
        filename = sprintf('./focus_diverging_sample_%d.png', n);
        case 8
        folderName = 'node_converging_field'; % Specify the name of the folder you want to check/create
        % Check if the folder exists
        if ~exist(folderName, 'dir')
        % Folder does not exist so create it
            mkdir(folderName);
        end
        % filename = sprintf('./node_converging_field/sample_%d.png', n);
        filename = sprintf('./node_converging_sample_%d.png', n);
        case 9
        folderName = 'node_diverging_field'; % Specify the name of the folder you want to check/create
        % Check if the folder exists
        if ~exist(folderName, 'dir')
        % Folder does not exist so create it
            mkdir(folderName);
        end
        % filename = sprintf('./node_diverging_field/sample_%d.png', n);
        filename = sprintf('./node_diverging_sample_%d.png', n);
        case 10
        folderName = 'improper_node_convergent_field'; % Specify the name of the folder you want to check/create
        % Check if the folder exists
        if ~exist(folderName, 'dir')
        % Folder does not exist so create it
            mkdir(folderName);
        end
        % filename = sprintf('./node_diverging_field/sample_%d.png', n);
        filename = sprintf('./improper_node_convergent_sample_%d.png', n);
        case 11
        folderName = 'improper_node_divergent_field'; % Specify the name of the folder you want to check/create
        % Check if the folder exists
        if ~exist(folderName, 'dir')
        % Folder does not exist so create it
            mkdir(folderName);
        end
        % filename = sprintf('./node_diverging_field/sample_%d.png', n);
        filename = sprintf('./improper_node_divergent_sample_%d.png', n);

        otherwise
        filename = sprintf('./sample_%d.png', n);

    end 

    saveas(gcf, filename);

end


end



function dYdt = spiral_sink_field(t, Y)
    % A spiral sink can be represented by a system such as:
    % dx/dt = a*x - b*y
    % dy/dt = b*x + a*y
    % where 'a' is negative (causing the spiral inwards) and 'b' causes rotation
    % a = -0.15;
    % b = -0.15;
    a = -0.15*2*rand()*2*rand();
    b = -0.15*2*rand()*2*rand();
    dYdt = [a*Y(1) - b*Y(2); b*Y(1) + a*Y(2)];
end

function dYdt = spiral_source_field(t, Y)
    % A spiral source can be represented by a system such as:
    % dx/dt = a*x + b*y
    % dy/dt = -b*x + a*y
    % where 'a' is positive (indicating source) and 'b' creates rotation
    a = 0.15*2*rand()*2*rand();
    b = 0.15*2*rand()*2*rand();
    dYdt = [a*Y(1) + b*Y(2); -b*Y(1) + a*Y(2)];
end


% Velocity field function
function dydt = center_velocity_field(~, Y)
    % This function returns the velocity at a point y
    % For a center, the velocity field is given by:
    % dx/dt = -y
    % dy/dt = x
    a = 0.15*2*rand()*2*rand();
    b = 0.15*2*rand()*2*rand();
    dydt = [-a*Y(2); b*Y(1)];
end


% Velocity field function for a saddle
function dYdt = saddle_field(t, Y)
    % A saddle can be represented by a system such as:
    % dx/dt = a*x
    % dy/dt = -b*y
    % where 'a' is positive (repelling along x) and 'b' is positive (attracting along y)
    a = -0.05*2*rand()*2*rand();  % repelling rate along x
    b = -0.05*2*rand()*2*rand();  % attracting rate along y
    dYdt = [a*Y(1); -b*Y(2)];
end

% Velocity field function for shear
function dYdt = shear_field(t, Y)
    % Define the shear rate constant
    k = 0.1*2*rand()*2*rand();
    % Velocity field equations
    y=Y(2);
    dYdt = [k*y; 0];
end


function dYdt = focus_converging_field(t, Y)
    % u = -x;
    % v = -y;
    a=0.15*2*rand()*2*rand();
    b=0.15*2*rand()*2*rand();
    dYdt = [-a*Y(1); -b*Y(2)];
end

function dYdt = focus_diverging_field(t, Y)
    % u = x;
    % v = y;
    a=0.15*2*rand()*2*rand();
    b=0.15*2*rand()*2*rand();
    dYdt = [a*Y(1); b*Y(2)];
end

function dYdt = node_converging_field(t, Y)
    x = Y(1);
    y = Y(2);
    r = sqrt(x^2 + y^2);
    % u = -5*x / r;
    % v = -y / r;
    a=0.3*2*rand();
    b=0.3*2*rand();
    dYdt = [-5*a*Y(1)/r; -b*Y(2)/r];
end

function dYdt = node_diverging_field(t, Y)
    x = Y(1);
    y = Y(2);
    r = sqrt(x^2 + y^2);
    % u = 5*x / (r);
    % v = y / (r);
    a=0.15*2*rand();
    b=0.15*2*rand();
    dYdt = [5*a*Y(1)/r; b*Y(2)/r];
end

function dYdt = improper_node_convergent(t, Y)
    % Improper node (convergent) velocity field
    % dx/dt = a*x + b*y
    % dy/dt = c*x + d*y
    % Here a and d are negative for convergence, b and c determine the asymmetry

    a = -0.5*2*rand();  % Negative for convergence
    b = 0.4*2*rand();   % Skewing effect
    c = -0.2*2*rand();  % Skewing effect
    d = -0.6*2*rand();  % Negative for convergence

    dYdt = [a*Y(1) + b*Y(2); c*Y(1) + d*Y(2)];
end


function dYdt = improper_node_divergent(t, Y)
    % Improper node (divergent) velocity field
    % dx/dt = a*x + b*y
    % dy/dt = c*x + d*y
    % Here a and d are positive for divergence, b and c determine the asymmetry

    a = 0.5*2*rand();   % Positive for divergence
    b = -0.4*2*rand();  % Skewing effect
    c = 0.2*2*rand();   % Skewing effect
    d = 0.6*2*rand();   % Positive for divergence

    dYdt = [a*Y(1) + b*Y(2); c*Y(1) + d*Y(2)];
end





% function dYdt = node_converging_field_plus_shear(t, Y)
%     x = Y(1);
%     y = Y(2);
%     r = sqrt(x^2 + y^2);
%     k = 0.05*2*rand()*2*rand();
%     % u = 5*x / (r);
%     % v = y / (r);
%     a=0.15*2*rand()*2*rand();
%     b=0.15*2*rand()*2*rand();
%     dYdt = [-5*a*Y(1)/r + k*Y(2); -b*Y(2)/r];
% end
