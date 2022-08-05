clc; clear all; close all;

%% compass gait walker with equations of motion of double pendulum
% all angles [theta; theta dot; phi; phi dot]
% theta = angle of stance leg with respect to the slope normal
% phi   = angle between stance and swing leg

% td = [0.02];
t = [0.2]; % t for gam = -0.009
td = [0.22];  % td for gam = -0.009
% t = [0.2]; % t for gam = -0.0145
% td = [0.25]; % td for gam = -0.0145
% t = [0.2]; % t for gam = -0.018
% td = [0.23]; % td for gam = -0.018
tdd = []; % theta double dot

% pd = [0.35];
p = [0.4]; % p for gam = -0.009
pd = [0.9]; % pd for gam = -0.009
% p = [0.3]; % p for gam = -0.0145
% pd = [0.9]; % pd for gam = -0.0145 905
% p = [0.3]; % p for gam = -0.018
% pd = [0.95]; % pd for game = -0.018
pdd = []; % phi double dot

M = 1; % hip mass
m = 1; % foot mass
l = 1; % leg length
gam = -0.009; % .11 ramp slope  (0.009, 0.0145, 0.018) should be 4 numbers
g = 9.81; % accel of grav
deltaT = 0.0001; % change in time for Eulers method
T = [0:deltaT:20]; % time 0-15 seconds

% matrix for changing velocities at healstrike
h = [-1 0 0 0; 0 0 0 0; -2 0 0 0; 0 0 0 0];

% booleans to turn off and on the search for a heelstike
boolSwitch = 0;
buffer = 0;

% need for plotting
originx = [0];
originy = [0];
origincount = 1;    % first origin is 0,0

% these will store the theta and phi that the switch happens at to use in
% the plotting loop
switchTH = [];
switchP = [];

% this will tell us how many switches happened and how many elements in
% switchTH and switchP (useful for plotting loop)
switchcount = 0;

for i=1:length(T)

    tdd(i)  = sin(t(i)-gam);

    td(i+1) = td(i) + (tdd(i) * deltaT);

    t(i+1)  =  t(i) + (td(i)  * deltaT);

    pdd(i)  = tdd(i) + (td(i)^2 * sin (p(i))) - (cos (t(i) - gam) * sin (p(i)));

    pd(i+1) = pd(i) + (pdd(i) * deltaT);

    p(i+1)  = p(i)  + (pd(i)  * deltaT);
    
    % first condition is heelstrike condition (phi = 2theta), second
    % condition is to ensure its past the vertical at the origin or base
    % of stance foot
    % third condition is to ensure it is the first point of contact with
    % the ramp
    if ((-0.01 <= p(i) - 2*t(i) && p(i) - 2*t(i) <= 0.01) || p(i) == 2*(t(i))) && ((originx(origincount) + (l * sin (t(i)))) + (l * sin (p(i)-t(i))) >= originx(origincount)+0.35) && boolSwitch == 0

        % point of contact found so turn on
        boolSwitch = 1;
        % buffer so that we dont include the other close values 
        buffer = 1;
        % found a hit so they switched
        switchcount = switchcount + 1;

        % the theta and phi that the switch happened at
        switchTH(switchcount) = t(i);
        switchP(switchcount)  = p(i);

        % with the switch found and new origin
        % add to origin x and y vectors
        origincount = origincount+1;
        originx(origincount) = (originx(origincount-1) + l * sin (t(i))) + l * sin (p(i)-t(i));
        originy(origincount) = (originy(origincount-1) + l * cos (t(i))) + l * -cos (p(i)-t(i));

        %fill in h (the matrix used in the switch equation) the rest of h is constant
        h(2,2) = cos(2*t(i));
        h(4,2) = cos(2*t(i))*(1-cos(2*t(i)));
           
        % angles just before healstrike
        aminus = [t(i); td(i); p(i); pd(i)];
 
        % calculate angles and velocities after healstike
        aplus  = h * aminus;

        % now add the new angle and velocity to the vector for each
        t(i+1)  = aplus(1);
        td(i+1) = aplus(2);
        p(i+1)  = aplus(3);
        pd(i+1) = aplus(4);

    end
    
    % determine when its far enough that the same strike wont make /\ that
    % if statement go off and turn it on again to find a new strike
    if (~((-0.05 <= p(i) - 2*t(i) && p(i) - 2*t(i) <= 0.05) || p(i) == 2*(t(i)))) && boolSwitch == 1 && buffer==1
        disp('in turn off');
        buffer = 0;
        boolSwitch = 0;
    end
end


%% start of plotting code
f1 = figure;
FPS = 20; % frames per second
t_anim = 0:1/FPS:max(T);
th_anim = interp1(T, t(1:numel(T)), t_anim); % thetas for plotting
p_anim  = interp1(T, p(1:numel(T)), t_anim); % phis for plotting
boolswitchPlot = 0;     % need so we can see if a switch happened 
switchcountPlot = 1;    % need so we can increase the theta and phi that the switch happened at
origincountPlot = 1;    % need so we can increase when the switch is found in plotting loop

% for plotting
leftX = originx(1)-1.5;
rightX = originx(1)+1.5;

for i = 1: numel(th_anim)


    clf
    %position of hip
    hipx =  originx(origincountPlot) + l * sin (th_anim(i));
    hipy =  originy(origincountPlot) + l * cos (th_anim(i));

    % position of swing foot
    swingx =  hipx + l * sin (p_anim(i) - th_anim(i));
    swingy =  hipy + l * -cos (p_anim(i)- th_anim(i));

    % thin stance leg
    plot([originx(origincountPlot) hipx], [originy(origincountPlot) hipy], 'Color', '#00132C', 'LineWidth', 3);
    hold on
    % thick swing leg
    plot([hipx swingx], [hipy swingy], 'Color', '#00132C', 'LineWidth', 5);
    plot(hipx, hipy, '.r', 'Color', '#F9CFC3', 'markersize', 30);

    yline(0, 'LineWidth',1.25, 'Color', '#00132C');
    axis equal
    axis([hipx-1.5 hipx+1.5 -1.5 1.5]); 


    % first condition is determining if there was a change from the first
    % loop (so there is no error for switchP and switchTH []
    % second condition is locating the heelstrike in the interp1 theta
    % vector
    % third condition is making sure its the first instance of it (so the
    % other close values don't make this go off)
    % fourth condition is making sure it doesn't pass how many strikes
    % there were
    if switchcount > 0 && (switchP(switchcountPlot) - 0.015 < p_anim(i) && p_anim(i) < switchP(switchcountPlot) + 0.015) && boolswitchPlot == 0 && switchcountPlot <= switchcount
        disp('here')
        boolswitchPlot = 1;
        if (switchcountPlot+1 <= switchcount)
        switchcountPlot = switchcountPlot+1;
        end
        if (origincountPlot+1 <= origincount)
        origincountPlot = origincountPlot+1;
        end
        disp(originx(origincountPlot));
        disp(originy(origincountPlot));
        disp(swingx);
        disp(swingy);
        disp('switch th is');
        disp(switchTH(switchcountPlot));
        
     
    end
    % find when to turn on boolswitchplot again 
    if boolswitchPlot == 1 && (~(switchTH(switchcountPlot) - 0.015 < th_anim(i) && th_anim(i) < switchTH(switchcountPlot) + 0.015)) 
        boolswitchPlot = 0;
    end
    drawnow
    movegui(f1, 'center');

%     if (i == 84)
%         return
%     end
    % if statement for calculating speed
    % beginning time for gam = -0.009 is 3.6, hipx = 1.4342
    % ending time for gam = -0.009 is 7.5, hipx = 2.2457
    % time difference = 7.5 - 3.6 = 3.9 s
    % x difference = 2.2457 - 1.4342 = 0.8115 ?
    % speed = 0.8115 / 3.9 = 0.20807692 = 0.2081 /s

    % beginning time for gam = -0.0145 is 4.5, hipx = 0.2152
    % ending time for gam = -0.0145 is 8.1, hipx = 0.2252
    % time difference = 8.1 - 4.5 = 3.6 s
    % x difference = 0.01
    % speed = 0.01 / 3.6 = 0.0027778 = 0.0028 /s

    % beginning time for gam = -0.018 is 4.1, hipx = 0.7490
    % ending time for gam = -0.018 is 8, hipx = 2.2814
    % time difference = 8 - 4.1 = 3.9
    % x difference = 2.2814 - 0.7490 = 1.5324
    % speed = 1.5324 / 3.9 = 0.39292308 = 0.3929 /s

end

%% plot data make graphs
f = figure;
plot (T, t(1:numel(T)));
hold on
plot (T, p(1:numel(T)));
plot (T, 2*t(1:numel(T)));
xlabel('Time (s)', 'FontName','Times New Roman', 'FontSize', 12);
ylabel('Leg Angles (Rad)', 'FontName','Times New Roman', 'FontSize', 12);
title ('Angular Position over Time', 'FontName','Times New Roman', 'FontSize', 14)
legend ('\theta', '\phi', '2\theta', 'FontSize', 12, 'FontName','Times New Roman');

movegui(f, "center");

f2 = figure;
plot (T,td(1:numel(T)));
hold on
plot (T, pd(1:numel(T)));
hold on
xlabel('Time (s)', 'FontName','Times New Roman', 'FontSize', 12);
ylabel ('Velocity (rad/s)', 'FontName','Times New Roman', 'FontSize', 12);
title ('Velocity over Time', 'FontName','Times New Roman', 'FontSize', 14)
legend ('$\dot{\theta}$', '$\dot{\phi}$', 'Interpreter', 'latex', 'FontSize', 12, 'FontName','Times New Roman')

