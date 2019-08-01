%simulation run time
runTime = 100000;

%Outline of the planetary system.
masses     = [5.9736.*10.^24, 7.3477.*10.^22];    % ['mass of body 1', 'mass of body 2', . . .]
names      = {['Earth'],['Moon']};                % ['name of body 1', 'name of body 2', . . .]
radii      = [6371000, 1737100];                  % ['radius of body 1', 'radius of body 2', . . .]
atmP       = [101.325, 0];                        %surface P, 0 to ignore drag.
positions  = [0, -5000000, 0; 0, 385000000, 0 ];  %initial ['X1','Y1','Z1'; 'X2','Y2','Z2', . . .]
velocities = [12.5, 0, 0; -1022, 0, 0];           %initial ['Vx1','Vy1','Vz1'; 'Vx2','Vy2','Vz2', . . .]
rotations  = [86400, 0, 0; 2592000, 0, 0];
%Spacecraft
scP = [0, (6371000-5000000), 0];
scV = [0, 0, 0];
scDD = [12.56636, 0.5]; %[area, drag coeff] command module had Cd of 0.5
stages   = [12500, 12500, 3200, 0;
        5000, 35000, 4400, 600;
        15000, 350000, 3100, 201];  %['stage fuel mass', 'stage dry mass', 'stage exhaust V', 'staging time'; . . . ]
scOnGround = [0, pi./2, 0];

%this is where things get tricky. The first column is the iteration number 
%of the beginning of the the burn. The second is the end iteration number.
%The third is the deltaV expended per iteration. The fourth Columns 5, 6 
%and 7 describe the equation of the angle(in radians) of the burn
%(theta = Ai^2 + Bi + C, where i is the number of burn iterations that have
%already been completed)(remember, "down" is the direction of gravitational
%force. so a positive thrust pushes you AWAY from the body.
%Burns MUST be at least two iterations apart.
burns    = [4.9, 199, 5850000, 1, -pi./2, -pi./160, 0;%pi./2, pi./160
            202, 349.35, 1002100, 1, -pi./2, -pi./160, 0;
            2164, 2177.5, 100000, 0, 0, 0, 0];


%[1=apogee, body numb, start time range, end time range, engine thrust]
apAutoBurns = [1, 1, 2800, 3200, 1002100];