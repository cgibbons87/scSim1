function [pStore] = scPlanets(config_file)
run(config_file);
t = 0;
h = 100;
pStore = cell(3,runTime./h);
for i = 1:(runTime./h)  
    accel = zeros(length(masses),3);
    for n = 1:length(masses)
        for m = 1:length(masses)
            %find the distance between bodies
            dist = (sum((positions(n,:)-positions(m,:)).^2)).^(1/2);
            if dist ~= 0
                %find the unit vector in the direction of the force
                uv = (positions(n,:)-positions(m,:))./(dist);
                %compute net acceleration of body n
                accel(n,:) = accel(n,:) - (6.67300.*10.^(-11)).*masses(m).*uv./(dist.^2);
            else
            end
        end
    end
    
    %store our position and velocity values
    pStore{1,i} = t;
    pStore{2,i} = positions;
    pStore{3,i} = velocities;
    
    %use Euler's method to compute change in position and update position array
    %for next iteration
    t = t + h;
    positions = positions + velocities.*h;
    velocities = velocities + accel.*h;
end