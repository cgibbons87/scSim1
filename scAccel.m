function [out alt] = scAccel(t, scVelocity, scPosition, isFuel)
global cf;
global A;
global store;
global int_plan;
global g_alt;

run(cf);


%use linear interpolation to find the postion and velocity for each planet
%at any arbitrary time, not just at an iteration number
%A = cell2mat(store(1,:));
[c d] = max(A(A<=t));

if (d+1) < length(A) && (A(d)~=A(d+1))
    L0 = (t-A(d));
    L1 = (t-A(d+1));
    L2 = (t-A(d+2));
    L01 = A(d)-A(d+1);
    L02 = A(d)-A(d+2);
    L12 = A(d+1)-A(d+2);
    W1 = L1.*L2./(L01.*L02);
    W2 = -L0.*L2./(L01.*L12);
    W3 = L0.*L1./(L02.*L12);
    pPositions = W1.*(store{2, d}) + W2.*(store{2, (d+1)}) + W3.*(store{2, (d+2)});
    pVelocities = W1.*(store{3, d}) + W2.*(store{3, (d+1)}) + W3.*(store{3, (d+2)});
else
    pPositions = store{2, d};
    pVelocities = store{3, d};
end


%acceleration due to gravity
gAccel = 0;
for i = 1:length(masses)
    %find the distance between spacecraft and the system's masses
    dist(i) = (sum((scPosition-pPositions(i,:)).^2)).^(1/2);
    if dist(i) ~= 0
        %find the unit vector in the direction of the force
        uv = (scPosition-pPositions(i,:))./(dist(i));
        %compute net acceleration of body n
        gAccel = gAccel - (6.67300e-11).*masses(i).*uv./(dist(i).^2);
        %gAccel(isinf(gAccel)) = 0;
    else
    end
end

%acceleration due to drag
[a b] = min(dist);
alt = (dist(b)-radii(b))./1000;
dAccel = 0;
if (atmP(b) ~= 0)&&(alt<100)
    %calculate the geopotential altitude
    gpAlt = radii(b) - (radii(b)^2)./dist(b);
    H = 6700; %atmopheric scale height
    rhoO = 1.752;
    rho = rhoO.*exp(-gpAlt./H);
    atm_uv = (pPositions(b,:)-scPosition)./((sum(((pPositions(b,:)-scPosition).^2)).^(1/2)));
    atm_uv_r = atm_uv*[cos(pi./2), -sin(pi./2), 0; sin(pi./2), cos(pi./2), 0; 0, 0, 1];
    atm_vel = (2.*pi.*radii(b))./rotations(b).*atm_uv_r;
    relVel = scVelocity - pVelocities(b,:) - atm_vel;
    uv = (relVel./((sum(relVel.^2)).^(1/2)));
    dAccel = -uv.*0.5.*rho.*scDD(1).*scDD(2).*(sum(relVel.^2))./sum(sum(stages(:,[1,2])));
    dAccel(isnan(dAccel)) = 0;
    dAccel(isinf(dAccel)) = 0;
end

%acceleration from thrusting
tAccel = 0;
if (any((burns(:,1) <= t) & (burns(:,2) >= t)))&&(isFuel)
    %figure out which pre-programed burn triggered the above conditional
    burn_index = (burns(:,1) <= t) & (burns(:,2) >= t);
    [e f] = max(int_plan(1,(int_plan(1,:)<=t)));
    curMass = int_plan(2,f) + (int_plan(2,f+1)-int_plan(2,f)).*(t-e)./(int_plan(1,f+1)-int_plan(1,f));
    %find thrust angle
    theta = burns(burn_index,6).*alt + burns(burn_index,7);
    if abs(theta) > abs(burns(burn_index,5))
        theta = burns(burn_index,5);
    end
    tm = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
    if burns(burn_index, 4) == 0
        uv = (scVelocity./((sum(scVelocity.^2)).^(1/2)))*tm;
    else
        uv = -(gAccel./((sum(gAccel.^2)).^(1/2)))*tm;
    end
    tAccel = burns(burn_index, 3).*uv./curMass;
end

g_alt = alt;
out = gAccel + dAccel + tAccel;
out(isnan(out)) = 0;
out(isinf(out)) = 0;
end

