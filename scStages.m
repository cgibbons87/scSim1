function [out] = scStages(sStore)
global cf;
global int_plan;
global g_alt;
g_alt = [];
run(cf);
out = {};
B = cell2mat(sStore(1,:));

stagingTimes = int_plan(1,logical(int_plan(4,:)));

for w = 1:length(stagingTimes)
    initialData = sStore(1:5,(stagingTimes(w)==B));
    
    optionset = odeset('Events',@events,'RelTol', 1e-8, 'AbsTol', 1e-7);
    tic
    [tArray,VelPosArray]=ode113(@(t,y) scODEstage(t,y), [stagingTimes(w),runTime], [initialData{4}, initialData{5}], optionset);
    toc
    stageStore = cell(3,length(tArray));
    
    for i=1:length(tArray)
        stageStore{1,i} = tArray(i);
        stageStore{2,i} = VelPosArray(i,1:3);
        stageStore{3,i} = VelPosArray(i,4:6);
    end
    
    if stageStore{1,end}<runTime
        stageStore = [stageStore onGround(stageStore(:,end), runTime)];
    end
    stageStore = {stageStore};
    out = [out stageStore];
end


end

function out = scODEstage(t, in)
scVi = in(1:3);
scPi = in(4:6);

dpdt = scVi';
dvdt = stageAccel(t, scVi', scPi');
dvdt = dvdt';
dpdt = dpdt';
out = [dvdt; dpdt];
end

function [value,isterminal,direction] = events(t, in)
global g_alt;
value = g_alt;
isterminal = 1;
direction = -1;
end

function out = stageAccel(t, stageVelocity, stagePosition)
global cf;
global A;
global store;
global g_alt;

run(cf);
%use interpolation to find the postion and velocity for each planet
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
    dist(i) = (sum((stagePosition-pPositions(i,:)).^2)).^(1/2);
    if dist(i) ~= 0
        %find the unit vector in the direction of the force
        uv = (stagePosition-pPositions(i,:))./(dist(i));
        %compute net acceleration of body n
        gAccel = gAccel - (6.67300e-11).*masses(i).*uv./(dist(i).^2);
        %gAccel(isinf(gAccel)) = 0;
    else
    end
end

%acceleration due to drag
[a b] = min(dist);
alt = (dist(b)-radii(b))./1000;
% dAccel = 0;
% if (atmP(b) ~= 0)&&(alt<100)
%     %calculate the geopotential altitude
%     gpAlt = radii(b) - (radii(b)^2)./dist(b);
%     H = 6700; %atmopheric scale height
%     rhoO = 1.752;
%     rho = rhoO.*exp(-gpAlt./H);
%     atm_uv = (pPositions(b,:)-scPosition)./((sum(((pPositions(b,:)-scPosition).^2)).^(1/2)));
%     atm_uv_r = atm_uv*[cos(pi./2), -sin(pi./2), 0; sin(pi./2), cos(pi./2), 0; 0, 0, 1];
%     atm_vel = (2.*pi.*radii(b))./rotations(b).*atm_uv_r;
%     relVel = scVelocity - pVelocities(b,:) - atm_vel;
%     uv = (relVel./((sum(relVel.^2)).^(1/2)));
%     dAccel = -uv.*0.5.*rho.*scDD(1).*scDD(2).*(sum(relVel.^2))./sum(sum(stages(:,[1,2])));
%     dAccel(isnan(dAccel)) = 0;
%     dAccel(isinf(dAccel)) = 0;
% end

g_alt = alt;
out = gAccel;% + dAccel;
out(isnan(out)) = 0;
out(isinf(out)) = 0;
end

function [out] = onGround(initialData, tf)
global store;
global A;
global cf;
run(cf);
out = cell(length(initialData),1000);

ti = initialData{1};
%since we will be using non-linear interpolation for the graph, we need to
%make the point of impact as stiff as possible
newT = linspace(ti, tf, 602);
newT = [linspace(newT(1), newT(1)+0.000001, 400) newT(3:end)];

for i = 1:1000
    out(1:3,i) = initialData(:);
    t = newT(i);
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
        pVelocities = W1.*(store{3, d}) + W2.*(store{3, (d+1)}) + W3.*(store{2, (d+2)});
    else
        pPositions = store{2, d};
        pVelocities = store{3, d};
    end
    
    out{1,i} = t;
    out{4,i} = pPositions;
    out{5,i} = pVelocities;
end
%find closest planet
for i = 1:length(masses)
    dist(i) = (sum((out{3,1}-out{4,1}(i,:)).^2)).^(1/2);
end
[x body_index] = min(dist);

initial_uv = (initialData{3}-out{4,1}(body_index,:))./((sum((initialData{3}-out{4,1}(body_index,:)).^2)).^(1/2));

for i = 1:1000
    r_shift = 2.*pi.*(newT(i) - initialData{1})./rotations(body_index);
    tm = [cos(-r_shift), -sin(-r_shift), 0; sin(-r_shift), cos(-r_shift), 0; 0, 0, 1];
    newUV = initial_uv*tm;
    magRotVel = (2.*pi.*radii(body_index))./rotations(body_index);
    newUVr = initial_uv*[cos(-r_shift-pi./2), -sin(-r_shift-pi./2), 0; sin(-r_shift-pi./2), cos(-r_shift-pi./2), 0; 0, 0, 1];
    
    out{2,i} = magRotVel.*newUVr;
    out{3,i} = out{4,i}(body_index,:) + newUV.*radii(body_index);
end
out = out(1:3,2:end);
end
