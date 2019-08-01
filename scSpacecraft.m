function [out] = scSpacecraft(scVi, scPi)
global cf;
global A;
global store;
global int_plan;
int_plan = [];
run(cf);
out = {};

for i = length(stages(:,1)):-1:1
    if i == length(stages(:,1))
        t = 0;
        isStaging = 0;
    else
        t = stages(i+1,4);
        isStaging = 1;
    end

    scMass = sum(sum(stages(1:i,[1,2])));
    isFuel = +(stages(i,2))>0;
    
    int_plan = [int_plan, [t;scMass;isFuel;isStaging]];
end
int_plan = sortrows(int_plan',1)';
    
for i = 1:length(burns(:,1))
    startBurnT = (burns(i,1));
    endBurnT = (burns(i,2));
    startMass = int_plan(2,(int_plan(1,:)<=startBurnT));
    startMass = startMass(end);
    %find number of prior staging events
    stagings = sum(int_plan(4,(int_plan(1,:)<=startBurnT)));
    %find the mass after the stage has run out of fuel
    newMass = sum(sum(stages((1:end-stagings),1:2)))-stages((end-stagings),2);
    dM = startMass - newMass;
    Ve = stages((1:end-stagings),3);
    outFuelTime = startBurnT+(Ve.*dM./abs(burns(i,3)));
    if outFuelTime < endBurnT
        endTime = outFuelTime;
        isFuel = 0;
    else
        endTime = endBurnT;
        newMass = startMass-(endBurnT - startBurnT).*abs(burns(i,3))./Ve(end);
        isFuel = 1;
    end
    isStaging = 0;
    int_plan = [int_plan, [startBurnT;startMass;isFuel;isStaging], [endTime;newMass;isFuel;isStaging]];
    int_plan = sortrows(int_plan',1)';
end
initialVPA = [scVi, scPi, 1];
int_plan = [int_plan, [runTime;int_plan(2,end);0;0]];

%optionset = odeset('RelTol', 1e-5, 'AbsTol', 1e-4);

for w = 1:(length(int_plan(1,:))-1)
    %find a time well inside the iteration
    %aveTime = (deltaT(w)+deltaT(w+1))./2;
    
    
    
    %(any((burns(:,1) <= (t + h)) & (burns(:,2) >= t))) && int_plan(3,w) ~= 1
    
    %scMass = [scMass [deltaT(w+1); afterIterMass]];
    
    optionset = odeset('Events',@events,'RelTol', 1e-8, 'AbsTol', 1e-7);
    tic
    [tArray,VelPosArray]=ode113(@(t,y) scODE113(t,y,int_plan(3,w)), int_plan(1,(w):(w+1)), initialVPA, optionset);
    toc
    sStore = cell(7,length(tArray));
    for i=1:length(tArray)
        sStore{1,i} = tArray(i);
        sStore{4,i} = VelPosArray(i,1:3);
        sStore{5,i} = VelPosArray(i,4:6);
        
        %use interpolating polynomials
        t = tArray(i);
        [a b] = max(A(A<=t));
        if (b+2 < length(A))
            L0 = (t-A(b));
            L1 = (t-A(b+1));
            L2 = (t-A(b+2));
            L01 = A(b)-A(b+1);
            L02 = A(b)-A(b+2);
            L12 = A(b+1)-A(b+2);
            W1 = L1.*L2./(L01.*L02);
            W2 = -L0.*L2./(L01.*L12);
            W3 = L0.*L1./(L02.*L12);
            sStore{2,i} = W1.*(store{2, b}) + W2.*(store{2, (b+1)}) + W3.*(store{2, (b+2)});
            sStore{3,i} = W1.*(store{3, b}) + W2.*(store{3, (b+1)}) + W3.*(store{3, (b+2)});
        else
            L0 = (t-A(b-2));
            L1 = (t-A(b-1));
            L2 = (t-A(b));
            L01 = A(b-2)-A(b-1);
            L02 = A(b-2)-A(b);
            L12 = A(b-1)-A(b);
            W1 = L1.*L2./(L01.*L02);
            W2 = -L0.*L2./(L01.*L12);
            W3 = L0.*L1./(L02.*L12);
            sStore{2,i} = W1.*(store{2, b-2}) + W2.*(store{2, (b-1)}) + W3.*(store{2, (b)});
            sStore{3,i} = W1.*(store{3, b-2}) + W2.*(store{3, (b-1)}) + W3.*(store{3, (b)});
        end
        dist = [0,0];
        for j=1:length(masses)
            dist(j) = (sum((sStore{5,i}-sStore{2,i}(j,:)).^2)).^(1/2);
        end
        [c d] = min(dist);
        alt = (dist(d)-radii(d))./1000;
        sStore{6, i} = dist;
        sStore{7, i} = alt;
        
        %linear interpolation to find mass. Mass change is linear so this is exact.
        [e f] = max(int_plan(1,(int_plan(1,:)<=t)));
        if t < max(int_plan(1,:))&&(int_plan(4,f+1)==0)
            sStore{8,i} = int_plan(2,f) + (int_plan(2,f+1)-int_plan(2,f)).*(t-e)./(int_plan(1,f+1)-int_plan(1,f));       
        else
            sStore{8,i} = int_plan(2,f);
        end
    end
    
    if sStore{1,end}<int_plan(1,w+1)
        sStore = [sStore onGround(sStore(:,end), int_plan(1,w+1))];
    end
    initialVPA = [sStore{4,end}, sStore{5,end}, sStore{7,end}];
    out = [out sStore(:,2:end)];
end

end

function [out alt] = scODE113(t, in, outFuelTime)
scVi = in(1:3);
scPi = in(4:6);

dpdt = scVi';
[dvdt alt] = scAccel(t, scVi', scPi',outFuelTime);
dvdt = dvdt';
dpdt = dpdt';
out = [dvdt; dpdt; alt];

end

function [value,isterminal,direction] = events(t, in)
global g_alt;
value = g_alt;
isterminal = 1;
direction = -1;
end

function [out] = onGround(initialData, tf)
global store;
global A;
global cf;
run(cf); 
out = cell(8,1000);
out(:,1) = initialData(:);
ti = initialData{1};
%since we will be using non-linear interpolation for the graph, we need to
%make the point of impact as stiff as possible
newT = linspace(ti, tf, 602);
newT = [linspace(newT(1), newT(1)+0.000001, 400) newT(3:end)];


[x body_index] = min(initialData{6});
initial_uv = (initialData{5}-initialData{2}(body_index,:))./((sum((initialData{5}-initialData{2}(body_index,:)).^2)).^(1/2));
for i = 1:1000
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
    out{2,i} = pPositions;
    out{3,i} = pVelocities;
    out{7,i} = 0;
    out{8,i} = initialData{8};
end
for i = 1:1000
    r_shift = 2.*pi.*(newT(i) - initialData{1})./rotations(body_index);
    tm = [cos(-r_shift), -sin(-r_shift), 0; sin(-r_shift), cos(-r_shift), 0; 0, 0, 1];
    newUV = initial_uv*tm;
    magRotVel = (2.*pi.*radii(body_index))./rotations(body_index);
    newUVr = initial_uv*[cos(-r_shift-pi./2), -sin(-r_shift-pi./2), 0; sin(-r_shift-pi./2), cos(-r_shift-pi./2), 0; 0, 0, 1];
    
    out{4,i} = magRotVel.*newUVr;
    out{5,i} = out{2,i}(body_index,:) + newUV.*radii(body_index);
    
    for j = 1:length(masses)
        dist(j) = (sum((out{5,i}-out{2,i}(j,:)).^2)).^(1/2);
    end
    
    out{6,i} = dist;
end

out = out(:,2:end);
end

% function [] = apAuto()
% 
% if (((cur_dist < old_dist&&(apAutoBurns(ap_burn_index,1)==1))||(cur_dist > old_dist&&(apAutoBurns(ap_burn_index,1)==0)))||(autoBurn(1) > 1))
%             uvSC = ((scVelocity-store{3, i}(body_ID,:))./((sum((scVelocity-store{3, i}(body_ID,:)).^2)).^(1/2)));     
%             autoBurn(1) = autoBurn(1) + 1;
%             if autoBurn(1)==2
%                 need_DV_vec = uvSC.*(6.67300.*10.^(-11).*masses(apAutoBurns(ap_burn_index,2))./cur_dist).^(1/2) - scVelocity + store{3, i}(body_ID,:);
%                 need_DV_mag = sum((need_DV_vec.^2)).^(1./2);
%             end
%             engineDVperIt = stages(end, 3).*(log(sum(sum(stages(:,[1,2])))) - log(sum(sum(stages(:,[1,2]))) - abs(apAutoBurns(ap_burn_index,5).*h./stages(end,3))));
%             totDVapplyed = ((sum((autoBurn([2,3,4])).^2)).^(1/2));
%             whole_iters = (need_DV_mag-totDVapplyed)./engineDVperIt;
%             if whole_iters >= 1
%                 bTime = h;
%                 deltaV = engineDVperIt + deltaV;
%                 oldmass = sum(sum(stages(:,[1,2])));
%                 stages(end,2) = stages(end,2) - abs(burns(burn_index,3).*bTime./stages(end,3));
%                 dDV = uvSC.*stages(end, 3).*log(oldmass./sum(sum(stages(:,[1,2]))));
%                 tDV = tDV + dDV;
%                 autoBurn([2,3,4]) = autoBurn([2,3,4]) + dDV;
%             else
%                 bTime = (stages(end, 3).*sum(sum(stages(:,[1,2]))./apAutoBurns(ap_burn_index,5)).*(1-1./(exp(engineDVperIt.*whole_iters./stages(end, 3)))));
%                 oldmass = sum(sum(stages(:,[1,2])));
%                 stages(end,2) = stages(end,2) - abs(burns(burn_index,3).*bTime./stages(end,3));
%                 dDV = uvSC.*stages(end, 3).*log(oldmass./sum(sum(stages(:,[1,2]))));
%                 tDV = tDV + dDV;
%                 autoBurn(:) = 0;
%             end
% 
% end
% 
%     %compute thrust from auto burns
%     if (any((apAutoBurns(:,3) <= (t)) & (apAutoBurns(:,4) >= t))) && Out_of_Fuel ~= 1 && (i>1)&&(autoBurn(1)~=0)
%         ap_burn_index = (apAutoBurns(:,3) <= (t) & (apAutoBurns(:,4) >= t));
%         body_ID = apAutoBurns(ap_burn_index,2);
%         cur_dist = sStore{3,i}(body_ID);
%         old_dist = sStore{3,(i-1)}(body_ID);
%         
%         
%         
%             %this resets the autoBurn conditional so that it can be
%             %triggered by later auto burns
%             if i == apAutoBurns(ap_burn_index,4)
%                 autoBurn(1) = 1;
%             end
%             
%         end
%     end

