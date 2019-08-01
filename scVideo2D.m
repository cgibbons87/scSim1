function [] = scVideo(sStore,stStore)
global cf;
global int_plan;
run(cf);

B = cell2mat(sStore(1,:));
C = [];
for i = 1:length(stStore) 
    C = [C(:) {cell2mat(stStore{i}(1,:))}];
end

t = 0;
vStore = cell(7,10000);

  for i=1:10000
    h = B(end)./10000;
    t = t + h;
    [a b] = max(B(B<=t));
    vStore{1, i} = t;
    if (b+2 < length(B))
        L0 = (t-B(b));
        L1 = (t-B(b+1));
        L2 = (t-B(b+2));
        L01 = B(b)-B(b+1);
        L02 = B(b)-B(b+2);
        L12 = B(b+1)-B(b+2);
        W1 = L1.*L2./(L01.*L02);
        W2 = -L0.*L2./(L01.*L12);
        W3 = L0.*L1./(L02.*L12);
        for j=2:length(sStore(:, 1))
            vStore{j, i} = W1.*(sStore{j, b}) + W2.*(sStore{j, (b+1)}) + W3.*(sStore{j, (b+2)});
        end
    else
        L0 = (t-B(b-2));
        L1 = (t-B(b-1));
        L2 = (t-B(b));
        L01 = B(b-2)-B(b-1);
        L02 = B(b-2)-B(b);
        L12 = B(b-1)-B(b);
        W1 = L1.*L2./(L01.*L02);
        W2 = -L0.*L2./(L01.*L12);
        W3 = L0.*L1./(L02.*L12);
        for j=2:length(sStore(:, 1))
            vStore{j, i} = W1.*(sStore{j, b-2}) + W2.*(sStore{j, (b-1)}) + W3.*(sStore{j, (b)});
        end
    end
    if (vStore{7, i} < -0.01)&&(b+1 < length(B))
        for j=2:length(sStore(:, 1))-1
            vStore{j, i} = sStore{j, b} + (sStore{j,b+1}-sStore{j,b}).*(t-sStore{1,b})./(sStore{1,b+1}-sStore{1,b});
        end
    elseif (vStore{7, i} < -0.01)&&(b+1 > length(B))
        for j=2:length(sStore(:, 1))-1
            vStore{j, i} = sStore{j, b} + (sStore{j,b}-sStore{j,b}).*(t-sStore{1,b})./(sStore{1,b}-sStore{1,b});
        end
    end
    
    [e f] = max(int_plan(1,(int_plan(1,:)<=t)));
    if t < max(int_plan(1,:))&&(int_plan(4,f+1)==0)
        vStore{8,i} = int_plan(2,f) + (int_plan(2,f+1)-int_plan(2,f)).*(t-e)./(int_plan(1,f+1)-int_plan(1,f));
    else
        vStore{8,i} = int_plan(2,f);
    end
    
    %spent stages.
    for j=1:length(stStore)
        [aa bb] = max(C{j}(C{j}<=t));
        dist = [];
        if t >= C{j}(1)&&bb+2<=length(C{j})
            L0 = (t-C{j}(bb));
            L1 = (t-C{j}(bb+1));
            L2 = (t-C{j}(bb+2));
            L01 = C{j}(bb)-C{j}(bb+1);
            L02 = C{j}(bb)-C{j}(bb+2);
            L12 = C{j}(bb+1)-C{j}(bb+2);
            W1 = L1.*L2./(L01.*L02);
            W2 = -L0.*L2./(L01.*L12);
            W3 = L0.*L1./(L02.*L12);
            vStore{8+j,i} = W1.*stStore{j}{3,bb} + W2.*stStore{j}{3,bb+1} + W3.*stStore{j}{3,bb+2};
            %yet another distance calc.
            for k = 1:length(masses)
                dist(k) = (sum((vStore{8+j,i}-vStore{2,i}(k,:)).^2)).^(1/2);
            end
        elseif t >= C{j}(1)
            vStore{8+j,i} = stStore{j}{3,bb};
            for k = 1:length(masses)
                dist(k) = (sum((vStore{8+j,i}-vStore{2,i}(k,:)).^2)).^(1/2);
            end
        end
        
        %use linear interp if altitude becomes negative
        [m mm] = min(dist);
        if (t >= C{j}(1))&&(m<radii(mm))&&(bb+1<=length(C{j}))
            vStore{8+j,i} = stStore{j}{3,bb} + (stStore{j}{3,bb+1}-stStore{j}{3,bb}).*(t-stStore{j}{1,bb})./(stStore{j}{1,bb+1}-stStore{j}{1,bb});
        elseif t >= C{j}(1)&&m<radii(mm)
            vStore{8+j,i} = stStore{j}{3,bb};
        end 
        
    end
  end


for i = 1:length(vStore)
    h = 10;
    t = vStore{1,i};
    positions = vStore{2, i};
    scVelocity = vStore{4, i};
    scPosition = vStore{5, i};
    dist = vStore{6, i};
    dV = 0;
    alt = vStore{7, i};
    [a b] = min(dist(:));
    scale = abs(1 + 5.*radii(b).*abs(alt)./dist(b) + abs(alt));
    
    axis equal
    axis([scPosition(1)./1000-scale,scPosition(1)./1000+scale,scPosition(2)./1000-scale,scPosition(2)./1000+scale]);
    
    hold on
    
    %for circle generation
    a = linspace(0,2.*pi,10000);
    cosA = cos(a);
    sinA = sin(a);
    
    %plot circles one at a time
    for j = 1:length(masses)
        if alt < 1000 && (j==b)
            x = radii(j).*cosA./1000 + positions(j,1)./1000;
            y = radii(j).*sinA./1000 + positions(j,2)./1000;
            plot(x,y, 'k')
        else
            rectangle('Position',[(positions(j,1)-radii(j))./1000,(positions(j,2)-radii(j))./1000,2.*radii(j)./1000,2.*radii(j)./1000],'Curvature',[1,1]);
        end
    end
    
    %add spacecraft and spacecraft data to graph
    plot(scPosition(1)./1000,scPosition(2)./1000, 'x')
    text(scPosition(1)./1000-scale+0.05.*scale,scPosition(2)./1000+scale-0.08.*scale, ['t(s)  = ', num2str(t)])
    text(scPosition(1)./1000-scale+0.05.*scale,scPosition(2)./1000+scale-0.16.*scale, ['dx = ', num2str(round(scVelocity(1)./10)./100)])
    text(scPosition(1)./1000-scale+0.05.*scale,scPosition(2)./1000+scale-0.24.*scale, ['dy = ', num2str(round(scVelocity(2)./10)./100)])
    text(scPosition(1)./1000-scale+0.05.*scale,scPosition(2)./1000+scale-0.32.*scale, ['dz = ', num2str(round(scVelocity(3)./10)./100)])
    text(scPosition(1)./1000-scale+0.05.*scale,scPosition(2)./1000+scale-0.40.*scale, ['magV = ', num2str(round((((scVelocity(1)^2)+(scVelocity(2).^2)+(scVelocity(3)^2)).^(1/2))./10)./100)])
    text(scPosition(1)./1000-scale+0.05.*scale,scPosition(2)./1000+scale-0.48.*scale, ['Mass(tons) = ', num2str((vStore{8,i}./100)./10)])
    text(scPosition(1)./1000-scale+0.05.*scale,scPosition(2)./1000+scale-0.56.*scale, ['dV  = ', num2str(round(dV./10)./100)])
    text(scPosition(1)./1000-scale+0.05.*scale,scPosition(2)./1000+scale-0.64.*scale, ['Alt = ', num2str(alt)])
    
    
    %add spent stages, if any
    if length(vStore(:,1))>8
        for j=9:length(vStore(:,1))
            if ~isempty(vStore{j,i})
                plot(vStore{j,i}(1)./1000,vStore{j,i}(2)./1000, 'o')
            end
        end
    end
    
    hold off
    %Grab a frame
    M(1) = getframe;
    
    %clear the graph
    clf
end
end
