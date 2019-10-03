clc
clear
close all
%import kml file
blub = kml2struct('george.kml');

%select first year vector to do skeletonzation for base centerline
lon = blub(1).Lon(1:end-1); lat = blub(1).Lat(1:end-1);
%convert to x and y (utm coordinates, wgs84), units are meters
[x,y] = ll2utm(lat,lon);
%make a skeleton
[v, e] = voronoiSkel([x/1e6 y/1e6]','boundary','fast',3); v = v.*1e6;
%plot land vector
figure(1);
plot(x,y), hold on
patch('faces',e,'vertices',v)
%find endnodes, by first finding the distance between the end-points of the
%graph.
ep = find(histc(e(:),unique(e(:)))==1);
[EP1,EP2] = meshgrid(v(ep,1),v(ep,2));
d = sqrt((EP1-EP1').^2+(EP2-EP2').^2);
%endnodes of the island are where end-points are furthest away from
%eachother
[~,idx,idy] = max2d(d); ep = [ep(idx) ep(idy)];
%find shortest path between end-points
g = (graph(e(:,1),e(:,2)));
path = shortestpath(g,ep(1),ep(2));
xp = v(path,1);
yp = v(path,2);
%smooth centerline, makes transects look better
[xp,yp] = smooth_contours(xp,yp,10);
plot(xp,yp), hold on

%resample along constant distance for transects
t = [0;cumsum(sqrt(diff(xp).^2+diff(yp).^2))];
ti = 0:100:t(end);
xl = interp1(t,xp,ti);
yl = interp1(t,yp,ti);
plot(xl,yl,'o')
%orthogonal lines along transect, increase length of lines by changing the
%15000 in the function for v.
v = [xl(2:end); yl(2:end)] - [xl(1:end-1); yl(1:end-1)];
v = [v, v(:,end)]; v = 15000*v / norm(v);
tr_x = [xl+v(2,:); xl-v(2,:)];
tr_y = [yl-v(1,:); yl+v(1,:)];
plot(tr_x,tr_y,'k'), axis equal
%
%get intersects of transects w/ barrier island to extract width
xint = zeros(size(tr_x)); yint = zeros(size(tr_y));
for ii=1:length(tr_x),
    [xint,yint] = polyxpoly(x,y,tr_x(:,ii),tr_y(:,ii));
    
    %get width by finding the largest distance between intersections.
    w(ii) = max2d(sqrt(bsxfun(@minus,xint,xint').^2+bsxfun(@minus,yint,yint').^2));
  

end
%plot for all years
figure(2);
yr = [1984:2018];
m=[1:35]
for i=m
    LON = blub(i).Lon(1:end-1); LAT = blub(i).Lat(1:end-1);
    [X,Y] = ll2utm(LAT,LON);
figure(2);
plot(X,Y,'DisplayName',num2str(yr(i))), hold on
    
    
    XINT = zeros(size(tr_x)); YINT = zeros(size(tr_y));
    for ii=1:length(tr_x),
    [XINT,YINT] = polyxpoly(X,Y,tr_x(:,ii),tr_y(:,ii));
    [Y_beach,beach_index]=min(YINT);
    X_beach=XINT(beach_index);
    [Y_back,back_index]=max(YINT);
    X_back=XINT(back_index);
    %get width by finding the largest distance between intersections.
    W(ii) = max2d(sqrt(bsxfun(@minus,XINT,XINT').^2+bsxfun(@minus,YINT,YINT').^2));
    
    %find width from centerline to the back and to the front of island
    W_beach(ii)=sqrt((X_beach-xl(ii))^2 + (Y_beach-yl(ii))^2);
    W_back(ii)=sqrt((X_back-xl(ii))^2 + (Y_back-yl(ii))^2);
    end
    
    CW(:,i)=[W i];
    CW_beach(:,i)=[W_beach i];
    CW_back(:,i)=[W_back i];
end
axis equal
% get rid of last row
CW(end,:) = [];
CW_beach(end,:) = [];
CW_back(end,:) = [];

%use cumulative maximum function to disregard noise,vegetation
CW = CW(:,m);
CW_beach = CW_beach(:,m);
CW_back = CW_back(:,m);
CW_back_cummax = cummax(CW_back,2);

%plot width
figure(3);
plot(CW);
%figure(4);
%plot(CW_beach);
%plot back width
figure(5);
plot(CW_back,'DisplayName',num2str(yr));
%calculate overwash and plot it
figure(6);
CW_back_change = diff(CW_back_cummax,1,2);
CW_beach_change = diff(CW_beach_cummin,1,2);
CW_change = diff(CW,1,2);
plot(CW_back_change);


CW_first26years = CW(:,1:(end-1));

%interpolate overwash and width 
figure(7)

for i=m
txt = ['Year ',num2str(1984+i)];
scatter(CW_first26years(:,i),CW_back_change(:,i),'filled','black','DisplayName',txt)

hold on

end


