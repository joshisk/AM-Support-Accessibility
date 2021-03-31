clc
clear all
clear vars
close all
rotate3d on
syms t;
[F,V,N]  = stlread('Support_Generation.stl');
l=length(F);

alpha = 30;
beta = 45;

Rx=[1     0           0        
    0 cosd(alpha)  sind(alpha) 
    0 -sind(alpha) cosd(alpha) ];

Ry=[cosd(beta) 0 -sind(beta)  
    0          1      0      
    sind(beta) 0 cosd(beta)];

V_t= V*Rx*Ry;
N_new = N*Rx*Ry;
angles_with_normals =acosd(N_new(:,3));
%% Plot part STL
figure(1);
 axis equal        
for i = 1:l
    f = [1 2 3];
    v = [V_t(F(i,1),1:3); V_t(F(i,2),1:3); V_t(F(i,3),1:3)];
    patch('Faces',f,'Vertices',v,'FaceColor','g')  

end

% for x=1:1:l
%  plot3(V_t((x*3)-2:(x*3),1), V_t((x*3)-2:(x*3),2), V_t((x*3)-2:(x*3),3),'g'); 
%  hold on 
%  plot3(V_t([(x*3)-2,(x*3)],1), V_t([(x*3)-2,(x*3)],2), V_t([(x*3)-2,(x*3)],3),'g'); 
%  hold on 
% end
hold on

%% plotting substrate 
Dx= (max(V_t(:,1))-min(V_t(:,1)))/20;
Dy= (max(V_t(:,2))-min(V_t(:,2)))/20;
Dz=(max(V_t(:,3))-min(V_t(:,3)));

Xv=[max(V_t(:,1))+Dx,min(V_t(:,1))-Dx];
Yv=[max(V_t(:,2))+Dy,min(V_t(:,2))-Dy];

% changed to vary the distance btwn substrate and part
%Zv=[min(V_t(:,3))-Dz/5,min(V_t(:,3))-Dz/2.5];
Zv=[min(V_t(:,3))-3,0];

Xb=[Xv(1),Xv(1),Xv(2),Xv(2)];
Yb=[Yv(1),Yv(2),Yv(1),Yv(2)];
Zb=[Zv(1),Zv(1),Zv(1),Zv(1)];

V_F= [Xv(2),Yv(2),Zv(1);Xv(2),Yv(1),Zv(1);Xv(1),Yv(1),Zv(1);Xv(1),Yv(2),Zv(1)];
Face_num=[1 2 3 4];
   patch('Faces',Face_num,'Vertices',V_F,'FaceColor','y')  
hold on


%% finding out the stl triangles requiring supports
s_num=0;d_num=0;

 for i=1:1:l    
    if angles_with_normals(i)>= 135       
       s_num=s_num+1;     
    end
    
    if  angles_with_normals(i)> 90 &&  angles_with_normals(i)< 135   
        d_num=d_num+1;
    end         
 end
 
support_requiring_faces=zeros(1,s_num);
 down_facing_surfaces = zeros(1,d_num);
 temp1=0;temp2=0;
 for i=1:1:l    
    if angles_with_normals(i)>= 135   
        temp1=temp1+1;
        support_requiring_faces(temp1) = i;      
    end
    
    if  (angles_with_normals(i)> 90 &&  angles_with_normals(i)< 135 )
        temp2=temp2+1;
        down_facing_surfaces(temp2) = i;
    end         
 end

Sf=support_requiring_faces;
Df= down_facing_surfaces;

 s=size(Sf);
d=size(Df);
areas=zeros(s(2),1);
edge1=zeros(s(2),3);  edge2=zeros(s(2),3);  edge3=zeros(s(2),3);

temp4=0;
temp3=zeros(1,s(2));
%points_Mat= zeros(s(2),1)
points_Mat=[];
if s(2)>0    
 for x=1:s(2) 
  plot3(V_t((Sf(x)*3)-2:(Sf(x)*3),1), V_t((Sf(x)*3)-2:(Sf(x)*3),2), V_t((Sf(x)*3)-2:(Sf(x)*3),3),'r'); 
 hold on 
 plot3(V_t([(Sf(x)*3)-2,(Sf(x)*3)],1), V_t([(Sf(x)*3)-2,(Sf(x)*3)],2), V_t([(Sf(x)*3)-2,(Sf(x)*3)],3),'r'); 
 hold on 
 
 a=sqrt((V_t((Sf(x)*3)-2,1)-V_t((Sf(x)*3)-1,1))^2  + (V_t((Sf(x)*3)-2,2)-V_t((Sf(x)*3)-1,2))^2  +  (V_t((Sf(x)*3)-2,3)-V_t((Sf(x)*3)-1,3))^2);
b=sqrt((V_t((Sf(x)*3)-1,1)-V_t((Sf(x)*3),1))^2  + (V_t((Sf(x)*3)-1,2)-V_t((Sf(x)*3),2))^2  +  (V_t((Sf(x)*3)-1,3)-V_t((Sf(x)*3),3))^2 );
c=sqrt((V_t((Sf(x)*3),1)-V_t((Sf(x)*3)-2,1))^2  + (V_t((Sf(x)*3),2)-V_t((Sf(x)*3)-2,2))^2  +  (V_t((Sf(x)*3),3)-V_t((Sf(x)*3)-2,3))^2);
abc=[a,b,c];
n= round(max(abc)/5);
half_S= 0.5*(a+b+c);
areas(x,1)= sqrt(half_S*(half_S-a)*(half_S-b)*(half_S-c));


for t=(linspace(0,1,n))
    if t~=1
temp3(x)=temp3(x)+3;
    end
end

points=zeros(temp3(x),3); % error lies here  ????
temp4=0;
for t=(linspace(0,1,n))
    if t~=1
edge1(i,:) = [V_t((Sf(x)*3)-2,1)+  t*(V_t((Sf(x)*3)-1,1)-V_t((Sf(x)*3)-2,1)),V_t((Sf(x)*3)-2,2)+  t*(V_t((Sf(x)*3)-1,2)-V_t((Sf(x)*3)-2,2)),V_t((Sf(x)*3)-2,3)+  t*(V_t((Sf(x)*3)-1,3)-V_t((Sf(x)*3)-2,3))];
edge2(i,:) = [V_t((Sf(x)*3)-1,1)+  t*(V_t((Sf(x)*3),1)-V_t((Sf(x)*3)-1,1)),V_t((Sf(x)*3)-1,2)+  t*(V_t((Sf(x)*3),2)-V_t((Sf(x)*3)-1,2)),V_t((Sf(x)*3)-1,3)+  t*(V_t((Sf(x)*3),3)-V_t((Sf(x)*3)-1,3))];
edge3(i,:) = [V_t((Sf(x)*3),1)+  t*(V_t((Sf(x)*3)-2,1)-V_t((Sf(x)*3),1)),V_t((Sf(x)*3),2)+  t*(V_t((Sf(x)*3)-2,2)-V_t((Sf(x)*3),2)),V_t((Sf(x)*3),3)+  t*(V_t((Sf(x)*3)-2,3)-V_t((Sf(x)*3),3))];

points(temp4+1,:) = edge1(i,:);
points(temp4+2,:) = edge2(i,:);
points(temp4+3,:) = edge3(i,:);

temp4=temp4+3;

    end
end
points_Mat = vertcat(points_Mat,points);

 end
end
hold on
 final_points=unique(points_Mat,'rows');
%final_points=points_Mat;

% if s(2)>0    
%  for x=1:s(2) 
%   plot3(V_t((Sf(x)*3)-2:(Sf(x)*3),1), V_t((Sf(x)*3)-2:(Sf(x)*3),2), V_t((Sf(x)*3)-2:(Sf(x)*3),3),'r'); 
%  hold on 
%  plot3(V_t([(Sf(x)*3)-2,(Sf(x)*3)],1), V_t([(Sf(x)*3)-2,(Sf(x)*3)],2), V_t([(Sf(x)*3)-2,(Sf(x)*3)],3),'r'); 
%  hold on 
%  
%  a=sqrt((V_t((Sf(x)*3)-2,1)-V_t((Sf(x)*3)-1,1))^2  + (V_t((Sf(x)*3)-2,2)-V_t((Sf(x)*3)-1,2))^2  +  (V_t((Sf(x)*3)-2,3)-V_t((Sf(x)*3)-1,3))^2);
% b=sqrt((V_t((Sf(x)*3)-1,1)-V_t((Sf(x)*3),1))^2  + (V_t((Sf(x)*3)-1,2)-V_t((Sf(x)*3),2))^2  +  (V_t((Sf(x)*3)-1,3)-V_t((Sf(x)*3),3))^2 );
% c=sqrt((V_t((Sf(x)*3),1)-V_t((Sf(x)*3)-2,1))^2  + (V_t((Sf(x)*3),2)-V_t((Sf(x)*3)-2,2))^2  +  (V_t((Sf(x)*3),3)-V_t((Sf(x)*3)-2,3))^2);
% abc=[a,b,c];
% n= round(max(abc)/5);
% half_S= 0.5*(a+b+c);
% areas(x,1)= sqrt(half_S*(half_S-a)*(half_S-b)*(half_S-c));
% 
% temp3=0;
% 
% for t=(linspace(0,1,n))
%     if t~=1
% temp3=temp3+3;
%     end
% end
% 
% points=zeros(temp3,3); % error lies here  ????
% 
% for t=(linspace(0,1,n))
%     if t~=1
%  edge1(i,:) = [V_t((Sf(x)*3)-2,1)+  t*(V_t((Sf(x)*3)-1,1)-V_t((Sf(x)*3)-2,1)),V_t((Sf(x)*3)-2,2)+  t*(V_t((Sf(x)*3)-1,2)-V_t((Sf(x)*3)-2,2)),V_t((Sf(x)*3)-2,3)+  t*(V_t((Sf(x)*3)-1,3)-V_t((Sf(x)*3)-2,3))];
% edge2(i,:) = [V_t((Sf(x)*3)-1,1)+  t*(V_t((Sf(x)*3),1)-V_t((Sf(x)*3)-1,1)),V_t((Sf(x)*3)-1,2)+  t*(V_t((Sf(x)*3),2)-V_t((Sf(x)*3)-1,2)),V_t((Sf(x)*3)-1,3)+  t*(V_t((Sf(x)*3),3)-V_t((Sf(x)*3)-1,3))];
% edge3(i,:) = [V_t((Sf(x)*3),1)+  t*(V_t((Sf(x)*3)-2,1)-V_t((Sf(x)*3),1)),V_t((Sf(x)*3),2)+  t*(V_t((Sf(x)*3)-2,2)-V_t((Sf(x)*3),2)),V_t((Sf(x)*3),3)+  t*(V_t((Sf(x)*3)-2,3)-V_t((Sf(x)*3),3))];
% 
% points(temp4+1,:) = edge1(i,:);
% points(temp4+2,:) = edge2(i,:);
% points(temp4+3,:) = edge3(i,:);
% 
% temp4=temp4+3;
% 
%     end
% end
% 
% 
%  end
% end

%%  Finding out the downfacing surfaces
% if d(2)>0    
%  for y=1:d(2)
%   plot3(V_t((Df(y)*3)-2:(Df(y)*3),1), V_t((Df(y)*3)-2:(Df(y)*3),2), V_t((Df(y)*3)-2:(Df(y)*3),3),'y'); 
%  hold on 
%  plot3(V_t([(Df(y)*3)-2,(Df(y)*3)],1), V_t([(Df(y)*3)-2,(Df(y)*3)],2), V_t([(Df(y)*3)-2,(Df(y)*3)],3),'y'); 
%  hold on
%  end
% end
% hold on

%% Plotting supports

 num_points=size(final_points);

% scatter3(final_points(:,1),final_points(:,2),final_points(:,3),'filled');
% hold on
for i=1:num_points(1)
   
     plot3([final_points(i,1),final_points(i,1)],[final_points(i,2),final_points(i,2)],[final_points(i,3),Zv(1)],'-.b');
    hold on
    
    
end
%% volume calculation 
sum_vol=0;
sum_area=0;

for i=1:s(2)
 
 v1 =  [V_t(Sf(i)*3-2,1),V_t(Sf(i)*3-2,2),V_t(Sf(i)*3-2,3)];
 v2 =  [V_t(Sf(i)*3-1,1),V_t(Sf(i)*3-1,2),V_t(Sf(i)*3-1,3)];
 v3 =  [V_t(Sf(i)*3,1),V_t(Sf(i)*3,2),V_t(Sf(i)*3,3)];
 
 v_1 = [v1(1),v1(2),Zv(1)];
 v_2 = [v2(1),v2(2),Zv(1)];
 v_3 = [v3(1),v3(2),Zv(1)];
 
 
area_v1v2v3 = 1/2*cross((v2-v1),(v3-v1));
vector1 = (v1-v1);
normal1 = cross((v3-v1),(v2-v1))  ;  
unit_normal1=normal1/norm(normal1);
% if isnan(unit_normal1);
%     unit_normal1=[0 0 0]
% else
%     unit_normal1=normal1/norm(normal1);
% end
vol_v1v1v2v3 = 1/3*dot(cross(area_v1v2v3, vector1),unit_normal1);

area_v1v2v_2v_1 = cross((v2-v1),(v_1-v1));
length2 = (v1-v1);
normal2 = cross((v2-v1),(v_1-v1));
unit_normal2=normal2/norm(normal2);
% if isnan(unit_normal2)
%     unit_normal2=[0 0 0]
% else
%     unit_normal2=normal2/norm(normal2);
% end
vol_v1v1v2v_2v_1= 1/3*dot(cross(area_v1v2v_2v_1, length2), unit_normal2);

 
 
area_v1v2v3 = norm(1/2*cross((v2-v1),(v3-v1)));
length1 = (v1-v1);
normal1 = (cross((v3-v1),(v2-v1)));
unit_normal1=normal1/norm(normal1);
vol_v1v1v2v3 = 1/3*(area_v1v2v3*dot(length1, unit_normal1));


area_v1v2v_2v_1 = norm(cross((v2-v1),(v_1-v1)));
length2 = (v1-v1);
normal2 = (cross((v2-v1),(v_1-v1)));
unit_normal2=normal2/norm(normal2);
vol_v1v1v2v_2v_1= 1/3*(area_v1v2v_2v_1*dot(length2, unit_normal2));

area_v1v_1v_3v3 = norm(cross((v_1-v1),(v3-v1)));
length3 = (v1-v1);
normal3 = cross((v_1-v1),(v3-v1));
unit_normal3=normal3/norm(normal3);
vol_v1v1v_1v_3v3 = 1/3*(area_v1v_1v_3v3*dot(length3, unit_normal3));

area_v_1v_2v_3 = norm(1/2*cross((v_2-v_1),(v_3-v_1)));
length4 = (v_1-v1);
normal4 = cross((v_2-v_1),(v_3-v_1));
unit_normal4=normal4/norm(normal4);
vol_v1v_1v_2v_3 = 1/3*(area_v_1v_2v_3*dot(length4, unit_normal4));

area_v2v3v_3v_2 = norm(cross((v3-v2),(v_2-v2)));
length5 = (v2-v1);
normal5 = cross((v3-v2),(v_2-v2));
unit_normal5=normal5/norm(normal5);
vol_v1v2v3v_3v_2 = 1/3*(area_v2v3v_3v_2*dot(length5, unit_normal5));

vol = abs(vol_v1v1v2v3+vol_v1v1v2v_2v_1+vol_v1v1v_1v_3v3+vol_v1v_1v_2v_3+vol_v1v2v3v_3v_2);

sum_vol=sum_vol +vol;
sum_area=sum_area + area_v1v2v3;
end

fprintf('Support contact area is %.2f: \n', sum_area);
fprintf('Support volume is %.2f: \n', sum_vol);

