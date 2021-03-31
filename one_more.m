clear all
close all
clc
syms t;
[F,V,N]  = stlread('Support_Generation.stl');
l=length(F);

alpha = 30;
beta = 30;

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


points=[]; % error lies here  ????
temp4=0;
for t=(linspace(0,1,5))
    if t~=1
vertex1(:,:) = [V_t((Sf(x)*3)-2,1)+  t*(V_t((Sf(x)*3)-2,1)-V_t((Sf(x)*3)-2,1)),V_t((Sf(x)*3)-2,2)+  t*(V_t((Sf(x)*3)-2,2)-V_t((Sf(x)*3)-2,2)),V_t((Sf(x)*3)-2,3)+  t*(Zv-V_t((Sf(x)*3)-2,3))];
vertex2(:,:) = [V_t((Sf(x)*3)-1,1)+  t*(V_t((Sf(x)*3)-1,1)-V_t((Sf(x)*3)-1,1)),V_t((Sf(x)*3)-1,2)+  t*(V_t((Sf(x)*3)-1,2)-V_t((Sf(x)*3)-1,2)),V_t((Sf(x)*3)-1,3)+  t*(Zv-V_t((Sf(x)*3)-1,3))];
vertex3(:,:) = [V_t((Sf(x)*3),1)+  t*(V_t((Sf(x)*3),1)-V_t((Sf(x)*3),1)),V_t((Sf(x)*3),2)+  t*(V_t((Sf(x)*3),2)-V_t((Sf(x)*3),2)),V_t((Sf(x)*3),3)+  t*(Zv-V_t((Sf(x)*3),3))];

points(temp4+1,:) = vertex1(1,1:3);
points(temp4+2,:) = vertex2(1,1:3);
points(temp4+3,:) = vertex3(1,1:3);

temp4=temp4+3;

    end
end
points_Mat = vertcat(points_Mat,points);

 end
end
hold on
final_points=unique(points_Mat,'rows');

%% Plotting supports
num_points=size(final_points);
for i=1:num_points(1)
scatter3([final_points(i,1),final_points(i,1)],[final_points(i,2),final_points(i,2)],[final_points(i,3),Zv(1)],'*b');
rotate3d on
hold on   
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')
end

%% Accessibility Analysis

%point_rayorigin = final_points;
direction = [1 2 3];
point_rayorigin = [50 48 97];
accessibility_direction = [];
Accessibility = [];

%for j = 1:length(point_rayorigin)
for i=1:length(F)
    
 triangle = [V_t((i*3)-2:(i*3),:);V_t((i*3)-2,:)];
 result = poi_accessibility_new(point_rayorigin, direction, triangle);
 accessibility_direction = [accessibility_direction; result];
 
end





























