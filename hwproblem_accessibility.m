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
%point_rayorigin = [-20 30 0];
%plot3(point_rayorigin(1),point_rayorigin(2),point_rayorigin(3),'*k')
%dir_ray = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];

point_rayorigin = final_points;

pos_x = [1 0 0];
neg_x = [-1 0 0];
pos_y = [0 1 0];
neg_y = [0 -1 0];
pos_z = [0 0 1];
neg_z = [0 0 -1];

final_result = [];
accessibility_neg_x = [];
accessibility_pos_x = [];
accessibility_neg_y = [];
accessibility_pos_y = [];
accessibility_neg_z = [];
accessibility_pos_z = [];
Accessibility = [];
for j = 1:length(point_rayorigin)
for i=1:length(F)
 triangle = [V_t((i*3)-2:(i*3),:);V_t((i*3)-2,:)];
 %plot3(triangle(:,1),triangle(:,2),triangle(:,3),'r')
 %face = [1 2 3];
 %patch('Faces',face,'Vertices',triangle,'FaceColor','c')

 result1 = poi_accessibility_new(point_rayorigin(j,:), pos_x, triangle);
 accessibility_neg_x = [accessibility_neg_x; result1];
 result2 = poi_accessibility_new(point_rayorigin(j,:), neg_x, triangle);
 accessibility_pos_x = [accessibility_pos_x; result2];
 result3 = poi_accessibility_new(point_rayorigin(j,:), pos_y, triangle);
 accessibility_neg_y = [accessibility_neg_y; result3];
 result4 = poi_accessibility_new(point_rayorigin(j,:), neg_y, triangle);
 accessibility_pos_y = [accessibility_pos_y; result4];
 result5 = poi_accessibility_new(point_rayorigin(j,:), pos_z, triangle);
 accessibility_neg_z = [accessibility_neg_z; result5]; 
 result6 = poi_accessibility_new(point_rayorigin(j,:), neg_z, triangle);
 accessibility_pos_z = [accessibility_pos_z; result6];

end


if  accessibility_neg_x(:,1) == -1
acc_neg_x = 1;
else
acc_neg_x = 0;
end

if  accessibility_pos_x(:,1) == 1
acc_pos_x = 1;
else
acc_pos_x = 0;
end

if  accessibility_neg_y(:,2) == -1
acc_neg_y = 1;
else
acc_neg_y = 0;
end

if  accessibility_pos_y(:,2) == 1
acc_pos_y = 1;
else
acc_pos_y = 0;
end

if  accessibility_neg_z(:,3) == -1
acc_neg_z = 1;
else
acc_neg_z = 0;
end

if  accessibility_pos_z(:,3) == 1
acc_pos_z = 1;
else
acc_pos_z = 0;
end
Acc_eachpoint = [acc_neg_x acc_pos_x acc_neg_y acc_pos_y acc_neg_z acc_pos_z];
Accessibility = [Accessibility; Acc_eachpoint];
end
negx_percent = (numel(find(Accessibility(:,1)==1))/numel(point_rayorigin(:,1)))*100
posx_percent = (numel(find(Accessibility(:,2)==1))/numel(point_rayorigin(:,1)))*100
negy_percent = (numel(find(Accessibility(:,3)==1))/numel(point_rayorigin(:,1)))*100
posy_percent = (numel(find(Accessibility(:,4)==1))/numel(point_rayorigin(:,1)))*100
negz_percent = (numel(find(Accessibility(:,5)==1))/numel(point_rayorigin(:,1)))*100
posz_percent = (numel(find(Accessibility(:,6)==1))/numel(point_rayorigin(:,1)))*100
all = [negx_percent; posx_percent; negy_percent; posy_percent; negz_percent; posz_percent];
sequence = sort(all,'descend')

























