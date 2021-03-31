clear all 
close all
clc
%% 
%point = [110 -10 10];
point = [25 0 10];

figure(1);
rotate3d on
hold on
axis equal
xlabel('X-Axis')
ylabel('Y-Axis')
zlabel('Z-Axis')
plot3(point(1),point(2),point(3),'*k')
hold on
%%
%[F,V,N]  = stlread('Z:\Shriyanka\Code\Matlab code\model2_bk.stl');
[F,V,N]  = stlread('Z:\Shriyanka\Modules\Accessibility\MATLAB code\Support_Generation.stl');
triangle = [];
pos_x_indices = [];
pos_y_indices = [];
pos_z_indices = [];
neg_x_indices = [];
neg_y_indices = [];
neg_z_indices = [];

for i=1:length(F)
 %%
 temp_triangle = [V((i*3)-2:(i*3),:);V((i*3)-2,:)];
 plot3(temp_triangle(:,1),temp_triangle(:,2),temp_triangle(:,3),'r')
 %face = [1 2 3];
 %patch('Faces',face,'Vertices',triangle,'FaceColor','c')
 triangle = [triangle; temp_triangle];
%  %%
%  temp_posx_index = (triangle(4*i,1)>= point(1));
%  temp_posy_index = (triangle(4*i,2)>= point(2));
%  temp_posz_index = (triangle(4*i,3)>= point(3));
%  
%  temp_negx_index = (triangle(4*i,1)< point(1));
%  temp_negy_index = (triangle(4*i,1)< point(2));
%  temp_negz_index = (triangle(4*i,1)< point(3));
%  
%  pos_x_indices = [pos_x_indices; temp_posx_index];
%  pos_y_indices = [pos_y_indices; temp_posy_index];
%  pos_z_indices = [pos_z_indices; temp_posz_index];
%  
%  neg_x_indices = [neg_x_indices; temp_negx_index];
%  neg_y_indices = [neg_y_indices; temp_negy_index];
%  neg_z_indices = [neg_z_indices; temp_negz_index];

%%
%  if pos_x_indices(i)==1
%  plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'g') 
%  end
 
% %%
%  if pos_y_indices(i)==1
%   plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'c')     
%  end
% 
% %%
%  if pos_z_indices(i)==1
%   plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'b')     
%  end
% 
% %%
% if neg_x_indices(i)==1
%  plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'g') 
% end
% 
% %%
% if neg_y_indices(i)==1
%  plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'c') 
% end
% 
% %%
% if neg_z_indices(i)==1
%  plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'b') 
% end


end

%%
figure(2);
rotate3d on
hold on
axis equal
xlabel('X-Axis')
ylabel('Y-Axis')
zlabel('Z-Axis')
plot3(point(1),point(2),point(3),'*k')
for i=1:length(F)
 %%
 temp_triangle = [V((i*3)-2:(i*3),:);V((i*3)-2,:)];
 plot3(temp_triangle(:,1),temp_triangle(:,2),temp_triangle(:,3),'r')
 %face = [1 2 3];
 %patch('Faces',face,'Vertices',triangle,'FaceColor','c')
 triangle = [triangle; temp_triangle];
 %%
 temp_posx_index = (triangle(4*i,1)>= point(1));
 pos_x_indices = [pos_x_indices; temp_posx_index];
 
%%
 if pos_x_indices(i)==1
 plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'g') 
 end
end

figure(3);
rotate3d on
hold on
axis equal
xlabel('X-Axis')
ylabel('Y-Axis')
zlabel('Z-Axis')
plot3(point(1),point(2),point(3),'*k')
for i=1:length(F)
 %%
 temp_triangle = [V((i*3)-2:(i*3),:);V((i*3)-2,:)];
 plot3(temp_triangle(:,1),temp_triangle(:,2),temp_triangle(:,3),'r')
 %face = [1 2 3];
 %patch('Faces',face,'Vertices',triangle,'FaceColor','c')
 triangle = [triangle; temp_triangle];
 %%
 temp_posy_index = (triangle(4*i,2)>= point(2));
 pos_y_indices = [pos_y_indices; temp_posy_index];
 
%%
 if pos_y_indices(i)==1
  plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'c')     
 end
end

%%
figure(4);
rotate3d on
hold on
axis equal
xlabel('X-Axis')
ylabel('Y-Axis')
zlabel('Z-Axis')
plot3(point(1),point(2),point(3),'*k')
for i=1:length(F)
 %%
 temp_triangle = [V((i*3)-2:(i*3),:);V((i*3)-2,:)];
 plot3(temp_triangle(:,1),temp_triangle(:,2),temp_triangle(:,3),'r')
 %face = [1 2 3];
 %patch('Faces',face,'Vertices',triangle,'FaceColor','c')
 triangle = [triangle; temp_triangle];
 %%
 temp_posz_index = (triangle(4*i,3)>= point(3));
 pos_z_indices = [pos_z_indices; temp_posz_index];
%%
 if pos_z_indices(i)==1
  plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'b')     
 end
end

%%
figure(5);
rotate3d on
hold on
axis equal
xlabel('X-Axis')
ylabel('Y-Axis')
zlabel('Z-Axis')
plot3(point(1),point(2),point(3),'*k')
for i=1:length(F)
 %%
 temp_triangle = [V((i*3)-2:(i*3),:);V((i*3)-2,:)];
 plot3(temp_triangle(:,1),temp_triangle(:,2),temp_triangle(:,3),'r')
 %face = [1 2 3];
 %patch('Faces',face,'Vertices',triangle,'FaceColor','c')
 triangle = [triangle; temp_triangle];
 %%
 temp_negx_index = (triangle(4*i,1)< point(1));
 neg_x_indices = [neg_x_indices; temp_negx_index];
%%
if neg_x_indices(i)==1
 plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'g') 
end
end

%%
figure(6);
rotate3d on
hold on
axis equal
xlabel('X-Axis')
ylabel('Y-Axis')
zlabel('Z-Axis')
plot3(point(1),point(2),point(3),'*k')
for i=1:length(F)
 %%
 temp_triangle = [V((i*3)-2:(i*3),:);V((i*3)-2,:)];
 plot3(temp_triangle(:,1),temp_triangle(:,2),temp_triangle(:,3),'r')
 %face = [1 2 3];
 %patch('Faces',face,'Vertices',triangle,'FaceColor','c')
 triangle = [triangle; temp_triangle];
 %%
 temp_negy_index = (triangle(4*i,2)< point(2));
 neg_y_indices = [neg_y_indices; temp_negy_index];
 
%%
if neg_y_indices(i)==1
 plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'c') 
end
end

%%
figure(7);
rotate3d on
hold on
axis equal
xlabel('X-Axis')
ylabel('Y-Axis')
zlabel('Z-Axis')
plot3(point(1),point(2),point(3),'*k')
for i=1:length(F)
 %%
 temp_triangle = [V((i*3)-2:(i*3),:);V((i*3)-2,:)];
 plot3(temp_triangle(:,1),temp_triangle(:,2),temp_triangle(:,3),'r')
 %face = [1 2 3];
 %patch('Faces',face,'Vertices',triangle,'FaceColor','c')
 triangle = [triangle; temp_triangle];
 %%
 temp_negz_index = (triangle(4*i,3)< point(3));
 neg_z_indices = [neg_z_indices; temp_negz_index];

%%
if neg_z_indices(i)==1
 plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'b') 
end
end,

