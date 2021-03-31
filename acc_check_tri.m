clear all 
close all
clc
%% 
support_point = [120 -50 10];
ray_dir = [1 0 0];

figure(1);
rotate3d on
hold on
axis equal
xlabel('X-Axis')
ylabel('Y-Axis')
zlabel('Z-Axis')
plot3(support_point(1),support_point(2),support_point(3),'*k')
hold on
quiver3(support_point(:,1),support_point(:,2),support_point(:,3),ray_dir(:,1),ray_dir(:,2),ray_dir(:,3))
hold on
%%
[F,V,N]  = stlread('Z:\Shriyanka\Code\Matlab code\model2_bk.stl');

triangle = [];
for i=1:length(F)
 
 temp_triangle = [V((i*3)-2:(i*3),:);V((i*3)-2,:)];
 plot3(temp_triangle(:,1),temp_triangle(:,2),temp_triangle(:,3),'r')
 %face = [1 2 3];
 %patch('Faces',face,'Vertices',triangle,'FaceColor','c')
 triangle = [triangle; temp_triangle];
 
end

length_stl_faces = length(F);
color = accessibility_check_triangle(support_point, ray_dir, triangle, length_stl_faces)