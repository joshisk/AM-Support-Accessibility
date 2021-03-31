clear all 
close all
clc
%% 

%[F,V,N]  = stlread('Z:\Shriyanka\Code\Matlab code\model2_bk.stl');
[F,V,N]  = stlread('Z:\Shriyanka\Modules\Accessibility\MATLAB code\Support_Generation.stl');

figure(1);
rotate3d on
hold on

%point_rayorigin = [120 -60 10];
point_rayorigin = [25 0 10];
%point_rayorigin = [110 -10 10];
%point_rayorigin = [40 10 15];
plot3(point_rayorigin(1),point_rayorigin(2),point_rayorigin(3),'*b')
dir_ray = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
pos_x = [1 0 0];
neg_x = [-1 0 0];
pos_y = [0 1 0];
neg_y = [0 -1 0];
pos_z = [0 0 1];
neg_z = [0 0 -1];
point_plane = [58.2957725524902,-80.1807022094727,20];
normal_plane = [-0.515038073062897,-0.857167303562164,0];

final_result = [];
accessibility_neg_x = [];
accessibility_pos_x = [];
accessibility_neg_y = [];
accessibility_pos_y = [];
accessibility_neg_z = [];
accessibility_pos_z = [];
for i=1:length(F)
 triangle = [V((i*3)-2:(i*3),:);V((i*3)-2,:)];
 plot3(triangle(:,1),triangle(:,2),triangle(:,3),'r')
 %face = [1 2 3];
 %patch('Faces',face,'Vertices',triangle,'FaceColor','c')

 result1 = poi_accessibility_new(point_rayorigin, pos_x, triangle);
 accessibility_neg_x = [accessibility_neg_x; result1];
 result2 = poi_accessibility_new(point_rayorigin, neg_x, triangle);
 accessibility_pos_x = [accessibility_pos_x; result2];
 result3 = poi_accessibility_new(point_rayorigin, pos_y, triangle);
 accessibility_neg_y = [accessibility_neg_y; result3];
 result4 = poi_accessibility_new(point_rayorigin, neg_y, triangle);
 accessibility_pos_y = [accessibility_pos_y; result4];
 result5 = poi_accessibility_new(point_rayorigin, pos_z, triangle);
 accessibility_neg_z = [accessibility_neg_z; result5]; 
 result6 = poi_accessibility_new(point_rayorigin, neg_z, triangle);
 accessibility_pos_z = [accessibility_pos_z; result6];

%  for j = 1:length(dir_ray)
% result = poi_accessibility_new(point_rayorigin, dir_ray(j,:), triangle);
% final_result = [final_result;result];
% end
end
% for j = 1:length(dir_ray)
% poi_accessibility_new(point_rayorigin, dir_ray(j,:), ver);
% end


if  accessibility_neg_x(:,1) == -1
acc_neg_x = 1;
Negative_X ='Point is Accessible'
else
acc_neg_x = 0;
Negative_X ='Point is Inaccessible'
end

if  accessibility_pos_x(:,1) == 1
acc_pos_x = 1;
Positive_X ='Point is Accessible'
else
acc_pos_x = 0;
Positive_X ='Point is Inaccessible'
end

if  accessibility_neg_y(:,2) == -1
acc_neg_y = 1;
Negative_Y ='Point is Accessible'
else
acc_neg_y = 0;
Negative_Y ='Point is Inaccessible'
end

if  accessibility_pos_y(:,2) == 1
acc_pos_y = 1;
Positive_Y ='Point is Accessible'
else
acc_pos_y = 0;
Positive_Y ='Point is Inaccessible'
end

if  accessibility_neg_z(:,3) == -1
acc_neg_z = 1;
Negative_Z ='Point is Accessible'
else
acc_neg_z = 0;
Negative_Z ='Point is Inaccessible'
end

if  accessibility_pos_z(:,3) == 1
acc_pos_z = 1;
Positive_Z ='Point is Accessible'
else
acc_pos_z = 0;
Positive_Z ='Point is Inaccessible'
end
















