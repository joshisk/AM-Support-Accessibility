function color = accessibility_check_triangle(support_point, ray_dir, triangle, length_stl_faces)
%%% , pos_x_index, pos_y_index, pos_z_index, neg_x_index, neg_y_index, neg_z_index

pos_x_indices = [];
neg_x_indices = [];
pos_y_indices = [];
neg_y_indices = [];
pos_z_indices = [];
neg_z_indices = [];
%temp_posx_index = (triangle(4*i,1)>= point(1))
for i = 1:length_stl_faces
%%
if ray_dir(1) >= 0
   pos_x_index = (triangle(4*i,1) >= support_point(:,1));
else
   neg_x_index = (triangle(4*i,1) < support_point(:,1));
end
%%
if ray_dir(2) >= 0
   pos_y_index = (triangle(4*i,2) >= support_point(:,2));
else
   neg_y_index = (triangle(4*i,2) < support_point(:,2));
end
%%
if ray_dir(3) >= 0
   pos_z_index = (triangle(4*i,3) >= support_point(:,3));
else
   neg_z_index = (triangle(4*i,3) < support_point(:,3));
end
%%
pos_x_indices = [pos_x_indices; pos_x_index];
pos_y_indices = [pos_y_indices; pos_y_index];
pos_z_indices = [pos_z_indices; pos_z_index];
neg_x_indices = [neg_x_indices; neg_x_index];
neg_y_indices = [neg_y_indices; neg_y_index];
neg_z_indices = [neg_z_indices; neg_z_index];
%%
color = [pos_x_indices pos_y_indices pos_z_indices neg_x_indices neg_y_indices neg_z_indices];
figure(1);
rotate3d on
hold on
axis equal
xlabel('X-Axis')
ylabel('Y-Axis')
zlabel('Z-Axis')
if color(i,1)==1
plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'g') 
end
if color(i,2)==1
plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'b') 
end
if color(i,3)==1
plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'c') 
end
if color(i,4)==1
plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'y') 
end
if color(i,5)==1
plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'m') 
end
if color(i,6)==1
plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'k') 
end
end

end