function plot_stl(V)

triangle = [];
for i=1:length(V)/3
    
temp_triangle = [V((i*3)-2:(i*3),:);V((i*3)-2,:)];
 %plot3(temp_triangle(:,1),temp_triangle(:,2),temp_triangle(:,3),'r');
 %face = [1 2 3];
 %patch('Faces',face,'Vertices',triangle,'FaceColor','c')
triangle = [triangle; temp_triangle];
plot3(triangle(((4*i-3):(4*i)),1),triangle(((4*i-3):(4*i)),2),triangle(((4*i-3):(4*i)),3),'r');

end
end