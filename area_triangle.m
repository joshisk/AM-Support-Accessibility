function area = area_triangle(triangle)

area = 1/2*norm(cross(triangle(1,:)-triangle(3,:), triangle(2,:)-triangle(3,:)));

end