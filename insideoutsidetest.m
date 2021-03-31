function decision = insideoutsidetest(triangle, point)

small_triangle1 = [triangle(1:2,:);point]; 
small_triangle2 = [triangle(2:3,:);point]; 
small_triangle3 = [triangle(1,:);triangle(3,:);point]; 

A1 = area_triangle(small_triangle1);
A2 = area_triangle(small_triangle2);
A3 = area_triangle(small_triangle3);

A = area_triangle(triangle);

if  A == A1+A2+A3
   decision = 1;
   'Point lies inside the triangle';
else
   decision = 0;
   'Point lies outside the triangle';
end

end