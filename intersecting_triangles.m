
function check = intersecting_triangles(vertex1, vertex2, vertex3, point_plane, normal_plane)

distance1 = distance_point_to_plane(vertex1, point_plane, normal_plane);
distance2 = distance_point_to_plane(vertex2, point_plane, normal_plane);
distance3 = distance_point_to_plane(vertex3, point_plane, normal_plane);

if(distance1 > 0 && distance2 > 0 && distance3 > 0)||(distance1 < 0 && distance2 < 0 && distance3 < 0)
    check = 0;
else
    check = 1;
end