
function distance = distance_point_to_plane(point_any, point_plane, normal_plane)

% normal_plane = cross(point2_plane-point1_plane, point3_plane-point1_plane)
% unit_normal = normal_plane/norm(normal_plane)
% distance = (dot(unit_normal, point1_plane) - dot(unit_normal, point_any))
distance = dot(normal_plane, (point_plane - point_any))/norm(normal_plane);

end