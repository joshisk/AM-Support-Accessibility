
% function [POI_v1v2, POI_v2v3, POI_v3v1] = points_of_intersection(first_vertex, second_vertex, point_plane, normal_plane)
function [POI_v1v2, POI_v2v3, POI_v3v1] = points_of_intersection(vertex1, vertex2, vertex3, const)

  %unitnormal_plane = normal_plane/norm(normal_plane);
  %POI_v1v2 = vertex1 + (dot(unitnormal_plane, (point_plane - vertex1))*(vertex2 - vertex1))/dot(unitnormal_plane, (vertex2 - vertex1));
  POI_v1v2(1) = vertex1(1) + (const-vertex1(3))*(vertex2(1)-vertex1(1))/(vertex2(3)-vertex1(3)); 
  POI_v1v2(2) = vertex1(2) + (const-vertex1(3))*(vertex2(2)-vertex1(2))/(vertex2(3)-vertex1(3)); 
  POI_v1v2(3) = vertex1(3) + (const-vertex1(3))*(vertex2(3)-vertex1(3))/(vertex2(3)-vertex1(3)); 
  
  %POI_v2v3 = vertex2 + (dot(unitnormal_plane, (point_plane - vertex2))*(vertex3 - vertex2))/dot(unitnormal_plane, (vertex3 - vertex2));
  POI_v2v3(1) = vertex2(1) + (const-vertex2(3))*(vertex3(1)-vertex2(1))/(vertex3(3)-vertex2(3));
  POI_v2v3(2) = vertex2(2) + (const-vertex2(3))*(vertex3(2)-vertex2(2))/(vertex3(3)-vertex2(3));
  POI_v2v3(3) = vertex2(3) + (const-vertex2(3))*(vertex3(3)-vertex2(3))/(vertex3(3)-vertex2(3));
  
  %POI_v3v1 = vertex3 + (dot(unitnormal_plane, (point_plane - vertex3))*(vertex1 - vertex3))/dot(unitnormal_plane, (vertex1 - vertex3));
  POI_v3v1(1) = vertex3(1) + (const-vertex3(3))*(vertex1(1)-vertex3(1))/(vertex1(3)-vertex3(3));
  POI_v3v1(2) = vertex3(2) + (const-vertex3(3))*(vertex1(2)-vertex3(2))/(vertex1(3)-vertex3(3));
  POI_v3v1(3) = vertex3(3) + (const-vertex3(3))*(vertex1(3)-vertex3(3))/(vertex1(3)-vertex3(3));
  
end

