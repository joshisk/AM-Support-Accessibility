
function result = poi_accessibility(point_rayorigin, dir_ray, point_plane, normal_plane)

 unitnormal_plane = normal_plane/norm(normal_plane);
 lambda = dot(unitnormal_plane, (point_plane - point_rayorigin))/dot(unitnormal_plane, dir_ray)
 if (lambda > 0)
 POI = point_rayorigin + lambda*dir_ray
 else
 POI = [nan nan nan]
 end
 
 if isnan(POI(1))||isnan(POI(2))||isnan(POI(3))||isinf(POI(1))||isinf(POI(2))||isinf(POI(3))
 result = -dir_ray
 else  
 result = [nan nan nan; -dir_ray]
 end

end