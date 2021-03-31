function result = poi_accessibility_new(point_rayorigin, dir_ray, triangle)

 point_plane = triangle(1,:);
 normal_plane = cross((triangle(2,:)-triangle(1,:)),(triangle(3,:)-triangle(1,:)));
 unitnormal_plane = normal_plane/norm(normal_plane);
 lambda = dot(unitnormal_plane, (point_plane - point_rayorigin))/dot(unitnormal_plane, dir_ray);
 
 if (lambda > 0)
 POI = point_rayorigin + lambda*dir_ray;
 else
 POI = [nan nan nan];
 end
 decision = insideoutsidetest(triangle, POI);
 %[in, on] = inpolygon(POI(:,1),POI(:,2),triangle(:,1),triangle(:,2));
 %if isnan(POI(1))||isnan(POI(2))||isnan(POI(3))||isinf(POI(1))||isinf(POI(2))||isinf(POI(3))||~in||~on
 if isnan(POI(1))||isnan(POI(2))||isnan(POI(3))||isinf(POI(1))||isinf(POI(2))||isinf(POI(3))||decision == 0
 result = -dir_ray;
 %fprintf('\nPOI is accessible in -dir_ray\n');
 else  
 result = [nan nan nan];
 %fprintf('\nPOI is inaccessible in -dir_ray\n');
 end

end