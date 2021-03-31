clc
clear all
close all

[F,V,N]  = stlread('Z:\Srikanth\line plot\6.STL file read\model2.stl');


figure(1);
axis equal
for x=1:length(F)
 plot3(V((x*3)-2:(x*3),1), V((x*3)-2:(x*3),2), V((x*3)-2:(x*3),3),'r'); 
 hold on 
 plot3(V([(x*3)-2,(x*3)],1), V([(x*3)-2,(x*3)],2), V([(x*3)-2,(x*3)],3),'r'); 
 hold on 
end
