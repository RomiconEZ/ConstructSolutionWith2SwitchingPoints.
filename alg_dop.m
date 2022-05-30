function [Dad] = ad(i,h,D)
D(i,:) = [];
D(:,h) = [];
Dad = det(D);
if mod((i+h),2) ~= 0
    Dad = Dad * (-1);
end
