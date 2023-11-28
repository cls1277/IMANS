function [GD]=convergence(obj,ref_point, fmin, fmax)%传入已经归一化的目标值和参考点
    [N,~]  = size(obj);
%     obj = (obj-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
    obj_size = N;
    [ref_size,~]=size(ref_point);
    distance=zeros(obj_size,ref_size);
    GD=0;
    for i=1:obj_size
        for j=1:ref_size
            distance(i,j)=(obj(i,1)-ref_point(j,1))^2+(obj(i,2)-ref_point(j,2))^2;
            distance(i,j)=sqrt(distance(i,j));
        end
        GD=GD+min(distance(i,:));
    end
    GD=GD/obj_size;
end