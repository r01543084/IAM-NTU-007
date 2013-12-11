%THOMAS 2D ALGORITHM
function [u] = thomas2d(aj,bj,cj,dj,u)

    for i = 2:length(u(1,:))
        cj(i,:) = cj(i,:)-bj(i-1)/dj(i-1)*cj(i-1,:);
    end
    
    u(end,:) = cj(end,:)/dj(end);%back substitution according
    
    for i = length(u(1,:))-1:-1:1
        u(i,:) = (cj(i,:)-aj(i)*u(i+1,:))/dj(i);
    end