function [total_B] = NV_cube(n, total_NV_centers)
% n is the number of configuratconfigurations
% total NV center is the number of NV centers (spin half electrons)
total_B = zeros(1,n);
for loop=1:n
    % simulating NV centers in a cube lattice x-y-z
    % grids = zeros(10,10,10);
    points = zeros(total_NV_centers,3);
    for i=1:total_NV_centers
        % 100x100 nm cube
        points(i,:) = randi([1 100],1,3);
    %     grids(point(1),point(2),point(3)) = 1;
    end
    
    % four tetrahedron direction
    dir = [1 1 1; 1 1 0; 1 0 1; 0 1 1];
    rand_dir = zeros(total_NV_centers,3);
    for i=1:total_NV_centers
        rand_dir(i,:) = dir(randi([1,4]),:);
    end
    % rand_dir = [dir(randi([1,4]),:); dir(randi([1,4]),:); dir(randi([1,4]),:); dir(randi([1,4]),:); dir(randi([1,4]),:); dir(randi([1,4]),:); dir(randi([1,4]),:); dir(randi([1,4]),:); dir(randi([1,4]),:); dir(randi([1,4]),:)];
    dis_mat = zeros(total_NV_centers,total_NV_centers);
    for i=1:total_NV_centers
        for j=i+1:total_NV_centers
            dis_mat(i,j) = norm(points(i,:) - points(j,:));
        end
    end
    
    total_omega = 0;
    for i=1:total_NV_centers
        for j=i+1:total_NV_centers
%             total_omega = total_omega + two_nv_Hamiltonian(dis_mat(i,j), rand_dir(i,:), rand_dir(j,:));
            total_omega = total_omega + Interaction(dis_mat(i,j), rand_dir(i,:), rand_dir(j,:), 0.5);
        end
    end
    total_B(loop) = total_omega;
end
return