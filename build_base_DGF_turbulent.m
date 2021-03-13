function [G_star] = build_base_DGF_turbulent(N,m)

% compute base DGF in dimensionless form
b1 = 1/9;
b2 = 9/10;
G_star = zeros(N);
for i = 1:N
    for j = 1:i
        if i == j
            G_star(i,j) = 1/((1-b1)*b2)*(j^b2-(j-1)^b2)^(1-b1);
        else
            G_star(i,j) = 1/((1-b1)*b2)*((i^b2-(j-1)^b2)^(1-b1) - ...
                                        ((i-1)^b2-(j-1)^b2)^(1-b1) - ...
                                        (i^b2-j^b2)^(1-b1) + ...
                                        ((i-1)^b2-j^b2)^(1-b1));
        end
    end
end

end