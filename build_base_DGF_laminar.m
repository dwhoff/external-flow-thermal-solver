function [G_star,b1,b2] = build_base_DGF_laminar(N,m)

% compute base DGF in dimensionless form
b1 = 0.042*m^0.26+1/3;
b2 = 3/4*(m+1)+(m+1)*(1-3*b1)/4/(b1-1);
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