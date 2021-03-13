function G = configure_DGF(G_lam,G_turb,ind_orig,ind_trans)

N = size(G_lam,1);

% stitch together laminar and turbulent DGFs according to the transition index
G_hybrid = G_lam;
if ind_trans<N
    G_hybrid(ind_trans+1:end,:) = G_turb(ind_trans+1:end,:);
end

% reindex base DGF according to specified origin
G = zeros(N);
G(ind_orig+1:end,ind_orig+1:end) = G_hybrid(1:N-ind_orig,1:N-ind_orig);
G(ind_orig:-1:1,ind_orig:-1:1) = G_hybrid(1:ind_orig,1:ind_orig);

end