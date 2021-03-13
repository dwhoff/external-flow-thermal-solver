function [rho,nu,cp,k] = get_fluid_props(fl_obj,species,P,T)

switch species
    case 'air'
        set(fl_obj,'T',T,'P',P);
        rho = density(fl_obj);
        mu = viscosity(fl_obj);
        nu = mu/rho;
        cp = cp_mass(fl_obj);
        k = thermalConductivity(fl_obj);
    case 'water'
        set(fl_obj,'T',T,'P',P);
        rho = density(fl_obj);
        cp = cp_mass(fl_obj);
        % NASA polynomials for transport properties
        % https://ntrs.nasa.gov/citations/19940013151
        % 300<T<1000
        mu = 1e-5*exp(0.7839*log(T)-0.3826e3/T+0.4904e5/T^2+0.8522);
        k = 1e-4*exp(0.1554e1*log(T)+0.6611e2/T+0.5597e4/T^2-0.3926); % uW/cm.K
        nu = mu/rho;
end

end