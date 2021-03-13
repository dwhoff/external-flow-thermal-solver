function Tnew = advance_RK4(f,T,dt)
% explicit RK4 time-stepping scheme for ODEs of the form --> dT/dt = f(T) 

k1 = dt*f(T);
k2 = dt*f(T+0.5*k1);
k3 = dt*f(T+0.5*k2);
k4 = dt*f(T+k3);
Tnew = T+(k1+2*k2+2*k3+k4)/6;

end