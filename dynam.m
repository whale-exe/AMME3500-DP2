function dxdt = dynam(t,x,A,B,Kgain,kr,r)

    u = -Kgain*x + kr*r;
    dxdt = A*x + B*u;
end