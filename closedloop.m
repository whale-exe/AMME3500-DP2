function dxdt = closedloop(t,xt,A,B,C,Kgain,kr,Lgain,r)

    x = xt(1:3);
    xh = xt(4:6);
    u = -Kgain*xh + kr*r;
    y = C*x;
    dxdt = [A*x + B*u; A*xh + B*u + Lgain*(y-C*xh)];

end
