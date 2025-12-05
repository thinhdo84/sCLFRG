function xd = vtol(x,u,epsilon)

xd = [x(2);
    - u(1)*sin(x(5)) + epsilon*u(2)*cos(x(5));
    x(4);
    u(1)*cos(x(5)) + epsilon*u(2)*sin(x(5))-1;
    x(6);
    u(2)    ];
end