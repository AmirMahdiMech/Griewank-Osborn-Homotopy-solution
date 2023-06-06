syms x y t
eps = 1;
a = 1 ; b= 0.1;
x_0 = [a;b];
H = [29*x^3/16-2*x*y-t*(29*a^3/16 - 2*a*b);y-x^2-t*(b-a^2)];
dHdxdy = jacobian(H,[x,y]);
dxdt = -dHdxdy\diff(H,t);
eps = dHdxdy\H;
dxdt = matlabFunction(dxdt,'Vars',{x y});
dHdxdy = matlabFunction(dHdxdy,'Vars',{x y});
H = matlabFunction(H,'Vars',{x y t});
eps = matlabFunction(eps,'Vars',{x y t});
It = [];
t_0 = 1;
delta = 1/2;counter1 = 0;counter2=0;
while norm(H(x_0(1),x_0(2),0)) > 1e-0 && counter1<5000
    x_tilda = x_0 - delta*dxdt(x_0(1),x_0(2));
    while norm(eps(x_tilda(1),x_tilda(2),t_0-delta*t_0)) > 1e-4 && counter2 < 20
        x_tilda = x_tilda - eps(x_tilda(1),x_tilda(2),t_0-delta);
        counter2 = counter2 +1;
        eps(x_tilda(1),x_tilda(2),t_0-delta*t_0);
    end
    if counter2 == 20 && norm(eps(x_tilda(1),x_tilda(2),t_0-delta)) > 1e-4 
        delta = delta/2;
        counter2=0;
    else
        t_0 = t_0 - delta;
        x_0 = x_tilda;
        It = [It [norm(x_0);delta]];
    end
    counter1 = counter1 +1;
end

    
 

    