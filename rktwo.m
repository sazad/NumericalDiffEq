% This M-File uses the second order Runge-Kutta method to
% solve the initial value problem 
%     Y' = f(x,y)
%     Y(a) = alpha
% on interval [a,b] with uniform mesh size h and returns
%     approx: (n+1) x 2 size matrix of points with the values of
%             x_i and approximations of y_i
%             
function approx = rktwo(f, alpha, a, b, h) 
% f: function that takes two arguments
% alpha: Y(a) = alpha = y0
% a: left endpoint = x0
% b: right endpoint = xn
% h: length of uniform mesh

n = (b-a)/h;
approx = zeros(n+1,2);
approx(1, :) = [a alpha]; 

for i = 1:n
    xi = approx(i,1);
    yi = approx(i,2);
    k1 = f(xi,yi);
    k2 = f(xi+h, yi+h*k1);
    ynext = yi + h*((k1+k2)/2);
    xnext = a + (i*h);
    approx(i+1,:) = [xnext ynext];
end
end
