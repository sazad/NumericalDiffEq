% This M-File uses the Second Order Taylor Series Method
% to solve the initial value problem 
%     Y' = f(x,y)
%     Y(a) = alpha
% on the interval [a,b] with uniform mesh size h where
%   dfx is the x partial derivative of f, and
%   dfy is the y partial derivative of f
% and returns 
%   approx:(n+1) x 2 size matrix that contains x_i and the approximations y_i

function approx = taylortwo(f, alpha, a, b, h, dfx, dfy)
% f: function with two arguments, f is C^2
% alpha: Y(a) = alpha
% a: left endpoint = x0
% b: right endpoint = xn
% h: uniform mesh size
% dfx: func. with 2 args
% dfy: func. with 2 args

n = (b-a)/h;
approx = zeros(n+1,2);
approx(1,:) = [a alpha];

for i = 1:n
    xi = approx(i,1);
    yi = approx(i,2);
    ynext = yi + (f(xi,yi)*h) + (((h^2)/2)*(dfx(xi,yi) + (dfy(xi,yi)*f(xi,yi))));
    xnext = a + (i*h);
    approx(i+1,:) = [xnext ynext];
end
end