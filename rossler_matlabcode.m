% Solution of an autonomous ODE y' = f(y) (Rössler attractor)
% using RK4 and RK5 methods where y(y1, y2, ..., yn)
% prints the solution to matrix y.

tic; clear;
f = inline ( '[-y(2) - y(3); y(1) + 0.1 * y(2); 0.1 + y(3) * (y(1) - 14)]', 'y' ); % equation
order = 5; % The order of the Runge-Kutta method must be 4 or 5
y0 = [1 1 0]'; % This is the initial point
k = 5000; % k+1 is the number of columns of the solution matrix y
dt = 0.05; % The step size of t for the numerical method
error = 10^(-4); % This is the maximum error of y
plotrange = [-20 15 -150 20 0 150]; % The region to be plotted
timelimit = 3000; % The maximum time allowed for the program to run.
dim = size(y0, 1); % The dimension of the system
m = size(y0,2); % The number of initial points
tmax = timelimit/4;

fprintf ('This code uses Runge-Kutta method of order %-2.0f\n', order)
fprintf ('starting with the initial point %2.2f %2.2f %2.2f\n', y0')
fprintf ('from time t=0 to t = %-3.0f', k * dt)
fprintf ('The stepsize is ? t = %2.12g\n', dt)
fprintf ('It will continue until the error falls below %2.12g', error)
fprintf ('or time limit %3.0f', timelimit)
fprintf ('seconds is exceeded\n\n\n')
for q = 1 : m
y = zeros (dim, k + 1);
y(:, 1) = y0 (:, q);
yprevious = repmat (inf, dim, k + 1);
ttt = 0;
e = inf;
p = 0;
while ttt<tmax & e>error
n = 2^p;
h = dt/n;
ta = 0;
ya = y0 (:, q);
t0 = clock;
for j = 1: n * k
if order = 4 % This is the 4th order RK method
k1 = h * f (ya);
k2 = h * f (ya + k1/2);
k3 = h * f (ya + k2/2);
k4 = h * f (ya + k3);
ya = ya + (k1 + 2 * k2 + 2 * k3 + k4)/6;
ta = j * h;
else if order = 5 % This is the 5th order RK method
k1 = h * f (ya);
k2 = h * f (ya + k1/2);
k3 = h * f (ya + (3 * k1 + k2)/16);
k4 = h * f (ya + k3/2);
k5 = h * f (ya+ (-3 * k2 + 6 * k3 + 9 * k4)/16);
k6 = h * f (ya + (k1 + 4 * k2 + 6 * k3-12 * k4 + 8 * k5)/7);
ya = ya + (7 * k1 + 32 * k3 + 12 * k4 + 32 * k5 + 7 * k6)/90;
ta = j * h;
else fprintf('The order of the method must be 4 or 5\n')
break
end
if mod (j, n)=0; 
i = j/n;
y (:, i+1) = ya;
end
end
e = max (max (abs (y-yprevious)));
yprevious = y;
fprintf ('The step size is dt/%-2.0f', 2^p)
fprintf (' = %2.12g\n', h)
fprintf ('The estimated error < %2.12g\n', e)
ttt = etime (clock, t0);
fprintf ('time in seconds = %2.1f\n\n', ttt)
if e<error
fprintf ('The error limit is satisfied\n\n')
elseif ttt>tmax
fprintf ('Time limit exceeded\n\n')
end
p = p + 1;
end
Figure (1)
plot 3 (y(1, :), y(2, :), y(3, :));
hold on
a = 0.1; b = 0.1; c = 14;A1 = c + sqrt (c^2-4 * a * b); A2 = -c-sqrt (c^2-4 * a * b);
A3 = c-sqrt (c^2-4 * a * b); A4 = -c + sqrt (c^2-4 * a * b);
plot 3 (A1/2, A2/(2 * a), A1/(2 * a), 'r*') % plots fixed points
hold on
plot 3 (A3/2, A4/(2 * a), A3/(2 * a), 'r*') % plots fixed points
axis (plotrange);
grid on;
end
fprintf ('The total time in seconds is %2.1f\n', toc)