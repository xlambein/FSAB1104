% Exercice 8 (APE2)

% 8.1

X = sin((0:3)*pi/6);
Y = cos((0:3)*pi/6);

x = linspace(X(1), X(end), 100);
plot(X, Y, 'r*', x, spline(X, Y, x));

% 8.2

T = (0:3)*pi/6;

t = linspace(T(1), T(end), 100);
plot(spline(T, X, t), spline(T, Y, t));

% 8.3

T = (0:12)*pi/6;
X = sin(T);
Y = cos(T);

t = linspace(T(1), T(end), 100);
plot(spline(T, X, t), spline(T, Y, t));