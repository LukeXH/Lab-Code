close all
clear

t = linspace(-pi+.000000001,pi);

map1 = @(x) x.*(x >= 0) + (2*pi + x).*(x < 0);

plot(t, map1(t))
des = pi*.75;
d1 = des - t;
d2 = map1(des) - map1(t);

%%% TEST ENCODER CALCULATOR
distCalc = GenEncoderCalculator([-pi, pi]);

figure
plot(t, d1)
hold on
plot(t,d2)
plot(t, distCalc(des,t))

