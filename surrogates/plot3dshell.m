

T = csvread('../data/DOE_air_p20_t20.csv', 1, 0);

% temperature,pressure,Re,hconv,v_in,v_max

t = T(:, 1);
p = T(:, 2);
re = T(:, 3);
h = T(:, 4);
v_in = T(:, 5);
v_max = T(:, 6);



figure();
plot3(t, p, h, '.');
xlabel('Temperature');
ylabel('Pressure');
zlabel('hconv');

figure();
plot3(t, p, re, '.');
xlabel('Temperature');
ylabel('Pressure');
zlabel('Re');

figure();
plot3(t, p, v_in, '.');
xlabel('Temperature');
ylabel('Pressure');
zlabel('V in');

figure();
plot3(t, p, v_max, '.');
xlabel('Temperature');
ylabel('Pressure');
zlabel('V max');
