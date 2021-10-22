T = csvread('../data/DOE_p30_t99.csv', 1, 0);

% temperature,pressure,hconv,dP_over_l,cp_mol,rho,Re,friction_factor

t = T(:, 1);
p = T(:, 2);
h = T(:, 3);
dp = T(:, 4);
cp = T(:, 5);
rho = T(:, 6);
re = T(:, 7);
f = T(:, 8);


figure();
plot3(t, 1e-6 * p, h, '.');
xlabel('Temperature (K)');
ylabel('Pressure (MPa)');
zlabel('H (w/m**2 K)');

figure();
plot3(cp, 1e-6 * p, h, '.');
xlabel('Cp (J / mol * K)');
ylabel('Pressure (MPa)');
zlabel('H (w/m**2 K)');
