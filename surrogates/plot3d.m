T = csvread('../data/DOE_p30_t99.csv', 1, 0);x = T(:, 1);y = T(:, 2);z = T(:, 3);plot3(x, y, z, '.');