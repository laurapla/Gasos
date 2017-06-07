x = load('x.dat');
y = load('y.dat');

%% Ra = 1e3

u3 = load('Matrixu3.dat');
v3 = load('Matrixv3.dat');

figure(3);
streamslice(x,y,u3,v3);
axis equal
xlim([0 1]);
ylim([0 1]);

%% Ra = 1e4

u4 = load('Matrixu4.dat');
v4 = load('Matrixv4.dat');

figure(4);
streamslice(x,y,u4,v4);
axis equal
xlim([0 1]);
ylim([0 1]);

%% Ra = 1e5

u5 = load('Matrixu5.dat');
v5 = load('Matrixv5.dat');

figure(5);
streamslice(x,y,u5,v5);
axis equal
xlim([0 1]);
ylim([0 1]);

%% Ra = 1e6

u6 = load('Matrixu6.dat');
v6 = load('Matrixv6.dat');

figure(6);
streamslice(x,y,u6,v6);
axis equal
xlim([0 1]);
ylim([0 1]);