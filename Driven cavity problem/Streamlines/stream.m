x = load('x.dat');
y = load('y.dat');

%% Re = 100

u100 = load('Resultats100.dat');
v100 = load('Resvltats100.dat');

figure(100);
streamslice(x,y,u100,v100);
axis equal
xlim([0 1]);
ylim([0 1]);

%% Re = 400

u400 = load('Resultats400.dat');
v400 = load('Resvltats400.dat');

figure(400);
streamslice(x,y,u400,v400);
axis equal
xlim([0 1]);
ylim([0 1]);

%% Re = 1000

u1000 = load('Resultats1000.dat');
v1000 = load('Resvltats1000.dat');

figure(1000);
streamslice(x,y,u1000,v1000);
axis equal
xlim([0 1]);
ylim([0 1]);

%% Re = 3200

u3200 = load('Resultats3200.dat');
v3200 = load('Resvltats3200.dat');

figure(3200);
streamslice(x,y,u3200,v3200);
axis equal
xlim([0 1]);
ylim([0 1]);

%% Re = 5000

u5000 = load('Resultats5000.dat');
v5000 = load('Resvltats5000.dat');

figure(5000);
streamslice(x,y,u5000,v5000);
axis equal
xlim([0 1]);
ylim([0 1]);

%% Re = 7500

u7500 = load('Resultats7500.dat');
v7500 = load('Resvltats7500.dat');

figure(7500);
streamslice(x,y,u7500,v7500);
axis equal
xlim([0 1]);
ylim([0 1]);

%% Re = 10000

u10000 = load('Resultats10000.dat');
v10000 = load('Resvltats10000.dat');

figure(10000);
streamslice(x,y,u10000,v10000);
axis equal
xlim([0 1]);
ylim([0 1]);