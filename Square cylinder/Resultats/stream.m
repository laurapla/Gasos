x = load('x.dat');
y = load('y.dat');

%% Re = 1
u1 = load('Matrixu1.dat');
v1 = load('Matrixv1.dat');

figure(1);
streamslice(x,y,u1,v1,27);
rectangle('Position',[12 3.5 1 1])
axis equal
xlim([11 17]);
ylim([2 6]);

%% Re = 3
u3 = load('Matrixu3.dat');
v3 = load('Matrixv3.dat');

figure(3);
streamslice(x,y,u3,v3,27);
rectangle('Position',[12 3.5 1 1])
axis equal
xlim([11 17]);
ylim([2 6]);

%% Re = 5
u5 = load('Matrixu5.dat');
v5 = load('Matrixv5.dat');

figure(5);
streamslice(x,y,u5,v5,28);
rectangle('Position',[12 3.5 1 1])
axis equal
xlim([11 17]);
ylim([2 6]);

%% Re = 10
u10 = load('Matrixu10.dat');
v10 = load('Matrixv10.dat');

figure(10);
streamslice(x,y,u10,v10,27);
rectangle('Position',[12 3.5 1 1])
axis equal
xlim([11 17]);
ylim([2 6]);

%% Re = 30
u30 = load('Matrixu30.dat');
v30 = load('Matrixv30.dat');

figure(30);
streamslice(x,y,u30,v30,25);
rectangle('Position',[12 3.5 1 1])
axis equal
xlim([11 17]);
ylim([2 6]);

%% Re = 50
u50 = load('Matrixu50.dat');
v50 = load('Matrixv50.dat');

figure(50);
streamslice(x,y,u50,v50,27);
rectangle('Position',[12 3.5 1 1])
axis equal
xlim([11 17]);
ylim([2 6]);

% %% Re = 60
% u60 = load('Matrixu60.dat');
% v60 = load('Matrixv60.dat');
% 
% figure(60);
% streamslice(x,y,u60,v60,27);
% rectangle('Position',[12 3.5 1 1])
% axis equal
% xlim([0 50]);
% ylim([0 8]);
% 
% %% Re = 100
% u100 = load('Matrixu100.dat');
% v100 = load('Matrixv100.dat');
% 
% figure(100);
% streamslice(x,y,u100,v100,27);
% rectangle('Position',[12 3.5 1 1])
% axis equal
% xlim([0 50]);
% ylim([0 8]);