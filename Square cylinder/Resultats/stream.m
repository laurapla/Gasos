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

%% Re = 30
u30 = load('Matrixu30.dat');
v30 = load('Matrixv30.dat');

figure(30);
streamslice(x,y,u30,v30,25);
rectangle('Position',[12 3.5 1 1])
axis equal
xlim([11 17]);
ylim([2 6]);