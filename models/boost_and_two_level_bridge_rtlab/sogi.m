k = 0.9;
wn = 50;
delta_T = 1/25e3;
x = 2 * wn * delta_T;
y = wn * wn * delta_T * delta_T;
a1 = 2 * (4-y)/(x+y+4);
a2 = (x-y-4)/(x+y+4);
b0 = x/(x+y+4);
b2 = -1*x/(x+y+4);
qb1 = 2*k*y/(x+y+4);
qb0 = k*y / (x+y+4);
qb2 = k*y / (x+y+4);


wn = 25e3/5/0.707;
td = 1/0.707/wn;
kp = 2 * 0.707 * wn
ki = wn^2
kp = 166.6;
ki= 27755.55;
b0_lf = 0.5 * (2*kp + delta_T*ki);
b1_lf = -0.5 * (2*kp - delta_T*ki);