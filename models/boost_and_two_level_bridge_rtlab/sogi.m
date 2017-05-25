clc
Q = 2^24;
k = 0.5;
wn = 50 * 2 * pi;
delta_T = 1/8e3;
x = 2 * wn * delta_T * k
int32(Q*x)
y = wn * wn * delta_T * delta_T
int32(Q*y)
a1 = 2 * (4-y)/(x+y+4)
int32(Q*a1)
a2 = (x-y-4)/(x+y+4)
int32(Q*a2)
b0 = x/(x+y+4)
int32(Q*b0)
b2 = -1*x/(x+y+4)
int32(Q*b2)
qb1 = 2*k*y/(x+y+4)
int32(Q*qb1)
qb0 = k*y / (x+y+4)
int32(Q*qb0)
qb2 = k*y / (x+y+4)
int32(Q*qb2)

%%
kp = 0.1;
ki = 2;
pr_a1 = -2 * ((delta_T * wn)^2 - 4)/((delta_T * wn)^2 + 4);
pr_b0 = (kp + 2 * ki * delta_T / (4 + delta_T^2 * wn^2));
pr_b1 = 2 * kp * (delta_T^2  * wn^2 - 4)/(delta_T^2 * wn^2 + 4);
pr_b2 = (kp - 2 * ki * delta_T / (4 + delta_T^2 * wn^2));
sys = tf([ki 0], [1 0 wn^2])+kp + tf([ki 0], [1 0 (3*wn)^2]);
bode(1/(1+sys))
dis_sys = c2d(sys,delta_T, 'tustin')
%%
wn = 25e3/5/0.707;
td = 1/0.707/wn;
kp = 2 * 0.707 * wn
ki = wn^2
kp = 166.6;
ki= 27755.55;
b0_lf = 0.5 * (2*kp + delta_T*ki);
b1_lf = -0.5 * (2*kp - delta_T*ki);

