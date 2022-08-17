function[xdot]=HROV3DOF_H(u,v,r,x,y,apsi,tauu,taur)
// 3DOF of hybrid underwater ROV
// M vdot+C(v)v+D(v)v+g(n)=tau
// coefficients in International Systems
// Demo: >> Vdot=HROV3DOF_H(0.5, 0.01, 0.001, 10, 2)

global m Iz xg yg yp1 yp2 xudot yvdot nrdot xuu yvv nrr

// Load the parameters
exec('quadrotorParameters.sce', -1);

V=[u v r]';
tau=[tauu 0 taur]';
// MRB due to the rigid body
MRB=[  m     0    -m*yg;
       0     m     m*xg;
     -m*yg  m*xg    Iz];
// MA due to the added mass
MA=[-xudot      0      0;
       0     -yvdot    0;
       0        0   -nrdot];
M=MRB+MA;
// CRB due to the rigid body 
CRB=[    0         0   -m*(xg*r+v);
         0         0        m*u;
      m*(xg*r+v) -m*u        0       ];
// CV coriolis forces due to added mass effects
CA=[     0              0        yvdot*v;
         0              0       -xudot*u;
    -yvdot*v         xudot*u        0       ];
CV=CRB+CA;
// Damping matrix
DV=[-xuu*abs(u)     0          0;
         0       -yvv*abs(v)     0;
         0           0      -nrr*abs(r)];
g=[0 0 0]';
// from M vdot+C(v)v+D(v)v+g(n)=tau
// 
Minv=inv(M);
Vdot=Minv*(-CV*V-DV*V-g+tau);
// here, the kinematic equations are added
// xdot=[Vdot ndot]', where:
// V=[u v r]' and
// n=[x y apsi]'.
ndot=[cos(apsi)*u;
      sin(apsi)*u+v;
      r];
xdot=[Vdot;
      ndot];
endfunction 