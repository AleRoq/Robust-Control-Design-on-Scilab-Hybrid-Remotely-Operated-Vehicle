function[xdot]=hrov3dof(u,v,r,x,y,apsi,tauu,taur)
// 3DOF model of an hybrid remotely operated
// vehicle (HROV)
// Mvdot+C(v)v+D(v)v+g(n)=tau
// SNAME notation

// Load the HROV parameters
exec('hrov3dofParameters.sce', -1);

// V is the velocity vector
V=[u v r]';
// n is the inertial frame vector
n=[x y apsi]';
// tau is the control vector
tau=[tauu 0 taur]';
// MRB is the matrix due to the rigid body mass
MRB=[  m     0    -m*yg;
       0     m     m*xg;
     -m*yg  m*xg    Iz];
// MA is the matrix due to the added mass
MA=[-xudot      0      0;
       0     -yvdot    0;
       0        0   -nrdot];
M=MRB+MA;
// CRB is the matrix due to the rigid body 
CRB=[    0         0   -m*(xg*r+v);
         0         0        m*u;
      m*(xg*r+v) -m*u        0    ];
// CV is matrix due to coriolis forces due to 
// added mass
CA=[     0         0      yvdot*v;
         0         0     -xudot*u;
    -yvdot*v    xudot*u      0    ];
CV=CRB+CA;
// DV is the damping matrix
DV=[-xuu*abs(u)    0         0;
         0    -yvv*abs(v)    0;
         0         0   -nrr*abs(r)];
g=[0 0 0]';
// from the dynamic equation, 
// Mvdot+C(v)v+D(v)v+g(n)=tau
Minv=inv(M);
Vdot=Minv*(-CV*V-DV*V-g+tau);
// ndot contains the kinematic equations
ndot=[cos(apsi)*u;
      sin(apsi)*u+v;
      r];
// Finally, the dynamics and kinematics in
// xdot=[udot vdot rdot xdot ydot apsidot]'
xdot=[Vdot;
      ndot];
endfunction
