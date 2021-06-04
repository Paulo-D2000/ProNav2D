import matplotlib.pyplot as plt
import numpy as np
import math

#Based on https://www.youtube.com/watch?v=CS9BL-0mdh8
#----------------------------------------------------------------------#

def nlinpronav_sim(t, y, HE_rad, Np, aT, VM, PN_type):
    #Pointers (indicies)
    sel_beta = 0;
    sel_RT1 = 1;
    sel_RT2 = 2;
    sel_RM1 = 3;
    sel_RM2 = 4;
    sel_VT1 = 5;
    sel_VT2 = 6;
    sel_VM1 = 7;
    sel_VM2 = 8;

    # alloc output vector
    dy = np.zeros(9);

    # target velocity magnitude
    VT = math.sqrt( pow(y[sel_VT1],2) + pow(y[sel_VT2],2) );

    # relative positions and velocities
    RTM1 = y[sel_RT1] - y[sel_RM1];
    RTM2 = y[sel_RT2] - y[sel_RM2];
    VTM1 = y[sel_VT1] - y[sel_VM1];
    VTM2 = y[sel_VT2] - y[sel_VM2];

    # relative distance
    RTM = math.sqrt( pow(RTM1, 2) + pow(RTM2, 2) );

    # line of sight angle and time derivative
    lambda_ = math.atan2( RTM2, RTM1 );
    lambda_dot = (RTM1*VTM2 - RTM2*VTM1)/pow(RTM, 2);
    
    # closing velocity
    VC = -(RTM1*VTM1 + RTM2*VTM2)/RTM;

    # DE  RHS computations y = [beta, RTx, RTz, RMx, RMz, VTx, VTz, VMx, VMz]
    dy[0] = aT/VT;
    dy[1] = VT*math.cos(y[sel_beta]);
    dy[2] = VT*math.sin(y[sel_beta]);
    dy[3] = y[sel_VM1];
    dy[4] = y[sel_VM2];
    dy[5] = aT*math.sin(y[sel_beta]);
    dy[6] = aT*math.cos(y[sel_beta]);

    # compute LHS of pursuer acceleration equation depending on PN type
    if(PN_type == 'True'):
        nc = Np*VC*lambda_dot;
        dy[7] = -nc*math.sin(lambda_);
        dy[8] =  nc*math.cos(lambda_);
    elif(PN_type == 'Pure'):
        Heading_pursuer = math.atan2(y[sel_VM2], y[sel_VM1]);
        nc = Np*VM*lambda_dot;
        dy[7] = -nc*math.sin(Heading_pursuer);
        dy[8] =  nc*math.cos(Heading_pursuer);
    else:
        print("ERROR: PN_type must be [string] with name = 'Pure' or 'True' ");
        exit(1);
    return dy;

# Program
#----------------------------------------------------------------------#

# ProNav Type
PN_type = 'Pure';
#PN_type = 'True';

# ProNav & Engagement Params
aT = 0;                 #Target Acceleration
HE_rad = -40*np.pi/180; #40 deg heading error
Np = 3;                 #ProNav Gain
tf = 22;                #Simulation time
h = 1e-2;               #0.01s timestep

# Engagement Initial Cond.
beta_rad = 0;           
RTx      = 40000;       #Target horizontal  range [ft]
RTz      = 10000;       #Target vertical    range [ft]
RMx      = 0;           #Missile horizontal range [ft]
RMz      = 10000;       #Missile vertical   range [ft]
VM       = 3000;        #Missile            speed [ft/s]
VT       = 1000;        #Target             speed [ft/s] 

# Initial Calcs.
#----------------------------------------------------------------------#

# resolve target vel. components inertial CS
VTx = VT*math.cos(beta_rad);
VTz = VT*math.sin(beta_rad);

# relative positions and velocities
RTMx = RTx - RMx;
RTMz = RTz - RMz;

# range
#RTM = math.sqrt(pow(RTMx,2) + pow(RTMz, 2));

# line of sight angle
lambda_ = math.atan2( RTMz, RTMx );

# missile lead angle
L = math.sin( VT*math.sin( beta_rad + lambda_ )/VM );

# missile velocity components
VMx = VM*math.cos(lambda_ * L + HE_rad);
VMz = VM*math.sin(lambda_ * L + HE_rad);

# Initial condition vector
#----------------------------------------------------------------------#

# IC vector
y0 = [beta_rad, RTx, RTz, RMx, RMz, VTx, VTz, VMx, VMz];

# Simulation - integrate nonlinear 2-D engagement with (RK4 ? )
#----------------------------------------------------------------------#

# Discretized time line
t = np.arange(0, tf, h);
nt = len(t)

# Prealloc solution matrix
yout = np.zeros((len(y0),nt));

# Assing initial condition
yout[:,0] = y0;

# Integrate with const dT, h [s]
y = y0;
for j in range(0,nt-1):
    s1 = nlinpronav_sim(t[j], y, HE_rad, Np, aT, VM, PN_type);
    s2 = nlinpronav_sim(t[j]+h/2, y+h*s1/2, HE_rad, Np, aT, VM, PN_type);
    s3 = nlinpronav_sim(t[j]+h/2, y+h*s2/2, HE_rad, Np, aT, VM, PN_type);
    s4 = nlinpronav_sim(t[j]+h, y+h*s3/2, HE_rad, Np, aT, VM, PN_type);
    y = y + h*(s1 + 2*s2 + 2*s2 + s4)/6;
    yout[:,j+1] = y;
y = yout;

# Post process & Display
#----------------------------------------------------------------------#

# Pointers (indicies)

RT1 = 1;
RT2 = 2;
RM1 = 3;
RM2 = 4;
VT1 = 5;
VT2 = 6;
VM1 = 7;
VM2 = 8;

# Plot range(Z) x time
plt.figure(1)
plt.title("Z Range(ft) x Time(s)")
plt.plot(t,y[RT2])
plt.plot(t,y[RM2])

# Plot range x time
plt.figure(2)
plt.title("Range(ft) x Time(s) ")
plt.plot(t,np.sqrt(np.power(y[RT1]-y[RM1],2)+np.power(y[RT2]-y[RM2],2)))
plt.show()
