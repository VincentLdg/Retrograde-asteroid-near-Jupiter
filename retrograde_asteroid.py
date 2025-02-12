# *******CIRI : Population of asteroids *******
#DUBOIS ETHAN
#LADUGUIE VINCENT

# %% Functions
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


def Derivee(U,t,mat_masse,mat_obj,planet):
    """
    This function return the acceleration of the asteroid at a given time t
    """
    Gsi=6.67384*10**-11
    day_sec=86400
    solarmass_kg=1.9891*10**30
    unitastro_m=1.496*10**11
    G=(Gsi*(1/unitastro_m)**3)*(solarmass_kg)*(day_sec)**2
    
    #First we calculate the acceleration due to the force exerted by the sun 
    Vx=U[1]
    Vy=U[3]
    Vz=U[5]
    dist_sun=((U[0])**2+(U[2])**2+(U[4])**2)**(3/2)
    acc_x=-G*mat_masse[0]*(U[0])/dist_sun
    acc_y=-G*mat_masse[0]*(U[2])/dist_sun
    acc_z=-G*mat_masse[0]*(U[4])/dist_sun
    
    #Then we do the same if there are other object
    for k in range(1,np.shape(mat_masse)[0]) : 
        autre_obj=mat_obj[k*6:k*6+6]
        sun=mat_obj[:6]
        dist_autre_obj=((U[0]-autre_obj[0])**2+(U[2]-autre_obj[2])**2+(U[4]-autre_obj[4])**2)**(3/2)
        dist_autre_obj_sol=((sun[0]-autre_obj[0])**2+(sun[2]-autre_obj[2])**2+(sun[4]-autre_obj[4])**2)**(3/2)
        acc_x= acc_x -G*mat_masse[k]*((U[0]-autre_obj[0])/dist_autre_obj+(autre_obj[0]/dist_autre_obj_sol))
        acc_y= acc_y -G*mat_masse[k]*((U[2]-autre_obj[2])/dist_autre_obj+(autre_obj[2]/dist_autre_obj_sol))
        acc_z= acc_z -G*mat_masse[k]*((U[4]-autre_obj[4])/dist_autre_obj+(autre_obj[4]/dist_autre_obj_sol))        
        
    return np.array([Vx , acc_x , Vy , acc_y, Vz , acc_z])

def RK4(U,t, mat_masse,mat_obj,h):
    """
    This function use the algorithm Runge Kutta 4 to calcultate the position at time t+1 using the previous position
    """
    k1=Derivee(U,t,mat_masse,mat_obj,planet)
    k2=Derivee(U+h*k1/2,t+0.5*h,mat_masse,mat_obj,planet)
    k3=Derivee(U+h*k2/2,t+0.5*h,mat_masse,mat_obj,planet)
    k4=Derivee(U+h*k3,t+h,mat_masse,mat_obj,planet)
    
    U=U+h*(k1+2*k2+2*k3+k4)/6
    return U

def Jupiter_circ(t):
    """
    Return the position of Jupiter at time t if we consider that the orbit is circular
    """
    Tj=4330.595 
    x=5.2*np.cos(2*np.pi*t/Tj)
    y=5.2*np.sin(2*np.pi*t/Tj) 
    return  np.array([x,0,y,0,0,0])

def Resolution(U0,t_start,t_end,h,mat_masse,planet):
    """
    This function calculates for a given duration :
        • U: the position of the asteroid
        • A: the semi major axis
        • E: the eccentricity
        • I: the inclination 
        • P: the position of other object
    """
    t=np.arange(t_start,t_end + 1,h)
    nb_pts=len(t)
    print(nb_pts)
    U=np.zeros((6,nb_pts))
    A=np.zeros((nb_pts))
    E=np.zeros((nb_pts))
    I=np.zeros((nb_pts))
    P=np.zeros((3*len(planet),nb_pts))
    
    U[:,0]=np.array(U0)
    pos=np.array([U[0,0],U[2,0],U[4,0]])
    vit=np.array([U[1,0],U[3,0],U[5,0]])
    A[0]=SemiMajorAxis(pos,vit)
    E[0]=Eccentricity(pos,vit)
    I[0]=Inclinaison(pos,vit)

    for j in range(len(planet)):
        r=planet[j](t_start)
        P[0+3*j,0]=r[0]
        P[1+3*j,0]=r[2]
        P[2+3*j,0]=r[4]
    
    for i in tqdm(range(1,nb_pts)):
        mat_obj=np.array([0,0,0,0,0,0])
        for j in range(len(planet)):
            r=planet[j](t[i])
            P[0+3*j,i]=r[0]
            P[1+3*j,i]=r[2]
            P[2+3*j,i]=r[4]
            mat_obj=np.hstack((mat_obj,r))

        U[:,i]=np.array(RK4(U[:,i-1],t[i],mat_masse,mat_obj,h))

        pos=np.array([U[0,i],U[2,i],U[4,i]])
        vit=np.array([U[1,i],U[3,i],U[5,i]])
        A[i]=SemiMajorAxis(pos,vit)
        E[i]=Eccentricity(pos,vit)
        I[i]=Inclinaison(pos,vit)

    return U,A,E,I,P

def SemiMajorAxis(pos,vit):
     Gsi=6.67384*10**-11
     day_sec=86400
     solarmass_kg=1.9891*10**30
     unitastro_m=1.496*10**11
     G=(Gsi*(1/unitastro_m)**3)*(solarmass_kg)*(day_sec)**2
     mu=G*1
     p=np.linalg.norm(pos)
     v=np.linalg.norm(vit)
     a=((2/p)-v**2/mu)**(-1)
     return a

def Eccentricity(pos,vit):
     Gsi=6.67384*10**-11
     day_sec=86400
     solarmass_kg=1.9891*10**30
     unitastro_m=1.496*10**11
     G=(Gsi*(1/unitastro_m)**3)*(solarmass_kg)*(day_sec)**2
     mu=G*1
     return np.linalg.norm((np.cross(vit,np.cross(pos,vit))/mu)-pos/np.linalg.norm(pos))

def Inclinaison(pos,vit):
    p=np.linalg.norm(pos)
    v=np.linalg.norm(vit)
    K=np.cross(pos,vit)/(p*v)
    return np.rad2deg(np.arccos(K[2]))


# Functions to convert orbital parameters in cartesian coordonates
def mean_motion(G,a):
    return np.sqrt(G/a**3)

def mean_anomaly(M0,n,t):
    return M0+n*t

def Ecc(E,M,e):
    return E-e*np.sin(E)-M

def Ecc_dot(E,e):
    return 1-e*np.cos(E)

def eccentric_anomaly(e,M0,M):
    E1=M0
    epsilon=10**-4
    E2=E1-Ecc(E1,M,e)/Ecc_dot(E1,e)
    while(E2-E1>epsilon):
        E1=E2
        E2=E1-Ecc(E1,M,e)/Ecc_dot(E1,e)
    
    return E1

def R1(theta):
    return np.array([[1,0,0],[0,np.cos(theta),np.sin(theta)],[0,-np.sin(theta),np.cos(theta)]])
def R3(theta):
    return np.array([[np.cos(theta),np.sin(theta),0],[-np.sin(theta),np.cos(theta),0],[0,0,1]])

def conversion_to_cartesian(a,e,i,longi,w,M0,t):
    Gsi=6.67384*10**-11
    day_sec=86400
    solarmass_kg=1.9891*10**30
    unitastro_m=1.496*10**11
    G=(Gsi*(1/unitastro_m)**3)*(solarmass_kg)*(day_sec)**2

    n=mean_motion(G,a)
    M=mean_anomaly(M0,n,t)
    E=eccentric_anomaly(e,M0,M)
    r=a*(1-e*np.cos(E))
    X=a*(np.cos(E)-e)
    Y=a*np.sqrt(1-e**2)*np.sin(E)
    Xdot=-np.sin(E)*n*a**2/r
    Ydot=np.sqrt(1-e**2)*np.cos(E)*n*a**2/r
    pos=R3(-longi)@R1(-i)@R3(-w)@np.array([[X],[Y],[0]])
    vit=R3(-longi)@R1(-i)@R3(-w)@np.array([[Xdot],[Ydot],[0]])

    return np.array([pos[0,0],vit[0,0],pos[1,0],vit[1,0],pos[2,0],vit[2,0]])


def Jupiter_ell(t):
    """
    Convert orbital paramters into cartesian coordinates for a given time t
    Epoch : 4 november 2013
    """
    a=5.202575
    e=0.048908
    i=1.3038*np.pi/180
    longi=100.5145*np.pi/180
    w=273.8752*np.pi/180
    M0=80.0392*np.pi/180
    
    return conversion_to_cartesian(a,e,i,longi,w,M0,t)

# This functions return the position and velocity of Mars and Earth
# Epoch J2000 so we shift the time t with the number of day between J2000 and the 4 november 2013 (5056 days)
def Mars(t):
    a=1.52366231
    e=0.09341233
    i=1.85061*np.pi/180
    longi=49.57854*np.pi/180
    w=286.46230*np.pi/180
    M0=19.41248*np.pi/180
    return conversion_to_cartesian(a,e,i,longi,w,M0,t+5056)

def Earth(t):
    a=1.0000010178
    e=0.0167086
    i=1.578690*np.pi/180
    longi=174.9*np.pi/180
    w=288.1*np.pi/180
    M0=357.51716*np.pi/180
    return conversion_to_cartesian(a,e,i,longi,w,M0,t+5056)

def Saturn(t) :
    a = 9.5820172
    e = 0.05555
    i = 2.485 * np.pi / 180
    longi = 113.715 * np.pi / 180
    w = 339.392 * np.pi / 180
    M0 = 49.944 * np.pi / 180
    return conversion_to_cartesian(a, e, i, longi, w, M0, t + 5056)

def Plot_trajectory(U,P,planet,name):
    c=['orange','red','green','magenta','black']
    plt.figure()
    plt.scatter(0,0,color='yellow',s=100)
    plt.scatter(U[0,-1],U[2,-1],color='blue',s=50,zorder=2)
    if len(planet)>0:
        for j in range(len(planet)):
            plt.scatter(P[0+3*j,-1],P[1+3*j,-1],color=c[j],s=50,zorder=2)
    plt.legend(['Sun','Asteroid']+name)

    plt.plot(U[0,:],U[2,:],'b-.',linewidth=1,zorder=1)
    if len(planet)>0:
        for j in range(len(planet)):
            plt.plot(P[0+3*j,:],P[1+3*j,:],color=c[j],linestyle='-.',zorder=2)
    plt.title("Trajectory of the asteroid")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis("equal")
    plt.show()

def Plot_trajectory3d(U,P,planet,name):
    ax = plt.figure().add_subplot(projection='3d')
    c=['orange','red','green','magenta','black']

    ax.scatter(0,0,0,color='yellow',s=100)
    ax.scatter(U[0,-1],U[2,-1],U[4,-1],color='blue',s=50,zorder=2)
    if len(planet)>0:
        for j in range(len(planet)):
            ax.scatter(P[0+3*j,-1],P[1+3*j,-1],P[2+3*j,-1],color=c[j],s=50,zorder=2)
    ax.legend(['Sun','Asteroid']+name)

    ax.plot(U[0,:],U[2,:],U[4,:],'b-.',linewidth=1,zorder=1)
    if len(planet)>0:
        for j in range(len(planet)):
            ax.plot(P[0+3*j,:],P[1+3*j,:],P[2+3*j,:],color=c[j],linestyle='-.',zorder=2)
    ax.set_title("Trajectory of the asteroid")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

def Plot_orbital_param(A,E,I,t_start,t_end,h):

    print(t_start)
    time_vector = np.arange(t_start,t_end + 1,h)/365.25
    fig,(ax1,ax2,ax3)=plt.subplots(3)
    fig.tight_layout(pad=3.3)
    ax1.plot(time_vector,A,'k')
    ax1.set_title("Semi-major axis")
    ax1.set_xlabel("years")
    ax1.set_ylabel("a [AU]")
    ax1.set_xlim([0,2000])

    ax2.plot(time_vector,E,'r')
    ax2.set_title("Eccentricity")
    ax2.set_xlabel("years")
    ax2.set_ylabel("e")
    ax2.set_xlim([0,2000])

    ax3.plot(time_vector,I,'b')
    ax3.set_title("Inclination")
    ax3.set_xlabel("years")
    ax3.set_ylabel("i [°]")
    ax3.set_xlim([0,2000])
    plt.show()

#%% Test with an asteroid at 2 UA from the Sun
h=1   #step 1 day

nb_its=36500   # duration of the simulation 
nb_pts=int(np.round(nb_its/h))
Gsi=6.67384*10**-11
day_sec=86400
solarmass_kg=1.9891*10**30
unitastro_m=1.496*10**11
G=(Gsi*(1/unitastro_m)**3)*(solarmass_kg)*(day_sec)**2


ms=1

m=0
x=2
Vx=0
y=0
Vy=np.sqrt(G/x)
z=0
Vz=0
U0=[x,Vx,y,Vy,z,Vz]
t0=0

planet=[]
mat_masse=np.array([ms])

# %%
U,A,E,I,P=Resolution(U0,t0,mat_masse,h,nb_pts,planet)

# %%
Plot_trajectory(U,P,planet,[])

# %%

fig,(ax1,ax2,ax3)=plt.subplots(3)
fig.tight_layout(pad=3.3)
ax1.plot(np.arange(t0,nb_its,h)/365.25,A,'k')
ax1.set_title("Semimajor axis over time")
ax1.set_xlabel("days")
ax1.set_ylabel("a")
ax1.set_ylim(1,3)

ax2.plot(np.arange(t0,nb_its,h)/365.25,E,'r')
ax2.set_title("Eccentricity over time")
ax2.set_xlabel("days")
ax2.set_ylabel("e")
ax2.set_ylim(-1,1)

ax3.plot(np.arange(t0,nb_its,h)/365.25,I,'b')
ax3.set_title("Inclination over time")
ax3.set_xlabel("days")
ax3.set_ylabel("i")
ax3.set_ylim(-1,1)
plt.show()


# %%
Gsi=6.67384*10**-11
day_sec=86400
solarmass_kg=1.9891*10**30
unitastro_m=1.496*10**11
G=(Gsi*(1/unitastro_m)**3)*(solarmass_kg)*(day_sec)**2

ms=1.988e30
mj=1.898e27

h=1
nb_its=36500
nb_pts=int(np.round(nb_its/h))

# asteroid's initial position and velocity
m=0
x=3.27
Vx=0
y=0
Vy=np.sqrt(G/x)
z=0
Vz=0
U0=[x,Vx,y,Vy,z,Vz]

m1=1
m2=mj/ms

mat_masse=np.array([m1 , m2])
planet=[Jupiter_circ]

t0=0

# %%
U,A,E,I,P=Resolution(U0,t0,mat_masse,h,nb_its,planet)

# %%
Plot_trajectory(U,P,planet,["Jupiter"])

# %%
Plot_orbital_param(A,E,I,nb_its,h,t0)

# %%
Gsi=6.67384*10**-11
day_sec=86400
solarmass_kg=1.9891*10**30
unitastro_m=1.496*10**11
G=(Gsi*(1/unitastro_m)**3)*(solarmass_kg)*(day_sec)**2

au=150e6
rjs=780e6
ms=1.988e30
mj=1.898e27
muj=G*mj/ms

nb_corps =3
h=1
nb_its=36500
nb_pts=int(np.round(nb_its/h))

a=3.27
e=0
i=0

m=0
x=a
Vx=0
y=0
Vy=np.sqrt(G/a)
z=0
Vz=0
U0=[x,Vx,y,Vy,z,Vz]

m1=1
m2=mj/ms
t0=0
mat_masse=np.array([m1 , m2])
planet=[Jupiter_ell]


# %%
U,A,E,I,P=Resolution(U0,t0,mat_masse,h,nb_pts,planet)

# %%
Plot_trajectory(U,P,planet,["Jupiter"])

# %%
Plot_orbital_param(A,E,I,nb_its,h,t0)


# %%
Gsi=6.67384*10**-11
day_sec=86400
solarmass_kg=1.9891*10**30
unitastro_m=1.496*10**11
G=(Gsi*(1/unitastro_m)**3)*(solarmass_kg)*(day_sec)**2

ms=1.988e30
mj=1.898e27

t_start= 0 * 365.25 # start year
t_end = 9000 * 365.25
h=10 # step time in days

# Asteroid's orbital parameters on the 4 of november 2013
a=5.454
e=0.3896
i=108.358*np.pi/180
longi=276.509*np.pi/180
w=226.107*np.pi/180
M0=146.88*np.pi/180
r=conversion_to_cartesian(a,e,i,longi,w,M0,0)
U0=[list(r)]

m1=1
m2=mj/ms

mat_masse=np.array([m1 , m2])

planet=[Jupiter_ell,Mars,Earth]

# %%
U,A,E,I,P=Resolution(U0,t_start,t_end,h,mat_masse,planet)

# %%
Plot_trajectory(U,P,planet,["Jupiter","Mars","Earth"])

# %%
Plot_orbital_param(A,E,I,t_start,t_end,h)

# %%
Plot_trajectory3d(U,P,planet,["Jupiter","Mars","Earth"])


#%% Incertitudes
Gsi=6.67384*10**-11
day_sec=86400
solarmass_kg=1.9891*10**30
unitastro_m=1.496*10**11
G=(Gsi*(1/unitastro_m)**3)*(solarmass_kg)*(day_sec)**2

ms=1.988e30
mj=1.898e27

t0=-1000
h=90
nb_its=(10000-t0)*365
nb_pts=int(np.round(nb_its/h))

# Incertitudes
ea=0.0156
ee=0.00170
ei=0.0261
elongi=0.00114
ew=0.0501
eM0=0.604

#asteroid 1
a=5.454-ea
e=0.3896-ee
i=(108.358-ei)*np.pi/180
longi=(276.509-elongi)*np.pi/180
w=(226.107-ew)*np.pi/180
M0=(146.88-eM0)*np.pi/180
r=conversion_to_cartesian(a,e,i,longi,w,M0,0)
U01=[list(r)]

#asteroid 2
a=5.454
e=0.3896
i=108.358*np.pi/180
longi=276.509*np.pi/180
w=226.107*np.pi/180
M0=146.88*np.pi/180
r=conversion_to_cartesian(a,e,i,longi,w,M0,0)
U02=[list(r)]

#asteroide 3
a=5.454+ea
e=0.3896+ee
i=(108.358+ei)*np.pi/180
longi=(276.509+elongi)*np.pi/180
w=(226.107+ew)*np.pi/180
M0=(146.88+eM0)*np.pi/180
r=conversion_to_cartesian(a,e,i,longi,w,M0,0)
U03=[list(r)]

m1=1
m2=mj/ms

mat_masse=np.array([m1 , m2])

planet=[Jupiter_ell,Mars,Earth]

#%%
U1,A1,E1,I1,P=Resolution(U01,t0,mat_masse,h,nb_its,planet)
U2,A2,E2,I2,P=Resolution(U02,t0,mat_masse,h,nb_its,planet)
U3,A3,E3,I3,P=Resolution(U03,t0,mat_masse,h,nb_its,planet)

#%% Trajectory
plt.figure()
plt.scatter(0,0,color='yellow',s=100)
plt.scatter(U1[0,-1],U1[2,-1],color='blue',s=50,zorder=2)
plt.scatter(U2[0,-1],U2[2,-1],color='red',s=50,zorder=2)
plt.scatter(U3[0,-1],U3[2,-1],color='green',s=50,zorder=2)

if len(planet)>0:
    for j in range(len(planet)):
        plt.scatter(P[0+3*j,-1],P[1+3*j,-1],s=50,zorder=2)
        
plt.legend(['Sun','Asteroid 1','Asteroid 2','Asteroid 3','Jupiter','Mars','Earth'])
plt.plot(U1[0,:],U1[2,:],'b-.',linewidth=1,zorder=1)
plt.plot(U2[0,:],U2[2,:],'r-.',linewidth=1,zorder=1)
plt.plot(U3[0,:],U3[2,:],'g-.',linewidth=1,zorder=1)

if len(planet)>0:
    for j in range(len(planet)):
        plt.plot(P[0+3*j,:],P[1+3*j,:],linestyle='-.',zorder=2)
plt.title("Trajectory of the asteroid")
plt.xlabel("X")
plt.ylabel("Y")
plt.axis("equal")
plt.show()

#%% Plot orbital parameters 

fig,(ax1,ax2,ax3)=plt.subplots(3)
fig.tight_layout(pad=3.3)
ax1.plot(np.arange(t0,nb_its,h)/365.25,A1,'b')
ax1.plot(np.arange(t0,nb_its,h)/365.25,A2,'r')
ax1.plot(np.arange(t0,nb_its,h)/365.25,A3,'g')
ax1.legend(['Asteroid 1','Asteroid 2','Asteroid 3'])
ax1.set_title("Semimajor axis over time")
ax1.set_xlabel("years")
ax1.set_ylabel("a")

ax2.plot(np.arange(t0,nb_its,h)/365.25,E1,'b')
ax2.plot(np.arange(t0,nb_its,h)/365.25,E2,'r')
ax2.plot(np.arange(t0,nb_its,h)/365.25,E3,'g')
ax2.legend(['Asteroid 1','Asteroid 2','Asteroid 3'])
ax2.set_title("Eccentricity over time")
ax2.set_xlabel("years")
ax2.set_ylabel("e")

ax3.plot(np.arange(t0,nb_its,h)/365.25,I1,'b')
ax3.plot(np.arange(t0,nb_its,h)/365.25,I2,'r')
ax3.plot(np.arange(t0,nb_its,h)/365.25,I3,'g')
ax3.legend(['Asteroid 1','Asteroid 2','Asteroid 3'])
ax3.set_title("Inclination over time")
ax3.set_xlabel("years")
ax3.set_ylabel("i")
plt.show()


#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime
# import matplotlib as mpl
# mpl.use('Qt5Agg')
# Simulation parameters
start_date = datetime.datetime(2013, 11, 4)
time_step_days = h

# Initialize figure
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-6, 6)
ax.set_ylim(-6, 6)
ax.set_xlabel("X (AU)")
ax.set_ylabel("Y (AU)")

# Colors
colors = {
    "earth": "#3b7ddd",  # Blue
    "mars": "#c1440e",   # Red-orange
    "jupiter": "#daa520",  # Golden-brown
    "saturn": "#f4a460",  # Sandy yellow
    "asteroid": "#555555"  # Dark gray
}

# Plot elements
asteroid_trail, = ax.plot([], [], '-', color=colors["asteroid"], alpha=0.6, label="Asteroid 2007 W266")
earth_trail, = ax.plot([], [], '-', color=colors["earth"], alpha=0.6, label="Earth")
mars_trail, = ax.plot([], [], '-', color=colors["mars"], alpha=0.6, label="Mars")
jupiter_trail, = ax.plot([], [], '-', color=colors["jupiter"], alpha=0.6, label="Jupiter")
#saturn_trail, = ax.plot([], [], '-', color=colors["saturn"], alpha=0.6, label="Saturn")

asteroid_dot, = ax.plot([], [], 'o', color=colors["asteroid"], markersize=6)
earth_dot, = ax.plot([], [], 'o', color=colors["earth"], markersize=6)
mars_dot, = ax.plot([], [], 'o', color=colors["mars"], markersize=6)
jupiter_dot, = ax.plot([], [], 'o', color=colors["jupiter"], markersize=6)
#saturn_dot, = ax.plot([], [], 'o', color=colors["saturn"], markersize=6)
sun_dot = ax.scatter(0,0,color="y",marker="o", s=60, label="Sun")
# Title
title = ax.set_title("")
plt.xlim([-7.5,7.5])
plt.ylim([-7.5,7.5])
# Update function for animation
def update(frame):
    # Update trails
    frame=int(frame*15)
    print(frame)
    jupiter_trail.set_data(P[0, :frame], P[1, :frame])
    mars_trail.set_data(P[3, :frame], P[4, :frame])
    earth_trail.set_data(P[6, :frame], P[7, :frame])
    #saturn_trail.set_data(P[9, :frame], P[10, :frame])
    asteroid_trail.set_data(U[0, :frame], U[2, :frame])

    jupiter_dot.set_data(P[0, frame], P[1, frame])
    mars_dot.set_data(P[3, frame], P[4, frame])
    earth_dot.set_data(P[6, frame], P[7, frame])
    #saturn_dot.set_data(P[9, frame], P[10, frame])
    asteroid_dot.set_data(U[0, frame], U[2, frame])
    # # Update moving points
    # asteroid_dot.set_data(positions_asteroid[0, frame], positions_asteroid[1, frame])
    # earth_dot.set_data(earth_positions[0, frame], earth_positions[1, frame])
    # mars_dot.set_data(mars_positions[0, frame], mars_positions[1, frame])
    # jupiter_dot.set_data(jupiter_positions[0, frame], jupiter_positions[1, frame])
    # saturn_dot.set_data(saturn_positions[0, frame], saturn_positions[1, frame])
    # Update title with month and year
    current_date = start_date + datetime.timedelta(days=frame * time_step_days)
    title.set_text(current_date.strftime("%B %Y"))

    return (asteroid_trail, earth_trail, mars_trail, jupiter_trail,
            asteroid_dot, earth_dot, mars_dot, jupiter_dot, title)

# Create animation
ani = animation.FuncAnimation(fig, update, frames=300, interval=0.1, repeat=True)

# Show animation
plt.legend()
plt.grid()
#plt.show()

# Save animation as GIF (PillowWriter)
gif_filename = "orbit_animation2.gif"
ani.save(gif_filename, writer=animation.PillowWriter(fps=240))

print(f"Animation saved as {gif_filename}")

#%%
dist_jup = np.linalg.norm(U[[0,2,4],:]-P[:3,:], axis=0)

plt.figure()
plt.plot(np.arange(t0,nb_its,h)/365.25, dist_jup)
plt.show()