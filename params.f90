
!Time step information
&step_params
iti=0          !Time step to begin
itmax=300000   !Max number of time steps
imean=100    !Number of time steps between writing mean profiles to output log file
ihst=10     !Number of time steps between writing profiles to history file
itape=10000  !Number of time steps before opening new history file and writing full 3D volumes
itstr=1
it_his=1     !Time step to begin recording history files
it_viz=1     !Time step to begin recording viz files
i_viz=10000   !Time steps between writing to viz file
itn = 0  !Index for "u.le.cou001" files
/

!Flags for various features
&flags
ismlt=0, ifree=0, iradup=0, iupwnd=0
ibuoy=0, ifilt=0, itcut=1, isubs=0, ibrcl=0, iocean=0
method=3, idebug=1, iz_space=1, ivis0=0, ifix_dt=0, new_vis=-1
isfc = 0, 1  ! isfc = 0 sets surface flux (qstar), isfc = 1 sets surface condition (Dirichlet through tsfcc)
iDNS = 1
ifields = 1  ! if set to 1, reads the initial condition from path_ran below; otherwise calls subroutine random
ilin = 0   !1 means linear interpolation for particle, 6th order otherwise
icouple = 0
iTcouple=1
iHcouple=1  !iHcouple also controls TE couple
ievap=1
ineighbor=0
icoalesce=1
ipart_method=2   !1 = same RK3 as flow
ipartdiff = 1
isfs = 2   !isfs=1 is simple random walk, isfs=2 is Weil et al. (2004) particle subgrid
iexner = 1  !iexner = 1 applies a correction to part%Tf to account for conversion b/w potential temp and temp. (need in PBL)
/

!Grid and domain parameters
&grid_params
ncpu_s=8

!Time step
cfl = 0.4
dt_new = 0.1

!Use for DNS:
Uo = 0.0  !Sets the magnitude of the top and bottom plate velocities (Couette flow -- assumes equal and opposite)
Ttop(1) = 300.15, 95.0  !Currently this should be (temperature, relative humidity) (used in DNS)
Tbot(1) = 300.15, 100.0  !Currently this should be (temperature, relative humidity) (used in DNS)

!Use for LES:
qstar = 0.5, 0.5  !Surface fluxes of (temperature, humidity) (used for LES and DNS)
tsfcc = 305.15, 100.0    !Surface conditions of (temp, humidity) (used for LES) -- make sure tsfcc is gt than t00 for both isfc=0 or 1

surf_RH = 100.0

psurf = 101325  !Pa -- nominal surface pressure, for use in the definition of potential temp

ugcont = 0.0   !The initial u-velocity in the field
vgcont = 0.0   !The initial v-velocity

dpdx = -0.69325  !The pressure gradient for channel flow (gets put into u_geo for DNS)

zi = 0.04  !This has to do with grid stretching; make equal to zl for Couette/channel
zl = 0.04 
xl = 0.251327  !2*pi*zl
yl = 0.125664  !2*pi*zl

zw1 = 0.00008  !The first grid point
/

!Set the paths for code I/O. Must be on the scratch directory, not AFS!
&path_names
path_seed="/scratch365/drichte2/tutorial/case1/"
path_part="/scratch365/drichte2/tutorial/case1/part.le.cou000"
path_res="/scratch365/drichte2/tutorial/case1/u.le.cou000"
path_sav="/scratch365/drichte2/tutorial/case1/"
path_his="/scratch365/drichte2/tutorial/case1/"
path_viz_xy="/scratch365/drichte2/tutorial/case1/"
path_viz_xz="/scratch365/drichte2/tutorial/case1/"
path_viz_yz="/scratch365/drichte2/tutorial/case1/"
path_stuf="/scratch365/drichte2/tutorial/case1/"
path_ran="/scratch365/drichte2/tutorial/case1/u.le.cou000"
path_histog="/scratch365/drichte2/tutorial/case1/"
/


!Material and particle properties and other constants
&constants
grav = 9.81
t00 = 273.0 !Reference temp for buoyancy
fcor = 0.0  !Coriolis parameter
zo = 0.1       !Roughness height (for LES)

!Air phase:
rhoa=1.0
nuf=1.57e-5  !m^2/s
Cpa=1006.0  !J/kg-K
Pra = 0.715
Sc = 0.615
Rd = 286.9     !J/kg-K  dry air gas constant

!Particles:
tnumpart = 1e5
mult_init = 1e8
rhow=1000.0  !kg/m^3
part_grav = 0,0,-9.81 
Cpp = 4179.0  !J/kg-K
Mw = 0.018015  !kg/mol
Ru = 8.3144    !J/mol-K Universal gas constant
Ms = 0.05844  !kg/mol: molecular weight of salt
Sal = 0.0   !INITIAL salinity
Gam = 7.28e-2
Ion = 2.0
Os = 1.093

!Particle initial conditions:
radius_init=22.8e-6
Tp_init = 300.0
vp_init = 0.0, 0.0, 0.0
qf_init = 0.01

/
