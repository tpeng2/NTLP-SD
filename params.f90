
!Time step information
&step_params
iti=0          !Time step to begin
itmax=10000   !Max number of time steps
imean=100    !Number of time steps between writing mean profiles to output log file
ihst=10      !Number of time steps between writing profiles to history file
itape=5000  !Number of time steps before opening new history file and writing full 3D volumes
itstr=1
it_his=1     !Time step to begin recording history files
it_viz=1     !Time step to begin recording viz files
i_viz=1000   !Time steps between writing to viz file
itn = 0  !Index for "u.le.cou001" files
itcplbegin=001
itclsbegin=001
/

!Flags for various features
&flags
ismlt=0, ifree=0, iradup=1, iupwnd=1
ibuoy=1, ifilt=0, itcut=1, isubs=0, ibrcl=0, iocean=0
method=3, idebug=0, iz_space=1, ivis0=1, ifix_dt=0, new_vis=-1
isfc = 0, 0  ! isfc = 0 sets surface flux (qstar), isfc = 1 sets surface condition (Dirichlet through tsfcc)
iDNS = 0
ifields = 0  ! if set to 1, reads the initial condition from path_ran below; otherwise calls subroutine random
isurfUlog=0 !0: constant lower boundary, 1: log profile for first 100 m
itempfld=1 ! 0: horizontally homogeneous, 1: Gaussian bubble
ilin = 0   !1 means linear interpolation for particle, 6th order otherwise
icouple=0
iTcouple=1
iHcouple=1  !iHcouple also controls TE couple
ievap=1
ineighbor=0
icoalesce=0
ipart_method=2   !1 = same RK3 as flow, 2 = backward Euler
ipartdiff = 1
isfs = 1   !isfs=1 is simple random walk, isfs=2 is Weil et al. (2004) particle subgrid
iexner = 1  !iexner = 1 applies a correction to part%Tf to account for conversion b/w potential temp and temp. (need in PBL)
iparsdefault=0 
ifixnumpart=1 
iSSGF=0 
ivpran=0 
irpsmp=0 
iFnumsmp=0 
iszcls=0 
impscl=0 
impmnl=0 
/

!Grid and domain parameters
&grid_params
ncpu_s=8
ichannel=0

!Time step
cfl = 0.5
dt_new = 1.0

!Use for DNS:
Uo = 0.0  !Sets the magnitude of the top and bottom plate velocities (Couette flow -- assumes equal and opposite)
Ttop(1) = 300.15, 95.0  !Currently this should be (temperature, relative humidity) (used in DNS)
Tbot(1) = 300.15, 100.0  !Currently this should be (temperature, relative humidity) (used in DNS)

!Use for LES:
qstar = 0.0, 0.0  !Surface fluxes of (temperature, humidity) (used for LES and DNS)
tsfcc = 300.0, 098.0    !Surface conditions of (temp, humidity) (used for LES) -- make sure tsfcc is gt than t00 for both isfc=0 or 1


psurf = 101325  !Pa -- nominal surface pressure, for use in the definition of potential temp

ugcont = 10.0   !The initial u-velocity in the field
vgcont = 0.0   !The initial v-velocity

dpdx = -0.69325  !The pressure gradient for channel flow (gets put into u_geo for DNS)

zi = 600.0  !This has to do with grid stretching; make equal to zl for Couette/channel
zl = 1000.0
xl = 1000.0
yl = 1000.0

zw1 = 1.25  !The first grid point
/

!Set the paths for code I/O. Must be on the scratch directory, not AFS!
&path_names
path_seed="/scratch365/tpeng2/NTLP-SD/ABLevap/"
path_part="/scratch365/tpeng2/NTLP-SD/ABLevap/part.le.cou000"
path_res="/scratch365/tpeng2/NTLP-SD/ABLevap/u.le.cou000"
path_sav="/scratch365/tpeng2/NTLP-SD/ABLevap/"
path_his="/scratch365/tpeng2/NTLP-SD/ABLevap/"
path_viz_xy="/scratch365/tpeng2/NTLP-SD/ABLevap/"
path_viz_xz="/scratch365/tpeng2/NTLP-SD/ABLevap/"
path_viz_yz="/scratch365/tpeng2/NTLP-SD/ABLevap/"
path_stuf="/scratch365/tpeng2/NTLP-SD/ABLevap/"
path_ran="/scratch365/tpeng2/NTLP-SD/ABLevap/u.le.cou000"
path_histog="/scratch365/tpeng2/NTLP-SD/ABLevap/"
/


!Material and particle properties and other constants
&constants
grav = 9.81
t00 = 273.0 !Reference temp for buoyancy
fcor = 1.0e-4  !Coriolis parameter
zo = 0.05       !Roughness height (for LES)

!Air phase:
rhoa=1.0
nuf=1.57e-5  !m^2/s
Cpa=1006.0  !J/kg-K
Pra = 0.715
Sc = 0.615
Rd = 286.9     !J/kg-K  dry air gas constant

!Particles:
tnumpart = 6e3
mult_init = 1e3
massfrac=0.05 !Before being multiped by mult_init
!tnumpart = 131072
!mult_init = 64e6
rhow=1000.0  !kg/m^3
part_grav = 0,0,-9.81
Vpmax = 0.9
Cpp = 4179.0  !J/kg-K
Mw = 0.018015  !kg/mol
Ru = 8.3144    !J/mol-K Universal gas constant
Ms = 0.05844  !kg/mol: molecular weight of salt
Sal = 34.0   !INITIAL salinity
Gam = 7.28e-2
Ion = 2.0
Os = 1.093

!Particle initial conditions:
radius_init=25e-6
Tp_init = 300.0
vp_init = 0.0, 0.0, 1.0
qf_init = 0.01

/
