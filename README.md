# NTLP (NCAR Turbulence with Lagrangian Particles)

## Compilation
To build on the CRC machines you must run the following commands:
```
module load cmake
module load mvapich2
module load intel

mkdir build
cd build 
cmake ..
make
```
## DEBUG (using -check bounds)
```
mkdir debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```

## SETUP AND RUNNING
To run, make a directory ("case1" or something) where les.run and params.in will go
(i.e., not out of the same directory as les.F)
Make sure all paths in these directories point to the proper locations

## Parameter files: params.f90
### Solver features
| Flag  | meaning | values |
| ------------- | ------------- | ------------- |
| ismlt  | Content Cell  | Content Cell  |
| ifree  | MO condition  | 0 (spatially averaged surface condition for MO, call lower); 1 (point-by-point conditions for MO free convection (call lower_free))  |
| iradup  | Content Cell  | Content Cell  |
| iupwnd  | Content Cell  | Content Cell  |
| ibuoy  | Content Cell  | Content Cell  |
| ifilt  | Content Cell  | Content Cell  |
| itcut  | Content Cell  | Content Cell  |
| isubs  | Content Cell  | Content Cell  |
| ibrcl  | Content Cell  | Content Cell  |
| icoean  | Content Cell  | Content Cell  |
| method  | Content Cell  | Content Cell  |
| idebug  | Content Cell  | Content Cell  |
| iz_space  | Content Cell  | Content Cell  |
| ivis0  | eddy viscosity model   | 0 (old); 1(new); will be passed to ivis  |
| ifix_dt  | fixed time step  | 0 (controled by CFL number); 1(fixed); set in sr, get_dt  |
| new_vis  | The iteration step for which the new model is turned on when ivis0=1;  | -1(always on);step;  |
| ```iDNS```  | DNS solver  | ```0```(DNS); ```1```(LES) |
| ifields  | Content Cell  | Content Cell  |
| ilin  | Content Cell  | Content Cell  |
| icouple  | Content Cell  | Content Cell  |
| iTcouple  | Content Cell  | Content Cell  |
| iHcouple  | Content Cell  | Content Cell  |
| ievap  | Content Cell  | Content Cell  |
| icoalesce  | Content Cell  | Content Cell  |
| ipart_method  | particle ODE intergration methods  | 1(RK3); 2(BE) |
| ipartdiff  | Content Cell  | Content Cell  |
| isfs  | Particle subgrid model  | 1(simple random walk);2(Weil et all (2004))  |
| iexner  | Content Cell  | Content Cell  |

### Time Stepping

| Flag  | meaning | values |
| ------------- | ------------- | ------------- |
| cfl  | CFL number  | 0.4 (default)  |
| dt_new  | default time steps; system takes the minimum between dt_new and CFL condition | 0.5 (default) |


### Initial and Boundary Conditions

| Flag  | meaning | values |
| ------------- | ------------- | ------------- |
| ```qstar```  | surface fluxes of (temperature,humidity) | ```(0.2,0.2)``` (default)  |
| ```dt_new```  | default time steps; system takes the minimum between dt_new and CFL condition | ```0.5``` (default) |
| ```ugcont```  | The initial u-velocity in the field | ```10``` (default) |
| ```dpdx```  | Geostrophic wind speed | ```10``` (default) |
| ```fcor```  | Coriolis force | ```1e-04``` (default) |
| ```zo```  | Roughness height (for LES) | ```0.05``` (default) |
