# pumpkyn / pumpkynPie signature catalog

**GENERATED — do not edit.** Regenerate with
`python3 gen_pumpkyn_catalog.py > pumpkyn_catalog.md`.
Judgement and campaign relevance live in `pumpkyn_reference.md`.

- `pumpkyn` @ **5f5ca31**  (`/Users/msc/Desktop/proj7/external/pumpkyn`)
- `pumpkynPie` @ **47be599**  (`/Users/msc/Desktop/proj7/external/pumpkynPie`)

## pumpkyn

### `+pumpkyn/+cr3bp` — 42 routines

| routine | signature | purpose |
|---|---|---|
| `cont` | `[x0,tau0,converged,nullDF,iter,convErr] = cont(x0,tau0,ds,nullDF,mu,tol)` | Given an initial state vector in the CR3BP system, this routine |
| `cont_np` | `[x0,iter] = cont_np(x0,tau0,mu,tol)` | Fixed Period Natural Parameter Continuation.  Given an initial guess |
| `directLambert` | `[v0,vf,converged] = directLambert(r0,rf,tof,pm,muStar,tStar,lStar,M,P,dim3)` | This routine will compute a single-revolution direct lambert-style |
| `dop` | `dop = dop(rRcv,rSat,maskIdx,dim3)` | this routine will compute various dillusion of precision metrics for |
| `eci2orb` | `oev = eci2orb(mu,r,v,dim3)` | convert eci state vector to six classical orbital |
| `elAng` | `phi = elAng(rObs,rTgt,dr,dim3)` | This routine will determine the elevation angle of rTgt relative to rObs. |
| `eom` | `etaDot = eom(tau,eta,mu,dim3)` | This code computes the equations of motion associated with the |
| `fromJ2K` | `[x0,tStar,lStar] = fromJ2K(jd0,rv,muStar,M,P,dim3)` | This routine will take J2000 Inertial states with respect to |
| `fromLLA` | `x = fromLLA(jd,lla,muStar,lStar,P,dim3)` | This routine will take observer latitude, longitude, and altitude |
| `fromOrb` | `rv = fromOrb(t,oev,M,muStar,tStar,lStar,P)` | This routine will take keplerian orbital elements with respect to |
| `fromPCI` | `rvND = fromPCI(t,rv,muStar,tStar,lStar,P)` | This routine will take Planet Ceneterd Inertial states with respect to |
| `getTulip` | `[tau0,x0,mu,tStar,lStar] = getTulip(tau0,Np,pm,tol)` | This routine will fetch previous cataloged initial states for all |
| `jacobi` | `[C,U] = jacobi(x0,muStar)` | This routine will compute the jacobi constant for a given |
| `lagrangePts` | `Lpts = lagrangePts(mu)` | This routine will compute the five lagrange points for any Circular |
| `manifolds` | `[x0_s,x0_u] = manifolds(tau0,x0,muStar,epsilon)` | Given an inital state of an unstable peridic three-body orbit, this |
| `maxGaps` | `maxGapTimes = maxGaps(tau,r,rPts,rPt0,minElAng,Ns,dTau)` | Computes the maximum surface gap times for a set of discretized points |
| `minDeltaV` | `dV = minDeltaV(x0,x1,muStar)` | This routine will compute the theoretical minimum cost associated with |
| `minLambert` | `[v0,vf,converged] = minLambert(rv0,rvf,tof,muStar,tStar,lStar,M,dim3)` | This routine will compute a single-revolution direct lambert-style |
| `monodromy` | `M = monodromy(x0,tau0,muStar)` | This routine will compute the monodromy matrix assoicated with |
| `occultationCalc` | `[tauDur,pIdx,oIdx] = occultationCalc(tau,x,rP1,rP2,muStar)` | This routine will determine the occultations of either primary body |
| `orb2eci` | `[r, v] = orb2eci(mu, oev, dim3)` | convert classical orbital elements to eci state vector |
| `orbitProperties` | `data = orbitProperties(x0,tau0,muStar,lStar)` | This routine efficiently computes and aggregates key diagnostic |
| `pointSphere` | `[r,lla,xyz] = pointSphere(N,dr,rP)` | The purpose of this routine is to create a sphere of points that could |
| `primary2PosVel` | `[r,v] = primary2PosVel(tau,muStar)` | This routine will compute the dimensionless position and velocity of |
| `prop` | `[tau,x] = prop(tau,x0,mu)` | This routine will take the initial states of a CR3BP and propagate them |
| `showEarth` | `[h,globe] = showEarth(jd0,lStar,muStar,hIn,varargin)` | This routine will properly place the Earth in dimensionless coordinates |
| `showEarthMoonSystem` | `[h,earthGlobe,moonGlobe,sunGlobe] = ...` |  |
| `showMoon` | `[h,globe,hIn] = showMoon(lStar,muStar,hIn,varargin)` | This routine will properly place the moon in dimensionless coordinates |
| `showSun` | `[h,globe] = showSun(jd0,lStar,muStar,M,hIn)` | Place or update the Sun in the Earth-Moon rotating barycentric frame. |
| `stabilityIndex` | `k = stabilityIndex(x0,tau0,muStar)` | This routine will compute the stability index of periodic orbit. |
| `stationKeeping_deltaV` | `[dV_LB,dV_UB] = stationKeeping_deltaV(x0,tau0,dVMaxErr,dTau,muStar,N)` | This routine will compute the upper bound associated with the station |
| `tfMin` | `sol_lambda0_tf = tfMin(rv0,rvf,lambda0_tf,Tmax,c,muStar,solverType)` | Solves a minimum-time transfer between two rotating-frame states in the |
| `tfMinEoM` | `[yDot,Ht,dHdy,aThrust] = tfMinEoM(tau,y,Tmax,c,muStar) %#ok<INUSD>` | This routine contains the low-thrust equations of motion in the CR3BP |
| `tfMinProp` | `[tau,y,Hf,dHdy] = tfMinProp(tf,y0,Tmax,c,muStar)` | This routine will use event detection to stop and restart integration |
| `threeBody` | `[rvDot,aP1,aP2] = threeBody(t,rv,muStar,M,tStar,lStar,P,dim3)` | This routine will propagate a satellite using the cartesian representation |
| `toJ2K` | `[rv,tStar,lStar] = toJ2K(jd0,x0,muStar,M,P,dim3)` | This routine will take CR3BP dimensionless states and convert them |
| `toLLA` | `lla = toLLA(jd,pos,muStar,lStar,P,dim3)` | This routine will take dimensionless states in the rotating barycentric |
| `toModEq` | `[meOEV,theta,psi] = toModEq(tau,rvND,aT,M,muStar,tStar,lStar,P)` | This routine will take the dimensionless states in the rotating |
| `toOrb` | `oev = toOrb(tau,rvND,M,muStar,tStar,lStar,P)` | This routine will take the dimensionless states in the rotating |
| `toPCI` | `rv = toPCI(tau,rvND,muStar,tStar,lStar,P)` | This routine will take CR3BP dimensionless states and convert them |
| `tulipConstellation` | `[tau,r,v,mu,lStar,tStar] = tulipConstellation(Np,tau0,Nr,pm,dtau)` | This routine will generate a constellation of satellites spaced equally |
| `twoBody` | `rvDot = twoBody(t,rv,mu,dim3)` | This routine will propagate a satellite using the two-body equations |

### `+pumpkyn/+pykep` — 1 routines

| routine | signature | purpose |
|---|---|---|
| `lambert2Body` | `[v0,vf] = lambert2Body(r0,rf,dt,mu,dir,nmax,dim3)` | This is a wrapper interface for the pykep implementation of the |

### `+pumpkyn/+util` — 35 routines

| routine | signature | purpose |
|---|---|---|
| `ECEF2LLA` | `[lla] = ECEF2LLA(posvececef,dim3)` | function [lla] = ECEF2LLA(posvececef, dim3) |
| `ECI2LLA` | `LLA = ECI2LLA(JD,r_ECI,dim)` | Convert ECI (CIS, Epoch J2000.0) Coordinates to Latitude/Longitude and |
| `JD2GAST` | `GAST = JD2GAST(JD)` | THETAm is the mean siderial time in degrees |
| `JD2GMST` | `GMST = JD2GMST(JD)` | Find the Julian Date of the previous midnight, JD0 |
| `LLA2ECEF` | `[xyz] = LLA2ECEF(lla,dim3)` | LLA2ECEF - convert latitude, longitude, and altitude to |
| `LLA2ECI` | `[r,v] = LLA2ECI(jd,lla,dim3)` | LLA2ECI - convert latitude, longitude, and altitude to |
| `Rx` | `R = Rx(theta)` | Generate a 3-D rotation matrix to rotate another DCM or vector about the |
| `Rz` | `R = Rz(theta)` | Generate a 3-D rotation matrix to rotate another DCM or vector about the |
| `addFigureLogo` | `[hLogoAxes,hLogo] = addFigureLogo( ...` |  |
| `arclength` | `[s,cs] = arclength(r,dim3)` | This routine computes the total length of a line given in vector form |
| `bsxAng` | `ang = bsxAng(A,B,dim3)` | Take two multi-dimensional arrays which contain vector data and compute |
| `bsxDot` | `C = bsxDot(A,B,dim3)` | Take two multi-dimensional arrays which contain vector data and compute |
| `chebyspace` | `y = chebyspace(a,b,n)` | chebyspace:  create a chebyshev node spaced vector of n points between |
| `eDim` | `eND = eDim(fND,fSeq)` | This routine will reconstruct an N-dimensional matrix from a flattened |
| `earth3D` | `[www,globe,hLight,cdata] = earth3D(options,www)` | Provides an interface to view the 3D earth in space.  This allows for |
| `fDim` | `[fND,fSeq] = fDim(ND,dim)` | This routine will flatten an N-dimensional matrix into a 2-D matrix |
| `getConst` | `tmpVar = getConst(varargin)` | Purpose; |
| `grs2rgb` | `res = grs2rgb(img, map)` | Convert grayscale images to RGB using specified colormap. |
| `juliandate` | `jd = juliandate( varargin )` | This routine is a replacement for the MATLAB juliandate routine which |
| `mean2True` | `TA = mean2True(MA,e)` | Convert mean anomaly to true anomaly |
| `mntimes` | `z = mntimes(x,y,xrowdim,xcoldim,yrowdim,ycoldim)` | matrix n-d times routine (mntimes.m) takes any n-d matricies and |
| `moon3D` | `[h,globe,hLight] = moon3D( ...` |  |
| `moonPosVel` | `[r,v] = moonPosVel(jd)` | High precision Earth Moon position and velocity routine adapted |
| `multiplyDCM` | `s = multiplyDCM(DCM,r,dim3)` | Take a L x M x N DCM and multiply it by a N x M vector |
| `planetPosVel` | `[r,v] = planetPosVel(jd0,center,target,varargin)` | This routine computes the position and velocity of one celestial body |
| `plotUnc` | `varargout = plotUnc(varargin)` | This routine will plot the 2D uncertainty bounds about a standard plot |
| `sDim` | `[sND1,sND2] = sDim(ND,dim)` | Given a single input matrix with it's singleton specifier defined by |
| `setOutwardNormals` | `normals = setOutwardNormals(globe,center)` | setOutwardNormals |
| `sphere3D` | `[x,y,z,xS,yS,zS] = sphere3D(N)` | This routine will compute a sphere using parametric equations which |
| `stars3D` | `hStarSphere = stars3D( ...` |  |
| `sun3D` | `[h,globe] = sun3D(posOffset,overlay,scale,h)` | Provides an interface to view the 3D Sun in space.  This allows for |
| `sunPosVel` | `[r,v,a] = sunPosVel(jd)` | Medium precision sun position and velocity routine adapted |
| `vmag` | `out = vmag(in,dim)` | Quickly compute the vector magnitude of any given N-D Vector. This is |

## pumpkynPie

### `+pumpkynPie/+bcr4b` — 3 routines

| routine | signature | purpose |
|---|---|---|
| `eom` | `etaDot = eom(tau,eta,opts)` | Bi-Circular Restricted Four-Body Equations of Motion and corresponding |
| `prop` | `[tau,x] = prop(tau,x0,opts)` | This routine will take the initial states of a BCR4B and propagate them |
| `showSun` | `[h,globe] = showSun(psi0,lStar,hIn)` | This routine will properly place the Sun in dimensionless coordinates |

### `+pumpkynPie/+cr3bp` — 28 routines

| routine | signature | purpose |
|---|---|---|
| `brouckeParams` | `[alpha,beta] = brouckeParams(x0,tau0,mu,dim3,sundman_flag)` | This routine will compute the stability parameters of an initial state |
| `cont` | `[x0,tau0,converged,nullDF,iter,convErr,ds] = cont(x0,tau0,ds,nullDF,mu,tol,sundman_flag)` | Given an initial state vector in the CR3BP system, this routine |
| `cont_np` | `[x0,tau0,iter] = cont_np(x0,t0,mu,tol,sundman_flag)` | Fixed Period Natural Parameter Continuation.  Given an initial guess |
| `eom` | `etaPrime = eom(~,eta,dim3,opts)` | This code computes the equations of motion associated with the |
| `geometricAccess` | `[losok, out] = geometricAccess(rObs,rTgt,params)` | This routine will compute geometric access between observer and target |
| `getAxial` | `[tau0,x0,mu,tStar,lStar] = getAxial(tau0,Lpt,pm)` | This routine will get Axial orbits about the Earth-Moon system |
| `getAxialData` | `data = getAxialData(tau,Lpt,pm)` | This routine will get a Axial orbit about the Earth-Moon system |
| `getCR3BParams` | `params = getCR3BParams()` | Earth-Moon Parameters: |
| `getCycler` | `[tau0,x0,mu,tStar,lStar] = getCycler(tau0,resStr,varargin)` | This routine will get cycler orbits about the |
| `getCyclerData` | `data = getCyclerData(tau,resRatio)` | This routine will get a cycler orbit about the Earth-Moon system |
| `getDPO` | `[tau0,x0,mu,tStar,lStar] = getDPO(tau0,varargin)` | This routine will get Distant Prograde orbits about the |
| `getDPOData` | `data = getDPOData(tau)` | This routine will get a Halo orbit about the Earth-Moon system |
| `getDRO` | `[tau0,x0,mu,tStar,lStar] = getDRO(tau0,varargin)` | This routine will get DRO orbits about the Earth-Moon system |
| `getDROData` | `data = getDROData(tau)` | This routine will get a Halo orbit about the Earth-Moon system |
| `getHalo` | `[tau0,x0,mu,tStar,lStar] = getHalo(tau0,Lpt,pm)` | This routine will get a Halo orbit about the Earth-Moon system |
| `getHaloData` | `data = getHaloData(tau,Lpt,pm)` | This routine will get a Halo orbit about the Earth-Moon system |
| `getLPO` | `[tau0,x0,mu,tStar,lStar] = getLPO(tau0,varargin)` | This routine will get Low Prograde orbits about the |
| `getLPOData` | `data = getLPOData(tau0,varargin)` | This routine will get the LPO data about the Earth-Moon system |
| `getLyapunov` | `[tau0,x0,mu,tStar,lStar] = getLyapunov(tau0,Lpt,pm)` | This routine will get a lyapunov orbit about the Earth-Moon system |
| `getLyapunovData` | `data = getLyapunovData(tau,Lpt,pm)` | This routine will get a lyapunov orbit about the Earth-Moon system |
| `getPumpkin` | `[tau0,x0,mu,tStar,lStar] = getPumpkin(Np,pm,varargin)` | This routine will get pumpkin orbits about the |
| `getTulipData` | `data = getTulipData(tau0,Np)` | This routine will get a cycler orbit about the Earth-Moon system |
| `orbitProperties` | `data = orbitProperties(x0,tau0,muStar,lStar,sundman_flag)` | This routine efficiently computes and aggregates key diagnostic |
| `phaseConstellation` | `x = phaseConstellation(tau,x0,tau0,N,muStar,phaseType)` | This routine will generate an N-satellite phased constellation from a |
| `plotCR3BParams` | `plotCR3BParams(data)` | This routine will view select plots provided by the data input structure |
| `prop` | `[tau,x,t] = prop(tau,x0,mu,sundman_flag)` | This routine will take the initial states of a CR3BP and propagate them |
| `propGPS` | `x = propGPS(tau,jd0,params,alm)` | This routine propagates the GPS constellation from a parsed YUMA |
| `visibility` | `Psi = visibility(r1,r2,mu,Re,Ne,dim3)` | Determine whether two or more satellites are visible to eachother by |

### `+pumpkynPie/+dyn` — 4 routines

| routine | signature | purpose |
|---|---|---|
| `multipleShooting` | `[t0c,X0c,exitflag] = multipleShooting(dynFun,t0,X0,odeOpts,sOpts,varargin)` | Multiple shooting method for transitioning ephemeris from the |
| `palc` | `[x0,tau0,converged,nullDF] = palc(dynFun,x0,tau0,ds,nullDF,tol,varargin)` | This routine will find a periodic orbit (and it's corresponding period |
| `propBetweenPatchPoints` | `[t,r,v] = propBetweenPatchPoints(dynFun,tPts,rPts,vPts,varargin)` | This routine will integrate (propagate) the equations of motion |
| `psdc` | `outputData = psdc(dynFun,tauDur,tauRef,xRef,opts,varargin)` | This uses the Principal Stretching Direction Control (PSDC) method, |

### `+pumpkynPie/+nav/+gnss` — 2 routines

| routine | signature | purpose |
|---|---|---|
| `readYuma` | `data = readYuma(source)` | This routine reads a GPS YUMA almanac from either a local file path or a |
| `yuma2ECEF` | `[rECEF,vECEF] = yuma2ECEF(jd,alm,wrapTK)` | This routine propagates GPS YUMA almanac records into WGS-84 |

### `+pumpkynPie/+plot` — 3 routines

| routine | signature | purpose |
|---|---|---|
| `outage` | `[h,s,fig,outputData] = outage(t,aOK,names)` | This routine will plot the outages associated with user provided |
| `polarcontourf` | `polarcontourf(lon_deg, lat_deg, Z, fill_levels, line_levels, varargin)` | POLARCONTOURF Creates a filled polar contour centered on the South Pole. |

### `+pumpkynPie/+util` — 10 routines

| routine | signature | purpose |
|---|---|---|
| `accessIntervals` | `[trueDur,tStart,tEnd] = accessIntervals(t,y)` | Compute the duration as well as the corresponding start time, and end time |
| `arclinspace` | `[tauPts,rvPts] = arclinspace(tau,rv,Npts,Nrev)` | equalArcTau |
| `bsxAng` | `ang = bsxAng(A,B,dim3)` | Take two multi-dimensional arrays which contain vector data and compute |
| `bsxDot` | `C = bsxDot(A,B,dim3)` | Take two multi-dimensional arrays which contain vector data and compute |
| `chebyspace` | `y = chebyspace(a,b,n)` | chebyspace:  create a chebyshev node spaced vector of n points between |
| `diagND` | `dND = diagND(ND,rowDim,colDim)` | Given an ND matrix containing two known dimensions (e.g. row and column |
| `maxMinAccessBlockBounds` | `[maxMinAccessDur,bestOffset,passFail] = maxMinAccessBlockBounds(t,idx,req)` | Find the bin-boundary offset that maximizes the minimum qualified access |
| `maxMinAccessBounds` | `[maxMinAccessFrac,bestOffset] = maxMinAccessBounds(t,idx,req)` | maxMinAccessBounds: find the bin-boundary offset that maximizes the |
| `plot3c` | `h = plot3c(varargin)` | This routine will plot a three-dimensional line with interpolated color |
| `updateSatTrailEffects` | `hTrail = updateSatTrailEffects( ...` |  |

