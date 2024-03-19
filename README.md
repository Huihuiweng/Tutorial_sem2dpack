
This code is used to simulate 2D and 2.5D dynamic rupture simulations, with slip-weakening, time-weakening, rate-and-state, and thermal pressurization friction laws. This tutorial aims to guide you to simulate the models and post process the results.

The structure of the tutorial document:
code              -- contains the original sem2dpack code
simulations       -- contains the model that can be run directly'
scripts_templates -- contains matlab and python scripts
post-process      -- choose the script you like and save it here


Current directory is defined as ${workdir}

Installation:
1. cd ${workdir}/code/sem2dpack/SRC
2. vi Makefile   ### Modify the executable dir and the type of compiler in Makefile
3. make

Running model:
1. cd ${workdir}/simulations
2. vi Par.inp    ### Modify the model parameters, refer to the manual of sem2dpack at the end of this file
3. ${path_to_bin}/sem2dsolve Par.inp

Post process results:
1. cd ${workdir}/post-process

! python version
2. cp ../scripts_templates/python/* .
3. python Example_figure.py

! matlab version
2. cp ../scripts_templates/matlab/* .  
3. open Example_figure.m in MATLAB and run the script directly.


# Manual for running sem2dpack (from sem2dpack)
 ----------------------------------------------------------------------------

 NAME   : BC_ABSORB
 GROUP  : BOUNDARY_CONDITION
 PURPOSE: Absorbing boundary
 SYNTAX : &BC_ABSORB stacey, let_wave /

  stacey   [log] [F] Apply Stacey absorbing conditions for P-SV.
                Higher order than Clayton-Engquist (the default).
  let_wave [log] [T] Allow incident waves across this boundary 
                if mechanism='WAVE' in &SRC_DEF

 NOTE   : Stacey conditions and incident waves are only implemented for 
          vertical and horizontal boundaries

 ----------------------------------------------------------------------------

 NAME   : BC_DIRNEU
 GROUP  : BOUNDARY_CONDITION
 PURPOSE: Dirichlet (null displacement) 
          and/or Neumann (null or time-dependent traction) 
          boundary conditions on vertical or horizontal boundaries
 SYNTAX : &BC_DIRNEU h, v, hsrc, vsrc /
          possibly followed by one or two STF_XXXX blocks

  h        [char]['N'] Boundary condition on the horizontal component
  v        [char]['N'] Boundary condition on the vertical component :
                       'N' : Neumann 
                       'D' : Dirichlet
  hsrc     [name]['none'] Name of the source time function for a
                time-dependent horizontal traction: 
                'RICKER', 'TAB', 'USER', etc  (see STF_XXXX input blocks)
  vsrc     [name]['none'] Same for the vertical component

 ----------------------------------------------------------------------------

 NAME   : BC_DYNFLT
 GROUP  : BOUNDARY_CONDITION, DYNAMIC_FAULT
 PURPOSE: Dynamic fault with friction
 SYNTAX : &BC_DYNFLT friction, cohesion|cohesionH, opening, Tn|TnH, Tt|TtH,
                     Sxx|SxxH, Sxy|SxyH, Sxz|SxzH, Syz|SyzH, Szz|SzzH
                     ot1, otd, oxi, osides /
          followed, in order, by:
          1. &DIST_XXX blocks (from the DISTRIBUTIONS group) for arguments
             with suffix H, if present, in the order listed above.
          2. &BC_DYNFLT_SWF, &BC_DYNFLT_TWF or &BC_DYNFLT_RSF block(s) 
             (if absent, default values are used)
          3. &BC_DYNFLT_NOR block (if absent, default values are used)

  friction [name(2)] ['SWF',''] Friction law type:
                  SWF = slip weakening friction
                  TWF = time weakening friction
                  RSF = rate and state dependent friction
                Some friction types can be combined. E.g. to set the 
                friction coefficient to the minimum of SWF and TWF, set 
                  friction='SWF','TWF'
  cohesion [dble] [0d0] part of the strength not proportional to normal stress
  opening  [log] [T] Allow fault opening instead of tensile normal stress
  Tn       [dble] [0d0] Initial normal traction (positive = tensile)
  Tt       [dble] [0d0] Initial tangent traction 
                 (positive antiplane: y>0; positive inplane: right-lateral slip)
  Sxx      [dble] [0d0] Initial stress sigma_xx
  Sxy      [dble] [0d0] Initial stress sigma_xy
  Sxz      [dble] [0d0] Initial stress sigma_xz
  Syz      [dble] [0d0] Initial stress sigma_yz
  Szz      [dble] [0d0] Initial stress sigma_zz
  otd      [dble] [0.d0] Time lag between outputs (in seconds)
                Internally adjusted to the nearest multiple of the timestep
                Its value can be found in the output file FltXX_sem2d.hdr
                The default internally resets otd = timestep
  ot1      [dble] [0.d0] Time of first output (in seconds)
                Internally adjusted to the nearest multiple of the timestep
                Its value can be found in the output file FltXX_sem2d.hdr
  oxi      [int(3)] [(1,huge,1)] First, last node and stride for output
                The default resets oxi(2) = last fault node
  osides   [log] [F] Export displacement and velocities on each side
                of the fault
  V        [dble] [1d-12] Initial velocity (needed for RSF)

 NOTE: The initial stress can be set as a stress tensor (Sxx,etc), as
       initial tractions on the fault plane (Tn and Tt) or as the sum of both.

 NOTE: We recommend to use dynamic faults with the leapfrog time scheme
       and a layer of Kelvin-Voigt damping material near the fault.

 ----------------------------------------------------------------------------

 NAME   : BC_DYNFLT_NOR
 GROUP  : DYNAMIC_FAULT
 PURPOSE: Normal stress response for dynamic faults.
 SYNTAX : &BC_DYNFLT_NOR kind, V, L, T /

  kind     [int] [1] Type of normal stress response:
                       0 = shear strength is independent of normal stress
                           (the cohesive strength is set as the product of
                           friction coefficient and initial normal stress)
                       1 = Coulomb 
                       2 = Prakash-Clifton with regularizing time scale
                       3 = Prakash-Clifton with regularizing length scale
  T        [dble] [1d0] Regularization time scale if kind=2
  V        [dble] [1d0] Characteristic velocity if kind=3
  L        [dble] [1d0] Regularization length scale if kind=3

 ----------------------------------------------------------------------------

 NAME   : BC_DYNFLT_RSF
 GROUP  : DYNAMIC_FAULT
 PURPOSE: Velocity and state dependent friction
 SYNTAX : &BC_DYNFLT_RSF kind, Dc | DcH, Mus | MusH , 
                         a | aH, b | bH, Vstar | VstarH /
          followed by &DIST_XXX blocks (from the DISTRIBUTIONS group) for
          arguments with suffix H, if present, in the order listed above.

  kind     [int] [1] Type of rate-and-state friction law:
                       1 = strong velocity-weakening at high speed
                           as in Ampuero and Ben-Zion (2008)
  Dc       [dble] [0.5d0] Critical slip 
  MuS      [dble] [0.6d0] Static friction coefficient
  a        [dble] [0.01d0] Direct effect coefficient
  b        [dble] [0.02d0] Evolution effect coefficient
  Vstar    [dble] [1d0] Characteristic or reference slip velocity
  theta    [dble] [1d0] State variable

 ----------------------------------------------------------------------------

 NAME   : BC_DYNFLT_SWF
 GROUP  : DYNAMIC_FAULT
 PURPOSE: Slip-weakening friction
 SYNTAX : &BC_DYNFLT_SWF Dc | DcH, MuS | MuSH , MuD | MuDH, healing /
          followed by &DIST_XXX blocks (from the DISTRIBUTIONS group) for
          arguments with suffix H, if present, in the order listed above.

  kind     [int] [1] Type of slip weakening function:
                       1 = linear
                       2 = exponential
  Dc       [dble] [0.5d0] Critical slip 
  MuS      [dble] [0.6d0] Static friction coefficient
  MuD      [dble] [0.5d0] Dynamic friction coefficient
  healing  [log] [F] Instantaneous healing upon slip arrest
               Healing is currently valid only with the leapfrog time scheme

 ----------------------------------------------------------------------------

 NAME   : BC_DYNFLT_TWF
 GROUP  : DYNAMIC_FAULT
 PURPOSE: Time weakening friction for dynamic faults 
          with prescribed rupture speed.
 SYNTAX : &BC_DYNFLT_TWF kind, MuS, MuD, Mu0, X, Z, V, L, T /

  kind     [int] [1] Type of time-weakening history:
               1 = expansion at constant speed V up to time T
               2 = expansion at decreasing speed then contraction
                   as in Andrews and Ben-Zion (JGR 1997, eqs 2 and 3)
  MuS      [dble] [0.6d0] Static friction coefficient
  MuD      [dble] [0.5d0] Dynamic friction coefficient
  Mu0      [dble] [0.6d0] Friction coefficient at the hypocenter at time=0
  X,Z      [dble] [0d0] Position of hypocenter
  V        [dble] [1d3] Rupture propagation speed (initial speed if kind=2)
  L        [dble] [1d0] Size of weakening zone
  T        [dble] [huge] Total duration

 NOTE   : Time-weakening is usually applied as an artificial nucleation procedure.
          The maximum size of the nucleation region is 2*V*T if kind=1, V*T/2 if kind=2

 ----------------------------------------------------------------------------

 NAME   : BC_DEF
 PURPOSE: Define a boundary condition
 SYNTAX : &BC_DEF tag, tags, kind /
          possibly followed by &BC_kind blocks 

   tag     [int] [none] A number assigned to the boundary. If you are
               using SEM2D built-in structured mesher the conventions are:
                       1       bottom
                       2       right
                       3       up
                       4       left
               If you are importing a mesh, you must use the tags assigned
               to the boundaries during the mesh construction.
   tags    [int(2)] [none] Two tags are needed for split-node interfaces (faults)
               and for periodic boundaries.
   kind    [char*6] [none] Type of boundary condition. The following are
               implemented:
               'DIRNEU', 'ABSORB', 'PERIOD', 'LISFLT', 'DYNFLT', 'KINFLT'

 NOTE   : Most of the boundary conditions need additional data, given
          in a BC_kind input block of the BOUNDARY_CONDITIONS group
          immediately following the BC_DEF block.

 ----------------------------------------------------------------------------

 NAME   : BC_KINFLT
 GROUP  : BOUNDARY_CONDITION
 PURPOSE: Prescribe a kinematic source (the spatio-temporal distribution of slip rate)
          on a finite fault
 STATUS : The current implementation has the following features and restrictions:
            . the fault is a flat horizontal boundary at the bottom of the model
            . prescribed (x,t)-dependent horizontal or vertical velocity
            . traction is free in the other component 
            . the source time function is the same everywhere (but variably shifted and scaled)
            . the rupture time, final slip and rise time can be spatially variable
 SYNTAX : &BC_KINFLT comp, trup|trupH, slipr|sliprH, tris|trisH, stf, ot1, otd, oxi /
          followed, in order, by:
          1. &DIST_XXX blocks (from the DISTRIBUTIONS group) for arguments
             with suffix H, if present, in the order listed above.
          2. one SOURCE TIME FUNCTION block (&STF_XXXX)
 
  comp     [int] [1] component to apply velocity, 1=horizontal, 2=vertical
  trup     [dble] [0d0] Rupture time (s)
  slipr    [dble] [0d0] Slip rate (m/s)
  tris     [dble] [0d0] Rise time (s)
  stf      [name] [none] Name of the source time function for slip rate:
                'RICKER', 'TAB', 'HARMONIC', 'BRUNE', 'USER', etc.
                See the SOURCE TIME FUNCTION blocks.
  otd      [dble] [0d0] Time lag between outputs (in seconds)
                Internally adjusted to the nearest multiple of the timestep
                Its value can be found in the output file FltXX_sem2d.hdr
                The default internally resets otd = timestep
  ot1      [dble] [0d0] Time of first output (in seconds)
                Internally adjusted to the nearest multiple of the timestep
                Its value can be found in the output file FltXX_sem2d.hdr
  oxi      [int(3)] [(1,huge,1)] First, last node and stride for output
                The default resets oxi(2) = last fault node
 

 ----------------------------------------------------------------------------

 NAME   : BC_LSF 
 GROUP  : BOUNDARY_CONDITION
 PURPOSE: Linear slip fault, a displacement discontinuity interface
          where stress and slip are linearly related
 SYNTAX : &BC_LSF Ktang | Ctang, Knorm | Cnorm /

  Ktang    [dble] [Inf] Tangential stiffness
  Ctang    [dble] [0d0] Tangential compliance
  Knorm    [dble] [Inf] Normal stiffness
  Cnorm    [dble] [0d0] Normal compliance

 NOTE: For each component:
       You can set K _or_ C, but _not_both_
       If C=0d0 or K=Inf then no discontinuity is allowed (transparent)
       If K=0d0 the fault is free stress boundary

 ----------------------------------------------------------------------------

 NAME   : DIST_GAUSSIAN
 GROUP  : DISTRIBUTIONS_2D
 PURPOSE: Bell shaped Gaussian 2D distribution 
 SYNTAX : &DIST_GAUSSIAN centered_at, length, offset, ampli, order /

  centered_at      [dble(2)] [0,0] Coordinates of the center point.
  length           [dble(2)] [1]   Characteristic lengths on each axis.
  offset           [dble] [0]      Background level.    
  ampli            [dble] [1]      Amplitude from background.
  order            [int] [1]       Exponent

 ----------------------------------------------------------------------------

 NAME   : DIST_GRADIENT
 GROUP  : DISTRIBUTIONS_2D
 PURPOSE: Constant gradient 2D distribution.
 SYNTAX : &DIST_GRADIENT file,valref ,grad,angle/

  file             [name] [none]    Name of the file containing the coordinates
                        of the points defining the reference line.
                        It is an ASCII file with 2 columns per line:
                        (1) X position (in m) and
                        (2) Z position (in m)
  valref           [dble] [none]    Value along the reference line
  grad             [dble >0] [none] Positive gradient (valref_units/meter)
  angle            [dble] [none]    Angle (degrees) between the vertical down 
                        and the grad+ direction. Anticlockwise convention (grad+
                        points down if 0, right if 90)

 NOTE   : Make sure the angle and ref-line are compatible. The code will
          abort if the ref-line is too short: some points of the domain
          cannot be projected to ref-line in the angle direction.

 ----------------------------------------------------------------------------

 NAME   : DIST_HETE1
 GROUP  : DISTRIBUTIONS_2D
 PURPOSE: Linear interpolation of values from a regular 2D grid.
 SYNTAX : &DIST_HETE1 file, col /

  file     [name] [none] Name of the file containing the definition
               of the regular grid and values at grid points.
               The format of this ASCII file is:
                 Line 1 :  ncol nx nz x0 z0 dx dz
                   ncol  = [int] number of data columns 
                   nx,nz = [2*int] number of nodes along x and z
                   x0,z0 = [2*dble] bottom-left corner 
                   dx,dz = [2*dble] spacing along x and z
                 Line 2 to nx*nz+1 : [ncol*dble] values at grid points
                   listed from left to right (x0 to x0+nx*dx), 
                   then from bottom to top (z0 to z0+nz*dx)
  col      [int] [1] Column of the file to be read

 NOTE   : The same file can contain values for (ncol) different properties,
          (e.g. rho, vp, vs) but each DIST_HETE1 block will read only one.

 NOTE   : Even if the original model domain has an irregular shape, 
          the regular grid where input values are defined must be rectangular
          and large enough to contain the whole model domain. 
          The regular grid possibly contains buffer areas with dummy values. 
          These dummy values should be assigned carefully (not random nor zero)
          because SEM2D might use them during nearest-neighbor interpolation.

 ----------------------------------------------------------------------------

 NAME   : DIST_LINEAR
 GROUP  : DISTRIBUTIONS_1D
 PURPOSE: Piecewise linear 1D distribution along X or Z.
 SYNTAX : &DIST_LINEAR n,dim,length /
            followed immediately by the interpolation data, 
            one line per point, two columns: 
            position (X or Z), value
          or
          &DIST_LINEAR file,dim,length / 
            and the interpolation data is read from a two-column file

  n        [int] [0] Number of points to be interpolated
  dim      [int] [1] Interpolate along X (dim=1) or along Z (dim=2)
  file     [name] [none] Name of the ASCII file containing the data
  length   [dble] [0] Smoothing length for sliding average window
               No smoothing if length=0

 NOTE: If the requested point is out of bounds we extrapolate linearly
       the two terminal values

 ----------------------------------------------------------------------------

 NAME   : DIST_ORDER0
 GROUP  : DISTRIBUTIONS_2D
 PURPOSE: Blockwise constant 2D distribution.
 SYNTAX : &DIST_ORDER0 xn, zn /
          x(1) ...  x(xn-1)
          z(1) ...  z(zn-1)
          v(1,1)  ... v(xn,1)          
            ...   ...   ...
          v(1,zn) ... v(xn,zn)          

  xn       [int] [none] Number of zones along X
  zn       [int] [none] Number of zones along Z
  x        [dble(xn-1)] [none] Boundaries of X-zones: 
                first zone  X < x(1), 
                second zone x(1) < X < x(2), ... 
                last zone   x(xn-1) < X
  z        [dble(zn-1)] [none] Boundaries of Z-zones
  v        [dble(xn,zn)] [none] Values inside each zone

 ----------------------------------------------------------------------------

 NAME   : DIST_PWCONR
 GROUP  : DISTRIBUTIONS_2D
 PURPOSE: Piecewise constant radial (2D) distribution.
          This distribution defines a set of annular zones, centered
          at an arbitrary reference point, and assigns constant values 
          within each zone.
 SYNTAX : &DIST_PWCONR num, ref /
             r(1)  ... ...  r(num-1)
          v(1) v(2) ... v(num-1) v(num)

  num      [int] [none] Number of annular zones (including inner and exterior)
  ref      [dble(2)] [(0d0,0d0)] Reference point: center of radial zones
  r        [dble(num-1)] [none] External radius of zones:
                first zone R <= r(1), 
                second r(1) < R <= r(2), ...
                last r(num-1) < R 
  v        [dble(num)] [none] Value inside each zone

 ----------------------------------------------------------------------------

 NAME   : DIST_SPLINE
 GROUP  : DISTRIBUTIONS_1D
 PURPOSE: Spline interpolated 1D distribution along X or Z.
 SYNTAX : &DIST_SPLINE file,dim /

  file     [name] [none] Name of the ASCII file containing
               the interpolation data, one line per point, two columns: 
               one line per point, two columns: 
               position (X or Z), value
  dim      [int] [1] Interpolate along X (dim=1) or along Z (dim=2)

 ----------------------------------------------------------------------------

 NAME   : GENERAL
 PURPOSE: General parameters
 SYNTAX : &GENERAL iexec, ngll, fmax, title, verbose, itInfo /

  iexec    [int] [0] Run level:
                       0 = just check
                       1 = solve
  ngll     [int] [9] Number of GLL nodes per edge on each spectral element
                ( polynomial order +1 ). Usually 5 to 9.
  fmax     [dble] [1.d0] The code checks if this maximum frequency is
                well resolved by the mesh and issues a warning if not. 
                This parameter is not used in computations, only for checking.
                To improve the resolution for a given fmax you must increase ngll 
                (but you will have to use shorter timesteps) or refine the mesh.
  ndof     [int] [2] Number of degrees of freedom per node
                       1 = SH waves, anti-plane
                       2 = P-SV waves, in-plane
  title    [word] [none] Title of the simulation
  verbose  [char(4)] ['1101'] Print progress information during each phase:
                       verbose(1) = input phase
                       verbose(2) = initialization phase
                       verbose(3) = check phase
                       verbose(4) = solver phase
                Example: '0001' is verbose only during solver.
  itInfo   [int] [100] Frequency (in number of timesteps) for printing
                progress information during the solver phase.
  abort_on_warnings        [log] [T] Abort if a severe warning occurs

 ----------------------------------------------------------------------------

 NAME   : MAT_DAMAGE
 GROUP  : MATERIALS
 PURPOSE: Set material properties for the damage rheology of 
          Lyakhovsky, Ben-Zion and Agnon (J. Geophys. Res. 1997) 
          and Hamiel et al (Geophys. J. Int. 2004)
 SYNTAX : &MAT_DAMAGE cp,cs,rho,phi,alpha,Cd,R,e0,ep /

  cp       [dble][0d0] P wave velocity
  cs       [dble][0d0] S wave velocity
  rho      [dble][0d0] density
  phi      [dble][0d0] internal friction angle
  alpha    [dble][0d0] initial value of damage variable
  Cd       [dble][0d0] damage evolution coefficient
  R        [dble][0d0] damage-related plasticity coefficient Cv
                 normalized by the inverse of the intact shear modulus
  e0       [dble(3)][0d0] initial total strain (11, 22 and 12)
  ep       [dble(3)][0d0] initial plastic strain (11, 22 and 12)

 ----------------------------------------------------------------------------

 NAME   : MAT_ELASTIC
 GROUP  : MATERIALS
 PURPOSE: Set material properties for a linear elastic medium
 SYNTAX : For isotropic material:
           &MAT_ELASTIC rho|rhoH, cp|cpH, cs|csH /
          For transverse anisotropy with vertical symmetry axis:
           &MAT_ELASTIC rho|rhoH, c11|c11H, c13|c13H, c33|c33H, c55|c55H, c66|c66H /
          Followed by one DIST_XXXX blocks for each argument present with suffix H,
          in the same order as listed above.

  cp       [dble][0d0] P wave velocity (m/s)
  cs       [dble][0d0] S wave velocity (m/s)
  rho      [dble][0d0] density (kg/m^3)
  c11,c13,c33,c55,c66  [dble][0d0] anisotropic elastic moduli (Pa)

 ----------------------------------------------------------------------------

 NAME   : MATERIAL
 PURPOSE: Define the material type of a tagged domain
 SYNTAX : &MATERIAL tag, kind /
          followed by one or two MAT_XXXX input blocks.

  tag      [int] [none]    Number identifying a mesh domain 
  kind     [name(2)] ['ELAST','']  Material types:
               'ELAST', 'DMG','PLAST', 'KV' 

 NOTE   : Some combinations of material kinds can be assigned to the same domain.
          Any material type can be combined with 'KV', for instance:
            &MATERIAL tag=1, kind='ELAST','KV' /
            followed by a &MAT_ELAST block and a &MAT_KV block
          sets an elastic material with Kelvin-Voigt damping. 

 ----------------------------------------------------------------------------

 NAME   : MAT_KV
 GROUP  : MATERIALS
 PURPOSE: Sets material properties for a Kelvin-Voigt viscous material.
          Adds a damping term C*v = K*eta*v, where eta is a viscous time.
          This produces attenuation with frequency-dependent quality factor
            Q(f) = 1/(eta*2*pi*f)
          Its main usage is for artificial damping of high-frequency 
          numerical artifacts generated by dynamic faults, which requires a
          thin layer of Kelvin-Voigt elements surrounding the fault
          with eta/dt = 0.1 to 0.3 and a layer thickness of 1 to 2 elements
          on each side of the fault.
 SYNTAX : &MAT_KV eta, ETAxDT /
          &MAT_KV etaH, ETAxDT / followed by a DIST_XXX input block

  eta      [dble][0d0] Viscosity coefficient
  ETAxDT   [log][T] If eta is given in units of dt (timestep)

 NOTE   : Kelvin-Voigt viscosity modifies the stability of time integration.
          The timestep (or the Courant number) must be set to a value
          smaller than usual. The critical timestep for a Kelvin-Voigt material 
          integrated with the leapfrog time scheme is
            dtc_kv = eta*( sqrt(1+dtc^2/eta^2)-1 )
          where dtc is the critical timestep for a purely elastic medium (eta=0).
          In terms of the normalized viscosity (if ETAxDT=T):
            dtc_kv = dtc / sqrt( 1+ 2*eta)

 ----------------------------------------------------------------------------

 NAME   : MAT_PLASTIC
 GROUP  : MATERIALS
 PURPOSE: Set material properties for elasto-plastic material
          with Mohr-Coulomb yield criterion
          and non-dilatant (null volumetric plastic strain)
 SYNTAX : &MAT_PLASTIC cp,cs,rho,phi,coh,Tv,e0 /

  cp       [dble][0d0] P wave velocity
  cs       [dble][0d0] S wave velocity
  rho      [dble][0d0] density
  phi      [dble][0d0] internal friction angle
  coh      [dble][0d0] cohesion
  Tv       [dble][0d0] visco-plastic relaxation time
  e0       [dble(3)][0d0] initial total strain (11, 22 and 12)

 ----------------------------------------------------------------------------

 NAME   : MAT_VISCO
 GROUP  : MATERIALS
 PURPOSE: Set material properties for viscoelastic medium with
          approximately constant quality factor in a prescribed frequency band
 SYNTAX : &MAT_VISCO cp,cs,rho,QP,QS,Nbody,fmin,fmax /

  cp       [dble][0d0] P wave velocity
  cs       [dble][0d0] S wave velocity
  rho      [dble][0d0] density
  QP       [dble][0d0] attenuation quality factor of P waves
  QS       [dble][0d0] attenuation quality factor of S waves
  Nbody    [int][0] number of viscoelastic mechanisms
  fmin     [dble][0d0] minimum frequency for constant Q
  fmax     [dble][0d0] maximum frequency for constant Q

 NOTE   : For Nbody=3, constant Q with less than 5% error can be achieved
          over a maximum bandwidth fmax/fmin ~ 100

 ----------------------------------------------------------------------------

 NAME   : MESH_CART
 GROUP  : MESH_DEF
 PURPOSE: Rectangular box with structured mesh.
 SYNTAX : &MESH_CART xlim, zlim, nelem, ezflt,fztag, FaultX /

  xlim     [dble(2)] [none] X limits of the box (min and max)
  zlim     [dble(2)] [none] Z limits of the box (min and max)
  nelem    [int(2)] [none]  Number of elements along each direction
  ezflt    [int][0] introduce a horizontal fault between the ezflt-th
                and the (ezflt+1)-th element rows. Rows are numbered from
                bottom to top, starting at ezflt=1.
                If ezflt=0, (default) no fault is introduced inside the box
                (for symmetric problems a fault can still be set at an external boundary)
                If ezflt=-1, a fault is introduced at/near the middle of the box
                (ezflt is reset to int[nelem(2)/2])
  fztag    [int][0] fault zone tag for elements close to the fault
                Useful to set a damping layer near the fault.
                If ezflt=0, a fault is assumed at the bottom boundary
  split    [log][F] splits the bottom boundary into two segments.
  splitD   [dble(2)][none] X distance that splits bottom boundary
  fznz     [int][1] vertical size (number of elements) of near-fault layer
  FaultX   [log][F] Same as ezflt=-1. Obsolete (will be deprecated) 

 NOTE: the following tags are automatically assigned to the boundaries: 
               1       Bottom or Bottom Right
               2       Right        
               3       Top  
               4       Left
               5       Fault, bottom side or Bottom Left
               6       Fault, top side

 ----------------------------------------------------------------------------

 NAME   : MESH_CART_DOMAIN
 PURPOSE: Define a subdomain within a structured meshed box.
 SYNTAX : &MESH_CART_DOMAIN tag,ex,ez /

  tag      [int] [none] Tag number assigned to this domain. 
  ex       [int(2)] [none]	Horizontal index of the first and last elements.
               The leftmost element column has horizontal index 1.
  ez       [int(2)] [none]	Vertical index of the first and last elements.
               The bottom element row has vertical index 1.

 NOTE   : If you ignore this input block a single domain (tag=1) will span 
          the whole box 

 ----------------------------------------------------------------------------

 NAME   : MESH_EMC2
 GROUP  : MESH_DEF
 PURPOSE: Imports a mesh from INRIA's EMC2 mesh generator in FTQ format
 SYNTAX : &MESH_EMC2 file /

  file     [name] [none] Name of the FTQ file, including suffix

 ----------------------------------------------------------------------------

 NAME   : MESH_DEF
 PURPOSE: Selects a method to import/generate a mesh.
 SYNTAX : &MESH_DEF method /
          followed by a &MESH_method input block

  method   [name] [none] Meshing method name:
               'CARTESIAN', 'LAYERED', 'EMC2', 'MESH2D'
               
 ----------------------------------------------------------------------------

 NAME   : MESH_LAYERED
 GROUP  : MESH_DEF
 PURPOSE: Structured mesh for layered medium 
          with surface and interface topography. 
 SYNTAX : &MESH_LAYERED xlim,zmin,nx,file,nlayer,ezflt,fztag /

  xlim     [dble(2)] [none] X limits of the box (min and max)
  zmin     [dble] [none] bottom Z limit of the box 
  nx       [int] [1]  Number of elements along the X direction.
                Not needed if ztopH='QSPLINE' in a &MESH_LAYER block.
  file     [string] [''] Only for flat layers,
                name of ASCII file containing layer parameters, 
                one line per layer, listed from top to bottom, 
                3 columns per line:
                (1) vertical position of top boundary,
                (2) number of elements along Z direction
                (3) material tag
  nlayer   [int] [none]  Number of layers
                If a file name is not given the layer parameters
                must be given immediately after the &MESH_LAYERED block
                by nlayer &MESH_LAYER input blocks,
                one for each layer, listed from top to bottom.
  ezflt    [int][0] introduce a fault between the ezflt-th and the
                (ezflt+1)-th element rows, numbered from bottom to top. 
                If ezflt=0 (default), no fault is introduced.
                If ezflt=-1, a horizontal fault is introduced at/near the 
                middle of the box: ezflt is reset to int[nelem(2)/2]
  fztag    [int][0] tag for elements near the fault
                Useful to set a damping layer near the fault.
  fznz     [int][1] vertical size of near-fault layer
                (half thickness in number of elements) 

 NOTE: the following tags are automatically assigned to the boundaries: 
               1       Bottom 
               2       Right        
               3       Top  
               4       Left
               5       Fault, lower side
               6       Fault, upper side

 ----------------------------------------------------------------------------

 NAME   : MESH_LAYER
 GROUP  : MESH_DEF
 PURPOSE: Define mesh parameters for one layer
 SYNTAX : &MESH_LAYER nz, ztop|ztopH, tag /
          followed by a DIST_XXXX block if ztopH is set

  nz       [int] [none]  Number of elements in layer along Z direction
  ztop     [dble] [none] Only for layers with flat top surface: 
                vertical position of top boundary
  ztopH    [string] ['none'] Name of the type of spatial distribution to 
                generate an irregular (non flat) top boundary. In general it is
                one of the 1D distribution available through a DIST_XXXX block: 
                  ztopH = 'LINEAR', or 
                  ztopH = 'SPLINE', etc. 
                There are two methods to generate a curve with a smooth normal, 
                typically to guarantee smooth boundary conditions on curved faults.
                The first method is based on quadratic splines and sometimes 
                produces degenerated elements:
                  ztopH='QSPLINE', followed by a &QC_SPLINE block
                The second method is based on cubic splines and is more robust:
                  ztopH='CSPLINE', followed by a &QC_SPLINE block
  tag      [int] [none]  Material tag
                 If not given, a tag is automatically assigned to the layer, 
                 sequentially numbered from top to bottom (top layer tag =1)

 NOTE: If ztopH='LINEAR' the mesh uses linearly deformed (Q4) elements,
       otherwise it uses quadratically deformed (Q9) elements

 ----------------------------------------------------------------------------

 NAME   : QC_SPLINE
 GROUP  : MESH_LAYER
 PURPOSE: Define the boundary of a layer using quadratic or cubic splines and
          enforcing smooth (continuous) normal between elements, for instance
          to guarantee smooth boundary conditions on curved faults.
          
 SYNTAX : &QC_SPLINE file /

  file     [string] [''] Name of ASCII file containing information of
                all the element vertex nodes lying on the boundary curve.
                One line per node, ordered by increasing x, 3 columns per line:
                  (1) x position
                  (2) z position
                  (3) derivative dz/dx of the curve at the node
                All QC_SPLINE curves in a mesh must have the same number of nodes.
                The parameter nx in &MESH_LAYERED is automatically reset
                (nx = number of nodes in QC_SPLINE - 1)

 ----------------------------------------------------------------------------

 NAME   : MESH_MESH2D
 GROUP  : MESH_DEF
 PURPOSE: Imports a mesh in mesh2d format 
	   as defined by the PRE/mesh2d mesh generator tools for Matlab
 SYNTAX : &MESH_MESH2D file /

  file     [name] [none] Name of the MESH2D file, including suffix.
                 The format of this file is:

  "NEL NPEL NNOD NBC"
  1 line with 4 integers: 
  nb of elements, nodes per element, total nb of nodes, nb of boundaries
  "NID X Y"
  NNOD lines, one per node, with 1 integer and 2 reals: 
  node id, x, y
  "EID NODES TAG"
  NEL lines, one per element, with NPEL+2 integers: 
  element id, NPEL node ids, tag. 
  "BCTAG NBEL"                                         |
  2 integers: boundary tag, nb of boundary elements    |
  "BID EID EDGE"                                       | repeat for each of
  NBEL lines, one per boundary element, 3 integers:    | the NBC boundaries
  boundary element id, bulk element id, edge id        |

 ----------------------------------------------------------------------------

 NAME   : SNAP_DEF
 GROUP  : SNAPSHOT_OUTPUTS
 PURPOSE: Set preferences for exporting snapshots
 SYNTAX : &SNAP_DEF it1, itd, fields, components, bin, visual3, avs, ps, gmt /
          Followed by a &SNAP_PS block if ps=T.

  it1      [int] [0]   Time step of first snapshot output
  itd      [int] [100] Number of timesteps between snapshots
  fields   [char*] ['V'] fields to export in snapshots (the prefix of the 
                output file names is given in parenthesis):
                 'D'     displacements (dx,dy,dz,da)
                 'V'     velocity (vx,vy,vz,va)
                 'A'     acceleration (ax,ay,az,aa)
                 'E'     strain (e11,e22,e12,e23,e13)
                 'S'     stress (s11,s22,s12,s33,e13,e23)
                 'd'     divergence rate (dvx/dx + dvz/dz)
                 'c'     curl rate (dvx/dz - dvz/dx)
  components [char*] ['ya'] components for PostScript outputs:
                 in P-SV: 'x','z' and/or 'a' (amplitude). 'y' is ignored 
                 in SH:   'y' only. Other values are ignored. 
  ps       [log] [T] PostScript (see &SNAP_PS input block)
  gmt      [log] [F] output triangulation file grid_sem2d.gmt
                 to be used in "pscontour -T" of the General Mapping Tool (GMT)
  avs      [log] [F] AVS (only for D,V and A fields)
  visual3  [log] [F] Visual3 (only for D,V and A fields)
  bin      [log] [T] binary
               
 NOTE   : E and S fields are exported only as binary.

 ----------------------------------------------------------------------------

 NAME   : SNAP_PS
 GROUP  : SNAPSHOT_OUTPUTS
 PURPOSE: Preferences for PostScript snapshots
 SYNTAX : &SNAP_PS vectors, mesh, background, color,
               isubsamp, boundaries, symbols, numbers, legend,
               ScaleField, Interpol, DisplayPts /

  vectors          [log] [F] Plots a vectorial field with arrows
  mesh             [log] [F] Plots the mesh on background
  background       [char] [''] Filled background, only for vector plots:
                                   ''   none 
                                   'P'  P-velocity model
                                   'S'  S-velocity model
                                   'T'  domains 
  isubsamp         [int] [3] Subsampling of the GLL nodes for the
                                 output of velocity model. 
                                 The default samples every 3 GLL points.
  boundaries       [log] [T] Colors every tagged boundary
  symbols          [log] [T] Plots symbols for sources and receivers
  numbers          [log] [F] Plots the element numbers
  legend           [log] [T] Writes legends
  color            [log] [T] Color output
  ScaleField       [dble] [0d0] Fixed amplitude scale (saturation),
                       convenient for comparing snapshots and making movies. 
                       The default scales each snapshot by its maximum amplitude
  Interpol         [log] [F] Interpolate field on a regular subgrid 
                       inside each element
  DisplayPts       [log] [3] Size of interpolation subgrid inside each 
                       element is DisplayPts*DisplayPts. The default plots at 
                       vertices, mid-edges and element center.
               
 ----------------------------------------------------------------------------

 NAME   : REC_LINE
 PURPOSE: Defines a line of receivers
 SYNTAX : If single receiver line: 
            &REC_LINE number,first,last,AtNode,isamp,field,irepr /
          If receiver locations from file:
            &REC_LINE file,AtNode,isamp,field,irepr /

  number   [int] [0] Number of stations in the line
  first    [dble(2)] Receivers can be located along a line,
                this is the position (x,z) of the first receiver
  last     [dble(2)] Position (x,z) of the last receiver,
                other receivers will be located with regular spacing
                between First and Last.
  file     [name] ['none'] Station positions can instead be read 
                from an ASCII file, with 2 columns: X and Z (in meters)
  AtNode   [log] [T] Relocate the stations at the nearest GLL node
  isamp    [int] [1] Sampling stride (in number of timesteps). Note that
                for stability reasons the timestep can be very small.
  field    [char] ['V'] The field in the seismogram:
                               'D'     displacement
                               'V'     velocity
                               'A'     acceleration
  irepr    [char] ['D'] Abscissa for the seismic multitrace plot:
                               'X' Horizontal position
                               'Z' Depth
                               'D' Distance to the first station

 NOTE   : to locate receivers at the free surface set their vertical position 
          above the free surface and AtNode=T

 ----------------------------------------------------------------------------

 NAME   : SRC_FORCE
 GROUP  : SOURCE MECHANISM
 PURPOSE: Point force source
 SYNTAX : &SRC_FORCE angle /

  angle    [dble] [0d0]	For P-SV, the angle of the applied force, 
                  in degrees, counterclockwise from Z-UP, e.g.: 
                  90 points left, 180 points down
                  For SH, angle is ignored and the SRC_FORCE block is not required.

 ----------------------------------------------------------------------------

 NAME   : SRC_DEF
 PURPOSE: Define the sources.
 SYNTAX : &SRC_DEF stf, mechanism, coord /
          &SRC_DEF stf, mechanism, file /
          followed by one SOURCE TIME FUNCTION block (STF_XXXX)
          and one SOURCE MECHANISM block (SRC_XXXX) 

  stf        [name] [none] Name of the source time function:
                  'RICKER', 'TAB', 'HARMONIC', 'BRUNE' or 'USER'
  mechanism  [name] [none] Name of the source mechanism:
                  'FORCE', 'EXPLOSION', 'DOUBLE_COUPLE', 'MOMENT' or 'WAVE'
  coord      [dble(2)] [huge] Location (x,z) of the source (m). 
  file       [name] ['none'] Name of file containing source parameters.
                  The file format is ASCII with one line per source and
                  2, 3 or 4 columns per line:
                    (1) X position (in m)
                    (2) Z position (in m)
                    (3) time delay (in seconds)
                    (4) relative amplitude
                  If column 4 is absent, amplitude = 1.
                  If columns 3 and 4 are absent, delay = 0 and amplitude = 1.

 ----------------------------------------------------------------------------

 NAME   : SRC_DOUBLE_COUPLE
 GROUP  : SOURCE MECHANISM
 PURPOSE: Define a double-couple source
 SYNTAX : &SRC_DOUBLE_COUPLE  dip /

  dip           [dble] [90] Dip angle, in degrees, clockwise 
                    from the positive X direction
                     
 NOTE   : Sign convention: if the source amplitude is positive the right block
          moves up (positive Z direction) in PSV and forward (positive Y 
          direction) in SH.

 NOTE   : The source time function gives the cumulative seismic moment Mo(t), 
          NOT the seismic moment rate.

 NOTE   : The seismic moment Mo must be rescaled because a 2D point source is 
          equivalent to a 3D line source. A proper scaling is obtained by 
          dividing the original 3D moment by the characteristic size of the 
          rupture area in the off-plane dimension. An approximate scaling for
          a fault area with aspect ratio close to unity is
            Mo_2D = (Mo_3D/dtau)^2/3 * dtau
          where dtau is the stress drop (typically a few MPa).

 ----------------------------------------------------------------------------

 NAME   : SRC_MOMENT
 GROUP  : SOURCE MECHANISM
 PURPOSE: Define a moment tensor source
 SYNTAX : &SRC_MOMENT Mxx,Mxz,Mzx,Mzz /
          &SRC_MOMENT Myx,Myz /

  Mxx,Mxz,Mzx,Mzz [dble] [0] Tensor components for PSV
  Myx,Myz         [dble] [0] Tensor components for SH
                     
 ----------------------------------------------------------------------------

 NAME   : SRC_WAVE
 GROUP  : SOURCE MECHANISM
 PURPOSE: Incident plane wave through the absorbing boundaries
 SYNTAX : &SRC_WAVE angle, phase /

  angle    [dble] [0d0]    Incidence angle in degrees within [-180,180]
                 counterclockwise from the positive Z (up) direction
                 to the wave vector direction:
                 Exs: incidence from below if angle in ]-90,90[
                      normal incidence from below if angle=0 
                      from bottom right if angle=+45 
                      from bottom left if angle=-45 
  phase    [char] ['S']    'S' or 'P' (only needed in PSV, ignored in SH)

 NOTE   : Incident waves enter through the absorbing boundaries.
          An incident wave is applied on every absorbing boundary
          unless "let_wave = F" in the respective BC_ABSO block.
          Incident waves are not implemented for "Stacey" absorbing boundaries.

 ----------------------------------------------------------------------------

 NAME   : STF_BRUNE
 GROUP  : SOURCE TIME FUNCTIONS
 PURPOSE: Brune (1970)'s model with omega-squared spectral fall-off:
            stf(t) = ampli*( 1 - (1+2*pi*fc*t)*exp(-2*pi*fc*t) )
 SYNTAX : &STF_BRUNE ampli, fc /

  ampli    [dble] [1d0] Amplitude (usually the seismic moment)
  fc       [dble] [1d0] Corner frequency (Hz)

 ----------------------------------------------------------------------------

 NAME   : STF_GAUSSIAN
 GROUP  : SOURCE TIME FUNCTIONS
 PURPOSE: Gaussian source time function  
          stf(t) = ampli*exp[-(pi*f0*(t-onset))^2]
 SYNTAX : &STF_GAUSSIAN ampli, f0, onset /

  ampli    [dble] [1d0] Amplitude
  onset    [dble] [0d0] Delay time (s)
  f0       [dble] [1d0] Characteristic frequency bandwidth (Hz)  

 ----------------------------------------------------------------------------

 NAME   : STF_HARMONIC
 GROUP  : SOURCE TIME FUNCTIONS
 PURPOSE: Harmonic source time function f(t) = ampli*sin(2*pi*t*f0)
 SYNTAX : &STF_HARMONIC ampli, f0 /

  ampli    [dble] [0d0] Amplitude
  f0       [dble] [0d0] Frequency

 ----------------------------------------------------------------------------

 NAME   : STF_RICKER
 GROUP  : SOURCE TIME FUNCTIONS
 PURPOSE: The Ricker wavelet is the second derivative of a gaussian.
 SYNTAX : &STF_RICKER ampli, f0, onset /

  ampli    [real] [1.] Signed amplitude of the central peak
  f0       [real >0] [0] Fundamental frequency (Hz).
                 distribution: it has a peak at f0 and an exponential
                 decay at high frequency. The cut-off high frequency is usually
                 taken as fmax = 2.5 x f0. 
  onset    [real >1/f0] [0] Delay time (secs) with respect to the peak value. 

 NOTE   : The spectrum has a peak at f0 and decays exponentially at high 
          frequencies. Beyond 2.5*f0 there is little energy, this is a 
          recommended value for fmax.
 NOTE   : onset>1/f0 is needed to avoid a strong jump at t=0, which can cause
          numerical oscillations. Ignore if using incident waves.

 ----------------------------------------------------------------------------

 NAME   : STF_TAB
 GROUP  : SOURCE TIME FUNCTIONS
 PURPOSE: Source time function spline-interpolated from values in a file 
 SYNTAX : &STF_TAB file /

  file     [string] ['stf.tab'] ASCII file containing the source time function,
               two columns: time and value. Time can be irregularly sampled and
               must increase monotonically.

 NOTE   : assumes value(t<min(time))=value(min(time)) 
          and value(t>max(time))=value(max(time))

 ----------------------------------------------------------------------------

 NAME   : STF_USER
 GROUP  : SOURCE TIME FUNCTIONS
 PURPOSE: A template for user-supplied source time function.
          File stf_user.f90 must be modified by the user to fit
          special needs.
 SYNTAX : &STF_USER ampli, onset, par1, par2, ipar1, ipar2 /

  ampli    [dble] [1.] Amplitude
  onset    [dble] [0]  Delay time (secs)
  par1     [dble] [0]  Example parameter
  par1     [dble] [0]  Example parameter
  par1     [int] [0]  Example parameter
  par1     [int] [0]  Example parameter


 ----------------------------------------------------------------------------

 NAME   : TIME
 PURPOSE: Defines time integration scheme
 SYNTAX : &TIME kind, {Dt or Courant}, {NbSteps or TotalTime} /
          Possibly followed by a TIME_XXXX block.

  kind      [char*10] ['leapfrog'] Type of scheme:
                'newmark'       Explicit Newmark
                'HHT-alpha'     Explicit HHT-alpha
                'leapfrog'      Central difference
                'symp_PV'       Position Verlet
                'symp_PFR'      Position Forest-Ruth (4th order)
                'symp_PEFRL'    Extended PFR (4th order)
		 'quasi-static'  Quasi-static (Kaneko, et al, 2011)
  Dt        [dble] [none] Timestep (in seconds)
  Courant   [dble] [0.5d0] the maximum value of the Courant-Friedrichs-Lewy 
                stability number (CFL), defined as
                  CFL = Dt*wave_velocity/dx 
                where dx is the distance between GLL nodes. Tipically CFL<= 0.5
  NbSteps   [int] [none] Total number of timesteps
  TotalTime [int] [none] Total duration (in seconds)

 NOTE   : The leap-frog scheme is recommended for dynamic faults. It is equivalent 
          to the default Newmark scheme (beta=0, gamma=1/2). However it is 
          faster and requires less memory.

 ----------------------------------------------------------------------------

 NAME   : TIME_NEWMARK
 GROUP  : TIME SCHEMES
 PURPOSE: Explicit Newmark time integration scheme 
 SYNTAX : &TIME_NEWMARK gamma, beta /

  beta     [dble] [0d0] First Newmark parameter.
               If beta=0 the scheme is fully explicit (the update of
               displacement depends only on the last value of acceleration),
               otherwise it is a single-predictor-corrector scheme
  gamma    [dble] [0.5d0] Second Newmark parameter.
               Second order requires gamma=1/2.

 ----------------------------------------------------------------------------

 NAME   : TIME_HHTA
 GROUP  : TIME SCHEMES
 PURPOSE: Explicit HHT-alpha time integration scheme, second order 
 SYNTAX : &TIME_HHTA alpha, rho /

  alpha    [dble] [0.5d0] Parameter in the HHT-alpha method. Values in [0,1].
               Defined here as 1 + HHT's original definition of alpha.
               When alpha=1 it reduces to second order explicit Newmark
               (beta=0, gamma=0.5).
  rho      [dble] [0.5d0] Minimum damping factor for high frequencies.
               Values in [0.5,1]. Rho=1 is non-dissipative.

 NOTE: We consider only second order schemes, for which alpha+gamma=3/2
       If  alpha<1, Newmark's beta is related to the HHT parameters by
         beta = 1 -alpha -rho^2*(rho-1)/[(1-alpha)*(1+rho)^3]
       If alpha=1, we set rho=1 (beta=0, gamma=0.5)
 
 NOTE: Dissipative schemes (rho<1) require slightly smaller Courant number
       (0.56 for rho=0.5, compared to 0.6 for rho=1)

 NOTE:	This is an explicit version of the HHT-alpha scheme of
         H.M. Hilber, T.J.R. Hughes and R.L. Taylor (1977) "Improved numerical 
         dissipation for time integration algorithms in structural dynamics" 
         Earthquake Engineering and Structural Dynamics, 5, 283-292
       implemented with a slightly different definition of alpha (1+original).
      	Its properties can be derived from the EG-alpha scheme of
         G.M. Hulbert and J. Chung (1996) "Explicit time integration 
         algorithms for structural dynamics with optimal numerical dissipation"
         Comp. Methods Appl. Mech. Engrg. 137, 175-188
       by setting alpha_m=0 and alpha=1-alpha_f.

