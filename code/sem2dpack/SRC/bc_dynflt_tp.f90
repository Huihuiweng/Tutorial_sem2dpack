! Thermal pressurization
! The main part of code was modified from MDSBI version 4.1.9 by Eric Dunham, https://pangea.stanford.edu/~edunham/codes/codes.html
! First implemented: 8 March 2024
! By Huihui Weng and Chao Liang

module bc_dynflt_tp

  use distribution_cd

  implicit none
  private

  integer,parameter :: pin = selected_int_kind(9)
  integer,parameter :: pr  = selected_real_kind(15,307)
  real(pr),parameter :: zero = 0._pr
  real(pr),parameter :: one = 1._pr
  real(pr),parameter :: two = 2._pr
  real(pr),parameter :: half = 0.5_pr
  real(pr),parameter :: pi = 3.14159265358979323846264338327950288419716939937510_pr
  real(pr),parameter :: twopi = 2._pr*pi

  type tp_input_type
    type(cd_type) :: rhoc,K,alphath,alphahy,Lambda,W,beta,deltaD,Phi 
  end type tp_input_type


  type tp_type
     private
     double precision, dimension(:,:), pointer :: T=>null(), dT=>null(), p=>null(), dp=>null()
     double precision, dimension(:), pointer   :: rhoc=>null(), K=>null(), alphath=>null(), alphahy=>null(), Lambda=>null()
     double precision, dimension(:), pointer   :: W=>null(), beta=>null(), deltaD=>null(), Phi=>null(), z=>null()
     real(pr)      :: T0, dz, dtime, sub_dt
     character(64) :: shape
     integer(pin)  :: nz

     type(tp_input_type) :: input
  end type tp_type

  public :: tp_type, tp_read, tp_init, thermpres_rate, getPorepressure, getTemperature 

  ! THERMPRES_TYPE is a derived type containing variables related to 
  ! thermal pressurization friction law
  !
  ! Parameters:
  ! RHOC = volumetric heat capacity
  ! K = thermal conductivity
  ! ALPHATH = thermal diffusivity
  ! ALPHAHY = hydraulic diffusivity
  ! LAMBDA = pore pressure rise per unit temprature rise
  ! SHAPE = shape of strain distribution
  ! W = width of deformation zone
  ! PHI = Maximum inelastic porosity change at large slip
  ! DELTAD = Characteristic slip distance associated with the change in porosity
  ! BETA =  volumetric pore fluid storage coefficient beta = n*(beta_n+beta_f)

contains

  subroutine tp_read(tp,ninput)
    ! READ_THERMPRES reads in thermpres variables from file

    use echo, only : echo_input,iout
    use stdio, only: IO_abort

    implicit none
    integer(pin),intent(in) :: ninput
    type(tp_type),intent(out) :: tp
    
    integer(pin)  :: nz
    real(pr)      :: rhoc,K,alphath,alphahy,Lambda,W,beta,deltaD,Phi
    real(pr)      :: T0, dz
    character(20) :: rhocH,KH,alphathH,alphahyH,LambdaH,WH
    character(20) :: betaH,deltaDH,PhiH
    character(64) :: shape

    ! make namelist of user input variables
    namelist /BC_DYNFLT_TP/ shape,T0,dz,nz,&
                            rhoc,K,alphath,alphahy,Lambda,W,beta,deltaD,Phi,&
                            rhocH,KH,alphathH,alphahyH,LambdaH,WH,betaH,deltaDH,PhiH

    ! assign default values
    shape   = 'plane'
    T0      = zero
    dz      = zero
    nz      = one
    rhoc    = zero
    K       = zero
    alphath = one
    alphahy = one
    Lambda  = one
    W       = half
    beta    = one
    deltaD  = one
    Phi     = zero
    rhocH   = ""
    KH      = ""
    alphathH= ""
    alphahyH= ""
    LambdaH = ""
    WH      = ""
    betaH   = ""
    deltaDH = ""
    PhiH    = ""

    ! read namelist from input file, call error routine if needed
    read(ninput,BC_DYNFLT_TP,END=300)
    300 continue

    tp%shape   = shape

    ! user only specifies one of rhoc and K (rhoc takes precedence if both given)
    if (rhoc/=zero) then
       K    = rhoc*alphath
    else
       rhoc = K/alphath
    end if

    call DIST_CD_Read(tp%input%rhoc,rhoc,rhocH,ninput,rhocH)
    call DIST_CD_Read(tp%input%K,      K,   KH,ninput,KH)
    call DIST_CD_Read(tp%input%alphath,alphath,alphathH,ninput,alphathH)
    call DIST_CD_Read(tp%input%alphahy,alphahy,alphahyH,ninput,alphahyH)
    call DIST_CD_Read(tp%input%Lambda,  Lambda, LambdaH,ninput,LambdaH)
    call DIST_CD_Read(tp%input%W,            W      ,WH,ninput,WH)
    call DIST_CD_Read(tp%input%beta,beta,betaH,ninput,betaH)
    call DIST_CD_Read(tp%input%deltaD,deltaD,deltaDH,ninput,deltaDH)
    call DIST_CD_Read(tp%input%Phi,Phi,PhiH,ninput,PhiH)

    tp%T0 = T0
    tp%dz = dz
    tp%nz = nz

  if (echo_input) write(iout,400) nz,dz,T0,rhocH,KH,alphathH,alphahyH,LambdaH,WH,betaH,deltaDH,PhiH,shape
  return

  400 format(5x,'Friction law  . . . . . . . . . . . . . .  = thermal pressurization', &
            /5x,'  Number of finite difference elements . (nz) = ',I,&
            /5x,'  Size of finite difference elements . . (dz) = ',F,&
            /5x,'  Initial temperature. . . . . . . . . . (T0) = ',F,&
            /5x,'  Volumetric heat capacity  . . . . .  (rhoc) = ',A,&
            /5x,'  Thermal conductivity . . . . . . . . . .(K) = ',A,&
            /5x,'  Thermal diffusivity . . . . . . . (alphath) = ',A,&
            /5x,'  Hydraulic diffusivity  . . . . . .(alphahy) = ',A,&
            /5x,'  Pore pressure rise  . . . . . . . .(lambda) = ',A,&
            /5x,'  Width of deformation zone  . . . . . . .(W) = ',A,&
            /5x,'  Maximum inelastic porosity  . . . . . (Phi) = ',A,&
            /5x,'  Characteristic slip distance . . . (DeltaD) = ',A,&
            /5x,'  Volumetric fluid storage coefficient (beta) = ',A,&
            /5x,'  Shape of strain distribution . . . .(shape) = ',A)

  end subroutine tp_read


  subroutine tp_init(tp,coord,dt)
    ! INIT_THERMPRES initializes thermpres variables

    implicit none
    
    type(tp_type),intent(inout)   :: tp
    double precision, intent(in)  :: coord(:,:)
    double precision, intent(in)  :: dt
    integer(pin) :: k,darray
    integer(pin) :: ny

    tp%dtime = dt

    ! Number of fault nodes
    ny = size(coord,2) 

    ! allocate and initialize diffusion grids
    allocate(tp%z( tp%nz))
    allocate(tp%T( tp%nz,ny))
    allocate(tp%dT(tp%nz,ny))
    allocate(tp%p( tp%nz,ny))
    allocate(tp%dp(tp%nz,ny))
    
    ! initialize temperature with linear profile and pressure as zero
    tp%z  = real( (/ (k, k=0,tp%nz-1) /) ,pr)*tp%dz
    tp%T  = zero
    tp%dT = zero
    tp%p  = zero
    tp%dp = zero

    ! initialize fault zone properties
    call DIST_CD_Init(tp%input%rhoc,   coord, tp%rhoc)
    call DIST_CD_Init(tp%input%K,      coord, tp%K)
    call DIST_CD_Init(tp%input%alphath,coord, tp%alphath)
    call DIST_CD_Init(tp%input%alphahy,coord, tp%alphahy)
    call DIST_CD_Init(tp%input%Lambda, coord, tp%Lambda)
    call DIST_CD_Init(tp%input%W,      coord, tp%W)
    call DIST_CD_Init(tp%input%beta,   coord, tp%beta)
    call DIST_CD_Init(tp%input%deltaD, coord, tp%deltaD)
    call DIST_CD_Init(tp%input%Phi,    coord, tp%Phi)

    ! add constant background
    do k = 1,tp%nz
       tp%T(k,:) = tp%T(k,:)+tp%T0
    end do

  end subroutine tp_init


!  subroutine destroy_thermpres(tp)
!    ! DESTROY_THERMPRES destroys derived type tp
!
!    implicit none
!    type(tp_type),intent(inout) :: tp
!
!    if (allocated(tp%T))  deallocate(tp%T)
!    if (allocated(tp%dT)) deallocate(tp%dT)
!    if (allocated(tp%p))  deallocate(tp%p)
!    if (allocated(tp%dp)) deallocate(tp%dp)
!
!  end subroutine destroy_thermpres


  subroutine thermpres_rate(tp,coord,V,Tau,D)
    ! THERMPRES_RATE solves the 1D diffusion equations for thermal pressurization
    ! using explicit finite differences

    implicit none
    type(tp_type),intent(inout) :: tp
    double precision, intent(in)    :: coord(:,:)
    double precision, dimension(:), intent(in) :: V,Tau,D

    double precision, dimension(size(V)) :: Q
    double precision  :: dt, sub_dt, FTCS
    integer(pin)      :: i, j, ny, sub_steps

    ny = size(coord,2)
    dt = tp%dtime

    ! Use FTCS scheme to stabilize the finite difference solution
    ! i.e., delta_t < delta_z**2 / (2*max(alphahy,alphath))
    FTCS = max(maxval(tp%alphath)*dt/tp%dz**2, maxval(tp%alphahy)*dt/tp%dz**2)
    if (FTCS < 0.2) then
       sub_steps = one
       sub_dt    = dt
    else
       sub_steps =  ceiling(FTCS/0.2)
       sub_dt    =  dt/sub_steps
    endif

    ! shear heating
    Q = V * Tau

    do j = 1,sub_steps
       do i = 1,ny
          call thermpres_rate_point(tp%dT(:,i),tp%dp(:,i),tp%T(:,i),tp%p(:,i),Q(i),V(i),D(i),i,tp)
       end do

       ! Update T and P
       tp%T = tp%T + tp%dT*sub_dt
       tp%p = tp%p + tp%dp*sub_dt
    end do

  end subroutine thermpres_rate


  subroutine thermpres_rate_point(dT,dp,T,p,Q,V,D,y_ind,tp)
    ! THERMPRES_RATE_POINT solves the 1D diffusion equations for thermal pressurization
    ! using explicit finite differences at one point on fault

    use echo, only : echo_input,iout
    use stdio, only: IO_abort

    implicit none

    real(pr),dimension(:),intent(out) :: dT,dp
    real(pr),dimension(:),intent(in) :: T,p
    real(pr),intent(in) :: Q
    real(pr),intent(in) :: V,D
    type(tp_type),intent(in) :: tp
    integer(pin),intent(in)  :: y_ind

    ! Internal Parameters:
    ! nz = length of T,p
    ! RT = alpha_th/dz**2
    ! RP = alpht_hy/dz**2
    ! DTDZ = derivative of T with respect to z at left boundary
    ! G = weighting function of heat source
    ! Z = distance normal to fault
    integer(pin) :: i,nz
    real(pr) :: dTdz,rt,rp,rd
    real(pr),dimension(size(T)) :: G

    ! The TP propertis at each fault node
    real(pr)      :: rhoc,K,alphath,alphahy,Lambda,W,beta,deltaD,Phi
    real(pr)      :: dz 

    rhoc     = tp%rhoc(y_ind)
    K        = tp%K(y_ind)
    alphath  = tp%alphath(y_ind)
    alphahy  = tp%alphahy(y_ind)
    Lambda   = tp%Lambda(y_ind)
    W        = tp%W(y_ind)
    beta     = tp%beta(y_ind)
    deltaD   = tp%deltaD(y_ind)
    Phi      = tp%Phi(y_ind)
    dz       = tp%dz

    nz   = size(T)
    ! Initialization
    dT   = zero
    dp   = zero
    dTdz = zero

    ! shear heating, either as boundary condition for slip-on-plane model or as
    ! volumetric heat source term for finite width shear zone model

    select case(tp%shape)
    case default
         call IO_abort('Invalid shear zone shape: thermpres_rate')
    case('plane')
        ! heat flux on the fault plane: no volumetric source term
        dTdz = -half*Q/K
        G    = 0

    case('box')
        ! volumetric source term
        ! calculate nondimensionalized strain rate, G=W*(strain rate)/V,
        ! note that integral over z of G(z) must equal W, in order that
        ! integral over z of strain rate equals V

        ! calculate coefficient to properly normalize strain rate distribution
        ! Accounting the symmetry, W and G are halved
        do i = 1,nz
           if (tp%z(i)<half*W-half*dz) then
               G(i) = half
           elseif (tp%z(i)>=half*W-half*dz.and.tp%z(i)<half*W+half*dz) then
               G(i) = (half*W+half*dz-tp%z(i))/dz
           else
               G(i) = zero
           end if
        end do

        ! Add the shear heating term
        dT = G*Q/(W*rhoc)

     case('gaussian')
        ! W is standard deviation of gaussian
        do i = 1,nz
            G(i) = one/sqrt(twopi) * exp(-tp%z(i)**2/(two*W**2)) 
        end do

        ! Add the shear heating term
        dT = G*Q/(W*rhoc)

    end select

    ! special case of adiabatic, undrained response (only works for finite width shear zone)
    if (nz==1) then
       dp = Lambda*dT
       return
    else
    end if

    ! add contribution from diffusion of heat and fluid mass
    ! parameter combinations
    rt = alphath/dz**2
    rp = alphahy/dz**2
    rd = (1/beta)*Phi*(V/deltaD)*exp(-D/deltaD)
        
    ! diffusion terms in interior
    do i = 2,nz-1
       dT(i) = dT(i)+rt*(T(i+1)-two*T(i)+T(i-1))
       dp(i) = dp(i)+rp*(p(i+1)-two*p(i)+p(i-1)) - rd*G(i)
    end do

    ! and at left boundary
    dT(1) = dT(1)+rt*two*(T(2)-T(1)  - dz*dTdz)
    dp(1) = dp(1)+rp*two*(p(2)-p(1)) - rd*G(1)  ! dp/dz = 0 on fault (no fluid flow across it by symmetry)

    ! thermal pressurization term in pressure rate
    dp = dp+Lambda*dT

    ! fix T and p at right boundary
    dT(nz) = zero
    dp(nz) = zero

  end subroutine thermpres_rate_point



!===============================================================
  function getPorepressure(tp) result(p)

  type(tp_type),intent(in)   :: tp
  real(pr),dimension(:),allocatable :: p

  allocate(p(size(tp%p,2)))
  p = tp%p(1,:)

  end function getPorepressure

!--------------------------------
  function getTemperature(tp) result(T)

  type(tp_type),intent(in)   :: tp
  real(pr),dimension(:),allocatable :: T

  allocate(T(size(tp%T,2)))
  T = tp%T(1,:)

  end function getTemperature


end module bc_dynflt_tp

