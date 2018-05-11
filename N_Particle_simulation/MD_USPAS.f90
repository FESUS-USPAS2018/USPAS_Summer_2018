program MDelectron
!unit setting: using electron@LAMMPS
!atomic mass units - Bohr radius - fs - Hartree
!!!velocity = [Bohr radius/atomic time unit] ~ 1.03275 fs
!Force = [Hartrees/Bohr radius]

  implicit none
  !system parameters
  real(8), parameter :: Pi=4*atan(1._8)
  real(8), parameter :: tt = 1.03275  !time-factor due to units
  real(8), parameter :: vc =  5.85E3  !speed of light
  real(8), parameter :: m = 5.4858E-4

  real(8), parameter :: mc2 = (m*vc)**2 !parameter for P to V
  real(8), parameter :: vc2 = vc**2 !parameter c^2 
  
  real(8), parameter :: cutoff2 = 400.0 !proximity limit
  integer, parameter :: N=8000, Ntime=10
  real(8) :: dt=1.0, realt = 0.0
  integer :: plotstride = 2
  !experiment: lz = 0.4 um, rx = ry = 100.0 um with 10^6 electrons
  !simulation: lz = 0.08 um, rx = ry = 20.0 um with 8000 electrons
  real(8), parameter :: L = 3.779E5  , Lz = 1511.8
  
  !using rest mass of electrondd
  !define the Electronic Field force coefficient in [Hartree/(a_0*electron_charge)]
  !to make F=eE agree with the Force unit
  real(8), parameter :: E_coeff=1.94E-12
  !Extraction Field in [Volts/meter]
  real(8), parameter :: EE = 5.0E6  ! for 5MV/m
  ! z_anode = 2cm to make 100keV electrons
  real(8), parameter :: z_anode_on = 1.89E8, z_anode_off = 3*z_anode_on
  real(8), parameter :: z_anode = 2*z_anode_on
  
  real(8) :: R(3,N), P(3,N), F(3,N)
  !R is for position of the atoms, P = gamma*m*v Momentum, F for Force
  real(8) :: V(3,N)
  !V array used to facilitate the p -> v -> r 
  real(8) :: PE, KE
  integer :: Time,i,j,k
  real(8) :: v2,gamma_i !v2= |v|^2 for calculating gamma
  integer :: flag(N),tout !flag: entry time for particle, tout for whole bunch
  integer :: tid,ptr,etime
  integer, dimension(:), allocatable :: emitted

   
!  open(UNIT=13, file="RandP_3D_Uniform.xyz", status="replace")
  open(UNIT=13, file="RandP_3D_Gaussian.xyz", status="replace")
  open(UNIT=11, file="init_xypxpypz.dat", action="read")
   
  !-----Initialization-----
  call init_R_Uniform(R,N)
!  call init_R_Gaussian(R,N)
!  R = 0.0
  P = 0.0
  call readin_mks(R,P,flag,dt,tout)
  call PtoV(P,V)
  F = 0.0
  
  !-----PhotoEmission-----
  !We are takke the advantage of the fact that the order of input particle number is
  !sorted by its emission time
  Allocate(emitted(tout-1))
  tid = 1
  do ptr = 1, N
    if (flag(ptr) > tid) then
      do i = tid, flag(ptr) - 1
        emitted(i) = ptr -1
      enddo
      tid = flag(ptr)
    endif
  enddo
  do etime = 1, tout-1
    call position_verlet_emitting(P,R,F,V,dt,emitted(etime))
!    print *, etime, emitted(etime)
    realt = realt + dt
  enddo
  Deallocate(emitted)

  !-----Simulation-----
  do Time=1,Ntime
    !call verlet_init(P,R,F,dt)
    call position_verlet(P,R,F,V,dt)
    realt = realt + dt
    if (modulo(Time,plotstride) == 1) then
      call getPE(R,N,PE)
      KE = 0.0
      do i = 1, N
        v2 = 0.0
        do j = 1, 3
          v2 = v2 + V(j,i)*V(j,i)
        enddo
         KE = KE + 1.0/sqrt(1.0-v2/vc2)
      enddo
      KE = (KE - N)*m*vc2 !KE = (gamma - 1)*m_0*c^2
      write(*,'(4F15.5)') realt, PE, KE, PE+KE
      write(13,'(i5)') N
      write(13,'(F15.5)') realt
      do i =1, N
        write(13,'(i5,3F15.3,3F15.7)') 1,R(:,i),P(:,i)
      end do
    end if  
    
    !change time step size for better resolution
    if (Time == 200) then
      dt = 2.5
      print *, 'dt = ', dt
    endif
    if (Time == 500) then
      dt = 5.0
      print *, 'dt = ', dt
    endif
    if (Time == 700) then
      dt = 10.0
      print *, 'dt = ', dt
    endif
    if (Time == 2000) then
      plotstride = 100
    endif
  end do

  close(13)
  close(11)

contains
 
  subroutine Init_R_Uniform(R,N)
    real(8), intent(inout) :: R(:,:)
    integer, intent(in) :: N
    integer :: numb,check,i
    real(8) :: r1,r2,r3,s1,s2,s3,rel(3)
    call random_seed()
    numb = 0
    do while (numb < N)
      call random_number(r1)
      call random_number(r2)
      call random_number(r3)
      !Uniform distribution
      s1 = 2.0*r1-1.0
      s2 = 2.0*r2-1.0
      s3 = 2.0*r3-1.0
      !check for proximity
      check = 0
      if (sqrt(s1*s1+s2*s2+s3*s3) < 1.0) then
        s1 = L*s1
        s2 = L*s2
        s3 = Lz*s3
        do i = 1, numb
          rel = R(:,i) - (/s1,s2,s3/)
          if (sum(rel**2) < cutoff2 ) then
            check = 1
            exit
          endif
        enddo
        if (check == 0) then
          numb = numb + 1
          R(:,numb) = (/s1,s2,s3/)
        endif
      endif
    enddo
  end subroutine Init_R_Uniform

  subroutine Init_R_Gaussian(R,N)
    real(8), intent(inout) :: R(:,:)
    integer, intent(in) :: N
    integer :: numb,check,i
    real(8) :: r1,r2,r3,r4,s1,s2,s3,rel(3)
    call random_seed()
    numb = 0
    do while (numb < N)
      call random_number(r1)
      call random_number(r2)
      call random_number(r3)
      call random_number(r4)
      !Box-Muller for Gaussian distribution
      s1 = L*sqrt(-2*log(r1))*cos(2*pi*r2)
      s2 = L*sqrt(-2*log(r1))*sin(2*pi*r2)
      s3 = Lz*sqrt(-2*log(r3))*sin(2*pi*r4)
      !check for proximity
      check = 0
      do i = 1, numb
        rel = R(:,i) - (/s1,s2,s3/)
        if (sum(rel**2) < cutoff2 ) then
          check = 1
          exit
        endif
      enddo
      if (check == 0) then
        numb = numb + 1
        R(:,numb) = (/s1,s2,s3/)
      endif
    enddo
  end subroutine Init_R_Gaussian

  subroutine readin_mks(R,P,flag,dt,tout)
    real(8), intent(inout) :: P(:,:), R(:,:)
    integer, intent(inout) :: flag(:),tout   
    real(8), intent(in) :: dt
    real(8) :: t,x,y,px,py,pz,t_fs,delta_t,tick_shift
    real(8) :: m_to_a0 = 1.89E10
    !for photoemission, p ~ mv
    !drifting delta_x[a_0] = px*drift_converter*dt = px/m_e[kg]*1.0E-15*m_to_a0,
    real(8) :: drift_converter = 2.0746E25
    !momentum converter,p[kg*m/s]*p_converter = p[amu*a0/tau]
    real(8):: p_converter = 1.175E22
    integer :: i
    
    do i=1,N
      read(11,*) t,x,y,px,py,pz
      t_fs = t*1.0E15
      flag(i) = ceiling(t_fs/dt)
      delta_t = dt*flag(i) - t_fs
      R(1,i) = x*m_to_a0 + px*drift_converter
      R(2,i) = y*m_to_a0 + py*drift_converter
      R(3,i) = pz*drift_converter
      P(:,i) = (/px,py,pz/)*p_converter
    enddo
    tick_shift = MINVAL(flag) - 1
    flag = flag - tick_shift
    tout = MAXVAL(flag)
  end subroutine readin_mks

  subroutine verlet_init(P,R,F,dt)
    real(8), intent(inout) :: P(:,:), R(:,:), F(:,:)
    real(8), intent(in) :: dt
    real(8) :: dtttm, hdttt
    dtttm = dt*tt/m
    hdttt = 0.5*dt*tt
      
    P=P+F*hdttt
    R=R+P*dtttm
    call Force(F,R)
    P=P+F*hdttt
  end subroutine verlet_init
 
  subroutine position_verlet(P,R,F,V,dt)
    real(8), intent(inout) :: P(:,:), R(:,:), F(:,:), V(:,:)
    real(8), intent(in) :: dt
    real(8) :: hdttt, dttt
    hdttt = 0.5*dt*tt
    dttt = dt*tt
      
    R=R+V*hdttt
    call Force(F,R)
    P=P+F*dttt
    call PtoV(P,V)
    R=R+V*hdttt
  end subroutine position_verlet

  subroutine position_verlet_emitting(P,R,F,V,dt,e_num)
    real(8), intent(inout) :: P(:,:), R(:,:), F(:,:), V(:,:)
    real(8), intent(in) :: dt
    integer, intent(in) :: e_num
    integer :: i
    real(8) :: hdttt, dttt
    hdttt = 0.5*dt*tt
    dttt = dt*tt
      
    do i = 1, e_num
      R(:,i)=R(:,i)+V(:,i)*hdttt
    enddo
    call Force_emitting(F,R,e_num)
    do i = 1, e_num
      P(:,i)=P(:,i)+F(:,i)*dttt
    enddo
    call PtoV_emitting(P,V,e_num)
    do i = 1, e_num
      R(:,i)=R(:,i)+V(:,i)*hdttt
    enddo
  end subroutine position_verlet_emitting

  subroutine PtoV(P,V)
    real(8), intent(in) :: P(:,:)
    real(8), intent(out):: V(:,:)
    integer :: i,j
    real(8) :: p2,pvc
    do i = 1, N
      p2 = 0.0
      do j= 1, 3
        p2 = p2 + p(j,i)**2
      enddo
      pvc = 1.0/(m*sqrt(1+p2/mc2))
      do j= 1, 3
        V(j,i) = P(j,i)*pvc
      enddo
    enddo
  end subroutine  PtoV

  subroutine PtoV_emitting(P,V,e_num)
    real(8), intent(in) :: P(:,:)
    real(8), intent(out):: V(:,:)
    integer, intent(in) :: e_num
    integer :: i,j
    real(8) :: p2,pvc
    do i = 1, e_num
      p2 = 0.0
      do j= 1, 3
        p2 = p2 + p(j,i)**2
      enddo
      pvc = 1.0/(m*sqrt(1+p2/mc2))
      do j= 1, 3
        V(j,i) = P(j,i)*pvc
      enddo
    enddo
  end subroutine  PtoV_emitting

  !Force Calculation
  subroutine Force(F,R)
    real(8), intent(inout) :: F(:,:)
    real(8), intent(in)   :: R(:,:)
    integer :: i,j 
    real(8) :: REL(3), RELlength, Fij(3), dz_ee
    real(8) :: zforee !this create a register for R(3,i) to save time
    F=0.0
    !$omp parallel do private(fij,rel,Rellength)
    do j=1,N-1
      do i=j+1,N
        REL=R(:,i)-R(:,j)
        !REL=REL-nint(Rel/L)*L
        RELlength=sqrt(sum(REL**2))
        !FC = 1.0 so we can skip it
        Fij=(1/(RELlength)**3)*REL
        !$omp atomic
        F(1,j)=F(1,j)-Fij(1)
        !$omp atomic
        F(2,j)=F(2,j)-Fij(2)
        !$omp atomic
        F(3,j)=F(3,j)-Fij(3)
        !$omp atomic
        F(1,i)=F(1,i)+Fij(1)
        !$omp atomic
        F(2,i)=F(2,i)+Fij(2)
        !$omp atomic
        F(3,i)=F(3,i)+Fij(3)
      end do 
    end do 
    !$omp end parallel do

    !the Extraction Field on z-dirction
!    do i = 1, N
!      zforee = R(3,i)
!      if ( zforee < z_anode_off) then 
!        if (zforee < z_anode_on) then
!          F(3,i) = F(3,i) + E_coeff*EE 
!        else
!          dz_ee = (zforee - z_anode)/z_anode_on
!          F(3,i) = F(3,i) + E_coeff*0.5*EE*(1.0-tanh(5.0*dz_ee))
!        endif
!      endif
!    enddo

  end subroutine Force

  subroutine Force_emitting(F,R,e_num)
    real(8), intent(inout) :: F(:,:)
    real(8), intent(in)   :: R(:,:)
    integer, intent(in) :: e_num
    integer :: i,j 
    real(8) :: REL(3), RELlength, Fij(3), dz_ee
    real(8) :: zforee !this create a register for R(3,i) to save time
    F=0.0
    !$omp parallel do private(fij,rel,Rellength)
    do j=1,e_num-1
      do i=j+1,e_num
        REL=R(:,i)-R(:,j)
        !REL=REL-nint(Rel/L)*L
        RELlength=sqrt(sum(REL**2))
        !FC = 1.0 so we can skip it
        Fij=(1/(RELlength)**3)*REL
        !$omp atomic
        F(1,j)=F(1,j)-Fij(1)
        !$omp atomic
        F(2,j)=F(2,j)-Fij(2)
        !$omp atomic
        F(3,j)=F(3,j)-Fij(3)
        !$omp atomic
        F(1,i)=F(1,i)+Fij(1)
        !$omp atomic
        F(2,i)=F(2,i)+Fij(2)
        !$omp atomic
        F(3,i)=F(3,i)+Fij(3)
      end do 
    end do 
    !$omp end parallel do

    !the Extraction Field on z-dirction
!    do i = 1, e_num
!      zforee = R(3,i)
!      if ( zforee < z_anode_off) then 
!        if (zforee < z_anode_on) then
!          F(3,i) = F(3,i) + E_coeff*EE 
!        else
!          dz_ee = (zforee - z_anode)/z_anode_on
!          F(3,i) = F(3,i) + E_coeff*0.5*EE*(1.0-tanh(5.0*dz_ee))
!        endif
!      endif
!    enddo
  end subroutine Force_emitting

  !get the potential energy of Atomic system
  subroutine getPE(R,N,PE)
  real(8), intent(in) :: R(:,:)
  integer, intent(in) :: N
  real(8), intent(out):: PE
  integer :: i,j
  real(8) :: rel_pe(3),rel_2,vij
  PE = 0.0
  do j=1,N-1
    do i=j+1,N
      rel_pe=R(:,i)-R(:,j)
      rel_2=sqrt(sum(rel_pe**2))
      vij = 1.0/(rel_2)
      PE = PE + vij
    end do
  end do
  end subroutine getPE

end program MDelectron
