c===================================================================
c
c  1. This UMAT subroutine implements a hybrid rate dependent 
c     polycrystal plasticity theory. It is applicable to materials
c     with strain rate sensitivity exponents < 0.1.
c  2. See the VERY IMPORTANT comments in the "flow_rule" subroutine.
c  3. There are no limitations on the isotropic and kinematic
c     hardening laws.
c  4. The calculations are performed in each grain's coordinate
c     system because that is more efficient.
c  5. Subroutines for the flow rule, isotropic hardening law, 
c     kinematic hardening law, and Jacobian are at the bottom
c     of this file.
c  6. By default, the following are preprogrammed into this UMAT:
c     FCC slip systems, power law flow rule, direct hardening /
c     dynamic recovery hardening laws for iso and back stresses.
c  7. There is a commented-out write statement that can be activated
c     to write out Euler angles in Roe convention around line 890.
c  8. Details of the algorithm are given in my Ga Tech PhD thesis
c     and in a submitted paper.
c
c     Bob McGinty
c     bmcginty@merc-mercer.org
c
c===================================================================

      subroutine umat (stress,  statev,  ddsdde,  sse,     spd,
     &			scd,     rpl,     ddsddt,  drplde,  drpldt,
     &			strain,  dstrain, time,    dtime,   temp,
     &			dtemp,   predef,  dpred,   cmname,  ndi,
     &			nshr,    ntens,   nstatv,  props,   nprops,
     &			coords,  drot,    pnewdt,  celent,  dfgrd0,
     &			dfgrd1,  noel,    npt,     layer,   kspt,
     &			kstep,   kinc)

      implicit double precision (a-h,o-z)

      parameter(Max_active_sys =  5,        ! 3-D => 5,   2-D => 2
     &          ISVs_per_grain = 51,
     &          Num_slip_sys   = 12,
     &          Tolerance      = 0.001,
     &          Tiny           = 1.D-30,
     &          Pi             = 3.14159265358979 )

c-------------------------------------------------------------------
c   Determine number of grains.  Divide the total number of ISVs
c   by the number of ISVs per grain.
c   3 - Xtal Orientation (Euler Angles)
c   6 - Iso  Stress in Crystal Orientation
c   6 - Back Stress in Crystal Orientation
c  12 - gamma_dot(i)     for 12 slip systems
c  12 - g(i) Iso  stress for 12 slip systems
c  12 - x(i) Back stress for 12 slip systems
c ----
c  51 
c-------------------------------------------------------------------

      character*8 cmname
      logical Converged, Out_of_time
      external transpose

c-------------------------------------------------------------------
c  Dimension arrays passed into the UMAT sub
c-------------------------------------------------------------------

      dimension
     &	coords(3),	! Coordinates of Gauss pt. being evaluated
     &	ddsdde(ntens,ntens), ! Tangent Stiffness Matrix
     &	ddsddt(ntens),	! Change in stress per change in temperature
     &	dfgrd0(3,3),	! Deformation gradient at beginning of step
     &	dfgrd1(3,3),	! Deformation gradient at end of step
     &	dpred(1),	! Change in predefined state variables
     &	drplde(ntens),	! Change in heat generation per change in strain
     &	drot(3,3),	! Rotation matrix
     &	dstrain(ntens),	! Strain increment tensor stored in vector form
     &	predef(1),	! Predefined state vars dependent on field variables
     &	props(nprops),	! Material properties passed in
     &	statev(nstatv),	! State Variables
     &	strain(ntens),	! Strain tensor stored in vector form
     &	stress(ntens),	! Cauchy stress tensor stored in vector form
     &	time(2)		! Step Time and Total Time
                                                                                
c-------------------------------------------------------------------            
c  Dimension other arrays used in this UMAT
c-------------------------------------------------------------------            

      dimension
     &	array1	(3,3),		! Dummy array
     &	array2	(3,3),		! Dummy array
     &	array3	(3,3),		! Dummy array
     &	array4	(6,6),		! Dummy array used in Voigt notation
     &	array5	(6,6),		! Inverse of array4().
     &	array6	(3,3,3,3),	! 4th rank dummy array
     &	array7	(3,3,3,3),	! Another 4th rank dummy array
     &  a	(3,3),		! Back stress tensor
     &	D	(3,3),		! Rate of deformation tensor
     &	D_el	(3,3),		! Elastic Rate of deformation tensor
     &	D_p	(3,3),		! Plastic Rate of deformation tensor
     &	D_xtal	(3,3),		! Rate of deformation tensor in xtal frame
     &	del	(3,3),		! Kronecker delta tensor
     &	delta_tau(num_slip_sys),! Change in tau_critical during step
     &	ddpdsig	(3,3,3,3),	! deriv of D_p wrt sig * dt
     &	ddsdde_4th(3,3,3,3),	! 4th rank tangent stiffness tensor
     &	dir_cos	(3,3),		! Direction cosines 
     &	g0	(num_slip_sys), ! Iso stress at beginning of step
     &	g	(num_slip_sys), ! Iso stress at end of step
     &	gamma_check(num_slip_sys),! shear strain rate on system
     &	gamma_dot(num_slip_sys),! shear strain rate on system
     &	icode	(num_slip_sys), ! = -1 if active, or >0 if simultaneous
     &	id	(5),		! ID array of active slip systems
     &	PdP	(5,5),		! P double dot P
     &	PdCdP	(5,5),		! P double dot C double dot P
     &	PdCdD	(5),		! P double dot C double dot D
     &	psi	(3),		! Euler angles
     &	p	(3,3,num_slip_sys),!Symmetric Schmid tensor
     &	pm	(3,3,num_slip_sys),!"merged"Symmetric Schmid tensor
     &	q	(3,3,num_slip_sys),!Antisymmetric part of Schmid tensor
     &	qm	(3,3,num_slip_sys),!"merged"Antisym part of Schmid tensor
     &	Schmid	(3,3,num_slip_sys),!Schmid tensor in xtal frame
     &	sig	(3,3),		! Stress tensor at end of step
     &	sig_avg	(3,3),		! Average stress
     &	sig_dot	(3,3),		! Rate of stress change for a grain
     &	tau	(num_slip_sys), ! Resolved shear stress 
     &	tau_c	(num_slip_sys), ! Critical tau 
     &	tau_dot	(num_slip_sys), ! Rate of resolved shear stress change
     &	time0	(num_slip_sys), ! Time required to activate slip sys
     &	vector	(num_slip_sys), ! Dummy working vector
     &	W_p	(3,3),		! Plastic spin in updated xtal frame
     &	W_el	(3,3),		! Substructure spin in xtal frame
     &	W_xtal	(3,3),		! Continuum Spin in xtal frame
     &	x	(num_slip_sys), ! Slip system back stress at end of step
     &	x0	(num_slip_sys), ! Slip system back stress at beginning
     &	xL_global(3,3),		! Eulerian velocity gradient
     &	xL_p	(3,3),		! Plastic part of Eulerian velocity gradient
     &	xL_xtal	(3,3),		! Velocity grad in crystal configuration
     &	xtal_rot_matrix(3,3),	! Xtal rotation during step
     &	y	(3,num_slip_sys),! Miller indices of slip directions 
     &	z	(3,num_slip_sys) ! Miller indices of slip plane normals


c-------------------------------------------------------------------
c  Divide #ISVs by #ISVs_per_grain to get number of grains
c-------------------------------------------------------------------

      num_grains = nstatv / ISVs_per_grain

c-------------------------------------------------------------------
c  Assign props() array to logical variable names
c-------------------------------------------------------------------

      C11        = props(1)   ! Elastic tensor component
      C12        = props(2)   ! Elastic tensor component
      C44        = props(3)   ! Elastic tensor component
      gamma_zero = props(4)   ! Gamma_dot_zero in flow rule
      flow_exp   = props(5)   ! Flow Exponent  in flow rule 
      G_zero     = props(6)   ! Iso hardening value at t=0
      G_dir      = props(7)   ! Direct  Hardening Constant (iso stress)
      G_dyn      = props(8)   ! Dynamic Recovery  Constant (iso stress)
      xLatent    = props(9)   ! Latent Hardening Ratio
      X_dir      = props(10)  ! Direct  Hardening Constant (back stress)
      X_dyn      = props(11)  ! Dynamic Recovery  Constant (back stress)

c-------------------------------------------------------------------
c  Example values to use as a guide for getting started
c-------------------------------------------------------------------

c      C11        = 150,000
c      C12        =  75,000
c      C44        =  37,500
c      gamma_zero =  0.0001
c      flow_exp   =  0.01
c      G_zero     =  13
c      G_dir      = 200
c      G_dyn      =   2
c      xLatent    = 1.0
c      X_dir      = 100
c      X_dyn      =  10

c-------------------------------------------------------------------
c  Initialize Kronnecker delta tensor
c-------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          del(i,j) = 0.0
        end do
        del(i,i) = 1.0
      end do

c-------------------------------------------------------------------
c  Assign slip system normals and slip directions for an FCC.
c-------------------------------------------------------------------

c     !   plane      dir
c     /  1, 1, 1, -1, 1, 0 /
c     /  1, 1, 1,  0,-1, 1 /
c     /  1, 1, 1,  1, 0,-1 /
c     / -1, 1, 1, -1,-1, 0 /
c     / -1, 1, 1,  1, 0, 1 /
c     / -1, 1, 1,  0, 1,-1 /
c     / -1,-1, 1,  1,-1, 0 /
c     / -1,-1, 1,  0, 1, 1 /
c     / -1,-1, 1, -1, 0,-1 /
c     /  1,-1, 1,  1, 1, 0 /
c     /  1,-1, 1, -1, 0, 1 /
c     /  1,-1, 1,  0,-1,-1 /


      do j = 1,num_slip_sys
         do i = 1,3
           y(i,j) = 0
           z(i,j) = 1
         end do
      end do
      do j = 4,9
         z(1,j) = -1
      end do
      do j = 7,12
         z(2,j) = -1
      end do
      y(2,1) = 1
      y(3,2) = 1
      y(1,3) = 1
      y(1,5) = 1
      y(3,5) = 1
      y(2,6) = 1
      y(2,8) = 1
      y(3,8) = 1
      y(1,7) = 1
      y(1,10) = 1
      y(2,10) = 1
      y(3,11) = 1
      y(1,1) = -1
      y(2,2) = -1
      y(3,3) = -1
      y(1,4) = -1
      y(2,4) = -1
      y(3,6) = -1
      y(1,9) = -1
      y(3,9) = -1
      y(2,7) = -1
      y(1,11) = -1
      y(2,12) = -1
      y(3,12) = -1

c-------------------------------------------------------------------
c  Normalise the Miller indices to length one.
c-------------------------------------------------------------------

      do i = 1,num_slip_sys
        call normalize_vector( y(1,i), y(2,i), y(3,i) )
        call normalize_vector( z(1,i), z(2,i), z(3,i) )
      end do

c--------------------------------------------------------------------
c  Calculate Schmid tensor in xtal orientation.
c--------------------------------------------------------------------

      do n = 1,num_slip_sys
        do i = 1,3
          do j = 1,3 
            Schmid(i,j,n) = y(i,n) * z(j,n)
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Calculate symmetric and antisymmetric parts of Schmid tensor
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
        do i = 1,3
          do j = 1,3
            p(i,j,k) = 0.5 * (Schmid(i,j,k) + Schmid(j,i,k))
            q(i,j,k) = 0.5 * (Schmid(i,j,k) - Schmid(j,i,k))
          end do
        end do
      end do

c===================================================================
c  INITIALIZE INTERNAL VARIABLES FOR INITIAL TIME STEP
c===================================================================

      if (time(2) .eq. 0.0) then

c-------------------------------------------------------------------
c  Check for normality of Miller indices.
c-------------------------------------------------------------------

        do k = 1,num_slip_sys
           sum = 0.0
          do i = 1,3
            sum = sum + y(i,k) * z(i,k)
          end do
          if (abs(sum) .gt. tolerance) then
            print*,'The Miller indices are WRONG!!!'
            print*,'on slip system # ',k
            STOP
          end if
        end do

c-------------------------------------------------------------------
c  Begin initialization loop over grains
c-------------------------------------------------------------------

        n = 0
        idum = 1
        do i = 1,num_grains

c-------------------------------------------------------------------
c  Initialize grain orientations in Euler angles
c-------------------------------------------------------------------

           psi1 =  Pi * ran0(idum)
           psi2 =  acos(ran0(idum))
           psi3 =  Pi * ran0(idum)

c-------------------------------------------------------------------
c  Save Euler Angles
c-------------------------------------------------------------------

           n = n + 3
           statev(n-2) = psi1
           statev(n-1) = psi2
           statev(n)   = psi3

c-------------------------------------------------------------------
c  Initialize sig, back stress, and gamma_dots to zero.
c-------------------------------------------------------------------

           do j = 1,12+num_slip_sys
             n = n + 1
             statev(n) = 0.0
           end do

c-------------------------------------------------------------------
c  Initialize iso shear stress for each grain's slip system
c-------------------------------------------------------------------

           do j = 1,num_slip_sys
             n = n + 1
             statev(n) = g_zero
           end do

c-------------------------------------------------------------------
c  Initialize back stress for each grain's slip system to zero.
c-------------------------------------------------------------------

           do j = 1,num_slip_sys
             n = n + 1
             statev(n) = 0
           end do

         end do ! num_grains
            
      end if  ! time<>0

c===================================================================
c  Calculate anisotropic compliance tensor components
c===================================================================

      det =  C11**3 + 2 * C12**3 - 3 * C11 * C12**2
      B11 = (C11 * C11 - C12 * C12) / det
      B12 = (C12 * C12 - C11 * C12) / det
      B44 =  0.25 / C44

c--------------------------------------------------------------------
c  Calculate L(i,j)
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          array1(i,j) = (dfgrd1(i,j) - dfgrd0(i,j)) / (dtime + Tiny)
          array2(i,j) = (dfgrd1(i,j) + dfgrd0(i,j)) / 2
        end do
      end do
      call inverse_3x3(array2,array3)
      call aa_dot_bb(3,array1,array3,xL_global)

c--------------------------------------------------------------------
c  Calculate D(i,j)
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          D(i,j) = 0.5 * (xL_global(i,j) + xL_global(j,i))
        end do
      end do

c--------------------------------------------------------------------
c  Calculate hydrostatic stress rate increase
c--------------------------------------------------------------------

      trace_D = D(1,1) + D(2,2) + D(3,3)
      sig_hydro_dot = (C11 + 2 * C12) * trace_D / 3

c--------------------------------------------------------------------
c  Calculate D_eff
c--------------------------------------------------------------------

      call aa_dot_dot_bb(3,D,D,sum)
      D_eff = sqrt(2*sum/3)

c--------------------------------------------------------------------
c  Initialize sig_avg(i,j) to zero.
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          sig_avg(i,j) = 0.0
        end do
      end do

c--------------------------------------------------------------------
c  Initialize ddsdde_4th(i,j,k,l) to zero.
c--------------------------------------------------------------------

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          ddsdde_4th(i,j,k,l) = 0.0
         end do
        end do
       end do
      end do

c-------------------------------------------------------------------
c  Initialize ISV counter
c-------------------------------------------------------------------

      Nisv  = 0	! for inputs  here
      Nisv2 = 0	! for outputs later

c====================================================================
c  BEGIN LOOP OVER GRAINS !!!
c====================================================================

      do m = 1,num_grains

c-------------------------------------------------------------------
c  Read in Euler angles defining current orientation of xtals.
c  Create direction cosine matrix
c-------------------------------------------------------------------

      do i = 1,3
        Nisv = Nisv + 1
        psi(i) = statev(Nisv)
      end do
      call euler_to_matrix(psi(1),psi(2),psi(3),dir_cos)
            
c-------------------------------------------------------------------
c  Read xtal iso stress.  Calc sig_hydro for later use
c-------------------------------------------------------------------

      do j = 1,3
        do i = 1,j
          Nisv = Nisv + 1
          sig(i,j) = statev(Nisv)
          sig(j,i) = sig(i,j)
        end do
      end do
      trace_sig     = sig(1,1) + sig(2,2) + sig(3,3)
      sig_hydro_ref = trace_sig / 3

c-------------------------------------------------------------------
c  Read xtal back stress tensor
c-------------------------------------------------------------------

      do j = 1,3
        do i = 1,j
          Nisv = Nisv + 1
          a(i,j) = statev(Nisv)
          a(j,i) = a(i,j)
        end do
      end do

c-------------------------------------------------------------------
c  Read in shear strain rates.
c-------------------------------------------------------------------

      do i = 1,num_slip_sys
        Nisv = Nisv + 1
        gamma_dot(i) = statev(Nisv)
      end do

c-------------------------------------------------------------------
c  Read iso shear values
c-------------------------------------------------------------------

      do i = 1,num_slip_sys
        Nisv = Nisv + 1
        g0(i) = statev(Nisv)
      end do

c-------------------------------------------------------------------
c  Read back stress values
c-------------------------------------------------------------------

      do i = 1,num_slip_sys
        Nisv = Nisv + 1
        x0(i) = statev(Nisv)
      end do

c--------------------------------------------------------------------
c  Rotate velocity gradient BACKWARDS to xtal orientation.
c--------------------------------------------------------------------

      call transpose(3,dir_cos,array2)
      call aa_dot_bb(3,array2,xL_global,array3)
      call aa_dot_bb(3,array3,dir_cos,xL_xtal)

c--------------------------------------------------------------------
c  Calculate symmetric and antisymmetric parts of velocity gradient
c  in local xtal coordinates
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          D_xtal(i,j) = 0.5 * (xL_xtal(i,j) + xL_xtal(j,i))
          W_xtal(i,j) = 0.5 * (xL_xtal(i,j) - xL_xtal(j,i))
        end do
      end do

c-------------------------------------------------------------------
c  Calculate tau_cr() at beginning of this step
c-------------------------------------------------------------------

      call flow_rule(num_slip_sys,       g0,  gamma_dot,
     &                 gamma_zero, flow_exp,      tau_c)

c--------------------------------------------------------------------
c  Calculate resolved shear stresses, tau()
c--------------------------------------------------------------------
 
      do k = 1,num_slip_sys
        call aa_dot_dot_bb(3,sig,p(1,1,k),tau(k))
      end do
 
c--------------------------------------------------------------------
c  Temporarily scale tau_c up.  This does NOT affect final answer.
c--------------------------------------------------------------------
 
      factor = 1.0
      do i = 1,num_slip_sys
        factor = max(factor,abs(tau(i)/tau_c(i)))
      end do
      factor = 1.000001 * factor
      do i = 1,num_slip_sys
        tau_c(i) = tau_c(i) * factor
      end do
 
c--------------------------------------------------------------------
c  Reset slip system bookkeeping stuff
c--------------------------------------------------------------------

      N_active_sys = 0
      do i = 1,num_slip_sys
        id(i) = 0
        icode(i) = 0
        gamma_dot(i) = 0
      end do

c--------------------------------------------------------------------
c  Reinitilize "merged" Schmid tensors
c--------------------------------------------------------------------

      do k = 1,num_slip_sys
        do i = 1,3
         do j = 1,3
           pm(i,j,k) = p(i,j,k)
           qm(i,j,k) = q(i,j,k)
         end do
        end do
      end do

c====================================================================
c  TOP OF !!! LOOP !!! TO DETERMINE ADDITIONAL ACTIVE SLIP SYSTEMS
c====================================================================

      Nctr = 0
      time_so_far = 0
  100 continue

c--------------------------------------------------------------------
c  Calculate D_p 
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          D_p(i,j) = 0
          do k = 1,N_active_sys
            D_p(i,j) = D_p(i,j) + gamma_dot(id(k)) * pm(i,j,id(k))
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Calculate D_el
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          D_el(i,j) = D_xtal(i,j) - D_p(i,j)
        end do
      end do

c--------------------------------------------------------------------
c  Calculate sig_dot & tau_dot
c--------------------------------------------------------------------

      call C_dot_dot_bb(C11,C12,C44,D_el,sig_dot)

      do k = 1,num_slip_sys
        if (icode(k) .eq. 0) then
          call aa_dot_dot_bb(3,sig_dot,p(1,1,k),tau_dot(k))
        end if
      end do

c--------------------------------------------------------------------
c  Determine next active slip system and time required
c  to get to it.
c--------------------------------------------------------------------

      i_min = 0
      t_min = dtime
      do i = 1,num_slip_sys
        if (icode(i) .eq. 0) then
          time0(i) = abs(tau_c(i) * sgn(tau_dot(i)) - tau(i))
     &               /abs(tau_dot(i) + Tiny)
          if (time0(i) .lt. t_min) then
            t_min = time0(i)
            i_min = i
          end if
        end if
      end do

c--------------------------------------------------------------------
c  Determine if out of time for current time step.
c--------------------------------------------------------------------

      out_of_time = .false.
      if (time_so_far + t_min .ge. dtime) then
        t_min = dtime - time_so_far
        out_of_time = .true.
      end if

c--------------------------------------------------------------------
c  Update stress
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          sig(i,j) = sig(i,j) + sig_dot(i,j) * t_min
        end do
      end do

c--------------------------------------------------------------------
c  Update tau(i)
c--------------------------------------------------------------------

      do i = 1,num_slip_sys
        if (icode(i).eq.0) then
          tau(i) = tau(i) + tau_dot(i) * t_min
        end if
      end do

c--------------------------------------------------------------------
c  Activate next active slip system if not out of time
c--------------------------------------------------------------------

      if (Out_of_time) go to 200

      N_active_sys = N_active_sys + 1
      id(N_active_sys) = i_min
      icode(i_min) = -1

c--------------------------------------------------------------------
c  Identify any "simultaneous" active slip systems
c--------------------------------------------------------------------

        do n = 1,num_slip_sys
          if (icode(n) .eq. 0) then
            test = abs((time0(n) - t_min)) / (dtime + Tiny)
     &             +(abs(tau_dot(n)) - abs(tau_dot(i_min))) / 
     &             abs(tau_dot(i_min)+Tiny)
            if (abs(test) .lt. 1.E-4) icode(n) = id(N_active_sys)
          end if
        end do

c--------------------------------------------------------------------
c  Add (merge) any necessary Schmid tensors
c--------------------------------------------------------------------

        do n = 1,num_slip_sys
          if (icode(n) .eq. id(N_active_sys)) then
            k = icode(n)
            ksign = sgn(tau_dot(n)) * sgn(tau_dot(k))
            do i = 1,3
              do j = 1,3
                pm(i,j,k) = pm(i,j,k) + p(i,j,n) * ksign
                qm(i,j,k) = qm(i,j,k) + q(i,j,n) * ksign
              end do
            end do
          end if
        end do

c--------------------------------------------------------------------
c  Calculate lower triangular matrix for simultaneous eqs
c--------------------------------------------------------------------

      k = N_active_sys
      call C_dot_dot_bb(C11,C12,C44,pm(1,1,id(k)),array1)
      do j = 1,N_active_sys
        call aa_dot_dot_bb(3,array1,pm(1,1,id(j)),PdCdP(k,j))
      end do
      call symmetric_decomp(5,PdCdP,N_active_sys,N_active_sys)

c--------------------------------------------------------------------
c  Calculate right hand side of system of equations
c--------------------------------------------------------------------

      k = N_active_sys
      call C_dot_dot_bb(C11,C12,C44,pm(1,1,id(k)),array1)
      call aa_dot_dot_bb(3,D_xtal,array1,PdCdD(k))
      do i = 1,N_active_sys
        vector(i) = PdCdD(i)
      end do

c--------------------------------------------------------------------
c  Do double back-substitution to solve for gamma_dot(i)
c--------------------------------------------------------------------

      call symmetric_backsub(5,N_active_sys,PdCdP,vector)

c--------------------------------------------------------------------
c  Store solution in gamma_dot array
c--------------------------------------------------------------------

      do i = 1,N_active_sys
        gamma_dot(id(i)) = vector(i)
      end do

      do i = 1,num_slip_sys
        if (icode(i) .gt. 0) then
          ksign = sgn(tau_dot(i)) * sgn(tau_dot(icode(i)))
          gamma_dot(i) = gamma_dot(icode(i)) * ksign
        end if
      end do

c--------------------------------------------------------------------
c  If not Out_of_time then loop back up to top
c--------------------------------------------------------------------

      time_so_far = time_so_far + t_min

      if (N_active_sys .lt. Max_active_sys) go to 100

c====================================================================
c  BEGIN ITERATION FOR IMPLICIT VERSION
c====================================================================

  200 continue
      Nctr = Nctr + 1

c--------------------------------------------------------------------
c  Calculate xL_p
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          xL_p(i,j) = 0
          do k = 1,N_active_sys
            xL_p(i,j) = xL_p(i,j)+gamma_dot(id(k))*(pm(i,j,id(k))
     &                + qm(i,j,id(k)))
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Calculate D_p and W_p
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          D_p(i,j) = 0.5 * (xL_p(i,j) + xL_p(j,i))
          W_p(i,j) = 0.5 * (xL_p(i,j) - xL_p(j,i))
        end do
      end do

c--------------------------------------------------------------------
c  Calculate Lattice Spin
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          W_el(i,j) = W_xtal(i,j) - W_p(i,j)
        end do
      end do

c--------------------------------------------------------------------
c  Convert spin tensor to a rotation matrix
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          array1(i,j) = W_el(i,j) * dtime
        end do
      end do

      call spin_to_matrix(array1,xtal_rot_matrix)

c--------------------------------------------------------------------
c  Rotate D_xtal(i,j) BACKWARDS by (W_el * dtime)
c--------------------------------------------------------------------

      call transpose(3,xtal_rot_matrix,array2)
      call aa_dot_bb(3,array2,D_xtal,array3)
      call aa_dot_bb(3,array3,xtal_rot_matrix,D)

c--------------------------------------------------------------------
c  Calculate right hand side of system of equations
c--------------------------------------------------------------------

      do k = 1,N_active_sys
        call C_dot_dot_bb(C11,C12,C44,pm(1,1,id(k)),array1)
        call aa_dot_dot_bb(3,array1,D,PdCdD(k))
        vector(k) = PdCdD(k)
      end do

c--------------------------------------------------------------------
c  Do double back-substitution to solve for gamma_dot(i)
c--------------------------------------------------------------------

      call symmetric_backsub(5,N_active_sys,PdCdP,vector)

c--------------------------------------------------------------------
c  Store previous iteration for convergence check later.
c--------------------------------------------------------------------

      do i = 1,num_slip_sys
        gamma_check(i) = gamma_dot(i)
      end do

c--------------------------------------------------------------------
c  Store solution in gamma_dot array
c--------------------------------------------------------------------

      do i = 1,N_active_sys
        gamma_dot(id(i)) = vector(i)
      end do

      do i = 1,num_slip_sys
        if (icode(i) .gt. 0) then
          ksign = sgn(tau_dot(i)) * sgn(tau_dot(icode(i)))
          gamma_dot(i) = gamma_dot(icode(i)) * ksign
        end if
      end do

c--------------------------------------------------------------------
c  Check for convergence
c--------------------------------------------------------------------

      Converged = .true.
      do i = 1,num_slip_sys
        test = (gamma_dot(i) - gamma_check(i)) / (D_eff + Tiny)
        if (abs(test) .gt. Tolerance) Converged = .false.
      end do

c--------------------------------------------------------------------
c  If not converged, go back and do it all over again
c--------------------------------------------------------------------

      if (.not.Converged .and. Nctr.lt.10) go to 200

c====================================================================
c  END OF ITERATIONS
c====================================================================
c--------------------------------------------------------------------
c  Update xtal orientation.
c--------------------------------------------------------------------

      call aa_dot_bb(3,dir_cos,xtal_rot_matrix,array1)
      call matrix_to_euler(array1,psi(1),psi(2),psi(3))
      call euler_to_matrix(psi(1),psi(2),psi(3),dir_cos)

c--------------------------------------------------------------------
c  Write out Euler angles in Roe convention
c--------------------------------------------------------------------

c      write(7,'(e10.2,2i5,3f12.2)')time(2),m,kinc,(psi(i)*180/Pi,i=1,3)

c--------------------------------------------------------------------
c  Calculate hardening law for isotropic part
c--------------------------------------------------------------------

      call Iso_Hardening_Law(num_slip_sys,   G_dir,   G_dyn,  dtime, 
     &                       gamma_dot,    xLatent,      g0,      g)

c--------------------------------------------------------------------
c  Calculate hardening law for back stress
c--------------------------------------------------------------------

      call Kinematic_Law(num_slip_sys,  X_dir,  X_dyn,  dtime, 
     &                      gamma_dot,     x0,      x)

c====================================================================
c  MOVE STRESS OUT TO YIELD SURFACE
c====================================================================
c--------------------------------------------------------------------
c  Calculate resolved shear stresses
c--------------------------------------------------------------------

      do i = 1,num_slip_sys
        call aa_dot_dot_bb(3,sig,p(1,1,i),tau(i))
      end do

c--------------------------------------------------------------------
c  Calculate critical resolved shear stresses
c--------------------------------------------------------------------

      call flow_rule(num_slip_sys,        g,  gamma_dot,
     &                 gamma_zero, flow_exp,      tau_c)

c--------------------------------------------------------------------
c  Calculate delta_tau values with appropriate sign
c--------------------------------------------------------------------

      do i = 1,N_active_sys
        delta_tau(i) = tau_c(id(i))*sgn(tau_dot(id(i)))-tau(id(i))
        do j = 1,num_slip_sys
          if (id(i) .eq. icode(j)) then
            ksign = sgn(tau_dot(id(i))) * sgn(tau_dot(j))
            delta_tau(i) = delta_tau(i) + ksign *
     &                       (tau_c(j) * sgn(tau_dot(j)) - tau(j))
          end if
        end do
      end do

c--------------------------------------------------------------------
c  Create lower triangular matrix for simultaneous eqs
c  Note that nothing happens if N_active_sys=0 !!!
c--------------------------------------------------------------------

      do i = 1,N_active_sys
        do j = 1,i
          call aa_dot_dot_bb(3,pm(1,1,id(i)),pm(1,1,id(j)),PdP(i,j))
        end do
      end do

c--------------------------------------------------------------------
c  Perform Cholesky decomposition
c--------------------------------------------------------------------

      call symmetric_decomp(5,PdP,1,N_active_sys)

c--------------------------------------------------------------------
c  Do double back-substitution to solve for lambda's
c--------------------------------------------------------------------

      call symmetric_backsub(5,N_active_sys,PdP,delta_tau)

c--------------------------------------------------------------------
c  Update stress
c--------------------------------------------------------------------

      do k = 1,N_active_sys
        do i = 1,3
          do j = 1,3
            sig(i,j) = sig(i,j) + delta_tau(k) * pm(i,j,id(k))
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Add in hydrostatic part
c--------------------------------------------------------------------

      sig_tr = sig(1,1) + sig(2,2) + sig(3,3)
      do i = 1,3
        sig(i,i) = sig(i,i)-sig_tr/3+sig_hydro_ref+sig_hydro_dot*dtime
      end do

c--------------------------------------------------------------------
c  Calculate delta_tau values for slip system back stress
c--------------------------------------------------------------------

      do i = 1,N_active_sys
        delta_tau(i) = x(id(i)) - x0(id(i))
        do j = 1,num_slip_sys
          if (id(i) .eq. icode(j)) then
            ksign = sgn(gamma_dot(id(i))) * sgn(gamma_dot(j))
            delta_tau(i) = delta_tau(i) + ksign * (x(j) - x0(j))
          end if
        end do
      end do

c--------------------------------------------------------------------
c  Do double back-substitution to solve for lambda's
c--------------------------------------------------------------------

      call symmetric_backsub(5,N_active_sys,PdP,delta_tau)

c--------------------------------------------------------------------
c  Update back stress
c--------------------------------------------------------------------

      do k = 1,N_active_sys
        do i = 1,3
          do j = 1,3
            a(i,j) = a(i,j) + delta_tau(k) * pm(i,j,id(k))
          end do
        end do
      end do

c--------------------------------------------------------------------
c  Add iso and back stresses together
c--------------------------------------------------------------------

        do i = 1,3
          do j = 1,3
            array1(i,j) = sig(i,j) + a(i,j)
          end do
        end do

c--------------------------------------------------------------------
c  Rotate stress tensor forward to current orientation.
c--------------------------------------------------------------------

      call transpose(3,dir_cos,array2)
      call aa_dot_bb(3,dir_cos,array1,array3)
      call aa_dot_bb(3,array3,array2,array1)

c--------------------------------------------------------------------
c  Add xtal stress to sig_avg
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          sig_avg(i,j) = sig_avg(i,j) + array1(i,j)
        end do
      end do

c====================================================================
c  BEGIN CALCULATION OF THE JACOBIAN (TANGENT STIFFNESS MATRIX).
c====================================================================
c--------------------------------------------------------------------
c  Calculate the derivative of the plastic part of the rate of
c  deformation tensor in the current configuration wrt sigma.
c--------------------------------------------------------------------

      call Dp_wrt_sig(num_slip_sys,  dtime,  flow_exp,    p,     g,
     &                gamma_dot,   ddpdsig)

c--------------------------------------------------------------------
c  Need to take the inverse of Array4.  Since it relates two 2nd
c  rank tensors that are both symmetric, Array4 can be xformed to 
c  Voigt notation style in order to do the inverse, then xformed back.
c--------------------------------------------------------------------

      call forth_to_Voigt (ddpdsig,Array4)

c--------------------------------------------------------------------
c  Add compliance tensor to ddpdsig(i,j,k,l)
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          if (i.eq.j) then
            array4(i,j) = array4(i,j) + B11
          else
            array4(i,j) = array4(i,j) + B12
          end if
        end do
        k = i + 3
        array4(k,k) = array4(k,k) + B44
      end do

c--------------------------------------------------------------------
c  Take inverse of matrix and put back in 4th rank tensor
c--------------------------------------------------------------------

      call symmetric_inverse(6,6,Array4,Array5)
      call Voigt_to_forth   (Array5,Array6)

c--------------------------------------------------------------------
c  Rotate Jacobian from xtal frame to global orientation
c--------------------------------------------------------------------

      call rotate_4th(dir_cos,Array6,Array7)

c--------------------------------------------------------------------
c  Add xtal Jacobian to total Jacobian
c--------------------------------------------------------------------

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          ddsdde_4th(i,j,k,l) = ddsdde_4th(i,j,k,l) + Array7(i,j,k,l)
         end do ! l
        end do ! k
       end do ! j
      end do ! i

c===================================================================
c  Store the internal variables in the statev() array
c===================================================================
c-------------------------------------------------------------------
c  Store xtal orienations
c-------------------------------------------------------------------

      do i = 1,3
        Nisv2 = Nisv2 + 1
        statev(Nisv2) = psi(i)
      end do
            
c-------------------------------------------------------------------
c  Store each xtal iso stress
c-------------------------------------------------------------------

      do j = 1,3
        do i = 1,j
          Nisv2 = Nisv2 + 1
          statev(Nisv2) = sig(i,j)
        end do
      end do

c-------------------------------------------------------------------
c  Store each xtal back stress
c-------------------------------------------------------------------

      do j = 1,3
        do i = 1,j
          Nisv2 = Nisv2 + 1
          statev(Nisv2) = a(i,j)
        end do
      end do

c-------------------------------------------------------------------
c  Store shear strain rates.
c-------------------------------------------------------------------

      do i = 1,num_slip_sys
       Nisv2 = Nisv2 + 1
       statev(Nisv2) = gamma_dot(i)
      end do

c-------------------------------------------------------------------
c  Store the reference shear stresses.
c-------------------------------------------------------------------

      do i = 1,num_slip_sys
        Nisv2 = Nisv2 + 1
        statev(Nisv2) = g(i)
      end do

c-------------------------------------------------------------------
c  Store the slip system back stress
c-------------------------------------------------------------------

      do i = 1,num_slip_sys
        Nisv2 = Nisv2 + 1
        statev(Nisv2) = x(i)
      end do

c===================================================================
c  End loop over all the grains
c===================================================================

      end do ! m over all grains

c--------------------------------------------------------------------
c  Divide by number of grains to get avg
c--------------------------------------------------------------------

      do i = 1,3
        do j = 1,3
          sig_avg(i,j) = sig_avg(i,j) / num_grains
        end do
      end do

c--------------------------------------------------------------------
c  Average ddsdde_4th(i,j,k,l)
c--------------------------------------------------------------------

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          ddsdde_4th(i,j,k,l) = ddsdde_4th(i,j,k,l) / num_grains
         end do ! l
        end do ! k
       end do ! j
      end do ! i

c--------------------------------------------------------------------
c  Store the stress tensor in the ABAQUS stress 'vector'
c--------------------------------------------------------------------

      do i = 1,ndi
         stress(i) = sig_avg(i,i)
      end do
      if (nshr .eq. 1) stress(ndi+1) = sig_avg(1,2)
      if (nshr .eq. 3) then
         stress(4) = sig_avg(1,2)
         stress(5) = sig_avg(1,3)
         stress(6) = sig_avg(2,3)
      end if

c--------------------------------------------------------------------
c  Store the Jocbian in Voigt notation form.
c  One can try tinkering with the commented
c  line to improve convergence.
c--------------------------------------------------------------------

      do i = 1,3
       do j = i,3   ! not 1 !!!
        ia = i
        if (i.ne.j) ia=i+j+1
        do k = 1,3
         do l = k,3 ! not 1 !!!
          ib = k
          if (k.ne.l) ib=k+l+1
          array4(ia,ib) = ddsdde_4th(i,j,k,l)
c           if (ib.ge.4) array4(ia,ib) = 2 * array4(ia,ib)
         end do
        end do
       end do
      end do

c--------------------------------------------------------------------
c  1-D analysis
c--------------------------------------------------------------------

      if (ndi .eq. 1) then		! 1-D
         ddsdde(1,1) = array4(1,1)

c--------------------------------------------------------------------
c  2-D or axisymmetric analysis
c--------------------------------------------------------------------

      else if (ndi .eq. 2) then		! 2-D plane stress & axi
         do i = 1,2
            do j = 1,2
               ddsdde(i,j) = array4(i,j)
            end do
         end do
         ddsdde(1,3) = array4(1,4)
         ddsdde(2,3) = array4(2,4)
         ddsdde(3,1) = array4(4,1)
         ddsdde(3,2) = array4(4,2)
         ddsdde(3,3) = array4(4,4)

c--------------------------------------------------------------------
c  2-D plain strain analysis
c--------------------------------------------------------------------

      else if (ndi .eq. 3 .and. nshr .eq. 1) then ! plane strain
         do i = 1,4
            do j = 1,4
               ddsdde(i,j) = array4(i,j)
            end do
         end do

c--------------------------------------------------------------------
c  Full 3-D analysis
c--------------------------------------------------------------------

      else				! Full 3-D
         do i = 1,6
            do j = 1,6
               ddsdde(i,j) = array4(i,j)
            end do
         end do
      end if

c-------------------------------------------------------------------
c  THE END !!!
c-------------------------------------------------------------------

      return
      end

c====================================================================
c====================================================================
c====================== S U B R O U T I N E S =======================
c====================================================================
c====================================================================
c
c  Normalize the length of a vector to one.
c
c--------------------------------------------------------------------

      subroutine normalize_vector(x,y,z)

      implicit double precision (a-h,o-z)

      xlength = sqrt(x*x+y*y+z*z)
      if (xlength.eq.0) xlength = 1.D-300
      x = x / xlength
      y = y / xlength
      z = z / xlength

      return
      end

c====================================================================
c====================================================================
c
c  Transpose an ( n x n ) tensor.
c
c--------------------------------------------------------------------

      subroutine transpose(n,a,b)

      implicit double precision (a-h,o-z)

      dimension a(n,n), b(n,n)

      do i = 1,n
         do j = 1,n
            b(i,j) = a(j,i)
         end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the dot product of two 2nd rank tensors.
c  Result is stored in cc(i,j)
c
c--------------------------------------------------------------------

      subroutine aa_dot_bb(n,a,b,c)

      implicit double precision (a-h,o-z)

      dimension a(n,n), b(n,n), c(n,n)

      do i = 1,n
         do j = 1,n
            c(i,j) = 0
            do k = 1,n
               c(i,j) = c(i,j) + a(i,k) * b(k,j)
            end do
         end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the double dot product of two 2nd rank tensors.
c
c--------------------------------------------------------------------

      subroutine aa_dot_dot_bb(n,a,b,sum)

      implicit double precision (a-h,o-z)

      dimension a(n,n), b(n,n)

      sum = 0.0
      do i = 1,n
         do j = 1,n
            sum = sum + a(i,j) * b(i,j)
         end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the double dot product of the 4th rank elasticity
c  tensor and a symmetric 2nd rank tensor.  The elasticity
c  tensor MUST, MUST, MUST be in its principle orientation!!!!
c
c--------------------------------------------------------------------

      subroutine C_dot_dot_bb(C11,C12,C44,b,c)

      implicit double precision (a-h,o-z)

      dimension b(3,3), c(3,3)

      C11mC12 = C11 - C12
      twoC44  = 2 * C44
      trace = b(1,1) + b(2,2) + b(3,3)

      do i = 1,3
        c(i,i) = C11mC12 * b(i,i) + C12 * trace
      end do

      c(1,2) = twoC44 * b(1,2)
      c(1,3) = twoC44 * b(1,3)
      c(2,3) = twoC44 * b(2,3)

      c(2,1) = c(1,2)
      c(3,1) = c(1,3)
      c(3,2) = c(2,3)
      
      return
      end

c====================================================================
c====================================================================
c
c  Rotates any 3x3x3x3 tensor by a rotation matrix.
c
c  c(i,j,k,l) = a(i,m) * a(j,n) * a(k,p) * a(l,q) * b(m,n,p,q)
c
c--------------------------------------------------------------------

      subroutine rotate_4th(a,b,c)

      implicit double precision (a-h,o-z)

      dimension a(3,3), b(3,3,3,3), c(3,3,3,3), d(3,3,3,3)

      do m = 1,3
       do n = 1,3
        do k = 1,3
         do l = 1,3
          d(m,n,k,l) = a(k,1) * (a(l,1) * b(m,n,1,1) + 
     &		a(l,2) * b(m,n,1,2) + a(l,3) * b(m,n,1,3)) +
     &		a(k,2) * (a(l,1) * b(m,n,2,1) + 
     &		a(l,2) * b(m,n,2,2) + a(l,3) * b(m,n,2,3)) +
     &		a(k,3) * (a(l,1) * b(m,n,3,1) + 
     &		a(l,2) * b(m,n,3,2) + a(l,3) * b(m,n,3,3))
         end do
        end do
       end do
      end do

      do i = 1,3
       do j = 1,3
        do k = 1,3
         do l = 1,3
          c(i,j,k,l) = a(i,1) * (a(j,1) * d(1,1,k,l) + 
     &		a(j,2) * d(1,2,k,l) + a(j,3) * d(1,3,k,l)) +
     &		a(i,2) * (a(j,1) * d(2,1,k,l) + 
     &		a(j,2) * d(2,2,k,l) + a(j,3) * d(2,3,k,l)) +
     &		a(i,3) * (a(j,1) * d(3,1,k,l) + 
     &		a(j,2) * d(3,2,k,l) + a(j,3) * d(3,3,k,l))
         end do
        end do
       end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the inverse of a 3 x 3 matrix.
c
c--------------------------------------------------------------------

      subroutine inverse_3x3(a,b)

      implicit double precision (a-h,o-z)

      dimension a(3,3), b(3,3)

      b(1,1) = a(2,2) * a(3,3) - a(3,2) * a(2,3)
      b(1,2) = a(3,2) * a(1,3) - a(1,2) * a(3,3)
      b(1,3) = a(1,2) * a(2,3) - a(2,2) * a(1,3)
      b(2,1) = a(3,1) * a(2,3) - a(2,1) * a(3,3)
      b(2,2) = a(1,1) * a(3,3) - a(3,1) * a(1,3)
      b(2,3) = a(2,1) * a(1,3) - a(1,1) * a(2,3)
      b(3,1) = a(2,1) * a(3,2) - a(3,1) * a(2,2)
      b(3,2) = a(3,1) * a(1,2) - a(1,1) * a(3,2)
      b(3,3) = a(1,1) * a(2,2) - a(2,1) * a(1,2)

      det = a(1,1) * b(1,1) + a(1,2) * b(2,1) + a(1,3) * b(3,1)

      do i = 1,3
         do j = 1,3
            b(i,j) = b(i,j) / det
         end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate the inverse of a matrix using 
c  LU decomposition (Crout's method)
c
c  Reference: "Numerical Recipes" Section 2.3  p. 31
c
c--------------------------------------------------------------------

      subroutine symmetric_inverse(NP,n,a,b)

      implicit double precision (a-h,o-z)

      dimension a(NP,NP), b(NP,NP), c(NP,NP)

      do i = 1,n
         do j = 1,n
            c(i,j) = a(i,j)
            b(i,j) = 0.0
         end do
         b(i,i) = 1.0
      end do

      call symmetric_decomp(NP,c,1,n)
      do j = 1,n
         call symmetric_backsub(NP,n,c,b(1,j))
      end do

      return
      end

c====================================================================
c====================================================================
c           T
c  Do an L.L  (Cholesky) decomposition on a symmetric matrix.
c  Only the LOWER TRIANGULAR part of the stiffness matrix is used.
c  istart & istop determine which rows are decomposed.
c
c--------------------------------------------------------------------

      subroutine symmetric_decomp(NP,a,istart,istop)

      implicit double precision (a-h,o-z)

      dimension a(NP,NP)

      do k = istart,istop
        do i = 1,k-1
          sum = 0
          do j = 1,i-1
            sum = sum + a(i,j) * a(k,j)
          end do
          a(k,i) = (a(k,i) - sum) / a(i,i)
        end do
        sum = 0
        do j = 1,k-1
          sum = sum + a(k,j) * a(k,j)
        end do
        a(k,k) = sqrt(a(k,k) - sum)
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Do the back substitution to solve for 'x'.  The original
c  coefficient matrix must have been symmetric.
c  Only the LOWER TRIANGULAR part of the stiffness matrix is used.
c
c--------------------------------------------------------------------

      subroutine symmetric_backsub(NP,n,a,b)

      implicit double precision (a-h,o-z)

      dimension a(NP,NP), b(NP)

      do k = 1,n
        sum = 0
        do i = 1,k-1
          sum = sum + a(k,i) * b(i)
        end do
        b(k) = (b(k) - sum) / a(k,k)
      end do

      do k = n,1,-1
        sum = 0
        do i = k+1,n
          sum = sum + a(i,k) * b(i)
        end do
        b(k) = (b(k) - sum) / a(k,k)
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Restore a symmetric 4th rank tensor stored in Voigt notation 
c  back to its 4th rank form.
c
c--------------------------------------------------------------------

      subroutine Voigt_to_forth(b,a)

      implicit double precision (a-h,o-z)

      dimension a(3,3,3,3), b(6,6)

      do i = 1,3
       do j = 1,3
        ia = i
        if (i.ne.j) ia=9-i-j
        do k = 1,3
         do l = 1,3
          ib = k
          if (k.ne.l) ib=9-k-l
          a(i,j,k,l) = b(ia,ib)
          if (ia.gt.3) a(i,j,k,l) = a(i,j,k,l) / 2
          if (ib.gt.3) a(i,j,k,l) = a(i,j,k,l) / 2
         end do
        end do
       end do
      end do

      return
      end


c====================================================================
c====================================================================
c
c  Store a SYMMETRIC 4th rank tensor in Voigt notation.
c
c--------------------------------------------------------------------

      subroutine forth_to_Voigt(a,b)

      implicit double precision (a-h,o-z)

      dimension a(3,3,3,3), b(6,6)

      do i = 1,3
       do j = i,3   ! not 1 !!!
        ia = i
        if (i.ne.j) ia=9-i-j
        do k = 1,3
         do l = k,3 ! not 1 !!!
          ib = k
          if (k.ne.l) ib=9-k-l
          b(ia,ib) = a(i,j,k,l)
         end do
        end do
       end do
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Return the sign of a number.
c
c--------------------------------------------------------------------

      function sgn(a)

      implicit double precision (a-h,o-z)

      sgn = 1.0
      if (a .lt. 0.0) sgn = -1.0

      return
      end 

c====================================================================
c====================================================================
c
c  Calculate the determinant of a 3 x 3 matrix.
c
c--------------------------------------------------------------------

      function determinant(a)

      implicit double precision (a-h,o-z)

      dimension a(3,3)

      b1 = a(2,2) * a(3,3) - a(3,2) * a(2,3)
      b2 = a(3,1) * a(2,3) - a(2,1) * a(3,3)
      b3 = a(2,1) * a(3,2) - a(3,1) * a(2,2)

      determinant = a(1,1) * b1 + a(1,2) * b2 + a(1,3) * b3

      return
      end

c===================================================================
c===================================================================
c
c  Print out an array.
c
c-------------------------------------------------------------------

      subroutine print_array(n,a)

      implicit double precision (a-h,o-z)
      dimension a(n,n)

      do i = 1,n
         write(6,'(10f12.5)')(a(i,j),j=1,n)
      end do
      print*,' '

      return
      end

c===================================================================
c===================================================================
c
c  Print out Euler angles in Roe notation.
c
c-------------------------------------------------------------------

      subroutine matrix_to_euler(a,psi,theta,phi)

      implicit double precision (a-h,o-z)

      dimension a(3,3)

      pi = 4 * atan(1.0d0)

      if (abs(a(3,3)) .gt. 0.99999) then
        psi   = 0.0
        theta = 0.0
        phi   = atan2(a(2,1),a(1,1))
      else
        psi   = atan2(a(2,3),a(1,3))
        theta = acos(a(3,3))
        phi   = atan2(a(3,2),-a(3,1))
      end if

      return
      end

c====================================================================
c====================================================================
c
c  Convert Euler angles to a rotation matrix.  Euler angles are
c  in Roe convention.
c  1) psi   around z-axis
c  2) theta around y'-axis
c  3) phi   around z"-axis
c
c--------------------------------------------------------------------

      subroutine euler_to_matrix(psi,theta,phi,a)

      implicit double precision (a-h,o-z)

      dimension a(3,3)

      s1 = sin(psi)
      c1 = cos(psi)
      s2 = sin(theta)
      c2 = cos(theta)
      s3 = sin(phi)
      c3 = cos(phi)
            
      a(1,1) = c1*c2*c3-s1*s3
      a(2,1) = c3*c2*s1+s3*c1
      a(3,1) = -c3*s2
      a(1,2) = -s3*c2*c1-c3*s1
      a(2,2) = -s3*c2*s1+c3*c1
      a(3,2) = s3*s2
      a(1,3) = s2*c1
      a(2,3) = s2*s1
      a(3,3) = c2

      return
      end

c====================================================================
c====================================================================
c
c  Generate uniform random numbers between 0 and 1.
c  This is a fairly coarse generator in that it returns only
c  714,025 different values.
c
c  ***IMPORTANT***
c  Set idum equal to any positive number between 1 and 714,024
c  to initialize or reinitialize the sequence!!!
c
c  call with:   x = ran0(idum)   for example.
c
c  Reference: "Numerical Recipes" Section 7.1  p. 195
c
c--------------------------------------------------------------------

      function ran0(idum)

      implicit double precision (a-h,o-z)

      Parameter (M=714025, IA=1366, IC=150889, RM=1./M)

      if (idum.lt.1 .or. idum.ge.M) idum = 1
      idum = mod(ia*idum+ic,m)
      ran0 = float(idum)/float(m)

      return
      end

c====================================================================
c====================================================================
c
c  Convert a spin tensor to a rotation matrix 
c
c--------------------------------------------------------------------

      subroutine spin_to_matrix(a,b)

      implicit double precision (a-h,o-z)

      dimension a(3,3), b(3,3)

      Pi = 4*atan(1.D0)

c--------------------------------------------------------------------
c  Store spin tensor, a(i,j), as a rotation vector
c--------------------------------------------------------------------

      p1 = a(3,2)
      p2 = a(1,3)
      p3 = a(2,1)
      ang = sqrt(p1*p1+p2*p2+p3*p3)

      s = sin(ang)
      c = cos(ang)

c--------------------------------------------------------------------
c  Normalize vector
c--------------------------------------------------------------------

      if (ang .le. 1.D-300) then
        p1 = 0
        p2 = 0
        p3 = 1
      else
        p1 = p1 / ang
        p2 = p2 / ang
        p3 = p3 / ang
      end if

c--------------------------------------------------------------------
c  Calculate rotation matrix
c--------------------------------------------------------------------

      b(1,1) = c + (1 - c) * p1 * p1
      b(1,2) = (1 - c) * p1 * p2 - s * p3
      b(1,3) = (1 - c) * p1 * p3 + s * p2
      b(2,1) = (1 - c) * p2 * p1 + s * p3
      b(2,2) = c + (1 - c) * p2 * p2
      b(2,3) = (1 - c) * p2 * p3 - s * p1
      b(3,1) = (1 - c) * p3 * p1 - s * p2
      b(3,2) = (1 - c) * p3 * p2 + s * p1
      b(3,3) = c + (1 - c) * p3 * p3

      return
      end

c====================================================================
c====================================================================
c
c  This subroutine calculates the critical tau values, tau_c().
c  Note that it is strictly an isotropic value, even if
c  kinematic hardening is present.
c  Also, it is CRITICAL that the flow rule used IN THIS UMAT
c  be a function of max(abs(gamma_dot(i=1..N))) only.
c  If you do otherwise, then you will still get results,
c  but they will NOT be correct, and they will NOT even be
c  closer to some expected classical value.  They will just be wrong
c  because the wrong slip systems will become active.
c  This UMAT could be modified to accept tau_critical
c  as a function of each individual gamma_dot, however it
c  would add a huge amount of complexity, and accomplish very,
c  very, very little because this entire UMAT is limited to 
c  materials with little strain rate sensitivity (flow_exp<0.1),
c  so any gamma_dot raised to a small flow_exp will be very
c  close to 1 anyway.
c
c  tau_critical = g * abs(gamma_dot_max / gamma_zero) ^ flow_exp
c
c--------------------------------------------------------------------

      subroutine flow_rule(N,g,gamma_dot,gamma_zero,flow_exp,tau_c)

      implicit double precision (a-h,o-z)

      dimension g(N), gamma_dot(N), tau_c(N)

      gamma_max = 0.0
      do i = 1,N
        gamma_max = max(gamma_max,abs(gamma_dot(i)))
      end do

      do i = 1,N
        if (gamma_max .eq. 0.0) then
          tau_c(i) = g(i)
        else
          tau_c(i) = g(i) * (gamma_max / gamma_zero) ** flow_exp
        end if
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate isotropic hardening law at end of time step.
c  It's a direct hardening, dynamic recovery law in implicit form.
c  Latent hardening is in the direct hard & dyn recovery terms.
c
c  g_dot = G_dir * SUM(h_matrix * abs(gamma_dot)) -
c          G_dyn * g * SUM(h_matrix * abs(gamma_dot))
c
c--------------------------------------------------------------------

      subroutine Iso_Hardening_Law(N,  G_dir,   G_dyn,  dtime, 
     &                             gamma_dot, xLatent,  g0, g)

      implicit double precision (a-h,o-z)

      dimension g0(N), g(N), gamma_dot(N)

      sum = 0
      do i = 1,N
        sum = sum + abs(gamma_dot(i))
      end do
      do k = 1,N
        sum00 = xLatent * sum - (xLatent - 1) * abs(gamma_dot(k))
        g(k) = (g0(k) + G_dir*sum00*dtime) / (1 + G_dyn*sum00*dtime)
      end do

      return
      end

c====================================================================
c====================================================================
c
c  Calculate kinematic hardening law at end of time step.
c  It's a direct hardening, dynamic recovery law in implicit form.
c
c  x_dot = X_dir * gamma_dot - X_dyn * x * abs(gamma_dot)
c
c--------------------------------------------------------------------

      subroutine Kinematic_Law(N,  X_dir,   X_dyn,  dtime, 
     &                         gamma_dot,      x0,      x)

      implicit double precision (a-h,o-z)

      dimension x0(N), x(N), gamma_dot(N)

      do k = 1,N
        x(k) = (x0(k) + X_dir * dtime * gamma_dot(k))
     &			 / (1 + X_dyn * dtime * abs(gamma_dot(k)))
      end do

      return
      end

c====================================================================
c====================================================================
c  Calculate the derivative of the plastic part of the rate of
c  deformation tensor in the current configuration wrt sigma.
c  This assumes that a power law form of flow rule is defined above.
c--------------------------------------------------------------------

      subroutine Dp_wrt_sig(N,  dtime,  flow_exp,    p,     g,
     &                      gamma_dot,   ddpdsig)

      implicit double precision (a-h,o-z)

      dimension gamma_dot(N), g(N), p(3,3,N), ddpdsig(3,3,3,3)
      Tiny = 1.D-30

      do i = 1,3
       do j = i,3 ! not 1
        do k = 1,3
         do l = k,3 ! not 1
           ddpdsig(i,j,k,l) = 0.0
           do m = 1,N
             ddpdsig(i,j,k,l) = ddpdsig(i,j,k,l)+abs(gamma_dot(m))
     &                        * p(i,j,m) * p(k,l,m) / g(m)
           end do ! num_slip_sys
           ddpdsig(i,j,k,l) = ddpdsig(i,j,k,l) * dtime / 
     &                        (flow_exp + Tiny)
         end do ! l
        end do ! k
       end do ! j
      end do ! i

      return
      end

