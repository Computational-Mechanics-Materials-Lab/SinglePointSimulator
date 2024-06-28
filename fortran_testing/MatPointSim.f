c-------------------------------------------------------------------
c  Abaqus Simulator for Tension, Compression, Torsion, 
c  Plane-Strain Tens/Comp, and Simultaneous Tens/Comp w/ Torsion
c-------------------------------------------------------------------

c       include 'umats/kitchen_sink.umat'
c       include 'umats/wrapper.umat'
c       include 'umats/fast.umat'
c       include 'umats/intergrain.umat'
c       include 'umats/sarma.umat'
c       include 'umats/sherby.umat'
c       include 'umats/bamfast.umat'
c       include 'umats/bammacrosurf.umat'
c       include 'cuitino_ortiz.umat'
c       include 'umats/umat.porous.MARIN-HORST.june00.txt'
      include 'Fast_Xtal_Plasticity_UMAT'

      implicit double precision (a-h,o-z)

c-------------------------------------------------------------------

      parameter(    Max_isv   = 20000,
     &		        Max_props = 100,
     &		        Tolerance = .001,
     &		        Max_itr   = 4,
     &		        ntens	  = 6,
     & 		        ndi       = 3,
     &		        nshr	  = 3,
     &		        noel	  = 1,
     &		        Num_modes = 5)

      character*8   cmname
      character*132 text
      character*80  loading_mode
      character*80  loading_key(2,Num_modes)

      logical Out_of_time

c-------------------------------------------------------------------
c  Dimension arrays used only in this main program.
c-------------------------------------------------------------------

      dimension
     &	array1	(3,3),	    ! Dummy Array
     &	array2	(3,3),	    ! Dummy Array
     &	array3	(3,3),	    ! Dummy Array
     &	D_dt	(3,3),	    ! Rate of Def tensor * dtime
     &	W_dt	(3,3),	    ! Spin tensor * dtime
     &	sig	(3,3),		    ! Stress Tensor
     &	statev_ref(Max_isv)	! Reference ISVs

c-------------------------------------------------------------------
c  Dimension arrays passed into the UMAT sub
c-------------------------------------------------------------------

      dimension
     &	coords(3),	            ! Coordinates of Gauss pt. being evaluated
     &	ddsdde(ntens,ntens),    ! Tangent Stiffness Matrix
     &	ddsddt(ntens),	        ! Change in stress per change in temperature
     &	dfgrd0(3,3),	        ! Deformation gradient at beginning of step
     &	dfgrd1(3,3),	        ! Deformation gradient at end of step
     &	dpred(1),	            ! Change in predefined state variables
     &	drplde(ntens),	        ! Change in heat generation per change in strain
     &	drot(3,3),	            ! Rotation matrix
     &	dstrain(ntens),	        ! Strain increment tensor stored in vector form
     &	predef(1),	            ! Predefined state vars dependent on field variables
     &	props(Max_props),	    ! Material properties passed in
     &	statev(Max_isv),	    ! State Variables
     &	strain(ntens),	        ! Strain tensor stored in vector form
     &	stress(ntens),	        ! Cauchy stress tensor stored in vector form
     &	time(2)		            ! Step Time and Total Time

c-------------------------------------------------------------------
c  Initialize loading modes
c-------------------------------------------------------------------

      data
     &	loading_key(1,1),loading_key(2,1) / 'C', 'TensionComp' /
     &	loading_key(1,2),loading_key(2,2) / 'T', 'Torsion'     /
     &	loading_key(1,3),loading_key(2,3) / 'P', 'PlaneStrain' /
     &	loading_key(1,4),loading_key(2,4) / 'S', 'SimuComTors' /
     &	loading_key(1,5),loading_key(2,5) / 'B', 'BiaxialTens' /

c-------------------------------------------------------------------
c  Initialize Pi
c-------------------------------------------------------------------

      Pi = 4. * atan(1.d0)

c-------------------------------------------------------------------
c  Initialize time & deformation gradients
c-------------------------------------------------------------------

      time(1)	= 0.0
      time(2)	= 0.0
      E_eff     = 0.0

      do i = 1,6
        strain(i)  = 0.0
        dstrain(i) = 0.0
      end do

      do i = 1,3
       do j = 1,3
         dfgrd0(i,j) = 0.0
         dfgrd1(i,j) = 0.0
         drot  (i,j) = 0.0
       end do
       dfgrd0(i,i) = 1.0
       dfgrd1(i,i) = 1.0
       drot  (i,i) = 1.0
      end do

c-------------------------------------------------------------------
c  Open ABAQUS files
c-------------------------------------------------------------------

      open(1,file='sim.inp',status='old')
      open(7,file='sim.msg',status='unknown')
      open(8,file='sim.dat',status='unknown')

c-------------------------------------------------------------------
c  Print header in .dat file
c-------------------------------------------------------------------

      write(8,'(42a,38a)')'     TIME     E_eff    SIG(1,1)  SIG(2,2)',
     & '  SIG(3,3)  SIG(1,2)  SIG(1,3) SIG(2,3)'

c-------------------------------------------------------------------
c  Read in number of ISVs
c-------------------------------------------------------------------

  100 read(1,'(a132)') text
      n = index(text,'DEPVAR')
      if (n .eq. 0) go to 100

      read(1,*) nstatv
      if (nstatv .gt. Max_isv) then
        print*,'STOP!!!  The number of state variables in the input'
        print*,'         file exceeds the max value in the Max_isv'
        print*,'         variable in the program.'
        Print*,' '
        Print*,'Set Max_isv in program to at least: ',nstatv+1,'.'
        Print*,' '
        STOP
      end if

c-------------------------------------------------------------------
c  Read in number of material properties
c-------------------------------------------------------------------

  200 read(1,'(a132)') text
      n = index(text,'CONSTANTS')
      if (n .eq. 0) go to 200

      read(text(n+10:n+20),*) nprops

      print*,' '

      if (nprops .gt. Max_props) then
        print*,'STOP!!!  The number of mat properties in the input'
        print*,'         file exceeds the max value in the Max_props'
        print*,'         variable in the program.'
        Print*,' '
        Print*,'Set Max_props in program to at least: ',nprops+1,'.'
        Print*,' '
        STOP
      end if

c-------------------------------------------------------------------
c  Read in material properties
c-------------------------------------------------------------------

      do i = 1,nprops/8
        read(1,'(a132)') text
        read(text,*)(props(8*(i-1)+j),j=1,8)
      end do

      itest = nprops - nprops / 8 * 8
      print*,'itest = ',itest
      if (itest .gt. 0) then
        read(1,'(a132)') text
        read(text,*)(props(j),j=(nprops/8)*8+1,nprops)
      end if

c===================================================================
c  Start new time step
c  Read in data
c===================================================================

  300 time(1) = 0.0
      Out_of_time = .false.

      read(1,'(a132)',END=999) text
      n = index(text,'STEP')
      if (n .eq. 0) go to 300

      read(1,'(a132)') text
      n = index(text,'**')
      if (n.eq.1) go to 300

      read(text,*)dtime,time_max,dtime_max
      print*,' '
      print*,' dtime    time_max    dtime_max'
      write(6,'(f7.4,f12.2,f11.2/)')dtime,time_max,dtime_max

c-------------------------------------------------------------------
c  Read in loading mode
c-------------------------------------------------------------------

  350 loading_mode = 'none'
      read(1,'(a132)') text
      do i = 1,Num_modes
        if (text(1:1) .eq. loading_key(1,i)) then
          loading_mode = loading_key(2,i)
        end if
      end do
      if (loading_mode .eq. 'none') go to 350

c-------------------------------------------------------------------
c  Read in loading direction
c-------------------------------------------------------------------

      read(text(2:2),*)i
      read(text(3:3),*)j

c-------------------------------------------------------------------
c  Read in loading rate
c-------------------------------------------------------------------

      if (loading_mode .eq. 'SimuComTors') then
        read(text(5:132),*)velocity2,velocity
      else if (loading_mode .eq. 'BiaxialTens') then
        read(text(5:132),*)velocity,velocity2
      else
        read(text(5:132),*)velocity
      end if



c===================================================================
c  Set up loading directions
c===================================================================
c-------------------------------------------------------------------
c  Tension & Compression
c-------------------------------------------------------------------

      if (loading_mode .eq. 'TensionComp') then

        print*,'Tension/Compression in the ',j,' direction.'
        Print*,' '
        if (i .ne. j) then
          Print*,'WARNING!!! Loading direction not clear in'
          Print*,'           input file.'
          Print*,' '
        end if
        k = j
        i = mod(k  ,3)+1
        j = mod(k+1,3)+1

c-------------------------------------------------------------------
c  Torsion
c-------------------------------------------------------------------

      else if (loading_mode .eq. 'Torsion') then

        print*,'Torsion in the ',i,' -',j,' plane.'
        Print*,' '
        if (i .eq. j) then
          Print*,'STOP!!! Loading plane must be two DIFFERENT'
          Print*,'        integers in input file.'
          Print*,' '
          STOP
        end if
        k = 6 - i - j

c-------------------------------------------------------------------
c  Plane Strain Compression
c-------------------------------------------------------------------

      else if (loading_mode .eq. 'PlaneStrain') then

        print*,'Plane-strain in the ',i,' -',j,' plane.'
        Print*,' '
        if (i .eq. j) then
          Print*,'STOP!!! Loading plane must be two DIFFERENT'
          Print*,'        integers in input file.'
          Print*,' '
          STOP
        end if
        k = 6 - i - j

c-------------------------------------------------------------------
c  Simultaneous Compression & Torsion
c-------------------------------------------------------------------

      else if (loading_mode .eq. 'SimuComTors') then

        print*,'Both Comp & Torsion in the ',i,' -',j,' plane.'
        Print*,' '
        if (i .eq. j) then
          Print*,'STOP!!! Loading plane must be two DIFFERENT'
          Print*,'        integers in input file.'
          Print*,' '
          STOP
        end if
        k = j
        j = 6 - i - k

c-------------------------------------------------------------------
c  Biaxial Tension (and Compression)
c-------------------------------------------------------------------

      else if (loading_mode .eq. 'BiaxialTens') then

        print*,'Biaxial tension in the ',i,' -',j,' plane.'
        Print*,' '
        if (i .eq. j) then
          Print*,'STOP!!! Loading plane must be two DIFFERENT'
          Print*,'        integers in input file.'
          Print*,' '
          STOP
        end if
        k = 6 - i - j

      end if

c-------------------------------------------------------------------
c  Print header
c-------------------------------------------------------------------

      write(6,'(42a,38a)')'     TIME     E_eff    SIG(1,1)  SIG(2,2)',
     & '  SIG(3,3)  SIG(1,2)  SIG(1,3) SIG(2,3)'


c-------------------------------------------------------------------
c  Call umat in order to get tangent stiffness matrix
c-------------------------------------------------------------------
      
c     print *,"before stress: ", stress
c     print *,"before ddsdde: ", ddsdde
c     print *,"before strain: ", strain
c     print *,"before dstrain: ", dstrain
c     print *,"before dfgrd1: ", dfgrd1

      call umat        (stress,  statev,  ddsdde,  sse,     spd,
     &			scd,     rpl,     ddsddt,  drplde,  drpldt,
     &			strain,  dstrain, time,    dtime,   temp,
     &			dtemp,   predef,  dpred,   cmname,  ndi,
     &			nshr,    ntens,   nstatv,  props,   nprops,
     &			coords,  drot,    pnewdt,  celent,  dfgrd0,
     &			dfgrd1,  noel,    npt,     layer,   kspt,
     &			kstep,   kinc)

      do m = 1,nstatv
        statev_ref(m) = statev(m)
      end do

c     print *,"after stress: ", stress
c     print *,"after ddsdde: ", ddsdde
c     print *,"after strain: ", strain
c     print *,"after dstrain: ", dstrain
c     print *,"after dfgrd1: ", dfgrd1

c===================================================================
c  Increments one time step
c===================================================================

  400 time(1) = time(1) + dtime
      time(2) = time(2) + dtime

c-------------------------------------------------------------------
c  Tension & Compression in the 'k' direction
c-------------------------------------------------------------------

      if (loading_mode .eq. 'TensionComp') then

      dfgrd1(k,k) = dfgrd0(k,k) + velocity * dtime
      delta_Dkk = (dfgrd1(k,k) - dfgrd0(k,k)) / dfgrd1(k,k)
      delta_Dii = -(stress(i)+stress(j) + (ddsdde(i,k)+ddsdde(j,k))
     &		  * delta_Dkk) / (ddsdde(i,i) + ddsdde(i,j) 
     &		  + ddsdde(j,i) + ddsdde(j,j))

      dfgrd1(i,i) = dfgrd0(i,i) / (1 - delta_Dii)
      dfgrd1(j,j) = dfgrd0(j,j) / (1 - delta_Dii)

c-------------------------------------------------------------------
c  Torsion in the i-j plane
c-------------------------------------------------------------------

      else if (loading_mode .eq. 'Torsion') then

      m = i + j + 1
      dfgrd1(i,j) = dfgrd0(i,j) + velocity * dtime
      delta_Dij = velocity * dtime / 2 / dfgrd1(i,i)
      delta_Djj = -(stress(j)+ddsdde(j,m)*delta_Dij) / ddsdde(j,j)
      dfgrd1(j,j) = dfgrd0(j,j) / (1 - delta_Djj)

c-------------------------------------------------------------------
c  Plane Strain Tension/Compression in the i-j plane
c-------------------------------------------------------------------

      else if (loading_mode .eq. 'PlaneStrain') then

      dfgrd1(j,j) = dfgrd0(j,j) + velocity * dtime
      delta_Djj = (dfgrd1(j,j) - dfgrd0(j,j)) / dfgrd1(j,j)
      delta_Dii = -(stress(i)+ddsdde(i,j)*delta_Djj) / ddsdde(i,i)
      dfgrd1(i,i) = dfgrd0(i,i) / (1 - delta_Dii)

c-------------------------------------------------------------------
c  Simultaneous Compress & Torsion
c-------------------------------------------------------------------

      else if (loading_mode .eq. 'SimuComTors') then

      m = i + k + 1
      dfgrd1(k,k) = dfgrd0(k,k) + velocity * dtime
      dfgrd1(i,k) = dfgrd0(i,k) + velocity2* dtime
      delta_Dkk = (dfgrd1(k,k) - dfgrd0(k,k)) / dfgrd1(k,k)
      delta_Dik = velocity2 * dtime / 2 / dfgrd1(i,i)

      det = ddsdde(i,i) * ddsdde(j,j) - ddsdde(i,j) * ddsdde(j,i)
      b1 = -stress(i) - ddsdde(i,k) * delta_Dkk
     &		- ddsdde(i,m) * delta_Dik
      b2 = -stress(j) - ddsdde(j,k) * delta_Dkk
     &		- ddsdde(j,m) * delta_Dik
      delta_Dii = (b1 * ddsdde(j,j) - b2 * ddsdde(i,j)) / det
      delta_Djj = (ddsdde(i,i) * b2 - ddsdde(j,i) * b1) / det
      dfgrd1(i,i) = dfgrd0(i,i) / (1 - delta_Dii)
      dfgrd1(j,j) = dfgrd0(j,j) / (1 - delta_Djj)

c-------------------------------------------------------------------
c  Biaxial Tension in the i - j plane
c-------------------------------------------------------------------

      else if (loading_mode .eq. 'BiaxialTens') then

      dfgrd1(i,i) = dfgrd0(i,i) + velocity * dtime
      dfgrd1(j,j) = dfgrd0(j,j) + velocity2* dtime
      delta_Dii = (dfgrd1(i,i) - dfgrd0(i,i)) / dfgrd1(i,i)
      delta_Djj = (dfgrd1(j,j) - dfgrd0(j,j)) / dfgrd1(j,j)
      delta_Dkk = (stress(k) - ddsdde(i,i) * delta_Dii
     &		  - ddsdde(j,j) * delta_Djj)
     &		  / ddsdde(k,k)

      dfgrd1(k,k) = dfgrd0(k,k) / (1 - delta_Dkk)

      end if

c===================================================================
c  Start new convergence iteration.  Restart with time-step
c  cut in half if necessary
c===================================================================

      Kinc = 0
  500 if (Kinc .ge. Max_itr) then
        time(1) = time(1) - dtime
        time(2) = time(2) - dtime
        dtime = dtime / 2
        go to 400
      end if
      Kinc = Kinc + 1

c-------------------------------------------------------------------
c  Calculate F_dot * dtime
c-------------------------------------------------------------------

      do m = 1,3
        do n = 1,3
           array1(m,n) = dfgrd1(m,n) - dfgrd0(m,n)
           array3(m,n) =(dfgrd1(m,n) + dfgrd0(m,n))/2
        end do
      end do

c-------------------------------------------------------------------
c  multiply F_dot * F_inverse * dtime
c-------------------------------------------------------------------

      call Xinverse_3x3(array3,array2)
      call Xaa_dot_bb(3,array1,array2,array3)

c-------------------------------------------------------------------
c  Get D_dt and W_dt
c-------------------------------------------------------------------

      do m = 1,3
        do n = 1,3
          D_dt(m,n) = (array3(m,n) + array3(n,m)) / 2
          W_dt(m,n) = (array3(m,n) - array3(n,m)) / 2
        end do
      end do

c-------------------------------------------------------------------
c  Store D_dt in dstrain
c-------------------------------------------------------------------

      dstrain(1) = D_dt(1,1)
      dstrain(2) = D_dt(2,2)
      dstrain(3) = D_dt(3,3)
      dstrain(4) = D_dt(1,2) * 2
      dstrain(5) = D_dt(1,3) * 2
      dstrain(6) = D_dt(2,3) * 2

c-------------------------------------------------------------------
c  Convert spin to drot(i,j) array for the UMAT
c-------------------------------------------------------------------

      call Xspin_to_matrix(W_dt,drot)

c-------------------------------------------------------------------
c  Call umat in order to get stress
c-------------------------------------------------------------------
c      print *,"ddsdde: ", ddsdde

      call umat(stress,  statev,  ddsdde,  sse,     spd,
     &			scd,     rpl,     ddsddt,  drplde,  drpldt,
     &			strain,  dstrain, time,    dtime,   temp,
     &			dtemp,   predef,  dpred,   cmname,  ndi,
     &			nshr,    ntens,   nstatv,  props,   nprops,
     &			coords,  drot,    pnewdt,  celent,  dfgrd0,
     &			dfgrd1,  noel,    npt,     layer,   kspt,
     &			kstep,   kinc)

c     print *, "-------------------------------------------"
c     print *, "stress: ", stress
c     print *, "ddsdde: ", ddsdde
c     print *, "dfgrd0: ", dfgrd0
c     print *, "dfgrd1: ", dfgrd1

c-------------------------------------------------------------------
c  Check to see if need to iterate 
c  Tension & Compression are immediately below
c-------------------------------------------------------------------

      if (loading_mode .eq. 'TensionComp') then

        test = abs(stress(i)+stress(j))/abs(stress(k))

      if (test .gt. tolerance) then
        delta_Dii = -(stress(i)+stress(j))/(ddsdde(i,i)+ddsdde(i,j) 
     &		  + ddsdde(j,i) + ddsdde(j,j))
        dfgrd1(i,i) = dfgrd1(i,i) / (1 - delta_Dii)
        dfgrd1(j,j) = dfgrd1(j,j) / (1 - delta_Dii)
        do m = 1,nstatv
          statev(m) = statev_ref(m)
        end do
        go to 500
      end if

c-------------------------------------------------------------------
c  Torsion
c-------------------------------------------------------------------

      else if (loading_mode .eq. 'Torsion') then

      m = i + j + 1
      test = abs(stress(j)/stress(m))

      if (test .gt. tolerance) then
        delta_Djj = -stress(j) / ddsdde(j,j)
        dfgrd1(j,j) = dfgrd1(j,j) / (1 - delta_Djj)
        do m = 1,nstatv
          statev(m) = statev_ref(m)
        end do
        go to 500
      end if

c-------------------------------------------------------------------
c  Plane Strain Compresssion
c-------------------------------------------------------------------

      else if (loading_mode .eq. 'PlaneStrain') then

      test = abs(stress(i)/stress(j))

      if (test .gt. tolerance) then
        delta_Dii = -stress(i) / ddsdde(i,i)
        dfgrd1(i,i) = dfgrd1(i,i) / (1 - delta_Dii)
        do m = 1,nstatv
          statev(m) = statev_ref(m)
        end do
        go to 500
      end if

c-------------------------------------------------------------------
c  Simultaneous Compress & Torsion
c-------------------------------------------------------------------

      else if (loading_mode .eq. 'SimuComTors') then

      m = i + k + 1
      test = (abs(stress(i)) + abs(stress(j))) /
     &	 (abs(stress(k)) + abs(stress(m)))

      if (test .gt. tolerance) then
        det = ddsdde(i,i) * ddsdde(j,j) - ddsdde(i,j) * ddsdde(j,i)
        b1 = -stress(i)
        b2 = -stress(j)
        delta_Dii = (b1 * ddsdde(j,j) - b2 * ddsdde(i,j)) / det
        delta_Djj = (ddsdde(i,i) * b2 - ddsdde(j,i) * b1) / det
        dfgrd1(i,i) = dfgrd1(i,i) / (1 - delta_Dii)
        dfgrd1(j,j) = dfgrd1(j,j) / (1 - delta_Djj)
        do m = 1,nstatv
          statev(m) = statev_ref(m)
        end do
        go to 500
      end if

c-------------------------------------------------------------------
c  Biaxial Tension
c-------------------------------------------------------------------

      else if (loading_mode .eq. 'BiaxialTens') then

      test = 2 * abs(stress(k)) / (abs(stress(i)) + abs(stress(j)))

      if (test .gt. tolerance) then
        delta_Dkk = -stress(k) / ddsdde(k,k)
        dfgrd1(k,k) = dfgrd1(k,k) / (1 - delta_Dkk)
        do m = 1,nstatv
          statev(m) = statev_ref(m)
        end do
        go to 500
      end if

      end if

c===================================================================
c  Finished increment!
c  Calc effective delta strain and add it to E_eff
c===================================================================

      sum = 0
      do m = 1,3
        do n = 1,3
          sum = sum + D_dt(m,n) * D_dt(m,n)
        end do
      end do
      dE_eff = sqrt(2. * sum / 3.)
      E_eff = E_eff + dE_eff

c-------------------------------------------------------------------
c  Update strain gradient
c-------------------------------------------------------------------

      do m = 1,6
        strain(m) = strain(m) + dstrain(m)
      end do

c-------------------------------------------------------------------
c  Update deformation gradient
c-------------------------------------------------------------------

      do m = 1,3
        do n = 1,3
          dfgrd0(m,n) = dfgrd1(m,n)
        end do
      end do

c-------------------------------------------------------------------
c  Update statev_ref
c-------------------------------------------------------------------

      do m = 1,nstatv
        statev_ref(m) = statev(m)
      end do

c-------------------------------------------------------------------
c  Write out results
c-------------------------------------------------------------------

      write(6,'(f10.4,f10.5,6(f10.3))')time(2),E_eff,
     &	(stress(m),m=1,6)
c     & (dfgrd1(n,n),n=1,3)


      write(8,'(f10.4,f10.5,6(f10.3))')time(2),E_eff,
     &	(stress(m),m=1,6)
c     & (dfgrd1(n,n),n=1,3)

c-------------------------------------------------------------------
c  If not out of time, then loop back up
c-------------------------------------------------------------------

      if (.not.Out_of_time) then
        if (time(1).ne.dtime) dtime = dtime * 1.5
        if (dtime .gt. dtime_max) dtime = dtime_max
        if (time_max-time(1).lt.dtime) then
          dtime = time_max - time(1)
          Out_of_time = .true.
        end if
        go to 400
      end if

c-------------------------------------------------------------------
c  Loop back up and check for another STEP
c-------------------------------------------------------------------

      go to 300

c-------------------------------------------------------------------
c  THE END!!!
c-------------------------------------------------------------------

  999 continue

      close(1)
      close(7)
      close(8)
      stop
      end

c====================================================================
c====================================================================
c
c  Calculate the dot product of two 2nd rank tensors.
c  Result is stored in cc(i,j)
c
c--------------------------------------------------------------------

      subroutine Xaa_dot_bb(n,a,b,c)

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
c  Calculate the inverse of a 3 x 3 matrix.
c
c--------------------------------------------------------------------

      subroutine Xinverse_3x3(a,b)

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
c  Convert a spin tensor to a rotation matrix 
c
c--------------------------------------------------------------------

      subroutine Xspin_to_matrix(a,b)

      implicit double precision (a-h,o-z)

      dimension a(3,3), b(3,3)

      Pi = 4*atan(1D0)

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
