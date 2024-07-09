

c--------------------------------------------------------------------
c     This umat is called by ale3d and wraps around a regular
c     abaqus umat call. Main differences are cmname, and initialization
c     is handled differently.
c--------------------------------------------------------------------
      subroutine umat ( stress,  statev,  ddsdde,  sse,     spd,
     &                  scd,     rpl,     ddsddt,  drplde,  drpldt,
     &                  strain,  dstrain, time,    dtime,   temp,
     &                  dtemp,   predef,  dpred,   cmname,  ndi,
     &                  nshr,    ntens,   nstatv,  props,   nprops,
     &                  coords,  drot,    pnewdt,  celent,  dfgrd0,
     &                  dfgrd1,  noel,    npt,     layer,   kspt,
     &                  kstep,   kinc )
     
      implicit none

      ! Dimension variables passed into the UMAT sub (not all are used)
      integer ndi, nshr, ntens, nstatv, nprops, noel 
      integer layer, kspt, kstep, kinc, npt 
      double precision cmname, sse, celent, dtime, temp, dtemp, pnewdt
      double precision spd, scd, rpl, drpldt
      character(8) cmnameale ! Material name
      double precision coords(3), ddsdde(6,6), ddsddt(ntens)
      double precision dfgrd1(3,3), dfgrd0(3,3), dpred(1), drplde(ntens)
      double precision drot(3,3), dstrain(ntens), predef(1)
      double precision props(nprops), statev(nstatv), strain(ntens)
      double precision stress(ntens), time(2) 

c     Helper variables
      cmnameale = 'UMAT_JC'

      if (nstatv.ne.6) then
        print*, 'Wrong number of state variables in cumat'
      end if
      

      call cumat(stress,  statev,  ddsdde,  sse,     spd,
     &  scd,     rpl,     ddsddt,  drplde,  drpldt,
     &  strain,  dstrain, time,    dtime,   temp,
     &  dtemp,   predef,  dpred,   cmnameale,  ndi,
     &  nshr,    ntens,   nstatv,  props,   nprops,
     &  coords,  drot,    pnewdt,  celent,  dfgrd0,
     &  dfgrd1,  noel,    npt,     layer,   kspt,
     &  kstep,   kinc)          
      
      return
      end

c--------------------------------------------------------------------
c     ABAQUS STRESS - SIG11, SIG22, SIG33, SIG12, SIG13, SIG23 
c     ABAQUS STRAIN - EPS11, EPS22, EPS33, GAM12, GAM13, GAM23
c                     WHERE GAM12 = 2*EPS12
c--------------------------------------------------------------------

      subroutine cumat ( stress,  statev,  ddsdde,  sse,     spd,
     &                  scd,     rpl,     ddsddt,  drplde,  drpldt,
     &                  strain,  dstrain, time,    dtime,   temp,
     &                  dtemp,   predef,  dpred,   cmname,  ndi,
     &                  nshr,    ntens,   nstatv,  props,   nprops,
     &                  coords,  drot,    pnewdt,  celent,  dfgrd0,
     &                  dfgrd1,  noel,    npt,     layer,   kspt,
     &                  kstep,   kinc )

      IMPLICIT NONE

      ! loop variables
      integer i
      
      ! Dimension variables passed into the UMAT sub (not all are used)
      integer ndi      ! Number of direct stress components
      integer nshr     ! Number of shear stress components
      integer ntens    ! Size of stess or stran array (ndi + nshr)
      integer nstatv   ! Number of SDVs
      integer nprops   ! Number of material constants
      integer noel     ! Element number
      integer layer    ! Layer number (for composites)
      integer kspt     ! Section point number within layer
      integer kstep    ! Step number
      integer kinc     ! Increment number
      integer npt      ! Integration point number
      character(8) cmname ! Material name
      double precision sse ! Specific elastic stain energy
      
      double precision
     & celent,         ! Characteristic element length
     & dtime,          ! Time increment
     & temp,           ! Temperature at start of increment
     & dtemp,          ! Temperature increment
     & pnewdt,         ! Ratio of new time increment to time
                       ! increment being used
     & spd,            ! Specific plastic dissipation
     & scd,            ! Specific creep dissipation
     & rpl,            ! Volumetic heat generation per unit time
     & drpldt,         ! Varation of rpl with temperature
     & coords(3),      ! Coordinates of Gauss pt. being evaluated
     & ddsdde(ntens,ntens), ! Tangent Stiffness Matrix
     & ddsddt(ntens),  ! Change in stress per change in temperature
     & dfgrd1(3,3),    ! Deformation gradient at end of step
     & dfgrd0(3,3),    ! Deformation gradient at beginning of step
     & dpred(1),       ! Change in predefined state variables
     & drplde(ntens),  ! Change in heat generation per change in strain
     & drot(3,3),      ! Rotation matrix
     & dstrain(ntens), ! Strain increment tensor stored in vector form
     & predef(1),      ! Predefined state vars dependent on field
                       ! variables
     & props(nprops),  ! Material properties passed in
     & statev(nstatv), ! State Variables
     & strain(ntens),  ! Strain tensor stored in vector form
     & stress(ntens),  ! Cauchy stress tensor stored in vector form
     & time(2)         ! Step Time and Total Time

      !Variables fed in from props 
      double precision lam      ! Elastic lambda (MPa)
      double precision mu       ! Shear modulus (MPa)
      double precision temp0    ! Initial temperature
      double precision tempmelt ! Melt temperature
c
      integer hsltype           ! Slip hardening model used
      double precision hsl1     ! Slip hardening parameter 1
      double precision hsl2     ! Slip hardening parameter 2
      double precision hsl3     ! Slip hardening parameter 3
      double precision hsl4     ! Slip hardening parameter 4
      double precision hsl5     ! Slip hardening parameter 5
      double precision hsl6     ! Slip hardening parameter 6
      double precision hsl7     ! Slip hardening parameter 7
      double precision hsl8     ! Slip hardening parameter 8
c
      integer dmgtype           ! Dmg model used
      double precision d1       ! Dmg parameter 1
      double precision d2       ! Dmg parameter 2
      double precision d3       ! Dmg parameter 3
      double precision d4       ! Dmg parameter 4
      double precision d5       ! Dmg parameter 5
      double precision d6       ! Dmg parameter 6
      double precision d7       ! Dmg parameter 7
      double precision d8       ! Dmg parameter 8
      double precision d9       ! Dmg parameter 9
      double precision d10      ! Dmg parameter 10
c      
      integer eosflag           ! 0 lin elast     
c                                 1 murn eos with cold curve
      double precision b0       ! Bulk modulus, Guinan and Stein, 74
      double precision dbdp     ! Bulk mod deriv w.r.t p, "" ref 
      double precision cv       ! MPa/K, Lee, Int J Thermophys, 13
      double precision grun     ! Grun coeff 
      double precision gruna    ! Volume dependence of Gruneisen paramater
      double precision ecolds   ! compression shift for cold curve
      double precision ecold0   ! cold curve value at minimum
      double precision ecold1   ! linear dependence of cold curve on compression
      double precision ecold2   ! quadratic dependence of cold curve oncompression
c
      !Variables related to statev
      double precision epsl0       ! Nonbas effplastic strain prev step
      double precision epdsl       ! "      " rate 
      double precision epdsl0      ! "      " from prev step step
      double precision energy      ! Int energy / ref vol
      double precision tempsv      ! Temperature stored as state var
      double precision dmg         ! Dmg variable 0-1
      double precision dmgflag     ! Dmg flag 0 or 1
c      
      !Utility
      logical actsl !whether modes are active
      double precision gsl   !trial shear on these
      double precision yssl0, yssl !initial and cur ys
      double precision dyssldepd       !change in ys wrt epd
      double precision depsl       !change in strains
      double precision p0, dp, p     !initial pressure, and increment
      double precision bmod       !bulk modulus
      double precision epdslmax   !estimate of max slip strain rate 
c      
      double precision 
     & stressd0(3,3),    ! Deviatoric trial stress
     & stressd(3,3),     ! Deviatoric stress
     & wpdt(3,3),        ! Plastic spin times dt
     & wdt(3,3)          ! Spin times dt (from abaqus drot)

c     very solution oriented
      double precision actTOL

c--------------------------------------------------------------------
c     Solution related parameters
c--------------------------------------------------------------------
      actTOL = 1.0d-8

c--------------------------------------------------------------------
c     Read in material properties
c--------------------------------------------------------------------
! c       sigy = (A+B*epsl**n)(1+C*ln(epdsl))(1-((t-t0)/(tmelt-t0))**m)
! c       hsl1 = A
! c       hsl2 = B
! c       hsl3 = n
! c       hsl4 = C
! c       hsl5 = m
! c       hsl6-7 = not used
! c       hsl8 = cutoff strain rate

      lam = props(1)
      mu = props(2)
      temp0 = props(3)
      tempmelt = props(4)
      hsltype = nint(props(5))     
      hsl1 = props(6)
      hsl2 = props(7)
      hsl3 = props(8)
      hsl4 = props(9)
      hsl5 = props(10)
      hsl6 = props(11)
      hsl7 = props(12)
      hsl8 = props(13)
      dmgtype = nint(props(14))
      d1 = props(15)
      d2 = props(16)
      d3 = props(17)
      d4 = props(18)
      d5 = props(19)
      d6 = props(20)
      d7 = props(21)
      d8 = props(22)
      d9 = props(23)
      d10 = props(24) !min tensile pressure if dmgd, for -5MPa its 5.0 
      eosflag = nint(props(25))
      cv = props(26)
      b0 = props(27)
      dbdp = props(28)
      grun = props(29)
      gruna = props(30)
      ecolds = props(31)
      ecold0 = props(32)
      ecold1 = props(33)
      ecold2 = props(34)


c--------------------------------------------------------------------
c     Initialize state variables on first step
c--------------------------------------------------------------------

	  if ((time(2)).le.(0.5d+0*dtime)) then            
        statev(1) = 0.0d+0 !epsilon slip
        statev(2) = 0.0d+0 !epsilondot slip
        statev(3) = 0.0d+0 !energy per unit reference volume
        statev(4) = temp !temperature stored in isv       
        statev(5) = 0.0d+0 !damage parameter continuous 0-1
        statev(6) = 0.0d+0 !damage flag 0 or 1
		
		if (noel.eq.1) then
			print *, noel, time(2), dtime, 'initialized'
		endif
		if (noel.eq.1) then
			print *, noel, statev(4), 'starting temp'
		endif	  
c	  		
	  else
	    statev(4) = temp
c	    !statev(4) = props(3)		
		if (noel.eq.1) then
			print *, noel, time(1), 'over 300 time'
		endif
		if (noel.eq.1) then
			print *, noel, statev(4), 'over 300 temp'
		endif
c
	  end if
      ! if ((time(1)-dtime).le.(0.5d+0*dtime)) then            
        ! statev(1) = 0.0d+0 !epsilon slip
        ! statev(2) = 0.0d+0 !epsilondot slip
        ! statev(3) = 0.0d+0 !energy per unit reference volume
! !        statev(4) = props(3) !temperature stored in isv       
        ! statev(4) = temp !temperature stored in isv (changed by mr1751)    
        ! statev(5) = 0.0d+0 !damage parameter continuous 0-1
        ! statev(6) = 0.0d+0 !damage flag 0 or 1
      ! end if 
	  ! if (noel.eq.1) then
		  ! print *, noel, time(1), dtime, 'initialized'
      ! endif
	  ! if (noel.eq.1) then
		  ! print *, noel, statev(4), 'starting temp'
	  ! endif	

c--------------------------------------------------------------------
c     Read in state variables - things from prev step
c--------------------------------------------------------------------
      call init_statevs(statev, nstatv, epsl0, epdsl0, energy, tempsv,
     &   dmg, dmgflag)      
      temp = tempsv
	  
 	  if (noel.eq.1) then
		  print *, noel, epsl0, epdsl0, tempsv, dmg, dmgflag,
     1'values after init_statevs'
	  endif	
 	  
      p0 = - (stress(1) + stress(2) + stress(3))/3.0d+0
	  
	  if (noel.eq.1) then
		  print *, noel, p0, stress(1),stress(2), stress(3),
     1'po and stress after init_statevs'
	  endif	  
c--------------------------------------------------------------------
c     Check if fully damaged - if so skip all plasticity
c--------------------------------------------------------------------      
      if (dabs(dmgflag-1.0d0).le.0.1d0) then !damaged
        dmg = 1.0d0
        dmgflag = 1.0d0
        call skip_plasticity(yssl, depsl, epdsl, stressd)
      else !undamaged      
c--------------------------------------------------------------------
c     Determine whether or not to bypass plastic routine
c--------------------------------------------------------------------   
	  if (noel.eq.1) then
		  print *, noel, stress, dstrain, mu, 'pre calcuation'
	  endif
c     Calculate trial stress tensor
      call calc_trial_devstress(stress, dstrain, mu, stressd0)
	  
	  if (noel.eq.1) then
		  print *, noel, stress, dstrain, mu, stressd0,
     1'calc trial devstress'
	  endif
	  
c     Calculate trial stresses on slip
      call calc_taus(stressd0, gsl)
	  
	  if (noel.eq.1) then
		  print *, noel, stressd0, gsl, 'calc taus'
	  endif
	  
c     Calculate strength from previous step - yssl0
      call calc_str_sl(.true., epsl0, epdsl0, temp, tempmelt,
     &  hsltype, hsl1, hsl2, hsl3, hsl4, hsl5, hsl6, hsl7, hsl8,
     &  yssl0, dyssldepd)
	 
	  if (noel.eq.1) then
		  print *, noel, yssl0, dyssldepd, 'calc str_s1'
	  endif	

c     Determine if potentially activate based on the trial stress      
      call calc_potactive(gsl, yssl0, actTOL, actsl)
	  
	  if (noel.eq.1) then
		  print *, noel, gsl, yssl0, actTOL, actsl, 'calc potactive'
	  endif	
c====================================================================
c     Start the solve for increment in plastic strain 
c====================================================================

c     Initialization for loop
      depsl = 0.0d+0
      epdsl = 0.0d+0
      yssl = yssl0

	  if (noel.eq.1) then
		  print *, noel, depsl, gsl, epdsl, yssl, 'before plastic'
	  endif	
	  
c     Calculate epdsl, depsl, yssl if potentially active
      if (actsl) then
c       Calc upper bound for epdsl, hsl8 is min epdsl allowed        
        call calc_max_epdsl(dstrain,gsl,lam,mu,dtime,hsl8,epdslmax)
        if (noel.eq.1) then
		    print *, noel, dstrain, gsl, lam, mu, dtime, hsl8, epdslmax,
     1'calc max epds1'
	    endif	
c       Calculate epdsl              
        call calc_epdsl(epsl0, gsl, 
     &        hsltype, hsl1, hsl2, hsl3, hsl4, hsl5, hsl6, hsl7, hsl8,
     &        tempmelt, stressd0, mu, temp, dtime, epdslmax, epdsl)

        depsl = epdsl * dtime

c       Calculate yssl
        call calc_str_sl(.true., epsl0+depsl, epdsl, temp, tempmelt, 
     &    hsltype,hsl1,hsl2,hsl3,hsl4, hsl5, hsl6, hsl7, hsl8, 
     &    yssl, dyssldepd)

	    if (noel.eq.1) then
		    print *, noel, yssl, dyssldepd, depsl, epdsl,
     1'calc str sl flag'
	    endif
		
        if (epdsl.lt.hsl8) then
          epdsl = 0.0d+0
          depsl = 0.0d+0
        end if      
      else
        epdsl = 0.0d+0
      end if 

c====================================================================
c     End plastic strain loop 
c====================================================================
      
c     At the end: calculate the deviatoric part of the stress     
      call calc_return_dev_stress(stressd0, mu, depsl, yssl, stressd) 
      end if !end damage loop
	  
	  if (noel.eq.1) then
		  print *, noel, stressd0, mu, depsl, yssl, stressd,
     1'calc return dev stress'
	  endif	
c--------------------------------------------------------------------
c     End damage loop
c--------------------------------------------------------------------

c====================================================================
c     Equation of state 
c====================================================================
      
c     Isothermal, linear elastic eos      
      if (eosflag.eq.0) then
        bmod = lam+2.0d+0*mu/3.0d+0
        dp = - bmod*(dstrain(1) + dstrain(2) + dstrain(3))
        p = p0+dp
!        tempsv = temp0 
		tempsv = temp
		
	    if (noel.eq.1) then
		    print *, noel, bmod, dp, p, tempsv, 'eos flag'
	    endif		
		
      else
        p = p0
c       Input/output - p and energy. Output - tempsv    
        call eos(dfgrd0, dfgrd1, b0, dbdp, grun, gruna, cv,
     &    yssl, depsl, ecolds, ecold0, ecold1, ecold2, dmgflag, d10,
     &    p, energy, tempsv)  
      end if

c====================================================================
c     Updates - including new damage state
c====================================================================

c
      stress(1) = stressd(1,1) - p
      stress(2) = stressd(2,2) - p
      stress(3) = stressd(3,3) - p
      stress(4) = stressd(1,2)
      stress(5) = stressd(1,3)      
      stress(6) = stressd(2,3)  
	  
	  
	  if (noel.eq.1) then
		  print *, noel, stress, 'updated stress'
	  endif	
	  
c     not necessary for eos implementation
      rpl = energy - statev(3)

	  if (noel.eq.1) then
		  print *, noel, rpl, 'updated rpl'
	  endif	
      
c     Update state variables      
      call update_statevs(depsl, epdsl, energy, tempsv, nstatv, statev) 

	  if (noel.eq.1) then
		  print *, noel, depsl, epdsl, energy, tempsv, nstatv, statev,
     1'update_statevs'
	  endif
      
c     Calculate damage - only used on next time step      
      call calc_dmg(dtime, tempmelt, stress, dmgtype, 
     &  d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, nstatv, statev)

 	  if (noel.eq.1) then
		  print *, noel, dtime, tempmelt, stress, dmgtype, d1,
     1d2, d3, d4, d5, d6, d7, d8, d9, d10, nstatv, statev,
     2'calc damage'
	  endif
c     Return ddsdde for implicit - not correct
      call elastddsdde(lam, mu, ntens, ndi, ddsdde)
 
 	  if (noel.eq.1) then
		  print *, noel, lam, mu, ntens, ndi, ddsdde, 'elastddsdde'
	  endif
	  
      return
      end



      
c====================================================================
c====================================================================
c  Read in state variables array to other variables
c--------------------------------------------------------------------

      subroutine init_statevs(statev, nstatv, epsl0, epdsl0, energy, 
     &  tempsv, dmg, dmgflag)     
      implicit none
      
c     input
      integer nstatv
      double precision statev(nstatv)
      
c     output      
      double precision epsl0, epdsl0, energy, tempsv, dmg, dmgflag
c      
      !Read in statev
      epsl0 = statev(1)
      epdsl0 = statev(2)
      energy = statev(3)
      tempsv = statev(4)
      dmg = statev(5)
      dmgflag = statev(6)

      return
      end

c====================================================================
c====================================================================
c     Zero out all vars needed for eos calc if material is damaged
c--------------------------------------------------------------------
      subroutine skip_plasticity(yssl, depsl, epdsl, stressd)
      implicit none

c     input/output 
      double precision yssl, depsl, epdsl, stressd(3,3)      
      
      yssl = 0.0d0
      depsl = 0.0d0
      epdsl = 0.0d0
      stressd(1,1) = 0.0d0
      stressd(1,2) = 0.0d0
      stressd(1,3) = 0.0d0
      stressd(2,1) = 0.0d0
      stressd(2,2) = 0.0d0
      stressd(2,3) = 0.0d0
      stressd(3,1) = 0.0d0
      stressd(3,2) = 0.0d0
      stressd(3,3) = 0.0d0
      
      return
      end
      
c====================================================================
c====================================================================
c  Return deviatoric trial stress in indicial notation
c--------------------------------------------------------------------

      subroutine calc_trial_devstress(STRESS, DSTRAIN, MU, SIGTDEV)

      implicit none
      
c     input
      double precision STRESS(6), DSTRAIN(6)
      double precision MU
      
c     output
      double precision SIGTDEV(3,3)
      
c     util
      double precision DEPSH, SIGH0 
      double precision DEPSD1, DEPSD2, DEPSD3, DEPSD4, DEPSD5, DEPSD6
      
C DEVIATORIC PART OF THE STRAIN INCREMENT (ABQ stores eng instead tens)

      DEPSH=(DSTRAIN(1)+DSTRAIN(2)+DSTRAIN(3))/3.0d+0
      DEPSD1=DSTRAIN(1)-DEPSH
      DEPSD2=DSTRAIN(2)-DEPSH
      DEPSD3=DSTRAIN(3)-DEPSH
      DEPSD4=DSTRAIN(4)/2.0d+0
      DEPSD5=DSTRAIN(5)/2.0d+0
      DEPSD6=DSTRAIN(6)/2.0d+0

      SIGH0=(STRESS(1)+STRESS(2)+STRESS(3))/3.0d+0
      SIGTDEV(1,1)=STRESS(1)-SIGH0+2.0d+0*MU*DEPSD1
      SIGTDEV(2,2)=STRESS(2)-SIGH0+2.0d+0*MU*DEPSD2
      SIGTDEV(3,3)=STRESS(3)-SIGH0+2.0d+0*MU*DEPSD3
      SIGTDEV(1,2)=STRESS(4)+2.0d+0*MU*DEPSD4
      SIGTDEV(2,1)=STRESS(4)+2.0d+0*MU*DEPSD4      
      SIGTDEV(1,3)=STRESS(5)+2.0d+0*MU*DEPSD5
      SIGTDEV(3,1)=STRESS(5)+2.0d+0*MU*DEPSD5
      SIGTDEV(2,3)=STRESS(6)+2.0d+0*MU*DEPSD6
      SIGTDEV(3,2)=STRESS(6)+2.0d+0*MU*DEPSD6

c      SIGT=SIGDEV(1)**2+SIGDEV(2)**2+SIGDEV(3)**2+2.0d+0*(SIGDEV(4)**2
c     & + SIGDEV(5)**2 + SIGDEV(6)**2)
c      SIGT=DSQRT(3.0d+0/2.0d+0*SIGT)


      return
      end

c====================================================================
c====================================================================
c  Calculate driving forces for yield based on trial stress 
c--------------------------------------------------------------------

      subroutine calc_taus(sigd, gsl)
      
      implicit none

      !input
      double precision sigd(3,3)

      !output
      double precision gsl

      gsl = dsqrt(3.0d+0/2.0d+0*(sigd(1,1)**2+
     & sigd(2,2)**2+sigd(3,3)**2+2.0d+0*(sigd(1,2)**2 + 
     & sigd(1,3)**2 + sigd(2,3)**2)))      

      return
      end      
      
c====================================================================
c====================================================================
c  Calculate strength for slip
c  -  bis just means you are in a bisection solve and dont need tangent
c--------------------------------------------------------------------

      subroutine calc_str_sl(bis, epsl, epdsl, temp, tempmelt,
     & hsltype, hsl1, hsl2, hsl3, hsl4, hsl5, hsl6, hsl7, hsl8, 
     & yssl, dyysdepd)

      implicit none
      
      !input
      logical bis
      double precision epsl, epdsl, temp, tempmelt
      integer hsltype 
      double precision hsl1, hsl2, hsl3, hsl4, hsl5, hsl6, hsl7, hsl8
      
      !output - yield strength, and deriv of ys wrt strain rate
      double precision yssl, dyysdepd
      
      !util
      double precision thom
      
c     IF MELTED, SET YSSL = 0, OTHERWISE CALC THOM
	  ! if temp. ge. tempmelt
		
      ! if (temp.ge.tempmelt) then
        ! yssl = 0.01d+0
	  if (temp.ge.tempmelt) then
		temp = tempmelt + 1.0d+0
	  end if
	  
        thom = (temp-293d+0)/(tempmelt-293d+0)
		!mr1751 changed treference from 293 to 20C)
	  
        if (hsltype.eq.1) then
c--------------------------------------------------------------------
c       hsltype = 1: Johnson-Cook, 
c                    uses a min epdsl (str rate) of props(8)
c--------------------------------------------------------------------
c       sigy = (A+B*epsl**n)(1+C*ln(epdsl))(1-((t-t0)/(tmelt-t0))**m)
c       hsl1 = A
c       hsl2 = B
c       hsl3 = n
c       hsl4 = C
c       hsl5 = m
c       hsl6-7 = not used
c       hsl8 = cutoff strain rate

          if (epdsl.le.hsl8) then
            epdsl = hsl8     
          end if  
          yssl=(hsl1+hsl2*dabs(epsl)**hsl3)*(1.0d+0+hsl4*dlog(epdsl))
     &      *(1.0d+0-thom**hsl5) 
     
          if (bis) return

          dyysdepd=hsl4*(1.0d+0-thom**hsl5)*(hsl1+hsl2*dabs(epsl)**hsl3)
          dyysdepd = dyysdepd / epdsl
               
        else if (hsltype.eq.2) then
c--------------------------------------------------------------------
c       hsltype = 2: Chang and Kochmann
c--------------------------------------------------------------------
c       sigy = (tau0+sig0*(1-exp(-h1*eps/sig0))+h2*eps)*(edot/edot0)**m
c       hsl1 = tau0
c       hsl2 = sig0
c       hsl3 = epd0
c       hsl4 = m
c       hsl5 = h1
c       hsl6 = h2
c       hsl7 = not used
c       hsl8 = cutoff strain rate

          if (epdsl.le.hsl8) then
            epdsl = hsl8              
          end if
          yssl=(hsl1+hsl2*(1.0d+0-dexp(-hsl5*epsl/hsl2))+hsl6*epsl)*
     &      (epdsl/hsl3)**hsl4       

          if (bis) return
          dyysdepd = hsl4/hsl3*(epdsl/hsl3)**(hsl4-1.0d+0)*
     &      (hsl1+hsl2*(1.0d+0-dexp(-hsl5*epsl/hsl2))+hsl6*epsl) 

        else
          print*, 'NB SLIP STRENGTH TYPE ', hsltype, ' NOT IMPLEMENTED'
        end if 
      
      return
      end    

c====================================================================
c====================================================================
c  Calculate which deformation mechs are potentially active based on
c    the trial stress
c--------------------------------------------------------------------

      subroutine calc_potactive(gsl, yssl, actTOL, actsl)
        
      implicit none
      
c     input
      double precision gsl
      double precision yssl, actTOL
      
c     output
      logical actsl
            
c     util
      double precision diff
      
c     slip
      diff = dabs(gsl) - yssl
      if (diff.ge.actTOL) then
        actsl = .true. 
      else
        actsl = .false.
      endif   

      return
      end

c====================================================================
c====================================================================
c  Calculate conservative estimate of maximum von mises strain rate 
c--------------------------------------------------------------------

      subroutine calc_max_epdsl(DSTRAN, SIGT, LAM, MU, DTIME, 
     &  EPSDCUTOFF, EPSDMAX) 

      implicit none

c     input
      double precision DSTRAN(6), SIGT, LAM, MU, DTIME, EPSDCUTOFF

c     output
      double precision EPSDMAX


c     util
      double precision DEPSD1, DEPSD2, DEPSD3, DEPSD4, DEPSD5, DEPSD6
      double precision DEPSH, DEPSE 
      double precision EMOD
      
      EMOD = MU*(3.0D+0*LAM+2.0D+0*MU)/(LAM+MU)

      DEPSH=(DSTRAN(1)+DSTRAN(2)+DSTRAN(3))/3.0d+0
      DEPSD1=DSTRAN(1)-DEPSH
      DEPSD2=DSTRAN(2)-DEPSH
      DEPSD3=DSTRAN(3)-DEPSH
      DEPSD4=DSTRAN(4)/2.0d+0
      DEPSD5=DSTRAN(5)/2.0d+0
      DEPSD6=DSTRAN(6)/2.0d+0 
      DEPSE=DEPSD1**2+DEPSD2**2+DEPSD3**2+2.0d+0*(DEPSD4**2+DEPSD5**2 
     & + DEPSD6**2)        
      DEPSE=DSQRT(2.0d+0/3.0d+0*DEPSE+EPSDCUTOFF)
      EPSDMAX=2.0d+0*(DEPSE+SIGT/EMOD)/DTIME

      return
      end
      
c====================================================================
c====================================================================
c  Calculate epd, von mises stress for nb slip in absence of 
c  basal slip and twin increments (still uses init vals to calc sl str)
c--------------------------------------------------------------------

      subroutine calc_epdsl(epsl0, gsl, 
     &  hsltype, hsl1, hsl2, hsl3, hsl4, hsl5, hsl6, hsl7, hsl8,
     &  tempmelt, stressd0, mu, temp, dt, epdslmax, epdsl)

      implicit none

c     input      
      double precision epsl0, gsl
c     - slip params and thermodynamic vars
      integer hsltype
      double precision hsl1, hsl2, hsl3, hsl4, hsl5, hsl6, hsl7, hsl8
      double precision tempmelt
c     - general
      double precision stressd0(3,3), mu, temp
c     - solution related
      double precision dt, epdslmax
      
c     output 
      double precision epdsl

c     util - solution related things
      logical bis
      double precision tempvar
      double precision FTOL, RTOL
      double precision X1, X2, F, FL, FH, DF, XL, XH, DX, DXOLD, FACT
      integer MAXIT, J

      DATA FTOL,RTOL/1.D-12,1.D-12/
      MAXIT = 100

c     initial bisection check over strain rates of interest      
      X1=hsl8
      X2=epdslmax
      FACT=3.0d+0*mu*dt
      bis=.true.

      call ksr(bis, X1, gsl, FACT, dt, epsl0, temp, tempmelt, 
     &  hsltype, hsl1,hsl2,hsl3,hsl4, hsl5, hsl6, hsl7, hsl8, 
     &  FL, DF) 

c     if FL is positive, no slip will happen     
      if (FL.ge.0.0D+0) then
        epdsl = 0.0d+0
        return
      end if
     
      call ksr(bis, X2, gsl, FACT, dt,epsl0, temp, tempmelt, 
     &  hsltype, hsl1,hsl2,hsl3,hsl4, hsl5, hsl6, hsl7, hsl8, 
     &  FH, DF) 

      IF(DABS(FL) .LT. FTOL)THEN
        epdsl=X1
        RETURN
        END IF
      IF(DABS(FH) .LT. FTOL)THEN
        epdsl=X2
        RETURN
        END IF
C
      IF((FL.GT.0.d0.AND.FH.GT.0.d0).OR.
     &  (FL.LT.0.0.AND.FH.LT.0.d0)) THEN
        WRITE(6,19)X1,FL,X2,FH
   19   FORMAT(' SOLUTION NOT BOUNDED IN EZ EVAL',4G12.5)
        epdsl = 0.0d+0
        RETURN
        END IF

c     associate high and low of strain rate to high and low of function        
      IF(FL .LT. 0.0d+0)THEN
        XL=X1
        XH=X2
      ELSE
        XH=X1
        XL=X2
      ENDIF

      epdsl = 0.5d+0*(X1+X2)
      DXOLD=DABS(X2-X1)
      DX=DXOLD
      bis = .false.
      
      call ksr(bis, epdsl, gsl, FACT, dt, epsl0, temp, tempmelt, 
     &  hsltype, hsl1,hsl2,hsl3,hsl4, hsl5, hsl6, hsl7, hsl8, 
     &  F, DF)     
C
C --- BISECT IF SOLUTION EXCEEDS LIMIT OR IF SLOW CONVERGENCE
C       OTHERWISE USE NEWTON ITERATION
C --- CONVERGENCE CHECKS ON BOTH STRAIN RATE AND NORMALIZED FUNCTION
C
      DO 10 J=1,MAXIT
C
        IF(((epdsl-XH)*DF-F)*((epdsl-XL)*DF-F) .GE. 0.0d+0 .OR. 
     *    DABS(2.0d+0*F) .GT. DABS(DXOLD*DF) ) THEN
          DXOLD=DX
          DX=0.5d+0*(XH-XL)
          epdsl=XL+2.0d+0/3.0d+0*DX
          IF(DABS(XL-epdsl)*dt .LT. RTOL .AND. 
     &       DABS(F) .LT. FTOL)RETURN
        ELSE
          DXOLD=DX
          DX=F/DF
          tempvar=epdsl
          epdsl=epdsl-DX
          IF(DABS(tempvar-epdsl)*dt .LT. RTOL .AND. 
     &       DABS(F) .LT. FTOL)RETURN
        ENDIF
C
C ---   GET FUNCTION AND SLOPE FOR NEXT ITERATION
C
        call ksr(bis, epdsl, gsl, FACT, dt,epsl0,temp, tempmelt, 
     &   hsltype, hsl1,hsl2,hsl3,hsl4, hsl5, hsl6, hsl7, hsl8, 
     &   F, DF) 

        IF(DABS(DX) .LT. RTOL .AND. DABS(F) .LT. FTOL) RETURN
C
        IF(F .LT. 0.0d+0) THEN
          XL=epdsl
        ELSE
          XH=epdsl
        ENDIF
   10 CONTINUE
C
C --- CUT TIME STEP IF NO CONVERGENCE
C
c      NFAIL=.TRUE.
      print*, 'EDOT SOLUTION DID NOT CONVERGE IN 100 STEPS IN EZ'
      
      return
      end      
      
      
c====================================================================
c====================================================================
c  Form constants to be used in state and deriv equations for epsdsl 
c--------------------------------------------------------------------

      subroutine ksr(bis, edot, sigt, fact, dt,
     &  epsl0, temp, tempmelt, hsltype,
     &  hsl1, hsl2, hsl3, hsl4, hsl5, hsl6, hsl7, hsl8, 
     &  f, df) 
            
      implicit none

c     input
      logical bis
      double precision edot, sigt
c     - parameters for loading and interaction     
      double precision fact, dt
c     - used in strength call only
      double precision epsl0, temp, tempmelt
      integer hsltype
      double precision hsl1, hsl2, hsl3, hsl4, hsl5, hsl6, hsl7, hsl8 
      
c     output
      double precision f, df

c     util - temporary var1, 2, von mises stress, dvm / dedot
      double precision sigvm, dsded
      double precision epsl
      
      epsl = epsl0 + edot*dt

      call calc_str_sl(bis, epsl, edot, temp, 
     & tempmelt,hsltype, hsl1, hsl2, hsl3, hsl4, hsl5, hsl6, hsl7, hsl8,
     & sigvm, dsded)

      f = (sigvm + fact*edot-sigt)/sigt
      
      if (bis) return
      
      df = (dsded + fact) / sigt
      
      return
      end        

c====================================================================
c====================================================================
c  Calculate deviatoric stress based on trial stress and amount 
c    of plastic deformation 
c--------------------------------------------------------------------

      subroutine calc_return_dev_stress(sigdT, mu, depsl, yssl, sigd)
      
      implicit none
      
c     input
      double precision sigdT(3,3)
      double precision mu, depsl, yssl

c     output
      double precision sigd(3,3)

c     util     
      integer i,j

      do i=1,3
        do j=1,3
        sigd(i,j) = sigdT(i,j)/(1.0d+0+3.0d+0*mu*depsl/yssl)
        end do
      end do

      return
      end


c====================================================================
c  Construct state variable array from other variables
c--------------------------------------------------------------------
      subroutine update_statevs(depsl, epdsl, energy, tempsv, nstatv,
     &    statev)     
      implicit none

c     input      
      double precision depsl, epdsl
      double precision energy, tempsv, dtime
      integer nstatv
      
c     output
      double precision statev(nstatv)      
c
      statev(1) = statev(1) + depsl
      statev(2) = epdsl
      statev(3) = energy
      statev(4) = tempsv
      
      return
      end

c====================================================================
c====================================================================
c  Just return ddsdde for lin elast umat 
c--------------------------------------------------------------------

      subroutine elastddsdde(lam, mu, ntens, ndi, ddsdde)
      implicit none
      
      !input
      integer i,j, ntens, ndi
      double precision lam, mu
      
      !output
      double precision ddsdde(ntens,ntens)

      do i=1,ntens
        do j=1,ntens
          ddsdde(i,j) = 0.0d0
        end do
      end do

      do i=1, ndi
        do j=1, ndi
          ddsdde(j,i) = lam
        end do
        ddsdde(i,i) = lam+2.0d0*mu
      end do

      do i=ndi+1, ntens
        ddsdde(i,i) = mu
      end do

      return
      end

c====================================================================
c====================================================================
c  Calculate damage based on updated state variables
c--------------------------------------------------------------------
      subroutine calc_dmg(dt, tempmelt, stress, 
     &  dmgtype, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, nstatv,
     &  statev)

      implicit none
      
      !input
      double precision dt, tempmelt
      double precision stress(6)
      integer dmgtype 
      double precision d1, d2, d3, d4, d5, d6, d7, d8, d9, d10
      integer nstatv
      
      !input/output
      double precision statev(nstatv)
      
      !util
      double precision thom, sigvm, sigm, triax, srterm, ef
      double precision temp, ep, epd
      double precision dmg, dmgflag
      

      ep = statev(1)
      epd = statev(2)
      temp = statev(4)
      dmg = statev(5)
      dmgflag = statev(6)

c--------------------------------------------------------------------
c       dmgtype = 0: Skip dmg routine 
c--------------------------------------------------------------------      
      if (dmgtype.eq.0) then
        dmg = 0.0d0
        dmgflag = 0.0d0
      else 
      
c     if melted or already damaged, then set dmg = 1
      if (temp.ge.tempmelt) then
        dmg = 1.0d0
        dmgflag = 1.0d0
      else if (dabs(dmgflag-1.0d0).le.0.1d0) then
        dmg = 1.0d0
        dmgflag = 1.0d0
      else
        thom = (temp-293d+0)/(tempmelt-293d+0)
        sigvm = dsqrt(0.5d0*((stress(1)-stress(2))**2+
     &    (stress(2)-stress(3))**2+(stress(3)-stress(1))**2)+
     &    3.0d0*(stress(4)**2+stress(5)**2+stress(6)**2))   
        sigm = (stress(1)+stress(2)+stress(3))/3.0d0 
        if (dabs(sigvm).le.1.0d-8) then !triax would be NaN if sigvm=0 
          sigvm = 1.0d-8
        end if 
        triax = sigm/sigvm
        
        if (dmgtype.eq.1) then
c--------------------------------------------------------------------
c       dmgtype = 1: Johnson-Cook dmg 
c                  - Does not use linear decay to spall stress
c                  - Same as assuming efmin = ef(triax=1.5)
c                  - Strain rate less than epd0 does not degrade 
c                  - Should be similar to Alegra implementation
c--------------------------------------------------------------------        
c       ef = (d1+d2exp(d3*triax))*(1+d4ln(epd/epd0))(1+d5*thom)
c--------------------------------------------------------------------
c       d1 = d1
c       d2 = d2
c       d3 = d3
c       d4 = d4
c       d5 = d5
c       d6 = epd0
c       d7-d10 = not used

          if (epd.le.dexp(1.0d0)*d6) then !do not evaluate sr term
            srterm = 1.0d0
          else
            srterm = 1.0d0+d4*dlog(epd/d6)
          end if 
          if (triax.le.1.5d0) then 
            ef=(d1+d2*dexp(d3*triax))*srterm*(1.0d0+thom*d5)     
          else
            ef=(d1+d2*dexp(d3*1.5d0))*srterm*(1.0d0+thom*d5) 
          end if
          dmg = dmg + epd*dt/ef
        else
c--------------------------------------------------------------------
c       dmgtype = 2: Johnson-Cook dmg 
c                    Uses a linear decay to spall stress after triax>1.5
c                    Strain rate less than epd0 does not degrade 
c                    Should be equivalent to ale3d implementation
c--------------------------------------------------------------------        

          print*, 'NO OTHER DMG TYPES IMPLEMENTED'
          dmg = dmg 
        end if !end dmg types

c       check damage flag, clean up
        if (dmg.ge.1.0d0) then
          dmg = 1.0d0
          dmgflag = 1.0d0
        else
          dmgflag = 0.0
        end if
      end if !end initially undamaged iteration         
      end if !end whether or not to skip damage (dmgtype=0)

c     put updated dmg variables back into statev
      statev(5) = dmg
      statev(6) = dmgflag

      return
      end
      
c====================================================================
c====================================================================
c  Calculate energy and pressure based volumetric strain
c--------------------------------------------------------------------

      subroutine eos(dfgrd0, dfgrd1, b0, dbdp, gam, gruna, cv,
     &   yssl, depsl, ecolds, ecold0, ecold1, ecold2, dmgflag, d10,
     &   p, energy, tempsv)       
      implicit none

c     input
      double precision dfgrd0(3,3), dfgrd1(3,3)
      double precision b0, dbdp, gam, gruna, cv
      double precision depsl, yssl
      double precision ecolds, ecold0, ecold1, ecold2
      double precision dmgflag, d10

c     input/output
      double precision p, energy

c     output
      double precision tempsv

c     util
      double precision jnew, jold, p0, xmu
      double precision determinant, locv, locdv, loce, loce0, locg

      double precision ONE
      parameter (ONE=1.0D+0)

c     determine jacobian, initialize pressure and energy
      jnew = determinant(dfgrd1)
      jold = determinant(dfgrd0)
      p0 = p
      loce0 = energy
      loce = loce0

c     set bounds on compression and tension
      locv = jnew
      locdv = (jnew-jold)
      if ( jnew .lt. 0.40 ) then
        locv = 0.40d0
        locdv = 0.d0
      else if ( jnew .gt. 1.15 ) then
        locv = 1.15d0
        locdv = 0.d0
      end if

c     calculate pressure - stiff in compression, linear in tension
      xmu = ONE/locv - ONE
      locg = gam + gruna*xmu
      if (jnew .le. ONE) then
        p = b0*xmu*(ONE+(ONE-0.5d0*locg)*xmu)/(ONE-(dbdp-ONE)*xmu)**2
      else
        p = b0*xmu
      end if

c     half step energy central difference on gruneisen contribution
      loce = loce - 0.5d0*p0*locdv
      p = (p + locg*loce ) / (ONE + 0.5d0*locg*locdv)        
      loce = loce - 0.5d0*p*locdv
        
c     If p is tensile and matl damaged, set p=pmin and e=e0
      if ((p.lt.-d10).and.(dabs(dmgflag-1.0d0).le.0.1d0)) then
        p = -d10
        loce = loce0
      end if
        
c     temperature calculation using the cold curve
      tempsv = (loce - 
     &  (ecold0 + ecold1*(xmu-ecolds) + ecold2*(xmu-ecolds)**2)) / cv
c
c     set integrated energy per unit reference volume energy for sdv
      energy = loce + jnew*yssl*depsl
      return
      end

c====================================================================
c====================================================================
c   Calculate the determinant of a 3 x 3 matrix.
c--------------------------------------------------------------------
      double precision function determinant(a)
      implicit none

      double precision a(3,3)
      double precision b1, b2, b3

      b1 = a(2,2) * a(3,3) - a(3,2) * a(2,3)
      b2 = a(3,1) * a(2,3) - a(2,1) * a(3,3)
      b3 = a(2,1) * a(3,2) - a(3,1) * a(2,2)

      determinant = a(1,1) * b1 + a(1,2) * b2 + a(1,3) * b3
      if (dabs(determinant).le.1d-10) then
        print*, 'WARNING: JUST TOOK DET LESS THAN 1E-10'
      end if

      return
      end
