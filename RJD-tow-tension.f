      SUBROUTINE VUMAT_TOW(STATEV,DDSDDE,
     2 CMNAME,
     3 NDI, NSHR,NTENS,NSTATV,PROPS,NPROPS,
     4 CELENT, SIG, TIMEVUMAT)
     
      IMPLICIT NONE

      ! Abaqus Defined Variables
      CHARACTER*80 CMNAME
      INTEGER:: NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NTENS,NSTATV,NPROPS
      INTEGER:: NDI,NSHR
      REAL*8,DIMENSION(6):: STRESS,DSTRAN,STRAN,DDSDDT,DRPLDE
      REAL*8,DIMENSION(6,6):: DDSDDE
      REAL*8,DIMENSION(NSTATV):: STATEV
      REAL*8,DIMENSION(NPROPS):: PROPS
      REAL*8,DIMENSION(2):: TIME
      REAL*8:: DTIME,TEMP,DTEMP,CELENT,PNEWDT,DRPLDT
      REAL*8:: SSE, SPD, SCD,RPL
      REAL*8:: PREDEF(1)
      REAL*8:: DPRED(1)
      REAL*8,DIMENSION(3,3):: COORDS,DROT,DFGRD0,DFGRD1
      
      
      
      REAL*8,DIMENSION(NTENS):: STRAIN, OLDSTRESS, dSTRAIN_MECH
      REAL*8,DIMENSION(NTENS):: SIG, STRAN_DEFG, DSTRAN_DEFG
      REAL*8,DIMENSION(6,6):: SS, CC
      REAL*8,DIMENSION(3,3):: DFGRD1t, eye3, DFGRD0t, MM
      REAL*8,DIMENSION(3,3):: STRANt, DSTRANt, STRESSTEMP
      REAL*8:: Constant, E, nu, JAC, MM1, MM2, MM3
      REAL*8:: E11, E22, E33, G12, G13, G23
      REAL*8:: E110, E220, E330
      REAL*8:: v12, v13, v23
      REAL*8:: GIc, sigcr11, BZ11
      REAL*8:: D11, eps01, h1, epsf1
      REAL*8:: GdissA, GdissB, GdissRatio
      REAL*8:: GIc2, sigcr22, BZ22
      REAL*8:: D22, eps02, h2, epsf2
      REAL*8:: GdissA2, GdissB2, GdissRatio2
      REAL*8:: GIc3, sigcr33, BZ33
      REAL*8:: D33, eps03, h3, epsf3
      REAL*8:: GdissA3, GdissB3, GdissRatio3
      REAL*8,DIMENSION(8):: CCMinput
      REAL*8,DIMENSION(8):: CCMoutput
      REAL*8:: Em, vm
      REAL*8:: E11c, E22c, E33c, G12c, G13c, G23c
      REAL*8:: v12c, v13c, v23c, Volf, Volf_mean
      REAL*8:: Ef1, Ef2, vf12, vf23, Gf12, Gf23
      REAL*8:: D1LIM, D2LIM, D3LIM, D12LIM, D23LIM, r, erru
      REAL*8:: Cm(6,6), Cf(6,6)	  
      REAL*8:: eps_m(6), eps_m11(6), eps_m22(6), eps_m33(6)
      REAL*8:: eps_m12(6), eps_m13(6), eps_m23(6)
      REAL*8:: Em0, vm0, SIG_Y, k1, k2, EPS_Y, Dm, Vf
      REAL*8:: Em_s, vm_s, s_m, EPS_J2, SIG_J2, Em_s_prev, vm_s_prev
      INTEGER:: I
      REAL*8:: del_EPS_11, del_EPS_22, del_EPS_12, del_EPS_23, GIc1
      REAL*8:: D12, D23, GIc12, GIc23, sigcr12, sigcr23
      REAL*8:: HASHIN2
      REAL*8:: eps012, eps023, epsf12, epsf23, BZ12, BZ23
      REAL*8:: E120, E230, sigcr12_h, sigcr23_h	  

      real*8:: TIMEVUMAT, use_SDV, sigcr22_h, eps02_h, sigcr11_h, eps01_h	   
      real*8:: eps012_h, eps023_h

       


!C------------------------------------------
!C INPUTS AND STATE VARIABLES
!C------------------------------------------

      Ef1 = 233000.0d0        !Fiber elastic properties
      Ef2 = 15000.0d0
      vf12 = 0.2d0
      vf23 = 0.3d0
      Gf12 = 8963.0d0
      Gf23 = 8963.0d0 
      Volf = PROPS(4)
      Em = PROPS(1)
	  
	  IF (Em.lt.500.0d0) THEN       !@ RJD-REMOVE THIS-THIS IS "FIX" TO ENSURE Em is NOT LOW.... CHECK WHY IN THE CYCLE JUMPING AND MATLAB CODE
	      Em = 500.0d0
	  END IF
	  
	  
	  
      Em0 = Em
      vm = PROPS(2)
      vm0 = vm
      Dm = PROPS(3)
	  
	  SIG_Y = 15.66d0   ! Not activated
	  k1 = 3308.0d0      ! Not activated
	  k2 = 56.92d0       ! Not activated
	  EPS_Y = SIG_Y/Em0
	  



! TO DO 
! * initialization of statements in the master file
! * redundant Volf, Vf, Em, Em0, vm, vm0 etc...
!

      use_SDV = PROPS(5)
      sigcr11 = PROPS(14)
      sigcr22 = PROPS(6)
      sigcr12 = PROPS(18)
      sigcr23 = PROPS(22)	  
      GIc1 = PROPS(13)     
      GIc2 = PROPS(7)             
      STATEV(139) = PROPS(6)
	  STATEV(140) = PROPS(18)
	  STATEV(141) = PROPS(22)
      GIc3 = GIc2
      sigcr33 = sigcr22
	  GIc12 = PROPS(17)
	  GIc23 = PROPS(21)
	  
	  IF (PROPS(10).gt.0.5d0) THEN  !If already in Hashin CB
	     sigcr22 = PROPS(24)        ! these 3 sigmas do not change with N once in CB
		 sigcr12 = PROPS(25)
		 sigcr23 = PROPS(26)
	  END IF
      
	  
	  
      D1LIM = 0.0001d0
      D2LIM = 0.0001d0
      D3LIM = 0.0001d0
	  D12LIM = 0.0001d0
	  D23LIM = 0.0001d0

	  IF (STATEV(100).lt.0.5d0) THEN
           STATEV(9) = Em0            !Em: Young's modulus
           STATEV(10) = vm0           !vm: Poisson's ratio
           STATEV(44) = Dm           !Dm: Stiffness fraction in matrix
           STATEV(33) = NOEL
           STATEV(32) = Volf
           STATEV(106) = PROPS(10)         ! 22 crackband flag
           STATEV(101) = PROPS(11)         ! 11 crackband flag
	  	  !write(*,*) 'I m inside', STATEV(7), STATEV(9)
	  END IF	  

      STATEV(137) = Volf
      Vf = STATEV(137)

      !write(*,*) sigcr22, use_SDV, PROPS(10)
	 
!C----------------------------------------------------------
!C	  LIST AND DESCRIPTION OF STATE VARIABLES USED
!C----------------------------------------------------------


      !STATEV 1 - 6 : Stiffness from the CCM
      !STATEV 7 - 8 : Em and vm from previous increment
      !STATEV 9 - 10: Em and vm at the end of increment
      !STATEV 11: EPS_J in matrix
      !STATEV 12: EQ STRESS in matrix
      !STATEV 13: Nonlinear flag in the matrix
	  !STATEV 14: Element number

!C      STATEV(101) = Crack Flag 1-direction (0,1,2)
!C                  0: no failure, 1: softening
!C                  2: complete failure
!C      STATEV(102) = D11, Stiffness reduction factor 1-direction
!C      STATEV(103) = Strain in the 1-direction
!C      STATEV(104) = Stress in the 1-direction
!C      STATEV(105) = Energy dissipated 1-direction
!C      STATEV(106) = Crack flag 2-direction (0,1,2)
!C      STATEV(107) = D22, Stiffness reduction factor 2-direction
!C      STATEV(108) = Strain in the 2-direction
!C      STATEV(109) = Stress in the 2-direction
!C      STATEV(110) = Energy dissipated 2-direction
!C      STATEV(111) = Crack Flag 3-direction (0,1,2)
!C      STATEV(112) = D33, Stiffness reduction factor 3-direction
!C      STATEV(113) = Strain in the 3-direction
!C      STATEV(114) = Stress in the 3-direction
!C      STATEV(115) = Energy dissipated 3-direction
!C      STATEV(116) = Effective E in the 1-direction
!C      STATEV(117) = Effective E in the 2-direction
!C      STATEV(118) = Effective E in the 3-direction
!C      STATEV(120) = ELEMENT DELETION FLAG
!C      
!C      STATEV(121) = eps01    ! Crack strain at doorstep of crackband, 1-dir
!C      STATEV(122) = eps02    
!C      STATEV(123) = eps03    
!C      STATEV(124) = CELENT   ! Characteristic element length
!C      STATEV(125) = epsf1    ! Failure crack strain, tail of the softening 
!C      STATEV(126) = epsf2    
!C      STATEV(127) = epsf3         
!C      STATEV(131) = E11c     ! Composite stiffness (undamaged) 
!C      STATEV(132) = E22c     
!C      STATEV(133) = v12c     
!C      STATEV(134) = v23c     
!C      STATEV(135) = G12c     
!C      STATEV(136) = G23c
!C      STATEV(137) = Vf (random)   ! Vf computed from uniform distribution
!C      STATEV(138) = sig11_cr
!C      STATEV(139) = sig22_cr     
!
!C
!C      !!!!Applied strains during the loading step
!C      STATEV(151) = STRAIN(1)
!C      STATEV(152) = STRAIN(2)
!C      STATEV(153) = STRAIN(3)
!C      STATEV(154) = STRAIN(4)
!C      STATEV(155) = STRAIN(5)
!C      STATEV(156) = STRAIN(6)
!C
!C      STATEV(161) = STRESS(1)        ! Mechanical Strain component 1-direction, etc...
!C      STATEV(162) = STRESS(2)        
!C      STATEV(163) = STRESS(3)
!C      STATEV(164) = STRESS(4)
!C      STATEV(165) = STRESS(5)
!C      STATEV(166) = STRESS(6)
!C
!C      STATEV(193) = D11*E11c 
!C      STATEV(194) = D22*E22c      !RJD-CHECK-FOR-TRANSVERSE DAMGE FOR PARIS LAW     
!C      STATEV(199) = check flags in certain portions of the code
!
!
!     
      !Store Stress from the previous increment      
      STRESS(1) = STATEV(161)
      STRESS(2) = STATEV(162)
      STRESS(3) = STATEV(163)
      STRESS(4) = STATEV(164)
      STRESS(5) = STATEV(165)
      STRESS(6) = STATEV(166)

      !Applied strains during the loading step
      STRAIN(1) = STATEV(151)
      STRAIN(2) = STATEV(152)
      STRAIN(3) = STATEV(153)
      STRAIN(4) = STATEV(154)
      STRAIN(5) = STATEV(155)
      STRAIN(6) = STATEV(156)
	  
	  del_EPS_11 = STATEV(151) - STATEV(103)  
      del_EPS_22 = STATEV(152) - STATEV(108)
	  del_EPS_12 = abs(STATEV(154)) - abs(STATEV(170))        !@ fill this
	  del_EPS_23 = abs(STATEV(156)) - abs(STATEV(175))        !@ fill this
	  STATEV(195) = del_EPS_11
      STATEV(196) = del_EPS_22
	  STATEV(197) = del_EPS_12
      STATEV(198) = del_EPS_23	  
	  
      !-------------------------------------------------------
      !Compute effective composite stiffness
      !-------------------------------------------------------
	  
      CCMinput(1) = Ef1
      CCMinput(2) = Ef2
      CCMinput(3) = Gf12
      CCMinput(4) = Gf23
      CCMinput(5) = vf12
      CCMinput(6) = Volf
      CCMinput(7) = Em
      CCMinput(8) = vm
      
      CALL CCM_CALC_tow(CCMinput,CCMoutput,ntens,nstatv,statev)
	  
      E11c = CCMoutput(1)
      E22c = CCMoutput(2)
      v12c = CCMoutput(3)
      v23c = CCMoutput(4)
      G12c = CCMoutput(5)
      G23c = CCMoutput(6)
      E33c = E22c
      v13c = v12c
      G13c = G12c
        
      STATEV(131) = E11c
      STATEV(132) = E22c
      STATEV(133) = v12c
      STATEV(134) = v23c
      STATEV(135) = G12c
      STATEV(136) = G23c      

      !----------------------------------------------------

      !----------------------------------------------------


!Matrix secant modulus and secant poisson's ratio from previous inc
      
  	  !IF (STATEV(100).lt.0.5d0) THEN
	    !STATEV(9) = PROPS(1)
		!STATEV(10) = PROPS(2)
      !END IF
	  
	  STATEV(9) = Em0      !RJD-assumption is E_sec = Em0
	  Em_s_prev = STATEV(9)
      STATEV(10) = vm0     !RJD-assumption is v_sec = vm0
	  vm_s_prev = STATEV(10)
      STATEV(7) = Em_s_prev
	  STATEV(8) = vm_s_prev
      


      !Calculate stiffness terms for matrix and fiber
      CALL calc_stiffness_iso(Em_s_prev, vm_s_prev, Cm)
	                      STATEV(14) = Cm(1,1)  !RJD-CHECK
	                      STATEV(15) = Cm(2,2)  !RJD-CHECK
	                      STATEV(16) = Cm(3,3)  !RJD-CHECK
	                      STATEV(17) = Cm(4,4)  !RJD-CHECK
	                      STATEV(18) = Cm(5,5)  !RJD-CHECK
	                      STATEV(19) = Cm(6,6)  !RJD-CHECK
	  
      CALL calc_stiffness_transiso(Ef1, Ef2, vf12, 
     1                                 vf23, Gf12, Gf23, Cf)
	                      STATEV(20) = Cf(1,1)  !RJD-CHECK
	                      STATEV(21) = Cf(2,2)  !RJD-CHECK
	                      STATEV(22) = Cf(3,3)  !RJD-CHECK
	                      STATEV(23) = Cf(4,4)  !RJD-CHECK
	                      STATEV(24) = Cf(5,5)  !RJD-CHECK
	                      STATEV(25) = Cf(6,6)  !RJD-CHECK
	 
	  
      !Matrix constituent strains due to composite strain epsilon_11	  
	  eps_m11 = 0.0d0
      eps_m11(1) = STRAIN(1)
      eps_m11(2) = -1.0d0*STATEV(8)*STRAIN(1)
      eps_m11(3) = -1.0d0*STATEV(8)*STRAIN(1)

      !Matrix constituent strains due to composite strain epsilon_22
      eps_m22 = 0.0d0
      eps_m22(2) = STRAIN(2)*(1.0d0/(1.0d0 - DSQRT(Vf)))
     1                                  /(1.0d0 + (Cm(2,2)/Cf(2,2)))
      eps_m22(1) = -1.0d0*STATEV(8)*eps_m22(2)
      eps_m22(3) = -1.0d0*STATEV(8)*eps_m22(2)

      !Matrix constituent strains due to composite strain epsilon_33
      eps_m33 = 0.0d0
      eps_m33(3) = STRAIN(3)*(1.0d0/(1.0d0 - DSQRT(Vf)))
     1                                    /(1.0d0 + (Cm(3,3)/Cf(3,3)))
      eps_m33(1) = -1.0d0*STATEV(8)*eps_m33(3)
      eps_m33(2) = -1.0d0*STATEV(8)*eps_m33(3)

      !Matrix constituent strains due to composite shear strains
      eps_m12 = 0.0d0
      eps_m13 = 0.0d0
      eps_m23 = 0.0d0
      eps_m12(4) = STRAIN(4)*(1.0d0/(1.0d0 - DSQRT(Vf)))
     1                                  /(1.0d0 + (Cm(1,2)/Cf(1,2)))
      eps_m13(5) = STRAIN(5)*(1.0d0/(1.0d0 - DSQRT(Vf)))
     1                                  /(1.0d0 + (Cm(1,3)/Cf(1,3)))
      eps_m23(6) = STRAIN(6)*(1.0d0/(1.0d0 - DSQRT(Vf)))
     1                                  /(1.0d0 + (Cm(2,3)/Cf(2,3)))



      !Add constituent matrix strains using superposition principle
      eps_m = 0.0d0
      DO I=1,6
        eps_m(I) = eps_m11(I) + eps_m22(I) + eps_m33(I) + 
     1             eps_m12(I) + eps_m13(I) + eps_m23(I)
      END DO
	  
	                      STATEV(26) = eps_m(1)  !RJD-CHECK
	                      STATEV(27) = eps_m(2)  !RJD-CHECK
	                      STATEV(28) = eps_m(3)  !RJD-CHECK
	                      STATEV(29) = eps_m(4)  !RJD-CHECK
	                      STATEV(30) = eps_m(5)  !RJD-CHECK
	                      STATEV(31) = eps_m(6)  !RJD-CHECK


      !Calculate equivalent matrix strain from matrix constituent strains
      CALL calc_eqstrain(eps_m, EPS_J2)     
      !----------------------------------------------------------------
      !Equivalent matrix secant modulus and Poisson's ratio calculation
      !----------------------------------------------------------------
      !Flag for eq. strain in linear or nonlinear part in matrix curve

      !Compute equivalent stress
      !IF (EPS_J2.LE.(999999999.0d0*EPS_Y)) THEN
      SIG_J2 = Em0*EPS_J2      !Secant will be equal to initial modulus
      !ELSEIF (EPS_J2.GT.(999999999.0d0*EPS_Y)) THEN
      !STATEV(13) = 1.0
      !SIG_J2 = SIG_Y - (k1/k2)*(EXP(-k2*EPS_J2) - !EXP(-k2*SIG_Y/Em0))
      !ENDIF

      
      !Compute secant modulus and secant Poisson's ratio
      Em_s = SIG_J2/EPS_J2
      vm_s = 0.5d0 + (Em_s/Em0)*(vm0 - 0.5d0)

           STATEV(9) = Em_s
           STATEV(10) = vm_s
           STATEV(11) = EPS_J2
           STATEV(12) = SIG_J2
      !----------------------------------------------------------------
 
      !IF (KINC.EQ.1) THEN
	  IF (STATEV(100).lt.0.5d0) THEN
            Em_s = Em0
            vm_s = vm0
      END IF


!C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!C---------------------------------------------------------
!C  WHEN THE MATERIAL FIRST ENTERS CRACKBAND
!C---------------------------------------------------------
      !C     ! For 1-direction
      !When material first enters softening region
	  
      IF (use_SDV.lt.0.5d0) THEN   !!! RJD-EDIT
            IF ((STRESS(1).ge.sigcr11).and.(STATEV(101).eq.0.0D0)) THEN
               STATEV(101) = 1.0D0
               STATEV(116) = STRESS(1)/STRAIN(1)
            END IF 
      END IF
         
      IF (use_SDV.gt.0.5d0) THEN    !!! RJD-EDIT
          IF ((STRESS(1).ge.sigcr11).and.(PROPS(11).lt.0.5d0)) THEN
             STATEV(101) = 1.0D0
             STATEV(116) = STRESS(1)/STRAIN(1)
          END IF
          IF (PROPS(11).gt.0.5d0) THEN
             STATEV(101) = PROPS(11)      !Crackband flag type
             STATEV(116) = PROPS(12)      !Effective E in 1-direction
          END IF
      END IF
         
         
       eps01 = sigcr11/STATEV(116)
       STATEV(121) = eps01
 
         
          ! Crackband parameters stored once at crack initiation
          h1 = CELENT
          STATEV(124) = CELENT 
          epsf1 = (2.0D0*GIc1)/(h1*sigcr11)
          STATEV(125) = epsf1  
          !gg1 = GIc1/h1
          !STATEV(128) = gg1
               
          ! Check Bazant limit
          BZ11 = 2.0D0*E11c*GIc1/(sigcr11**2)  
          IF (h1.ge.0.95*BZ11) THEN
            write(*,*) 'Characteristic length of element is too large!!'
            write(*,*) 'Reduce mesh size'
            write(*,*) 'h1', h1
            write(*,*) 'BZ11', BZ11
		    write(*,*) 'GIc1', GIc1
		    write(*,*) 'sigcr11', sigcr11
            write(*,*) 'E11c', E11c
            write(*,*) 'Em', Em			
            CALL XPLB_EXIT
          END IF


	  
	  
	  
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	  
      ! Evaluate hashin criterion
      IF (STATEV(106).lt.0.5d0) THEN	  
	  HASHIN2 = (STRESS(2)/sigcr22)**2 + (STRESS(4)/sigcr12)**2
     1	                              + (STRESS(6)/sigcr23)**2
         IF (HASHIN2.ge.1.0d0) THEN
		    sigcr22 = STRESS(2)
			sigcr12 = STRESS(4)
			sigcr23 = STRESS(6)
			STATEV(184) = sigcr22
			STATEV(185) = sigcr12
			STATEV(186) = sigcr23
		 END IF

	  END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC





!C-------------------------------------------
!C     ! For 2-direction (Hashin criterion) 

      !When material first enters softening region	  
	  IF (use_SDV.lt.0.5d0) THEN
	     IF ((HASHIN2.ge.1.0d0).and.(STATEV(106).eq.0.0D0)) THEN
              STATEV(106) = 1.0D0
              STATEV(117) = STRESS(2)/STRAIN(2) 
			  STATEV(178) = DABS(STRESS(4)/STRAIN(4))  !@
			  STATEV(179) = DABS(STRESS(6)/STRAIN(6))  !@
         END IF
	  END IF
	  
	  IF (use_SDV.gt.0.5d0) THEN
	  	 IF ((HASHIN2.ge.1.0d0).and.(PROPS(10).lt.0.5d0)) THEN
	          STATEV(106) = 1.0D0
              STATEV(117) = STRESS(2)/STRAIN(2)
			  STATEV(178) = DABS(STRESS(4)/STRAIN(4)) !@
			  STATEV(179) = DABS(STRESS(6)/STRAIN(6)) !@
		 END IF
		 IF (PROPS(10).gt.0.5d0) THEN
		     STATEV(106) = PROPS(10)
             STATEV(117) = PROPS(8)
			 STATEV(178) = PROPS(16) !@
			 STATEV(179) = PROPS(20) !@
         END IF 
	  END IF
	  
	  
      eps02 = sigcr22/STATEV(117)
	  eps012 = sigcr12/STATEV(178)  !@
	  eps023 = sigcr23/STATEV(179)  !@
      STATEV(122) = eps02
	  STATEV(180) = eps012          !@
	  STATEV(181) = eps023          !@

	  
         ! Crackband parameters stored once at crack initiation
         h2 = CELENT
         STATEV(124) = CELENT 
         epsf2 = (2.0D0*GIc2)/(h2*sigcr22)
		 epsf12 = (2.0D0*GIc12)/(h2*sigcr12)  !@
		 epsf23 = (2.0D0*GIc23)/(h2*sigcr23)  !@	 	 
         STATEV(126) = epsf2 
         STATEV(182) = epsf12        !@
         STATEV(183) = epsf23 		 !@ 
         !gg2 = GIc2/h2
         !STATEV(129) = gg2
		 
	        
         ! Check Bazant limit
         !E22 = STATEV(32)                  !Equal to E22c
         BZ22 = 2.0D0*E22c*GIc2/(sigcr22**2)  
         BZ12 = 2.0D0*G12c*GIc12/(sigcr12**2)     !@
         BZ23 = 2.0D0*G23c*GIc23/(sigcr23**2)  	  !@	  
         IF ((h2.ge.(0.95*BZ22)).or.(h2.ge.(0.95*BZ12)).or.
     1                               (h2.ge.(0.95*BZ23))) THEN
           write(*,*) 'Characteristic length of element is too large!!'
           write(*,*) 'Reduce mesh size'
		   write(*,*) '********************'
           write(*,*) 'h2', h2
           write(*,*) 'BZ22', BZ22
		   write(*,*) 'BZ12', BZ12                
		   write(*,*) 'BZ23', BZ23                
		   write(*,*) '********************'
		   write(*,*) 'E22c', E22c
		   write(*,*) 'GIc2', GIc2
		   write(*,*) 'sigcr22', sigcr22
		   write(*,*) '********************'
		   write(*,*) 'G12c', G12c
		   write(*,*) 'GIc12', GIc12
		   write(*,*) 'sigcr12', sigcr12
		   write(*,*) '********************'
		   write(*,*) 'G23c', G23c
		   write(*,*) 'GIc23', GIc23
		   write(*,*) 'sigcr23', sigcr23
		   write(*,*) '********************'
		   write(*,*) 'Em', Em
           CALL XPLB_EXIT
         END IF



!C-------------------------------------------
!C     ! For 3-direction

      !When material first enters softening region
      IF ((STRESS(3).ge.1000.D0*sigcr33).and.(STATEV(111).eq.0.0D0)) 
     1 THEN
         STATEV(111) = 1.0D0
         STATEV(118) = STRESS(3)/STRAIN(3)
	  
         ! Crackband parameters stored once at crack initiation
         eps03 = sigcr33/STATEV(118)
         STATEV(123) = eps03
	  
         h3 = CELENT
         STATEV(124) = CELENT
	  
         epsf3 = (2.0D0*GIc3)/(h3*sigcr33)
         STATEV(127) = epsf3
	  
         !gg3 = GIc3/h3
         !STATEV(130) = gg3
	  
         ! Check Bazant limit
         !E33 = STATEV(32)                !Equal to E22c
         BZ33 = 2.0D0*E33c*GIc3/(sigcr33**2)  !RJD-artificially large for now
         IF (h3.ge.(0.8*BZ33)) THEN
             write(*,*) 'Characteristic length of element large!!'
             write(*,*) 'Reduce mesh size'
             write(*,*) 'h3', h3
             write(*,*) 'BZ33', BZ33
             CALL XPLB_EXIT
         END IF
	  
      END IF

      
      
      
!C-----------------------------------------------------------	  
!C-----------------------------------------------------------
!C STIFFNESS DEGRADATION FOR CRACKED/DAMAGED SOLID
!C-----------------------------------------------------------
!C-----------------------------------------------------------
      
	! For 1-direction ****************
      IF ((use_SDV).lt.0.5d0) THEN  ! RJD-EDIT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            IF (STATEV(101).eq.0.0D0) THEN   !If not yet in crackband
               D11 = 1.0D0
               STATEV(102) = D11
        ELSEIF ((STATEV(101).eq.1.0D0).and.(STRAIN(1).gt.STATEV(103))) 
     1           THEN
         E110 = STATEV(116)
         eps01 = STATEV(121)
         epsf1 = STATEV(125)
         D11=(sigcr11/E110)*(1.0D0/(epsf1-eps01))*((epsf1/STRAIN(1))-1)
          STATEV(102) = D11
                  IF (D11.lt.D1LIM) THEN
                    STATEV(102) = D1LIM
                    STATEV(101) = 2.0D0
                    STATEV(120) = 0.0D0   ! Delete element
                  ELSE IF (STATEV(102).gt.1.0D0) THEN
                    STATEV(102) = 1.0D0
                  END IF
        ELSEIF ((STATEV(101).eq.1.0D0).and.(STRAIN(1).le.STATEV(103))) 
     1         THEN
               D11 = STATEV(102)
            ELSEIF (STATEV(102).lt.D1LIM) THEN
               STATEV(102) = D1LIM
               STATEV(101) = 2.0D0
               STATEV(120) = 0.0d0
            END IF
           D11 = STATEV(102)
           END IF
     
     
           
           IF ((use_SDV).ge.0.5d0) THEN  !RJD-EDIT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            !write(*,*) 'This fork'
            IF (STATEV(101).eq.0.0D0) THEN   !If not yet in crackband
               D11 = 1.0D0
               STATEV(102) = D11
               STATEV(199) = 10.0d0
            ELSE IF ((STATEV(101).eq.1.0d0).and.(PROPS(11).lt.0.5d0)
     1              .and.(del_EPS_11.ge.0.0d0)) THEN                                   !First time entering crackband in this simulation
          STATEV(101) = 1.0d0
          E110 = STATEV(116)
          eps01 = STATEV(121)
          epsf1 = STATEV(125)
         D11=(sigcr11/E110)*(1.0D0/(epsf1-eps01))*((epsf1/STRAIN(1))-1)
          STATEV(199) = 20.0d0		  
          STATEV(102) = D11
              IF (D11.lt.D1LIM) THEN
                  D11 = D1LIM
                  STATEV(102) = D11
                  STATEV(101) = 2.0d0
                  STATEV(120) = 0.0d0
              ELSE IF (D11.gt.1.0d0) THEN
                  D11 = 1.0d0
                  STATEV(102) = D11
              END IF
        ELSE IF ((STATEV(101).eq.1.0d0).and.(PROPS(11).gt.0.5d0).and.
     1  (del_EPS_11.gt.0.0d0)) THEN                                   !Already had entered crackband in the previous simulation
          E110 = STATEV(116)
          eps01 = STATEV(121)
          epsf1 = STATEV(125)
          sigcr11_h = PROPS(15)*PROPS(12)*sigcr11*epsf1/(sigcr11  
     1               + PROPS(15)*PROPS(12)*(epsf1-eps01))   	 
          eps01_h = sigcr11_h/(PROPS(15)*PROPS(12))
          STATEV(199) = 30.0d0
          IF (STRAIN(1).ge.eps01_h) THEN             
              D11=(sigcr11_h/E110)*(1.0D0/(epsf1-eps01_h))*
     1                                         ((epsf1/STRAIN(1))-1)   
              STATEV(102) = D11
              STATEV(199) = 40.0d0			  
              IF (D11.le.D1LIM) THEN
                  D11 = D1LIM
                  STATEV(102) = D11
              END IF
          ELSE IF (STRAIN(1).lt.eps01_h) THEN
              D11 = PROPS(15)
              STATEV(102) = D11
              STATEV(199) = 50.0d0
          END IF
             
             IF (STATEV(102).lt.D1LIM) THEN
                STATEV(102) = D1LIM
                STATEV(101) = 2.0D0
                STATEV(120) = 0.0D0   ! Delete element
                STATEV(199) = 60.0d0
             ELSE IF (STATEV(102).gt.1.0D0) THEN
                STATEV(102) = 1.0D0
                STATEV(199) = 70.0d0
             END IF
       ELSEIF ((STATEV(101).eq.1.0D0).and.(del_EPS_11.lt.0.0d0)) 
     1  THEN
              D11 = STATEV(102)
              STATEV(199) = 80.0d0
            END IF
            D11 = STATEV(102)
           
           END IF !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!C-------------
      !For 2-direction   *************************
	  IF ((use_SDV).lt.0.5d0) THEN  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       IF (STATEV(106).eq.0.0D0) THEN   !If not yet in crackband
          D22 = 1.0D0
		  STATEV(107) = D22
       ELSEIF ((STATEV(106).eq.1.0D0).and.(STRAIN(2).gt.STATEV(108))) 
     1  THEN
         E220 = STATEV(117)
         eps02 = STATEV(122)
         epsf2 = STATEV(126)
         D22=(sigcr22/E220)*(1.0D0/(epsf2-eps02))*((epsf2/STRAIN(2))-1)
		  !write(*,*) 'D22 ==', D22
		  STATEV(107) = D22
             IF (D22.lt.D2LIM) THEN
               STATEV(107) = D2LIM
               STATEV(106) = 2.0D0
               STATEV(120) = 0.0D0   ! Delete element
             ELSE IF (STATEV(107).gt.1.0D0) THEN
               STATEV(107) = 1.0D0
             END IF
       ELSEIF ((STATEV(106).eq.1.0D0).and.(STRAIN(2).le.STATEV(108))) 
     1  THEN
          D22 = STATEV(107)
       ELSEIF (STATEV(107).lt.D2LIM) THEN
          STATEV(107) = D2LIM
          STATEV(106) = 2.0D0
		  STATEV(120) = 0.0d0
       END IF
      D22 = STATEV(107)
      END IF


	  
	  IF ((use_SDV).ge.0.5d0) THEN  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       !write(*,*) 'This fork'
	   IF (STATEV(106).eq.0.0D0) THEN   !If not yet in crackband
          D22 = 1.0D0
		  STATEV(107) = D22
		  STATEV(199) = 10.0d0
!       ELSE IF ((STATEV(106).eq.1.0d0).and.(PROPS(10).lt.0.5d0)) THEN		  
       ELSE IF ((STATEV(106).eq.1.0d0).and.(PROPS(10).lt.0.5d0).and.
     1              (del_EPS_22.ge.0.0d0)) THEN                                   !First time entering crackband in this simulation
         STATEV(106) = 1.0d0
         E220 = STATEV(117)
         eps02 = STATEV(122)
         epsf2 = STATEV(126)
         D22=(sigcr22/E220)*(1.0D0/(epsf2-eps02))*((epsf2/STRAIN(2))-1)
         STATEV(199) = 20.0d0  
         STATEV(107) = D22
		      IF (D22.lt.D2LIM) THEN
			      D22 = D2LIM
				STATEV(107) = D22
				STATEV(106) = 2.0d0
				STATEV(120) = 0.0d0
			  ELSE IF (D22.gt.1.0d0) THEN
			      D22 = 1.0d0
				  STATEV(107) = D22
			  END IF
       !ELSE IF ((STATEV(106).eq.1.0d0).and.(PROPS(10).gt.0.5d0)) THEN			  
       ELSE IF ((STATEV(106).eq.1.0d0).and.(PROPS(10).gt.0.5d0).and.
     1  (del_EPS_22.gt.0.0d0)) THEN                                   !Already had entered crackband in the previous simulation
		  E220 = STATEV(117)
          eps02 = STATEV(122)
          epsf2 = STATEV(126)
		  sigcr22_h = PROPS(9)*PROPS(8)*sigcr22*epsf2/(sigcr22  
     1               + PROPS(9)*PROPS(8)*(epsf2-eps02))   	 
		  eps02_h = sigcr22_h/(PROPS(9)*PROPS(8))
		  STATEV(199) = 30.0d0
		  !write(*,*) 'Inner 3 fork', sigcr22_h, eps02_h, STRAIN(2)
		  IF (STRAIN(2).ge.eps02_h) THEN             
			  D22=(sigcr22_h/E220)*(1.0D0/(epsf2-eps02_h))*
     1                                         ((epsf2/STRAIN(2))-1)   
		      STATEV(107) = D22
			  STATEV(199) = 40.0d0
			  IF (D22.le.D2LIM) THEN
			      D22 = D2LIM
				  STATEV(107) = D22
			  END IF
		  ELSE IF (STRAIN(2).lt.eps02_h) THEN
		      D22 = PROPS(9)
		      STATEV(107) = D22
			  STATEV(199) = 50.0d0
	      END IF
             
			 IF (STATEV(107).lt.D2LIM) THEN
               STATEV(107) = D2LIM
               STATEV(106) = 2.0D0
               STATEV(120) = 0.0D0   ! Delete element
			   STATEV(199) = 60.0d0
             ELSE IF (STATEV(107).gt.1.0D0) THEN
               STATEV(107) = 1.0D0
			   STATEV(199) = 70.0d0
             END IF
       ELSEIF ((STATEV(106).eq.1.0D0).and.(del_EPS_22.lt.0.0d0)) 
     1  THEN
         D22 = STATEV(107)
	     STATEV(199) = 80.0d0
       END IF
      D22 = STATEV(107)
      
	  END IF !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	  
	  





!C-------------
      !For 12-direction   *************************
	  IF ((use_SDV).lt.0.5d0) THEN  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       IF (STATEV(106).eq.0.0D0) THEN   !If not yet in crackband
          D12 = 1.0D0
		  STATEV(169) = D12
       ELSEIF ((STATEV(106).eq.1.0D0).and.(STRAIN(4).gt.STATEV(170))) 
     1  THEN
          E120 = STATEV(178)
          eps012 = STATEV(180)
          epsf12 = STATEV(182)
      D12=(sigcr12/E120)*(1.0D0/(epsf12-eps012))*((epsf12/STRAIN(4))-1)
		  !write(*,*) 'D12 ==', D12
		  STATEV(169) = D12
             IF (D12.lt.D12LIM) THEN
               STATEV(169) = D12LIM
               STATEV(106) = 2.0D0
               STATEV(120) = 0.0D0   ! Delete element
             ELSE IF (STATEV(169).gt.1.0D0) THEN
               STATEV(169) = 1.0D0
             END IF
       ELSEIF ((STATEV(106).eq.1.0D0).and.(STRAIN(4).le.STATEV(170))) 
     1  THEN
          D12 = STATEV(169)
       ELSEIF (STATEV(169).lt.D12LIM) THEN
          STATEV(169) = D12LIM
          STATEV(106) = 2.0D0
		  STATEV(120) = 0.0d0
       END IF
      D12 = STATEV(169)
      END IF


	  
	  IF ((use_SDV).ge.0.5d0) THEN  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       !write(*,*) 'This fork'
	   IF (STATEV(106).eq.0.0D0) THEN   !If not yet in crackband
          D12 = 1.0D0
		  STATEV(169) = D12
		  STATEV(199) = 10.0d0
!       ELSE IF ((STATEV(106).eq.1.0d0).and.(PROPS(10).lt.0.5d0)) THEN		  
       ELSE IF ((STATEV(106).eq.1.0d0).and.(PROPS(10).lt.0.5d0).and.
     1              (del_EPS_12.ge.0.0d0)) THEN                                   !First time entering crackband in this simulation
        STATEV(106) = 1.0d0
        E120 = STATEV(178)
        eps012 = STATEV(180)
        epsf12 = STATEV(182)
      D12=(sigcr12/E120)*(1.0D0/(epsf12-eps012))*((epsf12/STRAIN(4))-1)
        STATEV(199) = 20.0d0  
        STATEV(169) = D12
		    IF (D12.lt.D12LIM) THEN
		    D12 = D12LIM
			STATEV(169) = D12
			STATEV(106) = 2.0d0
			STATEV(120) = 0.0d0
		  ELSE IF (D12.gt.1.0d0) THEN
		      D12 = 1.0d0
			  STATEV(169) = D12
		  END IF
       !ELSE IF ((STATEV(106).eq.1.0d0).and.(PROPS(10).gt.0.5d0)) THEN			  
       ELSE IF ((STATEV(106).eq.1.0d0).and.(PROPS(10).gt.0.5d0).and.
     1  (del_EPS_12.gt.0.0d0)) THEN                                   !Already had entered crackband in the previous simulation
		  E120 = STATEV(178)
          eps012 = STATEV(180)
          epsf12 = STATEV(182)
		  sigcr12_h = PROPS(19)*PROPS(16)*sigcr12*epsf12/(sigcr12  
     1               + PROPS(19)*PROPS(16)*(epsf12-eps012))   	 
		  eps012_h = sigcr12_h/(PROPS(19)*PROPS(16))
		  STATEV(199) = 30.0d0
		  !write(*,*) 'Inner 3 fork', sigcr12_h, eps012_h, STRAIN(4)
		  IF (STRAIN(4).ge.eps012_h) THEN             
			  D12=(sigcr12_h/E120)*(1.0D0/(epsf12-eps012_h))*
     1                                         ((epsf12/STRAIN(4))-1)   
		      STATEV(169) = D12
			  STATEV(199) = 40.0d0
			  IF (D12.le.D12LIM) THEN
			      D12 = D12LIM
				  STATEV(169) = D12
			  END IF
		  ELSE IF (STRAIN(4).lt.eps012_h) THEN
		      D12 = PROPS(19)
		      STATEV(169) = D12
			  STATEV(199) = 50.0d0
	      END IF
             
			 IF (STATEV(169).lt.D12LIM) THEN
               STATEV(169) = D12LIM
               STATEV(106) = 2.0D0
               STATEV(120) = 0.0D0   ! Delete element
			   STATEV(199) = 60.0d0
             ELSE IF (STATEV(169).gt.1.0D0) THEN
               STATEV(169) = 1.0D0
			   STATEV(199) = 70.0d0
             END IF
       ELSEIF ((STATEV(106).eq.1.0D0).and.(del_EPS_12.lt.0.0d0)) THEN 
              D12 = STATEV(169)
	          STATEV(199) = 80.0d0
       END IF
       D12 = STATEV(169)
      
	  END IF !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	  







!C-------------
      !For 12-direction   *************************
	  IF ((use_SDV).lt.0.5d0) THEN  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       IF (STATEV(106).eq.0.0D0) THEN   !If not yet in crackband
          D23 = 1.0D0
		  STATEV(174) = D23
       ELSEIF ((STATEV(106).eq.1.0D0).and.(STRAIN(6).gt.STATEV(175))) 
     1  THEN
          E230 = STATEV(178)
          eps023 = STATEV(181)
          epsf23 = STATEV(183)
      D23=(sigcr23/E230)*(1.0D0/(epsf23-eps023))*((epsf23/STRAIN(6))-1)
		  !write(*,*) 'D23 ==', D23
		  STATEV(174) = D23
             IF (D23.lt.D23LIM) THEN
               STATEV(174) = D23LIM
               STATEV(106) = 2.0D0
               STATEV(120) = 0.0D0   ! Delete element
             ELSE IF (STATEV(174).gt.1.0D0) THEN
               STATEV(174) = 1.0D0
             END IF
       ELSEIF ((STATEV(106).eq.1.0D0).and.(STRAIN(6).le.STATEV(175))) 
     1  THEN
          D23 = STATEV(174)
       ELSEIF (STATEV(174).lt.D23LIM) THEN
          STATEV(174) = D23LIM
          STATEV(106) = 2.0D0
		  STATEV(120) = 0.0d0
       END IF
      D23 = STATEV(174)
      END IF


	  
	  IF ((use_SDV).ge.0.5d0) THEN  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       !write(*,*) 'This fork'
	   IF (STATEV(106).eq.0.0D0) THEN   !If not yet in crackband
          D23 = 1.0D0
		  STATEV(174) = D23
		  STATEV(199) = 10.0d0
!       ELSE IF ((STATEV(106).eq.1.0d0).and.(PROPS(10).lt.0.5d0)) THEN		  
       ELSE IF ((STATEV(106).eq.1.0d0).and.(PROPS(10).lt.0.5d0).and.
     1              (del_EPS_23.ge.0.0d0)) THEN                                   !First time entering crackband in this simulation
        STATEV(106) = 1.0d0
        E230 = STATEV(178)
        eps023 = STATEV(181)
        epsf23 = STATEV(183)
      D23=(sigcr23/E230)*(1.0D0/(epsf23-eps023))*((epsf23/STRAIN(6))-1)
        STATEV(199) = 20.0d0  
        STATEV(174) = D23
		    IF (D23.lt.D23LIM) THEN
		    D23 = D23LIM
			STATEV(174) = D23
			STATEV(106) = 2.0d0
			STATEV(120) = 0.0d0
		  ELSE IF (D23.gt.1.0d0) THEN
		      D23 = 1.0d0
			  STATEV(174) = D23
		  END IF
       !ELSE IF ((STATEV(106).eq.1.0d0).and.(PROPS(10).gt.0.5d0)) THEN			  
       ELSE IF ((STATEV(106).eq.1.0d0).and.(PROPS(10).gt.0.5d0).and.
     1  (del_EPS_23.gt.0.0d0)) THEN                                   !Already had entered crackband in the previous simulation
		  E230 = STATEV(178)
          eps023 = STATEV(181)
          epsf23 = STATEV(183)
		  sigcr23_h = PROPS(23)*PROPS(20)*sigcr23*epsf23/(sigcr23  
     1               + PROPS(23)*PROPS(20)*(epsf23-eps023))   	 
		  eps023_h = sigcr23_h/(PROPS(23)*PROPS(20))
		  STATEV(199) = 30.0d0
		  !write(*,*) 'Inner 3 fork', sigcr23_h, eps023_h, STRAIN(6)
		  IF (STRAIN(6).ge.eps023_h) THEN             
			  D23=(sigcr23_h/E230)*(1.0D0/(epsf23-eps023_h))*
     1                                       ((epsf23/STRAIN(6))-1)   
		      STATEV(174) = D23
			  STATEV(199) = 40.0d0
			  IF (D23.le.D23LIM) THEN
			      D23 = D23LIM
				  STATEV(174) = D23
			  END IF
		  ELSE IF (STRAIN(6).lt.eps023_h) THEN
		      D23 = PROPS(23)
		      STATEV(174) = D23
			  STATEV(199) = 50.0d0
	      END IF
             
			 IF (STATEV(174).lt.D23LIM) THEN
               STATEV(174) = D23LIM
               STATEV(106) = 2.0D0
               STATEV(120) = 0.0D0   ! Delete element
			   STATEV(199) = 60.0d0
             ELSE IF (STATEV(174).gt.1.0D0) THEN
               STATEV(174) = 1.0D0
			   STATEV(199) = 70.0d0
             END IF
       ELSEIF ((STATEV(106).eq.1.0D0).and.(del_EPS_23.lt.0.0d0)) THEN 
              D23 = STATEV(174)
	          STATEV(199) = 80.0d0
       END IF
       D23 = STATEV(174)
      
	  END IF !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	  







	  
	  
!C-------------
      ! TO DO: PUT D3LIM, same as D2LIM for 22
      ! For 3-direction
      IF (STATEV(111).eq.0.0D0) THEN
          D33 = 1.0D0
      ELSEIF ((STATEV(111).eq.1.0D0).and.(STRAIN(3).gt.STATEV(113))) 
     1 THEN
         E330 = STATEV(118)
         eps03 = STATEV(123)
         epsf3 = STATEV(127)
         D33=(sigcr33/E330)*(1.0D0/(epsf3-eps03))*((epsf3/STRAIN(3))-1)
      ELSEIF ((STATEV(111).eq.1.0D0).and.(STRAIN(3).le.STATEV(113))) 
     1 THEN
          D33 = STATEV(112)
      ELSEIF (D33.le.0.01D0) THEN
          D33 = 0.01D0
          STATEV(111) = 2.0D0
      END IF










!C-----------------------------------------------------------
!C UPDATE STATE VARIABLES
!C-----------------------------------------------------------
      
      !For 1-direction ***********************
      IF ((STATEV(102).lt.D1LIM).and.(STATEV(101).eq.1.0D0)) THEN
            STATEV(102) = D1LIM
      END IF	  
          D11 = STATEV(102)
          STATEV(103) = STRAIN(1)
          STATEV(104) = STRESS(1)
          IF (PROPS(11).gt.1.5d0) THEN
             STATEV(102) = D1LIM
             D11 = STATEV(102)
          END IF
	  
      !For 2-direction ***********************	  
      IF ((STATEV(107).lt.D2LIM).and.(STATEV(106).eq.1.0D0)) THEN
	  STATEV(107) = D2LIM
	  END IF	  
          D22 = STATEV(107)
          STATEV(108) = STRAIN(2)
          STATEV(109) = STRESS(2)
	        IF (PROPS(10).gt.1.5d0) THEN
	           STATEV(107) = D2LIM
	           D22 = STATEV(107)
	        END IF
  
      !@ For 12-plane ***********************	  
      IF ((STATEV(169).lt.D12LIM).and.(STATEV(106).eq.1.0D0)) THEN
	  STATEV(169) = D12LIM
	  END IF	  
          D12 = STATEV(169)
          STATEV(170) = STRAIN(4)
          STATEV(171) = STRESS(4)
	        IF (PROPS(10).gt.1.5d0) THEN
	           STATEV(169) = D12LIM
	           D12 = STATEV(169)
	        END IF  

      !@ For 23-plane ***********************	  
      IF ((STATEV(174).lt.D23LIM).and.(STATEV(106).eq.1.0D0)) THEN
	  STATEV(174) = D23LIM
	  END IF	  
          D23 = STATEV(174)
          STATEV(175) = STRAIN(6)
          STATEV(176) = STRESS(6)
	        IF (PROPS(10).gt.1.5d0) THEN
	           STATEV(174) = D23LIM
	           D23 = STATEV(174)
	        END IF 
  
	  
	  !For 3-direction	----- RJD: To be updated as needed
	  IF ((STATEV(112).lt.D3LIM).and.(STATEV(111).eq.1.0D0)) THEN
	  STATEV(112) = D3LIM
	  END IF	  
      D33 = STATEV(112)
      STATEV(113) = STRAIN(3)
      STATEV(114) = STRESS(3)

	  
      
	  
	  
	  
	  
	  
	  
	  !ASSEMBLE THE STIFFNESS MATRIX	  
      SS = 0.0D0
      SS(1,1) = 1.0D0/(D11*E11c)
      SS(2,2) = 1.0D0/(D22*E22c)
      SS(3,3) = 1.0D0/(D33*E33c)
      SS(4,4) = 1.0D0/(D12*G12c)
      SS(5,5) = 1.0D0/(G13c)
      SS(6,6) = 1.0D0/(D23*G23c)
      SS(2,1) = -v12c/(E11c)
      SS(1,2) = SS(2,1)
      SS(3,1) = -v13c/(E11c)
      SS(1,3) = SS(3,1)
      SS(3,2) = -v23c/E22c
      SS(2,3) = SS(3,2)



      STATEV(193) = D11*E11c      !RJD-CHECK-FOR-AXIAL DAMAGE FOR PARIS LAW
      STATEV(194) = D22*E22c      !RJD-CHECK-FOR-TRANSVERSE DAMAGE FOR PARIS LAW
	  STATEV(191) = D12*G12c
	  STATEV(192) = D23*G23c
      
      ! Invert the compliance matrix
      CALL t_matrixInverse(SS, CC, 6, 6)

!C-----------------------------------------------------------
!C STRESS UPDATE
!C-----------------------------------------------------------
      !!!!Stress Update
      STRESS = MATMUL(CC, STRAIN)
	  
              
      !Store stress at the end of the increment	   
      STATEV(161) = STRESS(1)
      STATEV(162) = STRESS(2)
      STATEV(163) = STRESS(3)
      STATEV(164) = STRESS(4)
      STATEV(165) = STRESS(5)
      STATEV(166) = STRESS(6)
	  
      !!!! Pass stress to SIG variable
      SIG = STRESS

	  
      !IF (STATEV(100).lt.0.5d0) THEN
      !   write(*,*) 'Last', STATEV(7), STATEV(9)
      !ENDIF
      
	  STATEV(100) = 1.0d0
      RETURN
      END







      
      
      














 
      











!C---------------------------------------------------------------------
!C           SUBROUTINES START HERE
!C---------------------------------------------------------------------

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE calc_stress(Q_IN, EPS_IN, SIG_OUT)
      REAL*8, INTENT(IN) :: Q_IN(6,6)
      REAL*8, INTENT(IN) :: EPS_IN(6)
      REAL*8, INTENT(OUT) :: SIG_OUT(6)

      SIG_OUT = MATMUL(Q_IN, EPS_IN)

      RETURN
      END





!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !Calculates isotropic stiffness matrix
      SUBROUTINE calc_stiffness_iso(Em, vm, CC)
      REAL*8, INTENT(OUT):: CC(6,6)
      REAL*8:: Constant
      REAL*8, INTENT(IN)::Em, vm


        Constant = Em/((1.0d0 + vm)*(1.0d0 - 2.0d0*vm))
        CC = 0.0d0
        CC(1,1) = Constant*(1.0d0 - vm)
        CC(1,2) = Constant*vm
        CC(2,1) = CC(1,2)
        CC(1,3) = Constant*vm
        CC(3,1) = CC(1,3)
        CC(2,3) = Constant*vm
        CC(3,2) = CC(2,3)
        CC(2,2) = Constant*(1.0d0 - vm)
        CC(3,3) = Constant*(1.0d0 - vm)
        CC(4,4) = Constant*(1.0d0 - 2.0d0*vm)/2.0d0
        CC(5,5) = Constant*(1.0d0 - 2.0d0*vm)/2.0d0
        CC(6,6) = Constant*(1.0d0 - 2.0d0*vm)/2.0d0   

      
      RETURN
      END


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      !Calculates transversely isotropic stiffness matrix
      SUBROUTINE calc_stiffness_transiso(E11c, E22c, v12c, 
     1                                 v23c, G12c, G23c, Q)
      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE
      REAL*8, INTENT(IN)::E11c, E22c, v12c, v23c, G12c, G23c
      REAL*8::Delta
      REAL*8, INTENT(OUT)::Q(6,6)


      Delta = 2.0d0*E22c*v12c**2 + E11c*(v23c - 1.0d0)
      Q = 0.0
      Q(1,1) = E11c**2*(v23c - 1.0d0)/Delta
      Q(1,2) = -E11c*E22c*v12c/Delta
      Q(1,3) = Q(1,2)
      Q(2,1) = Q(1,2)
      Q(2,2) = E22c*(E22c*v12c**2-E11c)/((1.0d0 + v23c)*Delta)
      Q(2,3) = -(E22c*(E22c*v12c**2 + E11c*v23c))/((1.0+v23c)*Delta)
      Q(3,1) = Q(1,3)
      Q(3,2) = Q(2,3)
      Q(3,3) = Q(2,2)
      Q(4,4) = G12c
      Q(5,5) = G12c
      Q(6,6) = G23c
      
      RETURN
      END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE calc_eqstrain(EP, j2eps)
      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE
      REAL*8, INTENT(IN):: EP(6)
      REAL*8, INTENT(OUT):: j2eps
      REAL*8:: exx, eyy, ezz, gxy, gxz, gyz
  
      exx = (2.0d0/3.0d0)*EP(1) - (1.0d0/3.0d0)*EP(2) 
     1                 - (1.0d0/3.0d0)*EP(3)
      eyy = (2.0d0/3.0d0)*EP(2) - (1.0d0/3.0d0)*EP(1) 
     1                 - (1.0d0/3.0d0)*EP(3)     
      ezz = (2.0d0/3.0d0)*EP(3) - (1.0d0/3.0d0)*EP(1) 
     1                 - (1.0d0/3.0d0)*EP(2)
      gxy = EP(4)
      gxz = EP(5)
      gyz = EP(6)


      j2eps = (2.0d0/3.0d0)*DSQRT( 1.5d0*(exx**2 + eyy**2 + ezz**2) 
     1            + 0.75d0*(gxy**2 + gxz**2 + gyz**2))
           
      RETURN
      END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE calc_ccm(E11f, E22f, v12f, v23f, G12f, G23f, 
     1    Em, vm, Vf, E11c, E22c, v12c, v23c, G12c, G23c)
      
      !IMPLICIT NONE
      !INCLUDE 'ABA_PARAM.INC'
      REAL*8, INTENT(IN) :: E11f, E22f, G12f, G23f
      REAL*8, INTENT(IN) :: v12f, v23f, Vf, Em, vm
      REAL*8, INTENT(OUT) :: E11c, E22c, G12c, G23c, v12c, v23c
      REAL*8 :: gm_s, del_s, etaf_s, etam_s, eta4_s, Gm
      REAL*8 :: v21f, gm_num, gm_den, del_num, del_den
      REAL*8 :: etaf, etaf_num, etaf_den, etam_num, etam_den
      REAL*8 :: v12c_num1, v12c_num2, v12c_den

      
      !Assigning additional required variables
      Gm = 0.5d0*Em/(1.0d0 + vm)
      v21f = v12f*E22f/E11f
      
      !Calculating partial portions of the expressions
      gm_num = 2.0d0*Em*v21f*(1.0d0-Vf)*(v12f - vm)
      gm_den = E22f*(1.0d0+vm)*(1.0d0 + Vf*(1.0d0 - 2.0d0*vm)) +
     1         Em*(1.0d0 - v23f - 2.0d0*v12f*v21f)*(1 - Vf)
      del_num = 2.0d0*E22f*vm*Vf*(vm - v12f)
      del_den = E22f*(1.0d0+vm)*(1.0d0 + Vf*(1.0d0 - 2.0d0*vm)) +
     1         Em*(1.0d0 - v23f - 2.0d0*v12f*v21f)*(1.0d0 - Vf)

      etaf_num = E11f*Vf + ((1.0d0 - v12f*v21f)*Em + 
     1           vm*v21f*E11f)*(1.0d0 - Vf)
      etaf_den = E11f*Vf + Em*(1.0d0 - Vf)

      etam_num = ((1.0d0 - vm*vm)*E11f + vm*v12f*Em)*Vf + Em*(1.0d0- Vf)
      etam_den = E11f*Vf + Em*(1.0d0 - Vf)
      gm_s = gm_num/gm_den
      del_s = del_num/del_den
      etaf_s = etaf_num/etaf_den
      etam_s = etam_num/etam_den
      eta4_s = (3.0d0 - 4.0d0*vm + (Gm/G23f))/(4.0d0*(1.0d0 - vm))
      v12c_num1 = (1.0d0 - Vf)*(1.0d0 - v23f - 2.0d0*v12f*v21f)*vm*Em
      v12c_num2 = (vm + Vf*(2.0d0*v12f - vm) + 
     1                   (vm*vm*(1.0d0 - 2.0d0*Vf*v12f - Vf)))*E22f
      v12c_den = (1.0d0 - Vf)*(1.0d0 - v23f - 2.0d0*v12f*v21f)*Em +
     1           (1.0d0 + Vf + (1.0d0 - Vf)*vm - 2.0d0*Vf*vm*vm)*E22f

      !Passing back variables to the main program
      E11c = E11f*(1.0d0 + gm_s)*Vf + Em*(1.0d0 + del_s)*(1.0d0 - Vf)
      E22c = 1.0d0/((etaf_s*Vf/E22f) + (etam_s*(1.0d0 - Vf)/Em))
      G12c = Gm*(((Gm+G12f) - Vf*(Gm-G12f))/((Gm+G12f) + Vf*(Gm-G12f)))
      G23c = (Vf + eta4_s*(1.0d0 - Vf))/((Vf/G23f) + 
     1       (eta4_s*(1.0d0-Vf)/Gm))
      v23c = (0.5d0*E22c/G23c) - 1.0d0
      v12c = (v12c_num1 + v12c_num2)/v12c_den


      RETURN
      END


!!--------------------------------------------------------------------
!!!    SUBROUTINE FOR CCM

      SUBROUTINE CCM_CALC_tow(CCMinput,CCMoutput,ntens,nstatv,statev)
        
        IMPLICIT NONE      
        INTEGER*4::NTENS,NSTATV,KINC
        REAL*8,DIMENSION(NSTATV):: statev
        REAL*8,DIMENSION(8):: CCMinput
        REAL*8,DIMENSION(8):: CCMoutput
        REAL*8:: E11c, E22c, E33c, G12c, G13c, G23c
        REAL*8:: Ef1, Ef2, Gf23, Gf12, vf12, Volf, Volm, Em, vm
        REAL*8:: vf23, vf21, Gm, gamma, delta, etaf, etam, eta4
        REAL*8:: v12c, v23c
      
      Ef1 = CCMinput(1)
      Ef2 = CCMinput(2)
      Gf12 = CCMinput(3)
      Gf23 = CCMinput(4)
      vf12 = CCMinput(5)
      Volf = CCMinput(6)
      Em = CCMinput(7)
      vm = CCMinput(8)
	  
      Volm = 1.0D0 - Volf
      vf23 = ((Ef2/(2.0D0*Gf23)) - 1.0D0)
      vf21 = Ef2*vf12/Ef1
      Gm = Em/(2.0D0*(1.0D0 + vm))
      gamma=2*vf21*Em*(1-Volf)*(vf12-vm)/(Ef2*(1+vm)*(1+Volf*(1-2*vm))
     1        + Em*(1 - vf23 -2.0D0*vf12*vf21)*(1-Volf))
      delta = 2*Ef2*vm*Volf*(vm-vf12)/(Ef2*(1+vm)*(1+Volf*(1-2*vm)) 
     2        + Em*(1-vf23-2*vf12*vf21)*(1-Volf))
      E11c = Ef1*(1+gamma)*Volf + Em*(1 + delta)*(1 - Volf)	 
      etaf = (Ef1*Volf+((1-vf12*vf21)*Em + vm*vf21*Ef1)*(1-Volf))/
     1       (Ef1*Volf + Em*(1-Volf))	 
      etam = (((1- (vm**2))*Ef1+ vm*vf12*Em)*Volf + Em*Volm)/
     2       (Ef1*Volf+Em*(1-Volf))	 
      E22c = 1/((etaf*Volf/Ef2)+(etam*(1.0D0 - Volf)/Em))
     
	 
      G12c =Gm*(((Gm+Gf12)-Volf*(Gm-Gf12))/((Gm+Gf12) + Volf*(Gm-Gf12)))
      eta4 = (3.0D0 - 4.0D0*vm + (Gm/Gf23))/(4.0D0*(1.0D0 - vm))
      G23c = (Volf+eta4*(1-Volf))/((Volf/Gf23)+(eta4/Gm)*(1-Volf))
	  
      v12c = (((1-Volf)*(1 - vf23 - 2*vf12*vf21))*vm*Em 
     1 +(vm+Volf*(2*vf12-vm)+((vm**2)*(1-2*Volf*vf12-Volf)))*Ef2)
     2 /(((1-Volf)*(1-vf23-2*vf12*vf21))*Em
     3 +(1+Volf+(1-Volf)*vm - 2*Volf*(vm**2))*Ef2)

      v23c = (E22c/(2.0D0*G23c)) - 1.0D0
	 
      !Output of the CCM calculation
      CCMoutput(1) = E11c
      CCMoutput(2) = E22c
      CCMoutput(3) = v12c
      CCMoutput(4) = v23c
      CCMoutput(5) = G12c
      CCMoutput(6) = G23c
      CCMoutput(7) = Gm
      CCMoutput(8) = eta4

      RETURN
      END

	  
	  
	  
	  
	  
!C----------------------------------------------------------------------	  
!C            INTERPOLATION FUNCTION
!C----------------------------------------------------------------------
C!      FUNCTION interpolt(N_INTERVALS,t_cure,TIME,INP_TIME,
C!     1     INP_FUNC,KINC)
C!
C!      IMPLICIT NONE
C!      INTEGER*4:: J, N_INTERVALS, KINC
C!      REAL*8,DIMENSION(N_INTERVALS):: INP_TIME, INP_FUNC
C!      REAL*8:: t_cure, TIME
C!      REAL*8:: t_1(N_INTERVALS), t_2(N_INTERVALS)
C!      REAL*8:: phi_1(N_INTERVALS), phi_2(N_INTERVALS), toll
C!      REAL*8:: t_1_in, t_2_in, phi_1_in, phi_2_in
C!      REAL*8:: interpolt
C!      
C!      
C!       toll = 1e-05			
C!            interpolt = 0.0D0
C!              
C!              ! STORE THE VALUES FROM THE FILE              
C!              do J = 1, (N_INTERVALS-1)
C!              
C!              t_1(J) = INP_TIME(J)
C!              t_2(J) = INP_TIME(J+1)
C!                
C!              phi_1(J) = INP_FUNC(J)
C!              phi_2(J) = INP_FUNC(J+1)
C!                
C!              !write(*,*) t_1(J), t_2(J)
C!              !write(*,*) phi_1(J), phi_2(J)
C!              ! possible errors in the input files:       
C!900   FORMAT (A40)              
C!              if (t_1(J) == t_2(J)) then
C!                write(*,900) ' * ----   ERROR IN <INP_TIME.txt>   ---*'
C!                write(*,900) ' * ----	two data point are equal  ----*'
C!		        write(*,'(A30,I5,A5)') '* ----	         J = ',J      ,'----*'
C!                  stop
C!              endif 
C!              
C!          
C!      IF ( ( TIME >=(t_1(J)-toll) ).AND.( TIME <= (t_2(J)+toll) ) )then
C!              
C!                t_1_in = t_1(J)
C!                t_2_in = t_2(J)
C!                
C!                phi_1_in = phi_1(J)
C!                phi_2_in = phi_2(J)
C!
C!       endif
C!              enddo
C! 
C!             !write(*,*) 't_1_in ',t_1_in,'t_2_in ',t_2_in 
C!             !write(*,*) 'phi_1_in',phi_1_in,'phi_2_in',phi_2_in
C!             IF  ( (TIME <= t_cure) .and. (TIME < t_2_in ) ) then
C!                  
C!             interpolt=(phi_2_in-phi_1_in)/(t_2_in-t_1_in)*(TIME-t_1_in)
C!     &				 + phi_1_in
C!				
C!                  !write(*,*) TIME, interpol    
C!                  ! pause                             
C!		    ENDIF
C!              
C!			
C!      
C!      return
C!      end	  
	  
	  
	  


!C-------------------------------------------------------------------
!C-------------------------------------------------------------------
!C SUBROUTINES TO INVERT THE MATRIX
!C-------------------------------------------------------------------
!C-------------------------------------------------------------------
      subroutine t_matrixInverse(a,y,n,np)

      !taken from http://www.mpi-hd.mpg.de/astrophysik/HEA/internal/Numerical_Recipes/f2-3.pdf
      !a goes in, y=inverse(a) goes out

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      INTEGER np,indx(np),n,flag
      REAL*8 a(np,np),y(np,np),d
      integer i,j

      do i=1,n !Set up identity matrix.
      do j=1,n
      y(i,j)=0.
      enddo
      y(i,i)=1.
      enddo
      call t_ludcmp(a,n,np,indx,d,flag) !Decompose the matrix just once.
      do  j=1,n                !Find inverse by columns.
      call t_lubksb(a,n,np,indx,y(1,j))
      !Note that FORTRAN stores two-dimensional matrices by column, so y(1,j) is the
      !address of the jth column of y.
      enddo

      return
      end

!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
! same as above, quad precision

      subroutine t_qmatrixInverse(a,y,n,np)

      !taken from http://www.mpi-hd.mpg.de/astrophysik/HEA/internal/Numerical_Recipes/f2-3.pdf
      !a goes in, y=inverse(a) goes out

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      INTEGER np,indx(np),n,flag
      REAL*16 a(np,np),y(np,np),d
      integer i,j

      do i=1,n !Set up identity matrix.
      do j=1,n
      y(i,j)=0.
      enddo
      y(i,i)=1.
      enddo
      call t_ludcmp(a,n,np,indx,d,flag) !Decompose the matrix just once.
      do  j=1,n                !Find inverse by columns.
      call t_lubksb(a,n,np,indx,y(1,j))
      !Note that FORTRAN stores two-dimensional matrices by column, so y(1,j) is the
      !address of the jth column of y.
      enddo

      return
      end

!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------

      subroutine t_solveLinSysLU(A,b,n,flag)
      !solves Ax=b via LU-decomposition
      !double precision
      !input: A,b
      !output: x in place of b
      !        successfull: flag=1, not s.: flag=0

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      integer::n
      integer::indx(n)
      integer::flag
      real*8::A(n,n),b(n),d

!      write(*,*) 'A='
!      write(*,*) A
!      write(*,*) 'b='
!      write(*,*) b
!      write(*,*) '---'


      call t_ludcmp(A,n,n,indx,d,flag)
      if (flag.eq.1) then
      call t_lubksb(A,n,n,indx,b)
      endif

      return
      end


!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------

      subroutine t_qSolveLinSysLU(A,b,n,flag)
      !solves Ax=b via LU-decomposition
      !quad precision
      !input: A,b
      !output: x in place of b
      !        successfull: flag=1, not s.: flag=0

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      integer::n
      integer::indx(n)
      integer::flag
      real*16::A(n,n),b(n),d

      call t_qludcmp(A,n,n,indx,d,flag)
      call t_qlubksb(A,n,n,indx,b)

      return
      end

!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------

      SUBROUTINE t_ludcmp(a,n,np,indx,d,flag)

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      INTEGER n,np,indx(n),NMAX,flag
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=10,TINY=1.0e-20) !Largest expected n, and a small number.
      !Given a matrix a(1:n,1:n), with physical dimension np by np, this routine replaces it by
      !the LU decomposition of a rowwise permutation of itself. a and n are input. a is output,
      !arranged as in equation (2.3.14) above; indx(1:n) is an output vector that records the
      !row permutation effected by the partial pivoting; d is output as �1 depending on whether
      !the number of row interchanges was even or odd, respectively. This routine is used in
      !combination with lubksb to solve linear equations or invert a matrix.
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX) !vv stores the implicit scaling of each row.
      flag=1 !So far successfull ...
      d=1. !No row interchanges yet.
      do i=1,n !Loop over rows to get the implicit scaling information
      aamax=0.
      do j=1,n
      if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
      enddo
      if (aamax.eq.0.) then
        !write(*,*) 'singular matrix in ludcmp' !No nonzero largest element.
        !write(*,*)
        !write(*,*) a
        !CALL myExit()
        flag=0
        exit
      endif
      vv(i)=1./aamax !Save the scaling.
      enddo
      if (flag.eq.1) then
      do j=1,n !This is the loop over columns of Crout�s method.
      do i=1,j-1 !This is equation (2.3.12) except for i = j.
      sum=a(i,j)
      do k=1,i-1
      sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
      enddo
      aamax=0. !Initialize for the search for largest pivot element.
      do i=j,n !This is i = j of equation (2.3.12) and i = j+1. . .N
      sum=a(i,j) !of equation (2.3.13).
      do k=1,j-1
      sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
      dum=vv(i)*abs(sum) !Figure of merit for the pivot.
      if (dum.ge.aamax) then !Is it better than the best so far?
      imax=i
      aamax=dum
      endif
      enddo
      if (j.ne.imax)then !Do we need to interchange rows?
      do k=1,n !Yes, do so...
      dum=a(imax,k)
      a(imax,k)=a(j,k)
      a(j,k)=dum
      enddo
      d=-d !...and change the parity of d.
      vv(imax)=vv(j) !Also interchange the scale factor.
      endif
      indx(j)=imax
      if(a(j,j).eq.0.)a(j,j)=TINY
      !If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
      !For some applications on singular matrices, it is desirable to substitute TINY
      !for zero.
      if(j.ne.n)then !Now, finally, divide by the pivot element.
      dum=1./a(j,j)
      do i=j+1,n
      a(i,j)=a(i,j)*dum
      enddo
      endif
      enddo !Go back for the next column in the reduction.
      endif
      return
      END

!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
!same as above with quad precision

      SUBROUTINE t_qludcmp(a,n,np,indx,d,flag)

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      INTEGER n,np,indx(n),NMAX,flag
      REAL*16 d,a(np,np),TINY
      PARAMETER (NMAX=10,TINY=1.0e-20) !Largest expected n, and a small number.
      !Given a matrix a(1:n,1:n), with physical dimension np by np, this routine replaces it by
      !the LU decomposition of a rowwise permutation of itself. a and n are input. a is output,
      !arranged as in equation (2.3.14) above; indx(1:n) is an output vector that records the
      !row permutation effected by the partial pivoting; d is output as �1 depending on whether
      !the number of row interchanges was even or odd, respectively. This routine is used in
      !combination with lubksb to solve linear equations or invert a matrix.
      INTEGER i,imax,j,k
      REAL*16 aamax,dum,sum,vv(NMAX) !vv stores the implicit scaling of each row.
      flag=1 !So far successfull ...
      d=1. !No row interchanges yet.
      do i=1,n !Loop over rows to get the implicit scaling information
      aamax=0.
      do j=1,n
      if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
      enddo
      if (aamax.eq.0.) then
        !write(*,*) 'singular matrix in ludcmp' !No nonzero largest element.
        !write(*,*)
        !write(*,*) a
        !CALL myExit()
        flag=0
      endif
      vv(i)=1./aamax !Save the scaling.
      enddo
      if (flag.eq.1) then
      do j=1,n !This is the loop over columns of Crout�s method.
      do i=1,j-1 !This is equation (2.3.12) except for i = j.
      sum=a(i,j)
      do k=1,i-1
      sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
      enddo
      aamax=0. !Initialize for the search for largest pivot element.
      do i=j,n !This is i = j of equation (2.3.12) and i = j+1. . .N
      sum=a(i,j) !of equation (2.3.13).
      do k=1,j-1
      sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
      dum=vv(i)*abs(sum) !Figure of merit for the pivot.
      if (dum.ge.aamax) then !Is it better than the best so far?
      imax=i
      aamax=dum
      endif
      enddo
      if (j.ne.imax)then !Do we need to interchange rows?
      do k=1,n !Yes, do so...
      dum=a(imax,k)
      a(imax,k)=a(j,k)
      a(j,k)=dum
      enddo
      d=-d !...and change the parity of d.
      vv(imax)=vv(j) !Also interchange the scale factor.
      endif
      indx(j)=imax
      if(a(j,j).eq.0.)a(j,j)=TINY
      !If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
      !For some applications on singular matrices, it is desirable to substitute TINY
      !for zero.
      if(j.ne.n)then !Now, finally, divide by the pivot element.
      dum=1./a(j,j)
      do i=j+1,n
      a(i,j)=a(i,j)*dum
      enddo
      endif
      enddo !Go back for the next column in the reduction.
      endif
      return
      END

!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------

      SUBROUTINE t_lubksb(a,n,np,indx,b)

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      !Solves the set of n linear equations A � X = B. Here a is input, not as the matrix A but
      !rather as its LU decomposition, determined by the routine ludcmp. indx is input as the
      !permutation vector returned by ludcmp. b(1:n) is input as the right-hand side vector B,
      !and returns with the solution vector X. a, n, np, and indx are not modified by this routine
      !and can be left in place for successive calls with different right-hand sides b. This routine
      !takes into account the possibility that b will begin with many zero elements, so it is efficient
      !for use in matrix inversion.
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0 !When ii is set to a positive value, it will become the index
      !of the first nonvanishing element of b. We now do
      !the forward substitution, equation (2.3.6). The only new
      !wrinkle is to unscramble the permutation as we go.
      do i=1,n
      ll=indx(i)
      sum=b(ll)
      b(ll)=b(i)
      if (ii.ne.0)then
      do j=ii,i-1
      sum=sum-a(i,j)*b(j)
      enddo
      else if (sum.ne.0.) then
      ii=i !A nonzero element was encountered, so from now on we will
      endif !have to do the sums in the loop above.
      b(i)=sum
      enddo
      do i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
      sum=b(i)
      do j=i+1,n
      sum=sum-a(i,j)*b(j)
      enddo
      b(i)=sum/a(i,i) !Store a component of the solution vector X.
      enddo
      return !All done!
      END


!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
!same as above with quad precision

      SUBROUTINE t_qlubksb(a,n,np,indx,b)

      !INCLUDE 'ABA_PARAM.INC'
      IMPLICIT NONE

      INTEGER n,np,indx(n)
      REAL*16 a(np,np),b(n)
      !Solves the set of n linear equations A � X = B. Here a is input, not as the matrix A but
      !rather as its LU decomposition, determined by the routine ludcmp. indx is input as the
      !permutation vector returned by ludcmp. b(1:n) is input as the right-hand side vector B,
      !and returns with the solution vector X. a, n, np, and indx are not modified by this routine
      !and can be left in place for successive calls with different right-hand sides b. This routine
      !takes into account the possibility that b will begin with many zero elements, so it is efficient
      !for use in matrix inversion.
      INTEGER i,ii,j,ll
      REAL*16 sum
      ii=0 !When ii is set to a positive value, it will become the index
      !of the first nonvanishing element of b. We now do
      !the forward substitution, equation (2.3.6). The only new
      !wrinkle is to unscramble the permutation as we go.
      do i=1,n
      ll=indx(i)
      sum=b(ll)
      b(ll)=b(i)
      if (ii.ne.0)then
      do j=ii,i-1
      sum=sum-a(i,j)*b(j)
      enddo
      else if (sum.ne.0.) then
      ii=i !A nonzero element was encountered, so from now on we will
      endif !have to do the sums in the loop above.
      b(i)=sum
      enddo
      do i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
      sum=b(i)
      do j=i+1,n
      sum=sum-a(i,j)*b(j)
      enddo
      b(i)=sum/a(i,i) !Store a component of the solution vector X.
      enddo
      return !All done!
      END
