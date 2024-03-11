!C!For: Tow model - nonlinear behavior with crackband (Hashin initiation for 22, 12, 23)
!C!C Author: Royan Dmello, Aerospace Engineering, University of Michigan
!C!C Updated: 06/20/2022
      SUBROUTINE VUMAT_TOW(STATEV,DDSDDE,
     2 CMNAME,
     3 NDI, NSHR,NTENS,NSTATV,PROPS,NPROPS,
     4 CELENT, SIG, TIMEVUMAT)
     
      IMPLICIT NONE

      CHARACTER*80 CMNAME
      INTEGER:: NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NTENS,NSTATV,NPROPS
      INTEGER:: NDI,NSHR
      REAL*8:: STRESS(6), DSTRAN(6), STRAN(6), DDSDDT(6), DRPLDE(6)
      REAL*8:: DDSDDE(6,6)
      REAL*8:: STATEV(NSTATV), PROPS(NPROPS)
      REAL*8:: TIME(2)
      REAL*8:: DTIME,TEMP,DTEMP,CELENT,PNEWDT,DRPLDT
      REAL*8:: SSE, SPD, SCD, RPL
      REAL*8:: PREDEF(1)
      REAL*8:: DPRED(1)
      REAL*8:: COORDS(3,3), DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)
      
      
      
      REAL*8:: STRAIN(NTENS), OLDSTRESS(NTENS)
      REAL*8:: SIG(NTENS), STRAN_DEFG(NTENS), DSTRAN_DEFG(NTENS)
      REAL*8:: SS(6,6), CC(6,6)
      REAL*8:: Constant, E, nu
      REAL*8:: E11, E22, E33, G12, G13, G23
      REAL*8:: E110, E220, E330
      REAL*8:: v12, v13, v23
      REAL*8:: GIc, sigcr11, BZ11
      REAL*8:: D11, eps01, h1, epsf1
      REAL*8:: GIc2, sigcr22, BZ22
      REAL*8:: D22, eps02, h2, epsf2
      REAL*8:: D33, eps03, h3, epsf3
      REAL*8:: CCMinput(8)
      REAL*8:: CCMoutput(8)
      REAL*8:: Em, vm
      REAL*8:: E11c, E22c, E33c, G12c, G13c, G23c
      REAL*8:: v12c, v13c, v23c, Volf, Volf_mean
      REAL*8:: Ef1, Ef2, vf12, vf23, Gf12, Gf23
      REAL*8:: D1LIM, D2LIM, D3LIM, D12LIM, D23LIM
      REAL*8:: Cm(6,6), Cf(6,6)	  
      REAL*8:: eps_m(6), eps_m11(6), eps_m22(6), eps_m33(6)
      REAL*8:: eps_m12(6), eps_m13(6), eps_m23(6)
      REAL*8:: Em0, vm0, SIG_Y, k1, k2, EPS_Y, Dm, Vf
      REAL*8:: Em_s, vm_s, EPS_J2, SIG_J2, Em_s_prev, vm_s_prev
      INTEGER:: I
      REAL*8:: GIc1
      REAL*8:: D12, D23, GIc12, GIc23, sigcr12, sigcr23
      REAL*8:: HASHIN2
      REAL*8:: eps012, eps023, epsf12, epsf23, BZ12, BZ23
      REAL*8:: E120, E230	  
      real*8:: TIMEVUMAT

       


!C------------------------------------------
!C INPUTS AND STATE VARIABLES
!C------------------------------------------



      !write(*,*) 'in Tow'
      D1LIM = 0.02d0
      D2LIM = 0.02d0
      D3LIM = 0.02d0
      D12LIM = 0.02d0
      D23LIM = 0.02d0

!C    !Fiber elastic properties - hardcoded here
      Ef1 = 265000.0d0        
      Ef2 = 20100.0d0
      vf12 = 0.2d0
      vf23 = 0.25d0
      Gf12 = 32200.0d0
      Gf23 = 8040.0d0	  

      
	  	  
      Em0 = PROPS(1)	  
      vm0 = PROPS(2)  
      SIG_Y = PROPS(3)   !15.66d0 matrix yield stress   
      k1 = PROPS(4)      ! k1 = 3308.0d0 hardening parameter for matrix
      k2 = PROPS(5)      ! k2 = 56.92d0 hardening parameter for matrix
      Volf = PROPS(6)	  
      sigcr11 = PROPS(7)
      GIc1 = PROPS(8)
      sigcr22 = PROPS(9)
      GIc2 = PROPS(10)
      sigcr12 = PROPS(11)
      GIc12 = PROPS(12)
      sigcr23 = PROPS(13)	  
      GIc23 = PROPS(14)    

      EPS_Y = SIG_Y/Em0


      STATEV(139) = PROPS(9)
      STATEV(140) = PROPS(11)
      STATEV(141) = PROPS(13)

     	  	  
      !Initialize state variables
	  IF (STATEV(100).lt.0.5d0) THEN
           STATEV(9) = Em0            !Em: Young's modulus
           STATEV(10) = vm0           !vm: Poisson's ratio
           STATEV(33) = NOEL
           STATEV(32) = Volf
           STATEV(101) = 0.0d0
		   STATEV(106) = 0.0d0
		   STATEV(161) = 0.0d0
		   STATEV(162) = 0.0d0
		   STATEV(163) = 0.0d0
		   STATEV(164) = 0.0d0
		   STATEV(165) = 0.0d0
		   STATEV(166) = 0.0d0
	  END IF
        
      STATEV(137) = Volf
      Vf = STATEV(137)
     
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
	
      

      !---------------------------------------------------
      !Multiscaling done if not yet in crackband
      IF ((STATEV(101).lt.0.5d0).and.(STATEV(106).lt.0.5d0)) THEN
	  !write(*,*) 't2'
      
	      Em_s_prev = STATEV(9)  
	      vm_s_prev = STATEV(10)
          STATEV(7) = Em_s_prev
	      STATEV(8) = vm_s_prev
      


        !Calculate stiffness terms for matrix and fiber
        CALL calc_stiffness_iso(Em_s_prev, vm_s_prev, Cm)
             STATEV(14) = Cm(1,1)  
             STATEV(15) = Cm(2,2)  
             STATEV(16) = Cm(3,3)  
             STATEV(17) = Cm(4,4)  
             STATEV(18) = Cm(5,5)  
             STATEV(19) = Cm(6,6)  
	    
        CALL calc_stiffness_transiso(Ef1, Ef2, vf12, 
     1                                 vf23, Gf12, Gf23, Cf)
             STATEV(20) = Cf(1,1)
             STATEV(21) = Cf(2,2)  
             STATEV(22) = Cf(3,3)  
             STATEV(23) = Cf(4,4) 
             STATEV(24) = Cf(5,5)  
             STATEV(25) = Cf(6,6)  
	 
	  
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
	   
      STATEV(26) = eps_m(1) 
      STATEV(27) = eps_m(2)  
      STATEV(28) = eps_m(3)  
      STATEV(29) = eps_m(4)  
      STATEV(30) = eps_m(5)  
      STATEV(31) = eps_m(6)  


      !Calculate equivalent matrix strain from matrix constituent strains
      CALL calc_eqstrain(eps_m, EPS_J2)     
      IF (EPS_J2.LE.EPS_Y) THEN
      SIG_J2 = Em0*EPS_J2      !Secant will be equal to initial modulus
      ELSEIF (EPS_J2.GT.EPS_Y) THEN
      STATEV(13) = 1.0d0
      SIG_J2 = SIG_Y - (k1/k2)*(DEXP(-k2*EPS_J2) - DEXP(-k2*SIG_Y/Em0))
      ENDIF

      
      !Compute secant modulus and secant Poisson's ratio
	  	  IF (dabs(EPS_J2).gt.0.0d0) THEN
		  Em_s = DABS(SIG_J2/EPS_J2)
		  ELSE
		  Em_s = STATEV(9)
		  END IF
	  

      vm_s = 0.5d0 + (Em_s/Em0)*(vm0 - 0.5d0)


      IF (STATEV(100).lt.0.5d0) THEN
          Em_s = Em0
          vm_s = vm0
      END IF

      STATEV(9) = Em_s
      STATEV(10) = vm_s
      STATEV(11) = EPS_J2
      STATEV(12) = SIG_J2

      !----------------------------------------------------------------
 
!-------------------------------------------------------
      !Compute effective composite stiffness
!-------------------------------------------------------
	  
      CCMinput(1) = Ef1
      CCMinput(2) = Ef2
      CCMinput(3) = Gf12
      CCMinput(4) = Gf23
      CCMinput(5) = vf12
      CCMinput(6) = Volf
      CCMinput(7) = STATEV(9)
      CCMinput(8) = STATEV(10)
      
      CALL CCM_CALC_tow(CCMinput,CCMoutput)
	  
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
	  END IF
	  

      E11c = STATEV(131)
      E22c = STATEV(132)
      v12c = STATEV(133)
      v23c = STATEV(134)
      G12c = STATEV(135)
      G23c = STATEV(136)
	  E33c = E22c
      v13c = v12c
      G13c = G12c





!C---------------------------------------------------------
!C  WHEN THE MATERIAL FIRST ENTERS CRACKBAND
!C---------------------------------------------------------
      !C     ! For 1-direction
      !When material first enters crackband

      IF ((STRESS(1).ge.sigcr11).and.(STATEV(101).lt.0.5D0)) THEN
      STATEV(101) = 1.0D0      
      STATEV(116) = STRESS(1)/STRAIN(1)
      
      !This will be used to track max strain yet in CB
      STATEV(103) = STRAIN(1) 
        
      eps01 = sigcr11/STATEV(116)
      STATEV(121) = eps01
 
         
      ! Crackband parameters stored once at crack initiation
      h1 = CELENT
      STATEV(124) = CELENT 
      epsf1 = (2.0D0*GIc1)/(h1*sigcr11)
      STATEV(125) = epsf1  
               
      ! Check Bazant limit
      BZ11 = 2.0D0*E11c*GIc1/(sigcr11**2)  
         IF (h1.ge.0.95*BZ11) THEN
           write(*,*) 'Char. element length is too large in 11!!'
           write(*,*) 'Reduce mesh size'
           write(*,*) 'h1', h1
           write(*,*) 'BZ11', BZ11
           write(*,*) 'GIc1', GIc1
           write(*,*) 'sigcr11', sigcr11
           write(*,*) 'E11c', E11c
           write(*,*) 'Em', Em			
           CALL XPLB_EXIT
         END IF
	  END IF


	  
	  
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	  
      ! EVALUATE HASHIN CONDITION FOR 22, 12 and 23
      IF (STATEV(106).lt.0.5d0) THEN	  
	   HASHIN2 = (STRESS(2)/sigcr22)**2 + (STRESS(4)/sigcr12)**2
     1	                              + (STRESS(6)/sigcr23)**2
		 
		 IF ((HASHIN2.ge.1.0d0).and.(STRESS(2).gt.0.0d0)) THEN
			STATEV(106) = 1.0d0        ! Turn on crackband flag for 22,12,23
            sigcr22 = STRESS(2)
            sigcr12 = dabs(STRESS(4))
            sigcr23 = dabs(STRESS(6))
            STATEV(184) = sigcr22      !Store initiation stress conditions for CB
            STATEV(185) = sigcr12      !                ''
            STATEV(186) = sigcr23      !                ''
            ! Store initial crackband stiffness
            
            !STATEV(117) = STATEV(132)
            STATEV(117) = DABS(STRESS(2)/STRAIN(2))
            STATEV(178) = DABS(STRESS(4)/STRAIN(4))
            STATEV(179) = DABS(STRESS(6)/STRAIN(6))
                 
			
			! Store crackband critical strainß
            eps02 = DABS(STATEV(184)/STATEV(117))
            eps012 = dabs(STATEV(185)/STATEV(178))  
            eps023 = dabs(STATEV(186)/STATEV(179))
            STATEV(122) = eps02
            STATEV(180) = eps012          
            STATEV(181) = eps023          
            
			
			! Crackband parameters stored once at crack initiation
            h2 = CELENT
            STATEV(124) = CELENT 
            epsf2 = (2.0D0*GIc2)/(h2*STATEV(184)) 
		    epsf12 = (2.0D0*GIc12)/(h2*STATEV(185))  
		    epsf23 = (2.0D0*GIc23)/(h2*STATEV(186))  	 	 
            STATEV(126) = epsf2 
            STATEV(182) = epsf12        
            STATEV(183) = epsf23


            !These SDV will be used to keep track to max strain yet in CB
            STATEV(108) = STRAIN(2)
            STATEV(170) = STRAIN(4)
            STATEV(175) = STRAIN(6)

            !Check Bazant limit
            BZ22 = 2.0D0*E22c*GIc2/(STATEV(184)**2)  
            BZ12 = 2.0D0*G12c*GIc12/(STATEV(185)**2)     
            BZ23 = 2.0D0*G23c*GIc23/(STATEV(186)**2)  	  	  
            IF ((h2.ge.(0.95*BZ22)).or.(h2.ge.(0.95*BZ12)).or.
     1                               (h2.ge.(0.95*BZ23))) THEN
                write(*,*) 'Char. element length is too large in 2!!'
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
         END IF
      END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      
      
      
!C-----------------------------------------------------------	  
!C-----------------------------------------------------------
!C STIFFNESS DEGRADATION FOR CRACKED/DAMAGED SOLID
!C-----------------------------------------------------------
!C-----------------------------------------------------------
   
!C	!For 11-direction ****************
        IF (STATEV(101).eq.0.0D0) THEN   !If not yet in crackband
			D11 = 1.0D0
            STATEV(102) = D11
        ELSEIF ((STATEV(101).gt.0.5D0).and.(STATEV(101).lt.1.5D0).and. 
     1         (STRAIN(1).ge.STATEV(103))) THEN
			STATEV(103) = DABS(STRAIN(1))
		    E110 = STATEV(116)
            eps01 = STATEV(121)
            epsf1 = STATEV(125)
            D11=(sigcr11/E110)*(1.0D0/(epsf1-eps01))*
     1          ((epsf1/STRAIN(1))-1.0d0)
            STATEV(102) = D11
            IF (D11.lt.D1LIM) THEN
                STATEV(102) = D1LIM
                STATEV(101) = 2.0D0
                STATEV(120) = 0.0D0   ! Delete element
            ELSE IF (STATEV(102).gt.1.0D0) THEN
                STATEV(102) = 1.0D0
            END IF
        ELSEIF ((STATEV(101).gt.0.5D0).and.(STATEV(101).lt.1.5D0).and. 
     1         (STRAIN(1).lt.STATEV(103))) THEN
                D11 = STATEV(102)
        ELSEIF (STATEV(101).gt.1.5d0) THEN
            D11 = D1LIM
            STATEV(102) = D1LIM
        END IF
		
		IF (STATEV(102).lt.D1LIM) THEN
               STATEV(102) = D1LIM
               STATEV(101) = 2.0D0
               STATEV(120) = 0.0d0
            END IF
            D11 = STATEV(102)

!C-------------------------------------------------
!C      !For 22-direction
       IF (STATEV(106).lt.0.5D0) THEN   !If not yet in crackband
         D22 = 1.0D0
         STATEV(107) = D22
       ELSEIF ((STATEV(106).gt.0.5D0).and.(STATEV(106).lt.1.5D0)
     1        .and.(STRAIN(2).ge.STATEV(108))) THEN
		 STATEV(108) = DABS(STRAIN(2))
         E220 = STATEV(117)
         eps02 = STATEV(122)
         epsf2 = STATEV(126)
         sigcr22 = STATEV(184)    ! From Hashin intiation stress, not input cr value
         D22=(sigcr22/E220)*(1.0D0/(epsf2-eps02))
     1                     *((epsf2/STRAIN(2))-1.0d0)
		 STATEV(107) = D22
             IF (D22.lt.D2LIM) THEN
				STATEV(107) = D2LIM
                STATEV(106) = 2.0D0
                STATEV(120) = 0.0D0   ! Delete element
             ELSE IF (STATEV(107).gt.1.0D0) THEN
				STATEV(107) = 1.0D0
             END IF
       ELSEIF ((STATEV(106).gt.0.5D0).and.(STATEV(106).lt.1.5D0)
     1        .and.(STRAIN(2).lt.STATEV(108))) THEN
		  D22 = STATEV(107)
       ELSEIF (STATEV(106).gt.1.5d0) THEN
		  D22 = D2LIM
          STATEV(107) = D22
       END IF

	   
	   IF (STATEV(107).lt.D2LIM) THEN
            STATEV(107) = D2LIM
            STATEV(106) = 2.0D0
            STATEV(120) = 0.0d0
         END IF
         D22 = STATEV(107)

	  
	  


!C-------------
!C      !For 12-direction   *************************

       IF (STATEV(106).lt.0.5D0) THEN   !If not yet in crackband
		  D12 = 1.0D0
		  STATEV(169) = D12
       ELSEIF ((STATEV(106).gt.0.5D0).and.(STATEV(106).lt.1.5D0)
     1        .and.(DABS(STRAIN(4)).ge.STATEV(170))) THEN
		      STATEV(170) = DABS(STRAIN(4))
              E120 = STATEV(178)
              eps012 = STATEV(180)
              epsf12 = STATEV(182)
		  sigcr12 = STATEV(185)    ! From Hashin intiation stress, not input cr value
          D12=(sigcr12/E120)*(1.0D0/(DABS(epsf12) - DABS(eps012)))
     1                           *(DABS(epsf12/STRAIN(4)) - 1.0d0)
		  STATEV(169) = D12
             IF (D12.lt.D12LIM) THEN
			   STATEV(169) = D12LIM
               STATEV(106) = 2.0D0
               STATEV(120) = 0.0D0   ! Delete element
             ELSE IF (STATEV(169).gt.1.0D0) THEN
			   STATEV(169) = 1.0D0
             END IF
       ELSEIF ((STATEV(106).gt.0.5D0).and.(STATEV(106).lt.1.5D0)
     1        .and.(DABS(STRAIN(4)).lt.STATEV(170))) THEN
		  D12 = STATEV(169)
       ELSEIF (STATEV(106).gt.1.5d0) THEN
		  D12 = D12LIM
          STATEV(169) = D12      
       END IF


	   IF (STATEV(169).lt.D12LIM) THEN
		STATEV(169) = D12LIM
        STATEV(106) = 2.0D0
		STATEV(120) = 0.0d0
       END IF
       D12 = STATEV(169)



	 	  
!C-------------
!C      !For 23-direction   *************************
       IF (STATEV(106).lt.0.5D0) THEN   !If not yet in crackband
          D23 = 1.0D0
          STATEV(174) = D23
       ELSEIF ((STATEV(106).gt.0.5D0).and.(STATEV(106).lt.1.5D0)
     1        .and.(DABS(STRAIN(6)).ge.STATEV(175))) THEN
          STATEV(175) = DABS(STRAIN(6))
          E230 = STATEV(178)
          eps023 = STATEV(181)
          epsf23 = STATEV(183)
          sigcr23 = STATEV(186)    ! From Hashin intiation stress, not input cr value
          D23=(sigcr23/E230)*((1.0D0/(DABS(epsf23) - DABS(eps023)))
     1                         *(DABS(epsf23/STRAIN(6)) - 1.0d0))
		  STATEV(174) = D23
             IF (D23.lt.D23LIM) THEN
               STATEV(174) = D23LIM
               STATEV(106) = 2.0D0
               STATEV(120) = 0.0D0   ! Delete element
             ELSE IF (STATEV(174).gt.1.0D0) THEN
               STATEV(174) = 1.0D0
             END IF
       ELSEIF ((STATEV(106).gt.0.5D0).and.(STATEV(106).lt.1.5D0)
     1        .and.(DABS(STRAIN(6)).lt.STATEV(175))) THEN
          D23 = STATEV(174)
       ELSEIF (STATEV(106).gt.1.5d0) THEN
          D23 = D23LIM
          STATEV(174) = D23       
       END IF
	   
	   IF (STATEV(174).lt.D23LIM) THEN
            STATEV(174) = D23LIM
            STATEV(106) = 2.0D0
            STATEV(120) = 0.0d0
         END IF
         D23 = STATEV(174)	  
	  
!C-------------

      D33 = 1.0D0







!C-----------------------------------------------------------
!C UPDATE STATE VARIABLES
!C-----------------------------------------------------------
      
      !For 1-direction ***********************
      IF ((STATEV(102).lt.D1LIM).and.(STATEV(101).eq.1.0D0)) THEN
          STATEV(102) = D1LIM
      END IF	  
          D11 = STATEV(102)
	  
      !For 2-direction ***********************	  
      IF ((STATEV(107).lt.D2LIM).and.(STATEV(106).eq.1.0D0)) THEN
	  STATEV(107) = D2LIM
	  END IF	  
          D22 = STATEV(107)
  
      !For 12-plane ***********************	  
      IF ((STATEV(169).lt.D12LIM).and.(STATEV(106).eq.1.0D0)) THEN
	  STATEV(169) = D12LIM
	  END IF	  
          D12 = STATEV(169)

      !For 23-plane ***********************	  
      IF ((STATEV(174).lt.D23LIM).and.(STATEV(106).eq.1.0D0)) THEN
	  STATEV(174) = D23LIM
	  END IF	  
          D23 = STATEV(174)

      
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
	  
!C     !!!! Pass stress to SIG variable
      SIG = STRESS     
      STATEV(100) = 1.0d0
      
	  RETURN
      END




!C!C---------------------------------------------------------------------
!C!C           SUBROUTINES START HERE
!C!C---------------------------------------------------------------------

!C!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE calc_stress(Q_IN, EPS_IN, SIG_OUT)
      REAL*8, INTENT(IN) :: Q_IN(6,6)
      REAL*8, INTENT(IN) :: EPS_IN(6)
      REAL*8, INTENT(OUT) :: SIG_OUT(6)

      SIG_OUT = MATMUL(Q_IN, EPS_IN)

      RETURN
      END





!C!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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


!C!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C!Calculates transversely isotropic stiffness matrix
      SUBROUTINE calc_stiffness_transiso(E11c, E22c, v12c, 
     1                                 v23c, G12c, G23c, Q)
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

!C!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE calc_eqstrain(EP, j2eps)
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
           
      !write(*,*) 'eps_j2_tow_fun', j2eps
	  RETURN
      END

!C!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE calc_ccm(E11f, E22f, v12f, v23f, G12f, G23f, 
     1    Em, vm, Vf, E11c, E22c, v12c, v23c, G12c, G23c)
      
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


!C!!--------------------------------------------------------------------
!C!!!    SUBROUTINE FOR CCM

      SUBROUTINE CCM_CALC_tow(CCMinput_fun,CCMoutput_fun)
        
      IMPLICIT NONE      
      REAL*8:: CCMinput_fun(8)
      REAL*8:: CCMoutput_fun(8)
      REAL*8:: E11c, E22c, E33c, G12c, G13c, G23c
      REAL*8:: Ef1, Ef2, Gf23, Gf12, vf12, Volf, Volm, Em, vm
      REAL*8:: vf23, vf21, Gm, gamma, delta, etaf, etam, eta4
      REAL*8:: v12c, v23c
      
      Ef1 = CCMinput_fun(1)
      Ef2 = CCMinput_fun(2)
      Gf12 = CCMinput_fun(3)
      Gf23 = CCMinput_fun(4)
      vf12 = CCMinput_fun(5)
      Volf = CCMinput_fun(6)
      Em = CCMinput_fun(7)
      vm = CCMinput_fun(8)
	  
      Volm = 1.0D0 - Volf
      vf23 = ((Ef2/(2.0D0*Gf23)) - 1.0D0)
      vf21 = Ef2*vf12/Ef1
      Gm = Em/(2.0D0*(1.0D0 + vm))
      gamma=2.0d0*vf21*Em*(1.0d0-Volf)*(vf12-vm)/
     1         (Ef2*(1.0d0+vm)*(1.0d0+Volf*(1.0d0-2.0d0*vm))
     1        + Em*(1.0d0 - vf23 -2.0D0*vf12*vf21)*(1.0d0-Volf))
      delta = 2.0d0*Ef2*vm*Volf*(vm-vf12)/(Ef2*(1.0d0+vm)
     1        *(1.0d0+Volf*(1.0d0-2.0d0*vm)) 
     2        + Em*(1.0d0-vf23-2.0D0*vf12*vf21)*(1.0D0-Volf))
      E11c = Ef1*(1.0D0+gamma)*Volf + Em*(1.0D0 + delta)*(1.0D0 - Volf)	 
      etaf = (Ef1*Volf+((1.0D0-vf12*vf21)*Em + vm*vf21*Ef1)
     1       *(1.0D0-Volf))/(Ef1*Volf + Em*(1.0D0-Volf))	 
      etam = (((1.0D0- (vm**2.0D0))*Ef1+ vm*vf12*Em)*Volf + Em*Volm)/
     2       (Ef1*Volf+Em*(1.0D0-Volf))	 
      E22c = 1.0D0/((etaf*Volf/Ef2)+(etam*(1.0D0 - Volf)/Em))
     
	 
      G12c =Gm*(((Gm+Gf12)-Volf*(Gm-Gf12))
     1      /((Gm+Gf12) + Volf*(Gm-Gf12)))
      
      eta4 = (3.0D0 - 4.0D0*vm + (Gm/Gf23))/(4.0D0*(1.0D0 - vm))
      
      
      G23c = (Volf+eta4*(1.0D0-Volf))
     1       /((Volf/Gf23)+(eta4/Gm)*(1.0D0-Volf))
	  
      v12c = (((1.0D0-Volf)*(1.0D0 - vf23 - 2.0D0*vf12*vf21))*vm*Em 
     1      +(vm+Volf*(2.0D0*vf12-vm)+((vm**2.0D0)*
     2      (1.0D0-2.0D0*Volf*vf12-Volf)))*Ef2)
     3      /(((1.0D0-Volf)*(1.0D0-vf23-2*vf12*vf21))*Em
     4      +(1.0D0+Volf+(1.0D0-Volf)*vm - 2.0D0*Volf*(vm**2.0D0))*Ef2)

      v23c = (E22c/(2.0D0*G23c)) - 1.0D0
	 

      CCMoutput_fun(1) = E11c
      CCMoutput_fun(2) = E22c
      CCMoutput_fun(3) = v12c
      CCMoutput_fun(4) = v23c
      CCMoutput_fun(5) = G12c
      CCMoutput_fun(6) = G23c
      CCMoutput_fun(7) = Gm
      CCMoutput_fun(8) = eta4

      RETURN
      END

	  
	  
	  
	  

!C!C-------------------------------------------------------------------
!C!C-------------------------------------------------------------------
!C!C SUBROUTINES TO INVERT THE MATRIX
!C!C-------------------------------------------------------------------
!C!C-------------------------------------------------------------------
      subroutine t_matrixInverse(a,y,n,np)

      !taken from http://www.mpi-hd.mpg.de/astrophysik/HEA/internal/Numerical_Recipes/f2-3.pdf
      !a goes in, y=inverse(a) goes out


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

!C!---------------------------------------------------------------------------------
!C!!
!C!---------------------------------------------------------------------------------
!C! same as above, quad precision

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

!C!--------------------------------------------------------------------------------------------------
!C!
!C!--------------------------------------------------------------------------------------------------

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

!C!      write(*,*) 'A='
!C!      write(*,*) A
!C!      write(*,*) 'b='
!C!      write(*,*) b
!C!      write(*,*) '---'


      call t_ludcmp(A,n,n,indx,d,flag)
      if (flag.eq.1) then
      call t_lubksb(A,n,n,indx,b)
      endif

      return
      end


!C!--------------------------------------------------------------------------------------------------
!C!
!C!--------------------------------------------------------------------------------------------------

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

!C!--------------------------------------------------------------------------------------------------
!C!
!C!--------------------------------------------------------------------------------------------------

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

!C!--------------------------------------------------------------------------------------------------
!C!
!C!--------------------------------------------------------------------------------------------------
!C!same as above with quad precision

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

!C!--------------------------------------------------------------------------------------------------
!C!
!C!--------------------------------------------------------------------------------------------------

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


!C!--------------------------------------------------------------------------------------------------
!C!
!C!--------------------------------------------------------------------------------------------------
!C!same as above with quad precision

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
