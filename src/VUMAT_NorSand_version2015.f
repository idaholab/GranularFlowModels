! coded by Wencheng Jin, wencheng.jin@inl.gov
! Sequence of input variable for material constant
!     xGamma  = props(1)
!     xlambda = props(2)
!     xMtc    = props(3)
!     xN      = props(4)
!     xH      = props(5)
!     chi     = props(6)
!     GIr     = props(7)
!     xnu     = props(8)
!     nSubsteps =  props(9) :: = 0, no error control, forward Euler
!                           :: >=1, error control, modefied forward Euler
!     initial void ratio = props(10)
!     Drained/iUndrained = props(11) :: = 0, drained
!                                    :: = 1, iUndrained
      subroutine vumat(
! Read only (unmodifiable)variables -
     &  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     &  stepTime, totalTime, dt, cmname, coordMp, charLength,
     &  props, density, strainInc, relSpinInc,
     &  tempOld, stretchOld, defgradOld, fieldOld,
     &  stressOld, stateOld, enerInternOld, enerInelasOld,
     &  tempNew, stretchNew, defgradNew, fieldNew,
! Write only (modifiable) variables -
     &  stressNew, stateNew, enerInternNew, enerInelasNew)
      include 'vaba_param.inc' !assign precision
      character*80 cmname

C       double precision::stepTime, totalTime
C       integer nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal
!      real stepTime, totTime, dt
      dimension::props(nprops),density(nblock),coordMp(nblock,3),
     &  charLength(nblock), strainInc(nblock,ndir+nshr),
     &  relSpinInc(nblock,nshr), tempOld(nblock),
     &  stretchOld(nblock,ndir+nshr),
     &  defgradOld(nblock,ndir+nshr+nshr),
     &  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     &  stateOld(nblock,nstatev), enerInternOld(nblock),
     &  enerInelasOld(nblock), tempNew(nblock),
     &  stretchNew(nblock,ndir+nshr),
     &  defgradNew(nblock,ndir+nshr+nshr),
     &  fieldNew(nblock,nfieldv),
     &  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     &  enerInternNew(nblock), enerInelasNew(nblock)

C       double precision, intent(in) :: props
      
      !implicit none
      dimension:: DE(6,6), subdStrain(6),dFdSig(6), dstrain(6),dEpsE(6),
     &  stress(6), dStress(6), STRESSTrial(6), stress1(6), stress2(6), 
     &  dStress2(6),DDSDDE(6,6)
C       double precision ::  xGamma, xlambda, xMtc, xN, xH, chi, GIr, xnu
C       double precision ::  xG, xK, void, pImage, FTOL, STOL,dEpsV,dEpsD
C       double precision :: T, dT, porosity, PoreP, xM_tc_image
C       double precision :: pTrial, qTrial, etaTrial, yield, p,q,eta
C       double precision :: p1,q1,eta1,pImage1, void1, xMImage1
C       double precision :: xM1,  dpImage1,dpImage2
C       double precision :: theta,xJ3,xJ2,cos3Theta,xM,yieldTrial
C       double precision :: psi_image, chi_image, xMImage,subdEpsV,subdEpsD
C       double precision :: stress_power
C       double precision :: temp1, temp2, temp3, cosBeta
C       double precision :: dStrainE(6),dStress1(6)
C       double precision :: A, b, e_min, p_ref
C       integer :: km, ic, iflag, iSubsteps, iUndrain
C       integer :: iflag_yield, iflag_csd,ITER,nSubsteps
C       control iteration and error
      PARAMETER (FTOL=1.D-6,ITER=5,STOL=1.D-4)

      ! get model parameters
      xGamma= props(1)
      xlambda= props(2)
      xMtc   = props(3)
      xN    = props(4)
      xH    = props(5)
      chi   = props(6)
      GIr   = props(7)
      xnu    = props(8)
      nSubsteps = props(9)
      iUndrain = props(11)

      !input chi is the result from triaxial
      chi_image = chi / ( 1-xlambda*chi/xMtc)

      !write(6,*) props

      do km = 1, nblock
      	do ic = 1,ndir
            if (stressOld(km,ic) > 0) then
                stress(ic)  = -10.00001
            else
                stress(ic)  = stressOld(km,ic)
            endif
            dstrain(ic) = strainInc(km,ic)
        enddo

        if (nshr==1) then
            stress(4)= stressOld(km,4)
            dstrain(4) = strainInc(km,4)
            stress(5) = 0
            stress(6) = 0
            dstrain(5) = 0
            dstrain(6) = 0
        elseif (nshr==3) then
            do ic=4,6
              stress(ic)  = stressOld(km,ic)
              dstrain(ic) = strainInc(km,ic)
            enddo
        endif

        void   = stateOld(km,1)
        pImage = stateOld(km,2)
        xMImage = stateOld(km,3)
        PoreP  = -stateOld(km,8)
        iflag_csd = stateOld(km,7)

        !write(6,*) stateOld(km,1), stateOld(km,2), stateOld(km,3)
        !write(6,*) stateOld(km,8), stateOld(km,7)

        if (void == 0) then
            void = props(10)
            pImage = 0
            !pImage = props(12)
            !xMImage = props(13)
        endif

        if (iUndrain ==1) then !for iUndrained case
          stress(1) = stress(1) - PoreP
          stress(2) = stress(2) - PoreP
          stress(3) = stress(3) - PoreP
          porosity = void/(1+void)
          call getDevVolStrain(dstrain,dEpsV,dEpsD)
          PoreP = PoreP + 2E9*dEpsV/porosity ! change of pore water pressure
        endif

        !write(6,*) stress, dstrain, pImage, xMImage

        ! initiation of substepping
        T = 0.
        dT = 1. !start with the whole step, reduce it if the relative error is larger than tolerance

        !substepping
        do while (T < 1.0 )
          if ( T+ dT > 1. ) dT = 1. - T
C           do ic = 1, 6
C             subdStrain(ic) = dT*dstrain(ic)
C           end do
          subdStrain(1) = dT*dstrain(1)
          subdStrain(2) = dT*dstrain(2)
          subdStrain(3) = dT*dstrain(3)
          subdStrain(4) = dT*dstrain(4)
          subdStrain(5) = dT*dstrain(5)
          subdStrain(6) = dT*dstrain(6)

          call getPandQ (stress,p,q,eta)
          !Get the elasticitity modulus
C           e_min = 0.29268433
C           A = 416.5729504
C           b = 0.463124486
C           p_ref = 100E3
C           xG = 0.8*A*p_ref/(void-e_min)*(p/p_ref)**b
C           xK = (2.*(1.+xnu)) / (3.*(1.-2.*xnu))*xG
          xG = GIr * p
          xK = (2.*(1.+xnu)) / (3.*(1.-2.*xnu))*xG

          call ZERO2(DE,6)
          call elasticStiffnessMatrixNorSand (xG, xK, DE)
          call MatVec (DE,6, subdStrain, 6, dStress)
          call AddVec (stress, dStress, 1.d0, 1.d0, 6, STRESSTrial)
          call getPandQ (STRESSTrial,pTrial,qTrial,etaTrial)

          if ( pImage == 0 ) then
            ! Get critical state stress ratio xM(theta) with initial stress state and void ratio
            call getMlode (stress,xMtc,theta,xJ3,xJ2,cos3Theta,xM)
            ! Get xMImage and pImage from the nonlinear equation
            call getP_M_Image (void, xN, xM, xMtc, chi, xGamma, xlambda,
     &                 p, q, pImage, psi_image, chi_image, xMImage)
            !write(6,*) pImage, xMImage
          endif

          call getYieldFunctionNorSand (p,q,pImage,xMImage,yield)
          call getYieldFunctionNorSand (pTrial,qTrial,pImage,
     &                              xMImage,yieldTrial)
          if ( yieldTrial < FTOL ) then ! elastic
            call getDevVolStrain(subdStrain,subdEpsV,subdEpsD)
            void1 = void + subdEpsV * ( 1. + void)

            psi_image = void1 - xGamma + xlambda*log(pImage/1000)

            xM_tc_image = xMtc*(1 - chi_image*xN*abs(psi_image)/xMtc)  !image critical stress ratio
C             p_cap = pImage/exp(-chi_image/xM_tc_image*psi_image)
            !write(6,*), pTrial, p_cap, pImage 
            !if ( pTrial > p_cap ) then
              !DDSDDE=DE
              void = void1
              stress = STRESSTrial
              p = pTrial
              q = qTrial
              eta = qTrial/pTrial
              iflag_yield = 0
C             else ! plastic unloading
C               !if (km == 247) then
C                 !write(6,*), km, totalTime, stepTime
C                 !write(6,*), pTrial, p_cap, pImage
C                 !write(6,*), psi_image
C               !endif:,1
C               iflag = 1
C               iflag_yield = 2
C               iSubsteps = 1
C               ! dT is adjusted, the whole step is re-iterated
C               do while (iflag == 1)
C                 iSubsteps = iSubsteps + 1
C                 call getDevVolStrain(subdStrain,subdEpsV,subdEpsD)
C                 void1 = void + subdEpsV * ( 1. + void)
C                 call getElastoPlasticStiffnessUnloading(void1, stress,
C      &               subdStrain, p, q, pImage,xN, xH, xMtc, chi, xGamma,
C      &               xlambda, DE, DDSDDE, dpImage1,dStrainE)
C                 call MatVec (DDSDDE,6, subdStrain, 6, dStress)
C                 call MatVec (DE,6, dStrainE, 6, dStress1)
C                 call AddVec (stress, dStress, 1.d0, 1.d0,6, stress1)
C                 call getPandQ (stress1,p1,q1,eta1)
C                 call getMlode (stress1,xMtc,theta,xJ3,xJ2,cos3Theta,xM1)
C                 !State of the stress does not stay on the outer yield surface
C                 !pImage1 = pImage + dpImage1
                

C                 pImage0 = pImage - 1
C                 pImage1 = pImage
C                 do while (abs(pImage1 - pImage0)>1e-4)
C                   pImage0 = pImage1
C                   psi_image = void1-xGamma + xlambda*log(pImage0/1000)
C                   if (psi_image<0) then
C                     sign_ = -1
C                   else
C                     sign_ = 1
C                   endif
C                   xM_tc_image = xMtc - chi_image*xN*abs(psi_image)  !image critical stress ratio
C                   f_x = pImage0-p1*exp(-chi_image/xM_tc_image*psi_image)
C                   df_x = 1-p1*exp(-chi_image/xM_tc_image*psi_image)
C      &                 *(-chi_image*xlambda/pImage0*xM_tc_image
C      &                 -chi_image**2*xN*sign_*xlambda/pImage0*
C      &                 psi_image)/xM_tc_image**2
C                   pImage1 = pImage0 - f_x/df_x
C                 enddo
C                 psi_image = void1 - xGamma + xlambda*log(pImage1/1000)
C                 xM_tc_image = xMtc*(1 - chi_image*xN*abs(psi_image)/xMtc)
C                 yield = p1-pImage1/exp(-chi_image/xM_tc_image*psi_image)

C                 !write(6,*), dpImage1, pImage1-pImage
C                 !WRITE (6,'(''yield-value'')')
C                 !write(6,*),yield

C                 call updateM_Image(void1, xN, xM1, xMtc, chi, xGamma,
C      &                xlambda, p1, q1, pImage1, xMImage1)
                
C                 if (nSubsteps /= 0) then
C                   call getElastoPlasticStiffnessUnloading(void1,stress1,
C      &               subdStrain, p1, q1, pImage1,xN, xH, xMtc,chi,xGamma,
C      &               xlambda, DE, DDSDDE, dpImage2,dStrainE)
C                   call MatVec (DDSDDE,6, subdStrain, 6, dStress2)
C                   !call MatVec (DE,6, dStrainE, 6, dStress2)
C                   call errorEuler(stress, dStress1, dStress2, pImage,
C      &                 dpImage1, dpImage2, stress2, dT, iflag, STOL)
C                   if (iflag==1 .and. iSubsteps<nSubsteps) then 
C                     do ic = 1, 6
C                       subdStrain = dT*dstrain(ic)
C                     end do
C                   else
C                     iflag = 0
C                   endif
C                 else if ( nSubsteps == 0)  then
C                   iflag = 0
C                 endif
C               end do
C               if ( nSubsteps == 0) then
C                 stress = stress1
C                 pImage = pImage1
C                 void = void1
C                 p = p1
C                 q = q1
C                 eta = eta1
C               else
C                 stress = stress1
C                 void = void1
C                 call getPandQ (stress,p,q,eta)
C                 pImage = pImage1! + (dpImage1+dpImage)/2
C                 call getMlode (stress,xMtc,theta,xJ3,xJ2,cos3Theta,xM)
C                 call updateM_Image(void, xN, xM, xMtc, chi, xGamma,
C      &                xlambda, p, q, pImage, xMImage)
C               end if
C             endif
          else  ! plastic
            !check whether this increment is from elastic regime or not
            !if it is, update stress and void ratio, make them stay on yield surfaces
c            if ( yield < -FTOL ) then
c                temp1 = 0
c                temp2 = 1
c              call getElasticPartPegasus(yield,yieldTrial, temp1, temp2,
c     &            stress, dStress, pImage, xMImage, DE, subdStrain, void)
c              call getPandQ (stress,p,q,eta)
cc            else if ( abs(yield) < FTOL ) then
c              call getdFdSig(void, stress, p, q, pImage, xN, xMtc, chi,
c     &            xGamma, xlambda, dFdSig)
c              temp1 = 0
c              temp2 = 0
c              do ic=1,6
c              temp3 = 0
c                temp1 = temp1 - dStress(ic)*dFdSig(ic) !This part has some issues
c                temp2 = temp2 + dStress(ic)**2
c                temp3 = temp3 + dFdSig(ic)**2
c              end do
c              cosBeta = temp1/temp2/temp3
c              if ( cosBeta <= -0.1 ) then
c                call getElasticPartUnload(yield, stress, dStress,
c     &           pImage, xMImage, DE, subdStrain, void, FTOL)
c                call getPandQ (stress,p,q,eta)
c              !else - > pure plastic loading, starts from yield surface and will end on yield surface
c              end if
c            end if
             !if ( iflag_csd == 2 ) then
             !WRITE (6,'(''stress'')')
                !write(6,*), stress
                !WRITE (6,'(''stateV'')')
                !write(6,*), stateOld(km,1),stateOld(km,2),stateOld(km,3)
                !write(6,*), stateOld(km,4),stateOld(km,5),stateOld(km,6)
                !write(6,*), stateOld(km,7),stateOld(km,8)
                !call XPLB_EXIT
            !endif
            iflag_yield = 1
            iflag = 1
            iSubsteps = 1
            ! modified Eeuler integration scheme with error control, if relative error is not within tolerance
            ! dT is adjusted, the whole step is re-iterated
        	  do while (iflag == 1)
              iSubsteps = iSubsteps + 1
              call getDevVolStrain(subdStrain,subdEpsV,subdEpsD)
              void1 = void + subdEpsV * ( 1. + void)

              call getElastoPlasticStiffnessMatrix (void1, stress,
     &            subdStrain, p, q, pImage,xN, xH, xMtc, chi, xGamma,
     &            xlambda, DE, DDSDDE, dpImage1)
              call MatVec (DDSDDE,6, subdStrain, 6, dStress)
              call AddVec (stress, dStress, 1.d0, 1.d0,6, stress1)
              call getPandQ (stress1,p1,q1,eta1)
              !pImage1 = pImage + dpImage1
              call getMlode (stress1,xMtc,theta,xJ3,xJ2,cos3Theta,xM1)
                
              pImage1 = p1*exp(eta1/xMImage-1)
          
              call updateM_Image(void1, xN, xM1, xMtc, chi, xGamma,
     &                xlambda, p1, q1, pImage1, xMImage1)

              if (nSubsteps /= 0) then
                call getElastoPlasticStiffnessMatrix (void1, stress1, 
     &           subdStrain, p1, q1, pImage1, xN, xH, xMtc, chi,
     &               xGamma, xlambda, DE, DDSDDE, dpImage2)
                call MatVec (DDSDDE,6, subdStrain, 6, dStress2)
                !pImage2 = pImage1 + dpImage2
                !call getImagePressure(p1, q1, xMImage1, pImage1)
                call errorEuler(stress, dStress, dStress2, pImage,
     &                 dpImage1, dpImage2, stress2, dT, iflag, STOL)
                if (iflag==1 .and. iSubsteps<nSubsteps) then 
                  do ic = 1, 6
                    subdStrain = dT*dstrain(ic)
                  end do
                else
                  iflag = 0
                endif
              else if ( nSubsteps == 0)  then
                iflag = 0
              endif
            end do
            if ( nSubsteps == 0) then
              stress = stress1
              pImage = pImage1
              void = void1
              p = p1
              q = q1
              eta = eta1
            else
              stress = stress1
              void = void1
              call getPandQ (stress,p,q,eta)
              pImage = pImage1
              xMImage = xMImage1

              ! With updated stress and void ratio, update pImage and xMImage
!              call getMlode (stress,xMtc,theta,xJ3,xJ2,cos3Theta,xM)             
!              call getP_M_Image (void, xN, xM, xMtc, chi, xGamma, xlambda,
!     &                p, q, pImage, psi_image, chi_image, xMImage)
            end if
            !WRITE(222,*) yieldsave, yield
          end if ! end of plastic phase
          T = T + dT
        end do ! end substepping
        !check whether the updated stress is on the yield surface or not, if not, perform correction
         call getYieldFunctionNorSand (p,q,pImage,xMImage,yield)
C             WRITE(*,*) 'yield=',yield
C             if ( abs(yieldsave) < abs(yield) ) then
C               icount = 1
C               do while ( abs(yield) > FTOL .and. icount < ITER) 
C                 icount = icount + 1
C                 call getdFdSig(void, stress, p, q, pImage, xN, xMtc, chi,
C      &                 xGamma, xlambda, dFdSig)
C                 temp1 = 0
C                 do ic = 1, 6
C                   temp1 = temp1 + dFdSig(ic)**2
C                 end do
C                 dlambda = yield/temp1
C                 do ic = 1, 6
C                   stress(ic) = stress(ic) - dlambda*dFdSig(ic)
C                 end do
C                 call getPandQ (stress,p,q,eta)
C                 !call getMlode(stress,xMtc,theta,xJ3,xJ2,cos3Theta,xM)
C                 !call getxMImage(void, xN, xM, xMtc, chi, xGamma, xlambda,
C !     &                   p, pImage, psi_image, chi_image, xMImage)
C                 call getYieldFunctionNorSand (p,q,pImage,xMImage,yield)
C               end do
C             else 
C               icount =1
C               do while ( abs(yield) > FTOL .and. icount < ITER) 
C                 icount = icount + 1
C                 !perform stress correction iteratively
C                 call StressCorrection(void, stress, p, q, pImage,
C      &               xN, xH, xMtc, chi, xGamma, xlambda, DE, yield)
C                 call getPandQ (stress,p,q,eta)
C                 call getMlode(stress,xMtc,theta,xJ3,xJ2,cos3Theta,xM)

C                 call getP_M_Image (void, xN, xM, xMtc, chi, xGamma, 
C      &             xlambda, p, q, pImage, psi_image, chi_image, xMImage)
C                 call getYieldFunctionNorSand (p,q,pImage,xMImage,yield)
C               end do
C             end if

      !updated vector of additional state variables to abaqus statev vector
        stateNew(km,1) = void
        stateNew(km,2) = pImage
        stateNew(km,3) = xMImage
        stateNew(km,4) = p
        stateNew(km,5) = q
        if (totalTime == 0.0) then
          stateNew(km,6) = dEpsV
        else
          stateNew(km,6) = stateOld(km,6)+dEpsV !xMImage - eta
        endif
        stateNew(km,7) = iflag_yield
        stateNew(km,8) = xG/1E6!-PoreP

        if (iUndrain ==1) then !for iUndrained case
          stress(1) = stress(1) + PoreP
          stress(2) = stress(2) + PoreP
          stress(3) = stress(3) + PoreP
        endif

        !write(6,*), stateNew(km,2),stateNew(km,3),stateNew(km,4)

        do ic = 1,ndir
            if (stress(ic) > 0) then
                stressNew(km,ic)  = -10.00001
            else
                stressNew(km,ic)  = stress(ic)
            endif
        end do
        if (nshr==1) then
            stressNew(km,4) = stress(4)
        elseif (nshr==3) then
          stressNew(km,4) = stress(4)
          stressNew(km,5) = stress(5)
          stressNew(km,6) = stress(6)
        endif
        ! update internal stress:

        stress_power = 0.5 * (
     &        ( stressOld(km,1)+stressNew(km,1) )*strainInc(km,1) +
     &        ( stressOld(km,2)+stressNew(km,2) )*strainInc(km,2) +
     &        ( stressOld(km,3)+stressNew(km,3) )*strainInc(km,3))
        if(nshr == 1) then
            stress_power = stress_power + 0.5*(
     &            2*(stressOld(km,4)+stressNew(km,4))*strainInc(km,4))
        else
            stress_power = stress_power + 0.5*(
     &            2*(stressOld(km,4)+stressNew(km,4))*strainInc(km,4) +
     &            2*(stressOld(km,5)+stressNew(km,5))*strainInc(km,5) +
     &            2*(stressOld(km,6)+stressNew(km,6))*strainInc(km,6) )
        endif
        enerInternNew(km) = enerInternOld(km) 
     &      + stress_power/density(km)
            
      enddo

      return
      
      end subroutine vumat
                                                  
 ! -----------------------------------------------------------------------
 !                          INTERNAL SUBROUTINES
 ! -----------------------------------------------------------------------

      subroutine getElastoPlasticStiffnessMatrix(void,stress,subdStrain,
     &            p, q, pImage, xN, xH, xMtc, chi_tc, xGamma,
     &            xlambda, DE, DDSDDE, dpImage)
        include 'vaba_param.inc'
C         implicit none
        dimension:: stress(6), subdStrain(6), DE(6,6), DDSDDE(6,6),
     &    dJdSig(6), dxJ3dSig(6),dPdSig(6),dQdSig(6),dThetadSig(6),
     &    dFdSig(6), C_dFdSig(6), dStress(6), DE_dFdSig(6), 
     &    DE_dFdSig_dFdSig(6,6), DP(6,6)

C         double precision, intent(inout) :: stress(6)
C         double precision, intent(in) :: void, subdStrain(6)
C         double precision, intent(in) :: p,q, pImage, xN, xH, xMtc,chi_tc
C         double precision, intent(in) :: xGamma, xlambda, DE(6,6)
C         double precision, intent(out):: DDSDDE(6,6), dpImage
C         double precision :: dJdSig(6), dxJ3dSig(6),dPdSig(6),dQdSig(6)
C         double precision :: dThetadSig(6),dFdSig(6)!, dZetadSig(6) 
C         double precision :: dFdpImage, chi_image, xMImage, p_imageMax
C         double precision :: xJ2, xJ3, cos3Theta, theta, xM, psi_image
C         double precision :: dFdP, dFdQ, dFdM, dMdTheta, Zeta, dMdpImage
C         double precision :: C_dFdSig(6), dFdSig_C_dFdSig, dpImagedEpsQ
C         double precision :: xM_image_tc, pi, denominater, phi, dStress(6)
C         double precision :: DE_dFdSig(6), DE_dFdSig_dFdSig(6,6), DP(6,6)
C         double precision :: dFdq_,psi
C         integer :: I, J, Z

        PARAMETER (pi=3.14159265359D0)

        call ZERO1(dJdSig,6)
        call ZERO1(dxJ3dSig,6)
        call ZERO1(dPdSig,6)
        call ZERO1(dQdSig,6)

        call ZERO1(dThetadSig,6)
C         call ZERO1(dZetadSig,6)

        xJ2 = 0.
        xJ3 = 0.
        cos3Theta = 0.
        theta  = 0.
        xM = 0.

        psi_image = 0.
        chi_image = 0.
        xMImage = 0.

        dFdP = 0.
        dFdQ = 1.
        dFdM = 0
        dMdTheta = 0
        dFdpImage = 0.
        dMdpImage = 0
        Zeta = 0
           
        call getInvariants (stress, xJ2, dJdSig, xJ3, dxJ3dSig)
        call getMlode (stress,xMtc,theta,xJ3,xJ2,cos3Theta,xM)
        call getMaximuxMImagePressure (void, xN, xM, xMtc, chi_tc, xGamma,
     &            xlambda, p, pImage, p_imageMax, xMImage, xM_image_tc, 
     &            psi_image, chi_image)

!---------calculate dFdSigma --------------

        dFdP=xMImage - q/p;
        do I = 1, 3
               dPdSig(I) = 1./3.
               if (q == 0 ) then
                 dQdSig(I)=0
               else
                 dQdSig(I) =  3. / 2. / q * (-stress(I) - p)
               end if
        end do
        do I = 4, 6
               dPdSig(I) = 0. 
               if (q == 0 ) then
                 dQdSig(I)=0
               else
                 dQdSig(I) =  3. / q * (-stress(I) )
               endif
        end do

        dFdM  = -q/xMImage
        dMdTheta = (1 - chi_image*xN*abs(psi_image)/xMtc)
     &            *(3./2.* xMtc**2/( 3 + xMtc )) * sin( 3*theta/2 + pi/4) 
        if ( -0.01 < sqrt(xJ2) < 0.01) then ! avoid dividing by 0
            do I = 1, 6
                dThetadSig(I) = 0.
            end do
        else
            if (-1.e-16 < cos3Theta < 1e-16) then  ! avoid dividing by 0
              if (cos3Theta == 0) then 
                  do I = 1, 6
                      dThetadSig(I) = 0.
                  end do
              else
                  do I = 1, 6
                  ! Need to confirm the correct sign + or - ? 
                    dThetadSig(I) = sqrt(3.)/2./ cos3Theta/sqrt(xJ2**3.)
     &                            * (dxJ3dSig(I) - 3./2.*xJ3/xJ2*dJdSig(I))
                  end do
              end if
            end if
        end if

c        do I = 1, 6
c          if ( psi_image > 0 ) then 
c              dZetadSig(I) =  - *chi_image*xN/xMtc*xlambda/p*dPdSig(I)
c          else
c              dZetadSig(I) =  + *chi_image*xN/xMtc*xlambda/p*dPdSig(I)
c          end if
c        end do

        do I = 1, 6
          dFdSig(I) = dFdP*dPdSig(I) + dFdQ*dQdSig(I) 
     &          + dFdM*dMdTheta*dThetadSig(I)
        end do
        !dFdSig(1) = -dFdSig(1)

!---------calculate dFdpImage--------------

        dFdpImage = -xMImage*p/pImage 
        if ( psi_image > 0 ) then
          dMdpImage = - xM*chi_image*xN/xMtc*xlambda/pImage
        else
          dMdpImage =  xM*chi_image*xN/xMtc*xlambda/pImage
        end if
        dFdpImage = dFdpImage + dFdM*dMdpImage

!---------calculate dPidEpsQ for plasticity (hardening/softening) --------------

        psi = void - xGamma + xlambda*log(p/1000)

        dpImagedEpsQ=(xH-50*psi)*xMImage/xM_image_tc
     &               *(p_imageMax-pImage)*(p/pImage)
        !dpImagedEpsQ = xH*xMImage/xM_image_tc*(p_imageMax-pImage)*(p/pImage)
        
        call getdFdq_(dFdSig, dFdq_)

        do I = 1, 6
          C_dFdSig(I) = 0.
          do J = 1, 6
            C_dFdSig(I) = C_dFdSig(I) + dFdSig(J)*DE(J,I)
          end do
        end do
                
c        call MatVec (DE, 6, dFdSig, 6, C_dFdSig)

        dFdSig_C_dFdSig = 0
        do I = 1,6
          dFdSig_C_dFdSig = dFdSig_C_dFdSig + dFdSig(I)*C_dFdSig(I)
        end do

        denominater = dFdSig_C_dFdSig - dFdpImage*dpImagedEpsQ*dFdq_

        do I = 1, 6
          DE_dFdSig(I) = 0.
          do J = 1, 6
            DE_dFdSig(I) = DE_dFdSig(I)+DE(I,J)*dFdSig(J)
          end do
        end do

        do I = 1, 6
          do J = 1, 6
            DE_dFdSig_dFdSig(I,J) = DE_dFdSig(I) * dFdSig(J)
          end do
        end do
        
        do I = 1, 6
          do J = 1, 6
            DE_dFdSig_dFdSig(I,J) = DE_dFdSig_dFdSig(I,J)/denominater
          end do
        end do

        do I = 1, 6
          do J = 1, 6
            DP(I,J) = 0.
            do Z = 1, 6
      DP(I,J) = DP(I,J) + DE_dFdSig_dFdSig(I,Z) * DE(Z,J)
            end do
          end do
        end do

        ! obtain elastoplasticity stiffness matrix 

        do I = 1, 6
          do J = 1, 6
            DDSDDE(I,J) = DE(I,J) - 0.5*(DP(I,J)+DP(J,I))
          end do
        end do

        ! update the hardening softening parameter pImage with multipler
        call MatVec (DE,6, subdStrain, 6, dStress)

        phi = 0
        do I = 1,6
          phi = phi + dFdSig(I)*dStress(I)/denominater
        end do

        !Note dFdQ = 1, dEpsQ = phi*dFdQ = phi

        dpImage = -phi*dpImagedEpsQ

      end subroutine getElastoPlasticStiffnessMatrix
      
C       subroutine getPlasticStiffnessMatrixNorSand (DElast,
C      &           dFdSig,dFdPi,dPidEpsD, dFdQ, DPlast)

C            implicit none
C            double precision, intent(in) :: DElast(6,6)
C            double precision, intent(in) ::dFdSig(6),dFdPi,dPidEpsD,dFdQ
C            double precision, intent(out) :: DPlast(6,6)
C            double precision :: DE_dFdSig(6), DE_dFdSig_dFdSig(6,6)
C            double precision :: dFdSig_DE(6), dFdSig_DE_dFdSig
C            integer :: i, j, z
C            ! Evaluate De a 
C            do i = 1, 6
C                DE_dFdSig(i) = 0.
C                do j = 1, 6
C                    DE_dFdSig(i) = DE_dFdSig(i)+DElast(i,j)*dFdSig(j)
C                end do
C            end do
 
C            ! Evaluate (De a) bt
C            do i = 1, 6
C                do j = 1, 6
C                    DE_dFdSig_dFdSig(i,j) = DE_dFdSig(i) * dFdSig(j)
C                end do
C            end do
 
C            ! Evaluate b De
C            do i = 1, 6
C                dFdSig_DE(i) = 0.
C                do j = 1, 6
C                    dFdSig_DE(i) = dFdSig_DE(i)+dFdSig(j)*DElast(j,i)
C                end do
C            end do
 
C            !  Evaluate (b De) a
C            dFdSig_DE_dFdSig = 0.
C            do i = 1, 6
C             dFdSig_DE_dFdSig=dFdSig_DE_dFdSig+dFdSig_DE(i)*dFdSig(i)
C            end do
          
 
C            ! get Dplast
C            do i = 1, 6
C                do j = 1, 6
C                    DPlast(i,j) = 0.
C                    do z = 1, 6
C       DPlast(i,j) = DPlast(i,j)+DE_dFdSig_dFdSig(i,z) * DElast(z,j)
C                    end do
C                end do
C            end do          
C       end subroutine getPlasticStiffnessMatrixNorSand
      

C       subroutine StressCorrection(void, stress, p, q, pImage,
C      &    xN, xH, xMtc, chi_tc, xGamma, xlambda, DE, yield )
C         implicit none
C         double precision, intent(inout) :: void, stress(6), pImage
C         double precision, intent(in) :: p,q,xN,xH, xMtc,chi_tc,xGamma
C         double precision, intent(in) :: DE(6,6),xlambda, yield
C         double precision :: dJdSig(6), dxJ3dSig(6),dPdSig(6),dQdSig(6)
C         double precision :: dThetadSig(6),dFdSig(6), dZetadSig(6) 
C         double precision :: dFdpImage, chi_image, xMImage, p_imageMax
C         double precision :: xJ2, xJ3, cos3Theta, theta, xM, psi_image
C         double precision :: dFdP, dFdQ, dFdM, dMdTheta, Zeta, dMdpImage
C         double precision :: C_dFdSig(6), dFdSig_C_dFdSig, dpImagedEpsQ
C         double precision :: xM_image_tc, pi,alam
C         integer :: I
C         PARAMETER (pi=3.14159265359D0)

C         call ZERO1(dJdSig,6)
C         call ZERO1(dxJ3dSig,6)
C         call ZERO1(dPdSig,6)
C         call ZERO1(dQdSig,6)

C         call ZERO1(dThetadSig,6)
C         call ZERO1(dZetadSig,6)

C         xJ2 = 0.
C         xJ3 = 0.
C         cos3Theta = 0.
C         theta  = 0.
C         xM = 0.

C         psi_image = 0.
C         chi_image = 0.
C         xMImage = 0.

C         dFdP = 0.
C         dFdQ = 1.
C         dFdM = 0
C         dMdTheta = 0
C         dFdpImage = 0.
C         dMdpImage = 0
C         Zeta = 0
           
C         call getInvariants (stress, xJ2, dJdSig, xJ3, dxJ3dSig)
C         call getMlode (stress,xMtc,theta,xJ3,xJ2,cos3Theta,xM)
C         call getMaximuxMImagePressure (void, xN, xM, xMtc, chi_tc, xGamma,
C      &            xlambda, p, pImage, p_imageMax, xMImage, xM_image_tc, 
C      &            psi_image, chi_image)

C !---------calculate dFdSigma --------------

C         dFdP=xMImage - q/p;
C         do I = 1, 3
C                dPdSig(I) = 1./3.   
C                if (q == 0 ) then
C                  dQdSig(I)=0
C                else
C                  dQdSig(I) =  3. / 2. / q * (-stress(I) - p)
C                end if
C         end do
C         do I = 4, 6
C                dPdSig(I) = 0. 
C                if (q == 0 ) then
C                  dQdSig(I)=0
C                else
C                  dQdSig(I) =  3. / q * (-stress(I) )
C                endif
C         end do

C         dFdM  = -q/xMImage
C         dMdTheta = (3./2.* xMtc**2/( 3 + xMtc )) * sin( 3*theta/2 + pi/4) 
C         if ( -0.01 < sqrt(xJ2) < 0.01) then ! avoid dividing by 0
C             do I = 1, 6
C                 dThetadSig(I) = 0.
C             end do
C         else
C             if (-1.e-16 < cos3Theta < 1e-16) then  ! avoid dividing by 0
C               if (cos3Theta == 0) then 
C                   do I = 1, 6
C                       dThetadSig(I) = 0.
C                   end do
C               else
C                   do I = 1, 6
C                   ! Need to confirm the correct sign + or - ? 
C                     dThetadSig(I) = sqrt(3.)/2./ cos3Theta/sqrt(xJ2**3.)
C      &                            * (dxJ3dSig(I) - 3./2.*xJ3/xJ2*dJdSig(I))
C                   end do
C               end if
C             end if
C         end if

C         Zeta = 1 - chi_image*xN*abs(psi_image)/xMtc

C         do I = 1, 6
C           if ( psi_image > 0 ) then 
C               dZetadSig(I) =  - 2*chi_image*xN/xMtc*xlambda/p*dPdSig(I)
C           else
C               dZetadSig(I) =  + 2*chi_image*xN/xMtc*xlambda/p*dPdSig(I)
C           end if
C         end do

C         do I = 1, 6
C           dFdSig(I) = dFdP*dPdSig(I) + dFdQ*dQdSig(I) 
C      &          + dFdM*(dMdTheta*dThetadSig(I)*Zeta + xM*dZetadSig(I))  
C         end do

C !---------calculate dFdpImage--------------

C         dFdpImage = -xM*p/pImage 
C         if ( psi_image > 0 ) then
C           dMdpImage =   xM*chi_image*xN/xMtc*xlambda/pImage
C         else
C           dMdpImage = - xM*chi_image*xN/xMtc*xlambda/pImage
C         end if
C         dFdpImage = dFdpImage + dFdM*dMdpImage

C !---------calculate dPidEpsQ for plasticity (hardening/softening) --------------

C         dpImagedEpsQ=xH*xMImage/xM_image_tc*(p_imageMax-pImage)*(p/pImage)

C         call MatVec (DE, 6, dFdSig, 6, C_dFdSig)

C         dFdSig_C_dFdSig = 0
C         do I = 1,6
C           dFdSig_C_dFdSig = dFdSig_C_dFdSig + dFdSig(I)*C_dFdSig(I)
C         end do

C         alam = yield / ( dFdSig_C_dFdSig - dFdpImage*dpImagedEpsQ );

C         do I = 1,6
C           stress(I) = stress(I) - alam*C_dFdSig(I)
C         end do

C         pImage  = pImage + alam*dpImagedEpsQ

C       end subroutine StressCorrection



      subroutine getMlode (stress, xMtc, theta, xJ3, xJ2, cos3Theta, xM)
          include 'vaba_param.inc'
          dimension:: stress(6), dJdSig(6), dxJ3dSig(6)
C           implicit none
C           double precision, intent(in) :: xMtc
C           double precision, intent(inout):: stress(6)
C           double precision, intent(out) :: theta, xJ3, xJ2, cos3Theta, xM
C           double precision :: dJdSig(6), dxJ3dSig(6)
C           double precision ::  xJ3AxJ3,  sin3Theta, pi
          PARAMETER (pi=3.14159265359D0)
          
          call ZERO1(dJdSig,6)
          call ZERO1(dxJ3dSig,6)
          xM = 0.
          theta = 0.
          xJ3 = 0.
          xJ2 = 0.
          cos3Theta = 0.
          xJ3AxJ3 = 0.
          sin3Theta = 0.
 
          call getInvariants(stress,xJ2,dJdSig,xJ3,dxJ3dSig)
          if (xJ2 == 0.) then 
 	            xJ3AxJ3 = 0.
          else
              xJ3AxJ3 = xJ3/sqrt(xJ2**3)
          end if

          sin3Theta = 3.*sqrt(3.)/2. * xJ3AxJ3

          if (sin3Theta > 0.99) sin3Theta = 1.
          if (sin3Theta < -0.99) sin3Theta = -1.
 
          theta = 1./3. * asin(sin3Theta)
          if (theta > 0.523598) theta = 0.523598
          if (theta < -0.523598) theta = -0.523598

          cos3Theta = cos(3.*theta)
          if (-1.E-8 < cos3Theta < 1.E-8) cos3Theta = 0
          
          ! Jefferies \& Shuttle  2011
          xM = xMtc - xMtc**2./(3. + xMtc) * cos(3.*theta /2. + pi / 4.)
  
      end subroutine getMlode
 	  
      subroutine elasticStiffnessMatrixNorSand ( xG, xK, DElast)
          include 'vaba_param.inc'
          dimension:: DElast(6,6)
C         implicit none
C  	      double precision, intent(in) :: xG, xK
C  	      double precision, intent(out) :: DElast(6,6)
          integer :: I, J
          call ZERO2(DElast,6)
          do I = 1, 3
               do J = 1, 3
                   if ( I == J) then 
                       DElast(I,J) = xK + 4./3.*xG     ! diagonal of first block
                   else
                       DElast(I,J) = xK - 2./3.*xG   ! other terms of first block
                   end if
               end do
           end do
           do I = 4, 6
               DElast(I,I) = xG ! diagonal of fourth block
           end do
      end subroutine elasticStiffnessMatrixNorSand

      subroutine getdFdq_ (dFdSig, dFdq_)
        include 'vaba_param.inc'
        dimension:: dFdSig(6)
C         implicit none
C         double precision, intent(in) :: dFdSig(6) 
C         double precision, intent(out) :: dFdq_
C         double precision :: dFdSig_ii
        dFdSig_ii = ( dFdSig(1) + dFdSig(2) + dFdSig(3) )     
        dFdq_ = sqrt(2./3.) * sqrt((dFdSig(1) - dFdSig_ii/3.)**2. 
     &   + (dFdSig(2) - dFdSig_ii/3)**2. + (dFdSig(3)- dFdSig_ii/3)**2.
     &   + 2.*(dFdSig(4))**2.+ 2.*(dFdSig(5))**2.+2.*(dFdSig(6))**2.) 
       end subroutine getdFdq_
 	  
      subroutine getPandQ (stress,p,q,eta)
        include 'vaba_param.inc'
        dimension:: stress(6),sigma(6)
C         implicit none
C         double precision, intent(in) :: stress(6)
C         double precision, intent(out) :: p, q, eta
C         double precision :: sigma(6)
C         integer :: I
        sigma = -stress
        q = 0.
c        do while (q < 0.1)
          if (sigma(1) < 0.1 ) sigma(1) = 0.1     ! tension cut-off in x-direction
          if (sigma(2) < 0.1)  sigma(2) = 0.1     ! tension cut-off in y-direction
          if (sigma(3) < 0.1)  sigma(3) = 0.1     ! tension cut-off in z-direction
          p  =(sigma(1) + sigma(2) + sigma(3))/3.    ! mean effective stress 
          q  =sqrt(0.5*  ( (sigma(1)-sigma(2))**2. 
     &        +(sigma(2)-sigma(3))**2. + (sigma(3)-sigma(1))**2. 
     &        + 6*( sigma(4)**2. + sigma(5)**2. + sigma(6)**2.)))  ! deviatoric stress
c          if (q < 0.1) sigma(1) = sigma(1) + 0.03
c        end do
        eta = q/p
c        stress = -sigma
      end subroutine getPandQ

      subroutine getDevVolStrain (Eps, EpsV, EpsD)
        include 'vaba_param.inc'
        dimension:: Eps(6)
C         implicit none
C         double precision, intent(in) :: Eps(6)
C         double precision, intent(out) :: EpsV, EpsD
        EpsV = Eps(1) + Eps(2) + Eps(3)  ! volumetric strain
        EpsD = sqrt(2./3.)*sqrt( (Eps(1)-EpsV/3.)**2.
     &          + (EpS(2)-EpsV/3.)**2. + (EpS(3)-EpsV/3.)**2. 
     &          + (Eps(4)**2.)/2. + (Eps(5)**2.)/2. + (Eps(6)**2.)/2.) ! deviatoric strain  
      end subroutine getDevVolStrain
 	  
      subroutine getInvariants (stress, xJ2, dJdSig,xJ3, dxJ3dSig)
        include 'vaba_param.inc'
        dimension:: stress(6),dJdSig(6),dxJ3dSig(6),S(6)
C         implicit none
C         double precision, intent(inout) :: stress(6)
C         double precision, intent(out) ::xJ2,dJdSig(6),xJ3,dxJ3dSig(6)
C         double precision :: p, q, S(6), eta
C         integer :: I
         
        call getPandQ (stress,p,q,eta)
 
        S(1) = -stress(1) - p
        S(2) = -stress(2) - p
        S(3) = -stress(3) - p
        S(4) = -stress(4)
        S(5) = -stress(5)
        S(6) = -stress(6)                                          
 	  
        xJ3 = S(1)*S(2)*S(3)-S(1)*S(6)**2-S(2)*S(5)**2
     &       -S(3)*S(4)**2+2*S(4)*S(6)*S(5)
 	 
        xJ2 = 1./6.*((stress(1)-stress(2))**2
     &       +(stress(2)-stress(3))**2 + (stress(3)-stress(1))**2 )
     &       + stress(4)**2 + stress(5)**2 + stress(6)**2
           !if (abs(xJ2) < 0.0001) xJ2  =  0.0
 	 
        dJdSig(1) = S(1)
        dJdSig(2) = S(2)
        dJdSig(3) = S(3)
        dJdSig(4) = -2.*stress(4)
        dJdSig(5) = -2.*stress(5)
        dJdSig(6) = -2.*stress(6)

        dxJ3dSig(1) = -1./3.*S(1)*S(2) -1./3.*S(1)*S(3) +2./3.*S(2)*S(3)
     &               -2./3.*S(6)**2 +1./3.*S(5)**2 +1./3.*S(4)**2
        dxJ3dSig(2) = -1./3.*S(1)*S(2) +2./3.*S(1)*S(3) -1./3.*S(2)*S(3) 
     &               +1./3.*S(6)**2 -2./3.*S(5)**2 +1./3.*S(4)**2
        dxJ3dSig(3) = 2./3.*S(1)*S(2) -1./3.*S(1)*S(3) -1./3.*S(2)*S(3)
     &               +1./3.*S(6)**2 +1./3.*S(5)**2 -2./3.*S(4)**2
        dxJ3dSig(4) = -2.*S(3)*S(4)+2.*S(6)*S(5)
        dxJ3dSig(5) = -2.*S(2)*S(5)+2.*S(4)*S(6)
        dxJ3dSig(6) = -2.*S(1)*S(6)+2.*S(4)*S(5)

       end subroutine getInvariants

       subroutine getMaximuxMImagePressure (void, xN, xM, xM_tc, chi_tc, 
     &            xGamma, xlambda, p, p_image, p_imageMax, xM_image,
     &            xM_image_tc, psi_image, chi_image)
         include 'vaba_param.inc'
C          implicit none
C          double precision, intent(in) :: void,xN,xM,xM_tc,chi_tc,xGamma
C          double precision, intent(in) :: xlambda, p, p_image
C          double precision, intent(out) :: p_imageMax, xM_image,xM_image_tc
C          double precision, intent(out) :: psi_image, chi_image
C          double precision ::  eCrit, psi
               
! get Maximum Image Pressure for internal cap and dilation limitation

         chi_image = chi_tc / ( 1-xlambda*chi_tc/xM_tc);
         if (p > 1.) then
              eCrit = xGamma-xlambda*log(p/1000)      ! critical  void ratio
         else
              eCrit = xGamma                    ! critical  void ratio
         end if
         !if (eCrit > eMax) eCrit = eMax        ! maximum critical void ratio is the maximum void ratio
         !if (eCrit < eMin) eCrit = eMin        ! minimum critical void ratio is the minimum void ratio
         psi = void - eCrit                     !state parameter

         if (p > 1.) then
               psi_image = psi +  xlambda*log(p_image/p)     ! image state parameter
         else
               psi_image = psi +  xlambda*log(p_image)       ! image state parameter
         end if

         xM_image = xM*(1 - chi_image*xN*abs(psi_image)/xM_tc)  !image critical stress ratio 

         xM_image_tc = xM_tc*(1 - chi_image*xN*abs(psi_image)/xM_tc)  !image critical stress ratio at triaxial compression 

         p_imageMax=p*exp(- chi_image*psi_image/xM_image_tc)
          
      end subroutine getMaximuxMImagePressure

C       subroutine updateM_Image_unloading(void, xN, xM, xMtc,chi_tc, 
C      &                xGamma,xlambda, p, q, p_image, xM_image)

C            implicit none
C            double precision, intent(in) :: void, xN, xM, xMtc, chi_tc
C            double precision, intent(in) :: xGamma, xlambda, p, q
C            double precision, intent(inout):: p_image, xM_image
C            double precision :: eCrit, psi
C            double precision :: chi_image, psi_image

               
C ! get Maximum Image Pressure for internal cap and dilation limitation
           
C c          eta = q/p  
C c          p_image = p*exp( eta/xM_image-1)

C           chi_image = chi_tc / ( 1-xlambda*chi_tc/xMtc)

C           if (p > 1.) then
C               eCrit = xGamma-xlambda*log(p/1000)     !unit in kPa
C           else
C               eCrit = xGamma         ! critical  void ratio
C           end if
C           psi = void - eCrit         ! state parameter

C           if (p > 1.) then
C                psi_image = psi +  xlambda*log(p_image/p)     ! image state parameter
C           else
C                psi_image = psi +  xlambda*log(p_image)       ! image state parameter
C           end if
C           xM_image = xM*(1 - chi_image*xN*abs(psi_image)/xMtc)  !image critical stress ratio
          
C       end subroutine updateM_Image_unloading


      subroutine updateM_Image(void, xN, xM, xMtc,chi_tc, xGamma,
     &                xlambda, p, q, p_image, xM_image)
           include 'vaba_param.inc'
C            implicit none
C            double precision, intent(in) :: void, xN, xM, xMtc, chi_tc
C            double precision, intent(in) :: xGamma, xlambda, p, q
C            double precision, intent(inout):: p_image, xM_image
C            double precision :: eCrit, psi
C            double precision :: chi_image, psi_image
               
! get Maximum Image Pressure for internal cap and dilation limitation
           
c          eta = q/p  
c          p_image = p*exp( eta/xM_image-1)

          chi_image = chi_tc / ( 1-xlambda*chi_tc/xMtc)

          if (p > 1.) then
              eCrit = xGamma-xlambda*log(p/1000)     !unit in kPa
          else
              eCrit = xGamma         ! critical  void ratio
          end if
          psi = void - eCrit         ! state parameter

          if (p > 1.) then
               psi_image = psi +  xlambda*log(p_image/p)     ! image state parameter
          else
               psi_image = psi +  xlambda*log(p_image)       ! image state parameter
          end if
          xM_image = xM*(1 - chi_image*xN*abs(psi_image)/xMtc)  !image critical stress ratio
          
      end subroutine updateM_Image
      
      
       subroutine getP_M_Image(void, xN, xM, xMtc,chi_tc, xGamma,xlambda,
     &                 p, q, p_image, psi_image, chi_image, xM_image)
           include 'vaba_param.inc'
C            implicit none
C            double precision, intent(in) :: void, xN, xM, xMtc, chi_tc
C            double precision, intent(in) :: xGamma, xlambda, p, q
C            double precision, intent(out):: p_image, psi_image
C            double precision, intent(out):: chi_image, xM_image
C            double precision :: eCrit, psi, temp1, temp2, a, b, c, delta
C            double precision :: eta, x1, x2, pImage1, pImage2,F
C            double precision :: psi_image1, psi_image2
C            integer :: iflag, isymbol,loop
               
! get Maximum Image Pressure for internal cap and dilation limitation

           chi_image = chi_tc / ( 1-xlambda*chi_tc/xMtc);

           if (p > 1.) then
              eCrit = xGamma-xlambda*log(p/1000)        ! critical  void ratio in unit of kPa
           else
              eCrit = xGamma         ! critical  void ratio
           end if
           !if (eCrit > eMax) eCrit = eMax        ! maximum critical void ratio is the maximum void ratio
           !if (eCrit < eMin) eCrit = eMin        ! minimum critical void ratio is the minimum void ratio
           psi = void - eCrit                    !state parameter

           eta = q/p

           iflag =1
           isymbol =1
           loop = 0
           
           do while (iflag == 1 .and. loop<=2)
             loop = loop + 1
             temp1 = isymbol*psi
             temp2 = isymbol*xlambda
             a = - chi_image*xN*temp2/xMtc
             b = (1-chi_image*xN*temp1/xMtc)- chi_image*xN*temp2/xMtc
             c = 1-chi_image*xN*temp1/xMtc - eta/xM
             delta = b**2-4*a*c
             if (delta > 0) then
               x1 = ( -b+sqrt(delta) ) / (2*a)
               x2 = ( -b-sqrt(delta) ) / (2*a)    
               pImage1 = p*exp(x1)
               pImage2 = p*exp(x2)
               psi_image1 = psi + xlambda*log(pImage1/p)
               psi_image2 = psi + xlambda*log(pImage2/p)
               if  (psi_image1*isymbol > 0 ) then
                 p_image = pImage1
                 psi_image = psi_image1
                 xM_image = xM*(1-chi_image*xN/xMtc*abs(psi_image))
                 F = 1-log(p/p_image)-eta/xM_image
                 if (abs(F)<1e-6 
     &              .and. xM_image>0.1*xM .and. xM_image<xM) then
                   iflag = 0
                   return
                 end if
               end if
               if (psi_image2*isymbol > 0) then
                 p_image = pImage2
                 psi_image = psi_image2
                 xM_image = xM*(1-chi_image*xN/xMtc*abs(psi_image))
                 F = 1-log(p/p_image)-eta/xM_image
                 if (abs(F)<1e-6 
     &              .and. xM_image>0.1*xM .and. xM_image<xM) then
                   iflag = 0
                   return
                 end if
               end if
               isymbol = -1*isymbol
             else if (abs(delta) < 1e-6 ) then
               p_image = -b/(2*a)
               psi_image = psi + xlambda*log(p_image/p)
               xM_image = xM*(1-chi_image*xN/xMtc*abs(psi_image))
               F = 1-log(p/p_image)-eta/xM_image
               if ( abs(F)<1e-6 
     &               .and. xM_image>0.1*xM .and. xM_image<xM ) then
                 iflag = 0
                 return
               end if
             else
               isymbol = -1*isymbol
             end if
           end do
          
       end subroutine getP_M_Image

  
       subroutine getdFdSig(void, stress, p, q, pImage, xN, xMtc, chi_tc,
     &            xGamma, xlambda, dFdSig )
        include 'vaba_param.inc'
        dimension:: stress(6), dFdSig(6),dJdSig(6), dxJ3dSig(6),
     &     dPdSig(6),dQdSig(6),dZetadSig(6),dThetadSig(6)
C         implicit none
C         double precision, intent(inout) :: void, stress(6), pImage
C         double precision, intent(in) :: p, q, xN, xMtc, chi_tc
C         double precision, intent(in) :: xGamma, xlambda
C         double precision, intent(out) :: dFdSig(6)
C         double precision :: dJdSig(6), dxJ3dSig(6),dPdSig(6),dQdSig(6)
C         double precision :: dZetadSig(6),dThetadSig(6)
C         double precision :: xJ2, xJ3, cos3Theta, theta, xM, psi_image
C         double precision :: chi_image, xMImage
C         double precision :: dFdP, dFdQ, dFdM, dMdTheta, Zeta,pi
C         integer :: I
        PARAMETER (pi=3.14159265359D0)

        call ZERO1(dJdSig,6)
        call ZERO1(dxJ3dSig,6)

        call ZERO1(dPdSig,6)
        call ZERO1(dQdSig,6)

        call ZERO1(dThetadSig,6)
        call ZERO1(dZetadSig,6)

        xJ2 = 0.
        xJ3 = 0.
        cos3Theta = 0.
        theta  = 0.
        xM = 0.

        psi_image = 0.
        chi_image = 0.
        xMImage = 0.
        dFdP = 0.
        dFdQ = 1.

        dFdM = 0
        dMdTheta = 0

        Zeta = 0
           
        call getInvariants (stress, xJ2, dJdSig, xJ3, dxJ3dSig)
        call getMlode (stress,xMtc,theta,xJ3,xJ2,cos3Theta,xM)
        call getP_M_Image (void, xN, xM, xMtc, chi_tc, xGamma, xlambda,
     &              p, q, pImage, psi_image, chi_image, xMImage)

        dFdP=xMImage*log(p/pImage)
        do I = 1, 3
               dPdSig(I) = 1./3.   
               if (q == 0 ) then
                 dQdSig(I)=0
               else
                 dQdSig(I) =  3. / 2. / q * (-stress(I) - p)
               end if
        end do
        do I = 4, 6
               dPdSig(I) = 0. 
               if (q == 0 ) then
                 dQdSig(I)=0
               else
                 dQdSig(I) =  3. / q * (-stress(I) )
               endif
        end do

        dFdM = -q/xMImage
        Zeta = 1 - chi_image*xN*abs(psi_image)/xMtc
        dMdTheta = Zeta*(3./2. * xMtc**2 /(3+ xMtc))*sin(3*theta/2+pi/4) 
        if ( -0.01 < sqrt(xJ2) < 0.01) then ! avoid dividing by 0
               do I = 1, 6
                   dThetadSig(I) = 0.
               end do
        else
               if (-1.e-16 < cos3Theta < 1e-16) then  ! avoid dividing by 0
                   if (cos3Theta == 0) then 
                       do I = 1, 6
                           dThetadSig(I) = 0.
                       end do
                   else
                       do I = 1, 6
                       ! Need to confirm the correct sign + or - ? 
         dThetadSig(I) = sqrt(3.)/2./ cos3Theta/sqrt(xJ2**3.)*(dxJ3dSig(I)
     &                   - 3./2.*xJ3/xJ2*dJdSig(I) ) 
                       end do
                   end if
               end if
        end if

        

        do I = 1, 6
               dFdSig(I) = dFdP*dPdSig(I) + dFdQ*dQdSig(I) 
     &          + dFdM*dMdTheta*dThetadSig(I) 
        end do

       end subroutine  getdFdSig


       subroutine getdFdPImage(void, stress, p, q, pImage, xN, xMtc, 
     &            chi_tc, xGamma, xlambda, dFdpImage)
         include 'vaba_param.inc'
         dimension:: stress(6),dxJ3dSig(6),dPdSig(6),dJdSig(6),dQdSig(6)
C          implicit none
C          double precision, intent(inout) :: void, stress(6),pImage
C          double precision, intent(in) :: p, q, xN, xMtc, chi_tc
C          double precision, intent(in) :: xGamma, xlambda
C          double precision, intent(out) :: dFdpImage
C          double precision :: xJ2, xJ3, cos3Theta, theta, xM, psi_image
C          double precision :: dFdM, dMdpImage, chi_image, xMImage
C          double precision :: dxJ3dSig(6),dPdSig(6),dJdSig(6),dQdSig(6)
C          integer :: I

           call ZERO1(dJdSig,6)
           call ZERO1(dxJ3dSig,6)

           call ZERO1(dPdSig,6)
           call ZERO1(dQdSig,6)

           xJ2 = 0.
           xJ3 = 0.
           cos3Theta = 0.
           theta  = 0.
           xM = 0.

           psi_image = 0.
           chi_image = 0.
           xMImage = 0.

           dFdpImage = 0.
           dFdM = 0
           dMdpImage = 0
   
           call getMlode (stress,xMtc,theta,xJ3,xJ2,cos3Theta,xM)

           call getP_M_Image (void, xN, xM, xMtc, chi_tc, xGamma, xlambda,
     &              p, q, pImage, psi_image, chi_image, xMImage)

           dFdpImage = -xM*p/pImage 

           dFdM  = -q/xMImage

           if ( psi_image > 0 ) then
               dMdpImage =   xM*chi_image*xN/xMtc*xlambda/pImage
           else
               dMdpImage = - xM*chi_image*xN/xMtc*xlambda/pImage
           end if

           dFdpImage = dFdpImage + dFdM*dMdpImage

       end subroutine  getdFdPImage
    
       subroutine getYieldFunctionNorSand (p,q,p_image,xM_image,yield)
           include 'vaba_param.inc'
C            implicit none
C            double precision, intent(in) :: p, q, xM_image, p_image
C            double precision, intent(out) :: yield

           yield  =  q/p - xM_image + xM_image*log(p/p_image)
             
       end subroutine getYieldFunctionNorSand


 
C       subroutine getElasticPartPegasus(yield0, yield1, alpha0, alpha1, 
C      &   stress, dStrTrial, pImage, xMImage, DE, subdStrain, void)
C         implicit none
C         double precision,intent(in):: dStrTrial(6)
C         double precision,intent(in):: pImage, xMImage, DE(6,6)
C         double precision,intent(inout)::stress(6), void, alpha0, alpha1
C         double precision,intent(inout)::yield0, yield1, subdStrain(6)
C         double precision :: yieldNew,SigNew(6),pNew,qNew,etaNew,VOL,DEV
C         double precision :: dElastStrain(6), dElastStress(6), alpha
C         integer :: check, iteration, I, Tol

C         Tol = 1e-5
C         check = 0
C         !get elstic contribution (pegasus method)
C         yieldNew= 1000
C         iteration = 0
C         do while ( abs(yieldNew) > Tol ) 
C           alpha = alpha1 - yield1*(alpha1-alpha0)/(yield1-yield0)
C           do I = 1,6
C             SigNew(I) = stress(I) + alpha*dStrTrial(I)
C           end do
C           call getPandQ (SigNew, pNew, qNew, etaNew)
C           call getYieldFunctionNorSand(pNew,qNew,pImage,xMImage,yieldNew)
C           if ( (yieldNew*yield1) < 0. ) then
C             alpha0 = alpha1
C             yield0 = yield1
C           else
C             yield0 = yield1*yield0/(yield1+yieldNew)
C           end if    
C           alpha1 = alpha
C           yield1 = yieldNew

C           iteration = iteration + 1
C           if (iteration > 10) then !control maximum iteration
C             yieldNew = 0.
C             alpha = 0.
C             check = 999 
C           else
C             check = 1
C           end if
C         end do 
              
C         if (alpha /= alpha) then
C             alpha = 0.
C             check = -1
C         elseif (alpha > 1.) then
C               alpha = 0.
C               check = -1
C         elseif (alpha < 0.) then 
C               alpha = 0.
C               check = -1
C         end if
C         if (check == 1.) then ! update stress for elastic part
C           do I = 1, 6
C             dElastStrain (I) = alpha * subdStrain(I)
C             subdStrain (I) = (1-alpha) * subdStrain(I)
C           end do
C           call MatVec (DE,6, dElastStrain, 6, dElastStress)
C           do I = 1, 6
C             stress(I) = stress(I) + dElastStress(I)
C           end do
C           call getDevVolStrain(dElastStrain,VOL,DEV)
C           void = void + VOL * ( 1. + void)
C         end if
C       end subroutine getElasticPartPegasus

C       subroutine getElasticPartUnload(yield0, stress,dStrTrial,
C      &           pImage, xMImage, DE, subdStrain, void, FTOL)
C         implicit none
C         double precision,intent(in)::dStrTrial(6), FTOL
C         double precision,intent(in)::pImage,xMImage,DE(6,6)
C         double precision,intent(inout)::stress(6), void, subdStrain(6)
C         double precision,intent(inout):: yield0
C         double precision :: yieldNew,SigNew(6),pNew,qNew,etaNew,VOL,DEV
C         double precision :: dElastStrain(6),dElastStress(6)
C         double precision :: alpha0, alpha1, yield1, yieldsave
C         double precision :: dalpha, alpha
C         integer :: check, iteration, I, id, ic, Tol, Nsub, MaxIts
 
C         Tol = 1e-5
C         check = 0
C         !get elstic contribution (pegasus method)
C         yieldNew= 1000
C         iteration = 0
C         Nsub = 10
C         MaxIts = 3

C         alpha0 = 0.
C         alpha1 = 1.
C         yieldsave = yield0
C         do id =1, MaxIts
C           dalpha = (alpha1 - alpha0)/Nsub
C           do ic = 1, Nsub    
C             alpha = alpha0 + dalpha 
C             do I = 1,6
C               SigNew(I) = stress(I) + alpha*dStrTrial(I)
C             end do
C             call getPandQ (SigNew, pNew, qNew, etaNew)
C             call getYieldFunctionNorSand(pNew,qNew,pImage,
C      &               xMImage,yieldNew)
C             if ( yieldNew > FTOL ) then
C               alpha1 = alpha
C               if ( yield0 < -FTOL ) then
C                 yield1 = yieldNew
C                 go to 200
C               else
C                 alpha0 = 0
C                 yield0 = yieldsave
C                 EXIT
C               end if
C             else
C               alpha0 =  alpha
C               yield0 = yieldNew
C             end if
C           end do
C         end do
C         !did not found the intersection with yield surface
C         stop 'Error- cannot find the intersection during unloading' 
C         !call XPLB_EXIT
C         !found the intersection, call Pegasus scheme to determine the exact value
C 200     call getElasticPartPegasus(yield0, yield1, alpha0, alpha1,
C      &    stress,dStrTrial,pImage, xMImage, DE, subdStrain, void)

C       end subroutine getElasticPartUnload


      subroutine errorEuler (SigInitial, dSig1, dSig2, pImage, dPi1, 
     &       dPi2, Sigma, timeInc,iflag_sub,STOL)
        include 'vaba_param.inc'
        dimension:: SigInitial(6),dSig1(6),dSig2(6), Sigma(6),E(6)
C         implicit none
C         double precision, intent(in) :: SigInitial(6),dSig1(6),dSig2(6)
C         double precision, intent(in) :: pImage, dPi1, dPi2
C         double precision, intent(inout) :: timeInc
C         double precision, intent(out) :: Sigma(6)
C         integer, intent(inout) :: iflag_sub
C         double precision :: error,reductionFactor,E(6)!,error1, error2 
C         double precision :: EnorSqrt, SigNorSqrt, STOL!, SigmaNorm
C         integer :: I
c        STOL = 0.001
        EnorSqrt = 0.
        SigNorSqrt = 0.
        do I = 1,6
          E(I) = 0.5*(-dSig1(I) + dSig2(I))
          EnorSqrt = EnorSqrt + E(I)**2.
          Sigma(I) = SigInitial(I) + 0.5*(dSig1(I) + dSig2(I))
          SigNorSqrt = SigNorSqrt + Sigma(I)**2.
        end do
c        error1 = abs(dPi1-dPi2)/2/(pImage+0.5*(dPi1+dPi2))
        error = sqrt(EnorSqrt / SigNorSqrt)/2
C         error  = max(error1, error2)
        if (error > STOL) then
          reductionFactor = 0.9*sqrt(STOL/error)
          if (reductionFactor < 0.1) reductionFactor = 0.1
          if (reductionFactor > 2.) reductionFactor = 2.
          timeInc = reductionFactor * timeInc
          iflag_sub = 1
        else
          reductionFactor = min(0.9*sqrt(STOL/error),1.1)
          if (reductionFactor < 0.1) reductionFactor = 0.1
          if (reductionFactor > 2.)  reductionFactor = 2.
          timeInc = reductionFactor * timeInc
          iflag_sub = 0
        end if  
      end subroutine errorEuler

 ! -----------------------------------------------------------------------
 !                     INTERNAL MATH SUBROUTINES
 ! -----------------------------------------------------------------------
      subroutine ZERO2(A,N)
        include 'vaba_param.inc'
        dimension:: A(N,N)
C         implicit none
C         integer, intent(in) :: N
C         double precision, intent(inout) :: A (N,N)
C         integer :: I,J
 	     do I = 1,N
 	         do J=1,N
 	             A(I,J) = 0.
 	         end do
 	     end do
      end subroutine ZERO2
 
      subroutine ZERO1(A,N)
        include 'vaba_param.inc'
        dimension:: A(N)
C           implicit none
C           integer, intent(in) :: N
C           double precision, intent(inout) :: A(N)
C           integer :: I
 	      do I = 1, N
 	        A(I) = 0.
 	      end do
      end subroutine ZERO1
 	          
      subroutine addMatComp (A, B, C)
        include 'vaba_param.inc'
        dimension:: A(6,6), B(6,6), C (6,6)
C            implicit none
C            double precision, intent(in) :: A(6,6), B(6,6)
C            double precision, intent(out) :: C (6,6)
C            integer :: I , J
 	      do I = 1, 6
              do J = 1, 6
                      C(I,J) = A(I,J) + B(I,J)
              end do
           end do
      end subroutine addMatComp

      subroutine AddVec (A, B, coef1, coef2, N, C)
        include 'vaba_param.inc'
        dimension:: A(N), B(N), C(N)
C           implicit none
C           double precision, intent(in) :: A(N), B(N)
C           double precision, intent(in) :: coef1, coef2
C           double precision, intent(out) :: C(N)
C           integer :: I, J, N

          call ZERO1(C,N)
          do I = 1, N
             C(I) = coef1*A(I) + coef2*B(I)
          end do
      end subroutine AddVec
      
      subroutine MatVec (A, M, B, N, C)
        include 'vaba_param.inc'
        dimension:: A(M,N), B(N), C(N)
C           implicit none
C           double precision, intent(in) :: A(M,N), B(N)
C           double precision, intent(out) :: C(N)
C           integer :: I, J, M, N
          call ZERO1(C,N)
          do I = 1, M
              do J = 1, N
              	C(I) = C(I) + A(I,J) * B(J)
              end do
          end do
      end subroutine MatVec

C       subroutine getElastoPlasticStiffnessUnloading(void,stress,
C      &            subdStrain, p, q, pImage, xN, xH, xMtc, chi_tc,
C      &            xGamma,xlambda, DE, DDSDDE, dpImage,dStrainE)
C         implicit none
C         double precision, intent(inout) :: stress(6)
C         double precision, intent(in) :: void, subdStrain(6)
C         double precision, intent(in) :: p,q, pImage, xN, xH, xMtc,chi_tc
C         double precision, intent(in) :: xGamma, xlambda, DE(6,6)
C         double precision, intent(out):: DDSDDE(6,6), dpImage,dStrainE(6)
C         double precision :: dJdSig(6), dxJ3dSig(6),dPdSig(6),dQdSig(6)
C         double precision :: dThetadSig(6),dFdSig(6)!, dZetadSig(6),dFdQ,dFdM,
C         double precision :: dFdpImage, chi_image, xMImage, p_imageMax
C         double precision :: xJ2, xJ3, cos3Theta, theta, xM, psi_image
C         double precision :: dFdP, dMdTheta, Zeta, dMdpImage
C         double precision :: C_dFdSig(6), dpImagedEpsQ
C         double precision :: xM_image_tc, pi, denominater, phi, dStress(6)
C         double precision :: DE_dFdSig(6), DE_dFdSig_dFdSig(6,6), DP(6,6)
C         double precision :: dFdq_, dPOTdP, dPOTdQ, dPOTdM, q_yield
C         double precision :: dFdSig_C_dPdSig,dPOTdSig(6)
C         double precision :: xM_tc_image

C         integer :: I, J, Z, sign_
C         PARAMETER (pi=3.14159265359D0)

C         call ZERO1(dJdSig,6)
C         call ZERO1(dxJ3dSig,6)
C         call ZERO1(dPdSig,6)
C         call ZERO1(dQdSig,6)

C         call ZERO1(dThetadSig,6)
C C         call ZERO1(dZetadSig,6)

C         xJ2 = 0.
C         xJ3 = 0.
C         cos3Theta = 0.
C         theta  = 0.
C         xM = 0.

C         psi_image = 0.
C         chi_image = 0.
C         xMImage = 0.

C         dFdP = 0.
        
C C         dFdM = 0.
C         dMdTheta = 0.
C         dFdpImage = 0.
C         dMdpImage = 0.
C         Zeta = 0.

C         call getInvariants (stress, xJ2, dJdSig, xJ3, dxJ3dSig)
C         call getMlode (stress,xMtc,theta,xJ3,xJ2,cos3Theta,xM)
C         call getMaximuxMImagePressure (void, xN, xM, xMtc, chi_tc, xGamma,
C      &            xlambda, p, pImage, p_imageMax, xMImage, xM_image_tc,
C      &            psi_image, chi_image)

C !---------calculate dPdSigma --------------

C         !First P represents the potential function

C         !Calculate the deviatoric stress q which make (p,q) stays on yield surfaces at the cap location

C         q_yield = p*xMImage*(1-log(p/pImage))

C         dPOTdP=(xMImage - q_yield/p)
C         dPOTdQ = 1.
C         do I = 1, 3
C                dPdSig(I) = 1./3.
C                if (q == 0 ) then
C                  dQdSig(I)=0
C                else
C                  dQdSig(I) =  3. / 2. / q * (-stress(I) - p)
C                end if
C         end do
C         do I = 4, 6
C                dPdSig(I) = 0.
C                if (q == 0 ) then
C                  dQdSig(I)=0
C                else
C                  dQdSig(I) =  3. / q * (-stress(I) )
C                endif
C         end do

C         dPOTdM  = q_yield/xMImage
C         dMdTheta = (1 - chi_image*xN*abs(psi_image)/xMtc)
C      &            *(3./2.* xMtc**2/( 3 + xMtc )) * sin( 3*theta/2 + pi/4)
C         if ( -0.01 < sqrt(xJ2) < 0.01) then ! avoid dividing by 0
C             do I = 1, 6
C                 dThetadSig(I) = 0.
C             end do
C         else
C             if (-1.e-16 < cos3Theta < 1e-16) then  ! avoid dividing by 0
C               if (cos3Theta == 0) then
C                   do I = 1, 6
C                       dThetadSig(I) = 0.
C                   end do
C               else
C                   do I = 1, 6
C                   ! Need to confirm the correct sign + or - ? 
C                     dThetadSig(I) = sqrt(3.)/2./ cos3Theta/sqrt(xJ2**3.)
C      &                            * (dxJ3dSig(I) - 3./2.*xJ3/xJ2*dJdSig(I))
C                   end do
C               end if
C             end if
C         end if

C         !Use outer yield surface as potential for unloading, the direction should toward inward
C         do I = 1, 6
C           dPOTdSig(I) = (dPOTdP*dPdSig(I) + dPOTdQ*dQdSig(I) 
C      &          + dPOTdM*dMdTheta*dThetadSig(I))
C         end do

C         !WRITE (6,'(''plastic strain increment'')')
C         !write(6,*) dPOTdSig(1), dPOTdSig(2)

C !---------calculate dFdSigma --------------
C         DFdP = - pImage/p**2
C         do I = 1, 6
C           dFdSig(I) = dFdP*dPdSig(I)
C         end do

C !---------calculate dFdpImage--------------

C         chi_image = chi_tc / ( 1-xlambda*chi_tc/xMtc)
C         psi_image = void - xGamma + xlambda*log(pImage/1000)
C         xM_tc_image = xMtc*(1 - chi_image*xN*abs(psi_image)/xMtc)

C         if (psi_image<0) then
C             sign_ = -1
C         else
C             sign_ = 1
C          endif
C         dFdpImage = 1/p-exp(-chi_image/xM_tc_image*psi_image)
C      &              *(-chi_image*xlambda/pImage*xM_tc_image
C      &               -chi_image**2*xN*sign_*xlambda/pImage*
C      &               psi_image)/xM_tc_image**2

C !---------calculate dPidEpsQ for plasticity (hardening/softening) --------------

C         dpImagedEpsQ=  xH*xMImage/xM_tc_image/2*pImage

C         call getdFdq_(dPOTdSig, dFdq_) !Need double check this place

C         do I = 1, 6
C           C_dFdSig(I) = 0.
C           do J = 1, 6
C             C_dFdSig(I) = C_dFdSig(I) + dFdSig(J)*DE(J,I)
C           end do
C         end do
                
C         dFdSig_C_dPdSig = 0
C         do I = 1,6
C           dFdSig_C_dPdSig = dFdSig_C_dPdSig + dPOTdSig(I)*C_dFdSig(I)
C         end do

C         !check this place for positive ad negative plastic deviatoric strain

C         denominater = -dFdSig_C_dPdSig + dFdpImage*dpImagedEpsQ*dFdq_

C         do I = 1, 6
C           DE_dFdSig(I) = 0.
C           do J = 1, 6
C             DE_dFdSig(I) = DE_dFdSig(I)+DE(I,J)*dPOTdSig(J)
C           end do
C         end do

C         do I = 1, 6
C           do J = 1, 6
C             DE_dFdSig_dFdSig(I,J) = DE_dFdSig(I)*dFdSig(J)/denominater
C           end do
C         end do
        
C         do I = 1, 6
C           do J = 1, 6
C             DP(I,J) = 0.
C            do Z = 1, 6
C       DP(I,J) = DP(I,J) + DE_dFdSig_dFdSig(I,Z) * DE(Z,J)
C             end do
C           end do
C         end do

C         ! obtain elastoplasticity stiffness matrix 
C         !write(6,*), dFdSig_C_dPdSig/denominater

C         do I = 1, 6
C           do J = 1, 6
C             !DDSDDE(I,J) = DE(I,J)*dFdSig_C_dPdSig/denominater! - 0.5*(DP(I,J)+DP(J,I))
C             DDSDDE(I,J) = DE(I,J)- 0.5*(DP(I,J)+DP(J,I))
C           end do
C         end do

C         ! update the hardening softening parameter pImage with multipler
C         call MatVec (DE,6, subdStrain, 6, dStress)

C         phi = 0
C         do I = 1,6
C           phi = phi + dFdSig(I)*dStress(I)/denominater
C         end do
        
C         do I = 1,6
C           dStrainE(I) = subdStrain(I) - phi*dPOTdSig(I)
C         end do

C         !write(6,*), phi

C         !Note dFdQ = 1, dEpsQ = phi*dFdQ = phi

C         dpImage = phi*dpImagedEpsQ

C         !should check the positive anbd negative?

C       end subroutine getElastoPlasticStiffnessUnloading
