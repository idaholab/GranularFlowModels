! coded by Wencheng Jin, wencheng.jin@inl.gov
! Sequence of input variable for material constant
!      ei0   = props(1)
!      ed0   = props(2)
!      ec0   = props(3)
!      phic  = props(4)*PI/180
!      hs    = props(5)
!      beta  = props(6)
!      en    = props(7)
!      alpha = props(8)
!      alc    = props(9)    internal lenght parameter
!      1/0    = props(10)  switch of nonlocal enhancement
!      switch_time =  props(11)  control when to switch
!      ile_elem = props(12)  control when to switch
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
! Define the dimension and type of varibles from abaqus, no change needed
      include 'vaba_param.inc'
      integer nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal
!
      dimension props(nprops), density(nblock), coordMp(nblock,3),
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
!  Declaration of local variables
      integer ic, i, j, km, im, nsub, nerror, ny, nBG,
     &  nrelop, naccst, nonlocal_switch
      parameter (tol = 1.d-4, tolabs= 1.d-3)
! Maximum number of elements can be used in the nonlocal simulation
! Need to confirm for each of the simulation
      parameter (ile_elem=17000)
      dimension stress(6), defgrad(9), Tv(6), Dv(6), dstrain(6), D(3,3),
     &  DnlList(nblock), rotDiff(3,3), TvR(6), R(3,3), RT(3,3),
     &  coords(3), stretchU(6), eye(3,3), sigma(3,3), 
     &  coordX(2,ile_elem), coordY(2,ile_elem), coordZ(2,ile_elem),
     &  dloc(2,ile_elem), ElArea(2,ile_elem)
!
      common /KBLOKI/ time, lpocz, lblok, nshift, ncol
      common /KNLOCI/ coordX, coordY, coordZ, ElArea, dloc
!      
      logical prinfo
      character var1, var2
      character*80 cmname
!      
! print informations about time integration, useful when problems occur
      prinfo =.true.!.false.!
! minimal time substep size
      hmin = 1.d-20
! setting nerror = 0, calculation OK      
      nerror = 0
! provent from divided by zero with infinite output  
      if (dt<=0.d0) dt = 1.d0  
! Internal lenght parameter      
      alc   = props(9)
! Prepared for Gauss weight functions
      Gcoeff = 1./(alc*sqrt(2*asin(1.0)))
      Gcutoff =(3*alc)**2
! nonlocal switch  
! props(10)=1 -> nonlocal on
! props(10)=0 -> nonlocal off
      nonlocal_switch = props(10)
      switch_time = props(11)
      i_elem = props(12)
!      WRITE(6,*),'stepTime, totalTime, dt', stepTime, totalTime, dt
! parameter of the numerical differentiation
      theta = 1.0d-7
!*************************************************************************
! nonlocal preprocessing      
!*************************************************************************
      if (nonlocal_switch == 1 .and. totalTime >=(switch_time-dt)) then
        DnlList=0.d0
        if (nshift.eq.0.0.or.ncol.eq.0.0) then
          nshift=1
          ncol=2
        endif
        if (time.ne.stepTime) then !new time increment
          time=stepTime 
          lpocz=1
	        lblok=nblock
	        nshift=1+mod(nshift,2)
	        ncol=1+mod(ncol,2)
	      else !within the time increment
	        lpocz=lpocz+lblok
	        lblok=nblock
        endif
        !write(6,*)'nshift, ncol', nshift, ncol
        do km = 1, nblock
          coordX(ncol,lpocz-1+km)=coordMp(km,1)
          coordY(ncol,lpocz-1+km)=coordMp(km,2)
          if (nshr > 1) then
              coordZ(ncol,lpocz-1+km)=coordMp(km,3)
              ElArea(ncol,lpocz-1+km)=charLength(km)**3
          else
              coordZ(ncol,lpocz-1+km)=0
              ElArea(ncol,lpocz-1+km)=charLength(km)**2
          endif 
          defgrad=0.d0
          stretchU=0.d0
          dstrain=0.d0
          do ic = 1,ndir
              defgrad(ic) = defgradNew(km,ic)
              stretchU(ic)= stretchNew(km,ic)
              dstrain(ic) = strainInc(km,ic)
          enddo
          if (nshr==1) then
              stretchU(4)= stretchNew(km,4)
              dstrain(4) = strainInc(km,4)
              defgrad(4) = defgradNew(km,4)
              defgrad(7) = defgradNew(km,5)
          elseif (nshr==3) then
              do ic=4,6
                  defgrad(ic) = defgradNew(km,ic)
                  stretchU(ic)= stretchNew(km,ic)
                  dstrain(ic) = strainInc(km,ic)
              enddo
              defgrad(7) = defgradNew(km,7)
              defgrad(8) = defgradNew(km,8)
              defgrad(9) = defgradNew(km,9)
          endif
          !call computeRotationR(stretchU,defgrad,ndir,nshr,R) ! obtain the rotation tensor R
          !call rotation(Dv,ndir,nshr,dstrain,R,dt,1) ! Rotate the strength rate to the fixed framework
          call getV6(Dv,ndir,nshr,dstrain,dt,1)
          Dlocal = 0
          do ic = 1,ndir
              Dlocal = Dlocal + Dv(ic)**2
          enddo
          do ic = 1,nshr
              Dlocal = Dlocal + 2*Dv(3+ic)**2
          enddo
          !write(6,*) 'Dlocal', Dlocal
          !bell shaped weight function is used here, correspondingly, alc=shear band size
          if (Dlocal<0) then
            dloc(ncol,lpocz-1+km)=0
          else
	          dloc(ncol,lpocz-1+km)=sqrt(Dlocal)
          endif
        enddo
        !write(6,*) 'coordX', coordX
        if (abs(totalTime-dt) <= dt) then 
          do km = 1, nblock
            DnlList(km)=dloc(nshift,lpocz-1+km)
          enddo
        else
!	************************************
          do km = 1,nblock
            totelWeight = 0.d0
            Dnl = 0.d0
            do j=1, i_elem
              distx= coordX(nshift,j)-coordX(nshift,lpocz-1+km)
              disty= coordY(nshift,j)-coordY(nshift,lpocz-1+km)
              distz= coordZ(nshift,j)-coordZ(nshift,lpocz-1+km)
              rdist=sqrt(distx*distx+disty*disty+distz*distz)
              !write(6,*) 'i,rdist', k,rdist
              if (rdist.le.Gcutoff) then
                weight = Gcoeff*exp(-(rdist/alc))  !Gauss weight function
!	            weight = (1 - rdist**2/alc**2)**2    ! Bell shaped weight function, no cut off
		            Dnl = Dnl+dloc(nshift,j)*weight*ElArea(nshift,j)
		            totelWeight=totelWeight+weight*ElArea(nshift,j)
	            endif
            enddo
            !write(6,*) 'i,Dnl,totelWeight', k,Dnl,totelWeight
	          !DnlList(k)=(1-nlbeta)*dloc(nshift,lpocz-1+k)+nlbeta*Dnl/totelWeight
            DnlList(km)=Dnl/totelWeight
          enddo
!	******************************
        endif
      endif
      
!      do i=1, ile_elem
!         write(6,*) 'i,dloc', i,dloc(1,i),dloc(2,i)
!      enddo
!      do i=1, nblock
!         write(6,*) 'i,dnonloc', i,DnlList(i)
!      enddo

!*************************************************************************
! Time integration
!*************************************************************************
      do 100 km = 1,nblock
        stress=0.d0
        defgrad=0.d0
        stretchU=0.d0
        dstrain=0.d0
        do ic = 1,ndir
            if (stressOld(km,ic).ge.0) then
                stress(ic)  = -10.00001
            else
                stress(ic)  = stressOld(km,ic)
            endif
            defgrad(ic) = defgradNew(km,ic)
            stretchU(ic)= stretchNew(km,ic)
            dstrain(ic) = strainInc(km,ic)
        enddo
        if (nshr==1) then
            stress(4)= stressOld(km,4)
            stretchU(4)= stretchNew(km,4)
            dstrain(4) = strainInc(km,4)
            defgrad(4) = defgradNew(km,4)
            defgrad(7) = defgradNew(km,5)
        elseif (nshr==3) then
            do ic=4,6
              stress(ic)  = stressOld(km,ic)
              defgrad(ic) = defgradNew(km,ic)
              stretchU(ic)= stretchNew(km,ic)
              dstrain(ic) = strainInc(km,ic)
            enddo
            defgrad(7) = defgradNew(km,7)
            defgrad(8) = defgradNew(km,8)
            defgrad(9) = defgradNew(km,9)
        endif
        
        ! obtain the rotation tensor R
        ! call computeRotationR(stretchU,defgrad,ndir,nshr,R)
        ! Rotate the strain rate to the fixed framework
        ! call rotation(Dv,ndir,nshr,dstrain,R,dt,1)
        ! Rotate the old stress to the fixed framework
        ! call rotation(Tv,ndir,nshr,stress,R,dt,0)
        
        call getV6(Tv,ndir,nshr,stress,dt,0)
        call getV6(Dv,ndir,nshr,dstrain,dt,1)

C        WRITE(*,*),'stressold=',stress
       
        
! void ratio from last/initial step 
        void = stateOld(km,1)
C         WRITE(*,*),'void ratio=',void
! suggested first time substep size of Richardson extrapolation
        !tsub = stateOld(km,2)
! get the coordinates
        !coords = (/coordMp(km,:)/)
        
        rotDiff = 0.d0
        if (nshr == 1) then
            rotDiff(2,1) = relSpinInc(km,1)/dt
            rotDiff(1,2) = -rotDiff(2,1)
        endif
        if (nshr == 3) then
            rotDiff(3,2) = relSpinInc(km,1)/dt
            rotDiff(2,3) = -rotDiff(3,2)
            rotDiff(1,3) = relSpinInc(km,2)/dt
            rotDiff(3,1) = -rotDiff(1,3)
            rotDiff(2,1) = relSpinInc(km,3)/dt
            rotDiff(1,2) = -rotDiff(2,1)
        endif
        
        !WRITE(6,*),'rotDiff=',rotDiff

! get vectors with 6 components (3D) of stress and strain vector
        !call getV6(Tv,ndir,nshr,stress,dt,0)
        ! get vectors with 6 components (3D) of strain rate vector
        !call computeStretchRate(Dv,ndir,nshr,defgrad,stretchU,dt)
        
        !WRITE(6,*),'Incremental=',relSpinInc
        !call getV6(Dv,ndir,nshr,dstrain,dt,1)
        !call getV6(Dv,ndir,nshr,stretchU,dt,1)
        !WRITE(6,*),'StretchRateTrue=',Dv

        call v2m(D,Dv)
        if (nonlocal_switch == 1 .and. totalTime > switch_time) then
          Dnonlocal = DnlList(km)
        else
          call Aij_Bij(D,D,TrD)
          Dnonlocal = sqrt(TrD)
        endif
        !call Aij_Bij(D,D,TrD)
        !Dlocal = sqrt(TrD)
        !WRITE(6,*),'Deformation Gradiant (F)=',defgrad
        !WRITE(6,*),'Dnonlocal=',Dnonlocal
        !WRITE(6,*),'Dlocal=',Dlocal
        !WRITE(6,*),'dstrain=',dstrain
        !WRITE(6,*),'new_stretchRate=',Dv
        !write(6,*),'void', void

!****************************************************
! local extrapolation based on forward Euler,
! variable substeps, consistent Jacobian 
! and nerror estimation
        !write(6,*),'dtsub', dtsub
        !if ((dtsub<=0.d0).or.(dtsub>dt)) dtsub = dt
        !write(6,*),'y=',y
        !write(6,*),'Tv=', Tv
        !write(6,*),'Dv=', Dv
        !write(6,*),'void, dtsub, Dnonlocal', void,dtsub,Dnonlocal
        !write(6,*),'B,G', B,G
        !write(6,*),'dt,tol,tolabs,hmin,theta', dt,tol,tolabs,hmin,theta
        !write(6,*),'props,nerror,naccst', props,nerror,naccst
        
        call rateT(TvR,Tv,Dv,voidR,void,props,nprops,tolabs,Dnonlocal,
     &             rotDiff, nerror, prinfo)
        if (nerror>=2) then
          return
        endif
        
        !WRITE(6,*),'Dv', Dv
        !WRITE(6,*),'Tv', Tv
        !WRITE(6,*),'TvR', TvR
        
        stress = Tv + TvR * dt
        void = void + voidR * dt
        
        ei0   = props(1)
        ed0   = props(2)
        hs    = props(5)
        en    = props(7)
        trT = Tv(1)+Tv(2)+Tv(3)
        edd = ed0*exp(-(-trT/hs)**en)*1.001
        eii = ei0*exp(-(-trT/hs)**en)*0.999
        if (void.le.edd) then 
          void=edd
        end if
        if (void.ge.eii) then
          void=eii
        end if

!        nsub=1
!        call euler(y,ny,Tv,void,Dv,B,G,nBG,dt,tolabs,Dnonlocal,
!&             rotDiff,nsub,theta,props,nprops,nerror,prinfo)
!        naccst = 0
!        call richardson(y,ny,Tv,Dv,void,B,G,nBG,Dnonlocal,dtsub,
!     &                      dt,tol,tolabs,hmin,theta,rotDiff,
!     &                      props,nprops,nerror,prinfo,naccst)
        
        if ((nerror>0).and.(prinfo)) call wrista(2,y,ny,Dv,dt,coords)    
        if (nerror==2 .or. nerror==3 ) then
          call XPLB_EXIT
          return
        elseif (nerror==4) then
          call XPLB_EXIT
        endif
        ! Rotate the old stress to the fixed framework
        !do ic = 1,ndir+nshr
        !    stress(ic)  = y(ic)
        !    !stressNew(km,ic)  = y(ic)
        !end do
        !call rotationBack(Tv,ndir,nshr,stress,R)   
        
! updated stresses
!        do ic = 1,ndir
!            if (Tv(ic).ge.0) then
!                Tv(ic) = -10.00001
!            endif
!        end do
        do ic = 1,ndir+nshr
            !stressNew(km,ic)  = Tv(ic)
            if (stress(ic).ge.0 .and. ic <= ndir) then
                stressNew(km,ic)  = -10.00001
            else
                stressNew(km,ic)  = stress(ic)
            endif
        end do
        
! updated vector of additional state variables to abaqus statev vector
              
        stateNew(km,1) = void
       !tateNew(km,2) = dtsub
       !stateNew(km,3) = naccst
        
! output variable to check error         
        if (0) then
          WRITE(6,*),'Guass Point Number', km
          WRITE(6,*),'Void ratio e=', void
          WRITE(6,*),'stress=', stress
          !call wrista(1,y,ny,Dv,dtsub,coords)
        endif

  100 continue

      return
      end
!-----------------------------------------------------------------------------
!---------------------Subroutines for VUMAT files
!-----------------------------------------------------------------------------
      subroutine richardson(y,ny,Tv,Dv,void,B,G,nBG,Dnonlocal,h,
     &                      time,tol,tolabs,hmin,theta,rotDiff,
     &                      props,nprops,nerror,prinfo,naccst)
!-----------------------------------------------------------------------------
!  numerical solution of y'=f(y)
!  forward and backward Euler method, with local extrapolation 
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      integer  nprops, nerror, ny, nBG, i, naccst, nrejst
!      real*8 void, h, time, hmin, err, h2, hsav, actt,  
!     &  theta, tol, tolabs, fhnew, void2, hout,Dnonlocal,
!     &  Tv, Dv, B, G, Tv2, props, y, y1, y2,
!     &  B2, G2, sc, errt, est, rotDiff
      
      dimension Tv(6), Dv(6), B(6,6), G(6), rotDiff(3,3),
     &  Tv2(6), props(nprops), y(ny), y1(ny), y2(ny),
     &  B2(6,6), G2(6), sc(6), errt(6), est(ny) 
      
      logical prinfo
      character(11), parameter :: subr = 'richardson'

      if ((h<0.d0).or.(h>time)) h = time 
      nerror = 0
      actt = 0
      
      if (prinfo) write(6,*) 'time integration starts with h = ',h
      if (prinfo) write(6,*) 'end time', time

      do while (actt<time)
        if (h<hmin) then
          call errorcheck(nerror,3,'h     ',
     &                    h,3,'hmin  ',hmin,subr,prinfo)
          write(6,*) 'error h',h
          return
        endif

        call euler(y1,ny,Tv,void,Dv,B,G,nBG,h,tolabs,Dnonlocal,
     &             rotDiff,1,theta,props,nprops,nerror,prinfo)

        if (nerror>=2) then
          write(6,*) 'error in euler step with stepsize h'
          return
        endif
      
! two Euler steps with step size h/2
        h2 = h/2
        call euler(y2,ny,Tv,void,Dv,B,G,nBG,h2,tolabs,Dnonlocal,
     &             rotDiff,1,theta,props,nprops,nerror,prinfo)

        if (nerror>=2) then
          write(6,*) 'error in first euler step with stepsize h/2'
          return
        endif

        call getTQBG(y2,ny,Tv2,void2,B2,G2,nBG)
     
        call euler(y2,ny,Tv2,void2,Dv,B2,G2,nBG,h2,tolabs,Dnonlocal,
     &             rotDiff,1,theta,props,nprops,nerror,prinfo)

        if (nerror>=2) then
          write(6,*) 'error in second euler step with stepsize h/2'
          return
        endif
        
!     only the stress components are used for error estimation         
        do i=1,ny
          est(i) = (y2(i)-y1(i))
          if (i<=6) then 
            sc(i) = tolabs + max( abs(y(i)), abs(2*y2(i)-y1(i)) )
            errt(i) = abs(est(i)/sc(i))
          endif
        enddo
        err = maxval(abs(errt))
        err = max(err,1.d-16)

! calculate factor of step change (times a safety factor)

        fhnew = 0.90d0*sqrt(tol/err)
      
        if (err>tol) then

! do not accept step size, reject step

          nrejst = nrejst + 1
          if (prinfo) write(6,*) 'reject: h, actt', h,actt
        
          h = max(0.2d0, fhnew)*h
        else
!
! update and try to calculate time rate at the end of the substep
! second oder update
!
          y = y2 + est
          call getTQBG(y,ny,Tv,void,B,G,nBG)
       
! updated state variables can be handled by constitutive law: accept
! step size and update state variables, reuse the above objective time
! rate of state variables for new step
          actt = actt + h
          naccst = naccst + 1
          if (prinfo) write(6,*) 'accept: h, actt', h,actt

! suggested new step size limited by a factor of 2
          hsav = h        
          h = min(2.d0,fhnew)*h
          if ((actt+h>time).and.(actt<time)) then
! substep size for next time step 
! first suggested substep for next increment, reduced by a safety factor
            hout = 0.8d0*hsav
            h = time - actt
            if (abs(h)<hmin) then
                exit
            endif
          endif
        endif  
      enddo
      
      h = hout
      if (prinfo) then 
            write(6,*) 'end of time integration'
            write(6,*) 'number of accepted substeps: ', naccst
            write(6,*) 'number of rejected substeps: ', nrejst
            write(6,*) 'time integration proposes ',
     &                           'new start h = ', hsav
      endif

      return
      end       

! ============================================================================
! Subroutines for the time integration process
! ============================================================================
! -----------------------------------------------------------------------------
      subroutine euler(y,ny,Tv,void,Dv,B,G,nBG,time,tolabs,Dnonlocal,
     &                 rotDiff,nsub,theta,props,nprops,nerror,prinfo)

!-----------------------------------------------------------------------------
!
!  numerical solution of y'=f(y) forward Euler method 
      include 'vaba_param.inc'
      
      integer ny, nprops, nerror, nBG, i, j, nsub, 
     &  E(6,6), nrelop
 
      dimension TR(6), Tv(6), TRB(6), TvB(6), Tv1(6), Dv(6), DvB(6), 
     &  B(6,6), G(6), Bp(6,6), Gp(6), B1(6,6), G1(6), props(nprops),
     &  rotDiff(3,3), y(ny)
!      real*8, intent(out) :: y(ny)
      logical prinfo
      character(6) var1, var2
      character(11), parameter :: subr = 'euler'
    
      nerror = 0
      E = 0.d0
      
      do i=1,6
        E(i,i) = 1
      enddo
       
      h = time/nsub

      do i=1,nsub
        call rateT(TR,Tv,Dv,voidR,void,props,nprops,tolabs,Dnonlocal,
     &             rotDiff, nerror, prinfo)
        
        if (nerror>=2) then
          return
        elseif (nerror==1) then
          if (prinfo) then
            air = real(i)
            call errorcheck(nerror,2,var1,air,nrelop,'stress',var2v,
     &                      subr, prinfo)      
          endif
        endif
        
!        do j=1,6
!	    TvB = Tv+theta*B(:,j)
!          voidB = void+theta*G(j)
!          if (j<=ntens) then
!            DvB = Dv+theta*E(:,j)
!          else
!            DvB = Dv
!          endif
!          call rateT(TRB,TvB,DvB,voidRB,voidB,props,nprops,tolabs,
!     &               Dnonlocal, rotDiff, nerror, prinfo)     
!          Bp(:,j) = (TRB - TR)/theta
!          Gp(j) = (voidRB - voidR)/theta
!        enddo
        
        Tv1 = Tv + TR * h
        void1 = void + voidR * h
        B1 = B !+ Bp * h
        G1 = G !+ Gp * h
        
        ei0   = props(1)
        ed0   = props(2)
        hs    = props(5)
        en    = props(7)
        trT = Tv1(1)+Tv1(2)+Tv1(3)
        edd = ed0*exp(-(-trT/hs)**en)*1.001
        eii = ei0*exp(-(-trT/hs)**en)*0.999
        if (void1.le.edd) then 
	    void1=edd
        end if
        if (void1.ge.eii) then 
	    void1=eii
        end if
        if (nsub>1) then
          Tv = Tv1
          Qv = Qv1
          B = B1
          G = G1        
        endif
      enddo
      
!      y = (/ Tv1, void1 /)
      call getY(y,ny,Tv1,void1,B1,G1,nBG)
      
      return
      end

! ============================================================================
! Following subroutines have to be adapted when changing the constitutive law
! ============================================================================
!-----------------------------------------------------------------------------
      subroutine rateT(TR,Tv,Dv,voidR,void,props,nprops,tol,Dnonlocal,
     &                 rotDiff, nerror, prinfo)
!-----------------------------------------------------------------------------
! calculate objective time rate of stresses and additional state variables
! abaqus requires Green-McInnis-Naghdi stress rate
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      integer nprops, nerror, i, j
      logical prinfo
      character(11), parameter :: subr = 'rateT'

!     &      voidR, void, tol, valD, trT, term1, term2,
!     &	  term3, term4, term5, phic, hs, en, ed0, ec0, ei0,
!     &      alpha, beta, alc, trTsv2, trTsv3, sinphi, sq2,
!     &      sq3, sq6, c3t, c1, c2, a1, trTs2, ed, ec, rotDiff,
!     &      ei, fs, fd, trD, Dnonlocal, TR, Tv, Dv, props,
!     &      Ts, Ts2,  Tsv, Tsv2, Tsv3, D, T, TR1,
!     &      rateTerm2, rateTerm3, Temp
	 
      dimension TR(6), Tv(6), Dv(6), props(nprops), 
     &      FLL(3,3), FNN(3,3),  Ts(3,3), rotDiff(3,3),
     &      Ts2(3,3),  Tsv(3,3), Tsv2(3,3), Tsv3(3,3),
     &      D(3,3), T(3,3), TR1(3,3), rateTerm2(3,3),
     &      rateTerm3(3,3), Temp(3,3)
      parameter ( PI=3.141592653589793d0, zero=1d-10 )
      
      nerror = 0

      ei0   = props(1)
      ed0   = props(2)
      ec0   = props(3)
      phic  = props(4)*PI/180
      hs    = props(5)
      beta  = props(6)
      en    = props(7)
      alpha = props(8)
      !lc    = props(9)
 
!      if (nprops<9) then
!        npr = real(nprops)
!        call errorcheck(nerror,4,'nprops',npr,3,'9       ',9.d0,subr,
!     &                  prinfo)   
!        return
!     endif
! time derivative of additional state parameter void ratio
      call v2m(D,Dv)
      !write(6,*),'Dv=',Dv
      trD = Dv(1)+Dv(2)+Dv(3)
      voidR = ( 1 + void )*trD

      call v2m(T,Tv)
      trT = Tv(1)+Tv(2)+Tv(3)

!     ed = ed0*exp(-(-trT/hs)**en)
!     ec = ec0*exp(-(-trT/hs)**en)
!     ei = ei0*exp(-(-trT/hs)**en)
      
! check user input on severe errors
! which would lead to a crash in evaluation of constitutive law 
      if (phic<=zero) 
     &  call errorcheck(nerror,4,'phic  ',phic,5,'0     ',zero,subr,
     $                  prinfo)
      if (hs<=zero) 
     &  call errorcheck(nerror,4,'hs   ',hs,5,'0     ',zero,subr,prinfo)
      if (ed0<=zero) 
     &  call errorcheck(nerror,4,'ed0   ',ed0,5,'0     ',zero,subr,
     &                  prinfo)
      if (ec0<=ed0) 
     &  call errorcheck(nerror,4,'ec0   ',
     &                  ec0,5,'ed0   ',ed0,subr,prinfo)
      if (ei0<=ec0) 
     &  call errorcheck(nerror,4,'ei0   ',
     &                  ei0,5,'ec0   ',ec0,subr,prinfo)
! check actual void ratio
      if (void<=0) 
     &  call errorcheck(nerror,4,'e     ',void,5,'0     ',zero,subr,
     &                  prinfo)
! time derivative of stresses: hypoplastic law
! calculate L and N as functions of T, e
      
! For zero stress --> hypoplastic compression law
      if (abs(trT)<tol) then
        TR = 0.d0
        if ((trD<0).and.(void<=ei0)) then
! derivative of compression law
            TR(1) = voidR/void/en/3 * tol / ( tol/hs )**en
            TR(2) = TR(1)
            TR(3) = TR(1)
        endif
        return
      endif
      
      !write(6,*),'substressT=',Tv
! undefined stress state
      if (trT>=0) then
        call errorcheck(nerror,1,'trT   ',trT,4,'0     ',zero,subr,
     &                  prinfo)   
        return
      endif
! calculate relative stress: Ts = T / trac(T)
      term1 = 1.d0/trT
      do i=1,3
        do j=1,3
            Ts(i,j) = T(i,j)*term1
        enddo
      enddo

! calculate square of relative stress: Ts2 = Ts * Ts
      call Aik_Bkj(Ts,Ts,Ts2)
      trTs2  = Ts2(1,1)+Ts2(2,2)+Ts2(3,3)
      
      !write(6,*),'Ts2=',Ts2

! calculate deviator: Tsv = Ts - 1/3 I
      do i=1,3
        do j=1,3
            Tsv(i,j) = Ts(i,j)
        enddo
        Tsv(i,i) = Tsv(i,i) - 1.0d0/3.0d0
      enddo

      call Aik_Bkj(Tsv,Tsv,Tsv2)
      call Aik_Bkj(Tsv2,Tsv,Tsv3)
      
      !write(6,*),'Tsv2=',Tsv2
      !write(6,*),'Tsv3=',Tsv3
      
! calculate norm of relative deviatoric stress, ||Tsv||
      trTsv2 = Tsv2(1,1)+Tsv2(2,2)+Tsv2(3,3)
      trTsv3 = Tsv3(1,1)+Tsv3(2,2)+Tsv3(3,3)
      sinphi = sin(phic)
      sq2    = sqrt(2.0d0)
      sq3    = sqrt(3.0d0)
      sq6    = sqrt(6.0d0)

      if (trTsv2<=1.d-10) then
        c3t = 1.0d0
      else
        c3t = -sq6*trTsv3/trTsv2**1.5d0
        if (c3t > 1.0d0)  c3t =  1.0d0
        if (c3t < -1.0d0) c3t = -1.0d0
      endif
      c1 = sq3*( 3.0d0-sinphi )/( 2.0d0*sq2*sinphi )
      c2 = 3*( 3.0d0+sinphi )/( 8.0d0*sinphi )
      a1 = 1.0d0/(c1+c2*sqrt(trTsv2)*(1+c3t))
      
      !write(6,*),'c1, c2, a1, c3t', c1, c2, a1, c3t

      ed = ed0*exp(-(-trT/hs)**en)
      ec = ec0*exp(-(-trT/hs)**en)
      ei = ei0*exp(-(-trT/hs)**en)

! check actual void ratio
      if (void<ed) call warning('e     ',void,3,'ed    ',ed,set2,subr,
     &                          prinfo)
      if (void>ei) call warning('e     ',void,1,'ei    ',ei,set2,subr,
     &                          prinfo)

      term2 = ( (ei0-ed0)/(ec0-ed0) )**alpha
      term3 = ( 1/c1**2.0d0 + 1/3.0d0 - term2/(c1*sq3))**(-1.d0)
      term4 = hs/en*( ei/void )**beta*( 1+ei )/ei*( -trT/hs )**(1.d0-en)

      fs  = term3*term4
      
      fd = ( void-ed )/( ec-ed )
      if (fd>0.0d0) then
        fd = fd**alpha
      else
        fd = 0.0d0
      endif

      call Aij_Bij(Ts,D,term5)

      !write(6,*),'term2,term3,term4,term5,fs',term2,term3,term4,term5,fs
      !write(6,*),'a1,D,Ts',a1,D,Ts
      do i=1,3
        do j=1,3
            FLL(i,j) = a1**2*D(i,j) + term5*Ts(i,j)
            FNN(i,j) = a1*(Ts(i,j)+ Tsv(i,j))
        enddo
      enddo

! calculate  TR_ij = L_ijkl*D_kl+N_ij*||D||
      do i=1,3
        do j=1,3
            TR1(i,j) = fs*(FLL(i,j)+fd*Dnonlocal*FNN(i,j))
        enddo
      enddo
      !write(6,*), 'fs=', fs
      !write(6,*), 'LL=', FLL
      !write(6,*), 'NN=', FNN
      !write(6,*), 'Jaumann rate=', TR1
      
      !write(6,*), 'rotDiff=', rotDiff
      !write(6,*), '-rotDiff=', -rotDiff
      !adjust the Jaumann rate to Green-Naghdi rates
      
      !call Aik_Bkj(rotDiff,T,rateTerm2)
      !call Aik_Bkj(T,-rotDiff,rateTerm3)
      !call Aij_plus_Bij(TR1,rateTerm2,Temp)
      !call Aij_plus_Bij(Temp,rateTerm3,TR1)
      
      !write(6,*), 'Green-Naghdi rate rate=', TR1
      
      call m2v(TR,TR1)
      !write(6,*), 'stress rate=', TR

      return
      end
      
!-----------------------------------------------------------------------------
      subroutine getTQBG(y,ny,Tv,void,B,G,nBG)
!-----------------------------------------------------------------------------
! split super vector y into vectors T and Q and matrices B and G 
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      integer ny, i, nBG
!      real*8 void, y, BG, B, G,Tv
      dimension y(ny), BG(nBG), B(6,6), G(6),Tv(6)
      
      Tv = y(1:6)
      void = y(7)
      BG = y(8:ny)
      
      do i=1,6
        B(:,i) = BG( (i-1)*7+1 : (i-1)*7+6 )
        G(i) = BG( i*7 )
      enddo

      return
      end      
    
!-----------------------------------------------------------------------------
      subroutine getY(y,ny,Tv,void,B,G,nBG)
!-----------------------------------------------------------------------------
! transform input vectors and matrices into supervector y form 
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      integer ny, i, nBG
      dimension y(ny), BG(nBG),B(6,6), G(6),Tv(6)
!      real*8 y, BG,B, G,Tv,void

      do i=1,6
        BG(1+(i-1)*7 : i*7) = (/B(:,i), G(i)/)
      enddo

      y = (/ Tv, void, BG /)

      return
      end

!-----------------------------------------------------------------------------
      subroutine  computeRotationR(stretchUv,defgradv,ndi,nshr,Rt)
!-----------------------------------------------------------------------------
! compute rotation tensor
!-----------------------------------------------------------------------------    
      include 'vaba_param.inc'
      dimension stretchUv(6),defgradv(9),U(3,3),F(3,3),Uinv(3,3),Rt(3,3)
!      real stretchUv,defgradv,F,Rt,U,Uinv
      integer ndi, nshr, i
      
      F=0.d0
      U=0.d0
      
      !write(6,*),'defgrad=',defgradv
      
      ! do i=1,3
          ! do j=1,3
              ! F(i,j)=0.0d0
              ! U(i,j)=0.0d0
          ! end do 
      ! end do
      do i=1,3
          F(i,i) = defgradv(i)
          U(i,i) = stretchUv(i)
      enddo

      if (nshr==1) then
          U(1,2) = stretchUv(4)
          U(2,1) = stretchUv(4)
          
          F(1,2) = defgradv(4)
          F(2,1) = defgradv(7)      
      endif
      if (nshr==3) then
          U(1,2) = stretchUv(4)
          U(2,1) = stretchUv(4)
          U(2,3) = stretchUv(5)
          U(3,2) = stretchUv(5)
          U(3,1) = stretchUv(6)
          U(1,3) = stretchUv(6)
          
          F(1,2) = defgradv(4)
          F(2,3) = defgradv(5)
          F(3,1) = defgradv(6)
          F(2,1) = defgradv(7)
          F(3,2) = defgradv(8)
          F(1,3) = defgradv(9)
      endif
      
      !write(6,*),'F=',F
      !write(6,*),'U=',U

      CALL INVERSE(U,3,3,Uinv)
      call Aik_Bkj(F,Uinv,Rt)
      
      !write(6,*),'Uinv=',Uinv
      !write(6,*),'R=',Rt
      
      end

!-----------------------------------------------------------------------------
      subroutine rotation(Dv,ndi,nshr,V,R,dtime,fl)
!-----------------------------------------------------------------------------
! compute stretch tensor
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      dimension Dv(6),V(6),vTensor(3,3),R(3,3),RT(3,3),
     &     rTensor(3,3),vTemp(3,3)
!      real*8 Dv,V,R,RT,dtime,vTensor,rTensor,vTemp
      integer ndi, nshr, i, fl
      
      vTensor = 0.d0
      !for the intial step
      if (dtime.le.1.0E-20) then 
          dtime = 1.0
      endif

      do i=1,ndi
          if (fl==1) then
! use small strain to calculate the rate of stretching tensor
            vTensor(i,i) = V(i)/dtime
          else
            vTensor(i,i) = V(i)
          endif
      enddo

      if (nshr==1) then
          if (fl==1) then
            vTensor(1,2) = V(4)/dtime
            vTensor(2,1) = V(4)/dtime
          else
            vTensor(1,2) = V(4)
            vTensor(2,1) = V(4)
          endif
      endif
      if (nshr==3) then
          if (fl==1) then
              vTensor(1,2) = V(4)/dtime
              vTensor(2,1) = V(4)/dtime
              vTensor(2,3) = V(5)/dtime
              vTensor(3,2) = V(5)/dtime
              vTensor(3,1) = V(6)/dtime
              vTensor(1,3) = V(6)/dtime
          else
              vTensor(1,2) = V(4)
              vTensor(2,1) = V(4)
              vTensor(2,3) = V(5)
              vTensor(3,2) = V(5)
              vTensor(3,1) = V(6)
              vTensor(1,3) = V(6)
          endif
      endif

      call getTrans(R,RT)
      call Aik_Bkj(R,vTensor,vTemp)
      call Aik_Bkj(vTemp,RT,rTensor)
      
   
      Dv = 0.d0
      Dv(1) = rTensor(1,1)
      Dv(2) = rTensor(2,2)
      Dv(3) = rTensor(3,3)
      Dv(4) = (rTensor(1,2)+rTensor(2,1))/2.0d0
      Dv(5) = (rTensor(2,3)+rTensor(3,2))/2.0d0
      Dv(6) = (rTensor(1,3)+rTensor(3,1))/2.0d0

      return
      end
!-----------------------------------------------------------------------------
      subroutine rotationBack(Dv,ndi,nshr,V,R)
!-----------------------------------------------------------------------------
! compute stretch tensor
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      dimension Dv(6),V(6),vTensor(3,3),R(3,3),RT(3,3),
     &     rTensor(3,3),vTemp(3,3)
!      real*8 Dv,V,R,RT,vTensor,rTensor,vTemp
      integer ndi, nshr, i
      
      vTensor = 0.d0
      do i=1,ndi
          vTensor(i,i) = V(i)
      enddo

      if (nshr==1) then
          vTensor(1,2) = V(4)
          vTensor(2,1) = V(4)
      endif
      if (nshr==3) then
          vTensor(1,2) = V(4)
          vTensor(2,1) = V(4)
          vTensor(2,3) = V(5)
          vTensor(3,2) = V(5)
          vTensor(3,1) = V(6)
          vTensor(1,3) = V(6)
      endif

      call getTrans(R,RT)
      call Aik_Bkj(RT,vTensor,vTemp)
      call Aik_Bkj(vTemp,R,rTensor)
      
      Dv = 0.d0
      Dv(1) = rTensor(1,1)
      Dv(2) = rTensor(2,2)
      Dv(3) = rTensor(3,3)
      Dv(4) = (rTensor(1,2)+rTensor(2,1))/2.0d0
      Dv(5) = (rTensor(2,3)+rTensor(3,2))/2.0d0
      Dv(6) = (rTensor(1,3)+rTensor(3,1))/2.0d0

      return
      end
!----------------------------------------------------------------------------- 
      subroutine getTrans(R,RT)
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
!      real*8 R,RT
      dimension R(3,3), RT(3,3)
      
      RT = 0.d0
      do i=1,3
          RT(i,i) = R(i,i)
      enddo
      
      if (nshr==1) then
           RT(1,2) = R(2,1)
           RT(2,1) = R(1,2)
      endif
      
      if (nshr==3) then
           RT(1,2) = R(2,1)
           RT(2,1) = R(1,2)
           RT(1,3) = R(3,1)
           RT(3,1) = R(1,3)
           RT(2,3) = R(3,2)
           RT(3,2) = R(2,3)  
      endif
      end
      
!-----------------------------------------------------------------------------
      subroutine getV6(V6,ndi,nshr,V,dtime,fl)
!-----------------------------------------------------------------------------
! transform input vectors into 6x1 vector form
! fl...flag: 1 ... strain vector, for stretch rate
!            0 ... arbitrary vector
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      integer ndi, nshr, fl, i, j
!      real*8 dtime, V, V6
      dimension V(ndi+nshr), V6(6)

      V6 = 0.d0
      !for the intial step
      if (dtime.le.1.0E-20) then 
          dtime = 1.0
      end if
      
      !WRITE(6,*),'check time=',dtime

      do i=1,ndi
          if (fl==1) then
! use small strain to calculate the rate of stretching tensor
            V6(i) = V(i)/dtime
        else
            V6(i) = V(i)
        endif
      enddo

      if (nshr>=1) then
        if (fl==1) then
            V6(4) = V(ndi+1)/dtime
        else
            V6(4) = V(ndi+1)
        endif
      endif
      if (nshr>=2) then
        if (fl==1) then
            V6(5) = V(ndi+2)/dtime
        else
            V6(5) = V(ndi+2)
        endif
      endif
      if (nshr>=3) then
        if (fl==1) then
            V6(6) = V(ndi+3)/dtime
        else
            V6(6) = V(ndi+3)
        endif
      endif

      return
      end

!-----------------------------------------------------------------------------
      subroutine v2m(T,V)
!-----------------------------------------------------------------------------
! transform 6x1 vector into 3x3 matrix form
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      dimension V(6), T(3,3)
!      real*8 V,M

      T(1,1) = V(1)
      T(2,2) = V(2)
      T(3,3) = V(3)
      T(1,2) = V(4)
      T(2,1) = V(4)
      T(2,3) = V(5)
      T(3,2) = V(5)
      T(1,3) = V(6)
      T(3,1) = V(6)

      return
      end

!-----------------------------------------------------------------------------
      subroutine m2v(V,T)
!-----------------------------------------------------------------------------
! transform 3x3 matrix into 6x1 vector form
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      dimension V(6), T(3,3)
!      real*8 V,T

      V(1) = T(1,1)
      V(2) = T(2,2)
      V(3) = T(3,3)
      V(4) = (T(1,2)+T(2,1))/2.0d0
      V(5) = (T(2,3)+T(2,3))/2.0d0
      V(6) = (T(3,1)+T(1,3))/2.0d0

      return
      end

!-----------------------------------------------------------------------------
      subroutine Aij_Bij(A,B,C)
!-----------------------------------------------------------------------------
! inner product of a second order tensor T with a second order tensor D
! T:D = T_ij D_ij
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      dimension A(3,3), B(3,3)
!      real*8 A, B, C
      integer i, j 

      C = 0.d0
      do i=1,3
        do j=1,3
            C = C + A(i,j)*B(i,j)
        end do
      end do

      return
      end

!-----------------------------------------------------------------------------
      subroutine Aik_Bkj(A,B,C)
!-----------------------------------------------------------------------------
! inner product of a second order tensor A with a second order tensor B
! (matrix multiplication)
! T.D = T_ij D_jk
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      dimension A(3,3), B(3,3), C(3,3)
      integer i, j, k
!      real*8 A, B, C

      do i=1,3
        do j=1,3
            C(i,j) = 0.d0
            do k=1,3
                C(i,j) = C(i,j) + A(i,k)*B(k,j)
            end do
        end do
      end do

      return
      end

      subroutine Aij_plus_Bij(A,B,C)
!-----------------------------------------------------------------------------
! inner product of a second order tensor A with a second order tensor B
! (matrix multiplication)
! T.D = T_ij D_jk
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      dimension A(3,3), B(3,3), C(3,3)
      integer i, j, k
!      real*8 A, B, C

      do i=1,3
        do j=1,3
            C(i,j) = A(i,j)+B(i,j)
        end do
      end do

      return
      end

!-----------------------------------------------------------------------------
      subroutine subroutine Aij_Bkl(A,B,C)
!-----------------------------------------------------------------------------
! tensorial product of a second order tensor T with a second order tensor D
! T x D = T_ij D_kl
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      dimension A(3,3), B(3,3), C(3,3,3,3)
      integer i, j, k, l
!      real*8 A, B, C

      do i=1,3
        do j=1,3
            do k=1,3
                do l=1,3
                    C(i,j,k,l) = A(i,j)*B(k,l)
                enddo
            enddo
        enddo
      enddo

      return
      end
      
      SUBROUTINE INVERSE(A,N,NP,AINV)
!========================================================================
!
!    CALCULATE THE SECOND ORDER TENSOR A'S INVERSE, AINV
!    A^{-1} = AINV    
!    this subroutine inverses a (n x n) A matrix
!	 following a Gauss-Jordan elimination process
!
!========================================================================
      INCLUDE 'vaba_param.inc'
      integer N,NP
!     &    A,IPIV,INDXR,INDXC,A0,AINV
!  
      dimension A(NP,NP),IPIV(NP),INDXR(NP),INDXC(NP),
     &  A0(NP,NP),AINV(NP,NP)
!
      DO J=1,N
        IPIV(J)=0
      END DO
!
!     storage of the original A matrix
!
      DO I=1,N
        DO J=1,N
          A0(I,J)=A(I,J)
        END DO
      END DO
!
!	find a pivot (largest absolute value) among the rows of A that have not already been reduced
!
      DO I=1,N
        BIG=0.0D0
        DO J=1,N
          IF(IPIV(J).NE.1)THEN
            DO K=1,N
                IF(IPIV(K).EQ.0)THEN
                  IF(ABS(A(J,K)).GE.BIG)THEN
                    BIG=ABS(A(J,K))
                    IROW=J
                    ICOL=K
                    PIV=A(J,K)
                  END IF
                ELSEIF(IPIV(K).GT.1)THEN
                  write (7,*) 'Singular Matrix'
                END IF
            END DO
          END IF
        END DO
!
        IPIV(ICOL)=IPIV(ICOL)+1
        INDXR(I)=IROW
        INDXC(I)=ICOL
!	  
!     interchange the rows to put the pivot on the diagonal
!
        IF(IROW.NE.ICOL)THEN
          DO L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
          END DO
        END IF
!
!     reduction of the row of the pivot
!       
        IF(PIV.EQ.0) write (7,*) 'Singular Matrix2'
!       
        PIVIN=1./PIV          ! numerical stabilization
!
        A(ICOL,ICOL)=1.       ! numerical stabilization
        DO L=1,N
           A(ICOL,L)=A(ICOL,L)*PIVIN
        END DO
!
!     reduction of the column of the pivot
!
        DO LL=1,N
           IF(LL.NE.ICOL)THEN
             DUM=A(LL,ICOL)
             A(LL,ICOL)=0.    ! numerical stabilization
             DO L=1,N
                A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
             END DO
           END IF
        END DO
      END DO
!
!     unscramble the columns to get A^{-1}
!		
      DO J=N,1,-1   ! reverse DO loop
        DO K=1,N
          DUM=A(K,INDXR(J))
          A(K,INDXR(J))=A(K,INDXC(J))
          A(K,INDXC(J))=DUM
        END DO
      END DO
!
!	restitution process of A and Ainv
!
      DO I=1,N
        DO J=1,N
          AINV(I,J)=A(I,J)
        END DO
      END DO
!
      DO I=1,N
        DO J=1,N
          A(I,J)=A0(I,J)
        END DO
      END DO
!     
      RETURN
      END
!
!-----------------------------------------------------------------------------
      subroutine errorcheck(nerror,mode,var1,var1v,nrelop,var2,var2v,
     & subr, prinfo)
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      integer mode, nrelop
      integer, intent(out) :: nerror    
!      real*8 var1v, var2v
      logical prinfo
      character(11)  subr
      character(6) var1, var2
      
      if (prinfo) then 
      write(6,*) '==================================================='
      write(6,*) 'ERROR: abaqus job failed during call of vumat'
      write(6,*) '==================================================='
      
      write(6,*) 'ERROR in subroutine ', subr
      
      if (nrelop==0) then
        write(6,*) 'cannot evaluate ', var1
      elseif (nrelop==1) then
        write(6,112) var1,'=',var1v,' is greater than     ',
     &               var2,'=',var2v
        write(6,*) var1,' must be less or equal ',var2
      elseif (nrelop==2) then
        write(6,112) var1,'=',var1v,'is equal to         ',
     &               var2,'=',var2v
        write(6,*) var1,' must be unequal to ',var2
      elseif (nrelop==3) then
        write(6,112) var1,'=',var1v,' is less than        ',
     &               var2,'=',var2v
        write(6,*) var1,' must be greater or equal ',var2
      elseif (nrelop==4) then
        write(6,112) var1,'=',var1v,' is greater or equal ',
     &               var2,'=',var2v
        write(6,*) var1,' must be less than ',var2
      elseif (nrelop==5) then
        write(6,112) var1,'=',var1v,' is less or equal    ',
     &               var2,'=',var2v
        write(6,*) var1,' must be greater than ',var2
      endif
      endif
      
      if (mode==1) then
      if (prinfo) then 
        write(6,*) 'ERROR 1'
        write(6,*) 'Problems in evaluation of the time rate ',
     &             '(e.g. undefined stress state)'
      endif   
        nerror = 1
      elseif (mode==2) then
      if (prinfo) then 
        write(6,*) 'ERROR 2'
        write(6,*) 'Problems in time integration: ',var2,' rate ',
     &             'undefined in substep',int(var1v)
      endif
        nerror = 2
      elseif (mode==3) then
      if (prinfo) then 
        write(6,*) 'ERROR 3'
        write(6,*) 'Problems in time integration, ', 
     &             'stepsize smaller than minimum specified'
      endif
        nerror = 3
      elseif (mode==4) then
        write(6,*) 'ERROR 4'
        write(6,*) 'fatal error'
        write(6,*) 'Program terminated.'
        nerror = 4
      endif
      
112   format(1X,a6,a2,E9.2,a17,a6,a2,E9.2)      
      
      return
      end      
      
!-----------------------------------------------------------------------------
      subroutine wrista(mode,y,ny,Dv,dtime,coords)
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      integer mode, ny    
      dimension y(ny), Dv(6), coords(3)
       
      if (mode==1) then
        write(6,*) '==================================================='
        write(6,*) 'Call of vumat:'
        write(6,*) '==================================================='

      elseif (mode==2) then
        write(6,*) '==================================================='
        write(6,*) 'ERROR: abaqus job failed during call of vumat'
        write(6,*) '==================================================='
        write(6,*) 'state dumb:'
        write(6,*) 
      endif

        write(6,*) 'Coordinates of material point:'
        write(6,104) 'x1 = ',coords(1),' x2 = ',coords(2),' x3 = ',
     &    coords(3)

        write(6,*) 'Stresses:'
        write(6,*) 
        write(6,106) 'T(1,1) = ',y(1),'T(1,2) = ',y(4),'T(1,3) = ',
     &    y(5)
        write(6,106) 'T(2,1) = ',y(4),'T(2,2) = ',y(2),'T(2,3) = ',
     &    y(6)
        write(6,106) 'T(3,1) = ',y(5),'T(3,2) = ',y(6),'T(3,3) = ',
     &    y(3)
        write(6,*) 
        write(6,*) 'Stretching rate:'
        write(6,*) 
        write(6,101) 'D(1,1) = ',Dv(1),'D(1,2) = ',Dv(4),'D(1,3) = ',
     &    Dv(5)
        write(6,101) 'D(2,1) = ',Dv(4),'D(2,2) = ',Dv(2),'D(2,3) = ',
     &    Dv(6)
        write(6,101) 'D(3,1) = ',Dv(5),'D(3,2) = ',Dv(6),'D(3,3) = ',
     &    Dv(3)
        write(6,*) 
        write(6,*) 'Time increment:'
        write(6,*) 
        write(6,108) 'dtime = ',dtime
        write(6,*) 
        write(6,*) 'Void ratio:'
        write(6,*) 
        write(6,109) 'e = ',y(7)
        write(6,*) 
        write(6,*) '==================================================='
  
101   format(1X,3(a9,e10.3,2X))
102   format(1X,a25,i1)
103   format(1X,a7,i5)
104   format(1X,3(a5,f10.4,2X))
105   format(1X,a5,i2,a4,f20.3)
106   format(1X,3(a9,E12.4,2X))
107   format(1X,3(a10,f12.4,2X))
108   format(1X,a8,e12.4)
109   format(1X,a4,f10.4)
110   format(1X,a5,f10.4)
111   format(1X,a6,i4,2X,a11,i4,2X,a9,i10,2X,a19,i4)
       
      return
      end
      
      
!-----------------------------------------------------------------------------
      subroutine warning(var1,var1v,nrelop,var2,var2v,set2,subr,prinfo)
!-----------------------------------------------------------------------------
      include 'vaba_param.inc'
      integer nrelop    
!      real*8 set2, var1v, var2v
      logical prinfo
      character(11) subr
      character(6) var1, var2
      
      if (prinfo) then 
      write(6,*) 'VUMAT: WARNING in subroutine ', subr
      if (nrelop==0) then
        write(6,*) 'cannot evaluate ', var1
        write(6,*) var1,'=',set2
      elseif (nrelop==1) then
        write(6,112) var1,'=',var1v,'is greater than ',var2,'=',var2v
      elseif (nrelop==2) then
        write(6,112) var1,'=',var1v,'is equal to     ',var2,'=',var2v
      elseif (nrelop==3) then
        write(6,112) var1,'=',var1v,'is less than    ',var2,'=',var2v
      endif
      endif
      
      if (nrelop==0) var1v = set2
      
112   format(1X,a6,a2,f9.6,a17,a6,a2,f9.6)      
      
      return
      end      
