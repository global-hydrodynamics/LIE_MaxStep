      program inrt_1d
! ================================================
! - to precisely estimate maximum acceptable timestep for stable simulation using Local Inertial Equation
! - Copyright (c) 2023 Dai Yamazaki & Tomohiro Tanaka
! - Distributed with MIT License
! - https://github.com/global-hydrodynamics/LIE_MaxStep/blob/master/LICENSE
! ================================================
      implicit none
!
      integer                ::  iseq, nseq              !! number of cells
      parameter                 (nseq=10)
!
      real*8                 ::  dx                      !! cell size
      real*8                 ::  dt0                     !! baseline dt for "true" simulation : dt0 = CFL*0.001
      real*8                 ::  dt,   dt_pre            !! dt used for simulation
      real*8                 ::  tnow                    !! time of current time step
      real*8                 ::  tend                    !! simulation end time (used to judge instability)

      real*8                 ::  cfl, froude
!
      real*8                 ::  elevtn(nseq+1)          !! river bed elevation
!
      real*8                 ::  rivdph(nseq+1)          !! river water depth
      real*8                 ::  sfcelv(nseq+1)          !! water surface elevation (sfcelv = elevtn + rivdph) 
      real*8                 ::  rivsto(nseq)            !! water storage (rivsto = rivdph * rivwth * dx)
      real*8                 ::  rivout(nseq)            !! river ourflow (to downstream cell)
      real*8                 ::  rivinf(nseq)            !! river inflow  (from upstream cell): rivinf(iseq+1)=rivout(iseq)

      real*8                 ::  outpre(nseq)            !! river flow        at previous time step
      real*8                 ::  dphpre(nseq+1)          !! river water depth at previous time step
!
!
      real*8                 ::  slope                   !! water surface slope
      real*8                 ::  fdph,  fare             !! flow deepth (Hflow), flow area (Hflow * width)
      real*8                 ::  fdph1, fare1, fdph2, fare2  !! used for Hflow scheme
      real*8                 ::  rivvel                  !! velocity

      real*8                 ::  change                  !! boundary conditon change rate
      parameter                 (change=0.0001)
      real*8                 ::  dphini                  !! Initial Depth
      real*8                 ::  dphend                  !! Simulation end depth:  dphend = dphini*(1+change)
                                                         !!  dphend is used as downstream boundary condition

      real*8                 ::  outini                  !! Initial Flow        :  outini = manval**(-1)*bedslp**0.5*dphini**(5/3)
      real*8                 ::  outend                  !! Simulation end flow :  outend = manval**(-1)*bedslp**0.5*dphend**(5/3)
      real*8                 ::  inflow                  !! Upstream inflow     :  inflow = outend


      real*8                 ::  bedslp                  !! river bed slope
      real*8                 ::  rivwth                  !! river width (set to 1: unit width)
      real*8                 ::  manval, gravity         !! Manning's n, gravity acceralation
      parameter                 (gravity=9.8)
      parameter                 (rivwth=1.0)
      !
      real*8                 ::  outerr, dpherr          !! Flow error, Depth error
      real*8                 ::  errthrs                 !! Error threshold to judge stability
      parameter                 (errthrs=0.001)

      integer                ::  istep
      integer                ::  iconv, nconv            !! # of time steps with error<thres, judged to be converged when iconv>nconv
      parameter                 (nconv=100)

      real*8                 ::  g, n, h, i, v, q, b, bi, gh  !! for easy result check (shorter name variables)

      real*8                 ::  ipow,  npow             !! parameters to increase dt
      real*8                 ::  npow1, npow2
      parameter                 (npow1=10)
      parameter                 (npow2=100)

      integer                ::  ifirst
! ##############################################
!!  model parameters

      dx     = 0.1
      dphini = 0.1
      bedslp = 1.0 / 1000.
      manval = 0.03

! ==============================================
      ifirst=1
 3000 continue             !!  restart simulation with modified condition

!!  change condition, and re-start simulation
      if( ifirst==0 )then

        dx=dx * 10.**(1./5.)
        if( dx>=26000. ) goto 7777    !!  simulaion end

!        bedslp=bedslp*(10.**0.05)
!        if( bedslp>0.2 ) goto 7777

!        manval=manval*(10.**0.2)
!        if( manval>0.5 ) goto 7777

!        dphini=dphini*(10.**0.1)
!        if( dphini>10.0 ) goto 7777
      endif
! ##############################################

!! calculate representative values (from initial condition)
      i = bedslp
      h = dphini
      g = gravity
      n = manval
      v = i**0.5 * n**(-1) * h**(2./3.)
      q = v * h

      b = n**2. / h**(10./3.)
      bi= (b*i)**0.5
      gh= g * h

      cfl = dx / (gh)**0.5 
      froude= v / (g*h)**0.5

! ================================================
! Set Topography

      do iseq=1, nseq+1
        elevtn(iseq)= bedslp * dx*(nseq-iseq+1)
      end do

! ##############################################
! Execute "True" Simulation (DT=0.001*CFL)

! Explicit form does not follow CFL. Fix dt to small value
!     dt0 = cfl * 0.001

! Fixed dt starts from a certain dx for saving calculation time
! Tanaka T. 2016/8/25
      if(dx > 1000) then
          dt0 = 1000 / (gh)**0.5 * 0.001
      else
          dt0 = cfl * 0.001
      end if
      
      dt  = dt0
      tend=1.e30

! ===
! Initial & Boundary Condition
      call set_init_cond(rivdph,sfcelv,rivsto,rivout,rivinf)

! ===
! time step loop

      tnow=0
      istep=0
      iconv=0
      do while( tnow<tend )
        tnow=tnow+dt
        istep=istep+1
        call calc_rivout(dt,rivout,rivinf,rivsto,rivdph,sfcelv,outpre,dphpre)
        call calc_error(rivout,rivdph,rivsto)

        if( outerr<=0.0001 .and. dpherr<=0.0001 )then
          iconv=iconv+1
          if( iconv>=nconv ) exit                !!  true simulation judged to be converged (err<thres contenuously for nconv steps)
        elseif( outerr>1000 )then
          print *, 'Not Converged', dx, dt0      !!  true simulation does not converge [something is wrong!]
          stop
        else
          iconv=0
        endif
      end do

      tend=tnow - dt*(nconv-1)                   !! simulation end time

!!      print '(a,10f20.6)', 'Calc.End Time (sec)', tend, outerr, dpherr
! ##############################################
! "Check" Simulation

      npow=npow1
      ipow=npow

 1000 continue

      dt=dt0 * 10.**(ipow/npow)      !!  increase dt

! ===
! Initial Condition
      call set_init_cond(rivdph,sfcelv,rivsto,rivout,rivinf)
! ===
! time step loop
      tnow=0
      istep=0
      iconv=0
      do while( tnow<tend*10 )
        tnow=tnow+dt
        istep=istep+1
        call calc_rivout(dt,rivout,rivinf,rivsto,rivdph,sfcelv,outpre,dphpre)
        call calc_error(rivout,rivdph,rivsto)

        if( outerr>1.e10 )then   !!  not converged
          exit
        endif

        if( outerr<=errthrs .and. dpherr<=errthrs )then
          iconv=iconv+1
          if( iconv>=nconv ) exit     !!  simulation judged to be converged, increase dt and re-do simulation
        else
          iconv=0
        endif
      end do
! ====================
! judge converge

      if( outerr<=errthrs .and. dpherr<=errthrs )then  !! simulation converged, increase dt and re-do simulation
        tnow=tnow - dt*(nconv-1)
!!        print '(a,10f20.6)', 'Calc.End (sec)', dt, cfl, dt/cfl, tnow, tnow, outerr, dpherr

      else                                !! simulation not converged
        if( npow==npow1 ) then               !! re-execute simulation with smaller dt increase step (to estimate precise max dt)
          ipow=ipow-1
          npow=npow2                         !!  change npow1 -> npow2 (to estimate exact dt+)
          ipow=ipow/npow1*npow2
          dt_pre=dt0 * 10.**(ipow/npow)
          ipow=ipow+1
          goto 1000        !!  repeat simulation with smaller dt step
        else
          goto 4000        !!  max acceptable time step dt+ is estimated, end simulation
        endif
      endif

      dt_pre=dt
      ipow=ipow+1
      goto 1000      !!  repeat simulation with updated dt

! ==============================================
!! simulation end, write results
 4000 continue 

      if( ifirst==1 )then  !! write header for first time
        print '(20a20)', 'dx', 'slope', 'depth', 'manning', 'max_dt+', 'CFL=dx(gh)^(-0.5)', 'dt+/CFL', 'Froude', 'Discharge'
      endif

      dt=dt_pre
      print '(20f20.6)', dx, i, h, n, dt, cfl, dt/cfl, froude, q

      ifirst=0

      goto 3000    !! repeat simulation with updated simulation setting
! =======================
! all simulation end
 7777 continue


      CONTAINS
! ##############################################


      subroutine  set_init_cond(dph,sfc,sto,out,inf)
! ===================================================
      real*8                 ::  dph(nseq+1)
      real*8                 ::  sfc(nseq+1)
      real*8                 ::  sto(nseq)
      real*8                 ::  out(nseq)
      real*8                 ::  inf(nseq)

      do iseq=1, nseq
        dph(iseq)=dphini
        sfc(iseq)=elevtn(iseq)+dph(iseq)
        sto(iseq)=dph(iseq)   *rivwth   *dx
      end do

      dph(nseq+1)=dphini                                       !! downstream boundary
      sfc(nseq+1)=elevtn(nseq+1)+dphini

      do iseq=1, nseq
        slope=(sfc(iseq)-sfc(iseq+1))/dx
        fdph= (dph(iseq)+dph(iseq+1))*0.5
        fare=fdph*rivwth
        rivvel = manval**(-1.) * slope**(0.5) * fdph**(2./3.)     !!  first guess by diffusion wave
        out(iseq)=fare*rivvel                                     !!  (preivous time step value is needed by LIE)
        if( iseq<nseq ) inf(iseq+1)=out(iseq)
      end do

  !!  boundary condition

      dphend = dphini * (1.+change)                            !! set downstream boundary condition to depth end
      dph(nseq+1) = dphend
      sfc(nseq+1) = elevtn(nseq+1) + dph(nseq+1)

      outini=out(1)                                            !  discharge at initial condition 
      outend=manval**(-1.) * bedslp**(0.5) * dphend**(5./3.) * rivwth  !! calculate simulation end flow
      inflow=outend
      inf(1)=inflow                                            !! set upstream boundary to simulation end flow

      end subroutine set_init_cond






      subroutine calc_rivout(dt1,out,inf,sto,dph,sfc,out0,dph0)
! ===================================================
      real*8                 ::  dt1
      real*8                 ::  out(nseq)
      real*8                 ::  inf(nseq)
      real*8                 ::  sto(nseq)
      real*8                 ::  dph(nseq+1)
      real*8                 ::  sfc(nseq+1)

      real*8                 ::  out0(nseq)
      real*8                 ::  dph0(nseq+1)


      !! update storage
      do iseq=1, nseq
        dph0(iseq)=dph(iseq)     !! keep previous time step depth and flow
        out0(iseq)=out(iseq)

        sto(iseq)=sto(iseq)+inf(iseq)*dt1-out(iseq)*dt1   !!  update river storage St+1=St + (Inflow-Ouflow)*dt
        dph(iseq)=sto(iseq)/dx/rivwth                     !!  river depth
        sfc(iseq)=elevtn(iseq)+dph(iseq)                  !!  water surface elevation
      end do

      !! calculate discharge
      do iseq=1, nseq
        slope=(sfc(iseq)-sfc(iseq+1))/dx
        fdph= (dph(iseq)+dph(iseq+1))*0.5
        fare=fdph*rivwth

! Hflow scheme (if needed)
        fdph1= (dph(iseq)*dph(iseq+1))**0.5
        fare1=fdph1*rivwth

        fdph2= (dph(iseq)*dph0(iseq))**0.5
        fare2=fdph2*rivwth

! choose hflow scheme (if needed)
!        fdph=fdph1
!        fare=fdph1
!        fdph=fdph2
!        fare=fare2

! Tanaka T. 2015/8/19
! Semi-implicit formation
!       if( fdph>0 )then
!         out(iseq)=(out0(iseq)+dt1*gravity*fare*slope) &
!                     / (1+(dt1*gravity*manval**2.*abs(out0(iseq)))/(fdph**(4./3.)*fare))
!       else
!         out(iseq)=0.
!       endif
! Explicit formation
        if( fdph>0 )then
          out(iseq)= out0(iseq)+dt1*gravity*fare*slope &
                      - (dt1*gravity*manval**2.*abs(out0(iseq))*out0(iseq))/(fdph**(4./3.)*fare)
        else
          out(iseq)=0.
        endif
        
        rivvel=out(iseq)/fare
        if( iseq<nseq ) inf(iseq+1)=out(iseq)
      end do

      end subroutine calc_rivout


      subroutine calc_error(out,dph,sto)
! ===================================================
      real*8                 ::  out(nseq)
      real*8                 ::  dph(nseq)
      real*8                 ::  sto(nseq)

!! both outerr and dpherr will be 0, when simulation is converged

      outerr=0
      do iseq=1, nseq
        outerr=outerr + abs(out(iseq)-inflow)/inflow
      end do
      outerr=outerr/dble(nseq)                          !!  averaged relative flow error

      dpherr=0
      do iseq=1, nseq
        dpherr=dpherr + abs(dph(iseq)-dphend)/dphend
      end do
      dpherr=dpherr/dble(nseq)                           !! averaged  relative depth error

      do iseq=1, nseq
        if( out(iseq)*dt > sto(iseq) )then      !!  if flow*dt>storage (overflow condition), set error to 1.e20
          outerr=1.e20
          dpherr=1.e20
        endif
      end do


      end subroutine calc_error




! ##############################################
      end program inrt_1d

