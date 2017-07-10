!  THIS IS THE CPD-NLG MODEL
      program cpd
      IMPLICIT real*4 (A-H,O-Z)
      save
!  Input parameters are found in PERKIN
!  Output parameters are found in PERKOT
!  Output parameters for nitrogen release found in perkot3
!  Output parameters for light gas species distribution in perkot4

      character*80 perkin,perkot,perkot2,perkot3,perkot4
      character*1 ans
      dimension y(4),yp(4),ypp(4),ypred(4),tim(50),tem(50),gasf(9,10)
      dimension heat(50)
      dimension yygas(5),ffgas(5)
      real*4 l0,l,ma,kb,kc,kg,kp,ft(35),mt(35),mgas,mtot,mw1,mdel,kn,mb
      integer*2 lib
!  intar = .true. calculates tar molecular weight distribution in subroutine perkp.
      logical intar,idiff,imolw,inside,yyyy
      common/cinit/l0,c0,g0,ma,rba,finf,sig,siginv,nmax,pstar
      common/rate/ab,eb0,ebsig,ac,ec0,ecsig,ag,eg0,egsig,rg
      common/tarmw/u0,beta,kkk,umax,yield,press
      common/timer/tms,time
      common/nitmod/fst,fnca0,fnca,an,en0,ensig,kn
      common/lgmod/yf,xoc,yhc
!  idiff = .true. for calculations on a differential basis rather than printing out integrated values.  This option worked on early versions,but has not been tested on the latest versions.
      data nmax/10/,nt/1/,intar/.false./,idiff/.false./
      data ftar0,fgas0,fchar0/0.,0.,1./,ip/30/,zero/0.0/
      data ftar/0.0/,fchar/0.0/,fmet/0.0/,fvol/0.0/,fvolold/0.0/
      data lib/0/
      data y/4*0./,yp/4*0./,ypp/4*0./ypred/4*0./,tim/50*0./,tem/50*0./,ft/35*0./,mt/35*0./,u0/0./,kkk/1/,umax/0./,yield/0./
      data perkot2/'cpdheat.out2'/,perkot3/'cpdheat.out3'/,perkot4/'cpdheat.out4'/
      data perkin/'cpdheat.inp'/,perkot/'cpdheat.out1'/
!  y1 = l	labile bridges
!  y2 = del	ends
!  y3 = c	char links
!  y4 = fnca    mass fraction of nitrogen in site
      data g0,g1,g2/0.,0.,0./
!   read in input data for kinetics and output file names
!      print*,' Enter input file name:'
!      read 5,cpdheat
!      print*,' Enter output file name for tar and gas release:'
!     read 5,perkot
!     print*,' Enter output file name for CPD model parameters:'
!      read 5,perkot2
!      print*,' Enter output file name for nitrogen release:'
!      read 5,perkot3
!      print*,' Enter output file name for light gas distribution:'
!      read 5,perkot4
5     format(a)
      open(unit=1,name='cpdheat.inp', type='old')
      open(unit=2,name='perkot1.txt', type='unknown')
      open(unit=20,name='perkot2.txt', type='unknown')
      open(unit=21,name='perkot3.txt', type='unknown')
      open(unit=22,name='perkot4.txt', type='unknown')
!  read in measured p0 from NMR data
      read(1,*)p0
!  estimate c0
      read(1,*)c0
!  read in sigma+1 (coordination number) from NMR data
      read(1,*)sigp1
!  read in real MW/Cluster from NMR data (includes side chains)
      read(1,*)mw1
!  read in real mdel from NMR data
      read(1,*)mdel
!  read in fnit (the percent nitrogen in parent coal)
      read(1,*)fnit
!  read in fst (the percent stable nitrogen in parent coal)
      read(1,*)fst
!  read in mass fractions of hydrogen, carbon, and oxygen
      print*,'p0,fst',p0,fst
      read(1,*)fhyd
      read(1,*)fcar
      read(1,*)foxy		
!  kinetic parameters
      read(1,*)ab
      read(1,*)eb0
      read(1,*)ebsig
      read(1,*)ac
      read(1,*)ec0
      read(1,*)ag
      read(1,*)eg0
      read(1,*)egsig
      read(1,*)acr
      read(1,*)ecr
!  kinetic parameters for nitrogen released from char
      read(1,*)an
      read(1,*)en0
      read(1,*)ensig
!  pressure in atmospheres
      read(1,*)press
      print*,'an,press', an,press
!  empirical correlation to allow a small portion of alpha-carbon to stay with the aromatic cluster
      mdel = mdel-7
!  now calculate other chemical structure coefficients
      l0 = p0 - c0
      mb = 2.*mdel
      ma = mw1-sigp1*mdel
      beta = ma
      sig = sigp1-1
      finf = 1./(2.*ma/(mb*sigp1*(1-c0))+1)
      rba = mb/ma
      print*,'rba=',rba
      fnca0 = fnit*mw1/(mw1-sigp1*mdel)
!   calculate O/C and H/C ratios for light gas model
      xoc = (foxy/16)/(fcar/12)
      yhc = fhyd/(fcar/12)
! test to see if p0 is above percolation threshhold
      pthresh = 1/sig
      if(p0.le.pthresh)then
        print*,'Error: value of p0 is below the percolation threshold'
	stop
      endif
!   read in input data for temperature profile  this version reads in heating rates and temperatures the time read in is valid only if there is a hold time.
!   This time is read in in seconds. 
      tim(1) = 0.0
      read(1,*)tem(1)
      read(1,*)ntim
      ntim=ntim+1
!      print*,tim(1),tem(1)
      do 10 i=2,ntim
        read(1,*)heat(i),tem(i),tim(i)
	if(heat(i).ne.0.)then
	  tim(i) = (tem(i)-tem(i-1))/heat(i)+tim(i-1)
	endif
!	write(2,*)heat(i),tem(i),tim(i)
10    continue
!   read time step for calculation
      read(1,*)dt0,iprint,dtmax
      read(1,*)timax
      read(1,*)nmax
      read(1,111)inside
111   format(L5)
!   initialize variables
      y(1) = l0
      y(2) = 2.*(1.-c0-l0)
      y(3) = c0
      y(4) = fnca0
      aind0 = l0 + (1.-c0-l0)
      siginv = 1./sig
      pstar = 0.5*siginv
      ynhcn = 0.0
      yntar = 0.0
      yytar = 0.0
      yf = 0.0
      yyyy=.false.
!   START calculation DO LOOP
      Rg = 1.987                       !CAL/GMOLE K
      write(20,203)
203   format('c time(ms)      l         c      del/2    g1/2    g2/2      gtot      p')
      write(21,204)
204   format('c time(ms)   temp(K)    mw        Nca      Nclus      fnchar     fntar     fnhcn     fntot')
      write(22,205)
205   format('c time(ms)      fgas     fh20     fco2     fch4        fco     fother    yh20     yco2     ych4      yco     yother    Xgas')
      iii = 0
      fgas0 = 0.
      fcross = 0.0
      fntar = 0.0
      fntd = 0.0
      fnt = 0.0
      cmw1 = 0.0
      smw1 = 0.0
      if(heat(2).ne.0.)dt0 = tim(2)/1000
      dt = min(dt0,dtmax)
      ntmax = 100*(timax/dt+1)
      DO 100 iii=1,ntmax
        time = time+dt
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C PREDICTOR
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!   CALCULATE particle temperature
1       IF(time.LE.tim(nt+1))then
          Tp = tem(nt) + (time-tim(nt))*(tem(nt+1)-tem(nt))/(tim(nt+1)-tim(nt))
        elseif(nt.lt.ntim)then
          nt = nt+1
	  if(heat(nt+1).ne.0.0)then
	    dt = (tem(nt+1)-tem(nt))/heat(nt+1)/10
!	    print*,'nt=',nt,'  tem(nt)=',tem(nt),'  new dt=',dt
!	    print*,'tem(nt+1)=',tem(nt+1),'  heat(nt+1)=',heat(nt+1)
	  endif
          go to 1
        ELSE
          PRINT*,'REACHED END OF GAS TEMPERATURE CORRELATION'
	   if(inside.eqv..false. .and. yyyy.eqv. .true.)then
	     print*,'!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!'
	     print*,'O/C and H/C ratios are outside the  '
	     print*,'bounds of the library coals. '
	     print*,'Estimation of Light gas distribution '
	     print*,'is based on library coal no.',lib
	   endif
          go to 400
        ENDIF
        if(tp.gt.4000.)then
           print*,'     >>>>   WARNING    <<<<<'
           print*,' gas temperature too high----',tp
        endif
!--    DEVOLATILIZATION RATES 
49      call perks(y,ypp,Tp)
! COMPONENT MASS CONSERVATION
        do 50 j=1,4
           ypred(j) = y(j) + dt*ypp(j)
           ypred(j) = max(ypred(j),zero)
50      continue
!        print*,' ypred and ypp=',ypred,ypp
!	 print*,' ypred(4) = ',ypred(4)
!ccccccccccccccccccccccccccccccccccccccccccccccc
!  CORRECTOR
!ccccccccccccccccccccccccccccccccccccccccccccccc
!--    DEVOLATILIZATION RATES 
        Call perks(ypred,yp,Tp)
! COMPONENT MASS CONSERVATION
!   adjust time scale
	      dtold = dt
	      dymax = 0.
	      dymin = 1.
           do 60 j=1,4
	     if(y(j).gt.5.e-3)then
                dy1 = abs(dt*0.5*(yp(j)+ypp(j)))
	        dymax=max(dymax,dy1)
             endif
60	   continue
	   dm = fvol-fvolold
	   bigdm = max(dymax,dm)
           if(bigdm.lt.0.0001)then
              itest = 1
	      dt = dt*1.2
!	      print*,'at time=',time'  tried to increase dt from',dtold,'     to ',dt
!              if(dt.lt.dtmax)print*,'at time=',time,'dt increased from',dtold,'     to ',dt
           elseif(bigdm.gt.0.0005)then
              itest = 2
	      dt = dtold/2
!              dt = 0.0005*dtold/bigdm
!              print*,'at time=',time, 'dt decreased from'dtold,' to ',dt
	      print*,'dy1=',dy1
           endif
           dt = min(dt,dtmax)
	   if((dtold.ne.dt).and.(itest.lt.1))go to 49
	   fvolold = fvol
!	   print*,'dt just before 70 do loop =',dt
           dtn = dt
! COMPONENT MASS CONSERVATION
!  y1 = bridges
!  y2 = side chains
!  y3 = cross links (mass fraction basis)
!  y4 = mass of nitrogen per aromatic site (N_site; see Perry dissertation, 1999
        do 70 j=1,4
           y(j) = y(j) + dt*0.5*(yp(j)+ypp(j))
           y(j) = max(zero,y(j))
70      continue
!        if(y(1).lt.1.e-5)dt = dt0*10.
        fracr = 1.
        if(fmet.gt.1.e-5.and.acr.gt.0.)then
!           print*,'<main>acr,ecr,rg,tp=',acr,ecr,rg,tp
           ratecr = acr*exp(-ecr/rg/tp)*fmet*dt
           fracr = 1.-ratecr/fmet
           fmet = fmet-ratecr
           fcross = fcross+ratecr
           if(fmet.lt.0.)then
             fcross = fcross-ratecr+fmet+ratecr
             fmet = 0.0
             fracr = 0.
           endif
         endif
!  print out data
        if(y(1).gt.1.e-5)intar = .true.
        call perkp(y,mgas,ftar,ftart,fgas,fchart,ft,mt,intar)
        intar = .false.
      tms = time*1.e3
      if(idiff)then
          fdift = ftar-ftar0
          ftar0 = ftar
          fdifg = fgas-fgas0
          fdifc = -(fchar-fchar0)
          fchar0 = fchar
        write(2,2010)tms,tp,fcross,y(1),fdift,fdifg,fdifc
        write(6,200)tms,tp,fcross,y(1),fdift,fdifg,fdifc
      else
          gasmw = rba*ma/2.
          if(fgas.ge.1.e-5)then
           fgasd = fgas-fgas0
           fgas0 = fgas
           imolw = .false.
           if(mod(iii,iprint).eq.0)imolw=.true.
           call flash(fgas,gasmw,ft,mt,fracr,ftar,fmet,tp,press,nmax,imolw,sumy)
          elseif(fgas.lt.1.e-5)then
           fmet = ftart
           ftar = 0.
          endif
          intar = .false.
          fvol = fgas+ftar
          fchar = 1.-fvol
          if(fvol.gt.0.995)stop
      if(iii.eq.1)then
      write(6,201)
      write(2,202)
201   format(' time(ms)      temp     fcross     labile    ftar    fgas   fsolid     ftot    fmet')
202   format('c time(ms)      temp     fcross     labile   ftar     fgas   fsolid   ftot   fmet')
      endif
        if(mod(iii-1,iprint).eq.0)then
          write(2,200)tms,tp,fcross,y(1),ftar,fgas,fchar,fvol,fmet
!          write(6,*)fcross,fmet
          write(6,200)tms,tp,fcross,y(1),ftar,fgas,fchar,fvol,fmet
        endif
      endif
200   format(' ',2(1pe10.3,2x),0pf6.4,2x,f8.5,2x,5(f6.4,2x))
2010  format(' ',2(1pe10.3,2x),0pf6.4,2x,f8.5,2x,4(1pe10.3,2x))
      l = y(1)
      del = y(2)
      c = y(3)
      fnca = y(4)
      p = l+c
      tarfac = 1.
      g1 = (2.*(1-p)-del)
      g2 = 2.*(c-c0)
      g = g1+g2
      del2 = del/2.
      aind = del2 + l
      g12 = g1/2.
      g22 = g2/2.
      gtot = g/2.
    
!  NITROGEN RELEASE CALCULATIONS

!  new average molecular weight
      cmw1 = ma + (c0 + l + del2)*sigp1*mdel
!  nitrogen content of char and tar 
      fnt = y(4)*ma/cmw1
!  differential mass of nitrogen released with tar
      yntd = fnt*sumy
!  total nitrogen released with tar
      yntar = yntar + yntd
!  nitrogen remaining in char
      ynchar = fchar*fnt
!  differential mass of nitrogen released as light gas     
      ynhcnd = fchar*dt*(-0.5)*(yp(4)+ypp(4))*ma/cmw1
!  total mass of nitrogen released as light gas 
      ynhcn = ynhcn + ynhcnd
!  mass of nitrogen released with tar by difference
      yntar = fnit - ynhcn - ynchar
!  fraction of original nitrogen remaining in char
      fnchar = ynchar/fnit
!  fraction of original nitrogen released as tar
      fntar = yntar/fnit
!  fraction of original nitrogen released as light gas
      fnhcn = ynhcn/fnit
!  total fractional release of nitrogen
      fntot = (fnit - fnt*fchar)/fnit

    
!  DISTRIBUTE LIGHT GAS INTO H20, CO2, CO, CH4, & other HC's yf is a CPD indicator of the fraction of total light gas that has been released. The look up table on light gas composition is based on yf. 
      if(inside)then
         yf = 1-aind/aind0
       
         call lightgas(yygas,inside,lib,yyyy)
      
!  calculate fraction of total mass release that is h2o, co2, ch4, co, and other light gases

          do 218 ik=1,5
            ffgas(ik)=fgas*yygas(ik)
218       continue 
      endif
      
!  PRINT RESULTS IN OUTPUT FILES

      if(mod(iii-1,iprint).eq.0)then
        write(20,220)tms,l,c,del2,g12,g22,gtot,p
        write(21,221)tms,tp,cmw1,y(4),fnt,fnchar,fntar,fnhcn,fntot
        if(inside)write(22,222)tms,fgas,ffgas(1),ffgas(2),ffgas(3),ffgas(4),ffgas(5),yygas(1),yygas(2),yygas(3),yygas(4),yygas(5),yf
      endif
220   format(' ',1pe10.3,2x,7(0pf7.5,2x),0pf12.9,0pf9.7,0pf10.5,0pf12.9)
221   format(' ',2(1pe10.3),0pf8.2,6(0pf10.6))
222   format(' ',1pe10.3,0pf9.4,11(0pf9.4))
100   continue
400   continue
          intar = .false.
          fvol = fgas+ftar
          fchar = 1.-fvol
          if(fvol.gt.0.995)stop
          write(6,201)
          write(2,200)tms,tp,fcross,y(1),ftar,fgas,fchar,fvol,fmet
          write(6,200)tms,tp,fcross,y(1),ftar,fgas,fchar,fvol,fmet
      stop
      end

!c------------------------------------------------------------------

      subroutine perks(y,yp,T)
      IMPLICIT real*4 (A-H,O-Z)
      save
!  this subroutine is the meat of the devolatilization model      
!  y1 = l	labile bridges
!  y2 = del	ends
!  y3 = c	char links
!  y4 = fnca    nitrogen release as HCN
!  yp(i) = derivative of y(i) in time
!     nmax = number of terms in expansion for mol. wt. distribution
      real*4 l0,l,ma,kb,kc,kg,kn,kp,ft(35),mt(35),y(4),yp(4)
      common/cinit/l0,c0,g0,ma,rba,finf,sig,siginv,nmax,pstar
      common/rate/ab,eb0,ebsig,ac,ec0,ecsig,ag,eg0,egsig,rg
      common/nitmod/fst,fnca0,fnca,an,en0,ensig,kn
      data fx/0./
      l = y(1)
      del = y(2)
      c = y(3)
      fnca = y(4)
      p = l+c
      g1 = 2.*(1-p)-del
      g2 = 2.*(c-c0)
      g = g1+g2
!  calculate current activation energy using error function solution
      if(c0.lt.1.)fx = g/(1.-c0)/2.
      call inverf(fx,x)
      eg = eg0 + x*egsig
      if (fnca0.gt.0)fx = 1.-(fnca-fst*fnca0)/(fst*fnca0)
      call inverf(fx,x)
      en = en0 + x*ensig
      if(l0.gt.0.)fx = 1.-l/l0
      call inverf(fx,x)
!      print*,'temperature=',t,'   eg=',eg,'  g/ginf=',fx
      eb = eb0+x*ebsig
      ec = ec0
!  calculate rate constants
      rt = rg*t
      kb = ab*exp(-eb/rt)
      rho = ac*exp(-ec/rt)
      kg = ag*exp(-eg/rt)
      kn = an*exp(-en/rt)
!  calculate rate of destruction of labile bridges
      yp(1) = -kb*l
!  calculate rate of formation of ends (danglers)
      yp(2) = 2.*rho*kb*l/(rho+1.) - kg*del
!  calculate rate of formation of char
      yp(3) = kb*l/(rho+1.)
!  calculate rate of nitrogen loss (as HCN)
      yp(4) = -kn*max((fnca-fst*fnca0),0.)
      return
      end
!------------------------------------------------------------------
      subroutine perkp(y,mgas,ftar,ftart,fgas,fchar,ft,mt,intar)
      IMPLICIT real*4 (A-H,O-Z)
      save
!   calculates fractions of tar, gas, and char from p, sig, l, and c
      real*4 l0,l,ma,mb,kb,kc,kg,kp,y(4),yp(4),ft(35),mt(35),mgas,mtot,sumu(200)
      logical intar
      common/cinit/l0,c0,g0,ma,rba,finf,sig,siginv,nmax,pstar
      common/rate/ab,eb0,ebsig,ac,ec0,ecsig,ag,eg0,egsig,rg
      common/tarmw/u0,beta,kkk,umax,yield,press
      data pi/3.141592/,iprint/0/

      l = y(1)
      del = y(2)
      c = y(3)
      p = l+c
!     print*,' <perk>l,del,c,g1,g2,g,p,mgas,mtot=',l,del,c,g1,g2,g,p,mgas,mtot
      if(intar)then
        if(p.gt.0.9999)then
          delfac = 1.
        else
          delfac = del/(1.-p)
        endif
        a = 1.+rba*(l/p + (sig-1.)/4. *delfac)
        b = (delfac/2. - l/p)
!  find pstar
        pstar0 = pstar
        pinv = siginv+1.e-4
        if(p.ge.0.9999)then
           pstar = 0.
        elseif(p.ge.pinv)then
           do 10 i=1,25
              f = pstar*(1-pstar)**(sig-1) - p*(1-p)**(sig-1)
              fp = (1-pstar)**(sig-1)-pstar*(sig-1)*(1-pstar)**(sig-2)
              ppstar = pstar - f/fp
              err= abs(1.-ppstar/pstar)
              if(err.le.1.e-4)go to 11
              pstar = ppstar
10         continue
           print*,' warning--pstar did not converge'
           print*,' p=',p,' sig=',sig,' pstar=',pstar,'  pstar0=',pstar0
11         continue
        else
           pstar = p
        endif
!  check to see if pstar is in the right range
      if(pstar.lt.0. .or. (p.ne.pstar .and. pstar.ge.siginv))then
           print*,' error--pstar out of acceptable ranges!'
           print*,'     pstar=',pstar,';  HOW DID YOU DO THAT?'
           print*,'p=',p,'  sig=',sig,'  pstar0=',pstar0
           stop
        endif
        sfac = (sig+1.)/(sig-1.)
        fp = (pstar/p)**sfac
        kp = (pstar/p)**sfac*(1.-(sig+1)/2. *pstar)
!  calculate wt fraction tar, gas, and char
        ftart = 2.*(a*fp+rba*b*kp)/(2.+rba*(1.-c0)*(sig+1.))
      endif
      tarfac = 1.-ftar
      g1 = (2.*(1.0 - p) - del)
      g2 = 2.*(c-c0)
      g = g1+g2
      mgas = rba*ma*g*(sig+1)/4.*tarfac
      mtot = ma + rba*ma*(sig+1)/2. *(1.-c0)
      fgas = mgas/mtot
      fchar = 1.-ftar-fgas
!      print*,'<perkp>ftar,fgas,fchar=',ftar,fgas,fchar
!  calculate tar molecular weight distribution
      if(.not.intar)return
      iprint = iprint+1
!      if(iprint.eq.1) open(unit=25,name='perktar.out',type='unknown')
      ftsum = 0.0
      do 20 n=1,nmax
         tn = n*(sig-1.)+2
         xm = n*sig+1.
         yk = n-1.
         xm1 = xm+1.
!  gamln is the solution to the gamma function in the Sandia Math Library
         fg1 = gamln(xm1)
         if(fg1.le.1.e-10) then
            fgam = 0.
         else
           yk1 = yk+1.
           fg2 = gamln(yk1)
           xmyk = xm-yk+1.
           fg3 = gamln(xmyk)
           fgam = exp(fg1-fg2-fg3)
         endif
         bnn = (sig+1.)/(n*sig+1.)*fgam
         qn = bnn*(p**(n-1))*((1-p)**tn)/n
!  ft(n) = weight fraction of each tar bin
         ft(n) = 2.*(n*a*qn+rba*b*qn)/(2.+rba*(1.-c0)*(sig+1.))
         ftsum = ftsum + ft(n)
!  check to not divide by zero
         if(p.le.1.e-9)then
            fac = 0
         else
            fac = l/p
         endif
         tst = 1.-p
         if(tst.le.1.e-9)then
            fac1 = 0.
         else
            fac1 = del/(1.-p)
         endif
!  mt(n) = molecular weight of each tar bin
         mt(n) = n*ma+(n-1)*rba*ma*fac+tn*rba*ma/4.*fac1
!         if(iprint.eq.1)then
!           if(n.eq.1)then
!             print*,'Initial Fragment size distribution'
!             print*,'n   fraction      Cumul.   MW'
!           Endif
!           print*,n,ft(n),ftsum,mt(n)
!           write(24,*)ft(n),ftsum,mt(n)
!         endif
20    continue
      return
      end
      SUBROUTINE INVERF(Y,X)
      IMPLICIT real*4 (A-H,O-Z)
      save
!  this program calculates the inverse of the area under the normal curve.
!  if y=area(x), then given y, this program will calculate x.
!  A table lookup is performed.
      dimension xx(18),yy(18)
      data xx/3.4,3.2,3.,2.8,2.6,2.4,2.2,2.,1.8,1.6,1.4,1.2,1.,.8,.6,.4,.2,0./
      data yy/.9997,.9993,.9987,.9974,.9953,.9918,.9861,.9772,.9641,.9452,.9192,.8849,.8413,.7881,.7257,.6554,.5793,.5/
      fac = 1.
!  check to see if y is within range
      if(y.lt.0.0228)then
         x = -2.0
         return
      elseif(y.lt.0.5)then
        yp = 1.-y
        fac = -1.
      elseif(y.gt.0.9997)then
        x = 3.5
        return
      else
        yp = y
      endif
!  search for range
      do 10 i=17,1,-1
        if(yp.le.yy(i+1))then
          x = xx(i) + (yp-yy(i))*(xx(i+1)-xx(i))/(yy(i+1)-yy(i))
          x = fac*x
          return
        endif
10    continue
      return
      end
      function gamln(x)
!   this is a program to calculate the ln of the gamma function,taken from Abramowitz, p. 257, 6.1.41
      data pi/3.141592/
         gamln = (x-.5)*alog(x)-x+.5*alog(2.*pi)+1./(12.*x)-1./(360.*x**3)+1./(1260.*x**5)-1./(1680.*x**7)
      return
      end

      subroutine flash(fgas,gasmw,ft,mt,fracr,ftar,fmet,temp,press,nmax,imolw,sumy)
      IMPLICIT real*4 (A-H,O-Z)
      save
!  flash distillation of metaplast to form liquid and tar vapor
      logical imolw
      real*4 k(36),l(36),Ltot,mt(35),metold(35)
      dimension f(36),ft(35),xmw(36),z(36),pv(36),x(36),y(36),v(36),ftold(35),tarold(35)
      common/timer/tms,time
      data small/1.e-3/,a,b,g/87058,299,0.5903/,x3,x2/.2,.3/,zero/0.0/,ip/30/,x/36*0./,y/36*0./
      data k/36*0/,l/36*0/,v/36*0/,metold/35*0./,ftold/35*0./,xmw/36*0.0/,tarold/35*0./,f/36*0./,z/36*0./,pv/36*0./,fgas0/0./

!  metold(i) = mass fraction of coal contained in metaplast of mer size i,  fracr = fraction to account for reduction of metaplast by crosslinking in latest time step
!  renormalize in tar fragment data to molar basis
!  f(i) = moles of mer i
!  ft(i) = mass of mer i
  
      Ftot = 0.0
      do 10 i=1,nmax
         i1=i+1
         xmw(i1) = mt(i)
         dif = ft(i)-ftold(i)
         dif = max(dif,zero)
         f(i1) = (dif+metold(i)*fracr)/mt(i)
         ftold(i) = ft(i)
         Ftot = Ftot + f(i1)
         if(f(i1).gt.1.0)print*,'f(i1),ft(i),ftold(i),metold(i)=',i1,f(i1),ft(i),ftold(i),metold(i)
10    continue
      ntot = nmax + 1
      f(1) = (fgas-fgas0)/gasmw
      f(1) = max(f(1),0.)
      fgas0 = max(fgas,fgas0)
      xmw(1) = gasmw
      Ftot = Ftot + f(1)
!  get mole fraction of components in the feed stream and compute equilibrium contants k(i) from vapor pressure and Raoults law expression
      sum = 0.0
      do 30 ii=1,ntot
        sum = sum + f(ii)
        pv(ii) = a*exp(-b*xmw(ii)**g/temp)
        k(ii) = pv(ii)/press
        if(k(ii).lt.0.001)k(ii) = 0.0
30    continue
      if(sum.le.1.e-8)return
      do 40 ii=1,ntot
        z(ii) = f(ii)/sum
40    continue
!  use the Rachford-Rice formulation for flash distillation
!  x = V/F, first guess
      x1 = x3
!  calculate sum (Eq. 11-24, Separation Processes, by King, 1971)
      f1 = 0.0
      do 50 ii=1,ntot
         f1 = f1 + z(ii)*(k(ii)-1)/((k(ii)-1)*x1+1)
50    continue
!  Secant Method for convergence
      test = x2-x1
      if(test.lt.0.005)x2=x1+0.005 
      do 70 iter=1,100
!  calculate sum (Eq. 11-24, Separation Processes, by King, 1971)
        f2 = 0.0
        do 60 ii=1,ntot
           f2 = f2 + z(ii)*(k(ii)-1)/((k(ii)-1)*x2+1)
60      continue
        if(abs(f2).le.small.or.abs(f2-f1).le.small**2)go to 100
        x3 = x2-f2*(x2-x1)/(f2-f1)
        if(x3.gt.1.0)x3 = 1.0-small**2
        if(x3.lt.0.)x3 = small**2
        if(x3.eq.x2)then
          print*,'problem---f(V/F) not converged, but x3=x2',x3,x2
          if(x2.ge.small)then
            x3=x2-small
          else
            x3 = x2+small
          endif
        endif
        if(x2.le.1.e-5.and.x1.le.1.e-5)then
           x2 = 1.e-7
           go to 100
        endif
        if(x2.ge..9999 .and. x1.ge..9999)then
           x2 = .9999
           go to 100
        endif
        f1 = f2
        x1 = x2
!  under-relax solution (especially near the V/F=1 limit
        x2 = 0.2*x2+0.8*x3
!        print*,x1,x2,x3
70      continue
        print*,'Convergence not achieved in Flash distillation'
        print*,' last two guesses for V/F were',x1,x3
        do 81 ii=1,ntot
          write(30,*)z(ii),k(ii)
81      continue
        stop
100     continue
!  now calculate molecular weight distributions on a light-gas free basis, wt fractions
        Vtot = Ftot*x2
        Ltot = Ftot-Vtot
        VoL = Vtot/Ltot
        sumx = 0.0
        sumy = 0.0
        xmwtot = 0.0
        ttot = 0.0
        do 200 ii=2,ntot
          i = ii-1
          l(ii) = f(ii)/(1.+k(ii)*VoL)
          v(ii) = f(ii)-l(ii)
          x(ii) = l(ii)*xmw(ii)
          y(ii) = v(ii)*xmw(ii)
          metold(i) = max(x(ii),zero)
          tarold(i) = tarold(i)+y(ii)
          xmwtot = xmwtot+tarold(i)*xmw(ii)
          ttot = ttot+tarold(i)
          sumx = sumx + x(ii)
          sumy = sumy + y(ii)
200     continue
        if(ttot.gt.0.)xmwtot = xmwtot/ttot
        do 250 ii=2,ntot
           if(sumx.gt.1.e-28)x(ii) = x(ii)/sumx
           if(sumy.gt.1.e-28)y(ii) = y(ii)/sumy
250     continue
        ftar = ftar+sumy
        fmet = sumx
        if(imolw)then
          ip = ip+1
!          print*,'Weight Av. Molecular Wt.=',xmwtot ,' Gas MW =',gasmw
        endif
210     format(' ',i2,2x,f5.0,4(1pe10.3,2x))
211     format(' ',f5.0,4(1pe10.3,2x))
        return
        end
!  ------------------------------------------------------------------
!  --------------LIGHT GAS DISTRIBUTION MODEL------------------------
!  ------------------------------------------------------------------
	subroutine lightgas(yygas,inside,lib,yyyy)
!  This program calculates the distribution of light gas species  based on a look up table of reference data

	IMPLICIT real*4 (A-H,O-Z)
	real*4 ind
	integer*2 lib
	logical inside,yyyy
	common/lgmod/yf,xoc,yhc
	save
	dimension xx(12),yy(12),a(23),b(23)
        integer*2 s1(12),s2(12),s3(12)
	integer*2 p1(12),p2(12),p3(12),n(12)
	dimension xz(12,12),fh2o(12,12),fco2(12,12)
	dimension out(5,14),fco(12,12),fch4(12,12),ygas(4,3),yygas(5)
	
! *******************************************************************
! ****************LIGHT GAS DISTRIBUTION REFERENCE LIBRARY***********
! *******************************************************************

!  This library can be moved outside of submodel as long as  it is linked to the light gas sub-model
!  xz = The index used in correlation.  In main program it  corresponds to variable "yf"
!  f*** = fraction of light gas that is species ***
!  the data is organized in the following order with 12 ordered pairs for each species (ordered with xz)
!  each table is organized in rows in the following order 
! 	1 Lower Kittaning (Chen)
!	2 Pocahontas #3 (ANL)
!	3 Upper Freeport (ANL)
!	4 Pittsburgh (Chen)
!	5 Lewis Stockton (ANL)	
!	6 Utah Blind Canyon (ANL)
!	7 Illinois #6 (ANL)
!	8 Illinois #6 (Chen)
!	9 Wyodak (ANL)
!	10 Beulah Zap (ANL)
!	11 Dietz (Chen)
!	12 PSOC 1448 (Serio)

!  reference data for xz = yf, the fractional light gas released

        data xz/0.,.04,.11,.14,.21,.27,.34,.675,.9,1.,.0,.0,.0,.161,.442,.663,.777,.874,.921,.967,1.,.0,.0,.0,.0,.022,.20,.430,.526,.64,.787,.875,.927,.955,1.,.0,.0,.04,.12,.15,.23,.29,.36,.68,.9,1.,.0,.0,.0,.018,.058,.21,.417,.572,.696,.778,.821,.883,.932,1.,.0,.052,.144,.291,.498,.639,.746,.859,.925,.949,.966,1.,.0,.063,.178,.33,.506,.612,.706,.813,.895,.94,1.,.0,.0,.04,.12,.15,.23,.29,.36,.68,.9,1.,.0,.0,.0,.061,.146,.374,.535,.622,.714,.8,.883,.931,.964,1.,.0,.034,.087,.179,.316,.472,.585,.694,.777,.872,.935,1.,.0,.04,.12,.16,.25,.31,.37,.68,.9,1.,.0,.0,.0,.02,.055,.17,.313,.434,.546,.716,.874,.935,.973,1./
     
!  fraction of light gas that is water
     
	  data fh2o/.772,.772,.738,.455,.371,.304,.290,.273,.218,.218,.0,.0,.699,.632,.299,.269,.247,.249,.236,.225,.226,.0,.0,.0,.0,.0,.35,.297,.301,.299,.284,.291,.306,.297,.283,.0,.636,.636,.646,.550,.436,.320,.186,.199,.195,.195,.0,.0,1.,.983,.754,.488,.413,.385,.373,.382,.377,0.362,.367,.348,.665,.636,.604,.508,.435,.409,.383,.362,.351,.343,.342,.339,.763,.737,.698,.572,.527,.470,.438,.411,.411,.396,.378,.0,.748,.748,.637,.704,.490,.446,.348,.268,.266,.266,.0,.0,.0,.0,.385,.461,.396,.369,.344,.323,.292,.277,.266,.257,.0,.0,.197,.267,.26,.333,.361,.369,.346,.306,.285,.267,.521,.521,.55,.523,.511,.46,.414,.388,.313,0.313,.0,.0,.0,.0,.291,.335,.264,.271,.261,.211,.171,.160,.153,.149/
     
!  fraction of light gas that is carbon dioxide

	  data fco2/.0,.0,.0,.174,.174,.167,.129,.102,.071,.071,.0,.0,.259,.234,.113,.086,.097,.109,.116,.118,.122,.0,.0,.0,.333,.327,.070,.052,.057,.06,.059,.062,.066,.08,0.115,.0,.194,.194,.152,.117,.116,.122,.081,.092,.065,.065,.0,.0,.0,.0,.0,.122,.103,.086,.083,.082,.085,.086,.093,.128,.332,.318,.165,.141,.120,.108,.105,.119,.120,.122,.125,.130,.229,.221,.125,.09,.07,.073,.083,.133,.132,.13,.147,.0,.111,.111,.142,.175,.149,.155,.136,.122,.133,.133,.0,.0,.98,.984,.55,.345,.317,.285,.286,.277,.273,.264,.254,.255,.993,.989,.786,.572,.519,.416,.375,.345,.335,.32,.303,.299,.363,.363,.353,.325,.321,.35,.318,.251,.249,.249,.0,.0,1.,.983,.448,.179,.104,.09,.104,.151,.166,.160,.158,.154/
     
!  fraction of light gas that is methane

	  data fch4/.203,.203,.078,.160,.180,.219,.258,.294,.320,.320,.0,.0,.041,.037,.388,.389,.359,.332,.323,.307,.299,.0,.0,.0,.667,.655,.42,.454,.444,.419,.382,.353,.331,.321,.306,.0,.055,.055,.073,.088,.116,.124,.170,.15,.189,.189,.0,.0,.0,.0,.188,.195,.234,.243,.224,.21,.2,.186,.177,.167,.0,.0,.11,.155,.176,.172,.185,.173,.163,.159,.156,.151,.0,.0,.075,.136,.159,.178,.174,.157,.143,.141,.132,.0,.02,.02,.026,.042,.045,.049,.064,.1,.128,.128,.0,.0,.0,.0,.0,.029,.048,.067,.069,.072,.069,.066,.063,.061,.0,.0,.0,.0,.035,.05,.061,0.058,.057,.053,.049,.046,.01,.01,.011,.016,.011,.021,.023,.035,.06,.06,.0,.0,.0,.0,.216,.262,.362,.327,.307,.25,.203,.189,.182,.177/
     
!  fraction of light gas that is carbon monoxide
     
	  data fco/.0,.0,.157,.121,.141,.112,.139,.085,.145,.145,.0,.0,.0,.0,.0,.057,.097,.109,.124,.15,.153,.0,.0,.0,.0,.0,.0,.0,.0,.024,.078,.097,.099,.104,.099,.0,.083,.083,.038,.066,.032,.168,.286,.324,.313,.313,.0,.0,.0,.0,.0,.0,.055,.091,.124,.131,.142,.171,.168,.162,.0,.0,.0,.028,.093,.129,.142,.162,.181,.191,.193,.195,.0,.0,.0,.075,.099,.122,.139,.133,.148,.167,.177,.0,.101,.101,.173,.054,.219,.247,.335,.349,.28,.280,.0,.0,.0,.0,.055,.115,.151,.168,.172,.2,.236,.264,.287,.298,.0,.0,.0,.133,.142,.150,.15,.173,.206,.265,.307,.331,.096,.096,.066,.113,.123,.13,.2,.281,.334,.334,.0,.0,.0,.0,.0,.084,.078,.115,.130,.191,.262,.294,.311,.322/
     
!	************************************************************
!	**********End of Reference library**************************
!	************************************************************
	
!  *********************************************************************	
!  *******DETERMINE THE APPROPRIATE TRIANGLE FOR INTERPOLATION**********
!  *********************************************************************

! define the equations of line, distance and area

	yyy(aa,xd,bb)=aa*xd+bb
	xxx(aa,yd,bb)=(yd-bb)/aa
	d(x2,x1,y2,y1)=((x2-x1)**2+(y2-y1)**2)**.5
	at(aa,bb,cc)=0.5*bb*cc*(1-((bb**2+cc**2-aa**2)/(2*bb*cc))**2)**.5

!  inititalize variables
	x = xoc
	y = yhc
	ind = yf
	
!  look up table of the reference points of the triangular mesh

	data xx/.0177734,.0203654,.0659401,.0773465,.0893623,.1077369,.1305803,.1357569,.1803479,.2093441,.2603201,.0687/
	data yy/.6717240,.5810955,.6550527,.8088697,.7575807,.8506428,.7671163,.8523190,.8499221,.7890888,.8572938,.863/
     
!  look up table for a and b of equations of all triangle sides

	data a/-34.965,1.6228,-.34612,2.3021,1.7337,1.1993,4.3774,-4.2685,.23134,5.0647,1.3746,-3.6565,.059818,16.459,1.6638,-.05375,.27897,-2.0979,.092179,1.3380,3.7559,-6.2604,-.31655/
	data b/1.2932,.54805,.67788,.63081,.54074,.65041,.36641,1.1390,.73691,.30499,.70255,1.2446,.84420,-1.3821,.54985,.85962,.73069,1.2283,.83330,.50899,.60497,1.2931,.88475/
     
!  look up table for the three sides that correspond to each triangle

	data s1/1,3,4,7,8,10,12,14,15,18,21,22/
	data s2/2,7,6,5,10,9,14,15,17,20,4,11/
	data s3/3,6,8,9,11,12,13,16,18,19,22,23/
	
!  do loop to find the appropriate triangle for interpolation
!	if(iii.eq.1)then
	m=0.0
	inside = .true.	  
	DO 10 i=1,12
	   z1=xxx(a(s1(i)),y,b(s1(i)))
	   z2=xxx(a(s2(i)),y,b(s2(i)))
	   z3=yyy(a(s3(i)),x,b(s3(i)))
	  if(x.GE.z1.AND.x.LE.z2.AND.y.LE.z3)THEN
	     m=i  
	  go to 30
	  endif
	  if(i.EQ.12.AND.m.EQ.0.0)Then
	     print*,'one or both ratios are out of bounds'
	     inside = .false.
	  endif    
10	continue
30	yyyy=.true.
!	endif
!  *****************************************************************	
!  *****************TRIANGULAR INTERPOLATION************************
!  *****************************************************************
	if(inside)then
!  This interpolation scheme is taken from Zhao et al., 25th Symp. on Comb. (Int.), 1994/pp. 553-560. 

!  look up table for points 1,2, and 3 for each triangle

	data p1/2,3,1,3,5,5,7,7,7,10,1,4/
	data p2/1,1,4,5,4,6,6,8,9,9,12,12/
	data p3/3,5,5,7,6,7,8,9,10,11,4,6/
	
!  calculate the length of each side
!	if(iii.eq.1)then
	ds1=d(xx(p1(i)),xx(p2(i)),yy(p1(i)),yy(p2(i)))
	ds2=d(xx(p3(i)),xx(p1(i)),yy(p3(i)),yy(p1(i)))
	ds3=d(xx(p3(i)),xx(p2(i)),yy(p3(i)),yy(p2(i)))
	ds4=d(x,xx(p2(i)),y,yy(p2(i)))
	ds5=d(xx(p1(i)),x,yy(p1(i)),y)
	ds6=d(xx(p3(i)),x,yy(p3(i)),y)
	
!	print*,' ds1,ds2,ds3 = ',ds1,ds2,ds3

!  calculate the area of each triangle used in interpolation scheme

	A1=at(ds1,ds2,ds3)
	A2=at(ds1,ds5,ds4)
	A3=at(ds5,ds2,ds6)
	
!	print*,' A1,A2,A3 = ',A1,A2,A3
	
!  calculate S and R, the weighted fraction of two of the points  the weighted fraction of other point will be 1-R-S

	S=A2/A1
	R=A3/A1
!	endif
!	************************************************************
!	*************CALCULATE LIGHT GAS DISTRIBUTION***************
!	************************************************************
!  n is the number of order pairs of data for coals 1-12

	data n/10,9,11,10,12,12,11,10,12,12,10,12/

!  do loop to calculate the light gas composition of each reference coal (point) of triangle.  this is accomplished using linear  interpolation between points in reference library
!  j specifies the point (1,2, or3) of the triangle selected above

	DO 35 j=1,3
	  if(j.EQ.1)THEN
	    lib=p1(i)
	  elseif(j.EQ.2)THEN
	    lib=p2(i)
	  else
	    lib=p3(i)
	  endif
	
!  do loop to find the two xz points that surround ind

	DO 40 ii=1,12
	  if(ind.GE.xz(ii,lib).AND.ind.LE.xz(ii+1,lib))THEN
	  go to 45
	  endif
40	continue

!	for ygas(k,j)
!		k=1;fh2o
!		k=2;fco2
!		k=3;fch4
!		k=4;fco
!	        k=5;other light gases (HC'2, H2, parafins, olifins)

!  linear interpolation to find reference values as a function of ind

45      ygas(1,j)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*fh2o(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-xz(ii,lib)))*fh2o(ii+1,lib)
        ygas(2,j)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*fco2(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-xz(ii,lib)))*fco2(ii+1,lib)
        ygas(3,j)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*fch4(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-xz(ii,lib)))*fch4(ii+1,lib)
        ygas(4,j)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*fco(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-xz(ii,lib)))*fco(ii+1,lib)
     
!	  print*,'lib=',lib
!	  print*,'ii=',ii
!	  print*,'xz=',xz(ii,lib)

35	continue
	endif
!	************************************************************
!	*******CALCULATE GAS COMPOSITION FROM LIBRARY COALS*********
!	************************************************************
	if(inside)then
!  yygas(k) = the fraction of light gas that is k,  1=h20, 2=co2, 3=ch4, 4=co

	DO 50 k=1,4	      
	  yygas(k)=(1-R-S)*ygas(k,1)+R*ygas(k,2)+S*ygas(k,3)
!	  print*,ygas(k,1),ygas(k,2),ygas(k,3)
!	  print*,'xz(5,1)',xz(5,1)
50	continue
	endif
	
!	************************************************************
!	*****Estimate composition for coals outside mesh************
!	************************************************************
	
	if(inside.eqv..false.)then
	data out/.0,.085,.835,1.5,12,.085,.12,.835,1.5,6,.12,.155,.835,1.5,8,.155,.222,.835,1.5,9,.222,.285,.75,1.5,11,.0,.089,.75,.835,4,.0,.05,.63,.75,1,.05,.175,.63,.69,3,.066,.175,.69,.835,7,.175,.222,.0,.835,10,.0,.175,.55,.63,2,.222,1.,.0,.75,10,.285,1,0.75,1.5,13,.0,.175,.0,.55,14/
     
        DO 60 kk=1,14
            if(x.GT.out(1,kk).AND.x.LE.out(2,kk))then
               if(y.GE.out(3,kk).AND.y.LT.out(4,kk))then
            lib=out(5,kk)
            if(lib.eq.13)then
              yygas(1)=0.24
              yygas(2)=0.37
              yygas(3)=0.06
              yygas(4)=0.28
              go to 70
            endif
            if(lib.eq.14)then
              yygas(1)=0.18
              yygas(2)=0.08
              yygas(3)=0.37
              yygas(4)=0.18
              go to 70
            endif
!  do loop to find the two xz points that surround ind

	DO 63 ii=1,12
	  if(ind.GE.xz(ii,lib).AND.ind.LE.xz(ii+1,lib))THEN
	  go to 64
	  endif
63	continue

!  linear interpolation to find reference values as a function of ind

64      yygas(1)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*fh2o(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-xz(ii,lib)))*fh2o(ii+1,lib)
        yygas(2)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*fco2(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-xz(ii,lib)))*fco2(ii+1,lib)
        yygas(3)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*fch4(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-xz(ii,lib)))*fch4(ii+1,lib)
        yygas(4)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*fco(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-xz(ii,lib)))*fco(ii+1,lib)
               print*,out(1,kk),out(2,kk),out(3,kk),out(4,kk),out(5,kk)
               go to 65
               endif
            endif
60	continue
65      print*,'Light gas distribution is based on ref. # ',lib
!	print*,out(1,6),out(2,6),out(3,6),out(4,6),out(5,6)
	print*,x,y
	endif
	
70	yygas(5)=1-yygas(1)-yygas(2)-yygas(3)-yygas(4)
!	print*,'R and S = ',R,S
!	print*,'fh2o, fco2, fch4, fco, & fother =',yygas(1),yygas(2),yygas(3),yygas(4),yygas(5)
!	print*,inside,m,x,y,ind
!	nn=ntmax-2
   	return    
90	end
