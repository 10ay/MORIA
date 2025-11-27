#define _LDIM_  20000
#define _MDIM_     20
#define _NDIM_     40
#define _LINUX_ .true.
c
c WARNING: Jay Anderson's routines to read and write fits files do not work
c          with Intel Fortran. Use gfortran instead!
c

      program uvp2tri_scon_fs_asym_mcmc
      implicit none
  
      integer NMAX, MA, MAp1
      real toler
      PARAMETER(NMAX=1500,toler=1.e-6,MA=35,MAp1=MA+1)
      REAL Aa(MAp1),V(MA,MA),uu(NMAX,MA),W(MA),BB(NMAX)
      real Aamin(MAp1)
      integer   L, Ls, ier, img, img1
      real      pp_psky(_LDIM_),pp_pskyu(_LDIM_)
      real      pp_pskyu_min(_LDIM_)
      real wthresh, WMAX
      real*8    ul(_LDIM_)
      real*8    vl(_LDIM_)
      real      pl(_LDIM_), ppl(_LDIM_)
      integer   il(_LDIM_)
      integer   jl(_LDIM_)
      real*8    xl(_LDIM_)
      real*8    yl(_LDIM_)
      real     dul(_LDIM_)
      real     dvl(_LDIM_)
      real*8    rl(_LDIM_)
      integer niml(_LDIM_)
      integer nstl(_LDIM_)
      integer count1,count2,i,j
      integer U, Us, Us0, nim_maxp1
      real*8 iu(_LDIM_), ju(_LDIM_), Io, Jo
      real*8 i0(_NDIM_), j0(_NDIM_)
      integer NIM, NIMs 

      real*8 AGN(_NDIM_), BGN(_NDIM_), CGN(_NDIM_), DGN(_NDIM_)
      real*8 GAN(_NDIM_), GBN(_NDIM_), GCN(_NDIM_), GDN(_NDIM_)
      real*8 DLTA

      integer Ss, Su 
      real*4 ps(999), pbar, psig
      real*4 po(999)
      double precision RN
      real    x, y,prob,prob_rand
      integer O, Os, Ou
      real*4 fo(999), fbar, fsig
      real*4 do(999), dbar, dsig

      real rpsf_photij
      real psf78u(101,101,7,8)
      common /psf78_/psf78u

      real SKYN(_NDIM_)

      real*8  ptot, ftot
      real*8  zo,zm
      real*8  chisq, chisq_min
      integer if2, ifit_f3
      real*8  f, f1, f2

      integer NST, NSTs
      real    sky_nm(_NDIM_,_MDIM_)

      real    rpsf_phot
      real    psfu(101,101)

      real    show_psf(501,1001)
      real    show_str(2501,2001)
      real    xpsf, ypsf, rpsf
      integer ipsf, jpsf
      real    z_n(_NDIM_)     
      real    z_mn(_MDIM_,_NDIM_), zbar, zsig      
      real   sz_mn(_MDIM_,_NDIM_)     
      integer u_mn(_MDIM_,_NDIM_)     
      real    e, find_error,sep,sep1 
      real    sep12x, sep12y, sig12x, sig12y, sep12xchk, sep12ychk
      real*8  eeo, eem

      character*80 STRING,rmpixfile
      character*80 pixoutfile,skyoutfile,puseoutfile,probeoutfile
      character*80 mcmcoutfile,finaloutfile,showoutfile,psffile,infile
      character*16 tag,tagin
      integer NIT, iseed

      real ix1o, iy1o,eem_chi2
      real ix2o, iy2o, ix3o,iy3o, ix3m, iy3m
      real    if, ifo,ifm, if1o, if1m,if1
      real    chi2,ix2m,iy2m,ix1m,iy1m
      integer nimu(_LDIM_)
      integer nstu(_LDIM_)
      real     duu(_LDIM_)
      real     dvu(_LDIM_)
      real      fu(_LDIM_),  fu_min(_LDIM_)
      real     ppu(_LDIM_), ppu_min(_LDIM_)
      real     sgu(_LDIM_)
      real x1, y1, x1min, y1min, x3, x3min
      real x2, y2, x2min, y2min, y3, y3min
      real dr1, dr2, dr3, df1, df2
      real x1o, y1o, x2o, y2o, x3o, y3o, f1o, f2o
      real x1in, y1in, x2in, y2in, x3in, y3in, f1in, f2in
      real  f1min, f2min, fudge, fudge2
      real zmin,emin,chi2min,sum,chi2scon,chi2scmin
      real dufitmn,dufitmx,dvfitmn,dvfitmx,chi2cut,chi2max,chi2tot,chi2o
      real chi2norm
      integer irenorm
      integer ichi2, idof
      integer Nmcmc, nrepeat, ltag, lenc, naccept, nreject
      integer ltagin,lpsff,lpixf
c       double precision rn
      real pix_show(4001,2001)
      COMMON/RN_COM/ISEED

       iseed = 45780
       count1 = 0
c---------------------------------------------------------------
c
c

      write(6,*) 'enter input PSF file:'
      read(5,'(a)') psffile
      lpsff = lenc(psffile)
      write(6,*) 'enter input file tag:'
      read(5,'(a)') tagin
      ltagin = lenc(tagin)
      write(6,*) 'enter output file tag:'
      read(5,'(a)') tag
      ltag = lenc(tag)
      pixoutfile = 'uvp2tri_scon_fsky_'//tag(1:ltag)//'.01.pix_all'
      lpixf = lenc(pixoutfile)
      rmpixfile = 'uvp2tri_scon_fsky_'//tag(1:ltag)//'.08.rm_pix'

      open(21,file=pixoutfile,status='new')
      infile = 'img2extract_wfc3uv_psflist_'//tagin(1:ltagin)//'.uvp'
      open(63,file=infile,status='old')
      NIMs = 0
      NSTs = 0
      Ls = 0

c     1st 12 parameteres in the img2extract_wfc3uv_psflist.uvp file:
c     (1) u position in the reference frame (u,v) are ref frame cords
c     (2) v position in the reference frame
c     (3) p, the pixel value for the pixel in question
c     (4) exposure number
c     (5) i coordinate of pixel 
c     (6) j coordinate of pixel
c     (7) x-coordinate of distortion-corrected location of the pixel
c     (8) y-coordinate of distortion-corrected location of the pixel
c     (9)  delta-u (offset from the center of this object)
c     (10) delta-v (offset from the center of this object)
c     (11) sqrt[(delta-u)^2+(delta-v)^2]
c     (12) star number in the NEARBY_SIM list (ignore last column)
c     Ls = the pixel number for a given star

  1   read(63,*,end=2) ul(Ls+1), vl(Ls+1), pl(Ls+1),
     .              niml(Ls+1), il(Ls+1), jl(Ls+1),
     .                xl(Ls+1), yl(Ls+1), 
     .               dul(Ls+1),dvl(Ls+1), rl(Ls+1),
     .              nstl(Ls+1) 
c     This next line selects the target star. For photomtry of other stars
c     pick a different number to set nstl(Ls+1) equal to.
      if (nstl(Ls+1).ne.1) goto 1
      Ls = Ls + 1
      if (Ls.gt._LDIM_) stop 'Ls.gt._LDIM_'
      if (niml(Ls).gt.NIMs) NIMs = niml(Ls)
      if (nstl(Ls).gt.NSTs) NSTs = nstl(Ls)
      if (NIMs.gt._MDIM_) stop 'NIMs.gt._MDIM_'
      if (NSTs.gt._NDIM_) stop 'NIMs.gt._NDIM_'
      write(21,121) Ls, dul(Ls),dvl(Ls),pl(Ls),niml(Ls),nstl(Ls)
  121 format(1x,i5.5,1x,f7.3,1x,f7.3,1x,f9.1,1x,i4,1x,i4)
      goto 1
    2 continue
      close(21)
      close(63)

      print*,'      '
      print*,'NSTs: ',NSTs
      print*,'NIMs: ',NIMs
      print*,'  Ls: ',Ls
      print*,'      '

c      input sep12x, sep12y, sig12x, sig12y
c      ------------------------------------
       write(6,*)
     &    'enter constraint parameters sep12x, sig12x, sep12y, sig12y:'
       read(5,*) sep12x, sig12x, sep12y, sig12y

c      input initial coordinate and flux ratio values
c      ----------------------------------------------
       write(6,*) 'enter initial x1, y1, x2, y2, x3, y3, f1, f2 values:'
       write(6,*) 'If abs(f1+f2-1.d0) < 1.d-6, source 3 flux f3 == 0'
       read(5,*) x1in, y1in, x2in, y2in, x3in, y3in, f1in, f2in
       if(abs(f1in+f2in-1.d0).le.1.d-6) then
         ifit_f3 = 0
       else
         ifit_f3 = 1
       endif
       write(6,*) 'enter dr1, dr2, dr3, df1, df2:'
       read(5,*) dr1, dr2, dr3, df1, df2
       write(6,*) 'enter error bar fudge factor:'
       read(5,*) fudge
       fudge2 = fudge**2
       write(6,*) 'enter Nmcmc:'
       read(5,*) Nmcmc
       write(6,*) 'enter dufitmn, dufitmx, dvfitmn, dvfitmx, chi2cut:'
       read(5,*) dufitmn, dufitmx, dvfitmn, dvfitmx, chi2cut

c-------------------------------------------------------------
cccc
cccc     This routine then determines a sky background value for each exposure,
cccc      and this is output in the second output file (uvp2dual_fit.02.sky_out).
ccc      skyoutfile = 'uvp2tri_scon_fsky_'//tag(1:ltag)//'.02.sky_all'
ccc      open(22,file=skyoutfile,status='new')
ccc      do NST = 1, NSTs
ccc      do NIM = 1, NIMs
ccc         Ss = 0
ccc         do L = 1, Ls
ccc            if (niml(L).eq.NIM.and.
ccc     .          nstl(L).eq.NST.and.
ccc     .          rl(L).gt.8.5.and.rl(L).lt.13.5) then
ccc                Ss = Ss + 1
ccc                ps(Ss) = pl(L)
ccc                endif
ccc            enddo
ccc         pbar = 0.
ccc         psig = 0.
ccc         Su = 0
ccc         if (Ss.gt.0) call barsig(ps,Ss,pbar,psig,Su)
ccc         SKYN(NIM) = pbar
ccc         write( *,122) NIM,NST,pbar,psig,Ss,Su
ccc         write(22,122) NIM,NST,pbar,psig,Ss,Su
ccc  122    format('SKY: ',1x,i2.2,3x,i2.2,5x,f8.3,1x,f8.3,1x,i3,1x,i3)
ccc         sky_nm(NST,NIM) = pbar
ccc         enddo
ccc         enddo
ccc      close(22)

c     read in the PSF model file
c     --------------------------
      print*,'READIN ', psffile(1:lpsff)
      call readfits_r4(psffile,psfu,101,101)

      print*,'             '
      print*,'psfu(51,51): ',psfu(51,51)
      print*,'             '

      do L = 1, Ls
         NST = nstl(L)
         NIM = niml(L)
ccc         ppl(L) = pl(L) - sky_nm(NST,NIM)
         pp_psky(L) = pl(L)
      enddo

c-------------------------------------------------------------
c
c
      puseoutfile = 'uvp2tri_scon_fsky_'//tag(1:ltag)//'.03.pix_use'
      open(23,file=puseoutfile,status='new')
      Us = 0
      do L = 1, Ls
ccc         if (abs(dul(L)).lt.dufitmx.and.
ccc     .       abs(dvl(L)).lt.dvfitmx) then
         if (dul(L).gt.dufitmn.and.dul(L).lt.dufitmx.and.
     .       dvl(L).gt.dvfitmn.and.dvl(L).lt.dvfitmx) then
             Us = Us + 1
             nimu(Us) = niml(L)
             nstu(Us) = nstl(L)
             duu(Us)  =  dul(L)
             dvu(Us)  =  dvl(L)
             ppu(Us)  =  ppl(L)
             pp_pskyu(Us) = pp_psky(L)
             sgu(Us)  =  sqrt(16.0
     .                       +max(pl(L),0.)
     .                       +(0.01*max(pl(L),0.))**2)
             write( *,177) L,Us,duu(Us),dvu(Us),
     .                          pp_pskyu(Us),nimu(Us),nstu(Us)
             write(23,177) L,Us,duu(Us),dvu(Us),
     .                          pp_pskyu(Us),nimu(Us),nstu(Us)
  177        format(1x,i5.5,1x,i4.1x,2f10.4,1x,f10.1,1x,2i4)
             endif
         enddo
      close(23)

      if(Us.gt.NMAX) then
        write(6,*) 'Us =',Us,' > NMAX =',NMAX
        stop
      endif
 
      Us0 = Us
      print*,'Ls: ',Ls
      print*,'Us: ',Us

      probeoutfile = 'uvp2tri_scon_fsky_'//tag(1:ltag)//'.04.probe_fit'
      mcmcoutfile = 'uvp2tri_scon_fsky_'//tag(1:ltag)//'.07.mcmc'

      chi2norm = 1.0
      irenorm = 0

 180  continue

      write(6,*) 'chi2norm =',chi2norm

c     reset initial conditions after removing a pixel
c     -----------------------------------------------
      x1 = x1in
      y1 = y1in
      x2 = x2in
      y2 = y2in
      x3 = x3in
      y3 = y3in
      f1 = f1in
      f2 = f2in
      x1o = x1
      y1o = y1
      x2o = x2
      y2o = y2
      x3o = x3
      y3o = y3
      f1o = f1
      f2o = f2

      open(24,file=probeoutfile,status='unknown')
      open(25,file=mcmcoutfile,status='unknown')
      open(26,file=rmpixfile,status='unknown')

c-----------------------------------------------------------------------------
c MCMC starting from this point
c---------------------------------------------------------------------------
c Aparna's initial conditions:
c      x1 = 0.600
c      y1 = 0.200
c      x2 = -0.400
c      y2 = -0.150
c      x3 = 0.200
c      y3 = -2.800
c      f1 = 0.400
c      f2 = 0.300
c      dr1 = 0.025
c      dr2 = 0.025
c      dr3 = 0.025
c      df1 = 0.025
c      df2 = 0.025

       idof = Us - 9 - NIMs

       do U = 1, Us
         fu(u) = (f1)*rpsf_phot(duu(u)-x1,dvu(u)-y1,psfu)
     .         + (f2)*rpsf_phot(duu(u)-x2,dvu(u)-y2,psfu)
     .         + (1-f1-f2)*rpsf_phot(duu(u)-x3,dvu(u)-y3,psfu)
       enddo

c      simultaneous fit sky method
c      ---------------------------
       if(nimu(Us)+1.gt.MA) then
         write(6,*) 'nimu(Us)+1',nimu(Us)+1,' > MA =',MA
         stop
       endif
       nim_maxp1 = NIMs + 1
       do U = 1,Us
         BB(U) = pp_pskyu(U)/sgu(U)
         uu(U,1) = fu(U)/sgu(U)
         do img = 1,NIMs
           uu(U,img+1) = 0.
         enddo
         uu(U,nimu(U)+1) = 1./sgu(U)
       enddo
       call SVDCMP(uu,Us,nim_maxp1,nmax,MA,W,V,ier)
       if(ier.ne.0) STOP 'SVDCMP failure'
       WMAX=0.
       do j = 1,MA
         if(W(j).gt.WMAX) WMAX=W(j)
       enddo
       wthresh=toler*WMAX
       do j = 1,MA
         if(W(j).LT.wthresh) W(j)=0.
       enddo
       CALL SVBKSB(uu,W,V,Us,nim_maxp1,nmax,MA,BB,Aa)

c      calculate chi^2
c      ---------------
       chi2tot = 0.
       do U = 1,Us
         sum = Aa(1)*fu(U) + Aa(nimu(U)+1)
         chi2 = ((pp_pskyu(U)-sum)/sgu(u))**2
         chi2tot = chi2tot + chi2
       enddo
       chi2tot = chi2tot/chi2norm

       chi2scon = ((sep12x-(x1-x2))/sig12x)**2
     &            + ((sep12y-(y1-y2))/sig12y)**2
       chi2tot = chi2tot + chi2scon

       chi2o = chi2tot

c      save initial values
c      -------------------
       x1min = x1
       y1min = y1
       x2min = x2
       y2min = y2
       x3min = x3
       y3min = y3
       f1min = f1
       f2min = f2
       chi2min = chi2tot
       chi2scmin = chi2scon
       do U = 1, Us
         fu_min(U)  =  fu(U)
         pp_pskyu_min(U) = pp_pskyu(U)
       enddo
       do img1 = 1,nim_maxp1
         Aamin(img1) = Aa(img1)
       enddo

c      write initial state to MCMC file
c      --------------------------------
       write(24,993) 0,x1,y1,x2,y2,x3,y3,f1,f2,chi2tot,
     &               Aa(1),(Aa(img+1), img = 1,NIMs)
       write(25,993) 0,x1,y1,x2,y2,x3,y3,f1,f2,chi2tot,
     &               Aa(1),(Aa(img+1), img = 1,NIMs)

c------------------------------------------------------------
       naccept = 0
       nreject = 0
       do i=1, Nmcmc
         x1 = x1o + dr1*(2.*rn()-1.)
         y1 = y1o + dr1*(2.*rn()-1.)
         x2 = x2o + dr2*(2.*rn()-1.)
         y2 = y2o + dr2*(2.*rn()-1.)
         x3 = x3o + dr3*(2.*rn()-1.)
         y3 = y3o + dr3*(2.*rn()-1.)
 22      f1 = f1o + df1*(2.*rn()-1.)
         if(ifit_f3.eq.0) then
           f2 = 1.-f1
           x3 = x3o
           y3 = y3o
         else
           f2 = f2o + df2*(2.*rn()-1.)
         endif
         if((f1+f2.gt.1.0).or.(f1.lt.0.).or.(f2.lt.0.)) go to 22

         do U = 1, Us
           fu(u) = (f1)*rpsf_phot(duu(u)-x1,dvu(u)-y1,psfu)
     .             + (f2)*rpsf_phot(duu(u)-x2,dvu(u)-y2,psfu)
     .             + (1-f1-f2)*rpsf_phot(duu(u)-x3,dvu(u)-y3,psfu)
         enddo
c        simultaneous fit sky method
c        ---------------------------
         if(nimu(Us)+1.gt.MA) then
           write(6,*) 'nimu(Us)+1',nimu(Us)+1,' > MA =',MA
           stop
         endif
         nim_maxp1 = NIMs + 1
         do U = 1,Us
ccc           BB(U) = pp_pskyu(U)/sgu(U)  ! already calculated above
           uu(U,1) = fu(U)/sgu(U)
           do img = 1,NIMs
             uu(U,img+1) = 0.
           enddo
           uu(U,nimu(U)+1) = 1./sgu(U)
         enddo
         call SVDCMP(uu,Us,nim_maxp1,nmax,MA,W,V,ier)
         if(ier.ne.0) STOP 'SVDCMP failure'
         WMAX=0.
         do j = 1,MA
           if(W(j).gt.WMAX) WMAX=W(j)
         enddo
         wthresh=toler*WMAX
         do j = 1,MA
           if(W(j).LT.wthresh) W(j)=0.
         enddo
         CALL SVBKSB(uu,W,V,Us,nim_maxp1,nmax,MA,BB,Aa)

c        calculate chi^2
c        ---------------
         chi2tot = 0.
         do U = 1,Us
           sum = Aa(1)*fu(U) + Aa(nimu(U)+1)
           chi2 = ((pp_pskyu(U)-sum)/sgu(u))**2
           chi2tot = chi2tot + chi2
         enddo
         chi2tot = chi2tot/chi2norm
         chi2scon = ((sep12x-(x1-x2))/sig12x)**2
     &            + ((sep12y-(y1-y2))/sig12y)**2
         chi2tot = chi2tot + chi2scon

         if (chi2tot .lt. chi2min) then
           x1min = x1
           y1min = y1
           x2min = x2
           y2min = y2
           x3min = x3
           y3min = y3
           f1min = f1
           f2min = f2
           chi2min = chi2tot
           chi2scmin = chi2scon
           do U = 1, Us
             fu_min(U) = fu(U)
             pp_pskyu_min(U) = pp_pskyu(U)
           enddo
           do img1 = 1,nim_maxp1
             Aamin(img1) = Aa(img1)
           enddo
         endif  
        
         if (chi2tot .le. chi2o) then
           naccept = naccept + 1
           x1o = x1
           x2o = x2
           y1o = y1
           y2o = y2
           x3o = x3
           y3o = y3
           f1o = f1
           f2o = f2
           chi2o = chi2tot
           nrepeat = 0
           write(25,993) nrepeat,x1,y1,x2,y2,x3,y3,f1,f2,chi2tot,
     &                   Aa(1),(Aa(img+1), img = 1,NIMs)
         else
           prob_rand = rn()
           prob = exp(-(chi2tot-chi2o)/2.00)
           if (prob .gt. prob_rand) then
             naccept = naccept + 1
             x1o = x1
             x2o = x2
             y1o = y1
             y2o = y2
             x3o = x3
             y3o = y3
             f1o = f1
             f2o = f2
             chi2o = chi2tot
             nrepeat = 0
             write(25,993) nrepeat,x1,y1,x2,y2,x3,y3,f1,f2,chi2tot,
     &                     Aa(1),(Aa(img+1), img = 1,NIMs)
           else
             nreject = nreject + 1
             nrepeat = nrepeat + 1
             write(25,992) nrepeat
           endif
           call flush(25)
         endif

ccc         write(*,129) x1,y1,x2,y2,ifm,zm,eem,
ccc     .                x1min,y1min,x2min,y2min,f1min,f2min,zmin,emin
         if((i/10000)*10000.eq.i) write(6,*) i,' steps completed'
       enddo
       write(6,*) 'accepted, rejected MCMC steps:',naccept,nreject

c      check if a pixel should be removed
c      ----------------------------------
       chi2max = 0.
       ichi2 = 0
       do U = 1,Us
         sum = Aamin(1)*fu_min(U) + Aamin(nimu(U)+1)
         chi2 = ((pp_pskyu_min(U)-sum)/sgu(u))**2
         chi2 = chi2/chi2norm
         if(chi2.gt.chi2max) then
           ichi2 = U
           chi2max = chi2
         endif
       enddo

c      possibly remove pixel
c      ---------------------
       if(ichi2.gt.0) then
         if(chi2max.ge.chi2cut) then
           write(6,994) nimu(ichi2),duu(ichi2),dvu(ichi2),
     &                  pp_pskyu(ichi2),chi2max
           write(26,994) nimu(ichi2),duu(ichi2),dvu(ichi2),
     &                  pp_pskyu(ichi2),chi2max
           do U = ichi2, Us-1
             nimu(U) = nimu(U+1)
             nstu(U) = nstu(U+1)
             duu(U) = duu(U+1)
             dvu(U) = dvu(U+1)
             ppu(U) = ppu(U+1)
             pp_pskyu(U) = pp_pskyu(U+1)
             sgu(U) = sgu(U+1)
           enddo
           Us = Us - 1
           idof = Us - 9 - NIMs
           chi2norm = (chi2min-chi2max)*chi2norm/idof
           close(24)
           close(25)
           irenorm = 1
           go to 180
         elseif(irenorm.eq.0) then
           chi2norm = chi2min*chi2norm/idof
           irenorm = 1
           go to 180
         endif
       endif

 129   format(2x,4f8.3,f9.4,f13.1,g14.5,
     .        9x,4f8.3,f13.4,f13.1,g14.5)
 992   format(i8)
 993   format(i4,8f10.5,2f13.4,30f10.4)
 994   format('removing pix at img, du, dv, val, chi2 =',i4,2f10.4,f10.2,
     &        g13.5)


36     continue
      close(24)
      close(25)
      close(26)
c  ALL min are changed o. ppumin and fumin remain same.
      finaloutfile = 'uvp2tri_scon_fsky_'//tag(1:ltag)//'.05.final_fit'
      open(44,file=finaloutfile,status='new')
      write(44,990) x1in,y1in,x2in,y2in,x3in,y3in,f1in,f2in,
     &      dr1,dr2,dr3,df1,df2,Nmcmc,sep12x,sig12x,
     &      sep12y,sig12y,fudge,dufitmn,dufitmx,dvfitmn,dvfitmx,chi2cut
 990  format('#input  x1   y1       x2       y2       x3       y3    ',
     &       '   f1        f2',/,'#',8f9.5,/,
     &       '#  dr1      dr2      dr3      df1      df2      Nmcmc',/,
     &       '#',5f9.5,i8,/,
     &       '#    sep12x         sig12x         sep12y         ',
     &       'sig12y      fudge',/,'#',5g15.7,/,
     &       '#   dufitmn   dufitmx   dvfitmn  dvfitmx  chi2cut',/,
     &       '#',5f10.5)
      write(44,'(''# '')')
      write(44,991)
 991  format('# input file:')
      write(44,'(''#  PSF  file ='',a)') psffile(1:lpsff)
      write(44,'(''# pixel output file ='',a)') pixoutfile(1:lpixf)
      write(44,'(''# '')')
      write(44,'(''# BEST-FIT MODEL '')')
      write(44,'(''# pix0  = '',i10)') Us0
      write(44,'(''# pix     '',i10)') Us
      write(44,'(''#   X1MIN = '',f10.5)') x1min
      write(44,'(''#   Y1MIN = '',f10.5)') y1min
      write(44,'(''#   X2MIN = '',f10.5)') x2min
      write(44,'(''#   Y2MIN = '',f10.5)') y2min
      write(44,'(''#   X3MIN = '',f10.5)') x3min
      write(44,'(''#   Y3MIN = '',f10.5)') y3min
      write(44,'(''#   F1MIN = '',f10.5)') f1min
      write(44,'(''#   F2MIN = '',f10.5)') f2min
      write(44,'(''# chi2min = '',f10.4)') chi2min
      write(44,'(''#chi2scmin= '',f10.4)') chi2scmin
      write(44,'(''#  A1_min = '',g15.7)') Aamin(1)
      do img = 1,NIMs
        write(44,114) img, Aamin(img+1)
 114    format('# Sky(',i2,')  =',f10.4)
      enddo
      write(44,'(''# '')')
      write(44,'(''# INDIVIDUAL PIXELS:'')')
      write(44,'(''# '')')
      chi2scon = ((sep12x-(x1-x2))/sig12x)**2
     &         + ((sep12y-(y1-y2))/sig12y)**2
      chi2tot = chi2scon
      do U = 1,Us
        sum = Aa(1)*fu(U) + Aa(nimu(U)+1)
        chi2 = ((pp_pskyu(U)-sum)/sgu(u))**2
        chi2 = chi2/chi2norm
        chi2tot = chi2tot + chi2
        write(44,113) U,duu(U),dvu(U),pp_pskyu_min(U),
     .                fu_min(U),Aa(1)*fu(U),sgu(U),
     .                rpsf_phot(duu(U),dvu(U),psfu),chi2
      enddo
  113 format(1x,i4,1x,f8.4,1x,f8.4,1x,f8.1,
     .          1x,f8.5,1x,f9.1,1x,f8.2,1x,f8.5,1x,g13.5)
      close(25)

      print*,'GENERATE PIX_SHOW...'
      do i = 0001, 2001
        if (i.eq.i/100*100) print*,'---> i: ',i
        do j = 0001, 2001
          x = 0.01*(i-1001)
          y = 0.01*(j-1001)
          Ss = 0
          do L = 1, Ls
            if (abs(dul(L)-x).le.0.5.and.
     .        abs(dvl(L)-y).le.0.5) then
              Ss = Ss + 1
              f = f1min*rpsf_phot(dul(l)-x1min,dvl(l)-y1min,psfu)
     .          + f2min*rpsf_phot(dul(l)-x2min,dvl(l)-y2min,psfu)
     .          + (1-f1min-f2min)
     &            *rpsf_phot(dul(l)-x3min,dvl(l)-y3min,psfu)
              po(Ss) = pp_psky(L) - Aamin(niml(L)+1)
              ps(Ss) = po(Ss) - f*Aamin(1)
            endif
          enddo
          call barsig(ps,Ss,pbar,psig,Su)
          pix_show(i+0000,j) = pbar
          call barsig(po,Ss,pbar,psig,Su)
          pix_show(i+2000,j) = pbar
        enddo
      enddo
c-----------------------------------------------------

      showoutfile = 'uvp2tri_scon_fsky_'//tag(1:ltag)//
     &     '.06.pix_show.fits'
      call writfits_r4(showoutfile,pix_show,4001,2001)

      STOP
      end




c--------------------------------------------------------------
c
c this routine does a least squares solution (with no
c data rejection or weighting) for the 6-param linear
c fit for: 
c
c   
c    x2 = A*(x1-x1o) + B*(y1-y1o) + x2o  
c    y2 = C*(x1-x1o) + D*(y1-y1o) + y2o  
c
c it may look like there are 8 params, but two of the
c offsets are arbitrary and are just set to the centroid
c of the distribution.
c
c

      subroutine glob_fit6nrDP(x1,y1,x2,y2,NUSE,
     .                         A,B,C,D,x1o,y1o,x2o,y2o)
      implicit none

      real*8 x1(1), y1(1)
      real*8 x2(1), y2(1)
      integer NUSE
      real*8 A, B, C, D
      real*8 x1o, y1o, x2o, y2o

      real*8 sx1o,sy1o
      real*8 sx2o,sy2o

      real*8 sxx, sx, swx, szx
      real*8 syy, sy, swy, szy
      real*8 sw,  sz, sxy
      real*8 dlta
      real*8 dsxx, dsyy, dsxy
      real*8 dswx, dswy, dszx, dszy

      integer n

      if (NUSE.lt.3) goto 999

      sx1o = 0
      sy1o = 0
      sx2o = 0
      sy2o = 0
      do n = 1, NUSE
         sx1o = sx1o + x1(n)
         sy1o = sy1o + y1(n)
         sx2o = sx2o + x2(n)
         sy2o = sy2o + y2(n)
         enddo

      x1o = sx1o/NUSE
      y1o = sy1o/NUSE
      x2o = sx2o/NUSE
      y2o = sy2o/NUSE
  
c     print*,'NUSE x1o: ',NUSE,x1o,sx1o
c     print*,'     y1o: ',NUSE,y1o,sy1o
c     print*,'     x2o: ',NUSE,x2o,sx2o
c     print*,'     y2o: ',NUSE,y2o,sy2o

      sxx = 0.0
      sx  = 0.0
      syy = 0.0
      sy  = 0.0
      swx = 0.0
      swy = 0.0
      szx = 0.0
      szy = 0.0
      sw  = 0.0
      sz  = 0.0
      sxy = 0.0
      do n = 1, NUSE
         sxy = sxy + (x1(n)-x1o)*(y1(n)-y1o)
         sxx = sxx + (x1(n)-x1o)*(x1(n)-x1o)
         sx  = sx  + (x1(n)-x1o)
         syy = syy + (y1(n)-y1o)*(y1(n)-y1o)
         sy  = sy  + (y1(n)-y1o)
         swx = swx + (x2(n)-x2o)*(x1(n)-x1o)
         swy = swy + (x2(n)-x2o)*(y1(n)-y1o)
         sw  = sw  + (x2(n)-x2o)
         szx = szx + (y2(n)-y2o)*(x1(n)-x1o)
         szy = szy + (y2(n)-y2o)*(y1(n)-y1o)
         sz  = sz  + (y2(n)-y2o)
         enddo


      dsxx = sx*sx - NUSE*sxx
      dsyy = sy*sy - NUSE*syy
      dsxy = sx*sy - NUSE*sxy
   
      dlta = dsxx*dsyy - dsxy*dsxy
   
      if (dlta.eq.0) goto 999

      dswx = sw*sx - NUSE*swx
      dswy = sw*sy - NUSE*swy
      dszx = sz*sx - NUSE*szx
      dszy = sz*sy - NUSE*szy

      A =   (dswx*dsyy-dswy*dsxy)/dlta
      B =   (dswy*dsxx-dswx*dsxy)/dlta
      C =   (dszx*dsyy-dszy*dsxy)/dlta
      D =   (dszy*dsxx-dszx*dsxy)/dlta

      return

 999  continue
      A  = 1
      B  = 0
      C  = 0
      D  = 1
      x1o = 0  
      y1o = 0
      x2o = 0
      y2o = 0
      return
      end




      subroutine barsig(xlist,NTOT,bar,sig,NUSE)
      implicit none
    
      integer NTOT
      real*4 xlist(NTOT)
      real*4 bar 
      real*4 sig
      integer NUSE

      integer n
      real*8    bsum, ssum
      integer nsum
      integer NIT


      bar = 0.e0
      sig = 9e9
      do NIT = 1, 30
         bsum = 0.
         ssum = 0.
         nsum = 0.
         do n = 1, NTOT
            if (abs(xlist(n)-bar).le.2.25*sig) then
               bsum = bsum + xlist(n)
               ssum = ssum + abs(xlist(n)-bar)
               nsum = nsum + 1
               endif
            enddo
         if (nsum.gt.0) bar = bsum / nsum 
         if (nsum.gt.1) sig = ssum/(nsum-1)
         enddo
      NUSE = nsum
      if (nsum.le.1) sig = 0.999 

      return
      end
 






c----------------------------------------------------
c
c this is the general function that evaluates the
c PSF for a given offset from the center and for
c a given location in the image.  If you have a pixel
c that is located at (ix,iy) and is (dx,dy) from the
c center of a star, this routine will tell you what
c fraction of the light should fall in that pixel.
c
      real function rpsf_photij(dx,dy,ix,iy)
      implicit none
      real    dx, dy
      integer ix, iy
 
      real psf78(101,101,7,8)
      common /psf78_/psf78
 
 
      real rpsf_phot
 
      real XMIN, XMAX
      real YMIN, YMAX
 
      real    fx, fy
      integer hx, hy
 
      hx = 1
      if (ix.gt.4096/6*1) hx = 2
      if (ix.gt.4096/6*2) hx = 3
      if (ix.gt.4096/6*3) hx = 4
      if (ix.gt.4096/6*4) hx = 5
      if (ix.gt.4096/6*5) hx = 6
 
      hy = 1
      if (iy.gt.4096/6*1) hy = 2
      if (iy.gt.4096/6*2) hy = 3
      if (iy.gt.4096/6*3) hy = 4
      if (iy.gt.4096/6*4) hy = 6
      if (iy.gt.4096/6*5) hy = 7
 
      XMIN = (hx-1)*4096/6
      XMAX = (hx  )*4096/6
      YMIN = (hy-1)*4096/6
      YMAX = (hy  )*4096/6
      if (hy.ge.4) then
         YMIN = (hy-2)*4096/6
         YMAX = (hy-1)*4096/6
         endif
 
      fx = (ix-XMIN)/(XMAX-XMIN)
      fy = (iy-YMIN)/(YMAX-YMIN)
 
 
c--------------------------------------------------
c
c linearly interpolate the value of the PSF at this
c (dx,dy) offset among the nearest 4 PSFs
c
      rpsf_photij
     .   = (1-fx)*(1-fy)*rpsf_phot(dx,dy,psf78(1,1,hx  ,hy  ))
     .   + (1-fx)*( fy )*rpsf_phot(dx,dy,psf78(1,1,hx  ,hy+1))
     .   + ( fx )*(1-fy)*rpsf_phot(dx,dy,psf78(1,1,hx+1,hy  ))
     .   + ( fx )*( fy )*rpsf_phot(dx,dy,psf78(1,1,hx+1,hy+1))
 
      return
      end
 
 
 
c------------------------------------------
c
c this will read-in the PSF from the file...
c
c
      subroutine readin_psf(FILENAME,psf78)
      implicit none
 
      character*80 FILENAME
      real psf78(101,101,7,8)
 
      real psfarr(1401,801)
      integer ipsf, jpsf
      integer i, j, jmax
      integer NIM
 
 
      call readfits_r4(FILENAME,psfarr,1401,801)
      do i = 1, 7
      do j = 1, 8
         do ipsf = 001, 101
         do jpsf = 001, 101
            psf78(ipsf,jpsf,i,j) = psfarr(ipsf+(i-1)*100,
     .                                    jpsf+(j-1)*100)
            enddo
            enddo
         enddo
         enddo
 
      print*,' '
      print*,' '
      do j = 8, 5, -1
         write(*,111) 2048+(j-5)*2048/3,(psf78(51,51,i,j),i=1,7)
         enddo
      do j = 4, 1, -1
         write(*,111) 0000+(j-1)*2048/3,(psf78(51,51,i,j),i=1,7)
  111    format(11x,i4.4,2x,7(f9.6,1x))
  112    format(11x,i4.2,2x,5(f9.6,1x),2x,i4.2)
         enddo
      print*,' '
      print*,' '
 
      return
      end
 

c--------------------------------------------------
c
c this just reads in an r4 fits image
c

      subroutine readfits_r4(FILE,pix,NDIMX,NDIMY)
      implicit none

      character*(*) FILE
      integer NDIMX,NDIMY
      real    pix(NDIMX,NDIMY)

      character*199 FILEU
      character*070 INFO(10)

      integer naxes
      integer laxis(3)

      character*8  field
      character*70 stream 
      
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer i,ios, k
      integer j

      character*2880 buffc
      byte   buffb(2880)
      equivalence (buffc,buffb)

      real*4 buffr(0720)
      integer ii,nn,nx,ny
      integer nxx, nyy
      integer ifirst, i1, i2

      integer np1, np2, npt
      integer nextend
      integer nread
      real bscale, bzero
      integer bitpix

      logical DIAG
      data DIAG /.false./

      character*70 HDR(25)
c     common/HDR/HDR

      logical extend_tf

      FILEU = FILE
      do i = 195,2,-1
         if (FILE(i:i+4).eq.'.fits') FILEU = FILE(1:i+4)
         enddo

      do i = 1, 25
         HDR(i) = ' '
         enddo


      if (DIAG) then
         print*,'enter readfits_r4...'
         print*,'FILE: ',FILE(1:60)
         endif
   
      open(10,file=FILEU,status='old',
     .     err=900,recl=2880,form='UNFORMATTED',
     .     access='DIRECT')

      if (DIAG) print*,'...opened'

      naxes = -1
      laxis(1) = 1
      laxis(2) = 1
      laxis(3) = 1 
      nextend = 0

      do i = 1, 10 
         INFO(i) = ' '
         enddo

      extend_tf = .false.

      i = 0
      nread = 0
 100  continue
      i = i + 1
      read(10,rec=i,iostat=ios) buffc
      if (DIAG) print*,'READREC: ',i
      do k = 0, 35, 1
         if (DIAG) write(*,'(i4,1x,i4,1x,a80)') 
     .                   i,k,buffc(k*80+1:k*80+80)
         do j = 1, 80 
            if (ichar(buffc(k*80+j:k*80+j)).eq.0)
     .                buffc(k*80+j:k*80+j) = ' '
            enddo
         field  = buffc(k*80+01:k*80+08)
         stream = buffc(k*80+10:k*80+79)
         if (DIAG) write(*,'(''   FIELD : '',8a)') field
         if (DIAG) write(*,'(''   STREAM: '',20a)') stream
         if (DIAG) write(*,'(''   STREAM: '',70a)') stream
         if (field.eq.'EXTEND  ') read(stream,*) extend_tf
         if (field.eq.'NAXIS   ') read(stream,*) naxes
         if (field.eq.'NAXIS1  ') read(stream,*) laxis(1)
         if (field.eq.'NAXIS2  ') read(stream,*) laxis(2)
         if (field.eq.'NAXIS3  ') read(stream,*) laxis(3)
         if (field.eq.'NEXTEND ') read(stream,*) nextend
         if (field.eq.'BITPIX  ') read(stream,*) bitpix
         if (field.eq.'BSCALE  ') read(stream,*) bscale
         if (field.eq.'BZERO   ') read(stream,*) bzero

         if (field.eq.'EXPTIME ') INFO(1) = stream
         if (field.eq.'FILTNAM1') INFO(2) = stream
         if (field.eq.'FILENAME') INFO(3) = stream
         if (field.eq.'DATE-OBS') INFO(4) = stream
         if (field.eq.'TIME-OBS') INFO(5) = stream
         if (field.eq.'DEC_TARG') INFO(6) = stream
         if (field.eq.'RA_TARG ') INFO(7) = stream
         if (field.eq.'PA_V3   ') INFO(8) = stream
         if (field.eq.'PROPOSID') INFO(9) = stream

         if (field.eq.'CRPIX1  ') HDR(01) = stream
         if (field.eq.'CRPIX2  ') HDR(02) = stream
         if (field.eq.'CRVAL1  ') HDR(03) = stream
         if (field.eq.'CRVAL2  ') HDR(04) = stream
         if (field.eq.'CTYPE1  ') HDR(05) = stream
         if (field.eq.'CTYPE2  ') HDR(06) = stream
         if (field.eq.'CD1_1   ') HDR(07) = stream
         if (field.eq.'CD1_2   ') HDR(08) = stream
         if (field.eq.'CD2_1   ') HDR(09) = stream
         if (field.eq.'CD2_2   ') HDR(10) = stream
         if (field.eq.'ORIENTAT') HDR(11) = stream
         if (field.eq.'PA_APER ') HDR(12) = stream
         if (field.eq.'PA_V3   ') HDR(13) = stream
         if (field.eq.'DATE-OBS') HDR(14) = stream
         if (field.eq.'TIME-OBS') HDR(15) = stream
         if (field.eq.'EXPTIME ') HDR(16) = stream
         if (field.eq.'ROOTNAME') HDR(17) = stream
         if (field.eq.'TARGNAME') HDR(18) = stream
         if (field.eq.'RA_TARG ') HDR(19) = stream
         if (field.eq.'DEC_TARG') HDR(20) = stream
         if (field.eq.'PROPOSID') HDR(21) = stream
         if (field.eq.'FILTER1 ') HDR(22) = stream
         if (field.eq.'FILTER2 ') HDR(23) = stream
         if (field.eq.'VAFACTOR') HDR(24) = stream

         if (field.eq.'END     ') goto 101
         enddo 
      goto 100
 101  continue

      nread = nread + 1

      if (extend_tf.and.NREAD.eq.1.and.NAXES.eq.0) goto 100 ! read the first extension...

      if (DIAG) then
         print*,'----------------------------------------'
         print*,'  NREAD: ',nread
         print*,'NEXTEND: ',nextend
         print*,'  NAXIS: ',naxes
         print*,'  LAXIS: ',laxis(1),laxis(2),laxis(3)
         print*,' BITPIX: ',bitpix
         print*,' BSCALE: ',bscale
         print*,'  BZERO: ',bzero
         endif

      ifirst = i+1
      i1 = i
      i2 = i

      nbper  = 4*laxis(1)*laxis(2)
      npt    = laxis(1)*laxis(2)
      nbyte1 = 1 
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880

      if (BITPIX.ne.-32) then
         print*,'   FILEU: ',FILEU(1:40)
         print*,'   FILE : ',FILE(1:40)
         print*,'BITPIX.ne.-32... unreal!'
         stop
         endif

      do i = 01, NDIMX
      do j = 01, NDIMY
         pix(i,j) = 0.
         enddo
         enddo

      do i = i1, i2, 1
         read(10,rec=i,iostat=ios) buffc
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/4 + 1
         np2 = (nbyteE-nbyte1)/4 + 1
         call buff2pix_r4(buffb,buffr,0001,0720)
         do ii = 001, 720
            nn = np1 + (ii-1)
            ny = 1 + (nn-1)/laxis(1)
            nx = nn-(ny-1)*laxis(1)
            if (nx.ge.001.and.nx.le.NDIMX.and.
     .          ny.ge.001.and.ny.le.NDIMY) then
                pix(nx,ny) = buffr(ii)
                nxx = nx
                nyy = ny
                endif
            enddo
         if (ny.gt.NDIMY) goto 899

c        if (DIAG) write(*,1115) i,np1,np2,npt,nxx,nyy
 1115    format(1x,i8,1x,i10,1x,i10,1x,i10,1x,2i6)
         enddo 

  899 close(10)
c     if (DIAG) write(*,1115) i,np1,np2,npt,nxx,nyy



      return

  900 continue
      print*,'                '
      print*,'READFITS ERROR: '
      print*,'                '
      write(*,'(''   could not read in file:   '',80a)') FILEU
      print*,'                '
      stop

      end

      subroutine buff2pix_r4(buff,pix,n1,nt)
      implicit none
      byte buff(2880)
      real pix(*)
      integer n1,nt

      byte b(4)
      real r
      equivalence(r,b)

      integer i, npu, nbu
      logical islinux

      do i = 1, 720
         npu = n1+i-1
         nbu = (i-1)*4
         if (.not.islinux()) then
            b(1) = buff(nbu+1)
            b(2) = buff(nbu+2)
            b(3) = buff(nbu+3)
            b(4) = buff(nbu+4)
            endif
         if (islinux()) then
            b(4) = buff(nbu+1)
            b(3) = buff(nbu+2)
            b(2) = buff(nbu+3)
            b(1) = buff(nbu+4)
            endif
         if (npu.ge.1.and.npu.le.nt) pix(npu) = r
         enddo

      return
      end

 
      logical function islinux()
      implicit none

      logical done
      logical islinux_save

      integer i
      data i/0/
      data islinux_save/.true./

      common /islinux_/i,islinux_save

      byte    b(2)
      equivalence(i,b)

      if (i.eq.1) goto 1
 
      i = 1
      islinux_save = .false.
      if (b(1).eq.1) islinux_save = .true.
 
    1 islinux = islinux_save
      return
      end

 
c--------------------------------------
c
c this will evaluate a PSF at a given (dx,dy)
c offset; for the regions within 4 pixels of
c the center, it uses bi-cubic interpolation
c
      real function rpsf_phot(x,y,psf)
      implicit none
      real x, y
      real psf(101,101)
 
 
      real    rx, ry
      integer ix, iy   !     3   4
      real    fx, fy   !    *1*  2
 
      real    dd
 
      real A1, B1, C1, D1, E1, F1, V1
      real A2, B2, C2, D2, E2, F2, V2
      real A3, B3, C3, D3, E3, F3, V3
      real A4, B4, C4, D4, E4, F4, V4
 
 
      rx = 51 + x*4
      ry = 51 + y*4
      ix = int(rx)
      iy = int(ry)
      fx = rx-ix
      fy = ry-iy
 
      dd = sqrt(x**2+y**2)
 
 
      rpsf_phot = 0.
      if (dd.gt.12.0) return
 
      if (dd.gt.4.0) then
         rpsf_phot = (1-fx)*(1-fy)*psf(ix  ,iy  )
     .             + ( fx )*(1-fy)*psf(ix+1,iy  )
     .             + (1-fx)*( fy )*psf(ix  ,iy+1)
     .             + ( fx )*( fy )*psf(ix+1,iy+1)
         return
         endif
 
 
      A1 =  psf(ix  ,iy  )
      B1 = (psf(ix+1,iy  )-psf(ix-1,iy  ))/2
      C1 = (psf(ix  ,iy+1)-psf(ix  ,iy-1))/2
      D1 = (psf(ix+1,iy  )+psf(ix-1,iy  )-2*A1)/2
      F1 = (psf(ix  ,iy+1)+psf(ix  ,iy-1)-2*A1)/2
      E1 = (psf(ix+1,iy+1)-A1)
 
      A2 =  psf(ix+1,iy  )
      B2 = (psf(ix+2,iy  )-psf(ix  ,iy  ))/2
      C2 = (psf(ix+1,iy+1)-psf(ix+1,iy-1))/2
      D2 = (psf(ix+2,iy  )+psf(ix  ,iy  )-2*A2)/2
      F2 = (psf(ix+1,iy+1)+psf(ix+1,iy-1)-2*A2)/2
      E2 =-(psf(ix  ,iy+1)-A2)
 
      A3 =  psf(ix  ,iy+1)
      B3 = (psf(ix+1,iy+1)-psf(ix-1,iy+1))/2
      C3 = (psf(ix  ,iy+2)-psf(ix  ,iy  ))/2
      D3 = (psf(ix+1,iy+1)+psf(ix-1,iy+1)-2*A3)/2
      F3 = (psf(ix  ,iy+2)+psf(ix  ,iy  )-2*A3)/2
      E3 =-(psf(ix+1,iy  )-A3)
 
      A4 =  psf(ix+1,iy+1)
      B4 = (psf(ix+2,iy+1)-psf(ix  ,iy+1))/2
      C4 = (psf(ix+1,iy+2)-psf(ix+1,iy  ))/2
      D4 = (psf(ix+2,iy+1)+psf(ix  ,iy+1)-2*A4)/2
      F4 = (psf(ix+1,iy+2)+psf(ix+1,iy  )-2*A4)/2
      E4 = (psf(ix  ,iy  )-A4)
 
 
      V1 = A1
     .   + B1*( fx )
     .   + C1*( fy )
     .   + D1*( fx )**2
     .   + E1*( fx )*( fy )
     .   + F1*( fy )**2
 
      V2 = A2
     .   + B2*(fx-1)
     .   + C2*( fy )
     .   + D2*(fx-1)**2
     .   + E2*(fx-1)*( fy )
     .   + F2*( fy )**2
 
      V3 = A3
     .   + B3*( fx )
     .   + C3*(fy-1)
     .   + D3*( fx )**2
     .   + E3*( fx )*(fy-1)
     .   + F3*(fy-1)**2
 
      V4 = A4
     .   + B4*(fx-1)
     .   + C4*(fy-1)
     .   + D4*(fx-1)**2
     .   + E4*(fx-1)*(fy-1)
     .   + F4*(fy-1)**2
 
      rpsf_phot = (1-fx)*(1-fy)*V1
     .          + ( fx )*(1-fy)*V2
     .          + (1-fx)*( fy )*V3
     .          + ( fx )*( fy )*V4
 
 
      return
      end




c------------------------------------------
c
c
c 
      subroutine mkimg(xu,yu,du,nu,NIMs,Us,pixu)
      implicit none

      integer Us
      real    xu(Us)
      real    yu(Us)
      real    du(Us)
      integer nu(Us)
      integer NIMs
      real    pixu(1001,1001)

      integer U
      integer i, j
      real    x, y
      real    d, dbar, dsig
      integer Ls, Lu
      real    dmin(99), dl(99)
      integer NIM, NIMu

      print*,'Us: ',Us
      do U = 1, Us
         print*,'U: ',U,xu(U),yu(U),du(U),nu(U)
         enddo

      do i = 0001, 1001
      do j = 0001, 1001
         x = (i-501)*0.01 
         y = (j-501)*0.01 
         Ls = 0
         do NIM = 1, NIMs
            dmin(NIM) = 9e9
            enddo
         do U = 1, Us
            d = sqrt((xu(U)-x)**2+(yu(U)-y)**2)
            if (d.lt.dmin(nu(U))) then
               dmin(nu(U)) = d
               dl(nu(U)) = du(U)
               endif
            enddo
         call barsig(dl,NIMs,dbar,dsig,NIMu)
         pixu(i,j) = dbar
         enddo
         enddo

      return
      end

  
      subroutine psf_norm(psfa,psfb)
      implicit none

      real psfa(101,101)
      real psfb(101,101)

      real*8 psftot

      integer ix, iy 
      integer ipsf, jpsf

      real strobe(5,5), strobeb
      real strobc(5,5), strobec

      real psfc(101,101)
      real xpsf, ypsf, rpsf, ff

      do ipsf = 001, 101
      do jpsf = 001, 101
         psfc(ipsf,jpsf) = 0.0
         do ix = -2, 2
         do iy = -2, 2
            ff = 1.00/16.00
            if (abs(ix).eq.2) ff = ff*0.5
            if (abs(iy).eq.2) ff = ff*0.5
            psfc(ipsf,jpsf) = psfc(ipsf,jpsf) + ff*psfa(ipsf+ix,jpsf+iy)
            enddo
            enddo
         enddo
         enddo

      do ix = 1, 4
      do iy = 1, 4
         strobe(ix,iy) = 0.0
         strobc(ix,iy) = 0.0
         enddo
         enddo
      do ipsf = 001, 101
      do jpsf = 001, 101
         ix = ipsf-(ipsf-1)/4*4
         iy = jpsf-(jpsf-1)/4*4
         xpsf = (ipsf-51)*0.25
         ypsf = (jpsf-51)*0.25
         rpsf = sqrt(xpsf**2+ypsf**2)
         if (rpsf.lt.5.00) then 
            strobe(ix,iy) = strobe(ix,iy) + psfa(ipsf,jpsf)
            strobc(ix,iy) = strobc(ix,iy) + psfc(ipsf,jpsf)
            endif
         enddo
         enddo

      if (.true.) then 
         do ix = 1, 4
         do iy = 1, 4
            write(*,'(2i2,1x,3f10.6)') ix,iy,
     .                                 strobe(ix,iy), 
     .                                 strobc(ix,iy),
     .                                 strobe(ix,iy)/strobc(ix,iy)
            strobe(ix,iy) = strobe(ix,iy)/strobc(ix,iy)
            enddo
            enddo
         endif

      strobeb = 0.
      do ix = 1, 4
      do iy = 1, 4
         strobeb = strobeb + strobe(ix,iy)/16.00
         enddo
         enddo

      print*,'          '
      print*,'PSF_NORM: ',strobeb
      write(*,'(10x,i1,5x,5f10.6,5x,5f10.6)') 5,
     .         (strobe(ix,1)             ,ix=1,4),
     .          strobe( 1,1),
     .         (strobe(ix,1)/strobeb-1.00,ix=1,4),
     .          strobe( 1,1)/strobeb-1.00
      do iy = 4, 1, -1
         write(*,'(10x,i1,5x,5f10.6,5x,5f10.6)') iy,
     .            (strobe(ix,iy)             ,ix=1,4),
     .             strobe( 1,iy),
     .            (strobe(ix,iy)/strobeb-1.00,ix=1,4),
     .             strobe( 1,iy)/strobeb-1.00
         enddo 
      print*,'          '

      do ipsf = 001, 101
      do jpsf = 001, 101
         ix = ipsf-(ipsf-1)/4*4
         iy = jpsf-(jpsf-1)/4*4
         psfb(ipsf,jpsf) = psfa(ipsf,jpsf)/strobe(ix,iy)
         enddo 
         enddo 

      return
      end


c-------------------------------------
c
c
      subroutine psf_smoo(mi,mo)
      implicit none

      real mi(101,101) 
      real mo(101,101) 

      real m0(101,101)
      real ma(101,101) 
      real mb(101,101) 
      real mc(101,101) 
      real md(101,101) 

      integer i, j
      real    r

      do i = 1, 101
      do j = 1, 101
         m0(i,j) = mi(i,j)
         enddo
         enddo

      call sm5quar(mi,ma,101,101)
      call sm5quad(mi,mb,101,101)
      call sm5plan(mi,mc,101,101)
      call sm7plan(mi,md,101,101)

      do i = 1, 101
      do j = 1, 101
         r = sqrt((i-51.)**2+(j-51.)**2)/4.0
         mo(i,j) = md(i,j)
         if (r.lt.08) mo(i,j) = (mc(i,j)+md(i,j))/2
         if (r.lt.07) mo(i,j) =  mc(i,j)
         if (r.lt.06) mo(i,j) = (mb(i,j)+mc(i,j))/2
         if (r.lt.05) mo(i,j) =  mb(i,j)
         if (r.lt.03) mo(i,j) = (ma(i,j)+mb(i,j))/2
         if (r.lt.02) mo(i,j) =  ma(i,j)
         enddo
         enddo

      return
      end





c------------------------------------
c
c
      subroutine sm7plan(r,s,NX,NY)
      implicit none

      integer NX, NY
      real r(NX,NY)
      real s(NX,NY)

      integer i, j
      integer iu, ju
      integer im, jm

      real dx, dy
      real A, B, C

      real AR,BR,CR
      real AA,BB,CC

      do i = 1, NX
      do j = 1, NY
         iu = max(min(i,NX-3),1+3)
         ju = max(min(j,NY-3),1+3)
         AR = 0.
         BR = 0.
         CR = 0.
         AA = 0.
         BB = 0.
         CC = 0.
         do im = 1, 7
         do jm = 1, 7
            AR = AR + r(iu-4+im,ju-4+jm)
            BR = BR + r(iu-4+im,ju-4+jm)*(im-4)
            CR = CR + r(iu-4+im,ju-4+jm)*(jm-4)
            AA = AA + 1
            BB = BB + (im-4)*(im-4)
            CC = CC + (jm-4)*(jm-4)
            enddo 
            enddo 
         A = AR/AA 
         B = BR/BB 
         C = CR/CC
         dx = i-iu 
         dy = j-ju
         s(i,j) = A + B*dx + C*dy
         enddo 
         enddo 

      return
      end



c------------------------------------
c
c
      subroutine sm5plan(r,s,NX,NY)
      implicit none

      integer NX, NY
      real r(NX,NY)
      real s(NX,NY)

      integer i, j
      integer iu, ju
      integer im, jm

      real dx, dy
      real A, B, C

      real AA( 5, 5)
      data AA / ! SUM:     16.0
     .      0.2500,   0.5000,   0.5000,   0.5000,   0.2500,
     .      0.5000,   1.0000,   1.0000,   1.0000,   0.5000,
     .      0.5000,   1.0000,   1.0000,   1.0000,   0.5000,
     .      0.5000,   1.0000,   1.0000,   1.0000,   0.5000,
     .      0.2500,   0.5000,   0.5000,   0.5000,   0.2500/

      real BB( 5, 5)
      data BB / ! SUM:     0.0
     .     -1.0000,  -0.5000,   0.0000,   0.5000,   1.0000,
     .     -1.0000,  -0.5000,   0.0000,   0.5000,   1.0000,
     .     -1.0000,  -0.5000,   0.0000,   0.5000,   1.0000,
     .     -1.0000,  -0.5000,   0.0000,   0.5000,   1.0000,
     .     -1.0000,  -0.5000,   0.0000,   0.5000,   1.0000/

      real CC( 5, 5)
      data CC / ! SUM:     0.0
     .     -1.0000,  -1.0000,  -1.0000,  -1.0000,  -1.0000,
     .     -0.5000,  -0.5000,  -0.5000,  -0.5000,  -0.5000,
     .      0.0000,   0.0000,   0.0000,   0.0000,   0.0000,
     .      0.5000,   0.5000,   0.5000,   0.5000,   0.5000,
     .      1.0000,   1.0000,   1.0000,   1.0000,   1.0000/


      do i = 1, NX
      do j = 1, NY
         iu = max(min(i,NX-2),3)
         ju = max(min(j,NY-2),3)
         A = 0.
         B = 0.
         C = 0.
         do im = 1, 5
         do jm = 1, 5
            A = A + AA(im,jm)*r(iu-3+im,ju-3+jm)/16
            B = B + BB(im,jm)*r(iu-3+im,ju-3+jm)/25
            C = C + CC(im,jm)*r(iu-3+im,ju-3+jm)/25
            enddo 
            enddo 
         dx = i-iu 
         dy = j-ju
         s(i,j) = A + B*dx + C*dy
         enddo 
         enddo 

      return
      end



c------------------------------------
c
c
      subroutine sm5quar(r,s,NX,NY)
      implicit none

      integer NX, NY
      real r(NX,NY)
      real s(NX,NY)

      integer i, j
      integer iu, ju
      integer im, jm

      real dx, dy
      real A

      real AA( 5, 5)
      data AA / ! SUM:     1.0
     .      1.0408,  -2.0204,   1.9592,  -2.0204,   1.0408,
     .     -2.0204,  -0.4898,   5.0204,  -0.4898,  -2.0204,
     .      1.9592,   5.0204,  11.0408,   5.0204,   1.9592,
     .     -2.0204,  -0.4898,   5.0204,  -0.4898,  -2.0204,
     .      1.0408,  -2.0204,   1.9592,  -2.0204,   1.0408/

      do i = 1, NX
      do j = 1, NY
         iu = max(min(i,NX-2),3)
         ju = max(min(j,NY-2),3)
         A = 0.
         do im = 1, 5
         do jm = 1, 5
            A = A + AA(im,jm)*r(iu-3+im,ju-3+jm)/25
            enddo 
            enddo 
         dx = i-iu 
         dy = j-ju
         s(i,j) = A 
         enddo 
         enddo 

      return
      end



c------------------------------------
c
c
      subroutine sm5quad(r,s,NX,NY)
      implicit none

      integer NX, NY
      real r(NX,NY)
      real s(NX,NY)

      integer i, j
      integer iu, ju
      integer im, jm

      real dx, dy
      real A, B, C, D, E, F

      real AA( 5, 5)
      data AA / ! SUM:     1.0
     .     -1.8571,   0.2857,   1.0000,   0.2857,  -1.8571,
     .      0.2857,   2.4286,   3.1429,   2.4286,   0.2857,
     .      1.0000,   3.1429,   3.8571,   3.1429,   1.0000,
     .      0.2857,   2.4286,   3.1429,   2.4286,   0.2857,
     .     -1.8571,   0.2857,   1.0000,   0.2857,  -1.8571/

      real BB( 5, 5)
      data BB / ! SUM:     0.0
     .     -1.0000,  -0.5000,   0.0000,   0.5000,   1.0000,
     .     -1.0000,  -0.5000,   0.0000,   0.5000,   1.0000,
     .     -1.0000,  -0.5000,   0.0000,   0.5000,   1.0000,
     .     -1.0000,  -0.5000,   0.0000,   0.5000,   1.0000,
     .     -1.0000,  -0.5000,   0.0000,   0.5000,   1.0000/

      real CC( 5, 5)
      data CC / ! SUM:     0.0
     .     -1.0000,  -1.0000,  -1.0000,  -1.0000,  -1.0000,
     .     -0.5000,  -0.5000,  -0.5000,  -0.5000,  -0.5000,
     .      0.0000,   0.0000,   0.0000,   0.0000,   0.0000,
     .      0.5000,   0.5000,   0.5000,   0.5000,   0.5000,
     .      1.0000,   1.0000,   1.0000,   1.0000,   1.0000/

      real DD( 5, 5)
      data DD / ! SUM:     0.0
     .      0.7143,  -0.3571,  -0.7143,  -0.3571,   0.7143,
     .      0.7143,  -0.3571,  -0.7143,  -0.3571,   0.7143,
     .      0.7143,  -0.3571,  -0.7143,  -0.3571,   0.7143,
     .      0.7143,  -0.3571,  -0.7143,  -0.3571,   0.7143,
     .      0.7143,  -0.3571,  -0.7143,  -0.3571,   0.7143/

      real EE( 5, 5)
      data EE / ! SUM:     0.0
     .      1.0000,   0.5000,   0.0000,  -0.5000,  -1.0000,
     .      0.5000,   0.2500,   0.0000,  -0.2500,  -0.5000,
     .      0.0000,   0.0000,   0.0000,   0.0000,   0.0000,
     .     -0.5000,  -0.2500,   0.0000,   0.2500,   0.5000,
     .     -1.0000,  -0.5000,   0.0000,   0.5000,   1.0000/

      real FF( 5, 5)
      data FF / ! SUM:     0.0
     .      0.7143,   0.7143,   0.7143,   0.7143,   0.7143,
     .     -0.3571,  -0.3571,  -0.3571,  -0.3571,  -0.3571,
     .     -0.7143,  -0.7143,  -0.7143,  -0.7143,  -0.7143,
     .     -0.3571,  -0.3571,  -0.3571,  -0.3571,  -0.3571,
     .      0.7143,   0.7143,   0.7143,   0.7143,   0.7143/


      do i = 1, NX
      do j = 1, NY
         iu = max(min(i,NX-2),3)
         ju = max(min(j,NY-2),3)
         A = 0.
         B = 0.
         C = 0.
         D = 0.
         E = 0.
         F = 0.
         do im = 1, 5
         do jm = 1, 5
            A = A + AA(im,jm)*r(iu-3+im,ju-3+jm)/25
            B = B + BB(im,jm)*r(iu-3+im,ju-3+jm)/25
            C = C + CC(im,jm)*r(iu-3+im,ju-3+jm)/25
            D = D + DD(im,jm)*r(iu-3+im,ju-3+jm)/25
            E = E + EE(im,jm)*r(iu-3+im,ju-3+jm)/25
            F = F + FF(im,jm)*r(iu-3+im,ju-3+jm)/25
            enddo 
            enddo 
         dx = i-iu 
         dy = j-ju
         s(i,j) = A + 
     .            B*dx    + C*dy    + 
     .            D*dx**2 + E*dx*dy + F*dy**2
         enddo 
         enddo 

      return
      end



c------------------------------------
c
c
      subroutine sm7quar(r,s,NX,NY)
      implicit none

      integer NX, NY
      real r(NX,NY)
      real s(NX,NY)

      integer i, j
      integer iu, ju
      integer im, jm

      real dx, dy
      real A

      real AA( 7, 7)
      data AA / !     1.0
     .  2.0808,-1.7576,-0.2424, 0.8990,-0.2424,-1.7576, 2.0808,
     . -1.7576,-2.8182, 0.3636, 2.0606, 0.3636,-2.8182,-1.7576,
     . -0.2424, 0.3636, 4.5455, 6.5758, 4.5455, 0.3636,-0.2424,
     .  0.8990, 2.0606, 6.5758, 8.7172, 6.5758, 2.0606, 0.8990,
     . -0.2424, 0.3636, 4.5455, 6.5758, 4.5455, 0.3636,-0.2424,
     . -1.7576,-2.8182, 0.3636, 2.0606, 0.3636,-2.8182,-1.7576,
     .  2.0808,-1.7576,-0.2424, 0.8990,-0.2424,-1.7576, 2.0808/
  
      do i = 1, NX
      do j = 1, NY
         A = 0
         if (i.ge.4.and.i.le.NX-3.and.j.ge.4.and.j.le.NX-3) then
            do im = 1, 7
            do jm = 1, 7
               A = A + AA(im,jm)*r(i-4+im,j-4+jm)/49
               enddo 
               enddo 
            endif
         s(i,j) = A 
         enddo 
         enddo 

      return
      end





c------------------------------------
c
c
      subroutine sm7quad(r,s,NX,NY)
      implicit none

      integer NX, NY
      real r(NX,NY)
      real s(NX,NY)

      integer i, j
      integer iu, ju
      integer im, jm

      real dx, dy
      real A

      real AA( 7, 7)
      data AA / ! SUM:     1.0
     . -2.3333,-0.6667, 0.3333, 0.6667, 0.3333,-0.6667,-2.3333,
     . -0.6667, 1.0000, 2.0000, 2.3333, 2.0000, 1.0000,-0.6667,
     .  0.3333, 2.0000, 3.0000, 3.3333, 3.0000, 2.0000, 0.3333,
     .  0.6667, 2.3333, 3.3333, 3.6667, 3.3333, 2.3333, 0.6667,
     .  0.3333, 2.0000, 3.0000, 3.3333, 3.0000, 2.0000, 0.3333,
     . -0.6667, 1.0000, 2.0000, 2.3333, 2.0000, 1.0000,-0.6667,
     . -2.3333,-0.6667, 0.3333, 0.6667, 0.3333,-0.6667,-2.3333/
  
      do i = 1, NX
      do j = 1, NY
         A = 0
         if (i.ge.4.and.i.le.NX-3.and.j.ge.4.and.j.le.NY-3) then
            do im = 1, 7
            do jm = 1, 7
               A = A + AA(im,jm)*r(i-4+im,j-4+jm)/49
               enddo 
               enddo 
            endif
         s(i,j) = A 
         enddo 
         enddo 

      return
      end



c-----------------------------------
c
c
      subroutine smoo_2d_9X(psf,pss)
      implicit none
      real psf(101,101,09,10)
      real pss(101,101,09,10)

      real a(9,10)
      real b(9,10)
      real c(9,10)

      real abar
      real renorm09x10
      integer ipsf, jpsf
      integer ir, jr
      real rpsf

      do ipsf = 001, 101
      do jpsf = 001, 101
         rpsf = sqrt((ipsf-51.)**2+(jpsf-51.)**2)/4
         do ir = 01, 09
         do jr = 01, 10
            a(ir,jr) = psf(ipsf,jpsf,ir,jr)
            enddo
            enddo
         call sm5quad(a,c,09,10)
         call sm3quad(a,b,09,10)
         do ir = 01, 09
         do jr = 01, 10
            pss(ipsf,jpsf,ir,jr) = a(ir,jr)
            if (ir.gt.01.and.ir.lt.09.and.
     .          jr.gt.01.and.jr.lt.10) then
                pss(ipsf,jpsf,ir,jr) =  a(ir,jr)
                if (rpsf.ge.3.00) pss(ipsf,jpsf,ir,jr) =  b(ir,jr)
                if (rpsf.ge.3.75) pss(ipsf,jpsf,ir,jr) = (b(ir,jr)+
     .                                                    c(ir,jr))/2
                if (rpsf.ge.4.50) pss(ipsf,jpsf,ir,jr) =  c(ir,jr)
                endif
            pss(ipsf,jpsf,ir,jr) = (psf(ipsf,jpsf,ir,jr)+
     .                              pss(ipsf,jpsf,ir,jr))/2
            enddo
            enddo
         if (ipsf.eq.51.and.jpsf.eq.51) then
            abar = 0.
            do ir = 1, 09
            do jr = 1, 10
               abar = abar + a(ir,jr)
               enddo
               enddo
            abar = abar/90
            write(*,*)
            write(*,'(10x,''smoo_2d_99: '',f8.6)') abar
            write(*,*)
            do jr = 10, 01, -1
               write(*,'(9x,i2,2x,9i6,2x,9i6,2x,9i6,2x,9i6)') jr,
     .            (int((psf(51,51,ir,jr))*1e5+0.5),ir=1,9),
     .            (int((pss(51,51,ir,jr)-
     .                  psf(51,51,ir,jr))*1e5+0.5),ir=1,9)
               enddo
            print*,'            '
            do jr = 10, 01, -1
               write(*,'(9x,i2,2x,9i6,2x,9i6,2x,9i6,2x,9i6)') jr,
     .            (int((c(ir,jr)-a(ir,jr))*1e5+0.5),ir=1,9),
     .            (int((b(ir,jr)-a(ir,jr))*1e5+0.5),ir=1,9)
               enddo
            print*,'            '
            do jr = 10, 01, -1
               write(*,'(9x,i2,2x,9i6,2x,9i6,2x,9i6,2x,9i6)') jr,
     .            (int((pss(51,51,ir,jr)     )*1e5+0.5),ir=1,9),
     .            (int((pss(51,51,ir,jr)-abar)*1e5+0.5),ir=1,9)
               enddo
            print*,'            '
            print*,'            '
            endif
         enddo
         enddo


      return
      end




c-----------------------------------
c
c
      subroutine sm3quad(r,s,NX,NY)
      implicit none

      integer NX, NY
      real r(NX,NY)
      real s(NX,NY)

      integer i, j
      integer iu, ju
      integer im, jm

      real dx, dy
      real A, B, C, D, E, F

      real AA( 3, 3)
      data AA / ! SUM:     1.0
     .     -1.0000,   2.0000,  -1.0000,
     .      2.0000,   5.0000,   2.0000,
     .     -1.0000,   2.0000,  -1.0000/

      real BB( 3, 3)
      data BB / ! SUM:     0.0
     .     -1.5000,   0.0000,   1.5000,
     .     -1.5000,   0.0000,   1.5000,
     .     -1.5000,   0.0000,   1.5000/

      real CC( 3, 3)
      data CC / ! SUM:     0.0
     .     -1.5000,  -1.5000,  -1.5000,
     .      0.0000,   0.0000,   0.0000,
     .      1.5000,   1.5000,   1.5000/

      real DD( 3, 3)
      data DD / ! SUM:     0.0
     .      1.5000,  -3.0000,   1.5000,
     .      1.5000,  -3.0000,   1.5000,
     .      1.5000,  -3.0000,   1.5000/

      real EE( 3, 3)
      data EE / ! SUM:     0.0
     .      2.2500,   0.0000,  -2.2500,
     .      0.0000,   0.0000,   0.0000,
     .     -2.2500,   0.0000,   2.2500/

      real FF( 3, 3)
      data FF / ! SUM:     0.0
     .      1.5000,   1.5000,   1.5000,
     .     -3.0000,  -3.0000,  -3.0000,
     .      1.5000,   1.5000,   1.5000/


      do i = 1, NX
      do j = 1, NY
         iu = max(min(i,NX-1),2)
         ju = max(min(j,NY-1),2)
         A = 0.
         B = 0.
         C = 0.
         D = 0.
         E = 0.
         F = 0.
         do im = 1, 3
         do jm = 1, 3
            A = A + AA(im,jm)*r(iu-2+im,ju-2+jm)/09
            B = B + BB(im,jm)*r(iu-2+im,ju-2+jm)/09
            C = C + CC(im,jm)*r(iu-2+im,ju-2+jm)/09
            D = D + DD(im,jm)*r(iu-2+im,ju-2+jm)/09
            E = E + EE(im,jm)*r(iu-2+im,ju-2+jm)/09
            F = F + FF(im,jm)*r(iu-2+im,ju-2+jm)/09
            enddo 
            enddo 
         dx = i-iu 
         dy = j-ju
         s(i,j) = A + 
     .            B*dx    + C*dy    + 
     .            D*dx**2 + E*dx*dy + F*dy**2
         enddo 
         enddo 

      return
      end



c-------------------------------------
c
c
      subroutine smoo_psfnim(mi,mo)
      implicit none

      real mi(101,101) 
      real mo(101,101) 

      real m0(101,101)
      real ma(101,101) 
      real mb(101,101) 

      integer i, j
      real    r

      do i = 1, 101
      do j = 1, 101
         m0(i,j) = mi(i,j)
         r = sqrt((i-51.)**2+(j-51.)**2)/4.0
         if (r.gt.5) m0(i,j) = 0. 
         enddo
         enddo

      call sm5quad(m0,ma,101,101)
      call sm7plan(m0,mb,101,101)

      do i = 1, 101
      do j = 1, 101
         r = sqrt((i-51.)**2+(j-51.)**2)/4.0
         mo(i,j) = mb(i,j)
         if (r.lt.04.0) mo(i,j) = (ma(i,j)+mb(i,j))/2
         if (r.lt.03.5) mo(i,j) =  ma(i,j)
         enddo
         enddo

c      print*,'---> SMOO_PSFNIM: ',mi(51,51),mo(51,51),
c     .                            ma(51,51),mb(51,51)

      return
      end

      real function find_error(pu,fu,Us) 
      implicit none

      integer Us
      real pu(Us)
      real fu(Us)

      integer U

      real*8 ftot, ptot
      real*8 z, ee 

      ftot = 0.
      ptot = 0.
      do U = 1, Us
         ftot = ftot + fu(u)
         ptot = ptot + pu(u)
         enddo

      z = ptot/ftot
 
      ee = 0.0
      do U = 1, Us
         ee = ee + abs(pu(u)-z*fu(u))/z
         enddo

      find_error = ee

      return
      end


c--------------------------------------------------
c
c determine whether any images are "much" different
c from the others for a given star... if so, reject
c that image!
c
      subroutine zm_anal(z_m,sz_m,u_m,zbar,zsig,Ms,Mu)
      implicit none

      integer Ms, Mu
      real    z_m(Ms)
      integer u_m(Ms)   
      real    zbar
      real    zsig
 
      real ztot, utot
      integer M, MM, MMAX
      real sz_m(_MDIM_)

      do M = 1, Ms
         u_m(M) = 1
         enddo

    1 MMAX = 0
      do M = 1, Ms
         if (u_m(M).eq.1) then
            ztot = 0.
            utot = 0
            do MM = 1, Ms
               if (M.ne.MM) then
                  ztot = ztot + u_m(MM)*z_m(MM)
                  utot = utot + u_m(MM)
                  endif
               enddo
            if (utot.le.2) stop 'utot.le.2'
            zbar = ztot/utot
            ztot = 0.
            do MM = 1, Ms
               if (M.ne.MM) then
                  ztot = ztot + u_m(MM)*(z_m(MM)-zbar)**2
                  endif
               enddo
            zsig = sqrt(ztot/(utot-1))
            sz_m(M) = abs(z_m(M)-zbar)/zsig
            if (MMAX.eq.0) MMAX = M
            if (sz_m(M).gt.sz_m(MMAX)) MMAX = M
            endif
         enddo

      if (sz_m(MMAX).gt.4.0) then
         u_m(MMAX) = 0
         if (utot.ge.4) goto 1
         endif

      ztot = 0.
      utot = 0
      do M = 1, Ms
         ztot = ztot + u_m(M)*z_m(M)
         utot = utot + u_m(M)
         enddo
      zbar = ztot/utot
      ztot = 0.
      do M = 1, Ms
         ztot = ztot + u_m(M)*(z_m(M)-zbar)**2
         enddo
      zsig = sqrt(ztot/utot)
      Mu = utot

      return
      end


c-----------------------------------------------------
c
c this just writes a real*4 fits image
c

      subroutine writfits_r4(FILE,pix,PXDIMX,PXDIMY)
      implicit none

      character*(*) FILE
      integer PXDIMX,PXDIMY
      real    pix(PXDIMX,PXDIMY)

      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer i,ios

      character*2880 buffc
      byte buffb(2880)
      equivalence (buffb,buffc)

      integer ifirst, i1, i2

      integer np1, np2, npt
      integer k

      character*80 FILEU
      character*70 HDR(25)
      common/HDR/HDR

      FILEU = FILE
      do i = 1, 195
         if (FILE(i:i+4).eq.'.fits') then
             FILEU = FILE(1:i+4)
             goto 1
             endif
         enddo
   
      print*,' '
      print*,'readfits_r4:    no .fits'
      print*,'      FILE:  ',FILE
      stop
    1 continue

      i = 1

      open(10,file=FILEU,status='unknown',
     .     err=900,recl=2880,form='UNFORMATTED',
     .     access='DIRECT')

      do k = 00, 34
         write(buffc(k*80+1:k*80+80),'(80x)') 
         enddo

      write(10,rec=i,iostat=ios) buffc

      write(buffc( 0*80+1: 1*80),'(''SIMPLE  =                    T'')')
      write(buffc( 1*80+1: 2*80),'(''BITPIX  =                  -32'')')
      write(buffc( 2*80+1: 3*80),'(''NAXIS   = '',8x,i12)') 2
      write(buffc( 3*80+1: 4*80),'(''NAXIS1  = '',8x,i12)') PXDIMX
      write(buffc( 4*80+1: 5*80),'(''NAXIS2  = '',8x,i12)') PXDIMY
      write(buffc( 5*80+1: 6*80),'(''DATATYPE= '',9a)') 
     .                          " 'REAL*4' "
      write(buffc(07*80+1:08*80),'(''COMMENT   '',a05)') '     '
      write(buffc(08*80+1:09*80),'(''COMMENT   '',a05)') '     '
      write(buffc(09*80+1:10*80),'(''CRPIX1  = '',a70)') HDR(01)
      write(buffc(10*80+1:11*80),'(''CRPIX2  = '',a70)') HDR(02)
      write(buffc(11*80+1:12*80),'(''CRVAL1  = '',a70)') HDR(03)
      write(buffc(12*80+1:13*80),'(''CRVAL2  = '',a70)') HDR(04)
      write(buffc(13*80+1:14*80),'(''CTYPE1  = '',a70)') HDR(05)
      write(buffc(14*80+1:15*80),'(''CTYPE2  = '',a70)') HDR(06)
      write(buffc(15*80+1:16*80),'(''CD1_1   = '',a70)') HDR(07)
      write(buffc(16*80+1:17*80),'(''CD1_2   = '',a70)') HDR(08)
      write(buffc(17*80+1:18*80),'(''CD2_1   = '',a70)') HDR(09)
      write(buffc(18*80+1:19*80),'(''CD2_2   = '',a70)') HDR(10)
      write(buffc(19*80+1:20*80),'(''ORIENTAT= '',a70)') HDR(11)
      write(buffc(20*80+1:21*80),'(''PA_APER = '',a70)') HDR(12)
      write(buffc(21*80+1:22*80),'(''PA_V3   = '',a70)') HDR(13)
      write(buffc(22*80+1:23*80),'(''DATE-OBS= '',a70)') HDR(14)
      write(buffc(23*80+1:24*80),'(''TIME-OBS= '',a70)') HDR(15)
      write(buffc(24*80+1:25*80),'(''EXPTIME = '',a70)') HDR(16)
      write(buffc(25*80+1:26*80),'(''ROOTNAME= '',a70)') HDR(17)
      write(buffc(26*80+1:27*80),'(''TARGNAME= '',a70)') HDR(18)
      write(buffc(27*80+1:28*80),'(''RA_TARG = '',a70)') HDR(19)
      write(buffc(28*80+1:29*80),'(''DEC_TARG= '',a70)') HDR(20)
      write(buffc(29*80+1:30*80),'(''PROPOSID= '',a70)') HDR(21)
      write(buffc(30*80+1:31*80),'(''FILTER1 = '',a70)') HDR(22)
      write(buffc(31*80+1:32*80),'(''FILTER2 = '',a70)') HDR(23)
      write(buffc(33*80+1:34*80),'(''VAFACTOR= '',a70)') HDR(24)
      write(buffc(32*80+1:33*80),'(''CCDGAIN = '',a70)') HDR(25)
      write(buffc(33*80+1:34*80),'(''COMMENT   '',a05)') '     '
      write(buffc(34*80+1:35*80),'(''COMMENT   '',a05)') '     '
      write(buffc(35*80+1:36*80),'(''END       '')')

      write(10,rec=i,iostat=ios) buffc


      ifirst = i+1
      i1 = i
      i2 = i

      nbper  = 4*PXDIMX*PXDIMY
      npt    =   PXDIMX*PXDIMY
      nbyte1 = 1 
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880

      do i = i1, i2, 1
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/4 + 1
         np2 = (nbyteE-nbyte1)/4 + 1
         call pix2buff_r4(buffb,pix,np1,npt)
         write(10,rec=i,iostat=ios) buffc
         enddo 

      close(10)

      return

 900  continue
      print*,'writfits_r4.f ERROR'
      print*,'   FILEU: ',FILEU
      stop

      end

c-------------------------------------------------------
c
c
      subroutine pix2buff_r4(buff,pix,n1,nt)
      implicit none
      byte buff(2880)
      real*4 pix(*)
      integer n1,nt

      byte b(4)
      real*4 r
      equivalence(r,b)

      integer i, npu, nbu
      logical islinux

      do i = 1, 720
         npu = n1+i-1
         nbu = (i-1)*4
         if (npu.ge.1.and.npu.le.nt) r = pix(npu)
         if (.not.islinux()) then
            buff(nbu+1) = b(1)
            buff(nbu+2) = b(2)
            buff(nbu+3) = b(3)
            buff(nbu+4) = b(4)
            endif
         if (islinux()) then
            buff(nbu+1) = b(4)
            buff(nbu+2) = b(3)
            buff(nbu+3) = b(2)
            buff(nbu+4) = b(1)
            endif
         enddo

      return
      end

c-------------------------------------------------------------
c RANDOM NUMBER GENERATOR
c----------------------------------------------------
      FUNCTION RN()
      IMPLICIT NONE
      INTEGER*4 ISEED
      LOGICAL INIT
      DOUBLE PRECISION    RN
      DOUBLE PRECISION    DS(2),    DM(2),    DSEED
      DOUBLE PRECISION    DX24,     DX48
      DOUBLE PRECISION    DL,       DC,       DU,       DR
      COMMON/RN_COM/ISEED
      DATA      DS     /  1665 1885., 286 8876.  /
      DATA      DM     /  1518 4245., 265 1554.  /
      DATA      DX24   /  1677 7216.  /
      DATA      DX48   /  281 4749 7671 0656.  /
      DATA INIT/.TRUE./
      IF(INIT)THEN
        INIT=.FALSE.
        DSEED=DFLOAT(ISEED)
        DS(2)  =  DINT(DSEED/DX24)
        DS(1)  =  DSEED - DS(2)*DX24
      END IF
      DL  =  DS(1) * DM(1)
      DC  =  DINT(DL/DX24)
      DL  =  DL - DC*DX24
      DU  =  DS(1)*DM(2) + DS(2)*DM(1) + DC
      DS(2)  =  DU - DINT(DU/DX24)*DX24
      DS(1)  =  DL
      RN=  (DS(2)*DX24 + DS(1)) / DX48
      RETURN
      END
c-----------------------------------------------------------

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      integer function lenc(s)
c
c  Get the length of a character variable.  Only necessary because of the
c  Unix lobotomy of Fortran.
c
      character*(*) s
      lenc = len(s)
      do 900 i = len(s),1,-1
          if(s(i:i).ne.' ') then
              lenc = i
              return
          end if
900   continue
      lenc = 0
      return
      end

c==============================================================================
C $Id: svdcmp.F,v 1.4 1992/10/21 12:40:05 bennett Exp $ Initial revision
c Initial revision
c
      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V,ier)

ccc       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=500)
      DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)

      ier=0

      G=0.0
      SCALE=0.0
      ANORM=0.0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.0
        S=0.0
        SCALE=0.0
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+ABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=0.0
                DO 13 K=I,M
                  S=S+A(K,I)*A(K,J)
13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0
        S=0.0
        SCALE=0.0
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+ABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=A(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=0.0
                DO 21 K=L,N
                  S=S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,N
                  A(J,K)=A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.0.0) THEN
            DO 26 J=L,N
              V(J,I)=(A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=0.0
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.0
            V(J,I)=0.0
31        CONTINUE
        ENDIF
        V(I,I)=1.0
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            A(I,J)=0.0
33        CONTINUE
        ENDIF
        IF (G.NE.0.0) THEN
          G=1.0/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=0.0
              DO 34 K=L,M
                S=S+A(K,I)*A(K,J)
34            CONTINUE
              F=(S/A(I,I))*G
              DO 35 K=I,M
                A(K,J)=A(K,J)+F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=0.0
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+1.0
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=0.0
          S=1.0
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((ABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=SQRT(F*F+G*G)
              W(I)=H
              H=1.0/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.0) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.30) then
             write(6,*) 'No convergence in 30 iterations of SVDCMP'
             ier=1
             RETURN
          endif
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
          G=SQRT(F*F+1.0)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C=1.0
          S=1.0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=SQRT(F*F+H*H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=SQRT(F*F+H*H)
            W(J)=Z
            IF (Z.NE.0.0) THEN
              Z=1.0/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 NM=1,M
              Y=A(NM,J)
              Z=A(NM,I)
              A(NM,J)= (Y*C)+(Z*S)
              A(NM,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=0.0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END
c==============================================================================
C $Id: svbksb.F,v 1.1 1993/03/20 17:55:15 bennett Exp $ $Log: svbksb.F,v $
c
      SUBROUTINE SVBKSB(U,W,V,M,N,MP,NP,B,X)

ccc       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=500)
      DIMENSION U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP),TMP(NMAX)
      DO 12 J=1,N
        S=0.
        IF(W(J).NE.0.)THEN
          DO 11 I=1,M
            S=S+U(I,J)*B(I)
11        CONTINUE
          S=S/W(J)
        ENDIF
        TMP(J)=S
12    CONTINUE
      DO 14 J=1,N
        S=0.
        DO 13 JJ=1,N
          S=S+V(J,JJ)*TMP(JJ)
13      CONTINUE
        X(J)=S
14    CONTINUE
      RETURN
      END
