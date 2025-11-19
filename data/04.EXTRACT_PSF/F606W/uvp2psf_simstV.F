 
#define _LDIM_ 500000
#define _MDIM_     20
#define _NDIM_     40
#define _LINUX_ .true.

c#include "/home/jayander/FORTRAN/MAC_OSX/FITS/writfits_r4.f"

      program uvp2anal
      implicit none
  
      integer   L, Ls
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

      integer U, Us
      real*8 uu(_LDIM_), vu(_LDIM_), Uo, Vo
      real*8 iu(_LDIM_), ju(_LDIM_), Io, Jo
      real*8 i0(_NDIM_), j0(_NDIM_)
      integer NIM, NIMs 

      real*8 AGN(_NDIM_), BGN(_NDIM_), CGN(_NDIM_), DGN(_NDIM_)
      real*8 GAN(_NDIM_), GBN(_NDIM_), GCN(_NDIM_), GDN(_NDIM_)
      real*8 DLTA

      integer Ss, Su 
      real*4 ps(999), pbar, psig

      integer O, Os, Ou
      real*4 fo(999), fbar, fsig
      real*4 po(999)
      real*4 do(999), dbar, dsig

      real rpsf_photij
      real psf78u(101,101,7,8)
      common /psf78_/psf78u

      real SKYN(_NDIM_)

      integer idx, idy
      integer idx1,idy1
      integer idx2,idy2
      real    dx, dy

      real   dxu(_LDIM_), dxu1, dxu2
      real   dyu(_LDIM_), dyu1, dyu2
      real    fu(_LDIM_), fu1(_LDIM_), fu2(_LDIM_)
      real    pu(_LDIM_)
      real    eu(_LDIM_)
      integer nu(_LDIM_)
      real    xu(_LDIM_)
      real    yu(_LDIM_)
      real    du(_LDIM_)
      real*8  ic, jc, ic1, jc1, ic2, jc2
      real*8  uc, vc, uc1, vc1, uc2, vc2
      real    pixu(1001,1001)

      real*8  ptot, ftot
      real*8  z
      real*8  chisq, chisq_min
      integer if2
      real*8  f1, f2

      integer NST, NSTs
      real    sky_nm(_NDIM_,_MDIM_)

      real    rpsf_phot
      real    psfo(101,101)
      real    psfa(101,101)
      real    psfe(101,101)
      real    psfn(101,101)
      real    psfs(101,101)
      real    psfu(101,101)

      real    show_psf(501,1001)
      real    show_str(4001,2001)
      real    xpsf, ypsf, rpsf
      integer ipsf, jpsf
      real    z_n(_NDIM_)     
      real    z_mn(_MDIM_,_NDIM_), zbar, zsig      
      real   sz_mn(_MDIM_,_NDIM_)     
      integer u_mn(_MDIM_,_NDIM_)     
      integer NIMu
      real    e, find_error 
      real*8  emin, xmin, ymin
      integer i, j
      integer iis,igood,is,ngood_psf

      character*80 STRING
      integer NIT

      logical psfuse_nstar(_NDIM_)

      real x1, y1
      real x2, y2

c
c
c---------------------------------------------------------------
c
c
ccc      psfuse_nstar(01) = .false.
ccc      psfuse_nstar(02) = .true.
ccc      psfuse_nstar(03) = .true.
ccc      psfuse_nstar(04) = .true.
ccc      psfuse_nstar(05) = .true.
ccc      psfuse_nstar(06) = .true.
ccc      psfuse_nstar(07) = .true.
ccc      psfuse_nstar(08) = .true.
ccc      psfuse_nstar(09) = .true.
ccc      psfuse_nstar(10) = .true.
ccc      psfuse_nstar(11) = .true.
ccc      psfuse_nstar(12) = .true.
ccc      psfuse_nstar(13) = .true.
ccc      psfuse_nstar(14) = .true.
ccc      psfuse_nstar(15) = .true.
ccc      psfuse_nstar(16) = .true.
ccc      psfuse_nstar(17) = .true.
ccc      psfuse_nstar(18) = .true.
ccc      psfuse_nstar(19) = .true.
ccc      psfuse_nstar(20) = .true.
ccc      psfuse_nstar(21) = .true.
ccc      psfuse_nstar(22) = .true.
ccc      psfuse_nstar(23) = .true.
ccc      psfuse_nstar(24) = .true.
ccc      psfuse_nstar(25) = .true.
ccc      psfuse_nstar(26) = .true.
ccc      psfuse_nstar(27) = .true.
ccc      psfuse_nstar(28) = .true.
ccc      psfuse_nstar(29) = .true.
ccc      psfuse_nstar(30) = .true.
ccc      psfuse_nstar(31) = .true.
ccc      psfuse_nstar(32) = .true.
ccc      psfuse_nstar(33) = .true.
ccc      psfuse_nstar(34) = .true.
ccc      psfuse_nstar(35) = .true.
ccc      psfuse_nstar(36) = .true.
ccc      psfuse_nstar(37) = .true.
ccc      psfuse_nstar(38) = .true.
ccc      psfuse_nstar(39) = .true.
ccc      psfuse_nstar(40) = .true.

c
c
c---------------------------------------------------------------
c
c

      NIMs = 0
      NSTs = 0
      Ls = 0
    1 read(*,*,end=2) ul(Ls+1), vl(Ls+1), pl(Ls+1),
     .              niml(Ls+1), il(Ls+1), jl(Ls+1),
     .                xl(Ls+1), yl(Ls+1), 
     .               dul(Ls+1),dvl(Ls+1), rl(Ls+1),
     .               nstl(Ls+1) 
      Ls = Ls + 1
      if (Ls.gt._LDIM_) stop 'Ls.gt._LDIM_'
      if (niml(Ls).gt.NIMs) NIMs = niml(Ls)
      if (nstl(Ls).gt.NSTs) NSTs = nstl(Ls)
      if (NIMs.gt._MDIM_) stop 'NIMs.gt._MDIM_'
      if (NSTs.gt._NDIM_) stop 'NIMs.gt._NDIM_'

      goto 1
    2 continue

c     Read list of good PSFs
c     ----------------------
      open(unit=3,file='IN.good_psf_list',status='old')
      ngood_psf = 0
      do is = 1,NSTs
        read(3,*) iis,igood
        if(iis.ne.is) then
          stop 'error in IN.good_psf_list'
        elseif(igood.eq.1) then
          psfuse_nstar(is) = .true.
          ngood_psf = ngood_psf + 1
        else
          psfuse_nstar(is) = .false.
        endif
      enddo
      psfuse_nstar(01) = .false.

      print*,'      '
      print*,'NSTs: ',NSTs
      print*,'NIMs: ',NIMs
      print*,'  Ls: ',Ls
      print*,'Ngood: ',ngood_psf
      print*,'      '

    
      print*,'      ' 
      print*,'  Ls: ',Ls
      print*,'NIMs: ',NIMs 
      print*,'      ' 

      write(6,99) (i, psfuse_nstar(i), i = 1,NSTs)
 99   format(i4,L3)

      open(12,file='uvp2psf_simstV_sky.out',status='unknown')
      do NST = 1, NSTs
      do NIM = 1, NIMs
         Ss = 0
         do L = 1, Ls
            if (niml(L).eq.NIM.and.
     .          nstl(L).eq.NST.and.
     .          rl(L).gt.8.5.and.rl(L).lt.13.5) then
                Ss = Ss + 1
                ps(Ss) = pl(L)
                endif
            enddo
         pbar = 0.
         psig = 0.
         Su = 0
         if (Ss.gt.0) call barsig(ps,Ss,pbar,psig,Su)
         SKYN(NIM) = pbar
         write( *,123) NIM,NST,pbar,psig,Ss,Su
         write(12,123) NIM,NST,pbar,psig,Ss,Su
  123    format(1x,i2.2,3x,i2.2,5x,f8.3,1x,f8.3,1x,i3,1x,i3)
         sky_nm(NST,NIM) = pbar
         enddo
         enddo
      close(12)

      do L = 1, Ls
         NST = nstl(L)
         NIM = niml(L)
         ppl(L) = pl(L) - sky_nm(NST,NIM)
         enddo


      do NST = 1, NSTs
         ptot = 0.
         do L = 1, Ls
            if (rl(L).lt.2.5.and.nstl(L).eq.NST) then
               ptot = ptot + ppl(L)
               endif        
            enddo 
         print*,'NST: ',NST,ptot/NIMs,ptot
         z_n(NST) = ptot/NIMs
         enddo 

      do NST = 1, NSTs
      do NIM = 1, NIMs
         u_mn(NIM,NST) = 1
         enddo
         enddo

      do ipsf = 001, 101
      do jpsf = 001, 101
         psfo(ipsf,jpsf) = 0.00 
         psfu(ipsf,jpsf) = 0.00 
         enddo
         enddo

      do i = 0001, 0501
      do j = 0001, 1001
         show_psf(i,j) = 0.
         enddo
         enddo

      do NIT = 1, 6
         do ipsf = 001, 101
         do jpsf = 001, 101
            psfo(ipsf,jpsf) = psfu(ipsf,jpsf)
            enddo
            enddo

         print*,'                             '
         print*,'*****************************'
         print*,'***                          '
         print*,'***  NIT: ',NIT
         print*,'***                          '
         print*,'*****************************'
         print*,'                             '

         print*,' '
         write(STRING,'(''LOGA_FNDPSF.NIT'',i2.2,''simstV.OUT'')') NIT
         open(38,file=STRING,status='unknown')
         do ipsf = 001, 101
         do jpsf = 001, 101
            xpsf = (ipsf-051)*0.25
            ypsf = (jpsf-051)*0.25
            Os = 0
            do L = 1, Ls
               if (abs(dul(L)-xpsf).lt.0.5.and.
     .             abs(dvl(L)-ypsf).lt.0.5.and.
     .             u_mn(niml(L),nstl(L)).eq.1.and.
     .             psfuse_nstar(nstl(L))) then
                   Os = Os + 1
                   fo(Os) = ppl(L)/z_n(nstl(L)) 
     .                    - rpsf_phot(dul(L),dvl(L),psfu)
                   endif
               enddo
            call barsig(fo,Os,fbar,fsig,Ou)
            psfa(ipsf,jpsf) = fbar
            psfe(ipsf,jpsf) = fsig
            write(38,138) ipsf,jpsf,fbar,fsig,Os,Ou,psfu(ipsf,jpsf)
            if (ipsf.eq.51.and.jpsf.eq.51)
     .      write( *,138) ipsf,jpsf,fbar,fsig,Os,Ou,psfu(ipsf,jpsf)
  138       format(1x,i3.3,1x,i3.3,1x,f10.6,
     .                1x,f10.6,1x,i4,1x,i4,1x,f8.5)
            enddo
            enddo
         close(38)
         print*,' '

         do ipsf = 001, 101
         do jpsf = 001, 101
            psfn(ipsf,jpsf) = psfn(ipsf,jpsf) + psfa(ipsf,jpsf)
            enddo
            enddo
         call psf_smoo(psfn,psfs)
         call psf_norm(psfs,psfu)

         write(STRING,'(''LOGB_SHOPSF.NIT'',i2.2,''simstV.OUT'')') NIT
         open(37,file=STRING,status='unknown')
         do i = -500, 500, 1
         do j = -500, 500, 1
            xpsf = i*0.01
            ypsf = j*0.01
            write(37,137) xpsf,ypsf,rpsf_phot(xpsf,ypsf,psfu)
  137       format(1x,f8.3,1x,f8.3,1x,f10.6)
            enddo
            enddo
         close(37)

c
c determine zbar, etc
c
         print*,'SHOW THE FLUX AS MEASD FOR EACH STAR/IMAGE'
         do NST = 1, NSTs
            do NIM = 1, NIMs
               ptot = 0.
               ftot = 0.
               do L = 1, Ls
                  if (rl(L).lt.2.5.and.nstl(L).eq.NST.and.
     .                                 niml(L).eq.NIM) then
                     ptot = ptot + ppl(L)
                     ftot = ftot + rpsf_phot(dul(L),dvl(L),psfu)
                     endif        
                  enddo 
               z_mn(NIM,NST) = ptot/ftot
               enddo 
            call zm_anal(z_mn(1,NST),sz_mn(1,NST),
     .                                u_mn(1,NST),zbar,zsig,NIMs,NIMu)
            write(*,119) NST, psfuse_nstar(NST),
     .                        zbar,zsig,zsig/zbar,NIMs,NIMu,
     .                       (min(sz_mn(NIM,NST),99.99),
     .                             u_mn(NIM,NST),NIM=1,NIMs)
  119       format(1x,i2.2,l1,4x,f8.1,1x,f8.2,1x,f8.4,3x,i2,1x,i2,
     .             9x,_MDIM_(f6.3,'_',i1,1x))
            z_n(NST) = zbar
            enddo 


         write(STRING,'(''LOGC_SHOSAM.NIT'',i2.2,''simstV.OUT'')') NIT
         open(39,file=STRING,status='unknown')
         do L = 1, Ls
            write(39,139) L,dul(L),dvl(L),ppl(L),z_n(nstl(L)),
     .                      rpsf_phot(dul(L),dvl(L),psfu),
     .                      nstl(L)
  139       format(1x,i6.6,1x,f8.3,1x,f8.3,1x,f10.1,1x,f12.1,
     .                     f10.6,1x,i2)
            enddo
         close(39)

         do ipsf = 001+1, 101-1
         do jpsf = 001+1, 101-1
            show_psf(000+ipsf,100*(NIT-1)+jpsf) = psfo(ipsf,jpsf)
            show_psf(100+ipsf,100*(NIT-1)+jpsf) = psfa(ipsf,jpsf)
            show_psf(200+ipsf,100*(NIT-1)+jpsf) = psfs(ipsf,jpsf)-
     .                                            psfo(ipsf,jpsf)
            show_psf(300+ipsf,100*(NIT-1)+jpsf) = psfu(ipsf,jpsf)-
     .                                            psfs(ipsf,jpsf)
            show_psf(400+ipsf,100*(NIT-1)+jpsf) = psfu(ipsf,jpsf)-
     .                                            psfo(ipsf,jpsf)
            enddo
            enddo
            
         call writfits_r4('show_psf_simstV.fits',show_psf,501,1001)

         enddo!NIT

      print*,'OUTPUT psfout_simstV.fits...'
      call writfits_r4('psfout_simstV.fits',psfu,101,101)

      do i = 0001, 4001
      do j = 0001, 2001
         show_str(i,j) = 0.
         enddo
         enddo

      do NST = 1, NSTs
         if (.not.psfuse_nstar(NST)) then
            do ipsf = 001+2, 101-2
            do jpsf = 001+2, 101-2
               show_str((NST-1)*100+ipsf,jpsf+000) = -0.10
               enddo
               enddo
            endif
         do ipsf = 001+6, 101-6
         do jpsf = 001+6, 101-6
            xpsf = (ipsf-051)*0.25  
            ypsf = (jpsf-051)*0.25  
            Os = 0
            do L = 1, Ls
               if (abs(dul(L)-xpsf).le.0.5.and.
     .             abs(dvl(L)-ypsf).le.0.5.and.
     .             u_mn(niml(L),nstl(L)).eq.1.and.
     .             nstl(L).eq.NST) then
                   Os = Os + 1
                   po(Os) = ppl(L)/z_n(nstl(L)) 
                   fo(Os) = rpsf_phot(dul(L),dvl(L),psfu)
                   do(Os) = po(Os)-fo(Os)
                   endif
               enddo
            call barsig(po,Os,pbar,psig,Ou)
            call barsig(fo,Os,fbar,fsig,Ou)
            call barsig(do,Os,dbar,dsig,Ou)
            show_str((NST-1)*100+ipsf,jpsf+200) = fbar
            show_str((NST-1)*100+ipsf,jpsf+100) = pbar
            show_str((NST-1)*100+ipsf,jpsf+000) = dbar
            enddo
            enddo
         Us = 0
         do L = 1, Ls
            if (rl(L).le.2.5.and.
     .          u_mn(niml(L),nstl(L)).eq.1.and.
     .          nstl(L).eq.NST) then
               Us = Us + 1
               pu(Us) = ppl(L)
               xu(Us) = dul(L)
               yu(Us) = dvl(L)
               endif
            enddo
         emin = 9e9
         xmin = 0.0
         ymin = 0.0
         do idx = -100, 100
         do idy = -100, 100
            do U = 1, Us
               fu(U) = rpsf_phot(xu(U)-idx*0.001,
     .                           yu(U)-idy*0.001,psfu)
               enddo
            e = find_error(pu,fu,Us) 
            if (e.lt.emin) then
               emin = e
               xmin = idx*0.001
               ymin = idy*0.001
               endif
            enddo
            enddo
         write(*,'(''SHOW_STR: '',i2.2,1x,2f8.4,1x,f8.4,1x,i4,1x,l1)') 
     .                    NST,xmin,ymin,emin,Us,psfuse_nstar(NST)
         call writfits_r4('show_str_simstV.fits',show_str,4001,2001)
         enddo

      print*,'OUTPUT psfout_simstV.fits...'
      call writfits_r4('show_str_simstV.fits',show_str,4001,2001)


      print*,'SHOW INDIVIDUAL STARS...'
      do NST = 1, NSTs
         print*,'---> NST: ',NST,' / ',NSTs
      do NIM = 1, NIMs
         do ipsf = 001+2, 101-2
         do jpsf = 001+2, 101-2
            show_str((NST-1)*100+ipsf,(NIM-1+3)*100+jpsf) =  0.00
            if (u_mn(NIM,NST).eq.0)
     .      show_str((NST-1)*100+ipsf,(NIM-1+3)*100+jpsf) = -0.10
            enddo
            enddo
         do ipsf = 001+6, 101-6
         do jpsf = 001+6, 101-6
            xpsf = (ipsf-051)*0.25  
            ypsf = (jpsf-051)*0.25  
            Os = 0
            do L = 1, Ls
               if (abs(dul(L)-xpsf).le.0.5.and.
     .             abs(dvl(L)-ypsf).le.0.5.and.
     .             nstl(L).eq.NST.and.
     .             niml(L).eq.NIM) then
                   Os = Os + 1
                   po(Os) = ppl(L)/z_n(nstl(L)) 
                   fo(Os) = rpsf_phot(dul(L),dvl(L),psfu)
                   do(Os) = po(Os)-fo(Os)
                   endif
               enddo
            call barsig(po,Os,pbar,psig,Ou)
            call barsig(fo,Os,fbar,fsig,Ou)
            call barsig(do,Os,dbar,dsig,Ou)
            show_str((NST-1)*100+ipsf,(NIM-1+3)*100+jpsf) = dbar
            enddo
            enddo
         enddo
         call writfits_r4('show_str_simstV.fits',show_str,4001,2001)
         enddo

      print*,'OUTPUT show_str_simstV.fits...'
      call writfits_r4('show_str_simstV.fits',show_str,4001,2001)

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
 
 


c
c
c
c
c

      subroutine readfits_r4(FILE,pix,NDIMX,NDIMY)
      implicit none

      character*80 FILE
      integer NDIMX,NDIMY
      real    pix(NDIMX,NDIMY)

      character*80 FILEU
      character*70 INFO(10)
      common / fitsinfo / INFO 

      integer naxes
      integer laxis(3)
      common/laxis3_/laxis

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
      common/HDR/HDR

      FILEU = FILE
      do i = 75,2,-1
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

      i = 0
      nread = 0
 100  continue
      i = i + 1
      read(10,rec=i,iostat=ios) buffc
      if (DIAG) print*,'READREC: ',i
      do k = 0, 35, 1
         if (DIAG) write(*,'(i4,1x,i4,1x,a80)') 
     .                   i,k,buffc(k*80+1:k*80+80)
         field  = buffc(k*80+01:k*80+08)
         stream = buffc(k*80+10:k*80+79)
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
         print*,'BITPIX.ne.-32... unreal!'
         stop
         endif

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

         if (DIAG) write(*,1115) i,np1,np2,npt,nxx,nyy
 1115    format(1x,i8,1x,i10,1x,i10,1x,i10,1x,2i6)
         enddo 

  899 close(10)
      if (DIAG) write(*,1115) i,np1,np2,npt,nxx,nyy



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

      do i = 1, 720
         npu = n1+i-1
         nbu = (i-1)*4
         if (.not.(_LINUX_)) then
            b(1) = buff(nbu+1)
            b(2) = buff(nbu+2)
            b(3) = buff(nbu+3)
            b(4) = buff(nbu+4)
            endif
         if ((_LINUX_)) then
            b(4) = buff(nbu+1)
            b(3) = buff(nbu+2)
            b(2) = buff(nbu+3)
            b(1) = buff(nbu+4)
            endif
         if (npu.ge.1.and.npu.le.nt) pix(npu) = r
         enddo

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

ccc      character*80 FILE
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
      do i = 75,2,-1
         if (FILE(i:i+4).eq.'.fits') FILEU = FILE(1:i+4)
         enddo

      open(10,file=FILEU,status='unknown',
     .     err=900,recl=2880,form='UNFORMATTED',
     .     access='DIRECT')

      write(buffc( 0*80+1: 1*80),'(''SIMPLE  =                    T'')')
      write(buffc( 1*80+1: 2*80),'(''BITPIX  =                  -32'')')
      write(buffc( 2*80+1: 3*80),'(''NAXIS   ='',8x,i12)') 2
      write(buffc( 3*80+1: 4*80),'(''NAXIS1  ='',8x,i12)') PXDIMX
      write(buffc( 4*80+1: 5*80),'(''NAXIS2  ='',8x,i12)') PXDIMY
      write(buffc( 5*80+1: 6*80),'(''DATATYPE='',9a)') 
     .                          " 'REAL*4' "
      write(buffc(07*80+1:08*80),'(''COMMENT  '',a05)') '     '
      write(buffc(08*80+1:09*80),'(''COMMENT  '',a05)') '     '
      write(buffc(09*80+1:10*80),'(''CRPIX1  ='',a70)') HDR(01)
      write(buffc(10*80+1:11*80),'(''CRPIX2  ='',a70)') HDR(02)
      write(buffc(11*80+1:12*80),'(''CRVAL1  ='',a70)') HDR(03)
      write(buffc(12*80+1:13*80),'(''CRVAL2  ='',a70)') HDR(04)
      write(buffc(13*80+1:14*80),'(''CTYPE1  ='',a70)') HDR(05)
      write(buffc(14*80+1:15*80),'(''CTYPE2  ='',a70)') HDR(06)
      write(buffc(15*80+1:16*80),'(''CD1_1   ='',a70)') HDR(07)
      write(buffc(16*80+1:17*80),'(''CD1_2   ='',a70)') HDR(08)
      write(buffc(17*80+1:18*80),'(''CD2_1   ='',a70)') HDR(09)
      write(buffc(18*80+1:19*80),'(''CD2_2   ='',a70)') HDR(10)
      write(buffc(19*80+1:20*80),'(''ORIENTAT='',a70)') HDR(11)
      write(buffc(20*80+1:21*80),'(''PA_APER ='',a70)') HDR(12)
      write(buffc(21*80+1:22*80),'(''PA_V3   ='',a70)') HDR(13)
      write(buffc(22*80+1:23*80),'(''DATE-OBS='',a70)') HDR(14)
      write(buffc(23*80+1:24*80),'(''TIME-OBS='',a70)') HDR(15)
      write(buffc(24*80+1:25*80),'(''EXPTIME ='',a70)') HDR(16)
      write(buffc(25*80+1:26*80),'(''ROOTNAME='',a70)') HDR(17)
      write(buffc(26*80+1:27*80),'(''TARGNAME='',a70)') HDR(18)
      write(buffc(27*80+1:28*80),'(''RA_TARG ='',a70)') HDR(19)
      write(buffc(28*80+1:29*80),'(''DEC_TARG='',a70)') HDR(20)
      write(buffc(29*80+1:30*80),'(''PROPOSID='',a70)') HDR(21)
      write(buffc(30*80+1:31*80),'(''FILTER1 ='',a70)') HDR(22)
      write(buffc(31*80+1:32*80),'(''FILTER2 ='',a70)') HDR(23)
      write(buffc(33*80+1:34*80),'(''VAFACTOR='',a70)') HDR(24)
      write(buffc(32*80+1:33*80),'(''CCDGAIN ='',a70)') HDR(25)
      write(buffc(33*80+1:34*80),'(''COMMENT  '',a05)') '     '
      write(buffc(34*80+1:35*80),'(''COMMENT  '',a05)') '     '
      write(buffc(35*80+1:36*80),'(''END      '')')

      do k = 00, 34
         if (buffc(k*80+11:k*80+30).eq.'                    ')
     .       write(buffc(k*80+01:k*80+80),'(80('' ''))')
         enddo

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

      do i = 1, 720
         npu = n1+i-1
         nbu = (i-1)*4
         if (npu.ge.1.and.npu.le.nt) r = pix(npu)
         if (.not.(_LINUX_)) then
            buff(nbu+1) = b(1)
            buff(nbu+2) = b(2)
            buff(nbu+3) = b(3)
            buff(nbu+4) = b(4)
            endif
         if ((_LINUX_)) then
            buff(nbu+1) = b(4)
            buff(nbu+2) = b(3)
            buff(nbu+3) = b(2)
            buff(nbu+4) = b(1)
            endif
         enddo

      return
      end


