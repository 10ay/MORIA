       program mcmc_expand_average
c
c  Expand MCMC chain from uvp2tri_scon *_07.mcmc (Python or legacy Fortran).
c  Lines may begin with optional headers:
c     # NSTAR N
c     # NIMS M
c  If absent, N=3 (legacy tripod).  Each full row: i4, (3N-1)*f10.5, 2f13.4, M*f10.4
c
       implicit double precision (a-h,o-z)
       integer, parameter :: nmax=1200000
       integer, parameter :: nstarm=20
       double precision x(nstarm), y(nstarm), fp(nstarm-1)
       double precision ff(nstarm)
       double precision xhist(nmax,nstarm), yhist(nmax,nstarm)
       double precision fphist(nmax,nstarm-1)
       double precision x2mx1a(nmax), y2my1a(nmax)
       double precision flux1a(nmax), flux2a(nmax)
       double precision x12flxwt(nmax), y12flxwt(nmax)
       double precision A0a(nmax)
       double precision xav(nstarm), yav(nstarm), fpav(nstarm-1)
       double precision xsig(nstarm), ysig(nstarm), fpsig(nstarm-1)
       double precision chi2tot, A0, fn, flux12
       character*120 infile
       character*600 string
       integer irep, jj, k, ls, len1d, n, nsaved, nlines, nims, i
       integer ndat, nskip_req, nskip_eff
       logical first

       call getarg(1,infile)
       if (lenc(infile).lt.1) then
         write(6,*) 'usage: mcmc_expand_average.xOg <.07.mcmc>'
         stop
       endif

       open(unit=1,file=infile,status='old')
       open(unit=2,file='expanded_mcmc.txt',status='unknown')

       nskip_req = 50000
       nstar = 0
       nims = 0
       len1d = 0

c  Skip comment / header lines; first non-# line is start of data
 5     read(1,'(a)',end=88) string
       ls = lenc(string)
       if (ls.lt.1) goto 5
       if (string(1:1).eq.'#') then
         if (ls.ge.8 .and. string(1:8).eq.'# NSTAR')
     &       read(string(9:),*,err=5) nstar
         if (ls.ge.7 .and. string(1:7).eq.'# NIMS')
     &       read(string(8:),*,err=5) nims
         goto 5
       endif
       if (nstar.eq.0) nstar = 3
       if (nstar.lt.1 .or. nstar.gt.nstarm) then
         write(6,*) 'Invalid NSTAR', nstar
         stop 1
       endif

       x = 0.d0
       y = 0.d0
       fp = 0.d0
       flux1av = 0.d0
       flux2av = 0.d0
       x2mx1av = 0.d0
       y2my1av = 0.d0
       x12flxwtav = 0.d0
       y12flxwtav = 0.d0
       A0av = 0.d0
       do k = 1, nstarm
         xav(k) = 0.d0
         yav(k) = 0.d0
       enddo
       do k = 1, nstarm-1
         fpav(k) = 0.d0
       enddo

       chi2tot = 0.d0
       A0 = 0.d0
       nlines = 0
       ndat = 0
       first = .true.

 301   continue
       if (.not.first) read(1,'(a)',end=11) string
       first = .false.
       nlines = nlines + 1
       ls = lenc(string)
       if (ls.lt.1) goto 301
       if (string(1:1).eq.'#') goto 301

       read(string,*) irep

         if (irep.ne.0) goto 301

           if (len1d.eq.0) len1d = ls
           if (nstar.eq.1) then
             if (len1d.gt.110) then
               read(string,*) jj,x(1),y(1),chi2tot,A0
             else
               read(string,*) jj,x(1),y(1),A0,chi2tot
             endif
           else
             if (len1d.gt.110) then
               read(string,*) jj,(x(k),y(k),k=1,nstar),
     &             (fp(k),k=1,nstar-1),chi2tot,A0
             else
               read(string,*) jj,(x(k),y(k),k=1,nstar),
     &             (fp(k),k=1,nstar-1),A0,chi2tot
             endif
           endif

         fn = 1.d0
         do k = 1, nstar-1
           fn = fn - fp(k)
           ff(k) = fp(k)
         enddo
         ff(nstar) = fn

         write(2,'(60f14.6)') (x(k),y(k),k=1,nstar),
     &     (fp(k),k=1,nstar-1), fn, chi2tot, A0,
     &     (ff(k)*A0,k=1,nstar)

           ndat = ndat + 1
           if (ndat.gt.nmax) then
             write(6,*) 'Too many irep=0 rows, increase nmax:', nmax
             stop 1
           endif
           do k = 1, nstar
             xhist(ndat,k) = x(k)
             yhist(ndat,k) = y(k)
           enddo
           do k = 1, nstar-1
             fphist(ndat,k) = fp(k)
           enddo
           A0a(ndat) = A0
           if (nstar.ge.2) then
             x2mx1a(ndat) = x(2) - x(1)
             y2my1a(ndat) = y(2) - y(1)
             flux1a(ndat) = fp(1)*A0
             flux2a(ndat) = fp(2)*A0
             flux12 = flux1a(ndat) + flux2a(ndat)
             if (flux12.gt.0.d0) then
               x12flxwt(ndat)=(flux1a(ndat)*x(1)+flux2a(ndat)
     &           *x(2))/flux12
               y12flxwt(ndat)=(flux1a(ndat)*y(1)+flux2a(ndat)
     &           *y(2))/flux12
             else
               x12flxwt(ndat) = 0.d0
               y12flxwt(ndat) = 0.d0
             endif
           endif

       goto 301

 11    continue
       close(1)
       close(2)

       n = ndat
       write(6,*) n,' data rows (irep=0);', nlines,' file lines read'
       if (n.lt.1) then
         write(6,*) 'No irep=0 data rows'
         stop 1
       endif
       if (n.gt.nskip_req) then
         nskip_eff = nskip_req
       else
         nskip_eff = max(0, min(nskip_req, n/10))
       endif
       nsaved = n - nskip_eff
       if (nsaved.lt.10) then
         nskip_eff = max(0, n - 10)
         nsaved = n - nskip_eff
       endif
       if (nsaved.lt.1) then
         write(6,*) 'Not enough data rows after burn-in:', n
         stop 1
       endif
       write(6,*) 'burn-in skip', nskip_eff,
     & ' retained', nsaved
       do i = nskip_eff+1, n
         do k = 1, nstar
           xav(k) = xav(k) + xhist(i,k)
           yav(k) = yav(k) + yhist(i,k)
         enddo
         do k = 1, nstar-1
           fpav(k) = fpav(k) + fphist(i,k)
         enddo
         A0av = A0av + A0a(i)
         if (nstar.ge.2) then
           x2mx1av = x2mx1av + x2mx1a(i)
           y2my1av = y2my1av + y2my1a(i)
           flux1av = flux1av + flux1a(i)
           flux2av = flux2av + flux2a(i)
           x12flxwtav = x12flxwtav + x12flxwt(i)
           y12flxwtav = y12flxwtav + y12flxwt(i)
         endif
       enddo
       do k = 1, nstar
         xav(k) = xav(k) / dble(nsaved)
         yav(k) = yav(k) / dble(nsaved)
       enddo
       do k = 1, nstar-1
         fpav(k) = fpav(k) / dble(nsaved)
       enddo
       A0av = A0av / dble(nsaved)
       if (nstar.ge.2) then
         x2mx1av = x2mx1av / dble(nsaved)
         y2my1av = y2my1av / dble(nsaved)
         flux1av = flux1av / dble(nsaved)
         flux2av = flux2av / dble(nsaved)
         x12flxwtav = x12flxwtav / dble(nsaved)
         y12flxwtav = y12flxwtav / dble(nsaved)
       endif

       do k = 1, nstarm
         xsig(k) = 0.d0
         ysig(k) = 0.d0
       enddo
       do k = 1, nstarm-1
         fpsig(k) = 0.d0
       enddo
       x2mx1sig = 0.d0
       y2my1sig = 0.d0
       flux1sig = 0.d0
       flux2sig = 0.d0
       A0sig = 0.d0
       x12flxwtsig = 0.d0
       y12flxwtsig = 0.d0

       do i = 1, n
         if (i.gt.nskip_eff) then
           do k = 1, nstar
             xsig(k) = xsig(k) + (xhist(i,k)-xav(k))**2
             ysig(k) = ysig(k) + (yhist(i,k)-yav(k))**2
           enddo
           do k = 1, nstar-1
             fpsig(k)=fpsig(k)+(fphist(i,k)-fpav(k))**2
           enddo
           A0sig = A0sig + (A0a(i)-A0av)**2
           if (nstar.ge.2) then
             x2mx1sig=x2mx1sig+(x2mx1a(i)-x2mx1av)**2
             y2my1sig=y2my1sig+(y2my1a(i)-y2my1av)**2
             flux1sig=flux1sig+(flux1a(i)-flux1av)**2
             flux2sig=flux2sig+(flux2a(i)-flux2av)**2
             x12flxwtsig=x12flxwtsig+(x12flxwt(i)-x12flxwtav)**2
             y12flxwtsig=y12flxwtsig+(y12flxwt(i)-y12flxwtav)**2
           endif
         endif
       enddo

       do k = 1, nstar
         xsig(k) = sqrt(xsig(k)/dble(nsaved))
         ysig(k) = sqrt(ysig(k)/dble(nsaved))
       enddo
       do k = 1, nstar-1
         fpsig(k) = sqrt(fpsig(k)/dble(nsaved))
       enddo
       A0sig = sqrt(A0sig/dble(nsaved))
       if (nstar.ge.2) then
         x2mx1sig = sqrt(x2mx1sig/dble(nsaved))
         y2my1sig = sqrt(y2my1sig/dble(nsaved))
         flux1sig = sqrt(flux1sig/dble(nsaved))
         flux2sig = sqrt(flux2sig/dble(nsaved))
         x12flxwtsig = sqrt(x12flxwtsig/dble(nsaved))
         y12flxwtsig = sqrt(y12flxwtsig/dble(nsaved))
       endif

       write(6,*) 'NSTAR', nstar
       write(6,*) 'means (x_k, y_k; k=1..N):'
       do k = 1, nstar
         write(6,'(a,i2,4f12.5)') '  k', k,
     &     xav(k), xsig(k), yav(k), ysig(k)
       enddo
       if (nstar.ge.2) then
         write(6,*) 'f1..f_{N-1} means / RMS:'
         do k = 1, nstar-1
           write(6,'(a,i2,2f12.5)') '  f', k, fpav(k), fpsig(k)
         enddo
         write(6,'(a,f12.5,a,f12.5)') '  A0 mean', A0av, '  RMS', A0sig
         write(6,*) 'flux1 flux2 (stars 1-2), sep:'
         write(6,'(4f12.5)') flux1av, flux2av, x2mx1av, y2my1av
         write(6,'(a,4f12.5)') 'RMS', flux1sig,
     &     flux2sig, x2mx1sig, y2my1sig
         write(6,*) 'Weighted average mass 1 and 2 positions:'
         write(6,*) '      x                        y'
         write(6,*) 'wtd avg x=', x12flxwtav, ' +-', x12flxwtsig
         write(6,*) 'wtd avg y=', y12flxwtav, ' +-', y12flxwtsig
       else
         write(6,'(a,f12.5,a,f12.5)') 'A0 mean', A0av, '  RMS', A0sig
       endif

       stop

 88    write(6,*) 'EOF before any data line'
       stop 1
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function lenc(s)
      character*(*) s
      lenc = len(s)
      do 900 i = len(s), 1, -1
          if (s(i:i).ne.' ') then
              lenc = i
              return
          end if
900   continue
      lenc = 0
      return
      end
