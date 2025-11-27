       program mcmc_expand_average

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       parameter (nmax=1200000)
       double precision x1a(nmax),y1a(nmax),x2a(nmax),y2a(nmax),
     &           x3a(nmax),y3a(nmax),f1a(nmax),f2a(nmax),A0a(nmax)
       double precision x12flxwt(nmax),y12flxwt(nmax),
     &                  flux1a(nmax),flux2a(nmax)
       double precision x2mx1a(nmax),y2my1a(nmax)
       character*120 infile
       character*200 string

       call getarg(1,infile)
       lin = lenc(infile)

       open(unit=1,file=infile,status='old')
       open(unit=2,file='expanded_mcmc.txt',status='unknown')

       nskip = 50000

       x1av = 0.d0
       y1av = 0.d0
       f1av = 0.d0
       x2av = 0.d0
       y2av = 0.d0
       f2av = 0.d0
       x3av = 0.d0
       y3av = 0.d0
       flux1av = 0.d0
       flux2av = 0.d0
       x2mx1av = 0.d0
       y2my1av = 0.d0
       A0av = 0.d0
       x12flxwtav = 0.d0
       y12flxwtav = 0.d0
       do i = 1,9999999
         read(1,'(a)',end=11) string
         if(i.eq.1) len1 = lenc(string)
         read(string,*) irep
         if(irep.eq.0) then
           if(len1.gt.110) then
             read(string,*) jj,x1,y1,x2,y2,x3,y3,f1,f2,chi2tot,A0
           else
             read(string,*) jj,x1,y1,x2,y2,x3,y3,f1,f2,A0,chi2tot
           endif
         endif
         write(2,96) x1,y1,x2,y2,x3,y3,f1,f2,A0,
     &               f1*A0,f2*A0,(1.d0-f1-f2)*A0,chi2tot
         if(i.gt.nskip) then
           x1a(i) = x1
           y1a(i) = y1
           f1a(i) = f1
           x2a(i) = x2
           y2a(i) = y2
           f2a(i) = f2
           x3a(i) = x3
           y3a(i) = y3
           A0a(i) = A0
           x2mx1a(i) = x2 - x1
           y2my1a(i) = y2 - y1
           flux1a(i) = f1*A0
           flux2a(i) = f2*A0
           flux12 = f1a(i) + f2a(i)
           x12flxwt(i) = (f1a(i)*x1a(i)+f2a(i)*x2a(i))/flux12
           y12flxwt(i) = (f1a(i)*y1a(i)+f2a(i)*y2a(i))/flux12
           x1av = x1av + x1a(i)
           y1av = y1av + y1a(i)
           f1av = f1av + f1a(i)
           x2av = x2av + x2a(i)
           y2av = y2av + y2a(i)
           f2av = f2av + f2a(i)
           x3av = x3av + x3a(i)
           y3av = y3av + y3a(i)
           x2mx1av = x2mx1av + x2mx1a(i)
           y2my1av = y2my1av + y2my1a(i)
           flux1av = flux1av + flux1a(i)
           flux2av = flux2av + flux2a(i)
           A0av = A0av + A0a(i)
           x12flxwtav = x12flxwtav + x12flxwt(i)
           y12flxwtav = y12flxwtav + y12flxwt(i)
         endif
       enddo

 11    continue
       close(1)
       close(2)

       n = i -1
       write(6,*) n,' lines expanded'
       nsaved = n - nskip
       x1av = x1av/nsaved
       y1av = y1av/nsaved
       f1av = f1av/nsaved
       x2av = x2av/nsaved
       y2av = y2av/nsaved
       f2av = f2av/nsaved
       x3av = x3av/nsaved
       x2mx1av = x2mx1av/nsaved
       y2my1av = y2my1av/nsaved
       y3av = y3av/nsaved
       flux1av = flux1av/nsaved
       flux2av = flux2av/nsaved
       A0av = A0av/nsaved
       x12flxwtav = x12flxwtav/nsaved
       y12flxwtav = y12flxwtav/nsaved

       x1sig = 0.d0
       y1sig = 0.d0
       f1sig = 0.d0
       x2sig = 0.d0
       y2sig = 0.d0
       f2sig = 0.d0
       x3sig = 0.d0
       y3sig = 0.d0
       x2mx1sig = 0.d0
       y2my1sig = 0.d0
       x1sig = 0.d0
       x1sig = 0.d0
       flux1sig = 0.d0
       flux2sig = 0.d0
       A0sig = 0.d0
       do i = 1,n
         if(i.gt.nskip) then
           x1sig = x1sig + (x1a(i)-x1av)**2
           y1sig = y1sig + (y1a(i)-y1av)**2
           f1sig = f1sig + (f1a(i)-f1av)**2
           x2sig = x2sig + (x2a(i)-x2av)**2
           y2sig = y2sig + (y2a(i)-y2av)**2
           f2sig = f2sig + (f2a(i)-f2av)**2
           x3sig = x3sig + (x3a(i)-x3av)**2
           x2mx1sig = x2mx1sig + (x2mx1a(i)-x2mx1av)**2
           y2my1sig = y2my1sig + (y2my1a(i)-y2my1av)**2
           y3sig = y3sig + (y3a(i)-y3av)**2
           flux1sig = flux1sig + (flux1a(i)-flux1av)**2
           flux2sig = flux2sig + (flux2a(i)-flux2av)**2
           A0sig = A0sig + (A0a(i)-A0av)**2
           x12flxwtsig = x12flxwtsig + (x12flxwt(i)-x12flxwtv)**2
           y12flxwtsig = y12flxwtsig + (y12flxwt(i)-y12flxwtv)**2
         endif
       enddo
       x1sig = sqrt(x1sig/nsaved)
       y1sig = sqrt(y1sig/nsaved)
       f1sig = sqrt(f1sig/nsaved)
       x2sig = sqrt(x2sig/nsaved)
       y2sig = sqrt(y2sig/nsaved)
       f2sig = sqrt(f2sig/nsaved)
       x3sig = sqrt(x3sig/nsaved)
       y3sig = sqrt(y3sig/nsaved)
       x2mx1sig = sqrt(x2mx1sig/nsaved)
       y2my1sig = sqrt(y2my1sig/nsaved)
       flux1sig = sqrt(flux1sig/nsaved)
       flux2sig = sqrt(flux2sig/nsaved)
       A0sig = sqrt(A0sig/nsaved)
       x12flxwtsig = sqrt(x12flxwtsig/nsaved)
       y12flxwtsig = sqrt(y12flxwtsig/nsaved)

       write(6,97) x1av,y1av,x2av,y2av,x3av,y3av,f1av,f2av,A0av,
     &             x1sig,y1sig,x2sig,y2sig,x3sig,y3sig,f1sig,f2sig,A0sig
       write(6,98) flux1av,flux2av,x2mx1av,y2my1av,
     &             flux1sig,flux2sig,x2mx1sig,y2my1sig
       write(6,99) x12flxwtav,x12flxwtsig,y12flxwtav,y12flxwtsig

 96    format(8f9.4,5f10.4)
 97    format('       x1      y1       x2       y2       x3       y3',
     &      '       f1       f2         A0',/,
     &      3x,8f9.4,f11.4,/,'RMS',8f9.4,f11.4)
 98    format('      flux1      flux2      x2-x1    y2-y1',/,
     &      3x,2f11.4,2f9.4,/,'RMS',2f11.4,2f9.4)
 99    format('  Weighted average mass 1 and 2 positions:',/,
     &        '      x                        y',/,f12.4,' +-',f9.4,
     &         f12.4,' +-',f9.4)

       stop
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      function lenc(s)
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
