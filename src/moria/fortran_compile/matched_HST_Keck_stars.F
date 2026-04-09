       program matched_HST_Keck_stars

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       parameter(maxKeck=3000,maxHST=10000)
       double precision xk(maxKeck),yk(maxKeck),mk(maxKeck)
       double precision mksig(maxKeck)
       double precision xh(maxHST),yh(maxHST),
     &            xhsig(maxHST),yhsig(maxHST),mh(maxHST),mhsig(maxHST)
       double precision Keck_x,Keck_y,Keck_aprx_err,Keck_est_err
       integer idk(maxKeck)

       character*120 Keckfile,HSTfile,outfile
       character*104 cline
       character*5 HST_fil,Keck_fil

c      hardwired parameter
c      -------------------
       hst_mag_lim = -7.0d0
       hst_aprx_err = 1.5d0
       Keck_aprx_err = 5.d0
       hst_est_err = 0.5d0
       Keck_est_err = 1.0d0

c      input files and output file
c      ---------------------------
       write(6,*) 'enter HST input file name:'
       read(5,'(a)') HSTfile
       write(6,*) 'enter Keck input file name (ALLSTAR *.als format:'
       read(5,'(a)') Keckfile
       write(6,*) 'enter HST filter (V, I, F555W, F606W, F814W, ...):'
       read(5,'(a)') HST_fil
       lHST_fil = lenc(HST_fil)
       if(lHST_fil.eq.1) then
         HST_fil = HST_fil(1:1)//'    '
       elseif(lHST_fil.eq.2) then
         HST_fil = HST_fil(1:2)//'   '
       elseif(lHST_fil.eq.3) then
         HST_fil = HST_fil(1:3)//'  '
       elseif(lHST_fil.eq.4) then
         HST_fil = HST_fil(1:4)//' '
       endif
       write(6,*) 'enter Keck filter (K, Ks, H, ...):'
       read(5,'(a)') Keck_fil
       lKeck_fil = lenc(Keck_fil)
       if(lKeck_fil.eq.1) then
         Keck_fil = Keck_fil(1:1)//'    '
       elseif(lKeck_fil.eq.2) then
         Keck_fil = Keck_fil(1:2)//'   '
       elseif(lKeck_fil.eq.3) then
         Keck_fil = Keck_fil(1:3)//'  '
       elseif(lKeck_fil.eq.4) then
         Keck_fil = Keck_fil(1:4)//' '
       endif
       write(6,*) 'enter output file name:'
       read(5,'(a)') outfile

c      read HST input file
c      -------------------
       nHST = 0
       open(unit=1,file=HSTfile,status='old')
       do i = 1,999999
         read(1,'(a)',end=11,err=11) cline
         if(cline(1:1).ne.'#') then
           read(cline,*) xin,yin,fmagin,xsig,ysig,fmagsig
           if(fmagin.le.hst_mag_lim) then
c            keep only stars that are brighter than hst_mag_lim
c            --------------------------------------------------
             nHST = nHST + 1
             if(nHST.gt.maxHST) then
               write(6,*) nHST,' is too large'
               stop
             endif
             xh(nHST) = xin
             yh(nHST) = yin
             xhsig(nHST) = xsig
             yhsig(nHST) = ysig
             mh(nHST) = fmagin
             mhsig(nHST) = fmagsig
           endif
         endif
       enddo
 11    continue
       close(1)
       write(6,*) 'HST stars =', nHST

c      read Keck input file
c      --------------------
       nKeck = 0
       open(unit=2,file=Keckfile,status='old')
       read(2,*)
       read(2,*)
       read(2,*)
       do i = 1,9999
         read(2,*,end=21,err=21) idkin,xin,yin,fmag,fmagsig
         nKeck = Nkeck + 1
         if(nKeck.gt.maxKeck) then
           write(6,*) nKeck,' is too large'
           stop
         endif
         idk(Nkeck) = idkin
         xk(Nkeck) = xin
         yk(Nkeck) = yin
         mk(Nkeck) = fmag
         mksig(Nkeck) = fmagsig
       enddo
 21    continue
       close(2)
       write(6,*) 'Keck stars =', nKeck

       open(unit=3,file=outfile,status='unknown')
       write(3,'(a)') 'HST '//HST_fil//
     &   '                                                      '//
     &   'Keck '//Keck_fil
       write(6,'(a)') 'HST '//HST_fil//
     &   '                                                      '//
     &   'Keck '//Keck_fil
       write(3,91)

c      read approximate matched coordinates
c      ------------------------------------
       open(unit=4,file='IN.approx_matched_HST_Keck_stars',status='old')
       write(6,*) 'reading matched coords HST_x, HST_y, Keck_x, Keck_y:'
       nHST_keck = 0
       do i = 1,9999
         read(4,*,end=31,err=31) HST_x, HST_y, Keck_x, Keck_y, itarg
         nHST_keck = nHST_keck + 1
c        find the closest HST star
c        -------------------------
         chi2minHST = 1.d+10
         ihbest = 0
         do ih = 1,nHST
           if(xh(ih).lt.HST_x-hst_aprx_err) then
           elseif(xh(ih).gt.HST_x+hst_aprx_err) then
           elseif(yh(ih).lt.HST_y-hst_aprx_err) then
           elseif(yh(ih).gt.HST_y+hst_aprx_err) then
           else
             chi2x = (xh(ih)-HST_x)**2/(xhsig(ih)**2+hst_est_err**2)
             chi2y = (yh(ih)-HST_y)**2/(yhsig(ih)**2+hst_est_err**2)
             chi2 = chi2x + chi2y
             if(chi2.lt.chi2minHST) then
               chi2minHST = chi2
               ihbest = ih
             endif
           endif
         enddo
c        find the closest Keck star
c        --------------------------
         chi2minKeck = 1.d+10
         ikbest = 0
         do ik = 1,nKeck
           if(xk(ik).lt.Keck_x-Keck_aprx_err) then
           elseif(xk(ik).gt.Keck_x+Keck_aprx_err) then
           elseif(yk(ik).lt.Keck_y-Keck_aprx_err) then
           elseif(yk(ik).gt.Keck_y+Keck_aprx_err) then
           else
             chi2x = (xk(ik)-Keck_x)**2/(Keck_est_err**2)
             chi2y = (yk(ik)-Keck_y)**2/(Keck_est_err**2)
             chi2 = chi2x + chi2y
             if(itarg.gt.0) then
               write(6,*) 'Target star matches:'
               write(6,92) xh(ihbest),yh(ihbest),mh(ihbest),
     &               xhsig(ihbest),yhsig(ihbest),mhsig(ihbest),
     &               idk(ik),xk(ik),yk(ik),mk(ik),mksig(ik),itarg
             endif
             if(chi2.lt.chi2minKeck) then
               chi2minKeck = chi2
               ikbest = ik
             endif
           endif
         enddo
         if(ihbest.eq.0.and.ikbest.eq.0) then
           write(6,*) 'No HST or Keck star found for input star:',i
           write(6,92) 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &               0,0.d0,0.d0,0.d0,0.d0,itarg
         elseif(ihbest.eq.0) then
           write(6,*) 'No HST star found for input star:',i
           write(6,92) 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,
     &               idk(ikbest),xk(ikbest),
     &               yk(ikbest),mk(ikbest),mksig(ikbest),itarg
         elseif(ikbest.eq.0) then
           write(6,*) 'No Keck star found for input star:',i
           write(6,92) xh(ihbest),yh(ihbest),mh(ihbest),xhsig(ihbest),
     &               yhsig(ihbest),mhsig(ihbest),0,0.d0,0.d0,0.d0,0.d0,
     &               itarg
         else
           write(3,92) xh(ihbest),yh(ihbest),mh(ihbest),xhsig(ihbest),
     &               yhsig(ihbest),mhsig(ihbest),idk(ikbest),xk(ikbest),
     &               yk(ikbest),mk(ikbest),mksig(ikbest),itarg
         endif
       enddo
 31    continue

 91    format('    xHST      yHST    HSTmag   x_err    y_err  ',
     &        'HSTmag_err  #    xKeck    yKeck  Keckmag Kmag_err')
 92    format(2f10.4,4f9.4,i5,2f9.3,2f8.4,i4)

cHST Ixxxx                                                      Keck Kxxxx
c    xHST      yHST    HSTmag   x_err    y_err  HSTmag_err  #    xKeck    yKeck  Keckmag Kmag_err
c  470.1146  591.0200  -8.8939   0.0204   0.0144   0.0251  319  479.923  594.873  13.983  0.0029
c12345678901234567890123456789123456789123456789123456789123451234567891234567891234567812345678

       stop
       END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      function lenc(s)
c
c  Get the length of a character variable.  Only necessary because of the
c  Unix lobotomy of Fortran.
c
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

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
