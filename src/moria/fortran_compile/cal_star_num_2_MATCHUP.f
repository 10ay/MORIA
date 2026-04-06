       program cal_star_num_2_MATCHUP

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       parameter(maxstars=30000,nbmax=100,Nfs=100)
       real*8 xuref(maxstars),yuref(maxstars),miuref(maxstars),
     &        mvuref(maxstars),upsf(maxstars)
       integer iMuref(maxstars)
       character*160 inline,inline2,inline3
       character*80 matchup_file_in,matchup_file_out,calstarfile,
     &              matchup_file_Cal_only_out

c      enter Calibration star list
c      ---------------------------
       write(6,*) 'enter Calibration star list (NOTFAR file):'
       read(5,'(a)') calstarfile
       open(unit=11,file=calstarfile,status='old')
       read(11,*)
       read(11,*)
       do ic = 1,99999
         read(11,*,end=2,err=2) xuref(ic),yuref(ic),miuref(ic),
     &      mvuref(ic),upsf(ic),iMuref(ic)
       enddo
 2     continue
       ncal = ic - 1

c      open and read Jay Anderson format MATCHUP.XYMEEE files (F814W & F555W)
c      ----------------------------------------------------------------------
       write(6,*) 'enter matchup_file_in and matchup_file_out names:'
       read(5,'(a)') matchup_file_in
       read(5,'(a)') matchup_file_out
       read(5,'(a)') matchup_file_Cal_only_out
       open(unit=2,file=matchup_file_in,status='old')
       open(unit=3,file=matchup_file_out,status='new')
       open(unit=4,file=matchup_file_Cal_only_out,status='new')

c      MATCHUP file lists for each star
c      x, y, mag, x_sig, y_sig, mag_sig, NIMf, NIMg, NIMm, NIMMINu,NPK, pki(NPK),pkj(NPK), 
c      new MATCHUP file adds
c        calibration star # (only for calibration stars)
       j = 0
       do i=1,99999
         read(2,'(a)',err=31,end=31) inline
         lin = lenc(inline)
         inline3 = inline
         if(inline(1:1).ne.'#') then
           lout = lin
           j = j + 1
           read(inline,*) xm, ym
           do ic = 1,ncal
             dr2 = (xuref(ic)-xm)**2 + (yuref(ic)-ym)**2
             if(dr2.lt.0.01d0) then
ccc             if(j+1.eq.iMuref(ic)) then
               lout = lin + 7
               write(inline(lin+1:lout),'(i7)') ic
               write(inline3(37:43),'(i7)') ic
               write(4,'(a)') inline3(1:43)
             endif
           enddo
         elseif(i.lt.10) then
           lout = lin
           write(4,'(a)') inline3(1:33)
         else
           lout = lin
         endif
         write(3,'(a)') inline(1:lout)
       enddo

 31    continue

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
