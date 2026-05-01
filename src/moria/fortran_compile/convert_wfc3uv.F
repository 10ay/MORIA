 
c
c this differs from _i4 because it
c properly accounts for zero <doh!>
c
 
#define _LINUX_ .true.
#define _PXDIMX_ 4096
#define _PXDIMY_ 4096
 
      program convert_WFC
      implicit none
 
      real      pixr(4096,4096)
      real      pixh(4096,4096)
      integer*2 pixi(4096,4096)
      integer*4 pixii
 
      integer i
      integer j
      integer NARGs, NARG, iargc
 
      character*80 FILEI
      character*80 FILEO
 
      character*70 HDR(25)
      common/HDR/HDR
      data HDR /
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    ',
     . '                    '/
 
      character*80 string_head(720)
 
      NARGs = iargc()
 
 144  format(1x,i4,1x,i4,1x,i7,1x,i7)
c      open(44,file='convert_WFC.SATN',status='unknown')
c      open(45,file='convert_WFC.USAT',status='unknown')
 
      do NARG = 1, NARGs
         call getarg(NARG,FILEI)
         write(*,'(''ENTER      WFCFLTREAD: '',a18)') FILEI
 
         call fitshead2stringe(FILEI,string_head)
 
         call WFC3UV_FLTREAD(FILEI,pixr)
         FILEO = FILEI(01:09) // '_WJ2.fits'
         if (FILEI(11:13).eq.'flc') FILEO = FILEI(01:09) // '_WJC.fits'
         print*,'FILEO: ',FILEO
 
         do i = 0001, 4096
            pixr(i,2046) = -750.0
            pixr(i,2047) = -750.0
            pixr(i,2048) = -750.0
            pixr(i,2049) = -750.0
            pixr(i,2050) = -750.0
            pixr(i,2051) = -750.0
            enddo
 
         do i = 1, 4096
         do j = 1, 4096
            pixii = int(1000+pixr(i,j)+0.5)
            pixii = pixii - 1000
            if (pixii.gt.55000) then
               pixii = 55000 + (pixii-55000+2)/5
               endif
            pixii = pixii - 32000
            if (pixii.lt.-32750) pixii = -32750
            if (pixii.gt. 32750) pixii =  32750
            pixi(i,j) = pixii
            enddo
            enddo
         write(*,'(''ENTER     writfits_i2_32K: '',a18)') FILEO
         call writfits_j2h(FILEO,pixi,4096,4096,string_head)
         print*,' '
         write(*,110) (i,i=0000+64,4096-64,128)
         print*,' '
         do j = 4096-64,2048+64,-128
            write(*,111) j,(pixi(i,j)+32000,i=0000+64,4096-64,128)
            enddo
         print*,' '
         write(*,110) (i,i=0000+64,4096-64,128)
         print*,' '
         do j = 2048-64,0000+64,-128
            write(*,111) j,(pixi(i,j)+32000,i=0000+64,4096-64,128)
            enddo
         print*,' '
         write(*,110) (i,i=0000+64,4096-64,128)
         print*,' '
         print*,' '
  110    format(4x,1x,33i6)
  111    format(i4,1x,33i6)
         enddo
 
      stop
      end
 
 
 
 
c----------------------------------------------
c
c this routine will search for a field in a fits header
c
 
      subroutine query_hdr(filename,FIELDX,streamx)
      implicit none
 
      character filename*80
 
      character*8  field
      character*20 stream
 
      character*8  fieldx
      character*20 streamx
 
      integer  i
      integer ios, k
 
      character*2880 buff
 
      integer nread
 
 
c-----------------------------------------------
 
      open(10,file=FILENAME,status='old',
     .  err=900,recl=2880,form='UNFORMATTED',
     .  access='DIRECT')
 
      streamx = '0'
      i = 0
      nread = 0
 100  continue
      i = i + 1
      read(10,rec=i,iostat=ios) buff
      if (ios.lt.0) goto 900
      do k = 0, 35, 1
         field  = buff(k*80+01:k*80+08)
         stream = buff(k*80+11:k*80+31)
         if (field.eq.fieldx) streamx = stream(1:20)
         if (field.eq.'END     ') goto 101
 109     continue
         enddo
      goto 100
 101  continue
      close(10)
 
      return
 
  900 continue
      print*,' '
      print*,'imginfo.e ERROR EXIT. '
      print*,' '
      print*,'ONE OF THE IMAGES WAS NOT IN STANDARD'
      print*,'HST FITS FORMAT.'
      print*,' '
      write(*,'(''PROBLEM FILE: '',a80)') FILENAME
      print*,' '
      stop
      end
 
 
c---------------------------------------------------
c
c This routine checks to see if the header of an image
c contains a particular keyword. If it does, it returns
c its value-string in STREAMX.
c
 
      subroutine query_hdre(filename,FIELDX,streamx)
      implicit none
 
      character*80 filename
 
      character*8  field
      character*20 stream
 
      character*8  fieldx
      character*20 streamx
 
      integer  i
      integer ios, k
 
      character*2880 buff
 
      integer nread
 
c-----------------------------------------------
 
      !print*,'query_hdr...',FILENAME
      !print*,'    fieldx: ',fieldx
      streamx = ' '
 
      close(10)
      open(10,file=FILENAME,status='old',
     .  err=900,recl=2880,form='UNFORMATTED',
     .  access='DIRECT')
 
      streamx = 'NULL'
      i = 0
      nread = 0
 100  continue
      i = i + 1
      read(10,rec=i,iostat=ios) buff
      if (ios.lt.0) goto 900
      do k = 0, 35, 1
         field  = buff(k*80+01:k*80+08)
         stream = buff(k*80+11:k*80+31)
         !print*,i,k,field,stream
         if (field.eq.fieldx) streamx = stream(1:20)
         if (field.eq.'END     ') goto 101
 109     continue
         enddo
      goto 100
 101  continue
      nread = nread + 1
      if (nread.le.1) goto 100
      close(10)
      !print*,'    streamx: ',streamx
 
      return
 
  900 continue
      print*,' '
      print*,'imginfo.e ERROR EXIT. '
      print*,' '
      print*,'ONE OF THE IMAGES WAS NOT IN STANDARD'
      print*,'HST FITS FORMAT.'
      print*,' '
      write(*,'(''PROBLEM FILE: '',a80)') FILENAME
      print*,' '
      stop
      end
 
 
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/FITSIO/GEN/readfits_r4e.f"
c****                                         
c*********************************************
                                              
 
c-----------------------------------------------------------------
c
c reads an r4 fits image, with extensions (reads in one extension)
c
 
      subroutine readfits_r4e(FILE,pix,NDIMX,NDIMY,NEXTENU)
      implicit none
 
      character*(*) FILE
      integer NDIMX,NDIMY
      real    pix(NDIMX,NDIMY)
      integer NEXTENU
 
      character*199 FILEU
      character*70 INFO(10)
      common / fitsinfo / INFO
 
      integer naxes
      integer laxis(3)
      common/laxis3_/laxis
 
      integer NXF
      integer NYF
 
      character*8  field
      character*40 stream
 
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer i,ios, k
      integer j
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffc,buffb)
 
      integer ifirst, i1, i2
 
      integer np1, np2, npt
      integer nextend
      integer nread
      real bscale, bzero
      integer bitpix
 
      character*70 HDR(25)
      common/HDR/HDR
 
      logical DIAG
      data DIAG /.false./
      integer NEND
      integer iend
 
      if (DIAG) then
         print*,'enter readfits_r4e...... '
         print*,'FILEi: ',FILE(1:60)
         endif
 
      FILEU = FILE
      iend = 0
      do i = 196,1,-1
         if (DIAG) print*,i,iend,FILE(i:i+4)
         if (FILE(i:i+4).eq.'.fits') iend = i+4
         enddo
      if (DIAG) then
         print*,'iend: ',iend
         endif
      if (iend.eq.0) stop 'NO .fits in FILENAME'
      FILEU = FILE(1:iend)
 
      if (DIAG) then
         print*,'purge...'
         print*,'FILEu: ',FILEU
         endif
 
      open(10,file=FILEU,status='old',
     .     err=900,recl=2880,form='UNFORMATTED',
     .     access='DIRECT')
 
      if (DIAG) print*,'...opened ',FILEU
 
      naxes = -1
      laxis(1) = 1
      laxis(2) = 1
      laxis(3) = 1
      nextend = 0
 
      do i = 1, 10
         INFO(i) = ' '
         enddo
 
      BSCALE = 1.0
      BZERO  = 0.0
 
      NEND = 0
      i = 0
      nread = 0
 100  continue
      i = i + 1
      if (DIAG) print*,'READREC: ',i
      read(10,rec=i,iostat=ios) buffc
      do k = 0, 35, 1
         if (DIAG) write(*,'(i4,1x,i4,1x,a80)')
     .                   i,k,buffc(k*80+1:k*80+80)
         field  = buffc(k*80+01:k*80+08)
         stream = buffc(k*80+11:k*80+51)
         if (field.eq.'NAXIS   ') read(stream,*) naxes
         if (field.eq.'NAXIS1  ') read(stream,*) laxis(1)
         if (field.eq.'NAXIS2  ') read(stream,*) laxis(2)
         if (field.eq.'NAXIS3  ') read(stream,*) laxis(3)
         if (field.eq.'NEXTEND ') read(stream,*) nextend
         if (field.eq.'BITPIX  ') read(stream,*) bitpix
         if (field.eq.'BSCALE  ') read(stream,*) bscale
         if (field.eq.'BZERO   ') read(stream,*) bzero
 
         if (field.eq.'EXPTIME ') INFO(1) = stream
         if (field.eq.'FILTER  ') INFO(2) = stream
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
         if (field.eq.'FILTER ')  HDR(22) = stream
         if (field.eq.'FILTER1 ') HDR(22) = stream
         if (field.eq.'FILTER2 ') HDR(23) = stream
         if (field.eq.'VAFACTOR') HDR(24) = stream
 
         if (field.eq.'END     ') then
           !print*,' ---> NEND: ',NEND,NEXTENU
            if (NEND.eq.NEXTENU) goto 101
            NEND = NEND + 1
            endif
         enddo
      goto 100
 101  continue
 
      nread = nread + 1
      if (DIAG) then
         print*,' '
         print*,'----------------------------------------'
         print*,'  NREAD: ',nread
         print*,'NEXTEND: ',nextend
         print*,'  NAXIS: ',naxes
         print*,'  LAXIS: ',laxis(1),laxis(2),laxis(3)
         print*,' BITPIX: ',bitpix
         print*,' BSCALE: ',bscale
         print*,'  BZERO: ',bzero
         print*,'  NDIMX: ',NDIMX
         print*,'  NDIMY: ',NDIMY
         print*,' '
         endif
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      NXF = laxis(1)
      NYF = laxis(2)
 
      nbper  = 4*laxis(1)*laxis(2)
      npt    = laxis(1)*laxis(2)
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
      if (BITPIX.ne.-32) then
         print*,'readfits_r4e...'
         print*,'BITPIX: ',BITPIX
         print*,'prob'
         stop
         endif
 
      do i = i1, i2, 1
         read(10,rec=i,iostat=ios) buffc
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/4 + 1
         np2 = (nbyteE-nbyte1)/4 + 1
c        call buff2pix_r4(buffb,pix,np1,npt)
         call buff2pix_r4_edge(buffb,pix,np1,npt,
     .                         NDIMX,NDIMY,laxis(1),laxis(2))
         if (DIAG) write(*,1115) i,np1,np2,npt
 1115    format(1x,i8,1x,i10,1x,i10,1x,i10)
         enddo
 
      close(10)
      if (DIAG) print*,'...closed ',FILEU
      return
 
 900  continue
      print*,'READFITS ERROR'
      stop
 
      end
 
 
c------------------------------------------------------
c
c
      subroutine buff2pix_r4_edge(buff,pix,n1,nt,
     .                            NXP,NYP,NXF,NYF)
      implicit none
c     character buff(2880)
      byte buff(2880)
      integer NXP,NYP
      real pix(NXP,NYP)
      integer n1,nt
      integer NXF,NYF
 
      real pbuff(720)
      common /sneaky/pbuff
      integer i
      integer npu, nbu
      integer NX, NY
 
      byte b(4)
      real r
      equivalence(r,b)
 
      do i = 1, 720
         npu = n1+i-1
         nbu = (i-1)*4
         if (npu.ge.1.and.npu.le.nt) then
            NX = npu - (npu-1)/NXF*NXF
            NY = 1 + (npu-1)/NXF
            if (NX.le.NXP.and.NY.le.NYP) then
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
               if (npu.ge.1.and.npu.le.nt) pix(NX,NY) = r
               endif
            endif
         enddo
 
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/FITSIO/GEN/readfits_i2e.f"
c****                                         
c*********************************************
                                              
 
c-----------------------------------------------------------
c
c read in an i2 image with extensions
c
      subroutine readfits_i2e(FILEI,pix,NDIMX,NDIMY,NEXTENU)
      implicit none
 
      character*(*) FILEI
      integer   NDIMX,NDIMY
      integer   NEXTENU
      integer*2 pix(NDIMX,NDIMY)
 
      character*199 FILEU
      character*070 INFO(10)
      common / fitsinfo / INFO
 
      integer naxes
      integer laxis(3)
      integer NXF
      integer NYF
 
      character*8  field
      character*40 stream
 
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer i,ios, k
      integer j
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffc,buffb)
 
      integer ifirst, i1, i2
 
      integer np1, np2, npt
      integer nextend
      integer nread
      real bscale, bzero
      integer bitpix
 
      character*70 HDR(25)
      common/HDR/HDR
 
      logical DIAG
      data DIAG /.false./
      integer NEND
      integer ii
 
      FILEU = FILEI
      do i = 195,2,-1
         if (FILEI(i:i+4).eq.'.fits') then
            FILEU = FILEI(1:i+4)
            do ii = i+5, 80
               FILEU(ii:ii) = ' '
               enddo
            endif
         enddo
 
      if (DIAG) then
         print*,'enter readfits_i2e...'
         write(*,'(''FILEI: '',80a)') FILEI
         write(*,'(''FILEU: '',80a)') FILEU
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
 
      BSCALE = 1.0
      BZERO  = 0.0
 
      NEXTEND = 0
      NEND = -1
      i = 0
      nread = 0
 100  continue
      i = i + 1
      if (DIAG) print*,'READREC: ',i
      read(10,rec=i,iostat=ios) buffc
      do k = 0, 35, 1
         if (DIAG) write(*,'(i4,1x,i4,1x,a80)')
     .                   i,k,buffc(k*80+1:k*80+80)
         field  = buffc(k*80+01:k*80+08)
         stream = buffc(k*80+11:k*80+51)
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
         if (field.eq.'DEC_DEG ') INFO(6) = stream
         if (field.eq.'RA_DEG  ') INFO(7) = stream
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
         if (field.eq.'RA_DEG  ') HDR(19) = stream
         if (field.eq.'DEC_DEG ') HDR(20) = stream
         if (field.eq.'PROPOSID') HDR(21) = stream
         if (field.eq.'FILTER1 ') HDR(22) = stream
         if (field.eq.'FILTER2 ') HDR(23) = stream
         if (field.eq.'VAFACTOR') HDR(24) = stream
         if (field.eq.'CCDGAIN ') HDR(25) = stream
         if (field.eq.'END     ') then
            if (NEXTENU.gt.NEXTEND.and.NEXTENU.gt.1) then
               print*,'                     '
               write(*,'(''readfits_i2e: '',80a)') FILEI
               print*,'     NEXTENU.lt.NEXTEND...'
               print*,'     ---> NEXTEND: ',NEXTEND
               print*,'     ---> NEXTENU: ',NEXTENU
               print*,'                          '
               stop
               endif
            NEND = NEND + 1
            if (NEND.ge.1) goto 101
            endif
         enddo
      goto 100
 101  continue
 
      nread = nread + 1
      if (DIAG) then
         print*,' '
         print*,'----------------------------------------'
         print*,'   NREAD: ',nread
         print*,' NEXTEND: ',nextend
         print*,' NEXTENU: ',nextenu
         print*,'   NAXIS: ',naxes
         print*,'   LAXIS: ',laxis(1),laxis(2),laxis(3)
         print*,'  BITPIX: ',bitpix
         print*,'  BSCALE: ',bscale
         print*,'   BZERO: ',bzero
         print*,'   NDIMX: ',NDIMX
         print*,'   NDIMY: ',NDIMY
         print*,' '
         endif
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      NXF = laxis(1)
      NYF = laxis(2)
 
      nbper  = 2*laxis(1)*laxis(2)
      npt    =   laxis(1)*laxis(2)
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
c     print*,'NEND: ',NEND,NEXTENU,BITPIX
      if (NEND.ne.NEXTENU) goto 100
 
      if (BITPIX.ne.16) then
         print*,'BITPIX: ',BITPIX
         print*,'prob'
         stop
         endif
 
      do i = i1, i2, 1
         read(10,rec=i,iostat=ios) buffc
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/2 + 1
         np2 = (nbyteE-nbyte1)/2 + 1
         call buff2pix_i2e(buffb,pix(1,1),np1,npt)
         if (DIAG) write(*,1115) i,np1,np2,npt,i/laxis(1)
 1115    format(1x,i8,1x,i10,1x,i10,1x,i10,1x,i6.6)
         enddo
 
 
      close(10)
      return
 
 900  continue
      print*,'READFITS_I2E ERROR'
      print*,'FILEU: ',FILEU
      stop
 
      end
 
 
 
 
      subroutine buff2pix_i2e(buff,pix,n1,nt)
      implicit none
      byte buff(2880)
      integer*2 pix(*)
      integer n1,nt
 
      byte b(2)
      integer ii
      equivalence(ii,b)
 
      integer i, npu, nbu
 
      do i = 1, 1440
         npu = n1+i-1
         nbu = (i-1)*2
         if (.not.(_LINUX_)) then
            b(1) = buff(nbu+1)
            b(2) = buff(nbu+2)
            endif
         if ((_LINUX_)) then
            b(2) = buff(nbu+1)
            b(1) = buff(nbu+2)
            endif
         if (npu.ge.1.and.npu.le.nt) pix(npu) = ii
         enddo
 
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/FITSIO/GEN/writfits_j2h.f"
c****                                         
c*********************************************
                                              
 
 
 
c
c this is equivalent to i2_32K
c
 
 
      subroutine writfits_j2h(FILE,pix,PXDIMX,PXDIMY,STRING_HEAD)
      implicit none
 
      integer PXDIMX,PXDIMY
      character*80 FILE
      integer*2 pix(PXDIMX,PXDIMY)
      character*80 STRING_HEAD(360)
 
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer i,ios
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffb,buffc)
 
      character*70 TEMP70
      integer ii, page
      integer ifirst, i1, i2
 
      integer np1, np2, npt
 
      character*70 HDR(25)
      common/HDR/HDR
 
      character*80 FILEU
      logical break
 
      FILEU = FILE
      do i = 75,2,-1
         if (FILE(i:i+4).eq.'.fits') FILEU = FILE(1:i+4)
         enddo
 
c     write(*,'(''ENTER   writfits_i2_32K: '',80a)') FILEU
 
      open(10,file=FILEU,
     .     status='unknown',
     .     err   =900,
     .     recl  =2880,
     .     form  ='UNFORMATTED',
     .     access='DIRECT')
 
      write(buffc( 0*80+1: 1*80),'(''SIMPLE  =                    T'')')
      write(buffc( 1*80+1: 2*80),'(''BITPIX  =                   16'')')
      write(buffc( 2*80+1: 3*80),'(''NAXIS   = '',8x,i12)') 2
      write(buffc( 3*80+1: 4*80),'(''NAXIS1  = '',8x,i12)') PXDIMX
      write(buffc( 4*80+1: 5*80),'(''NAXIS2  = '',8x,i12)') PXDIMY
      write(buffc( 5*80+1: 6*80),'(''BSCALE  = '',8x,i12)')  00001
      write(buffc( 6*80+1: 7*80),'(''BZERO   = '',8x,i12)')  32000
      do ii = 07, 35
         write(buffc(ii*80+1:ii*80+80),'(''COMMENT   '',a05)') '     '
         enddo
c     write(buffc(35*80+1:35*80+80),'(''END'')')
      i = 1
      write(10,rec=i,iostat=ios) buffc
 
      break = .false.
      do page = 01, 20
         if (.not.break) then
            do ii = 01, 36
               if (STRING_HEAD((page-1)*36+ii)(1:3).eq.'END')
     .             break = .true.
               if (STRING_HEAD((page-1)*36+ii)(1:6).eq.'SIMPLE'.or.
     .             STRING_HEAD((page-1)*36+ii)(1:6).eq.'BITPIX'.or.
     .             STRING_HEAD((page-1)*36+ii)(1:6).eq.'NAXIS '.or.
     .             STRING_HEAD((page-1)*36+ii)(1:6).eq.'NAXIS1'.or.
     .             STRING_HEAD((page-1)*36+ii)(1:6).eq.'NAXIS2'.or.
     .             STRING_HEAD((page-1)*36+ii)(1:6).eq.'NEXTEN'.or.
     .             STRING_HEAD((page-1)*36+ii)(1:6).eq.'EXTEND'.or.
     .             STRING_HEAD((page-1)*36+ii)(1:6).eq.'BZERO '.or.
     .             STRING_HEAD((page-1)*36+ii)(1:6).eq.'BSCALE'.or.
     .             STRING_HEAD((page-1)*36+ii)(1:6).eq.'GROUPS') then
                   TEMP70 = STRING_HEAD((page-1)*36+ii)(01:70)
                   write(STRING_HEAD((page-1)*36+ii),
     .                   '(''COMMENT   '',a70)') TEMP70
                   endif
               buffc((ii-1)*80+01:(ii-1)*80+80) =
     .               STRING_HEAD((page-1)*36+ii)
               enddo
            i = i + 1
            write(10,rec=i,iostat=ios) buffc
            endif
         enddo
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      nbper  = 2*PXDIMX*PXDIMY
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
         np1 = (nbyte0-nbyte1)/2 + 1
         np2 = (nbyteE-nbyte1)/2 + 1
        !if (i.lt.0010) print*,'i: ',i,i1,i2,np1,np2,nbyte0
        !if (i.gt.2900) print*,'i: ',i,i1,i2,np1,np2,nbyte0
         call pix2buff_j2h(buffb,pix,np1,npt)
         write(10,rec=i,iostat=ios) buffc
         enddo
 
      close(10)
 
      return
 
 900  continue
      print*,'WRITFITS.f ERROR'
      stop
 
      end
 
 
 
      subroutine pix2buff_j2h(buff,pix,n1,nt)
      implicit none
      byte buff(2880)
      integer*2 pix(*)
      integer n1,nt
 
      byte b(2)
      integer*2 ii
      equivalence(ii,b)
 
      integer i, npu, nbu
 
      do i = 1, 1440
         npu = n1+i-1
         nbu = (i-1)*2
         if (npu.ge.1.and.npu.le.nt) ii = pix(npu)
         if (.not.(_LINUX_)) then
            buff(nbu+1) = b(1)
            buff(nbu+2) = b(2)
            endif
         if ((_LINUX_)) then
            buff(nbu+1) = b(2)
            buff(nbu+2) = b(1)
            endif
         enddo
 
      return
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/FITSIO/GEN/fitshead2stringe.f"
c****                                         
c*********************************************
                                              
 
c---------------------------------------------------------
c
c this reads in the zero-header and saves it into a string
c ---> this routine
c
      subroutine fitshead2stringe(filename,string_head)
      implicit none
 
      character*80 filename
 
      character*80 string_head(720)
 
      character*2880 buff
 
      integer n, i, ios, k
      integer NEXTEND
      logical FIRST
 
c-----------------------------------------------
 
      do n = 001, 720
         string_head(n) = '                    '
         enddo
 
      open(10,file=FILENAME,status='old',
     .  err=900,recl=2880,form='UNFORMATTED',
     .  access='DIRECT')
 
      NEXTEND = 0
      i = 0
      n = 0
      FIRST = .true.
 100  continue
      i = i + 1
      read(10,rec=i,iostat=ios) buff
      if (ios.lt.0) goto 900
      do k = 0, 35, 1
         n = n + 1
         if (n.gt.720) stop 'fitshead2stringe:  n.gt.720'
         string_head(n) = buff(k*80+01:k*80+80)
         if (string_head(n)(1:3).eq.'END') goto 101
         if (string_head(n)(1:7).eq.'NEXTEND') then
            read(string_head(n)(11:40),*) NEXTEND
            endif
 109     continue
         enddo
      goto 100
 101  continue
      if (NEXTEND.gt.0.and.FIRST) then
         string_head(n)(1:6) = 'PREEND'
         FIRST = .false.
         goto 100
         endif
      close(10)
 
      return
 
  900 continue
      print*,' '
      print*,'fitshead2string.e ERROR EXIT. '
      print*,' '
      write(*,'(''PROBLEM FILE: '',a80)') FILENAME
      print*,' '
      stop
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "/user/jayander/FORTRAN/FITSIO/WFC3UV/SUBROUTINES/WFC3UV_FLTREAD.f"
c****                                         
c*********************************************
                                              
 
 
c
c check: UVIS-CENTER     --- full array --- ic9r16asq
c        UVIS1           --- udef sub   --- ichl08ujq
c        UVIS2           --- udef sub   --- ichl08ulq
c        UVIS2-M1K1C-SUB --- std sub    --- icdw04ljq_flt.fits
c        UVIS2-C1K1C-SUB --- std sub    --- icck01leq
c        UVIS2-C512C-SUB --- std sub    --- icjd04ezq_flt.fits
c        UVIS2-2K2C-SUB  --- std sub    --- icp601pyq_flt.fits
c        UVIS2-M512C-SUB --- std sub    --- ic8ea1eeq_flt.fits
c
c
 
c--------------------------------------------
c
c
      subroutine WFC3UV_FLTREAD(FILE,pix)
      implicit none
      character*(*) FILE
      real pix(4096,4096)
 
      character*20 NAMEAP
      real pix0512(0512,0512)
      real pix0513(0513,0512)
      real pix1024(1024,1024)
      real pix4096(4096,2051)
      real pix2048(2048,2048)
      real pix2047(2047,2050)
      integer i, j
      character*20 STREAM
 
      integer*2 pix0513i(0513,0512)
 
 
      character*20 STREAM_CENTERA1
      character*20 STREAM_CENTERA2
      character*20 STREAM_SIZAXIS1
      character*20 STREAM_SIZAXIS2
      character*20 STREAM_NAXIS1
      character*20 STREAM_NAXIS2
      character*20 STREAM_BITPIX
      character*20 STREAM_CCDAMP
 
      integer CENTERA1
      integer CENTERA2
      integer SIZAXIS1
      integer SIZAXIS2
      integer NAXIS1
      integer NAXIS2
      integer BITPIX
      character*4 CCDAMP
 
c
c----------------------------------------------------------
c
 
      logical ISSUBARRAY
 
      ISSUBARRAY = .true.
 
      call query_hdre(FILE,'APERTURE',NAMEAP)
      call query_hdre(FILE,'SUBARRAY',STREAM)
      call query_hdre(FILE,'CENTERA1',STREAM_CENTERA1)
      call query_hdre(FILE,'CENTERA2',STREAM_CENTERA2)
      call query_hdre(FILE,'SIZAXIS1',STREAM_SIZAXIS1)
      call query_hdre(FILE,'SIZAXIS2',STREAM_SIZAXIS2)
      call query_hdre(FILE,'NAXIS1  ',STREAM_NAXIS1)
      call query_hdre(FILE,'NAXIS2  ',STREAM_NAXIS2)
      call query_hdre(FILE,'BITPIX  ',STREAM_BITPIX)
      call query_hdre(FILE,'CCDAMP  ',STREAM_CCDAMP)
 
      read(STREAM_CENTERA1,*) CENTERA1
      read(STREAM_CENTERA2,*) CENTERA2
      read(STREAM_SIZAXIS1,*) SIZAXIS1
      read(STREAM_SIZAXIS2,*) SIZAXIS2
      read(STREAM_NAXIS1  ,*) NAXIS1
      read(STREAM_NAXIS2  ,*) NAXIS2
      read(STREAM_BITPIX  ,*) BITPIX
 
      CCDAMP = STREAM_CCDAMP(2:5)
 
      do i = 0001, 4096
      do j = 0001, 4096
         pix(i,j) = -750
         enddo
         enddo
 
 
      if (NAMEAP(01:01).eq.'''') NAMEAP = NAMEAP(02:20)
      do i = 1, 20
         if (NAMEAP(i:i).eq.'''') NAMEAP(i:i) = ' '
         enddo
 
      read(STREAM,*) ISSUBARRAY
 
      write(*,'(''WFC3UV_FLTREAD NAMEAP: '',20a,'' SUB? '',l1,5x,80a)')
     .            NAMEAP,ISSUBARRAY,FILE
 
      if (.not.ISSUBARRAY) then
         if (NAMEAP(01:06).eq.'UVIS  '.or.
     .       NAMEAP(01:06).eq.'UVIS1 '.or.
     .       NAMEAP(01:06).eq.'UVIS2 '.or.
     .       NAMEAP(01:09).eq.'UVIS-FIX '.or.
     .       NAMEAP(01:09).eq.'UVIS1-FIX'.or.
     .       NAMEAP(01:09).eq.'UVIS2-FIX'.or.
     .       NAMEAP(01:09).eq.'G280-REF '.or.
     .       NAMEAP(01:09).eq.'UVIS-QUAD'.or.
     .       NAMEAP(01:11).eq.'UVIS-CENTER'.or.
     .       NAMEAP(01:11).eq.'UVIS-IR-FIX'.or.
     .       NAMEAP(01:12).eq.'UVIS1-IR-FIX'.or.
     .       NAMEAP(01:12).eq.'UVIS2-IR-FIX'.or.
     .       NAMEAP(01:13).eq.'UVIS-QUAD-FIX') then
             call readfits_WFC3(FILE,pix4096,1)
             do i = 0001, 4096
             do j = 0001, 2048
                pix(i,j+0000) = pix4096(i,j)
                enddo
                enddo
             call readfits_WFC3(FILE,pix4096,4)
             do i = 0001, 4096
             do j = 0001, 2048
                pix(i,j+2048) = pix4096(i,j)
                enddo
                enddo
             do i = 0001, 4096
             do j = 2048-2, 2049+2
                pix(i,j) = -750
                enddo
                enddo
             return
             endif
         print*,' '
         print*,'WFC3UV_FLTREAD not yet designed to operate on'
         print*,'aperture: ',NAMEAP
         print*,'subarray: ',ISSUBARRAY
         print*,'AS NON-SUBARRAY...'
         print*,' '
         STOP   'HALT IN WFC3UV_FLTREAD'
         endif
 
      if (NAMEAP(01:15).eq.'UVIS2-C1K1C-SUB') then
          print*,'C1K1Ca: ',FILE
          call readfits_r4e(FILE,pix1024,1024,1024,1)
          do i = 0001, 1024
          do j = 0001, 1024
             pix(i+0000,j+0001) = pix1024(i,j)
             enddo
             enddo
          print*,'C1K1Cb: ',FILE
          return
          endif
 
      if (NAMEAP(01:15).eq.'UVIS2-M1K1C-SUB') then
          call readfits_r4e(FILE,pix1024,1024,1024,1)
          do i = 0001, 1024
          do j = 0001, 1024
             pix(i+1023,j+1027) = pix1024(i,j)
             enddo
             enddo
          return
          endif
 
      if (NAMEAP(01:15).eq.'UVIS2-C512C-SUB') then
          call readfits_r4e(FILE,pix0513 ,0513,0512,1)
          call readfits_i2e(FILE,pix0513i,0513,0512,3)
          do i = 0001, 0513
          do j = 0001, 0512
             pix(i+0000,j+0001) = pix0513(i,j)
             enddo
             enddo
          do i = 0001, 0513
          do j = 0002, 0511
	     if (iand(int(pix0513i(i,j)),256).ne.0) then
                pix(i+0000,j+0001+1) = -750
                pix(i+0000,j+0001  ) = -750
                pix(i+0000,j+0001-1) = -750
                endif
             enddo
             enddo
          return
          endif
 
 
      if (NAMEAP(01:15).eq.'UVIS2-M512C-SUB') then
          call readfits_r4e(FILE,pix0512,0512,0512,1)
          do i = 0001, 0512
          do j = 0001, 0512
             pix(i+1535,j+1539) = pix0512(i,j)
             enddo
             enddo
          return
          endif
 
 
      if (NAMEAP(01:15).eq.'UVIS2-2K2C-SUB ') then
          call readfits_r4e(FILE,pix2047,2047,2050,1)
          do i = 0001, 2047
          do j = 0001, 2050
             pix(i+0000,j+0001) = pix2047(i,j)
             enddo
             enddo
          return
          endif
 
 
      if (NAMEAP(01:15).eq.'UVIS2          ') then
         print*,'UVIS2... user defined subarray...'
         if (CENTERA1.ne.2073) stop 'CENTERA1.ne.2073'
         if (CENTERA2.ne.1027) stop 'CENTERA2.ne.1027'
         if (SIZAXIS1.ne.2048) stop 'SIZAXIS1.ne.2048'
         if (SIZAXIS2.ne.2048) stop 'SIZAXIS2.ne.2048'
         if (CCDAMP.ne.'C   ') stop 'CCDAMP.ne."C   "'
         call readfits_r4e(FILE,pix2048,2048,2048,1)
         do i = 0001, 2048
         do j = 0001, 2048
            pix(i+1024-1,j+0003-1) = pix2048(i,j)
            enddo
            enddo
         return
         endif
 
 
      if (NAMEAP(01:15).eq.'UVIS1          ') then
         print*,'UVIS1... user defined subarray...'
         if (CENTERA1.ne.2073) stop 'CENTERA1.ne.2073'
         if (CENTERA2.ne.1027) stop 'CENTERA2.ne.1027'
         if (SIZAXIS1.ne.2048) stop 'SIZAXIS1.ne.2048'
         if (SIZAXIS2.ne.2048) stop 'SIZAXIS2.ne.2048'
         if (CCDAMP.ne.'A   ') stop 'CCDAMP.ne."A   "'
         call readfits_r4e(FILE,pix2048,2048,2048,1)
         do i = 0001, 2048
         do j = 0001, 2048
            pix(i+1024-1,j+0003+2048-1) = pix2048(i,j)
            enddo
            enddo
         return
         endif
 
 
      print*,' '
      print*,'WFC3UV_FLTREAD not yet designed to operate on'
      print*,'aperture: ',NAMEAP
      print*,'subarray: ',ISSUBARRAY
      print*,'AS A SUBARRAY...'
      print*,' '
 
      stop
      end
 
 
 
c--------------------------------------------
c
c
      subroutine readfits_WFC3(FILE,pix,nimg)
      implicit none
 
      character*80 FILE
      real pix(4096,2051)
 
      character*70 INFO(10)
      common / fitsinfo / INFO
 
      integer nimg
      integer naxes
      integer laxis(3)
 
      character*8  field
      character*20 stream
 
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer i,ios, k
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffc,buffb)
 
      integer ifirst, i1, i2
 
      integer np1, np2, npt
      integer nextend
      integer nread
      real bscale, bzero
      integer bitpix
 
 
      logical LINUX
      data LINUX/.true./
 
 
      character*70 HDR(25)
      common/HDR/HDR
 
      open(10,file=FILE,status='old',
     .     err=900,recl=2880,form='UNFORMATTED',
     .     access='DIRECT')
 
 
      bscale = 1.0
      bzero  = 0.0
 
      naxes = -1
      laxis(1) = 1
      laxis(2) = 1
      laxis(3) = 1
      nextend = 0
 
      do i = 1, 10
         INFO(i) = ' '
         enddo
      do i = 1, 25
         HDR(i) = ' '
         enddo
 
      i = 0
      nread = 0
 100  continue
      i = i + 1
      read(10,rec=i,iostat=ios) buffc
      do k = 0, 35, 1
         field  = buffc(k*80+01:k*80+08)
         stream = buffc(k*80+11:k*80+31)
         if (field.eq.'NAXIS   ') read(stream,*) naxes
         if (field.eq.'NAXIS1  ') read(stream,*) laxis(1)
         if (field.eq.'NAXIS2  ') read(stream,*) laxis(2)
         if (field.eq.'NAXIS3  ') read(stream,*) laxis(3)
         if (field.eq.'NEXTEND ') read(stream,*) nextend
         if (field.eq.'BITPIX  ') read(stream,*) bitpix
         if (field.eq.'BSCALE  ') read(stream,*) bscale
         if (field.eq.'BZERO   ') read(stream,*) bzero
 
         if (field.eq.'EXPTIME ') INFO(01) = stream
         if (field.eq.'FILTER  ') INFO(02) = stream
         if (field.eq.'FILENAME') INFO(03) = stream
         if (field.eq.'DATE-OBS') INFO(04) = stream
         if (field.eq.'TIME-OBS') INFO(05) = stream
         if (field.eq.'DEC_TARG') INFO(06) = stream
         if (field.eq.'RA_TARG ') INFO(07) = stream
         if (field.eq.'PA_V3   ') INFO(08) = stream
         if (field.eq.'PROPOSID') INFO(09) = stream
         if (field.eq.'CCDGAIN ') INFO(10) = stream
 
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
         if (field.eq.'FILTER  ') HDR(22) = stream
         if (field.eq.'FILTER  ') HDR(23) = stream
         if (field.eq.'VAFACTOR') HDR(24) = stream
         if (field.eq.'CCDGAIN ') HDR(25) = stream
         if (field.eq.'END     ') goto 101
         enddo
      goto 100
 101  continue
 
 
      nread = nread + 1
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      if (naxes.eq.0) then    ! maybe multiple images stored as extensions
         if (nextend.eq.0) then
            print*,'THIS IS A NULL IMAGE: '
            print*,'NAXES: ',NAXES
            print*,'NEXND: ',NEXTEND
            stop
            endif
         endif
 
      if (nread.ne.nimg+1) then
         nbper  = abs(BITPIX/8)*laxis(1)*laxis(2)
         if (NAXES.eq.0) nbper = 0
         i = i + 1.0*nbper/2880 + 0.9999
         goto 100
         endif
 
      if (laxis(1).ne.4096.or.
     .    laxis(2).ne.2051) then
          print*,'  laxis1: ',laxis(1)
          print*,'  laxis2: ',laxis(2)
          print*,'    4096: ',4096
          print*,'    2051: ',2051
          stop
          endif
 
      if (naxes.eq.2) then   ! nimg is irrelevant; ignore
         nbper  = 4*laxis(1)*laxis(2)
         npt    = laxis(1)*laxis(2)
         nbyte1 = 1
         nbyte2 = nbper
         i1 = i+1 + nbyte1/2880
         i2 = i+1 + nbyte2/2880
         endif
 
      if (naxes.eq.3) then
         nbper  = 4*laxis(1)*laxis(2)
         nbyte1 = 1 + nbper*(nimg-1)
         nbyte2 =     nbper*(nimg  )
         i1 = i+1 + nbyte1/2880
         i2 = i+1 + nbyte2/2880
         npt = laxis(1)*laxis(2)
         endif
 
      do i = i1, i2, 1
         read(10,rec=i,iostat=ios) buffc
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/4 + 1
         np2 = (nbyteE-nbyte1)/4 + 1
         call buff2pix_r4q(buffb,pix,np1,npt)
         enddo
 
      return
 
 900  continue
      print*,'readfits_WFC: READFITS ERROR'
      stop
 
      end
 
 
      subroutine buff2pix_r4q(buff,pix,n1,nt)
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
 
 
 
