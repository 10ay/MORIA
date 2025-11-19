
#define _LINUX_ .true.

      program q
      implicit none

      real      pixr(1024,1024), pixu
      integer*2 pixi(4096,4096)
      integer   i, j
      character*80 FILENAME
      integer NARG, iargc

      do NARG = 1, iargc()

         call getarg(NARG,FILENAME)
         do i = 0001, 4096
         do j = 0001, 4096
            pixi(i,j) = -32750
            enddo
            enddo
      
         call readfits_r4e(FILENAME,pixr,1024,1024,1)
         do i = 0001, 1024
         do j = 0001, 1021
            pixu = int(pixr(i,j)+1000)-1000
            pixu = pixu - 32000
            if (pixu.lt.-32750) pixu = -32750
            if (pixu.gt. 32750) pixu =  32750
           !pixi(i+1023,j+1027) = pixu
            pixi(i+0000,j+0001) = pixu
            enddo
            enddo

         FILENAME(11:13) = 'WJ2'
         print*,'NARG: FILENAME: ',NARG,' ',FILENAME(1:20)
         call writfits_i2_32K(FILENAME,pixi,4096,4096)
         enddo

      stop
      end 

c*********************************************
c****
c**** #include "/home/jayander/FORTRAN/MAC_OSX/FITS/readfits_r4e.f"
c****
c*********************************************


c-----------------------------------------------------------------
c
c reads an r4 fits image, with extensions (reads in one extension)
c

      subroutine readfits_r4e(FILE,pix,NDIMX,NDIMY,NEXTENU)
      implicit none

      character*80 FILE
      integer NDIMX,NDIMY
      real    pix(NDIMX,NDIMY)
      integer NEXTENU

      character*80 FILEU
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
         print*,'enter readfits_r4...'
         print*,'FILE: ',FILE(1:60)
         endif

      FILEU = FILE
      iend = 0
      do i = 76,1,-1
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
         print*,'FILEU: ',FILEU(1:60)
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
c**** #include "/home/jayander/FORTRAN/MAC_OSX/FITS/writfits_i2_32K.f"
c****
c*********************************************


      subroutine writfits_i2_32K(FILE,pix,PXDIMX,PXDIMY)
      implicit none

      integer PXDIMX,PXDIMY
      character*80 FILE
      integer*2 pix(PXDIMX,PXDIMY)

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

      character*70 HDR(25)
      common/HDR/HDR

      character*80 FILEU

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


      i = 1
      write(buffc( 0*80+1: 1*80),'(''SIMPLE  = T    '')')
      write(buffc( 1*80+1: 2*80),'(''BITPIX  =  16  '')')
      write(buffc( 2*80+1: 3*80),'(''NAXIS   = '',i12)') 2
      write(buffc( 3*80+1: 4*80),'(''NAXIS1  = '',i12)') PXDIMX
      write(buffc( 4*80+1: 5*80),'(''NAXIS2  = '',i12)') PXDIMY
      write(buffc( 5*80+1: 6*80),'(''DATATYPE= '',9a)')
     .                          ' ''INTEGER*2 '''
      write(buffc( 6*80+1: 7*80),'(''DATE    = '',11a)')
     .                          ' ''28/01/00'''
      write(buffc( 7*80+1: 8*80),'(''BSCALE  = '',i12)') 00001
      write(buffc( 8*80+1: 9*80),'(''BZERO   = '',i12)') 32000
      write(buffc(09*80+1:10*80),'(''CRPIX1  = '',a20)') HDR(01)
      write(buffc(10*80+1:11*80),'(''CRPIX2  = '',a20)') HDR(02)
      write(buffc(11*80+1:12*80),'(''CRVAL1  = '',a20)') HDR(03)
      write(buffc(12*80+1:13*80),'(''CRVAL2  = '',a20)') HDR(04)
      write(buffc(13*80+1:14*80),'(''CTYPE1  = '',a20)') HDR(05)
      write(buffc(14*80+1:15*80),'(''CTYPE2  = '',a20)') HDR(06)
      write(buffc(15*80+1:16*80),'(''CD1_1   = '',a20)') HDR(07)
      write(buffc(16*80+1:17*80),'(''CD1_2   = '',a20)') HDR(08)
      write(buffc(17*80+1:18*80),'(''CD2_1   = '',a20)') HDR(09)
      write(buffc(18*80+1:19*80),'(''CD2_2   = '',a20)') HDR(10)
      write(buffc(19*80+1:20*80),'(''ORIENTAT= '',a20)') HDR(11)
      write(buffc(20*80+1:21*80),'(''PA_APER = '',a20)') HDR(12)
      write(buffc(21*80+1:22*80),'(''PA_V3   = '',a20)') HDR(13)
      write(buffc(22*80+1:23*80),'(''DATE-OBS= '',a20)') HDR(14)
      write(buffc(23*80+1:24*80),'(''TIME-OBS= '',a20)') HDR(15)
      write(buffc(24*80+1:25*80),'(''EXPTIME = '',a20)') HDR(16)
      write(buffc(25*80+1:26*80),'(''ROOTNAME= '',a20)') HDR(17)
      write(buffc(26*80+1:27*80),'(''TARGNAME= '',a20)') HDR(18)
      write(buffc(27*80+1:28*80),'(''RA_TARG = '',a20)') HDR(19)
      write(buffc(28*80+1:29*80),'(''DEC_TARG= '',a20)') HDR(20)
      write(buffc(29*80+1:30*80),'(''PROPOSID= '',a20)') HDR(21)
      write(buffc(30*80+1:31*80),'(''FILTER1 = '',a20)') HDR(22)
      write(buffc(31*80+1:32*80),'(''FILTER2 = '',a20)') HDR(23)
      write(buffc(32*80+1:33*80),'(''VAFACTOR= '',a20)') HDR(24)
      write(buffc(33*80+1:34*80),'(''COMMENT   '')')
      write(buffc(34*80+1:35*80),'(''COMMENT   '')')
      write(buffc(35*80+1:36*80),'(''END       '')')

      write(10,rec=i,iostat=ios) buffc

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
         call pix2buff_i2_32K(buffb,pix,np1,npt)
         write(10,rec=i,iostat=ios) buffc
         enddo

      close(10)

      return

 900  continue
      print*,'WRITFITS.f ERROR'
      stop

      end




      subroutine pix2buff_i2_32K(buff,pix,n1,nt)
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

