 
#define _LINUX_ .true.
 
c *****************************************************
c *****************************************************
c *****                                           *****
c *****                                           *****
c *****    img2sam                                *****
c *****                                           *****
c *****                                           *****
c *****************************************************
c *****************************************************
c
c  This routine will take a list of images, and a list
c  of MAT files and will combine the images into a single
c  stacked mosaic image.  The input file will be called
c  "IN.img2sam_UVIS" and will have 0 columns, like:
c
c      01 j90a01czq_WJ2.fits 6 0
c      02 j90a01d0q_WJ2.fits 6 0
c      03 j90a01dfq_WJ2.fits 6 0
c
c  The first column is the image number, the second the
c  4096x4096 int*4 WFC image, and the 3rd column tells
c  us which filter to use in the inverse transformation.
c  The fourth column is a correction to the sky. (It will
c  be subtracted from the exposure before stacking; must
c  be an integer.)
c
c  The way it makes the mosaic image is to go through the
c  destination image pixel by pixel.  It determines which
c  location in each of the NIMs images corresponds to that
c  location, then it lists the pixel values at those points.
c  It takes these NIMs "samples" for the destination-pixel
c  value and does some sigma clipping on them to remove
c  CRs and sometimes even to remove some diffraction spikes
c  and bled columns if enough of the images are taken at
c  different orientations.
c
c  Jay Anderson
c  jay@eeyore.rice.edu
c  Aug 17, 2004 ; updated for WFC3/UV Mar 4, 2011
c
c--------------------------------------------------------------
 
 
 
c------------------------------------
c max number of images to be stacked
c (at this point, the routine is not ready
c  for more than 20)
c
#define _NIMS_ 30
 
c------------------------------------
c
c
 
#define _PIXOUTX_   1250
#define _PIXOUTY_   1250
#define _F_            1
#define _XOFF_      0000
#define _YOFF_      0000
 
c
c
c   XMETA = (_F_)*(XREF-_XOFF_)   ! (_XOFF_,_YOFF_) is the
c   YMETA = (_F_)*(YREF-_YOFF_)   !  REF frame coord of the
c                                 !  lower left image coord (0,0)
c                                 !  (not actually in the image)
c
 
 
 
c------------------------------------
c max number of stars in each mat field
c
#define _NMATMX_ 59999
 
      program img2sam
      implicit none
 
      integer*2 pixnima(4096,4096,_NIMS_)
 
      real      pixr(4096,4096)
      integer*2 pixo(_PIXOUTX_,_PIXOUTY_)
      integer   pixh
 
      integer NIM, NIMs
      character*80 FILENAME(_NIMS_)
      character*80 FILEOUT
      integer      IFILT(_NIMS_)
      real*8 XOFF(_NIMS_), YOFF(_NIMS_)
      real*8 AG(_NIMS_), BG(_NIMS_), CG(_NIMS_), DG(_NIMS_)
      real*8 GA(_NIMS_), GB(_NIMS_), GC(_NIMS_), GD(_NIMS_)
      real*8 Xo(_NIMS_), Yo(_NIMS_), Wo(_NIMS_), Zo(_NIMS_)
      real*8 DELTA
      real*8 xc, yc, xgc, ygc
      real*8 xr, yr
      integer i, j
      integer ii,jj
      integer ir,jr
      real*8  plist(_NIMS_)
      character clist(99999)
      real*8 bar,sig
      integer NL,NLu
      integer NS
      integer n
      character*3 FIX
      real*8 xmap , ymap
      real*8 xmast, ymast
      real hmb, histmode_r4
      real*8 xref, yref
 
      real skyarr(17,17)
 
      real*8 x1(_NMATMX_), y1(_NMATMX_), m1(_NMATMX_)
      real*8 x2(_NMATMX_), y2(_NMATMX_), m2(_NMATMX_)
      integer NMATs, NMATu
      integer NSAT
      integer NCOR
      integer DSKYu, DSKY(_NIMS_)
      real*8  dm(_NMATMX_)
      real    FRAT(_NIMS_)
 
c--------------------------------------------------
c
c read in the basic info:  image names and filter info
c and the transformations from the indiv frame (x2,y2)
c into the master frame (x1,y1)
c
      print*,'MAT INFO: '
      NIMs = 0
      open(10,file='IN.img2sam_wfc3uv',status='old')
  2   continue
      read(10,*,end=1) NIM,FILEOUT,i,DSKYu
      DSKY(NIM) = DSKYu
      FRAT(NIM) = 1.00
      if (NIM.gt._NIMS_) then
         print*,' THIS ROUTINE WAS COMPILED ONLY TO ALLOW FOR '
         print*,' _NIMS_ = ',_NIMS_,' IMAGES.  THERE ARE AT LEAST'
         print*,'  NIMS  = ',NIMS,' IN THIS IN.img2sam FILE'
         stop
         endif
      if (NIM.gt.NIMs) NIMs = NIM
      FILENAME(NIM) = FILEOUT
      IFILT(NIM) = i
      write(FILEOUT,'(''MAT.'',i3.3)') NIM
      open(11,file=FILEOUT,status='old')
      NMATs = 0
  3   continue
      read(11,*,end=4) x1(NMATs+1),y1(NMATs+1),
     .                 x2(NMATs+1),y2(NMATs+1),
     .                 m1(NMATs+1),m2(NMATs+1)
      NMATs = NMATs + 1
      dm(NMATs) = m1(NMATs)-m2(NMATs)
      if (NMATs.gt._NMATMX_) then
         print*,'Increase _NMATMX_ ',(_NMATMX_)
         stop
         endif
      goto 3
  4   continue
      call glob_fit6nrDP(x1,y1,x2,y2,NMATs,
     .                   AG(NIM),BG(NIM),CG(NIM),DG(NIM),
     .                   Xo(NIM),Yo(NIM),Wo(NIM),Zo(NIM))
      DELTA = AG(NIM)*DG(NIM)-CG(NIM)*BG(NIM)
      GA(NIM) =  DG(NIM)/DELTA
      GB(NIM) = -BG(NIM)/DELTA
      GC(NIM) = -CG(NIM)/DELTA
      GD(NIM) =  AG(NIM)/DELTA
      XOFF(NIM) = GA(NIM)*(2048-Wo(NIM))
     .          + GB(NIM)*(2048-Zo(NIM))
     .          + Xo(NIM)-2048
      YOFF(NIM) = GC(NIM)*(2048-Wo(NIM))
     .          + GD(NIM)*(2048-Zo(NIM))
     .          + Yo(NIM)-2048
      call barsig(dm,NMATs,bar,sig,NMATu,clist)
      write(*,113) NIM,NMATs,AG(NIM),BG(NIM),CG(NIM),DG(NIM),
     .             XOFF(NIM),YOFF(NIM),IFILT(NIM),bar,sig,
     .             FILENAME(NIM)
 113  format(i3.3,4x,i5,4x,4f13.7,4x,2f14.2,4x,i3,1x,2f10.3,5x,60a)
 
      FRAT(NIM) = 10**(-bar/2.5)
 
      goto 2
  1   continue
      close(10)
 
 
      print*,' '
      print*,' NIMs: ',NIMs
      print*,' '
 
 
c--------------------------------------------------
c
c read in the fits image associated with each image number
c also correct each for a 8x8 panel model of the sky background,
c which we take as the 10% level  (pixel value at which 10%
c of the pixels in the region are lower)
c
      do NIM = 1, NIMs
         print*,'---> NIM: ',NIM,FILENAME(NIM)
         call readfits_i2r(FILENAME(NIM),pixr,4096,4096)
         do i = 0001, 4096
         do j = 0001, 4096
            if (pixr(i,j).lt.- 25) pixr(i,j) = -9999999
            pixr(i,j) = (pixr(i,j) - DSKY(NIM))*FRAT(NIM)
            if (pixr(i,j).lt.-500) pixr(i,j) = -750
            enddo
            enddo
        !print*,'-->  NIM: ',NIM,pixr(1024,1024),FILENAME(NIM)(1:25)
        !print*,'--> NCOR: ',NCOR
        !print*,'--> NO LOCAL SKY SUB!'
         if (.false.) then
            call find_skyarr(pixr,skyarr)
            do j = 17,01,-1
               write(15,115) j,(skyarr(i,j),i=01,17)
 115           format(i2,1x,17f6.1)
               enddo
            call subt_skyarr(pixr,skyarr)
            endif
         do i = 0001, 4096
         do j = 0001, 4096
            if (pixr(i,j).gt. 64000) pixr(i,j)  =  64000
            if (pixr(i,j).lt.  -750) pixr(i,j)  = -00750
            pixnima(i,j,NIM) = int(pixr(i,j)+1000.5)-1000-32000
            enddo
            enddo
         do i = 0001, 4096
            pixnima(i,2047,NIM) = -32750
            pixnima(i,2048,NIM) = -32750
            pixnima(i,2049,NIM) = -32750
            pixnima(i,2050,NIM) = -32750
            enddo
 
         hmb = histmode_r4(0001,4096,0001,4096,pixr)
         write(*,'(i2.2,1x,f8.2,1x,4f11.7,1x,2f10.2,2x,f10.4,1x,40a)')
     .         NIM,hmb,
     .         AG(NIM),BG(NIM),CG(NIM),DG(NIM),
     .         XOFF(NIM),YOFF(NIM),FRAT(NIM),FILENAME(NIM)(1:40)
         enddo
 
 
 
 
c----------------------------
c
c initiate the output image to zeros
c
      do i = 0001, _PIXOUTX_
      do j = 0001, _PIXOUTY_
         pixo(i,j) = -32000
         enddo
         enddo
c     do i = 0, 4
c        pixo(0001+i,3001) = 32000
c        pixo(3001,0001+i) = 32000
c        pixo(_PIXOUTX_-i,3001) = 32000
c        pixo(3001,_PIXOUTY_-i) = 32000
c        enddo
 
c--------------------------------------------------------
c
c 4096x4096 are the dimensions of the output image;
c sometimes it may be better to make a larger image,
c but that will be complicated as I have been using pixo
c for dual purposes (above).  Below, go through the
c pixels one by one, find out which one
c
      do j = 0001+0000,_PIXOUTY_-0000
      do i = 0001+0000,_PIXOUTX_-0000
         xref = 1.*i/(_F_)+ (_XOFF_)
         yref = 1.*j/(_F_)+ (_YOFF_)
         NL = 0
         NS = 0
         do NIM = 1, NIMs
            xc = 2048.00 + AG(NIM)*(xref-2048.00-XOFF(NIM))
     .                   + BG(NIM)*(yref-2048.00-YOFF(NIM))
            yc = 2048.00 + CG(NIM)*(xref-2048.00-XOFF(NIM))
     .                   + DG(NIM)*(yref-2048.00-YOFF(NIM))
            call wfc3uv_cgj(xc,yc,0,IFILT(NIM),xr,yr)
            ir = int(xr+0.5)  ! pixel location in image # NIM
            jr = int(yr+0.5)  ! (I don't do any interpolation here)
c            if (i.eq.(_PIXOUTX_)/2) then
c                write(*,114) i,j,NIM,ir,jr,NL
c                endif
            if (ir.ge.1.and.ir.le.4096.and.
     .          jr.ge.1.and.jr.le.4096) then
                if (pixnima(ir,jr,NIM).gt.-32700) then
                   NL = NL + 1
                   plist(NL) = pixnima(ir,jr,NIM   ) + 32000.0
                   clist(NL) = ' '
c                  if (i.eq.(_PIXOUTX_)/2) then
c                     write(*,114) i,j,NIM,ir,jr,NL,plist(NL)
c114                  format(10x,3x,2i5,3x,i2,1x,2i5,1x,i2,1x,f10.2)
c                     endif
                   endif
                endif
            enddo
 
         NLu = 0
         call dbubble(plist,NL)
 
         if (NL.le.10) then
             call barsig_sns(plist,NL,bar,sig,NLu)   ! if there are fewer than 10, then use
           else                                      !    specialized statistics...
             call barsig(plist,NL,bar,sig,NLu,clist) ! here's where I avg all the estimates together
           endif
 
         if (NL.lt.1) bar = -500
         FIX = '   '
         if (NL.ne.NLu) FIX = '***'
         if (i.eq.int(_PIXOUTX_/2))
     .   write( *,111) xref,yref,i,j,NLu,NL,int(bar),int(sig),
     .                 FIX,(int(plist(n)),clist(n),n=1,NL)
 111     format(f8.2,1x,f8.2,1x,i4.4,1x,i4.4,1x,
     .          i2,1x,i2,1x,i6,1x,i4,1x,a3,5x,50(i5,a))
         pixo(i,j) = int(bar)-32000
         if (NL.eq.0) pixo(i,j) = -32750
 
 137     continue
         enddo
         if (j.eq.j/250*250) then ! output early version...
            FILEOUT = 'outputq.fits'
            call writfits_i2_32K(FILEOUT,pixo,_PIXOUTX_,_PIXOUTY_)
            print*,'  writfits_i4: ',FILEOUT
            endif
         enddo
 
c---------------------------------------------
c
c output the final mosaicked image
c
      FILEOUT = 'outputq.fits'
      call writfits_i2_32K(FILEOUT,pixo,_PIXOUTX_,_PIXOUTY_)
      print*,'  writfits_i4: ',FILEOUT
 
 
      stop
      end
 
 
c--------------------------------------------
c
c bubble sort a real*8 into ascending order
c
      subroutine dbubble(r1,NTOT)
      implicit none
      real*8 r1(1)
      integer NTOT
      integer n
      real*8 temp
      logical change
 777  continue
      change = .false.
      do n = 1, NTOT-1
         if (r1(n).gt.r1(n+1)) then
            temp = r1(n)
            r1(n) = r1(n+1)
            r1(n+1) = temp
            change = .true.
            endif
         enddo
      if (change) goto 777
      end
 
 
c--------------------------------------------
c
c do an 8x8 region sky subtraction for an image
c
c      subroutine sky_sub_i2(pix)
c      implicit none
c      integer*2 pix(4096,4096)
c
c      real skyarr(8,8)
c      integer ix, iy
c      real histmode_i2
c      real rclip
c      real rx, ry
c      real fx, fy
c      integer i, j
c      real skyval
c
c      print*,' '
c      print*,'COMPUTE SKY-SUBTRACTION ARRAY: '
c      do ix = 01, 08
c      print*,'   ix(1-->8): ',ix
c      do iy = 01, 08
c         skyarr(ix,iy) = histmode_i2(1+(ix-1)*512,ix*512,
c     .                               1+(iy-1)*512,iy*512,pix)
c         enddo
c         enddo
c
c      print*,' '
c      print*,'THE SKY SUB ARRAY '
c      write(*,'(1x,3x,4x,1x,8i8)') (256+(ix-1)*512,ix=1,8)
c      do iy = 8,1,-1
c         write(*,'(i1,3x,i4,1x,8f8.2)') iy,(iy-1)*512+256,
c     .                     (skyarr(ix,iy),ix=1,8)
c         enddo
c      print*,' '
c      print*,'DO THE SKY SUBTRACION...'
c      print*,' '
c      do i = 0001, 4096
c      do j = 0001, 4096
c         rx = rclip(1+(i-256)/512.0,1.0000,7.9999)
c         ry = rclip(1+(j-256)/512.0,1.0000,7.9999)
c         ix = int(rx)
c         iy = int(ry)
c         fx = rx-ix
c         fy = ry-iy
c         skyval = (1-fx)*(1-fy)*skyarr(ix  ,iy  )
c     .          + (1-fx)*( fy )*skyarr(ix  ,iy+1)
c     .          + ( fx )*(1-fy)*skyarr(ix+1,iy  )
c     .          + ( fx )*( fy )*skyarr(ix+1,iy+1)
c         pix(i,j) = pix(i,j) - skyval
c         enddo
c         enddo
c
c      return
c      end
c
c
c-----------------------------------------
c
c
      real function rclip(rval,rlo,rhi)
      implicit none
 
      real rval
      real rlo
      real rhi
 
      rclip = rval
      if (rclip.gt.rhi) rclip = rhi
      if (rclip.lt.rlo) rclip = rlo
      if (.not.(rclip.lt.rhi).and.
     .    .not.(rclip.gt.rlo)) rclip = rhi
 
      return
      end
 
 
 
c-----------------------------------------
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
 
      subroutine glob_fit6nrDP(x1,y1,x2,y2,NUSE,
     .                         A,B,C,D,x1o,y1o,x2o,y2o)
      implicit none
 
      real*8 x1(1), y1(1)
      real*8 x2(1), y2(1)
      integer NUSE
      real*8 A, B, C, D
      real*8 x1o, y1o, x2o, y2o
 
      real*8 sxx, sx, swx, szx
      real*8 syy, sy, swy, szy
      real*8 sw,  sz, sxy
      real*8 dlta
      real*8 dsxx, dsyy, dsxy
      real*8 dswx, dswy, dszx, dszy
 
      integer n
 
      if (NUSE.lt.3) goto 999
 
      x1o = 0
      y1o = 0
      x2o = 0
      y2o = 0
      do n = 1, NUSE
         x1o = x1o + x1(n)
         y1o = y1o + y1(n)
         x2o = x2o + x2(n)
         y2o = y2o + y2(n)
         enddo
 
      x1o = x1o/NUSE
      y1o = y1o/NUSE
      x2o = x2o/NUSE
      y2o = y2o/NUSE
 
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
 
c      print*,'x1o: ',x1o,y1o
c      print*,'x2o: ',x2o,y2o
 
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
 
 
 
c---------------------------------------------------
c
c this is the routine that takes in several (ordered)
c estimates of the pixel value and computes a mode for
c them; it first finds the closest 50% of the points, then
c it goes out to the IQD from this location and takes
c an average of the pixels in this range
c
      subroutine modsig(xlist,NTOT,bar,sig,NUSE,clist)
      implicit none
 
      integer NTOT
      real*8 xlist(NTOT)
      real*8 bar
      real*8 sig
      integer NUSE
      character clist(NTOT)
 
      integer N50
      integer i50, i
      real    d50, d
 
      real    puse
 
      N50 = 0.5*NTOT
      d50 = xlist(1+N50)-xlist(1)
      i50 = 1
      do i = 2, NTOT-N50
         d = xlist(i+N50)-xlist(i)
         if (d.lt.d50) then
            i50 = i
            d50 = d
            endif
         enddo
      bar = (xlist(i50+N50)+xlist(i50))/2.0
      sig = (xlist(i50+N50)-xlist(i50))
 
      puse = 0.
      nuse = 0.
      do i = 1, NTOT
         clist(i) = '#'
         if (abs(xlist(i)-bar).le.2.0*sig) then
            puse = puse + xlist(i)
            nuse = nuse + 1
            clist(i) = ' '
            endif
         enddo
 
      bar = puse/nuse
      return
      end
 
 
c---------------------------------------------------
c
c this is the routine that takes in several estimates
c of the pixel value and computes *some* kind of
c sigma clipped average for them.
c
      subroutine barsig(xlist,NTOT,bar,sig,NUSE,clist)
      implicit none
 
      integer NTOT
      real*8 xlist(NTOT)
      real*8 bar
      real*8 sig
      integer NUSE
      character clist(NTOT)
 
      integer n
      real*8    bsum, ssum
      integer nsum
      integer NIT
 
 
      bar = 0.e0
      sig = 9e9
      do NIT = 01, 12
         bsum = 0.
         ssum = 0.
         nsum = 0.
         do n = 1, NTOT
            clist(n) = '*'
            if (abs(xlist(n)-bar).le.2.25*sig) then
               bsum = bsum + xlist(n)
               ssum = ssum + abs(xlist(n)-bar)
               nsum = nsum + 1
               clist(n) = ' '
               endif
            enddo
         if (nsum.gt.0) bar = bsum/ nsum
         if (nsum.gt.1) sig = ssum/(nsum-1)
         enddo
      NUSE = nsum
      if (nsum.le.1) sig = 0.999
 
      return
      end
 
 
 
 
 
c---------------------------------------------------
c
c this is the routine that takes in several estimates
c of the pixel value and computes *some* kind of
c sigma clipped average for them.
c
      subroutine barsug(xlist,NTOT,bar,sig,NUSE,clist)
      implicit none
 
      integer NTOT
      real*8 xlist(NTOT)
      real*8 bar
      real*8 sig
      integer NUSE
      character clist(NTOT)
 
      integer n
      real*8    bsum, ssum
      integer nsum
      integer NIT
 
      real*8 ylist(10)
      integer NY, NYs
      real*8 bart, sigt
 
      if (NTOT.eq.0) then
         bar = 0.
         sig = 0.
         return
         endif
 
 
 
      NYs = NTOT
      do NY = 1, NYs
         ylist(NY) = xlist(NY)
         enddo
 
      call dbubble(ylist,NYs)
 
 
 5    continue
 
      bar = 0.
      do n = 1, NYs
         bar = bar + ylist(n)
         enddo
      bar = bar/NYs
 
      sig = 0.
      do n = 1, NYs
         sig = sig + abs(ylist(n)-bar)
         enddo
      sig = sig/(NYs-1)
      if (NYs.le.1) sig = 9999
      NUSE = NYs
 
      if (NYs.eq.1) return
 
      if (NYs.eq.5) then
         bart = (ylist(1)+ylist(2)+ylist(3))/3.0
         sigt = min((ylist(3)-ylist(1))/2,
     .              (ylist(4)-ylist(2))/2)
         if (ylist(5).gt.bart+3*(sigt+5)) then
            NYs = NYs - 1
            goto 5
            endif
         bart = (ylist(2)+ylist(3)+ylist(4)+ylist(5))/4.0
         sigt = (ylist(4)-ylist(2))/2
         if (ylist(1).lt.bart-3*(sigt+5).and.ylist(1).lt.0) then
            bar = bart
            sig = sigt
            NUSE = 4
            return
            endif
         endif
      if (NYs.eq.4) then
         bart = (ylist(1)+ylist(2))/2.0
         sigt = min((ylist(3)-ylist(1))/2.,
     .              (ylist(2)-ylist(1)))
         if (ylist(4).gt.bart+3*(sigt+5)) then
            NYs = NYs - 1
            goto 5
            endif
         bart = (ylist(2)+ylist(3)+ylist(4))/3.0
         sigt = (ylist(4)-ylist(2))/2
         if (ylist(1).lt.bart-3*(sigt+5).and.ylist(1).lt.0) then
            bar = bart
            sig = sigt
            NUSE = 3
            return
            endif
         endif
      if (NYs.eq.3) then
         bart = (ylist(1)+ylist(2))/2.0
         sigt = (ylist(2)-ylist(1))
         if (ylist(4).gt.bart+3*(sigt+5)) then
            NYs = NYs - 1
            goto 5
            endif
         bart = (ylist(2)+ylist(3)+ylist(4))/3.0
         sigt = (ylist(4)-ylist(2))/2
         if (ylist(1).lt.bart-3*(sigt+5).and.ylist(1).lt.0) then
            bar = bart
            sig = sigt
            NUSE = 2
            return
            endif
         endif
      if (NYs.eq.2) then
         if (ylist(2).gt.ylist(1)+3*(5+sqrt(ylist(1)))) then
            bar = ylist(1)
            sig = 9999
            NUSE = 1
            endif
         endif
 
 
 
 
      return
      end
 
 
 
 
 
 
      real function root(n,list,NTOT)
      implicit none
      real    n
      integer list(1)
      integer NTOT
      integer nl
 
      do nl = 1, NTOT-2
         if (n.lt.list(nl+1)) goto 3
         enddo
 3    continue
 
      root = nl
      if (list(nl).ne.list(nl+1))
     .    root = nl + 1.*(n-list(nl))/(list(nl+1)-list(nl))
 
      return
      end
 
 
c------------------------------------------------------------
c
c the following code will output the raw pixel values of each
c contributing into the frame of the output image.  this way
c you can look at the many point-samplings that go into the
c master image
c
c
c      do i = 0001,0400
c      do j = 0001,0200
c         pixx(i,j) = 0
c         enddo
c         enddo
c
c      if (.false.) then  ! this little routine
c      i = 1931           ! will print out all the
c      j = 2192           ! exact position observations
c      i = 0553           ! about this (i,j) pixel value
c      j = 1041           ! in the master frame.
c      do NIM = 1, _NIMS_
c         xc = 2048.00 + AG(NIM)*(i-2048.00-XOFF(NIM))
c     .                + BG(NIM)*(j-2048.00-YOFF(NIM))
c         yc = 2048.00 + CG(NIM)*(i-2048.00-XOFF(NIM))
c     .                + DG(NIM)*(j-2048.00-YOFF(NIM))
c         call wfc_cgf(xc,yc,xr,yr,0,8)
c         do ii = xr-01,xr+01
c         do jj = yr-01,yr+01
c            call wfc_gcf(1.0d0*ii,1.0d0*jj,xgc,ygc,2,8)
c            xmap = GA(NIM)*(xgc-2048.00)
c     .           + GB(NIM)*(ygc-2048.00)
c     .           + 2048.00 + XOFF(NIM)
c            ymap = GC(NIM)*(xgc-2048.00)
c     .           + GD(NIM)*(ygc-2048.00)
c     .           + 2048.00 + YOFF(NIM)
c            write( *,33) NIM,xc,yc,xgc,ygc,xmap,ymap,
c     .                      ii,jj,pixnim(ii,jj,NIM)
c            write(33,33) NIM,xc,yc,xgc,ygc,xmap,ymap,
c     .                      ii,jj,pixnim(ii,jj,NIM)
c 33         format(i2,3x,2f10.3,6x,2f10.3,6x,2f10.3,6x,2i5,4x,i6)
c            enddo
c            enddo
c         enddo
c         if (.true.) stop
c         endif
 
 
 
c--------------------------------------------
c
c do an 8x8 region sky subtraction for an image
c
      subroutine find_skyarr(pix,skyarr)
      implicit none
      real pix(4096,4096)
      real skyarr(17,17)
 
      integer ix, iy
      real histmode_r4
      integer i1,i2
      integer j1,j2
      real hmb
 
      hmb =  histmode_r4(0001,4096,0001,4096,pix)
      print*,'---> HMB: ',hmb
      do ix = 1, 17
      do iy = 1, 17
         skyarr(ix,iy) = 0.
         enddo
         enddo
      if (.true.) return
 
      print*,' '
      print*,'COMPUTE SKY-SUBTRACTION ARRAY: '
      print*,'(MODIFIED TO DO ALL CHIP AS ONE...)'
      do iy = 17, 01, -1
      do ix = 01, 17, +1
         i1 = (ix-1.75)*256
         i2 = (ix-0.25)*256
         j1 = (iy-1.75)*256
         j2 = (iy-0.25)*256
         skyarr(ix,iy) = histmode_r4(i1,i2,j1,j2,pix)
         enddo
         write(*,'(25x,i2,3x,17i4)')
     .            iy,(int(skyarr(ix,iy)+0.5),ix=01,17)
         enddo
 
      return
      end
 
 
c--------------------------------------------
c
c do an 8x8 region sky subtraction for an image
c
      subroutine subt_skyarr(pix,skyarr)
      implicit none
      real pix(4096,4096)
      real skyarr(17,17)
 
      real rclip
      integer ix, iy
      real    rx, ry
      real    fx, fy
      integer i, j
      real    skyval
 
      print*,'    DO THE SKY SUBTRACION...'
      do j = 0001, 4096
      do i = 0001, 4096
         rx = rclip(1+i/256.0,1.00000,16.99999)
         ry = rclip(1+j/256.0,1.00000,16.99999)
         ix = int(rx)
         iy = int(ry)
         fx = rx-ix
         fy = ry-iy
         skyval = (1-fx)*(1-fy)*skyarr(ix  ,iy  )
     .          + (1-fx)*( fy )*skyarr(ix  ,iy+1)
     .          + ( fx )*(1-fy)*skyarr(ix+1,iy  )
     .          + ( fx )*( fy )*skyarr(ix+1,iy+1)
         pix(i,j) = pix(i,j) - skyval
         enddo
         enddo
 
      return
      end
 
 
 
c-----------------------------------------
c
c
      real*8 function dclip(dval,rlo,rhi)
      implicit none
 
      real*8 dval
      real   rlo
      real   rhi
 
      dclip = dval
      if (dclip.gt.rhi) dclip = rhi
      if (dclip.lt.rlo) dclip = rlo
      if (.not.(dclip.lt.rhi).and.
     .    .not.(dclip.gt.rlo)) dclip = rhi
 
      return
      end
 
 
 
 
 
c---------------------------------------------------
c
c find the most common value within a region of an
c image; this just finds the value below which 10%
c of the pixels land.
c
      real function histmode_r4(i1,i2,j1,j2,pixarr)
      implicit none
 
      integer   i1,i2,j1,j2
      real      pixarr(4096,4096)
 
      integer i, j, h
      integer mhi, mlo
      real hist_summ
 
      integer htot
      integer hist(500)
      integer hcum(500)
 
      real    tmin
      integer hmin, hh
      integer dmin
      integer dcum
      integer cmin
 
      integer ntot
      real    ptot
 
      do h = 1,500
         hist(h) = 0
         hcum(h) = 0
         enddo
 
      htot = 0
      do i = i1,i2
      do j = j1,j2
         if (i.lt.0001) goto 2
         if (j.lt.0001) goto 2
         if (i.gt.4096) goto 2
         if (j.gt.4096) goto 2
         if (pixarr(i,j).ge.10000) goto 2
         if (pixarr(i,j).le.-0900) goto 2
         h = pixarr(i,j)+100
         if (h.lt.  1) goto 2
         if (h.gt.500) goto 2
         hist(h) = hist(h) + 1
         htot = htot + 1
 2       continue
         enddo
         enddo
 
 
      histmode_r4 = 0
 
      hcum(1) = hist(1)
      do h = 2, 500
         hcum(h) = hcum(h-1) + hist(h)
         if (hcum(h).lt.0.1*htot) histmode_r4 = h-100
         enddo
 
 
      hmin = 1
      dmin = 499
      do h = 1, 350
         do hh = h+1,500
            dcum = hcum(hh)-hcum(h)
            if (dcum.gt.htot*0.33) goto 3
            enddo
  3      continue
         if (hh-h.lt.dmin.or.
     .       hh-h.eq.dmin.and.
     .       (dcum.gt.cmin)) then
            dmin = hh-h
            hmin = h
            cmin = dcum
c            print*,'---> h: ',h,dmin,cmin
            endif
         enddo
 
      ptot = 0.
      ntot = 0
      do h = hmin-1-dmin/3,hmin+dmin+dmin/3+1
         ptot = ptot + hist(h)*(h-100)
         ntot = ntot + hist(h)
c         write(*,'(i4,1x,i8)') h-100,hist(h)
         enddo
      histmode_r4 = ptot/ntot
 
      return
 
      end
 
 
c-----------------------------------------------------
c
c barsig for small number statistics... N<=6
c
c
      subroutine barsig_sns(xlist,NTOTT,bar,sig,NUSE)
      implicit none
 
      integer NTOTT
      real*8  xlist(10)
      real*8  bar
      real*8  sig
      integer NUSE
 
      integer NTOT
      integer n
      real bar1, sig1
      real bar2, sig2
      real bar3, sig3
      real*8 slist(10)
 
      NTOT = NTOTT
 
      if (NTOT.gt.10) then
         print*,'barsig_sns called with NTOT>6: ',NTOT
         stop
         endif
 
      if (NTOT.eq.0) then
         bar = 0.
         sig = 0.
         NUSE = 0
         return
         endif
 
      do n = 1, NTOT
         slist(n) = xlist(n)
         enddo
 
      NTOT = NTOTT
 777  continue
      call dbubble(slist,NTOT)
 
      if (NTOT.eq.10) then
         bar1 = (slist(1)+slist(2)+slist(3)+slist(4)
     .          +slist(5)+slist(6)+slist(7)+slist(8))/8
         sig1 = (slist(7)-slist(2))/2.5
         bar2 = (slist(2)+slist(3)+slist(4)+slist(5)
     .          +slist(6)+slist(7)+slist(8)+slist(9))/8
         sig2 = (slist(8)-slist(3))/2.5
         bar = bar1
         sig = sig1
         if (sig2.lt.sig1) then
            bar = bar2
            sig = sig2
            endif
         if (sig.lt.5) sig = 5.0
         if (sig.lt.0.10*abs(bar)) sig = 0.10*abs(bar)
         if (slist(10).gt.bar+3*sig) then
            NTOT = 9
            goto 777
            endif
         if (slist(01).lt.bar-3*sig) then
            NTOT = 9
            slist(1) = slist(2)
            slist(2) = slist(3)
            slist(3) = slist(4)
            slist(4) = slist(5)
            slist(5) = slist(6)
            slist(6) = slist(7)
            slist(7) = slist(8)
            slist(8) = slist(9)
            slist(9) = slist(10)
            goto 777
            endif
         bar = (slist(1)+slist(2)+slist(3)+slist(4)
     .         +slist(5)+slist(6)+slist(7)+slist(8)
     .         +slist(9)+slist(10))/10
         sig = (abs(slist(01)-bar)
     .         +abs(slist(02)-bar)
     .         +abs(slist(03)-bar)
     .         +abs(slist(04)-bar)
     .         +abs(slist(05)-bar)
     .         +abs(slist(06)-bar)
     .         +abs(slist(07)-bar)
     .         +abs(slist(08)-bar)
     .         +abs(slist(09)-bar)
     .         +abs(slist(10)-bar))/9
         NUSE = 10
         goto 666
         endif
 
      if (NTOT.eq.9) then
         bar1 = (slist(1)+slist(2)+slist(3)+slist(4)
     .          +slist(5)+slist(6)+slist(7))/7
         sig1 = (slist(7)-slist(2))/3.0
         bar2 = (slist(2)+slist(3)+slist(4)+slist(5)
     .          +slist(6)+slist(7)+slist(8))/7
         sig2 = (slist(8)-slist(3))/3.0
         bar = bar1
         sig = sig1
         if (sig2.lt.sig1) then
            bar = bar2
            sig = sig2
            endif
         if (sig.lt.5) sig = 5.0
         if (sig.lt.0.10*abs(bar)) sig = 0.10*abs(bar)
         if (slist(9).gt.bar+3*sig) then
            NTOT = 8
            goto 777
            endif
         if (slist(1).lt.bar-3*sig) then
            NTOT = 8
            slist(1) = slist(2)
            slist(2) = slist(3)
            slist(3) = slist(4)
            slist(4) = slist(5)
            slist(5) = slist(6)
            slist(6) = slist(7)
            slist(7) = slist(8)
            slist(8) = slist(9)
            goto 777
            endif
         bar = (slist(1)+slist(2)+slist(3)+slist(4)
     .         +slist(5)+slist(6)+slist(7)+slist(8)+slist(9))/9
         sig = (abs(slist(1)-bar)
     .         +abs(slist(2)-bar)
     .         +abs(slist(3)-bar)
     .         +abs(slist(4)-bar)
     .         +abs(slist(5)-bar)
     .         +abs(slist(6)-bar)
     .         +abs(slist(7)-bar)
     .         +abs(slist(8)-bar)
     .         +abs(slist(9)-bar))/8
         NUSE = 9
         goto 666
         endif
 
      if (NTOT.eq.8) then
         bar1 = (slist(1)+slist(2)+slist(3)+slist(4)+slist(5))/4
         sig1 = (slist(6)-slist(2))/3.0
         bar2 = (slist(2)+slist(3)+slist(4)+slist(5)+slist(6))/4
         sig2 = (slist(7)-slist(3))/3.0
         bar = bar1
         sig = sig1
         if (sig2.lt.sig1) then
            bar = bar2
            sig = sig2
            endif
         if (sig.lt.5) sig = 5.0
         if (sig.lt.0.10*abs(bar)) sig = 0.10*abs(bar)
         if (slist(8).gt.bar+3*sig) then
            NTOT = 7
            goto 777
            endif
         if (slist(1).lt.bar-3*sig) then
            NTOT = 7
            slist(1) = slist(2)
            slist(2) = slist(3)
            slist(3) = slist(4)
            slist(4) = slist(5)
            slist(5) = slist(6)
            slist(6) = slist(7)
            slist(7) = slist(8)
            goto 777
            endif
         bar = (slist(1)+slist(2)+slist(3)+slist(4)
     .         +slist(5)+slist(6)+slist(7)+slist(8))/8
         sig = (abs(slist(1)-bar)
     .         +abs(slist(2)-bar)
     .         +abs(slist(3)-bar)
     .         +abs(slist(4)-bar)
     .         +abs(slist(5)-bar)
     .         +abs(slist(6)-bar)
     .         +abs(slist(7)-bar)
     .         +abs(slist(8)-bar))/7
         NUSE = 8
         goto 666
         endif
 
      if (NTOT.eq.7) then
         bar1 = (slist(1)+slist(2)+slist(3)+slist(4))/4
         sig1 = (slist(5)-slist(1))/3.0
         bar2 = (slist(2)+slist(3)+slist(4)+slist(5))/4
         sig2 = (slist(6)-slist(2))/3.0
         bar = bar1
         sig = sig1
         if (sig2.lt.sig1) then
            bar = bar2
            sig = sig2
            endif
         if (sig.lt.5) sig = 5.0
         if (sig.lt.0.10*abs(bar)) sig = 0.10*abs(bar)
         if (slist(7).gt.bar+3*sig) then
            NTOT = 6
            goto 777
            endif
         if (slist(1).lt.bar-3*sig) then
            NTOT = 6
            slist(1) = slist(2)
            slist(2) = slist(3)
            slist(3) = slist(4)
            slist(4) = slist(5)
            slist(5) = slist(6)
            slist(6) = slist(7)
            goto 777
            endif
         bar = (slist(1)+slist(2)+slist(3)+
     .          slist(4)+slist(5)+slist(6)+slist(7))/7
         sig = (abs(slist(1)-bar)
     .         +abs(slist(2)-bar)
     .         +abs(slist(3)-bar)
     .         +abs(slist(4)-bar)
     .         +abs(slist(5)-bar)
     .         +abs(slist(6)-bar)
     .         +abs(slist(7)-bar))/6
         NUSE = 7
         goto 666
         endif
 
      if (NTOT.eq.6) then
         bar1 = (slist(1)+slist(2)+slist(3)+slist(4))/4
         sig1 = (slist(4)-slist(1))/2.5
         bar2 = (slist(2)+slist(3)+slist(4)+slist(5))/4
         sig2 = (slist(5)-slist(2))/2.5
         bar = bar1
         sig = sig1
         if (sig2.lt.sig1) then
            bar = bar2
            sig = sig2
            endif
         if (sig.lt.4) sig = 4.0
         if (sig.lt.0.10*abs(bar)) sig = 0.10*abs(bar)
         if (slist(6).gt.bar+3*sig) then
            NTOT = 5
            goto 777
            endif
         if (slist(1).lt.bar-3*sig) then
            NTOT = 5
            slist(1) = slist(2)
            slist(2) = slist(3)
            slist(3) = slist(4)
            slist(4) = slist(5)
            slist(5) = slist(6)
            goto 777
            endif
         bar = (slist(1)+slist(2)+slist(3)+
     .          slist(4)+slist(5)+slist(6))/6
         sig = (abs(slist(1)-bar)
     .         +abs(slist(2)-bar)
     .         +abs(slist(3)-bar)
     .         +abs(slist(4)-bar)
     .         +abs(slist(5)-bar)
     .         +abs(slist(6)-bar))/5
         NUSE = 6
         goto 666
         endif
 
      if (NTOT.eq.5) then
         if (slist(3).gt.35000.and.
     .       slist(2).lt.01000.and.
     .       slist(1).lt.01000) then
               NTOT = 2
               goto 777
               endif
         bar1 = (slist(1)+slist(2)+slist(3)+slist(4))/4
         sig1 = (slist(4)-slist(1))/2
         bar2 = (slist(2)+slist(3)+slist(4)+slist(5))/4
         sig2 = (slist(5)-slist(2))/2
         bar3 = (slist(1)+slist(2)+slist(3))/3
         sig3 = (slist(3)-slist(1))*1.5
         bar = bar1
         sig = sig1
         if (sig2.lt.sig1) then
            bar = bar2
            sig = sig2
            endif
         if (sig3.lt.sig1) then
            bar = bar3
            sig = sig3
            endif
         if (sig.lt.4) sig = 4.0
         if (sig.lt.0.10*abs(bar)) sig = 0.10*abs(bar)
         if (slist(5).gt.bar+3*sig) then
            NTOT = 4
            goto 777
            endif
         if (slist(1).lt.bar-3*sig) then
            NTOT = 4
            slist(1) = slist(2)
            slist(2) = slist(3)
            slist(3) = slist(4)
            slist(4) = slist(5)
            goto 777
            endif
         bar = (slist(1)+slist(2)+slist(3)+
     .          slist(4)+slist(5))/5
         sig = (abs(slist(1)-bar)
     .         +abs(slist(2)-bar)
     .         +abs(slist(3)-bar)
     .         +abs(slist(4)-bar)
     .         +abs(slist(5)-bar))/4
         NUSE = 5
         goto 666
         endif
      if (NTOT.eq.4) then
         if (slist(3).gt.35000.and.
     .       slist(2).lt.01000.and.
     .       slist(1).lt.01000) then
               NTOT = 2
               goto 777
               endif
         bar1 = (slist(1)+slist(2)+slist(3))/3
         sig1 = (slist(3)-slist(1))
         bar2 = (slist(2)+slist(3)+slist(4))/3
         sig2 = (slist(4)-slist(2))
         bar3 = (slist(1)+slist(2))/2
         sig3 = (slist(2)-slist(1))*3
         bar = bar1
         sig = sig1
         if (sig2.lt.sig1) then
            bar = bar2
            sig = sig2
            endif
         if (sig3.lt.sig1) then
            bar = bar3
            sig = sig3
            endif
         if (sig.lt.4) sig = 4.0
         if (sig.lt.0.10*abs(bar)) sig = 0.10*abs(bar)
         if (slist(4).gt.bar+3*sig) then
            NTOT = 3
            goto 777
            endif
         if (slist(1).lt.bar-3*sig) then
            NTOT = 3
            slist(1) = slist(2)
            slist(2) = slist(3)
            slist(3) = slist(4)
            goto 777
            endif
         bar = (slist(1)+slist(2)+slist(3)+slist(4))/4
         sig = (abs(slist(1)-bar)
     .         +abs(slist(2)-bar)
     .         +abs(slist(3)-bar)
     .         +abs(slist(4)-bar))/3
         NUSE = 4
         goto 666
         endif
      if (NTOT.eq.3) then
         bar1 = (slist(1)+slist(2))/2
         sig1 = (slist(2)-slist(1))
         bar2 = (slist(2)+slist(3))/2
         sig2 = (slist(3)-slist(2))
         bar = bar1
         sig = sig1
         if (sig2.lt.sig1) then
            bar = bar2
            sig = sig2
            endif
         if (sig.lt.4) sig = 4.0
         if (sig.lt.0.10*abs(bar)) sig = 0.10*abs(bar)
         if (slist(3).gt.bar+3*sig) then
            NTOT = 2
            goto 777
            endif
         if (slist(1).lt.bar-3*sig) then
            NTOT = 2
            slist(1) = slist(2)
            slist(2) = slist(3)
            goto 777
            endif
         bar = (slist(1)+slist(2)+slist(3))/3
         sig = (abs(slist(1)-bar)
     .         +abs(slist(2)-bar)
     .         +abs(slist(3)-bar))/2
         NUSE = 3
         goto 666
         endif
      if (NTOT.eq.2) then
         bar = slist(1)
         sig = 4
         if (sig.lt.4) sig = 4.0
         if (sig.lt.0.10*abs(bar)) sig = 0.10*abs(bar)
         if (slist(2).gt.bar+5*sig) then
            NUSE = 1
            bar = slist(1)
            sig = 999.9
            goto 666
            endif
         if (abs(slist(2)).lt.20.and.slist(1).lt.-50) then
            NUSE = 1
            bar = slist(2)
            sig = 999.9
            goto 666
            endif
         bar = (slist(1)+slist(2))/2
         sig = (slist(2)-slist(1))/2
         NUSE = 2
         goto 666
         endif
      if (NTOT.eq.1) then
         bar = slist(1)
         sig = 999.9
         NUSE = 1
         goto 666
         endif
 
      print*,' '
      print*,'---> STOP...'
      print*,' '
      print*,'---> SHOULD NOT BE HERE.'
      print*,' '
 
 666  continue
      return
 
      end
 
 
 
                                              
c*********************************************
c****                                         
c**** #include "/home/jayander/FORTRAN/CODE/JENA/wfc3uv_cg.f"
c****                                         
c*********************************************
                                              
c
c#include "JAY_JENA_WFC3_UVIS_gc_OFFICIAL.f"
c
c      program test
c      implicit none
c
c      real*8 xc ,yc
c      real*8 xr ,yr
c      real*8 xrr,yrr
c      real*8 xcc,ycc
c      integer kin
c      integer fin
c
c      kin = 1
c      fin = 6
c
c      xr = 1000.00
c      yr = 4000.00
c
c
c      do xr = 0050, 4050, 100.0
c      do yr = 0050, 4050, 100.0
c         call wfc3uv_gcj_meta(xr ,yr ,  fin,xc, yc )
c         call wfc3uv_cgj     (xc ,yc ,0,fin,xrr,yrr)
c         call wfc3uv_gcj_meta(xrr,yrr,  fin,xcc,ycc)
c         write(*,999) xr, yr, xc, yc,
c     .                xrr,yrr,xcc,ycc
c  999    format(8f10.2)
c         enddo
c         enddo
c
c
c      return
c      end
 
 
c----------------------------------------
c
c this will perform the inverse
c distortion correction...
c
      subroutine wfc3uv_cgj(xc,yc,kin,fin,xr,yr)
      implicit none
      real*8 xc,yc
      real*8 xr,yr
      integer kin, fin
 
      integer k
      real*8 xt,yt
 
      real*8 dx, dy, dd
      integer NIT
      real*8 xr1, yr1
      real*8 xr2, yr2
 
      k = kin
 
 
      xr1 = xc
      yr1 = yc
      NIT = 0
 1    continue
      NIT = NIT + 1
      call wfc3uv_gcj_meta_k(xr1,yr1,1,fin,xt,yt)
      dx = xc-xt
      dy = yc-yt
      dd = sqrt(dx**2+dy**2)
      xr1 = xr1 + (xc-xt)*0.5
      yr1 = yr1 + (yc-yt)*0.5
      if (dd.gt.0.001.and.NIT.lt.99) goto 1
 11   continue
 
      xr2 = xc
      yr2 = yc
      NIT = 0
 2    continue
      NIT = NIT + 1
      call wfc3uv_gcj_meta_k(xr2,yr2-2048,2,fin,xt,yt)
      dx = xc-xt
      dy = yc-yt
      dd = sqrt(dx**2+dy**2)
      xr2 = xr2 + (xc-xt)*0.5
      yr2 = yr2 + (yc-yt)*0.5
c     if (yr2.gt.2048) yr2 = 2048
      if (dd.gt.0.001.and.NIT.lt.99) goto 2
 22   continue
 
      if (k.eq.1)  then
         xr = xr1
         yr = yr1
         return
         endif
      if (k.eq.2)  then
         xr = xr2
         yr = yr2
         return
         endif
      if (k.le.0) then
         if (yr1.gt.0001.and.
     .       yr1.lt.2048) then
             xr = xr1
             yr = yr1
             if (xr.lt.0000.or.
     .           xr.gt.4096.or.
     .           yr.lt.0000.or.
     .           yr.gt.2048) then
                 xr =  -1.00
                 yr =  -1.00
                 endif
             return
             endif
         if (yr2.gt.2048.and.
     .       yr2.lt.4096) then
             xr = xr2
             yr = yr2
             if (xr.lt.0000.or.
     .           xr.gt.4096.or.
     .           yr.lt.2048.or.
     .           yr.gt.4096) then
                 xr =  -1.00
                 yr =  -1.00
                 endif
             return
             endif
         if (yr2.ge.2048.and.yr1.le.0001.and.k.eq.0) then
             xr = -1
             yr = -1
             return
             endif
         xr = -1
         yr = -1
         endif
 
 
 
      return
      end
 
 
 
                                              
c*********************************************
c****                                         
c**** #include "/home/jayander/FORTRAN/CODE/JENA/JAY_JENA_WFC3_UVIS_gc_OFFICIAL.f"
c****                                         
c*********************************************
                                              

c--------------------------------------------------------------
c 
c this file contains three subroutines:
c
c   1) subroutine wfc3uv_gcj(xr,yr,chip,filter,xc,yc)
c      ---> Jay's compact version of Jena's superb correction.
c
c   3) subroutine wfc3uv_gcj_meta_k(xr,yr,chip,filter,xc,yc)
c      ---> Jay's meta-frame wrapper
c
c   2) subroutine wfc3uv_gcj_meta(xr,yr,filter,xc,yc)
c      ---> Jay's meta-frame wrapper for abutted-chip system
c
c
c  JENA: do feel free to name these as you see fit; the "j"
c        at the end of #1 was only to help me distinguish
c        between my compact formulation and your original one
c        so that I could cross-compare them.  It shouldn't
c        be there on the final version, of course.
c--------------------------------------------------------------



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     This subroutine corrects UVIS individual-chip raw
C     positions
C     
C     USAGE: 
C     (xr,yr): raw sigle-chip positions
C     chip:    chip ID (1-2, as in Bellini & Bedin 2009, 
C              i.e. UVIS1=chip2, UVIS2=chip1)          
C     filter:  0 -> F225W
C              1 -> F275W
C              2 -> F336W
C              3 -> F390W
C              4 -> F438W
C              5 -> F555W
C              6 -> F606W
C              7 -> F775W
C              8 -> F814W
C              9 -> F850L 
C     (xc,yc): corrected single-chip positions
C     
C     NOTE: if (LOOKUP .eq. .false.) the lookup_table
C     correction will be skipped. 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C J This is Jay's version of Jena's correction ; it compiles
c   about 100x faster and uses integer*2 lookup tables that
c   preserve only 0.001 pixel precision.  The polynomials
c   also have a similar precision.  All the extra things are
c   now done in real*4, since there's no need for real*8
c   except at the end.
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine wfc3uv_gcj(xr,yr,chip,filter,xc,yc)
      implicit none
      
      integer i,chip,filter,j,jj

      real*8 xr, yr, xt, yt, xc, yc
      real*8 x_cen, y_cen
      real*8 x_tab(38), y_tab(18), y1_tab(18), y2_tab(18)

      real a(9), b(9)

      integer zx     ! zone number along x
      real    xr_zx  ! the start of each x-zone 
      real    rzx    ! the location within an x-zone 

      integer nx, ny ! table indices in x,y
      real    rx, ry ! rough mapping into tables ; not exact
      real    fx, fy ! the fractional offset between table entries

      real dxp, dyp
      real dxt, dyt
      character*1 char
      character*3 lab

c------> LOOKUP: if true -> perform lookup-table correction
C     if false -> jump to polynomial correction
      logical LOOKUP

      data x_tab/0000.00, 
     .      23.50,  192.25,  304.75,  417.25,  529.75,  698.50, 
     .     698.50,  867.25,  979.75, 1092.25, 1204.75, 1373.50, 
     .    1373.50, 1542.25, 1654.75, 1767.25, 1879.75, 2048.50,
     .    2048.50, 2217.25, 2329.75, 2442.25, 2554.75, 2723.50,
     .    2723.50, 2892.25, 3004.75, 3117.25, 3229.75, 3398.50,
     .    3398.50, 3567.25, 3679.75, 3792.25, 3904.75, 4073.50,
     .    4096.00/

      data y1_tab/0000, 170.8125, 284.6875, 398.5625, 512.4375,
     .     626.3125, 740.1875, 911,
     .     911,1082,  1196,  1310,  1424,  1538,  1652,  1766,
     .     1880,  2051/

      data y2_tab/0000,  171, 285,  399,  513,  627,  741,
     .     855,  969,  1140,
     .     1140, 1310.8125, 1424.6875, 1538.5625, 1652.4375,
     .     1766.3125, 1880.1875, 2051/

c----------------------------------------------------------------------
c
c these are Jena's polynomials ; one line for each coord/filter
c
c
      real apoly(9,2,10)
      data apoly /
     .0.0, 0.0, 12.1203,-6.2786, 0.0638, 0.1759,-0.0568, 0.0037, 0.0347,
     .0.0, 0.0, 12.0209,-6.0574, 0.0012, 0.1487, 0.0219, 0.0610, 0.0325,
     .0.0, 0.0, 12.1023,-6.2844, 0.0608, 0.1783,-0.0559, 0.0048, 0.0333,
     .0.0, 0.0, 12.0155,-6.0395, 0.0012, 0.1436, 0.0258, 0.0508, 0.0324,
     .0.0, 0.0, 12.0912,-6.1882, 0.0653,-0.0623,-0.0965, 0.0164, 0.0330,
     .0.0, 0.0, 11.9937,-6.1352, 0.0042,-0.1514, 0.0740, 0.0398, 0.0326,
     .0.0, 0.0, 12.0575,-6.1989, 0.0669, 0.0417,-0.0828, 0.0160, 0.0502,
     .0.0, 0.0, 11.9586,-6.1143,-0.0066, 0.0016, 0.0700, 0.0742, 0.0547,
     .0.0, 0.0, 12.0355,-6.2339, 0.0632, 0.1161,-0.0506, 0.0311, 0.0408,
     .0.0, 0.0, 11.9735,-6.0685, 0.0021, 0.0612, 0.0431, 0.0732, 0.0563,
     .0.0, 0.0, 12.1457,-6.2259, 0.0621, 0.3040,-0.0430, 0.0225, 0.0338,
     .0.0, 0.0, 12.0384,-6.0369,-0.0048, 0.2288,-0.0145, 0.0538, 0.0541,
     .0.0, 0.0, 12.0439,-6.2225, 0.0654, 0.1827,-0.0544, 0.0333, 0.0232,
     .0.0, 0.0, 11.9592,-6.0686, 0.0009, 0.0963, 0.0371, 0.0759, 0.0359,
     .0.0, 0.0, 12.0443,-6.2037, 0.0623, 0.2732,-0.0516, 0.0585, 0.0746,
     .0.0, 0.0, 12.0220,-6.0396, 0.0046, 0.1394, 0.1234, 0.0205, 0.0324,
     .0.0, 0.0, 12.0251,-6.2193, 0.0653, 0.1919,-0.0539, 0.0473, 0.0375,
     .0.0, 0.0, 12.0577,-6.0836, 0.0027, 0.0638, 0.0910, 0.0730, 0.0273,
     .0.0, 0.0, 12.0511,-6.2362, 0.0578, 0.2917,-0.0457, 0.0421, 0.0576,
     .0.0, 0.0, 11.9808,-5.9852,-0.0017, 0.1943, 0.0244, 0.0661, 0.0388/

      real bpoly(9,2,10)
      data bpoly /
     .129.2297, 1.9351,.5906,5.5526,-3.2274,.0291,0.0327,-0.0412,-.0230,
     .140.2701,-4.2153,.7732,5.4956,-3.0580,.1557,-.0093,-0.0264,0.0277,
     .129.2700, 1.9251,.5809,5.5472,-3.2410,.0333,0.0538,-0.0408,-.0118,
     .140.2850,-4.2207,.7807,5.4933,-3.0483,.1633,0.0072,-0.0253,0.0198,
     .129.4380, 1.7865,.6762,5.5652,-3.1553,.0038,0.0343,-0.0615,0.0163,
     .140.3147,-4.3218,.6715,5.4763,-3.1515,.1893,-.0271, 0.0050,0.0140,
     .129.3641, 1.8663,.6637,5.5509,-3.1606,.0137,0.0260,-0.0513,-.0139,
     .140.2162,-4.2584,.7200,5.4724,-3.1158,.1772,-.0166, 0.0014,0.0223,
     .129.3780, 1.8968,.6111,5.5392,-3.1606,.0303,0.0543,-0.0328,-.0219,
     .140.2763,-4.2583,.7713,5.4920,-3.0793,.1672,0.0269,-0.0219,0.0233,
     .129.3566, 1.9305,.5794,5.5224,-3.2361,.0322,0.0474,-0.0101,0.0053,
     .140.2924,-4.1853,.7695,5.4644,-3.0337,.1187,0.0147,-0.0192,0.0178,
     .129.3845, 1.8679,.6216,5.5261,-3.2111,.0259,0.0804,-0.0336,0.0093,
     .140.2800,-4.2506,.7619,5.4616,-3.0626,.1555,0.0086,-0.0196,0.0493,
     .129.3217, 1.9458,.5925,5.5377,-3.1707,.0017,0.0416, 0.0059,-.0007,
     .140.1670,-4.2114,.7804,5.5297,-2.9862,.2284,0.0016, 0.0373,-.0099,
     .129.3650, 1.8970,.6181,5.5162,-3.1900,.0373,0.0716,-0.0474,-.0039,
     .140.1682,-4.2747,.7646,5.5080,-3.0895,.2008,0.0225, 0.0417,0.0535,
     .129.3177, 1.9401,.5467,5.5222,-3.2432,.0243,0.0650,-0.0497,0.0007,
     .140.2988,-4.2084,.8176,5.5099,-2.9918,.1504,0.0267,-0.0299,0.0396/

c----------------------------------------------------------------------
c
c J these are the lookup-table entries for Jena's lookup-table
c   I multiply them by 0.001 to get the value to use.  These
c   have the filter as the fastest-moving element, then x then y ;
c   it probably would have been better to have filter move slowest,
c   but it wasn't convenient.
c
c
      integer*2 dx1(10,38,18)
      integer*2 dy1(10,38,18)
      integer*2 dx2(10,38,18)
      integer*2 dy2(10,38,18)

      data dx1 / 
     .-085,-076,-165,-094,-062,-030,-081,-051,-118,-055, ! dx1(1:10, 1, 1)
     .-021,-027,-147,-111,-038,-030,-061,-035,-104,-043, ! dx1(1:10, 2, 1)
     .-038,-041,-104,-070,-025,-025,-038,-009,-077,-002, ! dx1(1:10, 3, 1)
     .-011,-024,-086,-066,-045,-028,-045,-005,-070,-037, ! dx1(1:10, 4, 1)
     .-064,-044,-072,-045,-047,-046,-042,-013,-062,-026, ! dx1(1:10, 5, 1)
     .-031,-021,-032,-018,-021,-022,-021,-020,-030, 001, ! dx1(1:10, 6, 1)
     .-006,-004,-007,-005,-013,-038,-011,-023, 002,-016, ! dx1(1:10, 7, 1)
     .-056,-067,-061,-046,-061,-082,-061,-070,-041,-065, ! dx1(1:10, 8, 1)
     .-031,-018, 003,-002,-026,-030,-027,-047,-008,-012, ! dx1(1:10, 9, 1)
     .-041,-032,-001,-008,-029,-040,-032,-038,-019,-021, ! dx1(1:10,10, 1)
     .-002,-008, 010,-002,-025,-038,-031,-014,-011,-004, ! dx1(1:10,11, 1)
     .-009, 001, 026, 032,-010,-017,-013, 002, 009, 005, ! dx1(1:10,12, 1)
     .-009,-003, 018, 012,-008,-008,-008, 002, 024, 011, ! dx1(1:10,13, 1)
     . 022, 012, 032, 035, 006,-001,-002, 015, 036, 014, ! dx1(1:10,14, 1)
     . 014, 017, 032, 033, 024, 022, 016, 023, 052, 042, ! dx1(1:10,15, 1)
     . 008, 024, 027, 034, 022, 021, 023, 020, 053, 021, ! dx1(1:10,16, 1)
     . 030, 026, 036, 046, 037, 029, 034, 033, 053, 040, ! dx1(1:10,17, 1)
     . 041, 037, 044, 053, 053, 052, 037, 039, 053, 055, ! dx1(1:10,18, 1)
     . 047, 052, 047, 060, 056, 053, 052, 044, 038, 054, ! dx1(1:10,19, 1)
     . 020, 022, 027, 035, 035, 024, 023, 021, 020, 035, ! dx1(1:10,20, 1)
     . 045, 046, 017, 035, 048, 054, 033, 039, 017, 053, ! dx1(1:10,21, 1)
     . 005, 031, 015, 031, 034, 029, 016, 044, 004, 038, ! dx1(1:10,22, 1)
     . 029, 021, 016, 028, 036, 036, 023, 058, 002, 043, ! dx1(1:10,23, 1)
     . 021, 018, 000, 011, 025, 025, 017, 068, 003, 042, ! dx1(1:10,24, 1)
     . 000, 018,-003, 015, 020, 018, 008, 053, 012, 031, ! dx1(1:10,25, 1)
     . 000, 007,-031,-003, 003,-009,-018, 037,-021, 003, ! dx1(1:10,26, 1)
     .-008,-001,-037,-020,-013,-027,-020, 029,-015, 012, ! dx1(1:10,27, 1)
     . 005,-008,-022,-017,-006,-030,-018, 016,-020, 017, ! dx1(1:10,28, 1)
     .-001,-013,-029,-023,-011,-017,-028, 004,-025, 025, ! dx1(1:10,29, 1)
     .-011,-015,-017,-017,-013,-009,-022, 001,-031, 022, ! dx1(1:10,30, 1)
     .-007,-012,-010,-013,-002, 031,-020,-011,-014, 003, ! dx1(1:10,31, 1)
     .-007,-015,-031,-017,-012, 033,-024,-022,-020, 017, ! dx1(1:10,32, 1)
     . 015, 011, 012, 011, 008, 050,-015, 010, 014, 009, ! dx1(1:10,33, 1)
     . 007, 002, 019, 009, 002, 033,-021, 010, 024, 020, ! dx1(1:10,34, 1)
     . 004, 019, 036, 033, 018, 036,-025, 022, 034, 014, ! dx1(1:10,35, 1)
     . 006, 004, 058, 037, 007, 001,-032, 025, 023,-006, ! dx1(1:10,36, 1)
     . 013, 002, 094, 059, 013,-087,-034, 022,-013,-012, ! dx1(1:10,37, 1)
     . 060, 032, 121, 079, 030,-073,-010, 037, 014, 014, ! dx1(1:10,38, 1)
     .-044,-030,-120,-077,-016,-016,-035,-027,-067,-023, ! dx1(1:10, 1, 2)
     .-057,-040,-119,-083,-025,-015,-042,-022,-062,-015, ! dx1(1:10, 2, 2)
     .-041,-019,-064,-057,-020, 000,-027,-002,-050,-032, ! dx1(1:10, 3, 2)
     .-039,-034,-061,-050,-032,-025,-028,-012,-050,-027, ! dx1(1:10, 4, 2)
     .-031,-017,-040,-036,-020,-024,-026,-010,-036,-019, ! dx1(1:10, 5, 2)
     .-002,-006,-008,-007,-018,-026,-003,-013,-012,-011, ! dx1(1:10, 6, 2)
     .-023, 000, 013, 004, 000,-010, 006,-020, 004,-003, ! dx1(1:10, 7, 2)
     .-053,-027,-029,-033,-049,-052,-041,-066,-038,-046, ! dx1(1:10, 8, 2)
     .-025,-023, 001,-005,-028,-039,-023,-040,-014,-031, ! dx1(1:10, 9, 2)
     .-012,-016, 007,-007,-029,-031,-026,-038,-014,-026, ! dx1(1:10,10, 2)
     .-014,-008, 015,-001,-020,-034,-020,-025,-013,-031, ! dx1(1:10,11, 2)
     .-007,-005, 017, 000,-018,-022,-018,-018,-006,-003, ! dx1(1:10,12, 2)
     .-015,-001, 011, 002,-010,-021,-010,-021, 007,-020, ! dx1(1:10,13, 2)
     . 002, 016, 025, 016, 000,-005, 001,-008, 015, 000, ! dx1(1:10,14, 2)
     . 002, 012, 023, 023, 010, 003, 009,-007, 033, 000, ! dx1(1:10,15, 2)
     . 002, 010, 023, 023, 016, 009, 013,-001, 036, 016, ! dx1(1:10,16, 2)
     . 006, 016, 016, 033, 028, 023, 019, 002, 037, 019, ! dx1(1:10,17, 2)
     . 012, 022, 030, 036, 041, 031, 038, 011, 033, 030, ! dx1(1:10,18, 2)
     . 023, 019, 014, 024, 044, 040, 028, 019, 026, 025, ! dx1(1:10,19, 2)
     . 013, 000,-003, 009, 022, 027, 017, 001, 002, 004, ! dx1(1:10,20, 2)
     . 007, 011, 005, 018, 025, 026, 021, 014, 003, 023, ! dx1(1:10,21, 2)
     . 024, 022, 004, 026, 042, 036, 032, 026, 007, 024, ! dx1(1:10,22, 2)
     . 012, 033,-004, 022, 024, 016, 012, 035, 000, 022, ! dx1(1:10,23, 2)
     . 024, 034, 011, 010, 025, 021, 024, 046, 009, 023, ! dx1(1:10,24, 2)
     . 011, 002,-008, 003, 007,-006, 004, 031, 000, 018, ! dx1(1:10,25, 2)
     .-027,-018,-023,-020,-011,-017,-011, 015,-019,-003, ! dx1(1:10,26, 2)
     .-005,-016,-030,-027,-020,-026,-013, 006,-020, 002, ! dx1(1:10,27, 2)
     .-022,-012,-046,-038,-019,-039,-019,-009,-027,-006, ! dx1(1:10,28, 2)
     .-011,-009,-036,-035,-021,-036,-018,-013,-038,-003, ! dx1(1:10,29, 2)
     .-011,-011,-035,-031,-010,-008,-018,-015,-027, 002, ! dx1(1:10,30, 2)
     .-023,-018,-027,-028,-018, 033,-021,-025,-009, 001, ! dx1(1:10,31, 2)
     .-024,-029,-026,-031,-021, 032,-029,-037,-025,-003, ! dx1(1:10,32, 2)
     . 003,-003,-007,-015,-005, 041,-020,-014, 011, 010, ! dx1(1:10,33, 2)
     .-007, 001, 006, 005,-001, 038,-016,-006, 028, 000, ! dx1(1:10,34, 2)
     .-001,-003, 029, 019, 004, 018,-020, 001, 029,-004, ! dx1(1:10,35, 2)
     . 011, 000, 035, 027, 005, 002,-022, 002, 018,-007, ! dx1(1:10,36, 2)
     .-016, 000, 078, 039, 009,-088,-030, 004, 003,-020, ! dx1(1:10,37, 2)
     . 014, 031, 110, 079, 035,-072,-013, 033, 014,-014, ! dx1(1:10,38, 2)
     .-044,-046,-124,-066,-019,-019,-034,-028,-042,-044, ! dx1(1:10, 1, 3)
     .-023,-006,-082,-069,-010, 002,-030,-014,-045,-029, ! dx1(1:10, 2, 3)
     .-024,-028,-065,-043,-011,-011,-017, 000,-034,-013, ! dx1(1:10, 3, 3)
     .-022,-015,-039,-037,-019,-017,-022,-013,-033,-020, ! dx1(1:10, 4, 3)
     .-011,-010,-021,-022,-018,-019,-013,-017,-027,-022, ! dx1(1:10, 5, 3)
     .-002, 005, 000, 000,-005,-012, 002,-017, 000,-009, ! dx1(1:10, 6, 3)
     .-006,-002, 020, 003,-003,-021, 005,-022,-003,-022, ! dx1(1:10, 7, 3)
     .-013,-040,-024,-031,-038,-050,-030,-061,-043,-059, ! dx1(1:10, 8, 3)
     .-032,-013, 014,-006,-020,-035,-015,-034,-011,-029, ! dx1(1:10, 9, 3)
     .-035,-024, 010,-009,-024,-033,-024,-031,-015,-029, ! dx1(1:10,10, 3)
     .-015,-007, 017, 002,-023,-028,-016,-024,-008,-019, ! dx1(1:10,11, 3)
     .-023,-009, 021, 001,-020,-021,-011,-020,-009,-020, ! dx1(1:10,12, 3)
     .-010,-004,-002,-013,-023,-026,-020,-039,-012,-014, ! dx1(1:10,13, 3)
     . 006, 003, 014, 008,-003,-003,-001,-021, 007,-005, ! dx1(1:10,14, 3)
     . 006, 010, 020, 010, 011, 004, 011,-016, 021, 001, ! dx1(1:10,15, 3)
     .-004, 003, 014, 015, 014, 005, 016,-009, 026, 000, ! dx1(1:10,16, 3)
     . 004, 008, 016, 020, 026, 016, 022,-006, 028, 009, ! dx1(1:10,17, 3)
     . 006, 004, 006, 016, 027, 018, 023, 001, 023, 005, ! dx1(1:10,18, 3)
     . 001, 020, 010, 012, 026, 026, 027, 009, 010, 020, ! dx1(1:10,19, 3)
     .-013, 004,-013, 000, 017, 018, 013,-007,-001, 012, ! dx1(1:10,20, 3)
     . 012, 000, 000, 011, 023, 023, 013, 006, 006, 001, ! dx1(1:10,21, 3)
     . 008, 015,-008, 012, 020, 017, 017, 008, 001, 014, ! dx1(1:10,22, 3)
     . 018, 012,-001, 013, 026, 022, 012, 018, 002, 015, ! dx1(1:10,23, 3)
     . 013, 022,-004, 006, 022, 013, 012, 025,-001, 012, ! dx1(1:10,24, 3)
     . 012, 009, 004,-003, 000,-002, 007, 019,-002, 011, ! dx1(1:10,25, 3)
     .-009, 000,-020,-016,-007,-020,-015, 002,-013, 003, ! dx1(1:10,26, 3)
     .-006,-020,-029,-030,-012,-033,-010,-008,-025,-005, ! dx1(1:10,27, 3)
     .-013,-014,-043,-029,-026,-042,-020,-015,-032,-004, ! dx1(1:10,28, 3)
     .-011,-011,-038,-035,-016,-031,-022,-020,-040, 002, ! dx1(1:10,29, 3)
     . 004,-018,-031,-028,-013,-003,-006,-010,-020, 010, ! dx1(1:10,30, 3)
     .-007,-016,-025,-022,-007, 044,-021,-024,-010, 001, ! dx1(1:10,31, 3)
     .-014,-017,-039,-025,-011, 038,-028,-024,-015, 003, ! dx1(1:10,32, 3)
     .-005,-014,-003,-015,-002, 052,-014,-013, 011, 001, ! dx1(1:10,33, 3)
     . 003,-001, 000, 000,-004, 031,-015,-010, 027, 000, ! dx1(1:10,34, 3)
     . 001, 004, 014, 003, 000, 025,-013,-002, 031,-006, ! dx1(1:10,35, 3)
     .-002, 000, 037, 027, 009, 002,-017,-002, 025,-007, ! dx1(1:10,36, 3)
     . 029, 017, 068, 040, 007,-076,-014, 002, 006,-025, ! dx1(1:10,37, 3)
     . 031, 018, 105, 061, 036,-079, 004, 017, 012, 001, ! dx1(1:10,38, 3)
     .-025,-023,-075,-057,-005, 000,-023,-012,-043, 007, ! dx1(1:10, 1, 4)
     .-031,-007,-072,-037, 004, 001,-011,-009,-027, 000, ! dx1(1:10, 2, 4)
     .-007,-001,-033,-033,-006,-007,-020,-009,-024,-014, ! dx1(1:10, 3, 4)
     .-004,-001,-018,-018,-005,-007,-004,-003,-017,-007, ! dx1(1:10, 4, 4)
     .-014, 000,-010,-001,-005,-013,-003,-018,-007,-012, ! dx1(1:10, 5, 4)
     . 003, 016, 012, 004, 010,-012, 008,-010, 001,-012, ! dx1(1:10, 6, 4)
     .-013, 007, 017,-001, 000,-010, 007,-023, 000,-020, ! dx1(1:10, 7, 4)
     .-052,-021,-010,-017,-033,-045,-023,-053,-025,-039, ! dx1(1:10, 8, 4)
     .-016,-010, 001,-010,-025,-031,-016,-036,-017,-033, ! dx1(1:10, 9, 4)
     .-019, 000, 021, 000,-018,-016,-009,-018,-008,-023, ! dx1(1:10,10, 4)
     .-023,-007, 008,-007,-024,-025,-017,-024,-017,-025, ! dx1(1:10,11, 4)
     . 006, 010, 022, 010,-007,-008, 001,-014, 001,-015, ! dx1(1:10,12, 4)
     .-044,-018,-006,-016,-030,-030,-023,-046,-022,-033, ! dx1(1:10,13, 4)
     .-004, 001, 016, 003,-012,-006, 000,-022, 001,-012, ! dx1(1:10,14, 4)
     . 000, 013, 013, 015, 016, 007, 014,-017, 017,-003, ! dx1(1:10,15, 4)
     . 006,-008,-002, 011, 004, 000, 007,-026, 017,-005, ! dx1(1:10,16, 4)
     .-008, 012, 005, 013, 017, 004, 016,-011, 027,-006, ! dx1(1:10,17, 4)
     .-009, 000, 007, 010, 026, 018, 024, 007, 023, 008, ! dx1(1:10,18, 4)
     .-018, 000, 002, 006, 026, 016, 021, 007, 010,-004, ! dx1(1:10,19, 4)
     .-004,-004,-008, 000, 018, 008, 012,-004,-011,-012, ! dx1(1:10,20, 4)
     . 001, 011, 000, 011, 021, 017, 015, 002, 009, 015, ! dx1(1:10,21, 4)
     . 002, 013,-005, 005, 018, 016, 008,-002, 000, 001, ! dx1(1:10,22, 4)
     . 014, 017,-007, 003, 015, 011, 014, 000, 002, 004, ! dx1(1:10,23, 4)
     . 016, 013, 009, 009, 014, 012, 012, 009, 002, 012, ! dx1(1:10,24, 4)
     . 012, 009,-016,-008,-003,-004, 000, 001,-001, 001, ! dx1(1:10,25, 4)
     .-002,-001,-019,-015,-017,-016,-011,-024,-022,-008, ! dx1(1:10,26, 4)
     .-006,-011,-029,-026,-016,-036,-015,-009,-021, 002, ! dx1(1:10,27, 4)
     .-011,-012,-040,-034,-027,-041,-019,-022,-039,-008, ! dx1(1:10,28, 4)
     . 001,-016,-042,-032,-022,-026,-019,-021,-028,-003, ! dx1(1:10,29, 4)
     .-005, 005,-029,-026, 002, 005,-007,-002,-017, 004, ! dx1(1:10,30, 4)
     .-002,-029,-034,-045,-027, 030,-019,-019,-014,-011, ! dx1(1:10,31, 4)
     .-009,-024,-034,-040,-021, 035,-026,-015,-018,-005, ! dx1(1:10,32, 4)
     .-007,-011,-012,-021,-007, 045,-013,-004, 012,-002, ! dx1(1:10,33, 4)
     . 004,-013,-002,-016,-009, 035,-016,-004, 021,-006, ! dx1(1:10,34, 4)
     . 005,-006, 011, 006,-002, 020,-019,-006, 024,-008, ! dx1(1:10,35, 4)
     .-005, 009, 023, 020, 004, 014,-007,-006, 024,-007, ! dx1(1:10,36, 4)
     . 010,-001, 064, 021, 004,-087,-024,-024, 000,-007, ! dx1(1:10,37, 4)
     . 055, 034, 090, 052, 028,-067, 005, 006, 030, 002, ! dx1(1:10,38, 4)
     .-025, 000,-056,-040, 018,-021,-020,-006,-009,-008, ! dx1(1:10, 1, 5)
     .-027, 003,-047,-017, 016, 007,-004,-003,-013, 002, ! dx1(1:10, 2, 5)
     .-006, 005,-009,-003, 010, 005, 000,-006,-009, 000, ! dx1(1:10, 3, 5)
     . 009, 012,-006, 005, 004, 005, 000,-006,-005, 004, ! dx1(1:10, 4, 5)
     . 001, 019, 006,-001,-003,-010, 000,-016,-007,-014, ! dx1(1:10, 5, 5)
     . 024, 027, 028, 019, 012, 001, 019,-012, 011,-002, ! dx1(1:10, 6, 5)
     . 008, 024, 031, 030, 009, 001, 018,-006, 014, 001, ! dx1(1:10, 7, 5)
     .-028,-014,-009,-017,-035,-042,-024,-058,-035,-054, ! dx1(1:10, 8, 5)
     .-017,-002, 017, 004,-009,-020, 003,-016,-002,-019, ! dx1(1:10, 9, 5)
     .-004,-004, 013,-004,-020,-022,-012,-015,-014,-021, ! dx1(1:10,10, 5)
     .-005,-004, 015, 000,-021,-022,-008,-018,-014,-022, ! dx1(1:10,11, 5)
     .-003, 003, 024, 002,-013,-001, 000,-015,-003,-010, ! dx1(1:10,12, 5)
     .-018,-016, 000,-012,-020,-023,-012,-038,-024,-018, ! dx1(1:10,13, 5)
     .-002, 018, 020, 011, 000, 006, 012,-020, 003,-001, ! dx1(1:10,14, 5)
     . 013, 003, 011, 007, 007, 003, 010,-014, 014,-004, ! dx1(1:10,15, 5)
     .-005, 009,-003, 002, 007,-003, 007,-018, 019,-020, ! dx1(1:10,16, 5)
     .-011,-013,-006, 002, 010, 007, 015,-004, 017,-001, ! dx1(1:10,17, 5)
     .-009,-004, 000, 005, 027, 014, 018, 015, 023,-002, ! dx1(1:10,18, 5)
     . 012, 004,-004, 003, 016, 015, 019, 014, 014,-003, ! dx1(1:10,19, 5)
     .-016,-010,-008,-008, 009, 009, 004,-006, 003, 003, ! dx1(1:10,20, 5)
     . 001, 008,-002, 002, 015, 012, 016, 000, 005,-001, ! dx1(1:10,21, 5)
     .-005, 006,-004, 007, 012, 010, 008,-001, 001, 005, ! dx1(1:10,22, 5)
     . 005, 015, 000,-002, 013, 009, 009,-010, 005, 009, ! dx1(1:10,23, 5)
     . 009, 012, 002, 006, 007, 006, 010,-007, 004, 007, ! dx1(1:10,24, 5)
     .-004, 000,-006, 000,-004,-005,-003,-010,-009, 003, ! dx1(1:10,25, 5)
     .-021,-016,-027,-029,-029,-031,-020,-031,-018,-015, ! dx1(1:10,26, 5)
     .-008,-005,-018,-015,-014,-030,-010,-018,-021,-007, ! dx1(1:10,27, 5)
     .-008,-011,-033,-032,-022,-040,-017,-015,-035,-002, ! dx1(1:10,28, 5)
     .-016,-009,-031,-024,-018,-020,-012,-004,-026, 004, ! dx1(1:10,29, 5)
     . 003,-004,-024,-016,-002, 008, 003, 005,-017, 011, ! dx1(1:10,30, 5)
     .-012,-006,-032,-028,-008, 038,-026,-002,-013,-006, ! dx1(1:10,31, 5)
     .-017,-014,-035,-035,-016, 041,-012,-002,-018,-003, ! dx1(1:10,32, 5)
     . 013, 001,-005,-026,-002, 051,-008, 009, 016, 002, ! dx1(1:10,33, 5)
     . 006,-007,-016,-020,-010, 040,-015,-002, 006,-008, ! dx1(1:10,34, 5)
     . 000, 002, 007,-007,-008, 021,-012,-007, 020,-010, ! dx1(1:10,35, 5)
     . 019, 008, 028, 004, 001, 006,-007,-007, 020,-005, ! dx1(1:10,36, 5)
     .-003,-005, 044, 016, 000,-088,-011,-019,-003,-018, ! dx1(1:10,37, 5)
     . 010, 026, 075, 062, 020,-077, 027, 014, 027, 011, ! dx1(1:10,38, 5)
     .-011,-029,-052,-018,-011,-006,-018, 008,-016, 000, ! dx1(1:10, 1, 6)
     .-008, 006,-029,-002, 017, 015, 002, 014,-004, 005, ! dx1(1:10, 2, 6)
     . 021, 021, 006, 010, 019, 014, 004, 002, 006, 016, ! dx1(1:10, 3, 6)
     . 004, 018, 022, 007, 016, 002, 003,-004, 002, 008, ! dx1(1:10, 4, 6)
     . 008, 013, 012, 015, 002, 000, 005,-022, 006,-009, ! dx1(1:10, 5, 6)
     . 006, 026, 031, 021, 015, 003, 015,-007, 010, 005, ! dx1(1:10, 6, 6)
     . 033, 042, 042, 027, 014, 009, 023,-001, 013,-004, ! dx1(1:10, 7, 6)
     .-016, 003,-001,-001,-020,-024,-008,-038,-018,-023, ! dx1(1:10, 8, 6)
     .-007, 002, 015, 003,-008,-012, 000,-012,-008,-023, ! dx1(1:10, 9, 6)
     .-012,-012, 003,-005,-017,-013,-001,-009,-011,-024, ! dx1(1:10,10, 6)
     .-005, 005, 013, 000,-012,-008, 001,-011,-001,-014, ! dx1(1:10,11, 6)
     . 002,-006, 011,-001,-017,-005, 000,-019,-009,-014, ! dx1(1:10,12, 6)
     .-012,-002,-002,-010,-017,-014,-008,-036,-009,-026, ! dx1(1:10,13, 6)
     .-007, 001, 016, 007,-006, 002, 007,-024, 001,-012, ! dx1(1:10,14, 6)
     . 001, 016, 008, 006, 004, 005, 011,-013, 008,-003, ! dx1(1:10,15, 6)
     .-005,-007,-009,-002,-002, 000, 005,-002, 017,-018, ! dx1(1:10,16, 6)
     .-005,-001,-006, 001, 014, 004, 012, 012, 021,-002, ! dx1(1:10,17, 6)
     .-011,-007,-004,-005, 005, 003, 011, 019, 014,-004, ! dx1(1:10,18, 6)
     . 008, 022, 014, 014, 030, 022, 020, 028, 025, 007, ! dx1(1:10,19, 6)
     . 000, 004,-003, 003, 021, 016, 017, 016, 009,-006, ! dx1(1:10,20, 6)
     . 002, 011, 004, 013, 025, 017, 015, 004, 014, 015, ! dx1(1:10,21, 6)
     . 000, 011,-002, 003, 015, 014, 013,-009, 002, 000, ! dx1(1:10,22, 6)
     . 009, 007,-001, 011, 011, 009, 014,-017, 004, 013, ! dx1(1:10,23, 6)
     . 005, 012, 000, 009, 011, 005, 010,-014, 009, 008, ! dx1(1:10,24, 6)
     . 002, 004,-007,-016,-011,-025,-006,-021,-005,-003, ! dx1(1:10,25, 6)
     .-007,-004,-021,-014,-020,-026,-016,-035,-018,-002, ! dx1(1:10,26, 6)
     .-005,-015,-016,-015,-018,-030,-010,-012,-014, 001, ! dx1(1:10,27, 6)
     .-005,-007,-020,-018,-005,-032,-012,-001,-019,-003, ! dx1(1:10,28, 6)
     . 005,-006,-028,-020,-013,-020,-009, 001,-014, 001, ! dx1(1:10,29, 6)
     . 007,-003,-023,-021,-001, 016,-006, 016,-010, 004, ! dx1(1:10,30, 6)
     .-005,-014,-029,-028,-012, 041,-008, 008,-013,-007, ! dx1(1:10,31, 6)
     . 013,-022,-021,-028,-012, 043,-010, 017,-007, 000, ! dx1(1:10,32, 6)
     .-002, 003,-015,-024, 000, 057,-006, 027, 005, 004, ! dx1(1:10,33, 6)
     . 000,-019,-011,-035,-015, 036,-017,-002, 001,-015, ! dx1(1:10,34, 6)
     . 008,-006, 000,-019,-005, 022,-016,-009, 009, 001, ! dx1(1:10,35, 6)
     . 018, 001, 022, 013, 006, 013,-001,-003, 017,-007, ! dx1(1:10,36, 6)
     . 013, 003, 036, 011, 000,-075,-010,-018,-009,-007, ! dx1(1:10,37, 6)
     . 033, 035, 071, 051, 028,-066, 024, 002, 004, 022, ! dx1(1:10,38, 6)
     .-011, 021,-028, 003, 017, 009, 006, 029,-002,-006, ! dx1(1:10, 1, 7)
     . 016, 017,-021, 009, 030, 018, 009, 036, 005, 008, ! dx1(1:10, 2, 7)
     . 020, 026, 029, 024, 032, 027, 016, 017, 018, 020, ! dx1(1:10, 3, 7)
     . 032, 033, 030, 032, 024, 021, 020, 001, 021, 021, ! dx1(1:10, 4, 7)
     . 007, 015, 030, 018, 010, 009, 012,-014, 012, 003, ! dx1(1:10, 5, 7)
     . 026, 043, 042, 030, 026, 017, 026, 000, 023, 010, ! dx1(1:10, 6, 7)
     . 028, 027, 039, 036, 022, 010, 030, 007, 022, 008, ! dx1(1:10, 7, 7)
     .-008,-005,-003,-008,-020,-023,-013,-029,-019,-026, ! dx1(1:10, 8, 7)
     .-005, 004, 002, 000,-009,-012, 004,-012,-007,-021, ! dx1(1:10, 9, 7)
     .-015, 003, 002,-003,-013,-015,-001,-006,-009,-013, ! dx1(1:10,10, 7)
     . 004,-003, 001, 000,-017,-006, 000,-018,-011,-018, ! dx1(1:10,11, 7)
     . 005, 009, 014,-002,-006, 001, 007,-018, 000,-009, ! dx1(1:10,12, 7)
     .-021,-019,-016,-015,-029,-012,-012,-036,-026,-028, ! dx1(1:10,13, 7)
     . 025, 018, 008, 012,-002, 007, 011,-017, 002,-004, ! dx1(1:10,14, 7)
     . 000,-002, 000, 002, 001, 005, 008,-011, 007,-001, ! dx1(1:10,15, 7)
     .-007, 003,-011,-005,-003,-002, 003, 000, 012,-012, ! dx1(1:10,16, 7)
     .-002,-007,-003, 000, 011, 007, 010, 020, 019,-011, ! dx1(1:10,17, 7)
     .-002, 007, 001, 008, 024, 014, 021, 034, 028, 003, ! dx1(1:10,18, 7)
     . 003,-003, 002, 003, 017, 013, 009, 035, 015, 006, ! dx1(1:10,19, 7)
     . 000,-009,-004, 006, 011, 008, 008, 023, 012, 008, ! dx1(1:10,20, 7)
     . 011, 006, 009, 013, 023, 015, 015, 009, 010, 010, ! dx1(1:10,21, 7)
     .-005, 001, 000, 008, 014, 009, 012,-013, 009, 011, ! dx1(1:10,22, 7)
     . 008, 008, 003, 011, 009, 008, 010,-017, 007, 002, ! dx1(1:10,23, 7)
     . 006, 006, 008, 010, 008, 002, 017,-016, 014, 009, ! dx1(1:10,24, 7)
     . 006,-001,-003, 007,-001,-006, 002,-014, 012, 014, ! dx1(1:10,25, 7)
     .-013,-010,-019,-015,-009,-021,-005,-023, 000, 000, ! dx1(1:10,26, 7)
     .-010,-009,-012,-008,-012,-026,-004,-005,-005, 000, ! dx1(1:10,27, 7)
     . 001,-011,-014,-012,-007,-023,-006, 007,-007, 008, ! dx1(1:10,28, 7)
     .-004,-010,-018,-014,-007,-014,-005, 011,-012, 003, ! dx1(1:10,29, 7)
     .-009,-010,-012,-017, 000, 017,-004, 029,-004, 006, ! dx1(1:10,30, 7)
     . 000,-015,-030,-014,-005, 041,-006, 027,-006, 001, ! dx1(1:10,31, 7)
     .-016,-027,-024,-031,-004, 044,-012, 033,-004,-004, ! dx1(1:10,32, 7)
     .-001,-011,-029,-023,-005, 043,-013, 025,-007,-002, ! dx1(1:10,33, 7)
     .-004,-014,-012,-024,-006, 035,-006, 012, 001, 000, ! dx1(1:10,34, 7)
     . 026, 000,-005,-006,-003, 031,-008,-005, 002,-006, ! dx1(1:10,35, 7)
     . 006, 009, 021, 003, 004, 011, 000, 000, 007, 002, ! dx1(1:10,36, 7)
     .-007,-015, 030, 016, 000,-080,-004,-020,-010,-013, ! dx1(1:10,37, 7)
     . 042, 033, 065, 040, 022,-066, 016, 003, 012, 008, ! dx1(1:10,38, 7)
     . 005, 002,-008, 012, 043, 012,-010, 048, 008,-004, ! dx1(1:10, 1, 8)
     . 016, 010, 011, 031, 032, 016, 017, 040, 018, 032, ! dx1(1:10, 2, 8)
     . 038, 057, 049, 055, 055, 046, 037, 045, 041, 050, ! dx1(1:10, 3, 8)
     . 009, 026, 048, 032, 028, 025, 023, 011, 025, 027, ! dx1(1:10, 4, 8)
     . 035, 039, 048, 039, 030, 024, 019, 007, 025, 016, ! dx1(1:10, 5, 8)
     . 044, 043, 058, 047, 033, 022, 034, 015, 029, 032, ! dx1(1:10, 6, 8)
     . 018, 033, 044, 037, 028, 031, 041, 028, 030, 030, ! dx1(1:10, 7, 8)
     .-012,-002, 004, 000,-005,-012,-005,-009,-009,-019, ! dx1(1:10, 8, 8)
     .-003, 000, 002, 004,-007,-004, 005, 001, 002,-006, ! dx1(1:10, 9, 8)
     . 001,-006,-004,-003,-014, 001,-002,-011,-003,-016, ! dx1(1:10,10, 8)
     . 012, 008, 004,-001,-011,-002, 006,-024, 002,-005, ! dx1(1:10,11, 8)
     .-001, 006, 007, 008,-011, 009, 010,-025, 000, 004, ! dx1(1:10,12, 8)
     . 028, 004, 000,-008,-019,-005,-003,-030,-005,-014, ! dx1(1:10,13, 8)
     . 004, 004, 004, 008,-006, 014, 013,-012, 005, 006, ! dx1(1:10,14, 8)
     . 014, 016, 001, 003, 001, 005, 010,-005, 011, 003, ! dx1(1:10,15, 8)
     . 004, 004,-007,-005, 000,-002, 002, 003, 008,-008, ! dx1(1:10,16, 8)
     . 017, 014,-003, 009, 010,-002, 011, 025, 015, 011, ! dx1(1:10,17, 8)
     . 010, 003, 012, 012, 022, 012, 019, 052, 026, 006, ! dx1(1:10,18, 8)
     .-021,-003, 000, 012, 018, 010, 015, 030, 007, 006, ! dx1(1:10,19, 8)
     .-010, 001,-004, 004, 011,-002, 007, 023, 004,-006, ! dx1(1:10,20, 8)
     .-002,-003, 002, 003, 009, 004, 009, 004, 008, 003, ! dx1(1:10,21, 8)
     .-003,-007, 000, 005, 011, 007, 008, 000, 003, 010, ! dx1(1:10,22, 8)
     . 014, 002, 005, 021, 014, 009, 016,-008, 014, 016, ! dx1(1:10,23, 8)
     . 006,-008, 001, 016, 005, 008, 010,-006, 017, 015, ! dx1(1:10,24, 8)
     . 004, 000, 002, 011, 002,-021, 016,-002, 015, 011, ! dx1(1:10,25, 8)
     . 013,-003, 003, 020,-004,-011, 008,-002, 017, 016, ! dx1(1:10,26, 8)
     .-017,-035,-016,-015,-011,-028,-005, 002,-001, 000, ! dx1(1:10,27, 8)
     .-020,-022,-015,-013,-008,-030,-007, 002,-008, 000, ! dx1(1:10,28, 8)
     . 001,-012,-010,-017,-007,-015,-007, 021, 002, 007, ! dx1(1:10,29, 8)
     . 007, 004,-004,-009, 013, 017, 002, 044, 005, 015, ! dx1(1:10,30, 8)
     .-011,-012,-021,-027, 001, 044,-015, 038,-004, 003, ! dx1(1:10,31, 8)
     .-002,-016,-028,-019,-004, 039,-023, 030, 000,-003, ! dx1(1:10,32, 8)
     . 007,-002,-011,-017, 005, 066, 002, 040,-001, 014, ! dx1(1:10,33, 8)
     . 002,-006,-017,-012,-004, 037,-004, 018,-011, 001, ! dx1(1:10,34, 8)
     . 011,-008, 003,-013,-003, 031, 000, 011,-012, 003, ! dx1(1:10,35, 8)
     . 026, 008, 015, 016, 007, 024, 009, 005,-006, 009, ! dx1(1:10,36, 8)
     . 015, 020, 038, 011, 005,-075, 004,-012,-011, 007, ! dx1(1:10,37, 8)
     . 040, 017, 051, 040, 012,-057, 020, 002,-001, 026, ! dx1(1:10,38, 8)
     .-031,-047,-009, 008, 018,-009,-005, 030, 022, 001, ! dx1(1:10, 1, 9)
     .-006,-010,-045,-015,-010,-028,-027, 004,-029,-033, ! dx1(1:10, 2, 9)
     . 010, 000, 016, 013, 009, 002,-004,-008, 000, 003, ! dx1(1:10, 3, 9)
     . 000, 000, 004, 012, 002,-005,-011,-017,-010,-010, ! dx1(1:10, 4, 9)
     .-022,-003, 006,-001,-015,-012,-017,-041,-016,-014, ! dx1(1:10, 5, 9)
     . 011, 007, 027, 010,-002,-008, 000,-019,-002,-006, ! dx1(1:10, 6, 9)
     .-018, 009, 011, 014, 010,-003, 009,-004, 002,-006, ! dx1(1:10, 7, 9)
     .-018, 000,-004,-010,-014,-024,-009,-018,-019,-033, ! dx1(1:10, 8, 9)
     . 000,-006,-001,-003,-013,-012, 000,-006,-004,-009, ! dx1(1:10, 9, 9)
     .-006,-003,-005,-013,-018,-005,-001,-009,-003,-018, ! dx1(1:10,10, 9)
     . 006, 008, 000,-001,-015,-008, 004,-028,-007,-008, ! dx1(1:10,11, 9)
     . 022, 010, 012, 006,-008, 012, 012,-025, 011,-001, ! dx1(1:10,12, 9)
     .-005, 006,-007,-015,-017,-015,-003,-036,-016,-014, ! dx1(1:10,13, 9)
     .-011,-002,-015,-011,-025,-013,-014,-035,-018,-014, ! dx1(1:10,14, 9)
     .-012,-018,-020,-021,-026,-013,-015,-028,-017,-027, ! dx1(1:10,15, 9)
     .-013,-017,-026,-019,-017,-018,-015,-012,-010,-022, ! dx1(1:10,16, 9)
     .-001,-006,-013,-008, 000,-013,-004, 012,-002,-005, ! dx1(1:10,17, 9)
     .-008,-008,-006, 001, 008,-002, 004, 031, 004,-001, ! dx1(1:10,18, 9)
     .-023,-023,-018,-012, 000,-004, 000, 021, 000, 000, ! dx1(1:10,19, 9)
     .-012,-017,-023,-011,-006,-010,-006, 004,-013,-013, ! dx1(1:10,20, 9)
     .-022,-013,-004,-009, 000,-010,-001,-002,-002,-005, ! dx1(1:10,21, 9)
     .-012,-009,-006, 010, 005,-003, 006,-009, 003,-001, ! dx1(1:10,22, 9)
     .-007,-011, 000, 009, 004, 000, 003,-013, 005, 011, ! dx1(1:10,23, 9)
     .-011,-013,-004, 011, 003, 001, 014,-013, 008, 005, ! dx1(1:10,24, 9)
     . 008,-021,-003, 002,-009,-018, 003,-002, 009, 003, ! dx1(1:10,25, 9)
     .-060,-045,-050,-032,-050,-059,-037,-039,-024,-042, ! dx1(1:10,26, 9)
     .-021,-051,-035,-030,-035,-047,-029,-026,-025,-015, ! dx1(1:10,27, 9)
     .-052,-048,-046,-039,-033,-058,-029,-016,-030,-039, ! dx1(1:10,28, 9)
     .-017,-041,-033,-033,-028,-034,-030, 001,-030,-021, ! dx1(1:10,29, 9)
     .-021,-026,-026,-033,-015,-008,-028, 019,-017,-008, ! dx1(1:10,30, 9)
     .-027,-031,-038,-047,-018, 028,-029, 019,-022,-029, ! dx1(1:10,31, 9)
     . 000,-030,-038,-024,-005, 040,-013, 031,-007,-002, ! dx1(1:10,32, 9)
     . 010, 008,-004,-015, 015, 063, 002, 041, 002, 015, ! dx1(1:10,33, 9)
     . 006,-006,-003, 002, 013, 048, 001, 031, 000, 014, ! dx1(1:10,34, 9)
     . 018, 001, 007, 009, 003, 041, 006, 016,-007, 003, ! dx1(1:10,35, 9)
     . 032, 017, 034, 022, 021, 038, 020, 010, 007, 021, ! dx1(1:10,36, 9)
     . 004, 023, 061, 030, 024,-063, 020, 003,-001, 010, ! dx1(1:10,37, 9)
     . 027,-028, 009,-001,-020,-122,-022,-032,-041,-011, ! dx1(1:10,38, 9)
     . 004,-009,-014, 012, 009, 003,-002, 029,-009, 010, ! dx1(1:10, 1,10)
     .-022,-014,-012, 009, 010, 000,-006, 025,-005, 000, ! dx1(1:10, 2,10)
     .-021, 000, 010, 018, 013, 005,-006, 002,-004,-001, ! dx1(1:10, 3,10)
     . 010, 018, 018, 024, 011, 008, 000, 006, 005, 005, ! dx1(1:10, 4,10)
     .-013,-013, 012, 014, 000,-005,-008,-005,-007,-001, ! dx1(1:10, 5,10)
     . 000, 017, 024, 019, 018, 012, 013, 004, 012, 011, ! dx1(1:10, 6,10)
     . 027, 025, 034, 021, 017, 012, 026, 024, 027, 014, ! dx1(1:10, 7,10)
     .-001, 003, 006, 001, 001, 003, 002, 003, 009, 000, ! dx1(1:10, 8,10)
     . 009, 016, 004, 001,-001, 007, 010,-004, 011, 005, ! dx1(1:10, 9,10)
     . 015, 008, 001, 007,-005, 009, 007,-014, 008, 003, ! dx1(1:10,10,10)
     . 020, 022, 005, 006,-005, 011, 005,-021, 006, 002, ! dx1(1:10,11,10)
     . 023, 026, 013, 007,-001, 007, 011,-030, 010, 002, ! dx1(1:10,12,10)
     . 017, 006,-003, 002,-012, 003, 002,-033, 009,-006, ! dx1(1:10,13,10)
     . 008, 002,-013,-010,-021,-012,-005,-039,-007,-009, ! dx1(1:10,14,10)
     . 000,-004,-010,-007,-012,-006,-006,-021,-003,-003, ! dx1(1:10,15,10)
     . 000,-010,-014,-009,-010,-016,-008,-013,-010,-011, ! dx1(1:10,16,10)
     .-006,-008,-011, 001, 000,-006, 000, 012,-004,-002, ! dx1(1:10,17,10)
     . 002,-004,-006, 000, 002,-010,-002, 022,-007,-004, ! dx1(1:10,18,10)
     .-003,-006,-002, 008, 014, 000, 009, 019,-002, 002, ! dx1(1:10,19,10)
     .-020,-018,-011,-010,-003,-021,-013, 005,-007,-009, ! dx1(1:10,20,10)
     .-003,-005,-001, 007, 007, 006, 006, 003,-004, 003, ! dx1(1:10,21,10)
     .-007,-004, 009, 004, 004,-003, 008, 000, 002, 002, ! dx1(1:10,22,10)
     .-010,-015, 000, 006, 000,-004, 010,-002, 001, 004, ! dx1(1:10,23,10)
     . 001,-003, 015, 015, 000,-007, 014, 008, 013, 007, ! dx1(1:10,24,10)
     .-004,-008, 009, 011, 013,-010, 022, 022, 025, 018, ! dx1(1:10,25,10)
     .-044,-058,-040,-028,-041,-058,-032,-028,-028,-036, ! dx1(1:10,26,10)
     .-036,-051,-033,-025,-035,-053,-026,-006,-020,-027, ! dx1(1:10,27,10)
     .-025,-051,-033,-028,-029,-049,-026, 003,-021,-018, ! dx1(1:10,28,10)
     .-039,-043,-029,-024,-023,-030,-024, 008,-006,-017, ! dx1(1:10,29,10)
     .-022,-021,-018,-018,-010, 002,-019, 021, 004,-011, ! dx1(1:10,30,10)
     .-024,-037,-030,-032,-013, 014,-028, 008,-014,-013, ! dx1(1:10,31,10)
     . 003, 004,-002,-014, 017, 057,-005, 034, 018, 009, ! dx1(1:10,32,10)
     . 023, 019,-003, 000, 020, 056, 008, 044, 009, 008, ! dx1(1:10,33,10)
     . 011, 006, 001,-010, 012, 043, 000, 027,-007, 005, ! dx1(1:10,34,10)
     . 008, 024, 013, 006, 017, 034, 011, 022,-001, 013, ! dx1(1:10,35,10)
     . 036, 021, 032, 019, 023, 029, 023, 020, 000, 021, ! dx1(1:10,36,10)
     . 021, 024, 053, 031, 013,-057, 016, 007,-013, 008, ! dx1(1:10,37,10)
     .-029,-010, 014,-011,-018,-102,-017,-036,-056,-033, ! dx1(1:10,38,10)
     .-052,-004, 002, 021, 019, 004,-015, 035,-003, 002, ! dx1(1:10, 1,11)
     .-031,-022,-012, 000, 000,-003,-014, 020,-017,-011, ! dx1(1:10, 2,11)
     .-004, 000, 008, 012, 011, 013,-007, 009, 003, 006, ! dx1(1:10, 3,11)
     .-013, 002, 017, 019, 018, 014,-001, 012, 006, 009, ! dx1(1:10, 4,11)
     .-003, 010, 007, 002, 004, 000,-004, 010, 001,-001, ! dx1(1:10, 5,11)
     . 021, 012, 029, 025, 018, 016, 017, 025, 013, 023, ! dx1(1:10, 6,11)
     . 018, 031, 027, 021, 023, 025, 031, 036, 033, 024, ! dx1(1:10, 7,11)
     . 012, 003, 012, 006, 003, 007, 011, 015, 011, 005, ! dx1(1:10, 8,11)
     . 008, 009,-001, 000, 000, 010, 009,-004, 014, 000, ! dx1(1:10, 9,11)
     . 019, 017, 000, 000, 001, 014, 010,-015, 017, 008, ! dx1(1:10,10,11)
     . 022, 027, 000, 001,-006, 007, 003,-027, 009, 000, ! dx1(1:10,11,11)
     . 031, 022, 007, 010,-001, 015, 010,-028, 019, 012, ! dx1(1:10,12,11)
     . 006, 027, 006, 008, 002, 005, 010,-023, 011, 006, ! dx1(1:10,13,11)
     . 030, 008,-006,-001,-007,-003,-006,-035, 003, 000, ! dx1(1:10,14,11)
     . 007,-006,-006,-001,-007,-004,-003,-026,-003,-005, ! dx1(1:10,15,11)
     . 012,-005,-002,-007, 001,-007,-004,-005, 000, 007, ! dx1(1:10,16,11)
     . 003,-002,-007,-001, 000,-005, 000, 008,-011, 004, ! dx1(1:10,17,11)
     . 001,-002, 000, 003, 013,-004, 007, 022, 001, 006, ! dx1(1:10,18,11)
     .-008,-008,-002, 002, 007,-006, 000, 015,-008, 009, ! dx1(1:10,19,11)
     .-025,-010,-008,-009, 001,-011, 003, 005,-013,-011, ! dx1(1:10,20,11)
     . 001, 000, 004, 010, 007,-004, 009, 006, 000, 004, ! dx1(1:10,21,11)
     .-002,-006, 007, 010, 009, 001, 018, 006, 002, 007, ! dx1(1:10,22,11)
     .-009, 001, 015, 012, 004,-003, 018, 011, 010, 012, ! dx1(1:10,23,11)
     .-005, 005, 015, 015, 011,-008, 028, 016, 018, 010, ! dx1(1:10,24,11)
     . 006,-002, 026, 019, 005,-001, 021, 040, 015, 010, ! dx1(1:10,25,11)
     .-048,-039,-027,-027,-036,-052,-022,-004,-023,-041, ! dx1(1:10,26,11)
     .-041,-054,-027,-020,-032,-055,-023, 001,-022,-026, ! dx1(1:10,27,11)
     .-043,-048,-031,-023,-026,-055,-024, 002,-014,-030, ! dx1(1:10,28,11)
     .-018,-032,-017,-018,-015,-026,-021, 014,-001,-014, ! dx1(1:10,29,11)
     .-014,-031,-011,-015,-005, 001,-012, 017, 008,-010, ! dx1(1:10,30,11)
     .-019,-016,-024,-027,-006, 030,-026, 006, 000,-005, ! dx1(1:10,31,11)
     .-019,-003,-001,-008, 015, 049,-005, 030, 023, 003, ! dx1(1:10,32,11)
     . 020, 019, 005, 003, 030, 063, 016, 038, 016, 023, ! dx1(1:10,33,11)
     . 022, 018, 005, 005, 024, 048, 011, 034, 006, 014, ! dx1(1:10,34,11)
     . 021, 007, 012, 011, 020, 028, 011, 029, 002, 017, ! dx1(1:10,35,11)
     . 028, 026, 037, 022, 023, 026, 022, 026, 000, 015, ! dx1(1:10,36,11)
     . 037, 038, 047, 023, 011,-059, 013, 000,-021, 008, ! dx1(1:10,37,11)
     .-025,-027, 020,-011,-026,-115,-024,-032,-055,-027, ! dx1(1:10,38,11)
     .-014,-029,-011, 005, 000,-007,-008, 038,-007, 000, ! dx1(1:10, 1,12)
     .-005,-030,-008, 003, 002, 000,-016, 026,-011,-003, ! dx1(1:10, 2,12)
     .-013, 001, 008, 020, 010, 012, 000, 016, 002, 005, ! dx1(1:10, 3,12)
     . 005, 012, 022, 027, 015, 022, 006, 037, 014, 016, ! dx1(1:10, 4,12)
     . 002, 001, 009, 016, 008, 010, 000, 028,-009, 014, ! dx1(1:10, 5,12)
     . 033, 029, 028, 028, 023, 021, 019, 034, 025, 015, ! dx1(1:10, 6,12)
     . 013, 010, 009, 009, 022, 018, 019, 028, 024, 024, ! dx1(1:10, 7,12)
     . 016, 011, 001, 007, 010, 007, 013, 024, 007, 014, ! dx1(1:10, 8,12)
     . 010, 019, 008, 002, 000, 012, 009, 002, 015, 019, ! dx1(1:10, 9,12)
     . 011, 020, 000, 002,-001, 012, 009,-017, 010, 010, ! dx1(1:10,10,12)
     . 018, 015,-012,-008,-011, 003,-004,-040, 001, 000, ! dx1(1:10,11,12)
     . 009, 012,-005, 002,-006, 005, 005,-030, 019, 002, ! dx1(1:10,12,12)
     . 026, 012, 000, 008,-003, 001,-002,-033, 013, 010, ! dx1(1:10,13,12)
     . 004, 003,-003,-001,-008,-009,-008,-038, 003,-005, ! dx1(1:10,14,12)
     . 011, 000,-010, 000,-003,-004,-007,-027, 003, 008, ! dx1(1:10,15,12)
     . 006, 003,-001, 000, 003,-001,-003,-008,-008,-001, ! dx1(1:10,16,12)
     . 005,-010,-007, 000, 005,-005,-003, 004,-010, 011, ! dx1(1:10,17,12)
     . 005,-002,-002, 004, 007,-002, 001, 013,-004, 008, ! dx1(1:10,18,12)
     . 000, 004, 009, 011, 017,-001, 007, 012,-004, 008, ! dx1(1:10,19,12)
     . 001,-002,-005, 005, 003,-010, 000, 005,-018, 004, ! dx1(1:10,20,12)
     . 003,-003, 008, 013, 012, 003, 017, 010, 000, 008, ! dx1(1:10,21,12)
     . 000, 010, 012, 013, 005,-005, 015, 011, 000, 001, ! dx1(1:10,22,12)
     . 004, 002, 025, 021, 015, 003, 028, 018, 005, 014, ! dx1(1:10,23,12)
     . 010, 009, 025, 024, 011, 002, 030, 033, 013, 017, ! dx1(1:10,24,12)
     . 020, 000, 035, 023, 011,-005, 026, 048, 018, 012, ! dx1(1:10,25,12)
     .-032,-045,-014,-016,-036,-046,-019, 005,-022,-030, ! dx1(1:10,26,12)
     .-029,-047,-014,-010,-026,-044,-014, 016,-010,-019, ! dx1(1:10,27,12)
     .-025,-043,-012,-009,-017,-038,-018, 013,-009,-018, ! dx1(1:10,28,12)
     .-019,-033,-021,-006,-017,-022,-014, 010,-001,-012, ! dx1(1:10,29,12)
     .-011,-022,-006,-004, 000, 003,-013, 009, 016,-010, ! dx1(1:10,30,12)
     .-021,-026,-023,-024,-008, 024,-025,-002, 000,-013, ! dx1(1:10,31,12)
     . 012,-005, 000, 000, 015, 055,-002, 019, 024, 018, ! dx1(1:10,32,12)
     . 012, 019, 008, 011, 032, 053, 010, 035, 020, 017, ! dx1(1:10,33,12)
     . 041, 025, 010, 006, 029, 044, 013, 027, 013, 020, ! dx1(1:10,34,12)
     . 013, 028, 015, 009, 021, 027, 010, 027, 000, 015, ! dx1(1:10,35,12)
     . 035, 031, 030, 015, 018, 017, 015, 027,-002, 022, ! dx1(1:10,36,12)
     . 009, 022, 031, 012, 001,-068, 004, 002,-033,-003, ! dx1(1:10,37,12)
     .-012,-005,-009,-018,-042,-117,-033,-039,-066,-024, ! dx1(1:10,38,12)
     .-053,-031, 000, 003, 000, 003,-012, 036,-018,-007, ! dx1(1:10, 1,13)
     .-068,-028,-002, 011,-001,-002,-014, 030,-004,-017, ! dx1(1:10, 2,13)
     .-005,-019, 008, 017, 008, 008,-007, 027, 001, 018, ! dx1(1:10, 3,13)
     . 011, 024, 023, 034, 025, 022, 013, 045, 010, 022, ! dx1(1:10, 4,13)
     . 000, 001, 015, 015, 013, 016, 001, 036, 004, 007, ! dx1(1:10, 5,13)
     . 016, 030, 019, 022, 030, 032, 021, 049, 021, 032, ! dx1(1:10, 6,13)
     . 023, 031, 013, 017, 018, 020, 022, 035, 019, 016, ! dx1(1:10, 7,13)
     .-005,-001,-009,-004,-005, 000,-003, 007, 000, 001, ! dx1(1:10, 8,13)
     . 027, 024,-002, 001, 005, 015, 015, 007, 013, 011, ! dx1(1:10, 9,13)
     . 022, 026, 000, 003,-001, 017, 009,-013, 015, 012, ! dx1(1:10,10,13)
     . 012, 021,-014, 002,-009, 006,-004,-030, 007, 011, ! dx1(1:10,11,13)
     . 023, 006,-008,-002,-008, 010,-002,-037, 014, 013, ! dx1(1:10,12,13)
     . 002, 014,-002, 005,-003, 003,-002,-032, 008, 008, ! dx1(1:10,13,13)
     . 030, 001,-010, 001,-010, 004,-012,-038, 008, 010, ! dx1(1:10,14,13)
     . 020, 009,-001, 004, 005,-001,-003,-021, 009, 009, ! dx1(1:10,15,13)
     . 010, 000,-009,-002,-005,-012,-013,-014,-005, 011, ! dx1(1:10,16,13)
     . 016, 002,-005, 005, 006,-010,-007, 000,-009, 007, ! dx1(1:10,17,13)
     .-001,-006, 000, 008, 012,-007, 003, 010,-005, 016, ! dx1(1:10,18,13)
     .-012, 004, 004, 008, 015,-002, 004, 013,-006, 012, ! dx1(1:10,19,13)
     .-001,-008,-005, 003, 004,-013,-004, 001,-020, 005, ! dx1(1:10,20,13)
     . 007, 017, 023, 023, 018, 002, 022, 023, 007, 017, ! dx1(1:10,21,13)
     .-002, 005, 018, 010, 011,-001, 027, 012,-001, 009, ! dx1(1:10,22,13)
     . 004, 007, 026, 017, 009, 001, 026, 027, 011, 008, ! dx1(1:10,23,13)
     . 009, 011, 038, 022, 017, 006, 040, 041, 018, 011, ! dx1(1:10,24,13)
     . 000, 004, 024, 020, 004,-002, 025, 053, 017, 000, ! dx1(1:10,25,13)
     .-024,-029,-009,-019,-035,-047,-018, 012,-024,-030, ! dx1(1:10,26,13)
     .-035,-029,-004,-002,-014,-038,-006, 018,-012,-019, ! dx1(1:10,27,13)
     .-023,-034,-015,-009,-022,-040,-014, 010,-008,-019, ! dx1(1:10,28,13)
     .-008,-028,-006, 003,-013,-025,-017, 003, 001,-011, ! dx1(1:10,29,13)
     .-006,-015,-009,-001, 001, 003,-009, 008, 014,-005, ! dx1(1:10,30,13)
     .-018,-017,-016,-024,-006, 022,-024,-017, 004,-005, ! dx1(1:10,31,13)
     . 022, 015, 011, 009, 028, 053, 002, 019, 039, 014, ! dx1(1:10,32,13)
     . 018, 024, 015, 009, 030, 057, 011, 027, 030, 025, ! dx1(1:10,33,13)
     . 027, 021, 011, 006, 024, 042, 011, 024, 016, 012, ! dx1(1:10,34,13)
     . 021, 024, 005, 007, 024, 025, 009, 018, 001, 008, ! dx1(1:10,35,13)
     . 033, 033, 021, 011, 016, 011, 018, 023,-002, 018, ! dx1(1:10,36,13)
     . 019, 029, 033, 018,-003,-075, 000, 007,-034,-008, ! dx1(1:10,37,13)
     .-038,-020,-012,-029,-035,-127,-037,-033,-076,-032, ! dx1(1:10,38,13)
     .-037,-044,-009, 003,-012,-013,-022, 026,-010,-011, ! dx1(1:10, 1,14)
     .-049,-043,-005, 013,-003,-003,-013, 030,-008,-013, ! dx1(1:10, 2,14)
     .-021,-002, 022, 024, 009, 015, 005, 038, 012, 007, ! dx1(1:10, 3,14)
     . 012, 006, 021, 040, 024, 028, 015, 050, 015, 020, ! dx1(1:10, 4,14)
     . 011, 025, 021, 029, 027, 022, 017, 051, 008, 024, ! dx1(1:10, 5,14)
     . 015, 038, 029, 036, 033, 034, 028, 063, 026, 026, ! dx1(1:10, 6,14)
     . 013, 023, 010, 016, 020, 025, 024, 038, 017, 028, ! dx1(1:10, 7,14)
     .-001, 018, 001, 002, 003, 011, 011, 024, 004, 010, ! dx1(1:10, 8,14)
     . 008, 017,-011,-003,-003, 008, 007,-003, 004, 016, ! dx1(1:10, 9,14)
     . 010, 013,-009,-007,-004, 018, 002,-010, 010, 015, ! dx1(1:10,10,14)
     . 027, 016,-006, 003,-003, 017, 002,-024, 013, 022, ! dx1(1:10,11,14)
     . 012, 013,-013,-006,-015, 007,-011,-039, 008, 005, ! dx1(1:10,12,14)
     . 021,-003,-010, 007,-001, 004,-007,-035, 018, 020, ! dx1(1:10,13,14)
     .-005, 005,-010, 004,-006, 005,-007,-037, 012, 010, ! dx1(1:10,14,14)
     . 021, 001,-011, 003,-001,-007,-015,-022, 009, 007, ! dx1(1:10,15,14)
     . 012, 002,-010,-001, 009,-002,-007,-002,-004, 007, ! dx1(1:10,16,14)
     .-001, 001,-008, 003, 005,-013,-016,-004,-012, 011, ! dx1(1:10,17,14)
     . 015, 001, 001, 003, 012,-008, 001, 005,-012, 009, ! dx1(1:10,18,14)
     . 002, 007, 006, 011, 007,-002, 006, 014,-002, 014, ! dx1(1:10,19,14)
     . 009, 007, 009, 011, 016, 002, 010, 013,-003, 015, ! dx1(1:10,20,14)
     . 005, 008, 018, 014, 013,-001, 022, 023, 000, 007, ! dx1(1:10,21,14)
     . 015, 007, 025, 018, 013, 003, 023, 024, 014, 016, ! dx1(1:10,22,14)
     . 016, 025, 047, 029, 021, 014, 047, 035, 018, 017, ! dx1(1:10,23,14)
     . 011, 020, 045, 024, 015, 010, 040, 045, 022, 018, ! dx1(1:10,24,14)
     .-013,-011, 013, 004,-018,-021, 005, 035,-008,-004, ! dx1(1:10,25,14)
     .-054,-042,-015,-021,-044,-049,-017, 003,-040,-042, ! dx1(1:10,26,14)
     .-021,-026,-001,-004,-028,-040,-013, 019,-016,-018, ! dx1(1:10,27,14)
     .-030,-025,-010,-008,-024,-041,-013, 005,-013,-023, ! dx1(1:10,28,14)
     .-024,-039,-018,-009,-020,-029,-020,-008,-008,-022, ! dx1(1:10,29,14)
     .-012,-015,-002,-001,-005, 004,-010,-004, 011,-006, ! dx1(1:10,30,14)
     .-018,-022,-013,-017, 007, 022,-013,-018, 006,-005, ! dx1(1:10,31,14)
     .-001,-003, 009, 009, 022, 046, 004, 000, 029, 020, ! dx1(1:10,32,14)
     . 017, 029, 008, 002, 032, 051, 006, 009, 022, 015, ! dx1(1:10,33,14)
     . 011, 020, 001,-005, 014, 034, 002, 010, 015, 015, ! dx1(1:10,34,14)
     . 016, 032, 008, 002, 017, 021, 007, 008, 005, 014, ! dx1(1:10,35,14)
     . 029, 033, 014, 010, 006, 003, 004, 017,-004, 005, ! dx1(1:10,36,14)
     . 025, 021, 020, 000,-010,-073,-003,-016,-039,-012, ! dx1(1:10,37,14)
     .-003, 011,-007,-027,-065,-132,-048,-060,-079,-074, ! dx1(1:10,38,14)
     .-068,-054,-013,-006,-017,-003,-018, 040,-008,-002, ! dx1(1:10, 1,15)
     .-051,-048,-009, 000,-016,-019,-026, 025,-012,-025, ! dx1(1:10, 2,15)
     .-022,-024, 015, 026, 007, 008,-003, 037, 008, 014, ! dx1(1:10, 3,15)
     .-008, 004, 025, 026, 027, 021, 013, 049, 015, 014, ! dx1(1:10, 4,15)
     . 013, 013, 022, 032, 022, 029, 016, 054, 022, 027, ! dx1(1:10, 5,15)
     . 020, 026, 023, 025, 033, 032, 028, 057, 021, 032, ! dx1(1:10, 6,15)
     . 034, 027, 012, 022, 019, 034, 031, 055, 022, 033, ! dx1(1:10, 7,15)
     . 007, 009,-002, 007, 007, 016, 016, 028, 009, 030, ! dx1(1:10, 8,15)
     . 016, 029, 000, 000, 009, 017, 010, 014, 011, 018, ! dx1(1:10, 9,15)
     . 000, 017,-013, 000,-003, 022, 006,-007, 012, 020, ! dx1(1:10,10,15)
     .-002, 006,-017,-004,-012, 013,-006,-020, 006, 016, ! dx1(1:10,11,15)
     . 010,-001,-018,-007,-012, 011,-014,-034, 009, 007, ! dx1(1:10,12,15)
     . 005,-002,-022,-010,-016,-001,-024,-046, 001, 007, ! dx1(1:10,13,15)
     . 006,-007,-022,-011,-017,-004,-024,-047, 006, 007, ! dx1(1:10,14,15)
     . 010, 002,-010, 003, 003, 000,-008,-012, 006, 013, ! dx1(1:10,15,15)
     . 006, 001,-006, 008, 007,-009,-010, 000, 000, 016, ! dx1(1:10,16,15)
     . 003,-001,-008, 000, 004,-010,-013, 003,-014, 009, ! dx1(1:10,17,15)
     . 011, 013, 006, 021, 021, 007, 004, 016, 001, 026, ! dx1(1:10,18,15)
     . 013, 011, 008, 011, 021,-003, 012, 011,-008, 005, ! dx1(1:10,19,15)
     . 010, 010, 005, 011, 016,-001, 006, 010,-008, 014, ! dx1(1:10,20,15)
     . 013, 018, 024, 019, 019, 010, 029, 034, 013, 012, ! dx1(1:10,21,15)
     . 015, 021, 037, 018, 022, 005, 036, 035, 015, 013, ! dx1(1:10,22,15)
     . 040, 039, 056, 042, 038, 027, 056, 056, 038, 030, ! dx1(1:10,23,15)
     . 020, 029, 057, 034, 026, 018, 053, 051, 030, 017, ! dx1(1:10,24,15)
     .-028,-022, 002,-008,-027,-032, 000, 014,-020,-018, ! dx1(1:10,25,15)
     .-036,-034,-016,-033,-051,-050,-021,-003,-045,-039, ! dx1(1:10,26,15)
     .-023,-030,-002,-004,-022,-039,-002, 014,-023,-026, ! dx1(1:10,27,15)
     .-023,-043,-010,-013,-032,-045,-021,-003,-020,-024, ! dx1(1:10,28,15)
     .-030,-036,-021,-019,-025,-036,-025,-017,-015,-021, ! dx1(1:10,29,15)
     .-016,-018,-007,-010,-007,-009,-017,-018, 005,-010, ! dx1(1:10,30,15)
     .-019,-031,-021,-027,-006, 014,-027,-044,-002,-026, ! dx1(1:10,31,15)
     . 004, 007, 009, 001, 019, 041,-002,-015, 029, 011, ! dx1(1:10,32,15)
     . 017, 018, 012, 000, 021, 045, 007,-006, 024, 010, ! dx1(1:10,33,15)
     . 010, 016,-002, 001, 015, 033,-007,-006, 013, 003, ! dx1(1:10,34,15)
     . 021, 025,-002,-007, 008, 017, 001, 000, 003,-007, ! dx1(1:10,35,15)
     . 030, 024, 000, 004, 001,-003,-003,-009,-010, 000, ! dx1(1:10,36,15)
     . 002, 022,-007,-020,-043,-096,-028,-037,-062,-039, ! dx1(1:10,37,15)
     .-015,-020,-032,-049,-072,-136,-056,-067,-092,-084, ! dx1(1:10,38,15)
     .-066,-087,-019, 001,-033,-025,-031, 018,-024,-055, ! dx1(1:10, 1,16)
     .-094,-076,-017,-004,-034,-022,-034, 017,-028,-024, ! dx1(1:10, 2,16)
     .-045,-035, 006, 014,-001,-006,-005, 022,-001,-010, ! dx1(1:10, 3,16)
     .-023,-012, 013, 020, 008, 004, 005, 038, 008, 018, ! dx1(1:10, 4,16)
     . 000,-008, 010, 020, 014, 017, 015, 038, 012, 010, ! dx1(1:10, 5,16)
     . 007, 026, 019, 021, 030, 028, 020, 049, 024, 032, ! dx1(1:10, 6,16)
     . 019, 026, 011, 008, 023, 028, 028, 046, 019, 029, ! dx1(1:10, 7,16)
     . 013, 023,-004, 001, 013, 018, 018, 040, 012, 031, ! dx1(1:10, 8,16)
     . 002, 010,-004, 002, 004, 020, 019, 011, 015, 024, ! dx1(1:10, 9,16)
     . 034, 031, 000, 004, 006, 025, 008, 007, 022, 029, ! dx1(1:10,10,16)
     . 006, 010,-010,-001,-001, 028, 000,-013, 012, 020, ! dx1(1:10,11,16)
     . 001,-006,-026,-014,-012, 014,-018,-027, 000, 017, ! dx1(1:10,12,16)
     . 004, 004,-024,-011,-019,-005,-021,-040, 004, 005, ! dx1(1:10,13,16)
     . 012, 002,-031,-008,-013, 006,-017,-035, 005, 007, ! dx1(1:10,14,16)
     . 012,-002,-016,-005, 001,-014,-016,-016,-003, 010, ! dx1(1:10,15,16)
     . 018, 002,-007, 008, 009,-002,-004, 008, 000, 013, ! dx1(1:10,16,16)
     . 007,-003,-007, 008, 008,-010,-011, 007,-009, 019, ! dx1(1:10,17,16)
     . 025, 013, 012, 012, 024, 005, 008, 022, 000, 021, ! dx1(1:10,18,16)
     . 026, 017, 028, 029, 026, 012, 023, 036, 002, 023, ! dx1(1:10,19,16)
     . 018, 013, 014, 020, 022, 005, 014, 026, 000, 018, ! dx1(1:10,20,16)
     . 038, 039, 047, 028, 031, 011, 038, 041, 022, 027, ! dx1(1:10,21,16)
     . 029, 027, 044, 024, 023, 010, 038, 043, 024, 013, ! dx1(1:10,22,16)
     . 026, 042, 063, 034, 034, 022, 060, 050, 043, 029, ! dx1(1:10,23,16)
     . 022, 034, 054, 030, 024, 019, 058, 052, 040, 017, ! dx1(1:10,24,16)
     .-012,-012, 023, 004,-007,-012, 022, 023,-004,-006, ! dx1(1:10,25,16)
     .-024,-032,-001,-011,-022,-028, 006, 006,-022,-032, ! dx1(1:10,26,16)
     .-017,-021, 000,-009,-034,-041,-005, 006,-023,-028, ! dx1(1:10,27,16)
     .-038,-038,-020,-023,-033,-048,-016,-013,-032,-029, ! dx1(1:10,28,16)
     .-030,-048,-026,-025,-040,-046,-029,-032,-025,-036, ! dx1(1:10,29,16)
     .-031,-032,-022,-020,-026,-024,-029,-037,-014,-024, ! dx1(1:10,30,16)
     .-024,-022,-020,-030,-014, 005,-029,-061,-001,-030, ! dx1(1:10,31,16)
     . 015, 005, 010,-005, 021, 036,-001,-026, 025, 002, ! dx1(1:10,32,16)
     . 007, 016,-001,-004, 013, 032,-008,-035, 013,-003, ! dx1(1:10,33,16)
     . 013, 026, 000,-011, 008, 031,-001,-025, 010,-007, ! dx1(1:10,34,16)
     . 023, 026,-005,-002,-001, 018,-007,-022, 001,-002, ! dx1(1:10,35,16)
     .-002, 022,-011,-019,-013,-009,-018,-036,-022,-023, ! dx1(1:10,36,16)
     . 020, 007,-018,-020,-054,-110,-035,-078,-060,-048, ! dx1(1:10,37,16)
     .-070,-029,-054,-075,-085,-161,-067,-116,-103,-063, ! dx1(1:10,38,16)
     .-131,-102,-041,-039,-050,-061,-061, 002,-044,-058, ! dx1(1:10, 1,17)
     .-100,-102,-034,-024,-048,-055,-055,-009,-043,-061, ! dx1(1:10, 2,17)
     .-075,-053,-008,-007,-021,-020,-026, 015,-021,-021, ! dx1(1:10, 3,17)
     .-040,-036, 003, 001,-003,-007,-009, 011,-001,-014, ! dx1(1:10, 4,17)
     .-039,-025,-005,-006,-006,-005,-003, 015, 004, 001, ! dx1(1:10, 5,17)
     .-006,-001, 007, 004, 013, 019, 021, 038, 016, 019, ! dx1(1:10, 6,17)
     . 000, 020, 003, 004, 010, 015, 021, 039, 022, 025, ! dx1(1:10, 7,17)
     .-012, 016, 002,-002, 005, 015, 009, 029, 008, 023, ! dx1(1:10, 8,17)
     . 004, 013,-010,-008, 001, 014, 009, 010, 012, 022, ! dx1(1:10, 9,17)
     . 006, 006,-013,-009, 000, 023, 005, 002, 008, 019, ! dx1(1:10,10,17)
     . 020, 010,-002, 001, 001, 026, 002,-001, 017, 025, ! dx1(1:10,11,17)
     . 013, 004,-011, 000,-007, 020,-006,-013, 007, 023, ! dx1(1:10,12,17)
     . 015,-001,-019,-002,-002, 008,-016,-021, 003, 001, ! dx1(1:10,13,17)
     . 014,-004,-024,-001,-003, 005,-020,-022, 000, 017, ! dx1(1:10,14,17)
     . 020, 000,-010, 003, 003, 001,-010,-004, 000, 013, ! dx1(1:10,15,17)
     . 014, 002,-005, 010, 019, 002,-001, 007, 000, 016, ! dx1(1:10,16,17)
     . 031, 017, 001, 009, 022, 010, 002, 014,-002, 020, ! dx1(1:10,17,17)
     . 034, 018, 019, 026, 030, 013, 017, 035, 005, 024, ! dx1(1:10,18,17)
     . 032, 033, 033, 020, 032, 010, 025, 038, 005, 025, ! dx1(1:10,19,17)
     . 046, 040, 034, 018, 023, 008, 026, 035, 003, 021, ! dx1(1:10,20,17)
     . 024, 037, 044, 027, 031, 019, 036, 044, 019, 022, ! dx1(1:10,21,17)
     . 039, 044, 051, 025, 027, 018, 048, 038, 036, 019, ! dx1(1:10,22,17)
     . 031, 040, 055, 026, 027, 017, 050, 041, 030, 014, ! dx1(1:10,23,17)
     . 016, 024, 049, 023, 015, 012, 051, 042, 031, 010, ! dx1(1:10,24,17)
     . 026, 022, 041, 017, 010, 012, 052, 034, 018, 017, ! dx1(1:10,25,17)
     .-016,-013, 011,-007,-023,-021, 012, 006,-027,-026, ! dx1(1:10,26,17)
     .-020,-032,-012,-017,-034,-037,-003,-005,-030,-029, ! dx1(1:10,27,17)
     .-037,-038,-013,-022,-041,-046,-012,-018,-033,-032, ! dx1(1:10,28,17)
     .-035,-033,-028,-026,-037,-045,-024,-037,-030,-035, ! dx1(1:10,29,17)
     .-031,-042,-026,-038,-032,-037,-033,-053,-023,-033, ! dx1(1:10,30,17)
     .-019,-025,-030,-033,-018, 001,-032,-069,-013,-030, ! dx1(1:10,31,17)
     . 001,-002,-003,-002, 006, 027,-004,-044, 013,-007, ! dx1(1:10,32,17)
     . 005, 028,-005,-013, 007, 033,-010,-045, 011,-006, ! dx1(1:10,33,17)
     . 033, 030,-007,-011, 007, 029,-011,-053, 001,-012, ! dx1(1:10,34,17)
     . 001, 022,-011,-016,-006, 012,-017,-050,-005,-010, ! dx1(1:10,35,17)
     . 028, 019,-015,-017,-019,-009,-023,-065,-020,-019, ! dx1(1:10,36,17)
     .-006, 004,-025,-036,-053,-105,-045,-098,-071,-048, ! dx1(1:10,37,17)
     .-030,-018,-054,-064,-094,-156,-084,-132,-105,-086, ! dx1(1:10,38,17)
     .-161,-152,-066,-087,-086,-099,-097,-046,-088,-098, ! dx1(1:10, 1,18)
     .-168,-146,-074,-069,-090,-095,-093,-037,-080,-095, ! dx1(1:10, 2,18)
     .-119,-113,-042,-049,-056,-071,-059,-042,-052,-074, ! dx1(1:10, 3,18)
     .-094,-086,-043,-042,-052,-062,-047,-031,-033,-049, ! dx1(1:10, 4,18)
     .-083,-062,-024,-034,-030,-035,-031,-029,-012,-036, ! dx1(1:10, 5,18)
     .-066,-041,-025,-031,-025,-026,-019,-009, 000,-020, ! dx1(1:10, 6,18)
     .-034,-012,-016,-032,-003, 001, 001, 009, 006,-003, ! dx1(1:10, 7,18)
     .-029,-026,-027,-027,-016,-011, 001, 004, 003,-005, ! dx1(1:10, 8,18)
     .-022,-011,-029,-035,-020, 002,-012,-005,-004,-003, ! dx1(1:10, 9,18)
     . 006,-009,-019,-018,-017, 014,-012,-008, 001, 012, ! dx1(1:10,10,18)
     . 000, 003,-018,-016,-009, 023,-003,-008, 000, 010, ! dx1(1:10,11,18)
     . 003, 004,-010,-006, 002, 022,-006,-003, 009, 020, ! dx1(1:10,12,18)
     . 014, 014,-009,-007, 003, 014,-005,-013,-004, 011, ! dx1(1:10,13,18)
     . 039, 027, 012, 016, 016, 036, 014, 004, 019, 018, ! dx1(1:10,14,18)
     . 010,-010,-022,-009, 003, 005,-016,-006,-014, 010, ! dx1(1:10,15,18)
     . 028,-002,-008,-002, 004, 001,-010,-007,-015, 008, ! dx1(1:10,16,18)
     . 020, 007, 004, 010, 016, 000, 000, 006,-013, 015, ! dx1(1:10,17,18)
     . 035, 014, 012, 012, 026, 008, 015, 017,-008, 015, ! dx1(1:10,18,18)
     . 041, 033, 041, 024, 025, 027, 025, 041, 000, 022, ! dx1(1:10,19,18)
     . 039, 031, 034, 023, 027, 018, 024, 036, 000, 024, ! dx1(1:10,20,18)
     . 032, 042, 046, 017, 024, 009, 037, 013, 013, 005, ! dx1(1:10,21,18)
     . 020, 038, 042, 007, 015, 005, 039, 012, 017, 004, ! dx1(1:10,22,18)
     . 024, 030, 041, 015, 005, 010, 045, 016, 023, 000, ! dx1(1:10,23,18)
     . 007, 030, 049, 009, 006, 009, 050, 015, 022, 000, ! dx1(1:10,24,18)
     . 007,-002, 029, 017, 003,-006, 038, 014, 008,-011, ! dx1(1:10,25,18)
     .-034,-024,-005,-024,-045,-036, 010,-015,-019,-038, ! dx1(1:10,26,18)
     .-030,-023,-012,-036,-045,-043,-004,-028,-039,-042, ! dx1(1:10,27,18)
     .-028,-025,-020,-036,-046,-042,-011,-045,-039,-052, ! dx1(1:10,28,18)
     .-029,-039,-034,-054,-056,-055,-025,-059,-044,-048, ! dx1(1:10,29,18)
     .-038,-031,-041,-050,-049,-047,-033,-073,-043,-070, ! dx1(1:10,30,18)
     .-050,-039,-040,-065,-047,-022,-049,-096,-043,-070, ! dx1(1:10,31,18)
     .-013,-005,-024,-042,-012,-005,-015,-075,-012,-038, ! dx1(1:10,32,18)
     . 009, 006,-026,-041,-015, 014,-028,-090,-021,-043, ! dx1(1:10,33,18)
     .-017, 009,-035,-046,-021, 006,-027,-094,-019,-044, ! dx1(1:10,34,18)
     . 025, 020,-029,-039,-023, 001,-031,-112,-032,-050, ! dx1(1:10,35,18)
     .-010, 011,-046,-053,-049,-040,-046,-133,-046,-063, ! dx1(1:10,36,18)
     . 010, 027,-042,-070,-075,-114,-064,-182,-093,-080, ! dx1(1:10,37,18)
     .-027,-048,-075,-099,-120,-165,-098,-214,-139,-110/ ! dx1(1:10,38,18)

      data dx2 / 
     . 047, 016, 082, 134, 070, 021, 073, 051, 010, 008, ! dx2(1:10, 1, 1)
     . 038, 024, 081, 124, 067, 022, 072, 044, 000,-002, ! dx2(1:10, 2, 1)
     . 043, 040, 101, 142, 097, 056, 097, 054, 048, 032, ! dx2(1:10, 3, 1)
     . 089, 048, 093, 120, 083, 053, 096, 042, 062, 041, ! dx2(1:10, 4, 1)
     . 060, 052, 089, 112, 088, 058, 091, 054, 076, 048, ! dx2(1:10, 5, 1)
     . 037, 055, 065, 079, 072, 051, 090, 050, 070, 045, ! dx2(1:10, 6, 1)
     . 099, 079, 083, 101, 088, 076, 106, 079, 102, 076, ! dx2(1:10, 7, 1)
     . 038, 016, 012, 049, 039, 032, 051, 023, 044, 023, ! dx2(1:10, 8, 1)
     . 024, 025, 018, 030, 028, 036, 036, 022, 039, 033, ! dx2(1:10, 9, 1)
     . 024, 008, 002, 027, 018, 034, 026, 005, 035, 026, ! dx2(1:10,10, 1)
     . 024, 017, 011, 024, 027, 040, 025, 005, 037, 031, ! dx2(1:10,11, 1)
     . 031, 005,-007, 016, 021, 032, 010, 000, 033, 026, ! dx2(1:10,12, 1)
     .-004, 008, 000, 024, 018, 026, 008,-001, 027, 038, ! dx2(1:10,13, 1)
     .-018,-039,-053,-020,-028,-021,-040,-049,-021,-024, ! dx2(1:10,14, 1)
     .-029,-047,-059,-028,-025,-021,-043,-046,-022,-008, ! dx2(1:10,15, 1)
     .-027,-054,-054,-032,-026,-038,-044,-045,-040,-027, ! dx2(1:10,16, 1)
     .-013,-034,-052,-019,-020,-041,-047,-038,-038,-023, ! dx2(1:10,17, 1)
     .-016,-040,-037,-021,-014,-028,-037,-027,-047,-015, ! dx2(1:10,18, 1)
     .-001, 002, 000,-006, 011,-009,-009, 000,-006,-005, ! dx2(1:10,19, 1)
     . 011, 008, 001, 006, 014,-004,-008, 002,-010, 009, ! dx2(1:10,20, 1)
     . 032, 034, 033, 025, 032, 018, 025, 015, 020, 017, ! dx2(1:10,21, 1)
     . 028, 050, 044, 044, 036, 029, 040, 031, 038, 028, ! dx2(1:10,22, 1)
     . 033, 046, 054, 048, 045, 043, 055, 036, 050, 041, ! dx2(1:10,23, 1)
     . 055, 060, 075, 066, 060, 055, 076, 058, 060, 051, ! dx2(1:10,24, 1)
     . 025, 036, 050, 048, 036, 032, 057, 057, 040, 037, ! dx2(1:10,25, 1)
     .-001, 011, 021, 023, 018, 018, 030, 039, 017, 020, ! dx2(1:10,26, 1)
     . 003, 003, 014, 017, 005, 015, 025, 031, 002, 010, ! dx2(1:10,27, 1)
     .-008, 002, 004, 006,-002, 000, 005, 022,-003, 004, ! dx2(1:10,28, 1)
     . 001,-005, 002, 010, 000,-002,-001, 011,-014, 008, ! dx2(1:10,29, 1)
     .-008,-019,-006,-008,-007, 000,-021, 007,-017,-004, ! dx2(1:10,30, 1)
     .-007, 007,-016, 004, 005, 026,-025,-003,-010, 003, ! dx2(1:10,31, 1)
     .-015,-017,-015,-013,-005, 007,-039,-009,-018,-009, ! dx2(1:10,32, 1)
     .-001, 020,-013,-014, 011, 044,-026,-009,-013, 006, ! dx2(1:10,33, 1)
     .-014, 019,-006,-014, 006, 049,-032, 000,-016, 007, ! dx2(1:10,34, 1)
     . 024, 034, 006, 010, 019, 053,-018, 006,-009, 025, ! dx2(1:10,35, 1)
     .-014, 036, 017, 007, 007, 041,-012, 016,-006, 033, ! dx2(1:10,36, 1)
     . 064, 056, 055, 038, 019,-001, 005, 025,-022, 067, ! dx2(1:10,37, 1)
     .-013, 024, 024,-015,-017,-051,-027, 012,-029, 009, ! dx2(1:10,38, 1)
     . 007,-009, 063, 093, 047, 013, 051, 039,-016,-013, ! dx2(1:10, 1, 2)
     . 006, 006, 055, 097, 049, 008, 045, 030,-016,-018, ! dx2(1:10, 2, 2)
     . 022, 003, 051, 064, 035, 008, 045, 016, 006, 007, ! dx2(1:10, 3, 2)
     . 014, 014, 044, 070, 044, 011, 051, 015, 023, 006, ! dx2(1:10, 4, 2)
     . 022, 022, 028, 050, 038, 018, 052, 013, 035, 019, ! dx2(1:10, 5, 2)
     . 039, 029, 038, 046, 038, 024, 058, 021, 050, 030, ! dx2(1:10, 6, 2)
     . 038, 044, 040, 036, 048, 043, 066, 049, 068, 042, ! dx2(1:10, 7, 2)
     . 011, 012, 001,-003, 006,-001, 023, 005, 018, 005, ! dx2(1:10, 8, 2)
     . 002, 003,-011, 000, 001, 015, 013, 011, 024, 007, ! dx2(1:10, 9, 2)
     . 010, 010,-009,-001, 004, 026, 012, 012, 023, 015, ! dx2(1:10,10, 2)
     . 000, 004,-017, 002,-004, 021,-006, 013, 021, 017, ! dx2(1:10,11, 2)
     . 009, 001,-008,-005, 001, 019,-004, 010, 019, 014, ! dx2(1:10,12, 2)
     . 015,-004,-026,-007,-005, 004,-018, 011, 010, 008, ! dx2(1:10,13, 2)
     .-033,-033,-052,-045,-038,-021,-053,-023,-023,-016, ! dx2(1:10,14, 2)
     .-025,-056,-055,-046,-038,-043,-053,-036,-038,-035, ! dx2(1:10,15, 2)
     .-029,-040,-056,-039,-035,-044,-052,-040,-047,-030, ! dx2(1:10,16, 2)
     .-015,-046,-042,-035,-025,-036,-038,-034,-043,-020, ! dx2(1:10,17, 2)
     .-002,-013,-021,-017,-015,-037,-023,-021,-029,-014, ! dx2(1:10,18, 2)
     . 015, 001,-006,-001, 000,-024,-009,-013,-015,-001, ! dx2(1:10,19, 2)
     . 022, 007, 007, 010, 007,-017,-001,-007,-007, 001, ! dx2(1:10,20, 2)
     . 031, 028, 031, 030, 022, 007, 031, 005, 020, 021, ! dx2(1:10,21, 2)
     . 040, 035, 048, 036, 033, 027, 045, 016, 035, 030, ! dx2(1:10,22, 2)
     . 023, 035, 056, 043, 035, 025, 056, 032, 043, 035, ! dx2(1:10,23, 2)
     . 027, 038, 056, 045, 041, 035, 064, 042, 048, 036, ! dx2(1:10,24, 2)
     . 038, 043, 073, 055, 044, 049, 072, 064, 060, 055, ! dx2(1:10,25, 2)
     .-004,-003, 020, 009,-001, 006, 026, 015, 013, 003, ! dx2(1:10,26, 2)
     .-010,-006, 016, 012,-004,-003, 019, 027, 013, 001, ! dx2(1:10,27, 2)
     . 003,-011, 013, 007,-014,-008, 006, 020, 001,-003, ! dx2(1:10,28, 2)
     .-012,-013,-005, 001,-013,-010,-006, 020,-002,-005, ! dx2(1:10,29, 2)
     .-015,-016,-008,-001,-013,-015,-015, 015,-004,-012, ! dx2(1:10,30, 2)
     . 007,-007,-019,-010,-008, 011,-028, 002,-007,-013, ! dx2(1:10,31, 2)
     .-022,-013,-026,-020,-010, 005,-035,-002,-014,-011, ! dx2(1:10,32, 2)
     .-013,-009,-021,-015,-004, 021,-034,-007,-014,-009, ! dx2(1:10,33, 2)
     .-010, 006,-011,-011, 001, 036,-028,-013,-011,-001, ! dx2(1:10,34, 2)
     .-006, 004,-013,-015,-007, 029,-027,-018,-013, 002, ! dx2(1:10,35, 2)
     . 017, 015, 012,-003,-001, 027,-013,-022,-015, 027, ! dx2(1:10,36, 2)
     .-012, 029, 027, 002,-014,-026,-012,-034,-029, 026, ! dx2(1:10,37, 2)
     .-006,-010,-002,-012,-041,-062,-036,-062,-065, 015, ! dx2(1:10,38, 2)
     .-001,-015, 028, 058, 022,-017, 004, 020,-053,-021, ! dx2(1:10, 1, 3)
     .-015,-031, 014, 048, 003,-018, 008, 012,-051,-030, ! dx2(1:10, 2, 3)
     .-014,-007, 020, 035, 019,-014, 026, 000,-013,-015, ! dx2(1:10, 3, 3)
     . 014,-009, 016, 026, 005,-006, 022,-005, 007,-005, ! dx2(1:10, 4, 3)
     .-004,-008, 006, 017, 010,-009, 027,-008, 015, 005, ! dx2(1:10, 5, 3)
     . 020, 020, 022, 016, 025, 015, 046, 008, 044, 014, ! dx2(1:10, 6, 3)
     . 052, 033, 024, 016, 034, 023, 040, 029, 048, 031, ! dx2(1:10, 7, 3)
     .-006,-001,-023,-023,-010,-016,-001,-008, 007,-013, ! dx2(1:10, 8, 3)
     . 010, 003,-011,-011, 000, 008, 001, 016, 020, 000, ! dx2(1:10, 9, 3)
     .-001,-001,-024,-025,-014, 006,-015, 009, 004,-003, ! dx2(1:10,10, 3)
     . 002,-003,-029,-018,-016, 009,-018, 014, 004, 000, ! dx2(1:10,11, 3)
     . 007,-004,-023,-023,-001, 010,-018, 018, 009, 009, ! dx2(1:10,12, 3)
     . 009,-010,-024,-022,-017,-008,-029, 007,-005,-008, ! dx2(1:10,13, 3)
     .-055,-064,-075,-057,-059,-052,-075,-030,-049,-043, ! dx2(1:10,14, 3)
     .-009,-024,-043,-039,-035,-033,-049,-020,-036,-025, ! dx2(1:10,15, 3)
     .-018,-033,-040,-043,-034,-045,-051,-033,-051,-024, ! dx2(1:10,16, 3)
     .-007,-031,-040,-027,-031,-044,-039,-033,-040,-020, ! dx2(1:10,17, 3)
     .-007,-019,-020,-011,-010,-036,-021,-022,-028,-009, ! dx2(1:10,18, 3)
     . 024,-001,-002, 006,-009,-026,-005,-024,-014,-006, ! dx2(1:10,19, 3)
     . 030, 009,-001, 013, 000,-010, 002,-014,-005, 004, ! dx2(1:10,20, 3)
     . 012, 023, 030, 024, 019, 005, 031,-005, 015, 021, ! dx2(1:10,21, 3)
     . 032, 032, 044, 037, 022, 019, 048, 011, 030, 025, ! dx2(1:10,22, 3)
     . 031, 036, 052, 042, 034, 029, 058, 022, 042, 036, ! dx2(1:10,23, 3)
     . 032, 047, 062, 046, 037, 037, 067, 038, 045, 038, ! dx2(1:10,24, 3)
     . 039, 044, 067, 046, 036, 042, 074, 055, 058, 037, ! dx2(1:10,25, 3)
     .-005, 004, 021, 006,-010,-006, 022, 009, 010,-003, ! dx2(1:10,26, 3)
     . 001,-006, 017, 009,-008,-001, 019, 029, 019, 000, ! dx2(1:10,27, 3)
     .-015,-004, 009, 004,-011,-006, 013, 031, 011,-005, ! dx2(1:10,28, 3)
     . 007,-017, 000, 007,-020,-021,-009, 019, 008,-010, ! dx2(1:10,29, 3)
     .-019,-017,-012,-001,-015,-022,-017, 020, 006,-011, ! dx2(1:10,30, 3)
     . 000,-020,-017,-014,-016,-002,-030, 005,-009,-020, ! dx2(1:10,31, 3)
     .-011,-022,-021,-017,-018,-011,-038, 002,-013,-020, ! dx2(1:10,32, 3)
     .-001,-003,-011,-015,-004, 027,-032,-005,-006,-009, ! dx2(1:10,33, 3)
     .-004,-007,-019,-021,-010, 019,-039,-021,-021,-010, ! dx2(1:10,34, 3)
     .-007, 011,-012,-014,-004, 035,-028,-028,-014, 005, ! dx2(1:10,35, 3)
     .-017, 010,-007,-014,-011, 016,-021,-036,-015, 010, ! dx2(1:10,36, 3)
     . 003, 012, 011,-009,-019,-036,-019,-056,-035, 018, ! dx2(1:10,37, 3)
     .-038,-034,-014,-037,-057,-077,-053,-090,-072,-015, ! dx2(1:10,38, 3)
     .-038,-018, 016, 030, 007,-006, 008, 020,-038,-019, ! dx2(1:10, 1, 4)
     .-024,-029, 000, 028,-001,-024, 000, 013,-054,-034, ! dx2(1:10, 2, 4)
     .-012,-024, 010, 021, 001,-014, 012,-004,-018,-012, ! dx2(1:10, 3, 4)
     .-022,-021,-002, 005,-004,-018, 017,-024,-011,-013, ! dx2(1:10, 4, 4)
     .-003, 000,-005, 003, 001,-008, 021,-017, 009,-013, ! dx2(1:10, 5, 4)
     . 010, 006, 001, 001, 011,-002, 025,-009, 023, 003, ! dx2(1:10, 6, 4)
     . 008, 028, 019,-007, 015, 012, 028, 016, 040, 017, ! dx2(1:10, 7, 4)
     .-026,-026,-039,-047,-036,-034,-020,-030,-005,-031, ! dx2(1:10, 8, 4)
     . 007, 014,-018,-021,-003, 010,-003, 020, 011,-002, ! dx2(1:10, 9, 4)
     . 000, 000,-021,-027,-012, 000,-020, 018, 005,-004, ! dx2(1:10,10, 4)
     .-005,-004,-033,-039,-019,-007,-034, 011,-010,-017, ! dx2(1:10,11, 4)
     . 011, 000,-020,-024,-012, 005,-026, 028, 000,-003, ! dx2(1:10,12, 4)
     .-011,-010,-029,-030,-016,-011,-034, 014,-010,-013, ! dx2(1:10,13, 4)
     .-030,-051,-060,-070,-055,-051,-072,-021,-049,-046, ! dx2(1:10,14, 4)
     .-022,-036,-039,-047,-039,-041,-047,-022,-040,-035, ! dx2(1:10,15, 4)
     .-020,-034,-042,-035,-041,-045,-048,-028,-051,-037, ! dx2(1:10,16, 4)
     .-009,-017,-030,-021,-029,-041,-036,-027,-038,-014, ! dx2(1:10,17, 4)
     . 010,-008,-016,-004,-014,-024,-011,-021,-022,-008, ! dx2(1:10,18, 4)
     . 003,-005,-001, 000,-007,-026,-006,-032,-015, 001, ! dx2(1:10,19, 4)
     . 005, 000, 000, 005,-003,-026,-001,-032,-006,-001, ! dx2(1:10,20, 4)
     . 025, 017, 029, 021, 014, 007, 031,-013, 005, 019, ! dx2(1:10,21, 4)
     . 036, 028, 040, 023, 019, 013, 043,-007, 023, 023, ! dx2(1:10,22, 4)
     . 016, 022, 039, 022, 017, 018, 050, 000, 020, 016, ! dx2(1:10,23, 4)
     . 026, 040, 064, 039, 029, 032, 069, 021, 040, 040, ! dx2(1:10,24, 4)
     . 027, 032, 050, 035, 018, 025, 051, 040, 045, 024, ! dx2(1:10,25, 4)
     .-007,-003, 018, 001,-015,-008, 016, 005, 008,-008, ! dx2(1:10,26, 4)
     .-008, 000, 017, 009,-010,-005, 021, 025, 024,-008, ! dx2(1:10,27, 4)
     . 000,-012, 007, 005,-020,-016, 000, 025, 016,-014, ! dx2(1:10,28, 4)
     .-020,-025,-003, 000,-021,-023,-005, 026, 009,-020, ! dx2(1:10,29, 4)
     .-013,-016,-010,-010,-023,-031,-022, 020, 009,-024, ! dx2(1:10,30, 4)
     .-012,-020,-019,-017,-014,-010,-028, 020,-004,-018, ! dx2(1:10,31, 4)
     .-016,-020,-025,-029,-023,-009,-035, 007,-014,-028, ! dx2(1:10,32, 4)
     .-012,-002,-017,-007,-003, 009,-033, 000,-005,-011, ! dx2(1:10,33, 4)
     . 000,-009,-023,-019,-010, 021,-038,-016,-018,-008, ! dx2(1:10,34, 4)
     .-002, 008,-013,-012,-007, 023,-031,-025,-015, 000, ! dx2(1:10,35, 4)
     .-022,-001,-018,-019,-016, 021,-024,-047,-023, 007, ! dx2(1:10,36, 4)
     . 021, 033, 029, 003,-006,-021, 004,-058,-011, 030, ! dx2(1:10,37, 4)
     .-024,-011,-011,-031,-044,-074,-045,-106,-054,-012, ! dx2(1:10,38, 4)
     .-006,-011, 010, 028, 009,-041, 000, 023,-055,-020, ! dx2(1:10, 1, 5)
     .-030,-032,-006, 011,-013,-032,-015, 005,-057,-043, ! dx2(1:10, 2, 5)
     .-006,-027, 000, 002,-007,-021, 008,-011,-025,-013, ! dx2(1:10, 3, 5)
     .-030,-019,-004, 000,-006,-010, 010,-023,-008,-011, ! dx2(1:10, 4, 5)
     .-004,-007,-003,-008,-004,-016, 015,-029, 006,-005, ! dx2(1:10, 5, 5)
     . 013, 006, 008,-009, 008, 000, 024,-019, 023, 000, ! dx2(1:10, 6, 5)
     . 015, 011, 002,-011, 013, 009, 015, 009, 034,-001, ! dx2(1:10, 7, 5)
     .-037,-033,-047,-060,-042,-040,-035,-039,-020,-037, ! dx2(1:10, 8, 5)
     .-004, 003,-014,-033,-008,-005,-014, 018, 012,-016, ! dx2(1:10, 9, 5)
     .-007, 002,-027,-035,-013, 004,-025, 018,-002,-013, ! dx2(1:10,10, 5)
     . 000,-008,-038,-041,-029,-005,-043, 011,-016,-016, ! dx2(1:10,11, 5)
     .-002, 001,-020,-033,-012,-003,-033, 033,-008,-013, ! dx2(1:10,12, 5)
     .-001,-022,-022,-039,-025,-020,-039, 023,-015,-010, ! dx2(1:10,13, 5)
     .-048,-047,-066,-073,-058,-052,-073,-013,-055,-051, ! dx2(1:10,14, 5)
     .-029,-034,-037,-044,-049,-046,-058,-023,-046,-036, ! dx2(1:10,15, 5)
     .-017,-034,-040,-035,-043,-047,-046,-029,-047,-034, ! dx2(1:10,16, 5)
     . 003,-020,-025,-021,-034,-034,-032,-030,-039,-026, ! dx2(1:10,17, 5)
     . 000,-008,-015,-008,-015,-025,-012,-024,-020,-007, ! dx2(1:10,18, 5)
     . 026, 000, 004, 010,-003,-012, 012,-035,-011, 010, ! dx2(1:10,19, 5)
     . 021, 004, 008, 006,-006,-011, 008,-037,-011, 010, ! dx2(1:10,20, 5)
     . 029, 016, 025, 016, 006, 002, 028,-027, 006, 017, ! dx2(1:10,21, 5)
     . 024, 015, 026, 018, 009, 006, 038,-025, 001, 017, ! dx2(1:10,22, 5)
     . 023, 027, 038, 025, 010, 011, 044,-017, 014, 023, ! dx2(1:10,23, 5)
     . 023, 028, 050, 033, 018, 018, 051, 003, 024, 019, ! dx2(1:10,24, 5)
     . 042, 045, 055, 043, 024, 032, 060, 039, 038, 027, ! dx2(1:10,25, 5)
     . 012,-009, 019, 008,-021,-006, 018, 001, 000,-011, ! dx2(1:10,26, 5)
     .-008,-005, 016, 009,-009,-007, 011, 022, 018,-012, ! dx2(1:10,27, 5)
     .-011,-019, 002,-002,-027,-021,-002, 022, 015,-017, ! dx2(1:10,28, 5)
     .-034,-029,-016,-016,-028,-034,-018, 020, 013,-029, ! dx2(1:10,29, 5)
     .-009,-018,-008,-002,-017,-021,-017, 034, 019,-017, ! dx2(1:10,30, 5)
     .-012,-022,-024,-017,-021,-023,-028, 018,-002,-026, ! dx2(1:10,31, 5)
     .-008,-020,-038,-015,-021,-026,-039, 013,-007,-028, ! dx2(1:10,32, 5)
     .-009, 001,-013,-015,-006, 019,-025, 007,-005,-011, ! dx2(1:10,33, 5)
     .-011, 008,-019,-013,-003, 014,-029,-009,-014,-001, ! dx2(1:10,34, 5)
     .-018,-010,-016,-008,-010, 026,-026,-026,-023,-004, ! dx2(1:10,35, 5)
     . 010, 006,-001, 000,-001, 028,-008,-035,-003, 019, ! dx2(1:10,36, 5)
     . 009, 019, 023, 009,-003,-013, 001,-062,-015, 032, ! dx2(1:10,37, 5)
     .-048,-029,-003,-024,-042,-054,-020,-093,-042, 006, ! dx2(1:10,38, 5)
     .-005,-001, 006, 034, 008,-009, 006, 036,-037,-009, ! dx2(1:10, 1, 6)
     .-015,-032,-005, 007,-008,-027, 000, 015,-050,-018, ! dx2(1:10, 2, 6)
     .-032,-028,-009,-006,-013,-019,-001,-019,-023,-021, ! dx2(1:10, 3, 6)
     .-013,-016,-007,-012,-011,-021, 005,-029,-012,-008, ! dx2(1:10, 4, 6)
     .-009,-013,-003,-008,-005,-010, 015,-033, 003,-006, ! dx2(1:10, 5, 6)
     . 006, 020, 012,-001, 015, 003, 024,-010, 032, 012, ! dx2(1:10, 6, 6)
     . 008,-002,-011,-022, 000, 004, 001,-002, 018, 000, ! dx2(1:10, 7, 6)
     .-033,-031,-040,-060,-033,-037,-031,-041,-013,-044, ! dx2(1:10, 8, 6)
     .-002, 000,-029,-036,-012,-005,-027, 018, 002,-013, ! dx2(1:10, 9, 6)
     .-001, 002,-025,-037,-014,-008,-032, 017,-008,-014, ! dx2(1:10,10, 6)
     . 001,-007,-041,-053,-028,-016,-043, 018,-017,-028, ! dx2(1:10,11, 6)
     . 003,-003,-018,-036,-007, 001,-029, 038,-007,-015, ! dx2(1:10,12, 6)
     .-004,-011,-019,-023,-021,-012,-035, 028,-019,-005, ! dx2(1:10,13, 6)
     .-025,-056,-056,-059,-061,-049,-070,-009,-056,-052, ! dx2(1:10,14, 6)
     .-014,-026,-039,-035,-039,-041,-050,-015,-047,-029, ! dx2(1:10,15, 6)
     .-025,-036,-032,-034,-038,-041,-041,-028,-045,-031, ! dx2(1:10,16, 6)
     .-014,-024,-028,-022,-031,-037,-031,-033,-035,-012, ! dx2(1:10,17, 6)
     . 000,-012,-015,-012,-022,-021,-008,-032,-022,-003, ! dx2(1:10,18, 6)
     . 017, 002, 004,-002,-006,-018, 006,-045,-016, 001, ! dx2(1:10,19, 6)
     . 028, 004, 000, 008, 000,-011, 013,-039,-008, 007, ! dx2(1:10,20, 6)
     . 004, 017, 025, 018, 006, 004, 029,-040,-003, 013, ! dx2(1:10,21, 6)
     . 022, 012, 028, 022, 006, 001, 036,-033, 000, 014, ! dx2(1:10,22, 6)
     . 025, 028, 029, 027, 008, 013, 041,-024, 005, 018, ! dx2(1:10,23, 6)
     . 024, 035, 049, 038, 017, 019, 056, 001, 020, 018, ! dx2(1:10,24, 6)
     . 026, 022, 051, 039, 015, 022, 046, 024, 035, 015, ! dx2(1:10,25, 6)
     .-017,-007, 010,-003,-025,-017, 004,-010,-001,-017, ! dx2(1:10,26, 6)
     .-002,-012, 014, 004,-014,-014, 009, 021, 016,-018, ! dx2(1:10,27, 6)
     .-003,-016,-004,-005,-020,-015,-005, 025, 022,-022, ! dx2(1:10,28, 6)
     . 000,-024,-013,-010,-026,-030,-018, 028, 018,-023, ! dx2(1:10,29, 6)
     .-008,-008,-007, 002,-012,-029,-010, 037, 023,-014, ! dx2(1:10,30, 6)
     .-012,-017,-022,-013,-018,-032,-030, 029, 000,-020, ! dx2(1:10,31, 6)
     .-003,-022,-037,-019,-023,-031,-034, 016,-004,-035, ! dx2(1:10,32, 6)
     .-008,-010,-016,-008,-002, 007,-026, 020,-005,-011, ! dx2(1:10,33, 6)
     .-003,-010,-026,-018,-005, 016,-029,-003,-017,-012, ! dx2(1:10,34, 6)
     .-004, 006,-014,-004, 000, 031,-020,-010,-007, 005, ! dx2(1:10,35, 6)
     . 006, 019, 014, 011, 011, 038,-002,-021, 000, 031, ! dx2(1:10,36, 6)
     . 014, 017, 040, 028, 011, 005, 021,-052,-004, 044, ! dx2(1:10,37, 6)
     . 003, 000, 015,-009,-018,-049,-021,-096,-036, 011, ! dx2(1:10,38, 6)
     . 007,-005, 011, 026, 014, 000, 010, 042,-018, 010, ! dx2(1:10, 1, 7)
     .-029,-013,-007, 015,-002, 001, 000, 026,-034,-002, ! dx2(1:10, 2, 7)
     .-019,-023,-002,-012,-007,-019, 007,-017,-021,-014, ! dx2(1:10, 3, 7)
     .-025,-022,-008,-012,-014,-018, 003,-037,-004,-010, ! dx2(1:10, 4, 7)
     .-020,-006,-005,-018,-006,-009, 007,-036, 008,-004, ! dx2(1:10, 5, 7)
     . 016, 011, 006,-008, 012, 003, 019,-016, 032, 005, ! dx2(1:10, 6, 7)
     . 016, 018, 010,-006, 021, 011, 015, 017, 036, 007, ! dx2(1:10, 7, 7)
     .-012,-016,-040,-054,-033,-040,-037,-032,-016,-035, ! dx2(1:10, 8, 7)
     .-001, 005,-014,-038, 000,-003,-019, 021, 004,-016, ! dx2(1:10, 9, 7)
     .-002,-002,-018,-049,-011,-001,-027, 029,-001,-024, ! dx2(1:10,10, 7)
     .-010,-002,-028,-045,-019,-010,-042, 025,-015,-025, ! dx2(1:10,11, 7)
     . 009, 004,-013,-030,-012,-005,-028, 038,-008,-011, ! dx2(1:10,12, 7)
     . 008, 001,-013,-010,-012,-004,-029, 033,-014,-009, ! dx2(1:10,13, 7)
     .-029,-041,-054,-051,-054,-046,-066,-008,-053,-049, ! dx2(1:10,14, 7)
     .-027,-033,-030,-033,-040,-036,-043,-010,-044,-027, ! dx2(1:10,15, 7)
     .-012,-036,-034,-042,-035,-037,-029,-026,-037,-027, ! dx2(1:10,16, 7)
     .-016,-019,-023,-028,-030,-027,-023,-031,-038,-012, ! dx2(1:10,17, 7)
     . 001,-006,-003,-013,-015,-018,-001,-037,-016,-005, ! dx2(1:10,18, 7)
     . 017,-009, 000, 004,-011,-017, 008,-052,-017, 005, ! dx2(1:10,19, 7)
     . 011, 007, 008, 010,-005,-006, 022,-047,-013, 009, ! dx2(1:10,20, 7)
     . 029, 017, 023, 029, 009, 000, 033,-044, 000, 021, ! dx2(1:10,21, 7)
     . 014, 024, 026, 023, 007, 010, 039,-043, 000, 016, ! dx2(1:10,22, 7)
     . 036, 028, 036, 035, 013, 013, 044,-025, 009, 023, ! dx2(1:10,23, 7)
     . 030, 029, 050, 035, 017, 025, 046,-003, 019, 024, ! dx2(1:10,24, 7)
     . 018, 025, 041, 034, 016, 019, 037, 022, 029, 014, ! dx2(1:10,25, 7)
     .-010,-006, 009, 005,-013,-009, 008,-006,-002,-011, ! dx2(1:10,26, 7)
     .-002,-018, 011, 002,-015,-009, 000, 017, 017,-014, ! dx2(1:10,27, 7)
     .-004,-019, 000, 001,-013,-012,-008, 026, 024,-015, ! dx2(1:10,28, 7)
     .-006,-010,-009,-001,-017,-029,-016, 033, 023,-029, ! dx2(1:10,29, 7)
     . 000,-016,-016, 000,-007,-021,-012, 041, 029,-016, ! dx2(1:10,30, 7)
     .-025,-021,-026,-011,-014,-036,-028, 027, 003,-023, ! dx2(1:10,31, 7)
     .-008,-019,-033,-022,-018,-038,-034, 024,-007,-026, ! dx2(1:10,32, 7)
     .-009,-016,-024,-009,-005, 003,-018, 029,-002,-011, ! dx2(1:10,33, 7)
     .-006,-007,-016,-004,-001, 022,-022, 009,-003, 000, ! dx2(1:10,34, 7)
     . 012, 007,-002, 011, 009, 040,-009, 007, 004, 014, ! dx2(1:10,35, 7)
     . 031, 017, 013, 020, 019, 054, 017,-008, 013, 036, ! dx2(1:10,36, 7)
     . 022, 026, 050, 038, 018, 015, 031,-041, 007, 054, ! dx2(1:10,37, 7)
     .-020,-014, 011,-007,-020,-034,-010,-090,-030, 022, ! dx2(1:10,38, 7)
     . 023, 006, 008, 030, 015, 024, 019, 026,-015, 009, ! dx2(1:10, 1, 8)
     .-013,-015, 008, 015, 006, 002, 013, 024,-025, 013, ! dx2(1:10, 2, 8)
     .-018,-003, 011, 001,-001, 000, 009,-009, 000, 000, ! dx2(1:10, 3, 8)
     .-030,-011,-001,-016,-007,-016, 007,-037,-004,-009, ! dx2(1:10, 4, 8)
     . 004,-013, 000,-009,-001,-006, 012,-033, 017, 000, ! dx2(1:10, 5, 8)
     . 008, 013, 017,-004, 019, 011, 020,-005, 031, 008, ! dx2(1:10, 6, 8)
     . 026, 038, 011,-002, 026, 016, 024, 019, 044, 011, ! dx2(1:10, 7, 8)
     .-015,-019,-028,-041,-015,-025,-024,-023,-006,-024, ! dx2(1:10, 8, 8)
     . 004, 010,-005,-032, 001,-004,-018, 028, 010,-013, ! dx2(1:10, 9, 8)
     .-013,-001,-020,-030,-014,-004,-025, 037, 000,-014, ! dx2(1:10,10, 8)
     . 002,-008,-020,-032,-017,-011,-036, 037,-011,-021, ! dx2(1:10,11, 8)
     .-008, 013,-014,-016,-006, 000,-030, 041,-009,-010, ! dx2(1:10,12, 8)
     . 000,-010,-013,-025,-018,-015,-025, 027,-015,-014, ! dx2(1:10,13, 8)
     .-035,-031,-038,-049,-038,-033,-056,-003,-046,-035, ! dx2(1:10,14, 8)
     .-005,-025,-029,-041,-031,-030,-031,-010,-037,-023, ! dx2(1:10,15, 8)
     .-021,-023,-022,-039,-029,-029,-031,-024,-036,-022, ! dx2(1:10,16, 8)
     .-025,-020,-021,-027,-028,-026,-015,-038,-033,-006, ! dx2(1:10,17, 8)
     . 007,-012,-002,-006,-011,-010, 003,-040,-015, 005, ! dx2(1:10,18, 8)
     . 002, 003,-004, 016,-014,-007, 017,-057,-021, 009, ! dx2(1:10,19, 8)
     . 009, 004, 003, 020,-006,-009, 018,-058,-017, 010, ! dx2(1:10,20, 8)
     . 035, 025, 026, 029, 007, 009, 038,-044,-002, 022, ! dx2(1:10,21, 8)
     . 024, 021, 031, 032, 010, 013, 041,-040, 001, 023, ! dx2(1:10,22, 8)
     . 034, 027, 037, 026, 014, 017, 042,-030, 009, 020, ! dx2(1:10,23, 8)
     . 030, 029, 042, 035, 018, 025, 046,-003, 018, 032, ! dx2(1:10,24, 8)
     . 041, 034, 057, 046, 025, 035, 046, 024, 033, 016, ! dx2(1:10,25, 8)
     .-007,-016, 010, 001,-019,-003, 001,-015,-006,-014, ! dx2(1:10,26, 8)
     . 000,-009, 006, 009,-004,-004, 004, 023, 018,-005, ! dx2(1:10,27, 8)
     . 000,-012, 000, 004,-010,-005,-007, 028, 023,-015, ! dx2(1:10,28, 8)
     . 002,-013,-015, 002,-006,-013,-011, 034, 025,-015, ! dx2(1:10,29, 8)
     .-003,-006,-005, 010, 000,-020,-011, 045, 025,-009, ! dx2(1:10,30, 8)
     . 007,-014,-033,-008,-008,-023,-025, 037, 005,-009, ! dx2(1:10,31, 8)
     .-020,-031,-037,-021,-019,-039,-031, 031,-007,-030, ! dx2(1:10,32, 8)
     . 015,-003,-013, 003, 007, 008,-013, 043, 011, 006, ! dx2(1:10,33, 8)
     . 008,-001,-015, 000, 004, 023,-006, 027, 008, 003, ! dx2(1:10,34, 8)
     . 023, 012, 012, 022, 017, 051, 010, 022, 016, 028, ! dx2(1:10,35, 8)
     . 026, 028, 029, 033, 027, 061, 025, 010, 020, 046, ! dx2(1:10,36, 8)
     . 032, 030, 067, 049, 036, 024, 035,-031, 020, 054, ! dx2(1:10,37, 8)
     . 015, 000, 031, 010, 004,-019, 006,-056,-017, 034, ! dx2(1:10,38, 8)
     .-008, 006, 010, 014, 016, 023, 022, 036,-018, 017, ! dx2(1:10, 1, 9)
     . 008, 006, 014, 017, 007, 009, 017, 023,-014, 014, ! dx2(1:10, 2, 9)
     .-012,-009, 019, 013, 010, 009, 024,-008, 006, 017, ! dx2(1:10, 3, 9)
     .-009,-016, 001,-004,-001,-003, 014,-025, 013, 001, ! dx2(1:10, 4, 9)
     .-020, 007, 007,-001, 005,-003, 014,-028, 017,-002, ! dx2(1:10, 5, 9)
     . 016, 030, 022, 007, 030, 022, 030, 006, 049, 024, ! dx2(1:10, 6, 9)
     . 029, 018, 017,-002, 029, 021, 013, 034, 039, 016, ! dx2(1:10, 7, 9)
     .-030,-010,-032,-054,-016,-031,-025,-016,-003,-026, ! dx2(1:10, 8, 9)
     . 012, 006,-001,-022, 007, 002,-009, 044, 017, 004, ! dx2(1:10, 9, 9)
     . 003, 001,-006,-026,-003,-001,-020, 051, 005,-012, ! dx2(1:10,10, 9)
     . 000, 009,-010,-019,-008,-001,-026, 046,-004,-012, ! dx2(1:10,11, 9)
     . 017, 006, 000,-012,-002, 008,-018, 050, 000,-003, ! dx2(1:10,12, 9)
     . 006, 018,-006,-017, 000, 008,-007, 034,-006, 000, ! dx2(1:10,13, 9)
     .-028,-040,-045,-062,-044,-043,-053,-010,-050,-052, ! dx2(1:10,14, 9)
     .-022,-010,-016,-036,-021,-018,-026,-006,-030,-006, ! dx2(1:10,15, 9)
     .-012,-021,-024,-030,-026,-022,-017,-028,-036,-016, ! dx2(1:10,16, 9)
     .-013,-020,-007,-016,-016,-015,-005,-037,-027,-002, ! dx2(1:10,17, 9)
     .-004,-001, 004, 007,-006,-001, 010,-038,-019, 017, ! dx2(1:10,18, 9)
     . 013,-011, 000, 013,-004, 000, 023,-057,-019, 009, ! dx2(1:10,19, 9)
     . 010, 000, 004, 016,-001, 000, 022,-057,-017, 009, ! dx2(1:10,20, 9)
     . 023, 014, 023, 023, 010, 012, 033,-049,-005, 021, ! dx2(1:10,21, 9)
     . 022, 019, 030, 023, 011, 018, 040,-043,-001, 018, ! dx2(1:10,22, 9)
     . 016, 025, 036, 034, 020, 022, 043,-027, 009, 027, ! dx2(1:10,23, 9)
     . 036, 039, 052, 043, 026, 029, 049,-003, 022, 027, ! dx2(1:10,24, 9)
     . 037, 030, 047, 053, 033, 034, 045, 031, 030, 034, ! dx2(1:10,25, 9)
     . 005, 000, 010, 014,-009,-007, 002,-011,-003,-011, ! dx2(1:10,26, 9)
     . 009, 005, 015, 019, 007, 009, 004, 026, 017, 003, ! dx2(1:10,27, 9)
     . 008,-006,-003, 010, 000, 000,-007, 032, 021,-012, ! dx2(1:10,28, 9)
     . 007,-008,-009, 009, 000,-017,-014, 039, 024,-012, ! dx2(1:10,29, 9)
     . 000,-006,-016, 002, 003,-016,-015, 043, 026,-014, ! dx2(1:10,30, 9)
     .-024,-019,-022, 004,-002,-030,-020, 042, 012,-012, ! dx2(1:10,31, 9)
     .-010,-015,-037,-008,-007,-039,-030, 037, 002,-014, ! dx2(1:10,32, 9)
     .-005,-002,-015, 003, 011, 010, 000, 050, 016,-001, ! dx2(1:10,33, 9)
     . 000, 005,-008, 003, 006, 023, 002, 037, 017, 015, ! dx2(1:10,34, 9)
     . 021, 008, 003, 023, 015, 049, 014, 029, 019, 024, ! dx2(1:10,35, 9)
     . 027, 039, 041, 037, 035, 064, 036, 013, 030, 047, ! dx2(1:10,36, 9)
     . 032, 041, 078, 054, 046, 044, 056,-023, 025, 065, ! dx2(1:10,37, 9)
     . 015, 024, 054, 019, 015,-009, 030,-063,-015, 032, ! dx2(1:10,38, 9)
     . 028, 021, 020, 041, 029, 038, 033, 026, 005, 029, ! dx2(1:10, 1,10)
     . 004, 012, 008, 032, 019, 032, 017, 014,-015, 011, ! dx2(1:10, 2,10)
     .-007, 008, 028, 013, 024, 021, 029, 000, 018, 027, ! dx2(1:10, 3,10)
     .-002, 022, 017, 011, 014, 019, 026,-012, 021, 014, ! dx2(1:10, 4,10)
     .-003,-001, 006,-001, 009, 009, 014,-011, 030, 005, ! dx2(1:10, 5,10)
     . 019, 014, 019, 010, 028, 021, 024, 012, 044, 021, ! dx2(1:10, 6,10)
     . 008, 031, 028, 010, 041, 026, 036, 057, 058, 022, ! dx2(1:10, 7,10)
     .-012,-030,-023,-033,-016,-019,-026, 000, 000,-021, ! dx2(1:10, 8,10)
     . 004, 004, 006, 000, 011, 004,-005, 061, 025,-006, ! dx2(1:10, 9,10)
     . 000, 005, 008, 010, 008, 011,-007, 068, 018,-008, ! dx2(1:10,10,10)
     . 000,-005, 000,-014, 000, 000,-018, 055, 007, 001, ! dx2(1:10,11,10)
     . 015, 023, 012,-015, 013, 017,-002, 057, 012, 006, ! dx2(1:10,12,10)
     . 004, 011, 011,-012, 009, 009, 001, 044, 002, 017, ! dx2(1:10,13,10)
     .-020,-005,-022,-054,-023,-015,-032, 009,-033,-014, ! dx2(1:10,14,10)
     .-001,-005,-009,-016,-011,-008,-004,-001,-018,-013, ! dx2(1:10,15,10)
     .-023, 001, 002,-006,-002,-002, 000,-017,-019, 000, ! dx2(1:10,16,10)
     .-003,-004,-004, 010,-008, 003, 004,-035,-019,-002, ! dx2(1:10,17,10)
     . 011, 005, 005, 015, 005, 015, 022,-037,-008, 015, ! dx2(1:10,18,10)
     . 017, 010, 021, 011, 005, 014, 026,-051,-012, 020, ! dx2(1:10,19,10)
     . 028, 003, 019, 010, 007, 012, 027,-046,-010, 036, ! dx2(1:10,20,10)
     . 025, 022, 024, 020, 013, 018, 042,-047, 002, 023, ! dx2(1:10,21,10)
     . 028, 015, 033, 035, 018, 022, 040,-038, 007, 029, ! dx2(1:10,22,10)
     . 034, 034, 028, 040, 018, 021, 034,-025, 002, 028, ! dx2(1:10,23,10)
     . 032, 041, 035, 051, 031, 030, 045,-001, 017, 030, ! dx2(1:10,24,10)
     . 031, 052, 059, 057, 047, 045, 042, 029, 031, 037, ! dx2(1:10,25,10)
     .-003, 003, 008, 018, 008, 005, 007,-010,-010,-004, ! dx2(1:10,26,10)
     . 021, 014, 014, 025, 021, 013, 009, 022, 019, 000, ! dx2(1:10,27,10)
     . 015, 007, 003, 015, 013, 002,-002, 025, 018, 000, ! dx2(1:10,28,10)
     . 009, 000,-008, 007, 011, 007,-008, 036, 024,-006, ! dx2(1:10,29,10)
     . 004, 007,-007, 007, 010,-011,-010, 049, 032,-005, ! dx2(1:10,30,10)
     . 023, 000,-024,-002, 009,-009,-004, 057, 030, 004, ! dx2(1:10,31,10)
     .-008,-023,-029,-015,-004,-033,-024, 051, 016,-019, ! dx2(1:10,32,10)
     . 019, 008,-014, 006, 012,-003, 011, 061, 024, 022, ! dx2(1:10,33,10)
     . 010, 009, 004, 014, 021, 021, 019, 049, 022, 015, ! dx2(1:10,34,10)
     . 032, 022, 022, 023, 031, 046, 036, 034, 034, 027, ! dx2(1:10,35,10)
     . 035, 021, 052, 044, 040, 082, 046, 016, 038, 045, ! dx2(1:10,36,10)
     . 069, 074, 106, 079, 073, 066, 075,-014, 044, 071, ! dx2(1:10,37,10)
     . 058, 010, 063, 052, 024, 021, 030,-054,-001, 031, ! dx2(1:10,38,10)
     . 028, 005, 029, 009, 005, 015, 015, 009,-022, 007, ! dx2(1:10, 1,11)
     .-019,-011, 003, 004, 009, 008, 004, 001,-028, 007, ! dx2(1:10, 2,11)
     .-017, 000, 006,-005, 003, 009, 016,-014, 001, 001, ! dx2(1:10, 3,11)
     .-024,-002, 001,-003, 005,-002, 010,-026, 005, 002, ! dx2(1:10, 4,11)
     .-027, 010, 001,-009, 008, 000, 004,-017, 023, 002, ! dx2(1:10, 5,11)
     .-008, 002, 005,-011, 010, 004, 010, 002, 032, 009, ! dx2(1:10, 6,11)
     . 001, 006, 008, 000, 022, 009, 008, 029, 039, 004, ! dx2(1:10, 7,11)
     .-018,-019,-029,-031,-015,-028,-025,-006,-001,-036, ! dx2(1:10, 8,11)
     . 001,-008, 000, 001, 005, 000,-015, 050, 013,-004, ! dx2(1:10, 9,11)
     .-010, 008, 000,-009, 003,-001,-012, 060, 015,-008, ! dx2(1:10,10,11)
     .-008,-005,-008,-023,-006,-002,-020, 053,-002,-005, ! dx2(1:10,11,11)
     . 020, 018, 005,-016, 008, 011,-006, 053, 008, 002, ! dx2(1:10,12,11)
     .-006, 014, 006,-027, 001, 012,-004, 031,-005, 005, ! dx2(1:10,13,11)
     .-028,-014,-023,-050,-025,-025,-027, 007,-034,-021, ! dx2(1:10,14,11)
     .-013,-006,-011,-017,-015,-016,-012,-003,-019,-009, ! dx2(1:10,15,11)
     .-010,-008,-005,-008,-013,-010,-005,-028,-026,-008, ! dx2(1:10,16,11)
     .-004,-010,-006, 000,-010,-006, 001,-037,-026, 003, ! dx2(1:10,17,11)
     .-001, 004, 004, 002,-005, 003, 016,-044,-009, 009, ! dx2(1:10,18,11)
     . 016, 005, 013,-002, 009, 017, 027,-051,-014, 013, ! dx2(1:10,19,11)
     .-007,-018,-013,-020,-016,-011, 004,-078,-031,-011, ! dx2(1:10,20,11)
     .-018,-013,-005,-007,-018,-016, 006,-081,-032,-003, ! dx2(1:10,21,11)
     .-004,-004,-005, 000,-016,-017, 002,-074,-030,-007, ! dx2(1:10,22,11)
     .-017,-006, 000, 010,-013,-006, 006,-053,-025, 004, ! dx2(1:10,23,11)
     .-005, 008, 001, 013, 000,-010, 010,-033,-019,-008, ! dx2(1:10,24,11)
     . 022, 014, 018, 027, 010, 016, 009, 000, 000,-005, ! dx2(1:10,25,11)
     . 000,-006,-003, 005,-007,-011,-007,-023,-021,-021, ! dx2(1:10,26,11)
     .-012,-005,-002, 004, 006,-002,-011, 001,-003,-018, ! dx2(1:10,27,11)
     . 001,-008,-013,-005,-001,-005,-020, 007, 002,-027, ! dx2(1:10,28,11)
     .-008,-019,-029,-007,-013,-018,-028, 021, 004,-024, ! dx2(1:10,29,11)
     . 001,-016,-034,-010, 000,-031,-025, 028, 012,-023, ! dx2(1:10,30,11)
     .-006,-020,-038,-014,-009,-031,-027, 038, 006,-008, ! dx2(1:10,31,11)
     .-028,-049,-076,-050,-041,-060,-057, 008,-020,-047, ! dx2(1:10,32,11)
     .-033,-047,-060,-040,-029,-046,-039, 016,-022,-037, ! dx2(1:10,33,11)
     .-021,-035,-039,-027,-031,-019,-025, 005,-018,-022, ! dx2(1:10,34,11)
     .-029,-019,-021,-014,-013, 007,-009,-005,-011,-018, ! dx2(1:10,35,11)
     . 016,-011, 009, 002,-004, 027, 003,-026,-002, 011, ! dx2(1:10,36,11)
     .-018, 013, 051, 023, 025, 017, 024,-063,-009, 025, ! dx2(1:10,37,11)
     . 046,-001, 056, 018, 023, 023, 021,-044,-011, 012, ! dx2(1:10,38,11)
     . 000, 027,-022, 021, 013, 031, 014, 000,-015, 052, ! dx2(1:10, 1,12)
     .-002,-004,-011, 006,-003, 019, 006,-012,-027, 013, ! dx2(1:10, 2,12)
     .-024,-003, 001, 003, 005, 014, 014,-014,-001, 022, ! dx2(1:10, 3,12)
     .-020,-005, 000,-006, 004, 007,-001,-016, 006, 003, ! dx2(1:10, 4,12)
     .-013,-007,-004,-003, 003, 002, 001,-007, 016, 002, ! dx2(1:10, 5,12)
     .-014, 018, 019, 015, 027, 018, 017, 016, 042, 012, ! dx2(1:10, 6,12)
     . 005, 014, 009, 020, 016, 008, 010, 042, 039, 009, ! dx2(1:10, 7,12)
     .-053,-040,-025,-010,-021,-033,-034, 000,-010,-039, ! dx2(1:10, 8,12)
     .-012,-006, 003,-006, 002,-003,-008, 053, 017,-007, ! dx2(1:10, 9,12)
     .-011,-006, 003,-018, 000, 004,-014, 059, 007,-014, ! dx2(1:10,10,12)
     .-011,-006,-006,-028,-007,-001,-020, 053,-002,-009, ! dx2(1:10,11,12)
     . 003, 016, 016,-011, 013, 011, 002, 057, 007, 008, ! dx2(1:10,12,12)
     . 027, 025, 022, 011, 016, 019, 014, 048, 018, 015, ! dx2(1:10,13,12)
     .-009,-013,-016,-022,-022,-018,-026, 006,-024,-020, ! dx2(1:10,14,12)
     . 014, 021, 017, 012, 010, 011, 019, 016, 000, 008, ! dx2(1:10,15,12)
     .-015, 006, 004,-009,-005, 002, 001,-017,-018, 003, ! dx2(1:10,16,12)
     .-004, 000, 002,-004,-001, 007, 012,-032,-011, 005, ! dx2(1:10,17,12)
     . 006, 005, 013,-007, 005, 019, 021,-040,-014, 013, ! dx2(1:10,18,12)
     . 014, 012, 012, 018, 008, 013, 024,-055,-014, 022, ! dx2(1:10,19,12)
     .-013,-015,-009,-009,-022,-010,-002,-076,-038, 000, ! dx2(1:10,20,12)
     . 002, 000,-002, 005,-009,-007, 006,-073,-030,-003, ! dx2(1:10,21,12)
     .-008,-013,-005, 007,-013,-003, 006,-068,-031,-001, ! dx2(1:10,22,12)
     . 010,-004,-010, 010,-004,-009, 003,-048,-031,-006, ! dx2(1:10,23,12)
     . 000,-001, 000, 021, 005, 005, 003,-020,-019,-003, ! dx2(1:10,24,12)
     . 020, 021, 000, 014, 020, 003, 004,-005,-014, 015, ! dx2(1:10,25,12)
     .-012,-002,-013, 000, 002,-013,-017,-026,-027,-013, ! dx2(1:10,26,12)
     . 025, 005,-010, 011, 014, 006,-012, 000, 000,-003, ! dx2(1:10,27,12)
     . 002,-006,-013,-001, 012, 007,-017, 013, 006,-014, ! dx2(1:10,28,12)
     . 000,-017,-030,-008,-003,-005,-028, 014, 006,-018, ! dx2(1:10,29,12)
     . 004, 000,-012, 005, 008, 001,-009, 042, 030,-004, ! dx2(1:10,30,12)
     .-005,-019,-036,-016,-007,-027,-022, 035, 014,-015, ! dx2(1:10,31,12)
     .-059,-056,-069,-051,-036,-068,-054, 005,-015,-046, ! dx2(1:10,32,12)
     .-022,-029,-048,-026,-026,-043,-023, 023,-009,-019, ! dx2(1:10,33,12)
     .-025,-041,-045,-023,-024,-029,-019, 005,-020,-034, ! dx2(1:10,34,12)
     .-010,-016,-005,-001,-012, 007,-004,-003,-010, 000, ! dx2(1:10,35,12)
     .-003, 001, 023, 010, 006, 042, 020,-020,-006, 002, ! dx2(1:10,36,12)
     . 017, 012, 071, 036, 016, 020, 029,-069,-012, 026, ! dx2(1:10,37,12)
     . 016, 029, 078, 038, 026, 023, 033,-078,-005, 039, ! dx2(1:10,38,12)
     .-019,-006,-018, 015, 003, 033, 001,-023,-042, 005, ! dx2(1:10, 1,13)
     .-005, 000,-024, 001,-003, 013,-004,-031,-043, 010, ! dx2(1:10, 2,13)
     .-026, 006, 002, 003, 009, 020, 012,-015,-002, 019, ! dx2(1:10, 3,13)
     .-015, 000,-008,-005, 000, 009, 000,-025, 000, 003, ! dx2(1:10, 4,13)
     .-025, 005, 002, 008, 008, 011, 000,-006, 018, 006, ! dx2(1:10, 5,13)
     . 004, 019, 016, 034, 027, 023, 021, 026, 042, 014, ! dx2(1:10, 6,13)
     .-001, 014, 032, 035, 031, 013, 010, 046, 040, 012, ! dx2(1:10, 7,13)
     .-045,-027,-014,-010,-017,-024,-026, 011, 000,-027, ! dx2(1:10, 8,13)
     .-021,-010, 003,-004, 009, 010,-002, 050, 018,-004, ! dx2(1:10, 9,13)
     .-019,-014, 000,-017, 000, 000,-012, 063, 008,-010, ! dx2(1:10,10,13)
     .-019,-010,-004,-018, 000, 000,-010, 055, 000,-008, ! dx2(1:10,11,13)
     . 022, 020, 024, 010, 014, 019, 007, 066, 021, 009, ! dx2(1:10,12,13)
     . 017, 038, 026, 029, 020, 019, 021, 050, 015, 026, ! dx2(1:10,13,13)
     .-012,-005, 000,-009,-009,-010,-011, 020,-010,-014, ! dx2(1:10,14,13)
     . 004, 015, 019, 001, 002, 008, 011, 008, 000, 009, ! dx2(1:10,15,13)
     . 001, 016, 010,-016, 000, 000, 016,-014,-011, 005, ! dx2(1:10,16,13)
     . 005, 020, 017,-007, 006, 017, 018,-023,-004, 019, ! dx2(1:10,17,13)
     . 020, 027, 021, 015, 015, 028, 031,-026,-002, 032, ! dx2(1:10,18,13)
     . 019, 005, 014, 027, 013, 021, 023,-049,-011, 015, ! dx2(1:10,19,13)
     .-015,-014,-012,-004,-016,-012,-002,-078,-038,-003, ! dx2(1:10,20,13)
     . 004,-006, 000, 018,-006,-001, 008,-066,-028, 005, ! dx2(1:10,21,13)
     . 002, 004, 000, 009,-005,-005, 000,-062,-027, 000, ! dx2(1:10,22,13)
     . 000, 006,-001, 019, 004, 007, 009,-036,-024,-003, ! dx2(1:10,23,13)
     . 008, 009,-004, 014, 010, 003, 007,-019,-020, 012, ! dx2(1:10,24,13)
     . 018, 019, 007, 026, 034, 022, 003, 006,-001, 000, ! dx2(1:10,25,13)
     . 010, 001,-020, 002, 003,-003,-018,-019,-024,-014, ! dx2(1:10,26,13)
     . 005, 002,-013, 002, 021, 013,-007, 002, 000,-005, ! dx2(1:10,27,13)
     . 012,-004,-013, 001, 014, 015,-016, 018, 014, 000, ! dx2(1:10,28,13)
     .-003,-005,-022,-011, 011, 006,-013, 025, 015,-014, ! dx2(1:10,29,13)
     . 005,-011,-028,-002, 007, 001,-009, 039, 030, 000, ! dx2(1:10,30,13)
     .-008,-003,-028,-014, 000,-020,-013, 042, 022,-010, ! dx2(1:10,31,13)
     .-060,-043,-071,-040,-037,-056,-045, 013,-011,-045, ! dx2(1:10,32,13)
     .-021,-022,-043,-024,-017,-036,-020, 020,-008,-015, ! dx2(1:10,33,13)
     .-026,-024,-026,-021,-023,-036,-018, 013,-016,-017, ! dx2(1:10,34,13)
     .-017,-012,-001,-004,-009, 007, 002,-003,-014,-011, ! dx2(1:10,35,13)
     . 023, 000, 030, 020, 016, 037, 021,-022,-006, 009, ! dx2(1:10,36,13)
     . 020, 039, 091, 034, 034, 030, 038,-070,-007, 032, ! dx2(1:10,37,13)
     . 010, 018, 096, 039, 029, 016, 021,-071,-027,-014, ! dx2(1:10,38,13)
     .-003, 040,-008, 005, 020, 022, 018,-032,-039, 034, ! dx2(1:10, 1,14)
     .-004,-002,-027,-010, 006, 026, 006,-038,-049, 020, ! dx2(1:10, 2,14)
     .-029, 009,-018,-006,-001, 011, 002,-033,-014, 007, ! dx2(1:10, 3,14)
     .-033,-002, 000,-004, 011, 010, 003,-016, 000, 005, ! dx2(1:10, 4,14)
     . 001, 006, 006, 022, 006, 017, 002,-001, 018, 006, ! dx2(1:10, 5,14)
     .-013, 013, 014, 035, 024, 014, 013, 019, 031, 009, ! dx2(1:10, 6,14)
     .-010, 009, 025, 025, 023, 025, 019, 044, 040, 013, ! dx2(1:10, 7,14)
     .-049,-026,-024,-024,-016,-016,-027, 001,-001,-023, ! dx2(1:10, 8,14)
     .-017,-006, 006,-008, 007, 006,-002, 051, 020,-002, ! dx2(1:10, 9,14)
     .-018,-006, 010,-011, 009, 011, 000, 062, 018, 000, ! dx2(1:10,10,14)
     .-007,-006, 008,-006, 000, 006,-005, 057, 010,-004, ! dx2(1:10,11,14)
     .-013, 020, 023, 025, 016, 018, 009, 070, 024, 014, ! dx2(1:10,12,14)
     . 026, 034, 040, 021, 023, 022, 025, 055, 024, 022, ! dx2(1:10,13,14)
     .-017,-004, 009,-002,-013,-007,-004, 026,-005,-011, ! dx2(1:10,14,14)
     . 003, 021, 023,-016, 005, 000, 011, 009,-002, 011, ! dx2(1:10,15,14)
     . 009, 014, 017,-015, 004, 007, 011,-007,-009, 005, ! dx2(1:10,16,14)
     . 009, 012, 023, 008, 006, 017, 016,-022,-002, 018, ! dx2(1:10,17,14)
     . 020, 021, 020, 029, 018, 025, 032,-025, 004, 024, ! dx2(1:10,18,14)
     . 017, 022, 029, 037, 023, 020, 030,-039,-012, 033, ! dx2(1:10,19,14)
     .-021,-012,-006, 010,-007,-009, 000,-064,-034, 007, ! dx2(1:10,20,14)
     . 002, 002, 000, 016, 000, 006, 006,-062,-023, 001, ! dx2(1:10,21,14)
     . 004, 009, 007, 016, 005,-002, 009,-048,-020, 014, ! dx2(1:10,22,14)
     . 007, 005,-003, 014, 010, 003, 007,-034,-023, 002, ! dx2(1:10,23,14)
     . 012, 011, 002, 024, 024, 013, 005,-012,-013, 008, ! dx2(1:10,24,14)
     . 006, 023,-008, 018, 032, 021, 004, 000, 000, 013, ! dx2(1:10,25,14)
     .-003,-008,-022, 005, 018, 000,-011,-016,-018,-003, ! dx2(1:10,26,14)
     . 003, 010,-016, 002, 021, 013,-016, 007, 006,-008, ! dx2(1:10,27,14)
     . 000,-001,-017,-002, 015, 014,-010, 023, 012,-004, ! dx2(1:10,28,14)
     . 011,-003,-029,-012, 008, 007,-014, 026, 023,-005, ! dx2(1:10,29,14)
     .-007,-011,-022,-009, 011, 011,-008, 041, 033,-002, ! dx2(1:10,30,14)
     . 006,-003,-027,-004, 005,-014, 000, 052, 031,-006, ! dx2(1:10,31,14)
     .-019,-051,-063,-043,-035,-049,-043, 010,-010,-051, ! dx2(1:10,32,14)
     .-010,-017,-029,-020,-015,-040,-015, 029,-007,-009, ! dx2(1:10,33,14)
     .-022,-013,-014,-009,-011,-026,-008, 010,-013,-024, ! dx2(1:10,34,14)
     . 003,-002, 007,-001,-005, 005, 003,-005,-018,-005, ! dx2(1:10,35,14)
     . 001, 015, 047, 020, 019, 034, 027,-024,-003, 019, ! dx2(1:10,36,14)
     . 024, 019, 078, 033, 019, 025, 026,-085,-029,-002, ! dx2(1:10,37,14)
     .-001, 016, 091, 027, 031, 013, 035,-095,-022, 033, ! dx2(1:10,38,14)
     . 012, 011,-039,-004, 016, 033, 002,-030,-032, 031, ! dx2(1:10, 1,15)
     .-022, 017,-043,-005, 002, 021, 001,-043,-046, 034, ! dx2(1:10, 2,15)
     .-005, 000,-019,-006, 000, 019, 006,-036,-012, 012, ! dx2(1:10, 3,15)
     .-044,-002,-019, 000,-006, 006,-011,-037,-011,-007, ! dx2(1:10, 4,15)
     .-020,-003, 000, 023, 005, 011, 003,-007, 011, 011, ! dx2(1:10, 5,15)
     .-003, 019, 027, 038, 021, 031, 013, 019, 036, 014, ! dx2(1:10, 6,15)
     .-010, 000, 011, 009, 024, 011, 009, 037, 032, 007, ! dx2(1:10, 7,15)
     .-049,-036,-017,-022,-021,-025,-026,-008,-005,-032, ! dx2(1:10, 8,15)
     .-024,-007, 000,-003, 010, 008,-005, 045, 019, 000, ! dx2(1:10, 9,15)
     .-018, 001, 013, 008, 010, 006, 001, 057, 022, 004, ! dx2(1:10,10,15)
     .-018,-004, 006, 006, 001, 006,-006, 065, 015, 004, ! dx2(1:10,11,15)
     . 009, 019, 033, 027, 016, 022, 024, 075, 030, 013, ! dx2(1:10,12,15)
     . 017, 017, 041, 010, 022, 019, 021, 059, 026, 020, ! dx2(1:10,13,15)
     .-007, 003, 012,-019,-006,-006, 001, 027,-005, 000, ! dx2(1:10,14,15)
     . 006, 018, 031,-013, 007, 001, 012, 017, 001, 005, ! dx2(1:10,15,15)
     . 000, 018, 023, 004, 004, 002, 015,-008,-003, 008, ! dx2(1:10,16,15)
     . 022, 017, 018, 017, 011, 014, 022,-017,-009, 017, ! dx2(1:10,17,15)
     . 019, 029, 029, 040, 019, 022, 029,-021, 004, 034, ! dx2(1:10,18,15)
     . 009, 023, 023, 039, 022, 025, 032,-029,-002, 028, ! dx2(1:10,19,15)
     . 005, 004,-001, 014, 000, 002, 000,-058,-020, 003, ! dx2(1:10,20,15)
     . 001, 009, 005, 004, 002,-008, 007,-052,-025, 009, ! dx2(1:10,21,15)
     . 009, 008, 003, 015, 015, 004, 010,-041,-015, 007, ! dx2(1:10,22,15)
     . 005, 007,-004, 016, 015, 009, 004,-030,-016, 015, ! dx2(1:10,23,15)
     . 009, 007, 000, 021, 022, 013, 008,-011,-006, 011, ! dx2(1:10,24,15)
     . 022, 017,-004, 017, 035, 020, 002, 004, 005, 012, ! dx2(1:10,25,15)
     .-001, 005,-031,-002, 008, 000,-018,-018,-015,-003, ! dx2(1:10,26,15)
     . 008,-006,-021,-005, 017, 019,-015, 010, 008,-004, ! dx2(1:10,27,15)
     . 000,-004,-031,-008, 014, 017,-019, 025, 017, 002, ! dx2(1:10,28,15)
     . 013,-007,-032,-013, 006, 007,-012, 028, 025,-008, ! dx2(1:10,29,15)
     . 006, 012,-014,-010, 016, 019, 003, 055, 044, 005, ! dx2(1:10,30,15)
     .-007,-010,-014,-002, 003,-012, 000, 051, 027, 003, ! dx2(1:10,31,15)
     .-034,-025,-056,-038,-029,-043,-040, 022, 000,-026, ! dx2(1:10,32,15)
     . 000,-017,-019,-014,-017,-038,-014, 028,-008,-019, ! dx2(1:10,33,15)
     .-009,-005,-002,-005,-011,-036, 000, 013,-013,-017, ! dx2(1:10,34,15)
     . 005, 000, 009, 000,-004, 003, 002,-008,-012,-002, ! dx2(1:10,35,15)
     . 003,-003, 046, 017, 012, 028, 015,-037,-018,-002, ! dx2(1:10,36,15)
     . 002, 002, 090, 030, 023, 016, 025,-098,-032, 004, ! dx2(1:10,37,15)
     .-001, 021, 092, 039, 028, 000, 026,-093,-031,-007, ! dx2(1:10,38,15)
     .-005,-003,-066,-026,-008, 035,-014,-066,-062, 013, ! dx2(1:10, 1,16)
     .-025, 007,-067,-048,-004, 013,-013,-071,-060, 020, ! dx2(1:10, 2,16)
     .-014, 013,-031, 006, 004, 022, 011,-036,-017, 016, ! dx2(1:10, 3,16)
     .-012, 002,-018, 008, 001, 013,-006,-046,-012, 007, ! dx2(1:10, 4,16)
     .-034,-016,-017, 009,-012, 000,-008,-034,-005,-005, ! dx2(1:10, 5,16)
     .-023, 015, 012, 016, 023, 019, 013, 005, 031, 006, ! dx2(1:10, 6,16)
     .-021, 006, 024, 017, 020, 018, 012, 030, 025, 006, ! dx2(1:10, 7,16)
     .-071,-040,-032,-031,-023,-028,-041,-017,-016,-042, ! dx2(1:10, 8,16)
     .-017,-002, 014, 000, 007, 009, 001, 041, 023,-001, ! dx2(1:10, 9,16)
     .-025,-005, 009, 012, 006, 003, 000, 054, 020,-009, ! dx2(1:10,10,16)
     .-015, 005, 017, 014, 007, 007, 002, 059, 024, 003, ! dx2(1:10,11,16)
     .-001, 017, 032, 008, 017, 021, 020, 073, 028, 017, ! dx2(1:10,12,16)
     . 007, 046, 055,-001, 022, 018, 030, 061, 032, 023, ! dx2(1:10,13,16)
     .-007, 001, 014,-034,-008,-008,-004, 030, 002,-003, ! dx2(1:10,14,16)
     . 000, 020, 037, 003, 002, 005, 016, 019, 007, 008, ! dx2(1:10,15,16)
     . 007, 014, 030, 014, 007, 002, 014, 005,-001, 010, ! dx2(1:10,16,16)
     .-006, 013, 025, 025, 009, 006, 018,-009, 000, 018, ! dx2(1:10,17,16)
     . 037, 031, 030, 038, 023, 022, 033,-011, 006, 029, ! dx2(1:10,18,16)
     . 013, 024, 024, 033, 025, 019, 027,-025, 000, 030, ! dx2(1:10,19,16)
     . 002,-003, 003, 005, 000,-011, 003,-050,-030, 008, ! dx2(1:10,20,16)
     .-001, 006, 002, 013, 003, 003, 002,-046,-012, 003, ! dx2(1:10,21,16)
     . 004, 016, 005, 015, 013, 006, 011,-034,-015, 011, ! dx2(1:10,22,16)
     . 016, 012,-008, 015, 012, 001,-003,-029,-012,-001, ! dx2(1:10,23,16)
     . 002, 009,-001, 012, 025, 008, 005,-014,-004, 011, ! dx2(1:10,24,16)
     . 017, 016,-010, 015, 025, 031,-001, 011, 009, 009, ! dx2(1:10,25,16)
     .-008,-006,-026,-005, 009, 003,-024,-016,-009,-008, ! dx2(1:10,26,16)
     .-002, 001,-029,-004, 017, 013,-014, 014, 013,-002, ! dx2(1:10,27,16)
     .-016,-003,-030,-015, 013, 019,-014, 032, 026, 000, ! dx2(1:10,28,16)
     . 008,-001,-032,-012, 005, 014,-015, 036, 028,-003, ! dx2(1:10,29,16)
     . 006, 002,-018,-003, 015, 013, 002, 053, 041, 014, ! dx2(1:10,30,16)
     .-004, 004,-010,-007,-001,-009,-007, 054, 029,-005, ! dx2(1:10,31,16)
     .-030,-030,-036,-034,-028,-034,-029, 033,-004,-026, ! dx2(1:10,32,16)
     .-019,-013,-014,-008,-011,-031,-006, 029,-008,-009, ! dx2(1:10,33,16)
     . 014,-009, 000,-014,-011,-036,-008, 006,-024,-016, ! dx2(1:10,34,16)
     .-005,-023, 020, 000,-006,-014,-004,-020,-029,-014, ! dx2(1:10,35,16)
     . 000, 000, 040, 012, 006, 020, 011,-046,-022, 002, ! dx2(1:10,36,16)
     . 014,-003, 097, 017, 011, 011, 010,-108,-043,-034, ! dx2(1:10,37,16)
     .-001, 000, 100, 020, 006, 023, 005,-112,-042,-004, ! dx2(1:10,38,16)
     .-011, 009,-088,-044,-002, 032,-012,-085,-070, 019, ! dx2(1:10, 1,17)
     .-015,-006,-103,-069,-021, 006,-026,-100,-078,-007, ! dx2(1:10, 2,17)
     .-020, 007,-032, 000,-002, 023, 000,-059,-029, 026, ! dx2(1:10, 3,17)
     .-021,-008,-042, 004,-010, 007,-004,-054,-024,-008, ! dx2(1:10, 4,17)
     .-015,-009,-012,-012, 002, 005, 000,-046, 000,-004, ! dx2(1:10, 5,17)
     .-028,-001, 002, 001, 010, 015, 000,-020, 014,-003, ! dx2(1:10, 6,17)
     .-036,-018,-001, 001, 007, 001,-004, 009, 025, 000, ! dx2(1:10, 7,17)
     .-065,-046,-020,-026,-025,-027,-037,-035,-016,-041, ! dx2(1:10, 8,17)
     .-033,-009, 000, 004, 006,-006,-009, 033, 019,-006, ! dx2(1:10, 9,17)
     .-028, 002, 007, 022, 000,-001, 000, 051, 024,-007, ! dx2(1:10,10,17)
     .-024,-002, 017,-006, 004, 000, 003, 058, 022,-005, ! dx2(1:10,11,17)
     .-012, 015, 024,-015, 004, 002, 010, 061, 027, 006, ! dx2(1:10,12,17)
     . 021, 021, 051, 000, 021, 018, 028, 057, 031, 018, ! dx2(1:10,13,17)
     .-025,-007, 020,-025,-017,-019,-008, 025,-004,-017, ! dx2(1:10,14,17)
     .-008, 011, 024, 002,-005,-013, 006, 014, 000, 001, ! dx2(1:10,15,17)
     . 008, 013, 030, 019, 000,-005, 013, 007,-001, 008, ! dx2(1:10,16,17)
     . 019, 030, 036, 028, 015, 007, 024, 001, 009, 027, ! dx2(1:10,17,17)
     . 015, 023, 025, 027, 021, 013, 022,-007, 007, 028, ! dx2(1:10,18,17)
     . 029, 033, 023, 026, 020, 021, 020,-015, 011, 034, ! dx2(1:10,19,17)
     . 001, 006,-003, 002, 000,-009,-006,-042,-013, 011, ! dx2(1:10,20,17)
     . 001, 008, 004, 009, 004,-001, 000,-040,-018, 007, ! dx2(1:10,21,17)
     . 018, 009, 003, 014, 005,-004, 000,-032,-005, 008, ! dx2(1:10,22,17)
     .-017,-011,-002, 001, 006, 004,-012,-026,-012, 006, ! dx2(1:10,23,17)
     .-005, 001,-012, 006, 013, 008,-003,-005, 003, 002, ! dx2(1:10,24,17)
     . 006, 000,-010,-004, 024, 020,-001, 004, 017, 012, ! dx2(1:10,25,17)
     .-027,-022,-039,-019, 000,-004,-029,-014,-011,-016, ! dx2(1:10,26,17)
     . 005,-001,-037,-013, 006, 021,-016, 019, 022, 002, ! dx2(1:10,27,17)
     . 005, 000,-030,-014, 011, 016,-018, 032, 030,-003, ! dx2(1:10,28,17)
     .-005,-001,-024,-010, 003, 019,-013, 043, 038, 002, ! dx2(1:10,29,17)
     . 015, 006,-016,-002, 017, 017, 003, 059, 045, 012, ! dx2(1:10,30,17)
     . 011, 008, 000,-003, 000, 004, 005, 055, 025, 003, ! dx2(1:10,31,17)
     .-022,-025,-042,-031,-031,-035,-030, 025,-009,-033, ! dx2(1:10,32,17)
     . 000, 002, 000,-016,-018,-038,-010, 022,-014,-015, ! dx2(1:10,33,17)
     .-017,-013, 003,-010,-023,-047,-013,-001,-029,-025, ! dx2(1:10,34,17)
     . 007, 000, 025, 000,-007,-016,-005,-027,-027,-014, ! dx2(1:10,35,17)
     .-023,-011, 050,-002,-004, 005,-003,-053,-037,-026, ! dx2(1:10,36,17)
     .-025,-003, 095, 009, 002, 001, 007,-117,-055,-021, ! dx2(1:10,37,17)
     .-001, 003, 093, 009,-008,-018,-001,-123,-059,-015, ! dx2(1:10,38,17)
     .-012,-020,-131,-098,-038, 007,-032,-136,-109,-007, ! dx2(1:10, 1,18)
     .-047,-042,-138,-090,-053,-017,-044,-138,-114, 000, ! dx2(1:10, 2,18)
     . 000, 001,-089,-022,-022, 005,-018,-098,-059,-006, ! dx2(1:10, 3,18)
     .-039,-025,-066,-051,-026,-010,-019,-098,-039,-015, ! dx2(1:10, 4,18)
     .-059,-012,-045,-042,-024,-012,-026,-081,-025,-016, ! dx2(1:10, 5,18)
     .-025, 000,-009,-003, 012, 006, 011,-040, 022, 007, ! dx2(1:10, 6,18)
     .-026, 003, 020, 013, 016, 000, 000,-010, 024, 005, ! dx2(1:10, 7,18)
     .-065,-065,-053,-042,-039,-043,-053,-061,-029,-032, ! dx2(1:10, 8,18)
     .-052,-012,-007, 002,-013,-021,-017, 007, 011,-031, ! dx2(1:10, 9,18)
     .-015,-013, 009,-015,-006,-008,-010, 037, 014,-007, ! dx2(1:10,10,18)
     .-024,-003, 014,-034,-013,-006,-002, 042, 017,-009, ! dx2(1:10,11,18)
     .-003, 004, 026,-008, 000,-003, 004, 045, 015,-007, ! dx2(1:10,12,18)
     .-008, 028, 030, 002,-001,-004, 012, 041, 023, 006, ! dx2(1:10,13,18)
     .-047,-016, 004,-028,-034,-037,-025, 008,-011,-023, ! dx2(1:10,14,18)
     .-023,-018, 017,-003,-017,-035,-007, 015,-008,-016, ! dx2(1:10,15,18)
     .-035,-006, 009,-009,-021,-031,-012, 005,-009,-004, ! dx2(1:10,16,18)
     .-009, 001, 015, 004,-009,-021, 000,-004,-009, 001, ! dx2(1:10,17,18)
     . 035, 014, 026, 019, 005, 000, 011, 010, 006, 015, ! dx2(1:10,18,18)
     . 003, 020, 018, 018, 011, 004, 009,-005, 004, 012, ! dx2(1:10,19,18)
     .-014,-015,-005, 001,-018,-010,-010,-027,-017, 002, ! dx2(1:10,20,18)
     .-009,-011,-011,-002,-010,-016,-020,-036,-015,-013, ! dx2(1:10,21,18)
     .-023,-016, 000,-005,-004, 002,-017,-029,-009, 009, ! dx2(1:10,22,18)
     .-012, 000,-020,-013,-008,-018,-022,-025, 000,-001, ! dx2(1:10,23,18)
     .-021,-013,-017,-021,-005,-003,-024,-020, 003,-002, ! dx2(1:10,24,18)
     .-022, 004,-026,-012,-003, 004,-019, 011, 019,-009, ! dx2(1:10,25,18)
     .-050,-017,-053,-040,-027,-015,-043,-013, 000,-027, ! dx2(1:10,26,18)
     .-014,-013,-025,-031,-009,-002,-025, 017, 021,-018, ! dx2(1:10,27,18)
     .-007, 000,-037,-014,-009, 008,-021, 029, 037,-008, ! dx2(1:10,28,18)
     . 021, 001,-021,-022,-007, 000,-025, 039, 030, 002, ! dx2(1:10,29,18)
     . 019, 005,-004,-012, 002, 005,-010, 052, 049,-009, ! dx2(1:10,30,18)
     .-006, 014,-002,-004,-008,-026,-009, 047, 005, 008, ! dx2(1:10,31,18)
     .-046,-009,-042,-051,-051,-064,-043, 006,-026,-043, ! dx2(1:10,32,18)
     .-021,-019, 002,-029,-018,-042,-016, 011,-024,-024, ! dx2(1:10,33,18)
     .-027,-017, 012,-028,-018,-052,-018,-019,-032,-034, ! dx2(1:10,34,18)
     .-037,-033, 026,-027,-026,-046,-025,-040,-052,-051, ! dx2(1:10,35,18)
     .-017,-012, 055,-018,-023,-021,-017,-091,-063,-035, ! dx2(1:10,36,18)
     .-026,-053, 092,-009,-026,-014,-028,-138,-074,-058, ! dx2(1:10,37,18)
     .-127,-055, 081,-014,-011,-012,-024,-158,-074,-079/ ! dx2(1:10,38,18)

      data dy1 / 
     .-005,-001,-098,-089,-044, 072,-002,-064,-066, 031, ! dy1(1:10, 1, 1)
     .-021,-004,-106,-078,-064, 049, 000,-093,-050,-019, ! dy1(1:10, 2, 1)
     .-008,-031,-083,-088,-073, 021,-024,-074,-048,-012, ! dy1(1:10, 3, 1)
     .-011,-010,-065,-053,-061, 023,-015,-060,-022,-029, ! dy1(1:10, 4, 1)
     . 002,-021,-052,-044,-054, 025,-010,-037,-001,-013, ! dy1(1:10, 5, 1)
     .-020,-005,-039,-036,-044, 006,-020,-020, 000,-021, ! dy1(1:10, 6, 1)
     . 024,-005,-017, 001,-033, 026, 005, 009, 030,-005, ! dy1(1:10, 7, 1)
     .-005,-028,-037,-025,-047, 009,-016,-011, 003,-031, ! dy1(1:10, 8, 1)
     .-033,-024,-030,-032,-057,-007,-034,-013,-007,-028, ! dy1(1:10, 9, 1)
     . 001,-013,-015,-012,-034, 000,-013, 006, 011,-027, ! dy1(1:10,10, 1)
     .-010,-006,-006,-005,-028,-004,-010, 010, 008,-017, ! dy1(1:10,11, 1)
     . 000,-021, 002,-019,-034,-007,-018, 000, 000,-025, ! dy1(1:10,12, 1)
     . 001,-016, 006,-006,-017, 007,-005, 018, 005,-004, ! dy1(1:10,13, 1)
     .-027,-023, 000,-019,-036,-014,-024,-001,-024,-039, ! dy1(1:10,14, 1)
     .-036,-036,-016,-021,-039,-024,-031,-018,-031,-040, ! dy1(1:10,15, 1)
     .-015,-014,-007,-022,-025,-014,-022,-012,-039,-041, ! dy1(1:10,16, 1)
     .-019,-039,-010,-014,-018,-014,-027,-021,-035,-019, ! dy1(1:10,17, 1)
     .-032,-027,-015,-023,-018,-022,-028,-034,-048,-033, ! dy1(1:10,18, 1)
     .-054,-042,-012,-023,-034,-038,-035,-044,-057,-037, ! dy1(1:10,19, 1)
     .-044,-039,-025,-024,-037,-039,-040,-046,-056,-049, ! dy1(1:10,20, 1)
     .-043,-038,-024,-040,-020,-034,-038,-044,-052,-046, ! dy1(1:10,21, 1)
     .-022,-040,-020,-017,-023,-029,-017,-037,-035,-028, ! dy1(1:10,22, 1)
     .-041,-028,-025,-030,-016,-032,-017,-043,-033,-020, ! dy1(1:10,23, 1)
     .-023,-027,-007,-010, 000,-020,-010,-030,-014,-016, ! dy1(1:10,24, 1)
     .-002, 005, 020, 003, 003,-005, 001,-013,-001,-003, ! dy1(1:10,25, 1)
     .-030,-015,-009,-015,-019,-015,-020,-033,-014,-006, ! dy1(1:10,26, 1)
     .-001,-017, 017,-013,-003,-016,-005,-017,-010,-017, ! dy1(1:10,27, 1)
     .-006, 016, 020, 000, 011, 000, 016,-003, 012, 002, ! dy1(1:10,28, 1)
     . 010, 030, 024, 013, 007,-003, 018,-010, 016, 015, ! dy1(1:10,29, 1)
     . 029, 021, 018, 008, 015, 008, 021,-009, 011, 014, ! dy1(1:10,30, 1)
     . 024, 055, 035, 020, 007, 023, 036,-006, 032, 044, ! dy1(1:10,31, 1)
     . 020, 042,-002,-002, 002, 008, 014,-031, 005, 014, ! dy1(1:10,32, 1)
     . 008, 007, 002,-022,-018,-014, 004,-049, 010, 012, ! dy1(1:10,33, 1)
     . 009, 059, 004,-006,-007, 007, 019,-052, 009, 021, ! dy1(1:10,34, 1)
     . 022, 020,-008,-015,-009, 003, 026,-065, 013, 026, ! dy1(1:10,35, 1)
     . 048, 043,-013,-029, 001, 023, 034,-061, 034, 036, ! dy1(1:10,36, 1)
     . 045, 059,-011,-024,-002, 014, 057,-064, 041, 040, ! dy1(1:10,37, 1)
     . 040, 037,-030,-073,-018, 002, 027,-104, 013, 008, ! dy1(1:10,38, 1)
     . 029, 043,-032,-014,-013, 055, 030,-034,-008, 005, ! dy1(1:10, 1, 2)
     . 035, 013,-043,-039,-009, 053, 024,-038,-019, 011, ! dy1(1:10, 2, 2)
     . 023, 023,-037,-015,-017, 030, 010,-028, 001, 010, ! dy1(1:10, 3, 2)
     . 023, 006,-028,-013,-017, 032, 011,-019, 007, 010, ! dy1(1:10, 4, 2)
     . 026, 025,-007, 000,-007, 026, 008,-018, 015, 012, ! dy1(1:10, 5, 2)
     . 034, 023, 002, 006, 000, 035, 017,-005, 036, 022, ! dy1(1:10, 6, 2)
     . 035, 041, 031, 014, 018, 035, 018, 012, 036, 026, ! dy1(1:10, 7, 2)
     . 010, 018,-006, 001,-016, 010,-012,-014, 009, 000, ! dy1(1:10, 8, 2)
     . 012, 005, 004,-008,-010,-002,-006,-003, 005,-014, ! dy1(1:10, 9, 2)
     . 006, 012, 006, 008,-008, 007,-008, 016, 002,-002, ! dy1(1:10,10, 2)
     . 003, 006, 017, 005, 002, 012, 000, 024, 016, 000, ! dy1(1:10,11, 2)
     . 010, 014, 020, 019, 007, 017, 001, 032, 017, 002, ! dy1(1:10,12, 2)
     . 004, 023, 024, 018, 012, 016, 007, 033, 014, 004, ! dy1(1:10,13, 2)
     .-021,-019,-008,-018,-017,-014,-025, 005,-007,-014, ! dy1(1:10,14, 2)
     .-002,-004, 001,-006,-009, 000,-013, 006,-019,-012, ! dy1(1:10,15, 2)
     .-021,-017,-002,-012,-007,-006,-012, 001,-020,-001, ! dy1(1:10,16, 2)
     . 002,-001, 000,-001,-004,-005,-012, 005,-025,-017, ! dy1(1:10,17, 2)
     .-018,-011,-002,-006,-003,-009,-015, 003,-021,-012, ! dy1(1:10,18, 2)
     .-010,-013,-013,-006, 003,-011,-008, 008,-022,-020, ! dy1(1:10,19, 2)
     .-025,-018,-010,-016, 002,-019,-015, 002,-028,-013, ! dy1(1:10,20, 2)
     .-022,-027,-017,-020,-010,-026,-017, 003,-023,-031, ! dy1(1:10,21, 2)
     .-036,-029,-013,-019,-007,-019,-024, 002,-023,-017, ! dy1(1:10,22, 2)
     .-015,-013,-006,-006, 000,-015,-011, 005,-007,-019, ! dy1(1:10,23, 2)
     .-021,-015,-004,-007, 000,-013,-013,-002,-001,-016, ! dy1(1:10,24, 2)
     .-005,-008, 017, 009, 016,-002, 006, 006, 010, 007, ! dy1(1:10,25, 2)
     .-015,-008,-007,-011,-010,-023,-013,-013,-005,-026, ! dy1(1:10,26, 2)
     .-001, 009, 018, 013, 015,-001, 003,-014, 003, 003, ! dy1(1:10,27, 2)
     . 013, 004, 031, 019, 009, 003, 005,-013, 010, 015, ! dy1(1:10,28, 2)
     . 017, 019, 037, 022, 026, 014, 022,-006, 017, 017, ! dy1(1:10,29, 2)
     . 026, 037, 038, 025, 021, 019, 026,-008, 017, 029, ! dy1(1:10,30, 2)
     . 034, 032, 031, 038, 031, 022, 027,-003, 021, 038, ! dy1(1:10,31, 2)
     . 016, 020, 029, 001, 008, 012, 016,-014, 008, 023, ! dy1(1:10,32, 2)
     . 017, 024, 009, 012, 010, 009, 013,-021,-001, 010, ! dy1(1:10,33, 2)
     . 021, 002, 005, 000, 010, 013, 016,-021, 009, 010, ! dy1(1:10,34, 2)
     . 028, 036, 004, 000, 010, 012, 022,-017, 005, 018, ! dy1(1:10,35, 2)
     . 023, 038, 000,-009, 003, 020, 026,-020, 007, 026, ! dy1(1:10,36, 2)
     . 021, 017,-023,-014,-002, 034, 031,-030, 001, 021, ! dy1(1:10,37, 2)
     .-022,-001,-047,-030,-026, 010,-002,-056,-021, 026, ! dy1(1:10,38, 2)
     . 029, 017,-033,-012,-004, 059, 016,-006, 007, 064, ! dy1(1:10, 1, 3)
     . 005, 029,-036,-027,-008, 042, 016,-013,-005, 010, ! dy1(1:10, 2, 3)
     . 018, 003,-032,-022,-013, 015, 002,-029,-004, 013, ! dy1(1:10, 3, 3)
     . 015, 011,-010,-005, 000, 005, 005,-017, 011, 002, ! dy1(1:10, 4, 3)
     . 034, 016, 000, 005, 005, 022, 014,-005, 023, 017, ! dy1(1:10, 5, 3)
     . 029, 023, 015, 018, 014, 011, 014,-012, 031, 020, ! dy1(1:10, 6, 3)
     . 020, 029, 028, 026, 018, 022, 013, 006, 026, 017, ! dy1(1:10, 7, 3)
     . 029, 010, 004, 000, 000,-004,-008,-016, 009,-009, ! dy1(1:10, 8, 3)
     . 006, 001, 000, 009,-003,-004,-017,-016,-001, 000, ! dy1(1:10, 9, 3)
     . 013, 001, 012, 010, 004, 000,-006, 000, 008, 000, ! dy1(1:10,10, 3)
     . 007, 016, 019, 022, 016, 011, 001, 011, 011, 004, ! dy1(1:10,11, 3)
     . 021, 015, 018, 019, 015, 011, 008, 014, 016, 008, ! dy1(1:10,12, 3)
     . 005, 005, 017, 016, 014, 013, 008, 022, 012, 000, ! dy1(1:10,13, 3)
     .-002,-002,-001,-001,-005,-002,-011, 005,-007,-011, ! dy1(1:10,14, 3)
     .-011,-015,-007,-007,-005,-011,-017,-001,-018,-022, ! dy1(1:10,15, 3)
     .-003,-004, 000, 002, 000,-001,-009, 007,-018,-012, ! dy1(1:10,16, 3)
     .-009,-004,-002, 002, 006, 000, 000, 009,-011,-008, ! dy1(1:10,17, 3)
     .-003, 000,-006,-002, 007, 000,-005, 021,-014, 000, ! dy1(1:10,18, 3)
     .-021,-025,-013,-014,-004,-017,-013, 024,-019,-021, ! dy1(1:10,19, 3)
     .-010,-012,-015,-020,-007,-015,-011, 024,-016,-013, ! dy1(1:10,20, 3)
     .-036,-017,-010,-009,-001,-011,-011, 024,-016,-013, ! dy1(1:10,21, 3)
     .-016,-009,-009,-014,-004,-018,-011, 028,-010,-022, ! dy1(1:10,22, 3)
     .-020,-014, 002,-016,-003,-019,-004, 024, 000,-009, ! dy1(1:10,23, 3)
     .-002, 006, 000, 001, 012,-008,-004, 018, 012,-005, ! dy1(1:10,24, 3)
     .-009,-009, 016, 006, 013, 004, 005, 000, 014, 003, ! dy1(1:10,25, 3)
     .-013,-023, 000,-005,-003,-015,-014,-013, 000,-010, ! dy1(1:10,26, 3)
     .-004, 004, 012, 004, 000,-003,-004,-021, 005, 002, ! dy1(1:10,27, 3)
     . 022, 014, 025, 022, 020, 004, 011,-010, 017, 001, ! dy1(1:10,28, 3)
     . 025, 018, 033, 027, 024, 011, 015,-008, 020, 025, ! dy1(1:10,29, 3)
     . 017, 020, 031, 028, 017, 018, 020,-010, 019, 022, ! dy1(1:10,30, 3)
     . 021, 028, 040, 023, 027, 023, 027, 004, 020, 024, ! dy1(1:10,31, 3)
     . 014, 006, 008, 010, 004, 010, 012,-010, 010, 011, ! dy1(1:10,32, 3)
     . 014, 034, 013, 007, 005, 005, 011,-005, 003, 008, ! dy1(1:10,33, 3)
     . 019, 032, 024, 011, 012, 010, 020, 010, 003, 014, ! dy1(1:10,34, 3)
     . 007, 014,-001,-005,-003, 015, 010, 007,-001, 012, ! dy1(1:10,35, 3)
     . 023, 012,-004,-004, 004, 014, 015, 013,-012, 019, ! dy1(1:10,36, 3)
     . 005, 008,-022,-029,-021, 011, 013,-019,-025, 017, ! dy1(1:10,37, 3)
     . 012,-008,-054,-039,-031,-010,-011,-036,-049,-015, ! dy1(1:10,38, 3)
     . 034, 014,-027,-036,-003, 006, 011,-014, 001,-005, ! dy1(1:10, 1, 4)
     . 009,-004,-038,-038,-007, 007, 000,-027, 004, 011, ! dy1(1:10, 2, 4)
     . 018, 007,-012,-017, 002, 008, 006,-025, 001, 002, ! dy1(1:10, 3, 4)
     . 024, 011, 000,-009,-001, 008, 001,-026, 014, 017, ! dy1(1:10, 4, 4)
     . 028, 031, 017, 017, 019, 012, 014,-011, 025, 022, ! dy1(1:10, 5, 4)
     . 029, 014, 031, 014, 021, 011, 009,-016, 023, 012, ! dy1(1:10, 6, 4)
     . 040, 029, 044, 039, 034, 015, 019,-001, 037, 025, ! dy1(1:10, 7, 4)
     .-011,-004, 009, 007,-004,-013,-005,-034, 000,-009, ! dy1(1:10, 8, 4)
     . 011, 003, 012, 014, 005,-007,-015,-016,-003, 000, ! dy1(1:10, 9, 4)
     . 011,-002, 008, 009, 002,-009,-011,-015, 000,-007, ! dy1(1:10,10, 4)
     . 016, 014, 016, 015, 015, 000,-002, 000, 008, 000, ! dy1(1:10,11, 4)
     . 006, 005, 010, 006, 012, 001, 000,-001,-002, 001, ! dy1(1:10,12, 4)
     . 014, 016, 022, 022, 020, 013, 007, 017, 014, 006, ! dy1(1:10,13, 4)
     . 000, 000,-005,-005, 003,-004,-011,-007,-016,-014, ! dy1(1:10,14, 4)
     .-009,-008, 000,-005, 000,-005,-005,-003,-016,-013, ! dy1(1:10,15, 4)
     .-006,-004,-009,-005, 004,-003, 000, 003,-015,-008, ! dy1(1:10,16, 4)
     . 003, 000,-002, 000, 005,-001,-004, 008,-013,-002, ! dy1(1:10,17, 4)
     . 006,-005,-013,-010, 000, 000,-002, 019,-015,-010, ! dy1(1:10,18, 4)
     . 000,-003,-011,-011, 005,-007,-002, 037,-007,-009, ! dy1(1:10,19, 4)
     .-017,-016,-017,-009,-001,-008,-005, 030,-015,-022, ! dy1(1:10,20, 4)
     . 001,-009,-008,-012,-004,-011,-009, 031, 000,-008, ! dy1(1:10,21, 4)
     . 000,-006,-001,-013, 004,-006,-005, 037, 002,-003, ! dy1(1:10,22, 4)
     .-004,-010,-003,-008, 003,-003,-002, 028, 006,-009, ! dy1(1:10,23, 4)
     . 000, 005, 005,-011, 003,-008, 001, 020, 011,-004, ! dy1(1:10,24, 4)
     .-012, 001, 013, 000, 005,-002, 009, 003, 025, 012, ! dy1(1:10,25, 4)
     .-031,-009,-006,-004,-008,-016,-006,-012, 007, 000, ! dy1(1:10,26, 4)
     . 008,-005, 009, 008, 002,-007,-003,-012, 010,-003, ! dy1(1:10,27, 4)
     . 004, 010, 021, 017, 010, 007, 009,-009, 022, 013, ! dy1(1:10,28, 4)
     . 018, 016, 033, 032, 020, 023, 021, 002, 029, 022, ! dy1(1:10,29, 4)
     . 014, 025, 025, 025, 014, 012, 015, 005, 028, 015, ! dy1(1:10,30, 4)
     . 024, 011, 045, 041, 025, 022, 030, 022, 030, 029, ! dy1(1:10,31, 4)
     .-004, 013, 020, 014, 016, 010, 008, 012, 013, 006, ! dy1(1:10,32, 4)
     . 011, 001, 017, 012, 003, 002, 007, 023, 007, 004, ! dy1(1:10,33, 4)
     . 001, 016, 006, 007, 006, 010, 008, 034,-003, 014, ! dy1(1:10,34, 4)
     . 012, 001, 000, 001, 001, 004, 010, 033,-014, 002, ! dy1(1:10,35, 4)
     . 011,-001,-014,-016,-008, 008, 003, 025,-026, 007, ! dy1(1:10,36, 4)
     .-005, 003,-024,-022,-012, 010, 000, 018,-039, 007, ! dy1(1:10,37, 4)
     .-080,-031,-052,-054,-047,-021,-029, 000,-066, 005, ! dy1(1:10,38, 4)
     . 034, 007,-007,-014, 011, 018, 023,-033, 002,-021, ! dy1(1:10, 1, 5)
     . 002, 015,-027,-024,-004, 000,-001,-034,-006,-010, ! dy1(1:10, 2, 5)
     . 018, 009,-005,-001, 001,-003,-001,-025, 004, 001, ! dy1(1:10, 3, 5)
     . 016, 017, 010, 006, 010,-003, 005,-020, 016, 014, ! dy1(1:10, 4, 5)
     . 041, 021, 031, 023, 026, 005, 009,-016, 026, 007, ! dy1(1:10, 5, 5)
     . 025, 024, 027, 023, 018, 004, 007,-023, 014, 010, ! dy1(1:10, 6, 5)
     . 021, 024, 039, 030, 027, 007, 008,-010, 021, 020, ! dy1(1:10, 7, 5)
     . 013,-002, 014, 008, 000,-024,-020,-037,-005,-004, ! dy1(1:10, 8, 5)
     .-008,-003, 000,-002, 000,-019,-023,-032,-013,-016, ! dy1(1:10, 9, 5)
     . 005, 000, 006, 007, 004,-011,-011,-019,-008,-001, ! dy1(1:10,10, 5)
     .-002, 005, 010, 005, 014,-008, 000,-010, 000, 000, ! dy1(1:10,11, 5)
     .-003,-009,-005, 003, 008,-009,-006,-007,-007,-008, ! dy1(1:10,12, 5)
     . 011, 005, 006, 000, 012, 004, 010,-004,-003, 006, ! dy1(1:10,13, 5)
     . 000,-012,-011,-012, 002,-012,-014,-022,-019,-016, ! dy1(1:10,14, 5)
     .-006,-011,-016,-015,-007,-014,-010,-028,-023,-020, ! dy1(1:10,15, 5)
     .-006,-002,-010,-005, 003,-001, 000,-005,-019,-005, ! dy1(1:10,16, 5)
     . 004,-009,-006,-005, 002, 002, 002, 009,-010,-013, ! dy1(1:10,17, 5)
     .-007,-007,-023,-017,-006,-009,-009, 010,-025,-019, ! dy1(1:10,18, 5)
     . 005,-003,-004,-008,-004, 000, 000, 036,-005, 006, ! dy1(1:10,19, 5)
     . 001, 000,-008,-015,-008, 000,-002, 032,-009,-009, ! dy1(1:10,20, 5)
     .-011,-006,-014, 000, 000,-005,-005, 035,-002,-009, ! dy1(1:10,21, 5)
     . 000, 001,-009,-006,-006,-002,-005, 029,-002,-003, ! dy1(1:10,22, 5)
     . 003, 004, 000,-010,-003,-004, 001, 025, 002,-005, ! dy1(1:10,23, 5)
     .-007,-009,-004,-006,-004,-005,-001, 016, 005, 000, ! dy1(1:10,24, 5)
     . 008, 001, 005, 007, 010,-001, 005, 006, 015, 011, ! dy1(1:10,25, 5)
     .-001,-008,-007,-010,-014,-009,-010,-005, 003,-004, ! dy1(1:10,26, 5)
     .-006,-011, 003,-002,-004,-005,-003,-014, 008, 000, ! dy1(1:10,27, 5)
     . 006, 010, 015, 013, 007, 004, 009, 003, 030, 016, ! dy1(1:10,28, 5)
     . 013, 016, 025, 018, 010, 008, 014, 012, 033, 015, ! dy1(1:10,29, 5)
     . 003,-009, 023, 004, 003, 010, 001, 005, 023, 003, ! dy1(1:10,30, 5)
     .-002, 024, 039, 033, 023, 025, 022, 042, 034, 023, ! dy1(1:10,31, 5)
     . 001, 003, 022, 004,-005, 001, 006, 015, 019, 015, ! dy1(1:10,32, 5)
     .-021,-018, 000, 003,-006,-002,-004, 031, 003,-008, ! dy1(1:10,33, 5)
     . 010,-009, 007, 008,-005,-001, 006, 041,-011, 001, ! dy1(1:10,34, 5)
     .-012,-012, 000,-008,-003, 005, 001, 043,-024, 001, ! dy1(1:10,35, 5)
     .-015,-018,-013,-019,-011,-011,-006, 036,-039, 000, ! dy1(1:10,36, 5)
     .-013, 000,-031,-039,-021, 015, 003, 037,-060, 000, ! dy1(1:10,37, 5)
     .-025,-020,-052,-076,-047,-021,-027, 000,-092,-039, ! dy1(1:10,38, 5)
     . 021, 003,-011,-011, 008,-002,-006,-035,-002, 029, ! dy1(1:10, 1, 6)
     . 018,-002,-009,-010,-006,-010,-003,-037,-006, 010, ! dy1(1:10, 2, 6)
     . 007, 007, 002,-010, 005, 000,-003,-023, 000, 001, ! dy1(1:10, 3, 6)
     . 023, 009, 013, 009, 014,-004, 001,-019, 008, 003, ! dy1(1:10, 4, 6)
     . 018, 018, 032, 015, 019, 007, 007,-011, 018, 006, ! dy1(1:10, 5, 6)
     . 028, 025, 039, 021, 030, 008, 010,-010, 021, 014, ! dy1(1:10, 6, 6)
     . 015, 017, 032, 020, 022, 003, 003,-015, 009, 010, ! dy1(1:10, 7, 6)
     .-012,-007, 006,-006, 001,-014,-020,-036,-012,-006, ! dy1(1:10, 8, 6)
     .-014,-019, 000,-010,-014,-023,-025,-035,-023,-016, ! dy1(1:10, 9, 6)
     .-014,-015, 002,-008,-002,-017,-022,-024,-017,-010, ! dy1(1:10,10, 6)
     .-003,-008, 001,-002, 001,-004,-008,-016,-010, 000, ! dy1(1:10,11, 6)
     .-009,-003,-006,-006, 003,-009,-006,-026,-008, 000, ! dy1(1:10,12, 6)
     . 007, 000, 001,-002, 014, 002, 007,-017,-004,-010, ! dy1(1:10,13, 6)
     .-010,-028,-020,-018,-006,-015,-011,-037,-021,-015, ! dy1(1:10,14, 6)
     .-008,-010,-028,-025,-013,-015,-013,-033,-025,-022, ! dy1(1:10,15, 6)
     .-006,-010,-017,-015,-004,-009,-006,-018,-017,-022, ! dy1(1:10,16, 6)
     .-008,-002,-017,-020,-008, 000,-007,-007,-019,-004, ! dy1(1:10,17, 6)
     . 004, 006,-013,-008, 000, 006, 003, 008,-008, 003, ! dy1(1:10,18, 6)
     .-015,-009,-023,-018,-010,-005,-011, 007,-013,-012, ! dy1(1:10,19, 6)
     .-003,-011,-027,-019,-015,-014,-019, 003,-021,-019, ! dy1(1:10,20, 6)
     . 004,-003,-021,-025,-019,-013,-016, 019,-018,-010, ! dy1(1:10,21, 6)
     .-004, 000,-012,-008,-009,-007,-005, 018,-005,-006, ! dy1(1:10,22, 6)
     . 000, 003,-013,-006,-011,-001,-001, 014,-004, 004, ! dy1(1:10,23, 6)
     . 007, 002, 001, 000,-005,-002, 000, 011, 007, 005, ! dy1(1:10,24, 6)
     . 016, 010, 005, 022,-002, 007, 015, 011, 009, 005, ! dy1(1:10,25, 6)
     . 003,-011,-006,-003,-014,-011,-004,-011, 002,-001, ! dy1(1:10,26, 6)
     . 000,-008,-010, 004,-011,-011, 000,-002, 009, 002, ! dy1(1:10,27, 6)
     .-001,-007, 003, 011,-006,-004, 003, 002, 015, 002, ! dy1(1:10,28, 6)
     .-010,-011, 019, 009,-001, 007, 009, 009, 028, 005, ! dy1(1:10,29, 6)
     . 022, 015, 030, 020, 010, 012, 016, 017, 032, 029, ! dy1(1:10,30, 6)
     . 000, 002, 027, 008, 009, 013, 014, 040, 036, 010, ! dy1(1:10,31, 6)
     .-027,-011, 015, 002, 000,-007, 001, 024, 017,-001, ! dy1(1:10,32, 6)
     .-006,-014, 007, 001,-012,-008,-005, 031, 000,-002, ! dy1(1:10,33, 6)
     .-028,-014, 004, 003,-004, 000, 000, 047,-011,-003, ! dy1(1:10,34, 6)
     . 000,-006, 004, 005,-010, 002, 001, 045,-025, 008, ! dy1(1:10,35, 6)
     .-016,-017,-017,-018,-018,-009,-009, 041,-049,-008, ! dy1(1:10,36, 6)
     .-042,-010,-023,-039,-032,-013,-012, 028,-077,-009, ! dy1(1:10,37, 6)
     .-077,-041,-050,-053,-053,-039,-033, 014,-094,-050, ! dy1(1:10,38, 6)
     . 023, 005,-002,-011,-002, 019, 008,-036, 002,-013, ! dy1(1:10, 1, 7)
     .-006,-004,-015,-018, 000,-003,-009,-043,-012,-006, ! dy1(1:10, 2, 7)
     . 019,-002, 005, 002, 006,-004, 003,-025,-002, 007, ! dy1(1:10, 3, 7)
     . 003, 009, 013, 002, 010,-001, 003,-029, 000, 000, ! dy1(1:10, 4, 7)
     . 029, 013, 025, 014, 018, 008, 004,-015, 005, 014, ! dy1(1:10, 5, 7)
     . 007, 003, 024, 018, 015, 001, 001,-015, 002, 009, ! dy1(1:10, 6, 7)
     . 007, 008, 034, 010, 018, 004, 003,-008, 004, 006, ! dy1(1:10, 7, 7)
     .-008,-014, 003,-007,-008,-020,-024,-033,-023,-009, ! dy1(1:10, 8, 7)
     .-019,-032,-007,-017,-010,-024,-026,-025,-024,-017, ! dy1(1:10, 9, 7)
     .-026,-015,-007,-011,-005,-012,-015,-015,-013,-010, ! dy1(1:10,10, 7)
     .-015,-015,-011,-014,-005,-011,-015,-024,-011,-018, ! dy1(1:10,11, 7)
     . 002,-008,-008,-008, 003,-001,-008,-021,-011,-013, ! dy1(1:10,12, 7)
     . 000, 006,-004,-002, 004, 006, 002,-017, 003, 000, ! dy1(1:10,13, 7)
     .-019,-014,-034,-025,-017,-016,-020,-037,-022,-017, ! dy1(1:10,14, 7)
     .-008,-016,-026,-023,-015,-004,-015,-030,-012,-019, ! dy1(1:10,15, 7)
     .-002, 000,-022,-013,-010, 001,-006,-025,-010,-010, ! dy1(1:10,16, 7)
     . 011, 001,-019,-012,-010, 002,-010,-021,-007,-009, ! dy1(1:10,17, 7)
     .-002,-005,-021,-018,-012,-007,-018,-022,-016,-006, ! dy1(1:10,18, 7)
     . 022, 016,-008,-006, 000, 004,-001,-003,-007,-006, ! dy1(1:10,19, 7)
     .-006, 006,-014,-008,-009, 000,-004,-003,-009, 000, ! dy1(1:10,20, 7)
     . 000, 000,-018,-018,-010,-011,-011,-010,-016,-004, ! dy1(1:10,21, 7)
     .-010, 006,-010,-007,-005,-006,-013, 000,-010, 006, ! dy1(1:10,22, 7)
     . 009, 003,-009,-009,-009,-002,-007, 002,-003, 004, ! dy1(1:10,23, 7)
     . 006, 005,-010,-007,-012,-010,-007,-002,-011, 000, ! dy1(1:10,24, 7)
     . 001, 000, 003, 005,-005, 008, 012, 008, 011, 016, ! dy1(1:10,25, 7)
     .-012,-012,-016,-007,-016,-013,-005,-008,-011,-008, ! dy1(1:10,26, 7)
     . 002,-014, 002, 002,-012,-006, 000, 000, 002, 009, ! dy1(1:10,27, 7)
     .-009,-009,-005,-004,-015,-013, 003, 004, 010,-005, ! dy1(1:10,28, 7)
     . 000, 007, 014, 016, 000, 001, 011, 017, 016, 014, ! dy1(1:10,29, 7)
     . 002,-006, 016, 009,-001, 004, 012, 024, 020, 011, ! dy1(1:10,30, 7)
     . 007, 000, 035, 011, 009, 015, 014, 047, 029, 007, ! dy1(1:10,31, 7)
     .-016,-015, 016,-001,-010,-006,-004, 026, 012, 005, ! dy1(1:10,32, 7)
     .-023,-013, 007,-006,-012,-013,-005, 031,-002,-015, ! dy1(1:10,33, 7)
     . 000,-011, 009, 003,-011,-012,-003, 030,-013,-008, ! dy1(1:10,34, 7)
     .-028, 000,-001,-004,-011,-014,-007, 036,-028,-010, ! dy1(1:10,35, 7)
     .-014,-032,-014,-023,-024,-021,-011, 023,-056,-016, ! dy1(1:10,36, 7)
     .-021, 009,-025,-037,-025,-023,-018, 021,-073,-018, ! dy1(1:10,37, 7)
     .-059,-061,-059,-078,-060,-028,-034,-005,-102,-044, ! dy1(1:10,38, 7)
     . 001, 023,-001, 005, 028, 005, 004,-028,-005, 109, ! dy1(1:10, 1, 8)
     . 011, 021, 001,-008, 001, 005,-006,-025,-010,-006, ! dy1(1:10, 2, 8)
     . 000,-009, 000,-013,-004,-007,-015,-034,-017,-014, ! dy1(1:10, 3, 8)
     . 003,-009, 001, 002, 002,-005,-010,-016,-008,-007, ! dy1(1:10, 4, 8)
     .-001, 004, 009, 009, 009, 003, 003,-002,-005, 009, ! dy1(1:10, 5, 8)
     . 006, 001, 013, 009, 007, 001, 000, 000,-001,-006, ! dy1(1:10, 6, 8)
     . 002, 003, 014, 023, 013, 012, 008, 019, 006, 017, ! dy1(1:10, 7, 8)
     .-017,-026,-013,-005,-004,-009,-011, 000,-014,-015, ! dy1(1:10, 8, 8)
     .-034,-036,-020,-017,-021,-018,-025,-012,-020,-021, ! dy1(1:10, 9, 8)
     . 000,-028,-012,-018,-013,-013,-025,-010,-021,-022, ! dy1(1:10,10, 8)
     .-017,-010,-014,-010,-001,-004,-003, 001,-001,-013, ! dy1(1:10,11, 8)
     .-024,-016,-020,-012,-011,-003,-019,-013,-006,-009, ! dy1(1:10,12, 8)
     . 010,-004,-001, 000, 006, 012, 000,-007, 018,-001, ! dy1(1:10,13, 8)
     .-023,-015,-026,-023,-010,-010,-018,-028, 002,-026, ! dy1(1:10,14, 8)
     .-012,-007,-032,-026,-017,-010,-017,-033,-007,-019, ! dy1(1:10,15, 8)
     .-007,-004,-024,-018,-006, 000,-011,-029, 002,-014, ! dy1(1:10,16, 8)
     .-004,-006,-025,-016,-005, 000,-011,-033,-002,-009, ! dy1(1:10,17, 8)
     . 015, 001,-021,-013,-009, 002,-008,-028,-005,-001, ! dy1(1:10,18, 8)
     . 003, 013,-008,-003, 000,-002,-005,-024,-007, 009, ! dy1(1:10,19, 8)
     . 014, 013,-015,-014,-001,-002,-014,-032,-016, 003, ! dy1(1:10,20, 8)
     . 014, 008,-022,-012,-010,-011,-014,-037,-015,-003, ! dy1(1:10,21, 8)
     . 014, 009,-010,-015,-010,-003,-003,-028,-010,-004, ! dy1(1:10,22, 8)
     . 007, 016,-006,-005,-004,-005,-002,-025,-010, 010, ! dy1(1:10,23, 8)
     . 003, 008,-005,-012,-008,-009,-001,-015,-008, 002, ! dy1(1:10,24, 8)
     . 002, 005, 004, 000,-009,-014,-001,-003,-011,-004, ! dy1(1:10,25, 8)
     .-020,-006,-023,-032,-026,-026,-015,-019,-025,-013, ! dy1(1:10,26, 8)
     .-021,-005,-007,-009,-023,-018,-007,-004,-017,-016, ! dy1(1:10,27, 8)
     .-009,-001, 010, 004,-001,-008, 011, 008,-008,-001, ! dy1(1:10,28, 8)
     .-020,-013, 007, 004,-009,-008, 012, 013, 000,-006, ! dy1(1:10,29, 8)
     .-005,-002, 020, 012, 000, 002, 015, 023, 012,-001, ! dy1(1:10,30, 8)
     .-014,-001, 014, 014, 000,-007, 017, 025, 017, 004, ! dy1(1:10,31, 8)
     .-009,-023,-002,-002,-012,-013,-001, 010, 000,-023, ! dy1(1:10,32, 8)
     .-026,-022, 000,-011,-019,-022,-015,-008,-008,-019, ! dy1(1:10,33, 8)
     .-015,-015,-011,-017,-017,-026,-022,-004,-013,-020, ! dy1(1:10,34, 8)
     .-019,-024,-005,-011,-013,-018,-009, 000,-018,-020, ! dy1(1:10,35, 8)
     .-027,-026,-016,-027,-028,-037,-033,-010,-042,-036, ! dy1(1:10,36, 8)
     .-039,-034,-013,-025,-017,-028,-009,-001,-059,-013, ! dy1(1:10,37, 8)
     .-044,-018,-027,-033,-043,-051,-046,-026,-087,-035, ! dy1(1:10,38, 8)
     . 048, 003, 042, 024, 040, 021, 023,-004, 021, 067, ! dy1(1:10, 1, 9)
     . 023, 008, 012, 000, 009, 007,-010,-014,-003, 012, ! dy1(1:10, 2, 9)
     . 007,-007, 000,-004, 005,-001,-001,-020,-005,-001, ! dy1(1:10, 3, 9)
     . 010, 004, 011, 006, 004, 000,-008,-008,-003,-001, ! dy1(1:10, 4, 9)
     . 026, 009, 025, 021, 022, 018, 014, 014, 010, 024, ! dy1(1:10, 5, 9)
     .-006,-002, 009, 007, 006, 000,-002, 006,-003, 003, ! dy1(1:10, 6, 9)
     . 010, 011, 016, 018, 014, 013, 007, 015,-002,-007, ! dy1(1:10, 7, 9)
     .-009,-016,-003, 007, 004, 004,-004, 012,-006,-013, ! dy1(1:10, 8, 9)
     . 001,-015, 006, 000, 006, 007,-005, 016, 001, 007, ! dy1(1:10, 9, 9)
     .-001,-005, 004, 006, 003, 007, 002, 011, 005, 003, ! dy1(1:10,10, 9)
     . 005, 003, 008, 006, 017, 021, 006, 015, 017, 004, ! dy1(1:10,11, 9)
     . 000, 006, 000, 000, 008, 013, 000, 005, 007, 001, ! dy1(1:10,12, 9)
     . 023, 028, 010, 012, 025, 032, 021, 006, 038, 017, ! dy1(1:10,13, 9)
     . 030, 035, 004, 017, 018, 026, 013, 003, 027, 016, ! dy1(1:10,14, 9)
     . 030, 033, 008, 013, 028, 032, 027, 008, 044, 014, ! dy1(1:10,15, 9)
     . 032, 030, 017, 023, 029, 038, 024, 010, 041, 037, ! dy1(1:10,16, 9)
     . 035, 048, 023, 025, 036, 045, 035, 012, 045, 032, ! dy1(1:10,17, 9)
     . 046, 031, 013, 023, 027, 028, 020, 004, 027, 033, ! dy1(1:10,18, 9)
     . 048, 044, 023, 024, 029, 034, 027, 006, 027, 038, ! dy1(1:10,19, 9)
     . 009, 019,-012,-006, 005, 002,-003,-022, 001, 002, ! dy1(1:10,20, 9)
     . 011, 012,-003, 000, 001, 002, 000,-026,-006, 014, ! dy1(1:10,21, 9)
     . 023, 027, 000, 002, 006, 013, 006,-020, 004, 019, ! dy1(1:10,22, 9)
     . 022, 019, 007, 007, 001, 005, 003,-014,-002, 014, ! dy1(1:10,23, 9)
     . 019, 023, 005,-002, 000, 000, 005,-013, 005, 013, ! dy1(1:10,24, 9)
     . 003, 014, 003, 000, 004, 003, 017, 010, 007, 014, ! dy1(1:10,25, 9)
     . 007, 000, 001,-006,-003,-001, 008,-001, 000, 007, ! dy1(1:10,26, 9)
     .-003,-003, 000,-005,-014,-013,-003, 000,-012,-007, ! dy1(1:10,27, 9)
     . 012, 016, 026, 015, 010, 004, 026, 028, 008, 016, ! dy1(1:10,28, 9)
     . 008, 014, 027, 029, 015, 015, 030, 039, 011, 015, ! dy1(1:10,29, 9)
     . 015, 009, 033, 021, 014, 015, 026, 045, 021, 023, ! dy1(1:10,30, 9)
     . 022, 010, 028, 032, 016, 013, 029, 037, 029, 014, ! dy1(1:10,31, 9)
     . 009, 012, 038, 036, 023, 022, 031, 040, 035, 025, ! dy1(1:10,32, 9)
     . 000, 000, 014, 022, 005, 007, 011, 021, 016,-001, ! dy1(1:10,33, 9)
     .-001, 012, 021, 018, 016, 004, 019, 026, 020, 007, ! dy1(1:10,34, 9)
     . 017, 015, 027, 030, 021, 012, 022, 032, 016, 018, ! dy1(1:10,35, 9)
     .-016, 000, 014, 004, 002,-013, 006, 021,-012, 000, ! dy1(1:10,36, 9)
     .-010, 023, 019, 001, 008, 001, 009, 022,-027, 003, ! dy1(1:10,37, 9)
     .-041,-005,-025,-025,-035,-009,-027, 000,-059,-019, ! dy1(1:10,38, 9)
     . 031, 045, 025, 026, 024, 029, 015, 012, 027, 038, ! dy1(1:10, 1,10)
     .-007, 003, 003, 003, 005, 004,-006,-018, 000, 000, ! dy1(1:10, 2,10)
     . 004, 000, 011, 010, 002, 007,-007,-005, 000, 001, ! dy1(1:10, 3,10)
     .-020,-009, 003, 004,-005, 001,-007,-007,-007,-005, ! dy1(1:10, 4,10)
     .-009, 000, 008, 014, 004, 009, 000, 019, 000, 004, ! dy1(1:10, 5,10)
     .-006,-010, 001, 012, 005, 004,-001, 029,-001,-008, ! dy1(1:10, 6,10)
     .-023,-020,-011, 000,-004, 010,-008, 039,-004, 007, ! dy1(1:10, 7,10)
     .-026,-018,-015, 000,-005, 001,-013, 037,-009,-011, ! dy1(1:10, 8,10)
     .-027,-027,-019,-007,-015,-002,-014, 023,-010,-015, ! dy1(1:10, 9,10)
     .-020,-014,-018,-009,-003, 010,-008, 024,-003,-010, ! dy1(1:10,10,10)
     . 000, 000,-003, 004, 005, 015, 004, 029, 011, 000, ! dy1(1:10,11,10)
     . 000, 003,-007, 004, 010, 017, 004, 028, 015, 005, ! dy1(1:10,12,10)
     . 006, 004,-013, 005, 007, 015, 000, 023, 023, 002, ! dy1(1:10,13,10)
     . 004, 007, 004, 010, 026, 027, 010, 030, 035, 010, ! dy1(1:10,14,10)
     . 023, 019,-001, 010, 017, 024, 010, 023, 034, 024, ! dy1(1:10,15,10)
     . 015, 018, 000, 009, 026, 019, 011, 018, 035, 013, ! dy1(1:10,16,10)
     . 024, 017, 000, 015, 023, 030, 013, 016, 034, 033, ! dy1(1:10,17,10)
     . 028, 022, 007, 018, 025, 031, 012, 002, 031, 028, ! dy1(1:10,18,10)
     . 035, 035, 011, 020, 031, 032, 024, 002, 035, 036, ! dy1(1:10,19,10)
     . 012, 001,-017,-009,-006,-001,-011,-034,-002, 006, ! dy1(1:10,20,10)
     . 011, 009,-018,-007, 002, 001,-006,-031, 004, 010, ! dy1(1:10,21,10)
     . 012, 011,-014,-002,-001, 001,-002,-030, 000, 008, ! dy1(1:10,22,10)
     . 012, 010,-012,-007,-002, 002, 001,-029, 000, 012, ! dy1(1:10,23,10)
     .-003, 004,-006,-006, 000, 000, 004,-023,-004, 009, ! dy1(1:10,24,10)
     . 015, 015, 014,-002, 008, 002, 009,-009,-007, 001, ! dy1(1:10,25,10)
     . 004, 006, 000,-007,-009,-002, 000,-017,-016, 001, ! dy1(1:10,26,10)
     . 004, 014, 015, 001, 008, 003, 019, 013,-007, 006, ! dy1(1:10,27,10)
     . 000, 010, 014, 009, 003, 008, 015, 019,-007, 000, ! dy1(1:10,28,10)
     . 017, 014, 023, 012, 010, 012, 021, 029, 005, 008, ! dy1(1:10,29,10)
     . 013, 007, 017, 017, 006, 006, 019, 021, 008, 003, ! dy1(1:10,30,10)
     .-018,-007, 023, 013, 004, 009, 017, 005, 017,-001, ! dy1(1:10,31,10)
     . 009, 007, 021, 027, 014, 014, 024, 013, 019, 008, ! dy1(1:10,32,10)
     .-009,-003, 021, 015, 007, 000, 014,-012, 014, 005, ! dy1(1:10,33,10)
     . 004, 005, 023, 020, 012, 003, 017,-007, 019, 003, ! dy1(1:10,34,10)
     . 010, 012, 027, 024, 019, 006, 020, 000, 019, 018, ! dy1(1:10,35,10)
     . 000, 007, 012, 008, 011,-008, 004,-005, 008, 000, ! dy1(1:10,36,10)
     . 020, 006, 006, 003, 012,-001, 014, 005,-007, 004, ! dy1(1:10,37,10)
     .-003,-010,-013,-015,-003,-028,-013,-029,-034,-021, ! dy1(1:10,38,10)
     . 054, 025, 037, 021, 033, 039, 030, 027, 032, 032, ! dy1(1:10, 1,11)
     . 001, 009, 011, 011, 009, 014, 001, 000, 021, 010, ! dy1(1:10, 2,11)
     .-002,-011, 007, 008, 001, 004,-007, 001, 002, 011, ! dy1(1:10, 3,11)
     .-013,-014, 006, 003,-004, 002,-012, 007, 001,-001, ! dy1(1:10, 4,11)
     . 000,-007, 008, 005, 001, 007,-001, 024, 003,-003, ! dy1(1:10, 5,11)
     .-026,-013, 004, 004, 000, 006,-002, 032, 003, 000, ! dy1(1:10, 6,11)
     .-008,-014,-007, 000,-008, 001,-009, 031,-001,-018, ! dy1(1:10, 7,11)
     .-029,-034,-024,-008,-024,-005,-018, 019,-008,-017, ! dy1(1:10, 8,11)
     .-013,-016,-013,-005,-005, 001,-008, 036,-006,-009, ! dy1(1:10, 9,11)
     .-024,-016,-017,-004,-007, 000,-005, 031,-001,-011, ! dy1(1:10,10,11)
     .-018, 002,-005,-003, 001, 006,-001, 032, 006, 000, ! dy1(1:10,11,11)
     .-013,-004,-015,-005,-001, 005,-005, 021, 004, 002, ! dy1(1:10,12,11)
     .-019,-012,-019,-008, 000, 006,-002, 021, 009,-012, ! dy1(1:10,13,11)
     .-009, 003,-012, 004, 004, 015, 002, 031, 017, 003, ! dy1(1:10,14,11)
     . 003, 006,-010, 007, 018, 017, 008, 021, 027, 012, ! dy1(1:10,15,11)
     . 006, 008,-007, 009, 012, 018, 008, 016, 029, 019, ! dy1(1:10,16,11)
     . 021, 018, 002, 009, 020, 026, 010, 014, 034, 025, ! dy1(1:10,17,11)
     . 017, 019, 003, 016, 026, 028, 009, 011, 032, 028, ! dy1(1:10,18,11)
     . 027, 018, 003, 024, 028, 033, 015, 002, 035, 036, ! dy1(1:10,19,11)
     . 004, 000,-012,-004, 004, 005,-010,-022, 013, 008, ! dy1(1:10,20,11)
     . 006,-004,-021,-008, 000,-002,-011,-033, 001, 012, ! dy1(1:10,21,11)
     .-006,-006,-019,-009,-003,-005,-011,-035, 000, 000, ! dy1(1:10,22,11)
     . 003,-006,-016,-010,-001,-009,-004,-028,-007, 004, ! dy1(1:10,23,11)
     . 000,-003,-008,-008,-003,-006,-003,-023,-005,-002, ! dy1(1:10,24,11)
     . 012, 003, 004,-002, 000, 000, 006,-007,-001,-002, ! dy1(1:10,25,11)
     . 001, 002,-006,-008,-007,-006, 000,-015,-010, 000, ! dy1(1:10,26,11)
     . 012, 001, 012, 000, 003, 006, 005, 009,-009,-006, ! dy1(1:10,27,11)
     . 009, 008, 013, 002, 007, 007, 013, 015,-003,-002, ! dy1(1:10,28,11)
     . 005, 010, 020, 007, 007, 006, 013, 010, 001, 002, ! dy1(1:10,29,11)
     . 003, 000, 011, 000,-002, 006, 005, 000,-002,-007, ! dy1(1:10,30,11)
     . 007,-005, 008, 013, 010, 008, 014,-007, 001, 004, ! dy1(1:10,31,11)
     . 011,-003, 022, 019, 004, 006, 014,-007, 011,-008, ! dy1(1:10,32,11)
     . 000,-009, 009, 005, 000, 001, 008,-023, 003, 002, ! dy1(1:10,33,11)
     .-011,-003, 018, 005, 003,-005, 007,-025, 011,-006, ! dy1(1:10,34,11)
     . 006,-008, 016, 012, 010, 006, 012,-022, 020, 006, ! dy1(1:10,35,11)
     .-001, 000, 008, 002, 001,-009, 007,-023, 009, 007, ! dy1(1:10,36,11)
     . 004, 004, 005, 006, 011, 003, 010, 004, 014, 001, ! dy1(1:10,37,11)
     .-034,-034,-032,-026,-019,-030,-018,-026,-020,-016, ! dy1(1:10,38,11)
     . 022, 016, 038, 036, 030, 045, 025, 038, 046, 037, ! dy1(1:10, 1,12)
     .-005, 014, 021, 014, 010, 012, 008, 018, 019, 019, ! dy1(1:10, 2,12)
     .-011,-013, 010, 007,-006, 000,-003, 008, 010,-001, ! dy1(1:10, 3,12)
     .-015,-013, 003,-001,-007,-004,-009, 009, 006,-004, ! dy1(1:10, 4,12)
     .-013,-017, 015, 010, 002, 003, 003, 025, 016, 000, ! dy1(1:10, 5,12)
     .-009,-003, 005, 001,-005,-003,-001, 024, 006,-012, ! dy1(1:10, 6,12)
     .-027,-034, 000,-003,-006,-003,-001, 029, 000,-009, ! dy1(1:10, 7,12)
     .-023,-034,-007,-010,-015,-014,-014, 018,-006,-024, ! dy1(1:10, 8,12)
     .-030,-013,-016,-009,-016,-010,-010, 019,-012,-021, ! dy1(1:10, 9,12)
     .-009,-008,-014,-010,-006,-007,-007, 024,-006,-012, ! dy1(1:10,10,12)
     .-011,-006, 001, 002, 001, 001, 007, 030, 001, 004, ! dy1(1:10,11,12)
     .-012,-010,-013,-003, 000, 003, 004, 023, 000,-007, ! dy1(1:10,12,12)
     .-005,-005,-011,-006,-002, 000,-001, 017, 005, 008, ! dy1(1:10,13,12)
     . 007,-003,-007, 004, 008, 007, 009, 027, 011, 005, ! dy1(1:10,14,12)
     . 003, 004,-003, 003, 008, 013, 007, 017, 020, 013, ! dy1(1:10,15,12)
     . 014, 012,-002, 009, 019, 021, 009, 019, 027, 017, ! dy1(1:10,16,12)
     . 017, 013, 002, 018, 024, 019, 012, 017, 037, 027, ! dy1(1:10,17,12)
     . 014, 015, 000, 016, 020, 022, 004, 009, 025, 032, ! dy1(1:10,18,12)
     . 027, 025, 008, 017, 029, 025, 017, 006, 035, 022, ! dy1(1:10,19,12)
     . 002,-001,-020,-008, 001,-002,-015,-022, 008, 011, ! dy1(1:10,20,12)
     .-010,-011,-018,-009,-001,-005,-017,-027,-003,-005, ! dy1(1:10,21,12)
     .-002,-010,-021,-015,-003,-009,-017,-028,-009, 003, ! dy1(1:10,22,12)
     . 006,-007,-016,-005,-002,-010,-011,-023, 000, 000, ! dy1(1:10,23,12)
     .-003,-004,-012,-010,-001,-010,-006,-015,-002, 000, ! dy1(1:10,24,12)
     . 008, 007, 000, 000, 008, 000, 001,-002,-002,-005, ! dy1(1:10,25,12)
     .-004,-010,-010,-010, 000,-007,-004,-008,-003,-003, ! dy1(1:10,26,12)
     . 003, 002, 003,-008,-002, 001, 000, 004,-005,-009, ! dy1(1:10,27,12)
     . 008, 000, 006,-006,-003,-002, 006, 006,-003,-010, ! dy1(1:10,28,12)
     . 001, 000, 014,-003, 001, 004, 007, 005, 000,-003, ! dy1(1:10,29,12)
     .-004, 001, 005,-002, 000, 004, 004,-001,-004,-004, ! dy1(1:10,30,12)
     .-009,-022, 003, 000,-007,-004, 009,-014,-003, 000, ! dy1(1:10,31,12)
     .-001,-006, 008, 009, 005, 007, 021,-010, 001, 011, ! dy1(1:10,32,12)
     . 000,-011, 012,-002,-003,-007, 003,-027, 000,-008, ! dy1(1:10,33,12)
     .-006,-012, 010, 007, 004,-006, 003,-027, 011,-003, ! dy1(1:10,34,12)
     . 004, 000, 017, 006, 007, 002, 012,-025, 013, 005, ! dy1(1:10,35,12)
     . 007,-011, 005, 000, 004, 000, 008,-022, 019,-003, ! dy1(1:10,36,12)
     . 006, 014, 002, 003, 013, 011, 012,-008, 027, 021, ! dy1(1:10,37,12)
     .-015,-009,-026,-012,-015,-010,-009,-021, 005,-012, ! dy1(1:10,38,12)
     . 031, 035, 032, 036, 033, 031, 019, 026, 051, 019, ! dy1(1:10, 1,13)
     .-003,-012, 004, 001, 003, 011, 002, 008, 015, 013, ! dy1(1:10, 2,13)
     .-017,-012, 004,-004,-005,-002,-009, 006, 011,-007, ! dy1(1:10, 3,13)
     .-027,-020,-001,-008,-012,-011,-009, 002, 003,-005, ! dy1(1:10, 4,13)
     .-017,-017, 008, 000,-009,-007, 000, 009, 009,-004, ! dy1(1:10, 5,13)
     .-013,-023, 002, 002,-003,-008,-005, 013, 008,-006, ! dy1(1:10, 6,13)
     .-013,-004, 000,-002,-007,-003, 000, 021,-003,-013, ! dy1(1:10, 7,13)
     .-012,-010,-001,-003,-008,-007,-001, 018,-007,-001, ! dy1(1:10, 8,13)
     .-024,-017,-013, 000,-011,-009,-008, 009,-015,-011, ! dy1(1:10, 9,13)
     .-009,-017,-013,-006,-009,-019,-005, 009,-017,-017, ! dy1(1:10,10,13)
     .-022, 007, 000, 002, 000,-003, 007, 022, 001, 000, ! dy1(1:10,11,13)
     .-012, 001,-001, 005, 004, 000, 003, 021,-002, 003, ! dy1(1:10,12,13)
     .-019,-022,-009,-016,-008, 000,-003, 005,-010,-009, ! dy1(1:10,13,13)
     .-020,-011,-011,-008, 000, 000, 003, 008,-007,-001, ! dy1(1:10,14,13)
     .-004, 005,-004, 002, 004, 011, 003, 012, 007, 008, ! dy1(1:10,15,13)
     . 000, 002, 002, 009, 011, 009, 008, 013, 012, 019, ! dy1(1:10,16,13)
     . 007, 028, 004, 024, 023, 021, 010, 017, 026, 023, ! dy1(1:10,17,13)
     . 012, 007, 000, 016, 018, 013, 006, 012, 016, 017, ! dy1(1:10,18,13)
     . 021, 011, 000, 020, 021, 008, 001, 002, 015, 027, ! dy1(1:10,19,13)
     .-006,-018,-026,-008,-004,-013,-020,-023,-007,-009, ! dy1(1:10,20,13)
     .-013,-004,-024,-007,-004,-017,-017,-021,-010,-005, ! dy1(1:10,21,13)
     . 002,-005,-012,-006,-004,-009,-016,-010,-006, 001, ! dy1(1:10,22,13)
     . 010, 006,-008,-011, 005,-007,-006,-008, 003, 000, ! dy1(1:10,23,13)
     . 013, 004,-002,-001, 006,-003,-002, 001, 006, 000, ! dy1(1:10,24,13)
     . 004,-004, 001,-003, 013, 000, 007, 013, 010,-004, ! dy1(1:10,25,13)
     .-006,-007, 002,-006, 003, 000,-002, 003, 001,-010, ! dy1(1:10,26,13)
     .-004,-019,-013,-010,-009,-010,-007, 003,-007,-012, ! dy1(1:10,27,13)
     .-005,-004,-003,-018,-004,-008,-003, 002,-007,-017, ! dy1(1:10,28,13)
     .-008,-002, 006,-007,-004, 005, 003, 012,-002,-003, ! dy1(1:10,29,13)
     .-025,-012,-010,-016,-016,-012,-002,-006,-011,-019, ! dy1(1:10,30,13)
     .-010,-021,-002,-016,-009,-010,-007,-019,-015,-010, ! dy1(1:10,31,13)
     .-007,-007, 000,-006,-008,-011, 006,-008,-001, 000, ! dy1(1:10,32,13)
     .-018,-022, 004,-008,-003,-008,-004,-021,-006,-016, ! dy1(1:10,33,13)
     .-007,-019, 003,-007,-007,-004, 001,-024, 002, 001, ! dy1(1:10,34,13)
     . 000,-006, 009,-005, 000,-003, 005,-022, 018, 002, ! dy1(1:10,35,13)
     . 004,-007, 002, 002, 002, 000, 005,-013, 021, 003, ! dy1(1:10,36,13)
     .-012, 008,-002,-004, 003, 008, 004,-015, 026, 005, ! dy1(1:10,37,13)
     .-032,-029,-029,-030,-024,-009,-009,-040,-007,-013, ! dy1(1:10,38,13)
     . 022, 019, 012, 022, 021, 038, 031, 042, 038, 047, ! dy1(1:10, 1,14)
     . 002, 009, 002, 003, 007, 000,-001, 009, 012, 012, ! dy1(1:10, 2,14)
     .-018,-008,-013,-007,-007,-003,-006, 002,-002,-001, ! dy1(1:10, 3,14)
     .-023,-020,-005,-013,-016,-013,-009, 000, 000,-004, ! dy1(1:10, 4,14)
     .-030,-005,-005,-009,-007,-007,-010,-002, 000,-013, ! dy1(1:10, 5,14)
     .-013,-013, 001,-005,-014,-008,-003,-004, 000,-012, ! dy1(1:10, 6,14)
     .-018,-016,-006,-007,-006,-016,-003,-009,-011,-012, ! dy1(1:10, 7,14)
     .-012,-010,-005,-001,-012,-015, 000,-004,-012,-017, ! dy1(1:10, 8,14)
     .-021,-016,-001,-004,-010,-015,-005,-005,-011,-020, ! dy1(1:10, 9,14)
     .-006, 006,-006,-005,-008,-017, 000,-007,-012,-008, ! dy1(1:10,10,14)
     .-025,-012, 000,-001,-005,-010, 000, 001,-013,-009, ! dy1(1:10,11,14)
     .-011, 000,-001, 000,-004,-007, 002, 001,-012, 005, ! dy1(1:10,12,14)
     .-021,-021,-010,-011,-005,-005,-013,-006,-018,-009, ! dy1(1:10,13,14)
     .-016,-016,-001,-006,-007,-004,-002,-006,-009, 000, ! dy1(1:10,14,14)
     .-012,-001, 001, 006, 006, 007, 000, 004, 000, 007, ! dy1(1:10,15,14)
     . 001, 011, 000, 011, 016, 009, 008, 016, 006, 011, ! dy1(1:10,16,14)
     . 022, 012, 011, 022, 019, 014, 010, 023, 009, 022, ! dy1(1:10,17,14)
     . 004, 012, 004, 024, 021, 012, 006, 020, 011, 021, ! dy1(1:10,18,14)
     . 021, 026, 008, 022, 025, 012, 013, 023, 015, 021, ! dy1(1:10,19,14)
     .-006,-005,-014,-002,-003,-016,-020,-009,-020,-009, ! dy1(1:10,20,14)
     .-010,-020,-023,-007,-005,-015,-019,-004,-014,-013, ! dy1(1:10,21,14)
     .-004,-009,-013,-009,-003,-016,-018, 000,-011,-014, ! dy1(1:10,22,14)
     . 005, 000, 000, 004, 010, 001,-004, 022, 005,-001, ! dy1(1:10,23,14)
     . 024, 016, 021, 017, 028, 013, 014, 040, 031, 006, ! dy1(1:10,24,14)
     . 021, 018, 013, 011, 016, 011, 011, 028, 019, 009, ! dy1(1:10,25,14)
     . 007,-003, 000, 002, 006,-001, 000, 016, 008,-007, ! dy1(1:10,26,14)
     .-001,-011,-006,-010,-006,-004,-009, 007,-005,-012, ! dy1(1:10,27,14)
     .-004,-015,-002,-009,-010,-007,-002, 008,-009,-016, ! dy1(1:10,28,14)
     .-002,-005, 000,-009,-006,-003, 001, 011,-002,-011, ! dy1(1:10,29,14)
     .-016,-023,-011,-017,-014,-012,-007, 001,-023,-019, ! dy1(1:10,30,14)
     .-016,-017,-010,-020,-019,-015,-006,-007,-018,-012, ! dy1(1:10,31,14)
     .-008,-024,-011,-015,-009,-010, 000, 002,-016,-006, ! dy1(1:10,32,14)
     .-002,-010, 001,-006,-010,-010, 001,-008,-005,-002, ! dy1(1:10,33,14)
     .-001,-008, 003, 000, 000,-001, 008,-011, 000,-007, ! dy1(1:10,34,14)
     .-005,-014, 006, 000,-003, 008, 003,-001, 005, 006, ! dy1(1:10,35,14)
     .-001,-008, 000,-006, 000, 010, 008,-009, 013, 003, ! dy1(1:10,36,14)
     . 020, 014,-005, 007, 005, 015, 006, 006, 017, 025, ! dy1(1:10,37,14)
     .-020, 005,-019,-032,-024, 003,-019,-023, 003,-017, ! dy1(1:10,38,14)
     . 044, 031, 027, 028, 023, 034, 020, 029, 023, 032, ! dy1(1:10, 1,15)
     . 000, 009,-013,-003, 001, 016, 006, 010,-001, 018, ! dy1(1:10, 2,15)
     . 012,-003, 000, 002,-002, 000,-003, 003,-003, 007, ! dy1(1:10, 3,15)
     .-007, 000,-011,-005,-006,-004,-004,-006,-019,-006, ! dy1(1:10, 4,15)
     .-014,-020,-010,-004,-011,-010,-006,-012,-014,-007, ! dy1(1:10, 5,15)
     .-019,-009,-005,-003,-006,-010, 001,-020,-014,-005, ! dy1(1:10, 6,15)
     .-016,-016,-009,-007,-019,-018,-009,-032,-013,-017, ! dy1(1:10, 7,15)
     .-026,-031,-003,-007,-020,-023,-012,-033,-013,-018, ! dy1(1:10, 8,15)
     .-020,-009,-008,-006,-013,-015,-004,-028,-018,-007, ! dy1(1:10, 9,15)
     .-018,-008,-002,-003,-010,-016, 000,-020,-014,-010, ! dy1(1:10,10,15)
     .-005,-011,-002,-001,-003,-010,-002,-016,-006,-001, ! dy1(1:10,11,15)
     .-012,-012, 005,-005,-003,-006,-001,-014,-011,-006, ! dy1(1:10,12,15)
     .-006, 004, 006, 007, 004,-003, 006, 000, 001, 008, ! dy1(1:10,13,15)
     .-014, 004, 004, 008, 005, 000, 003,-002,-006, 007, ! dy1(1:10,14,15)
     . 001, 009, 007, 012, 014, 007, 009, 010, 000, 012, ! dy1(1:10,15,15)
     .-009, 005, 008, 015, 013, 001, 009, 009,-001, 005, ! dy1(1:10,16,15)
     . 006, 006, 010, 022, 025, 021, 018, 024, 012, 022, ! dy1(1:10,17,15)
     . 003, 004, 010, 023, 020, 010, 010, 023, 003, 019, ! dy1(1:10,18,15)
     . 010, 015, 010, 024, 028, 016, 013, 034, 012, 015, ! dy1(1:10,19,15)
     .-016,-011,-009, 002, 001,-009,-008, 011,-013,-008, ! dy1(1:10,20,15)
     .-013,-016,-010,-005,-004,-019,-019, 003,-023,-012, ! dy1(1:10,21,15)
     .-004,-003,-007, 003, 001,-006,-014, 021,-008,-006, ! dy1(1:10,22,15)
     . 001,-001,-001,-002, 003,-006,-011, 025,-005,-005, ! dy1(1:10,23,15)
     . 005, 004,-003, 001,-002,-008,-011, 030, 000,-016, ! dy1(1:10,24,15)
     .-007,-011,-003,-001, 007,-001,-008, 023, 001,-014, ! dy1(1:10,25,15)
     .-006,-008,-004,-009,-007,-005,-007, 020, 005,-011, ! dy1(1:10,26,15)
     . 003,-015,-001,-005,-002,-007,-008, 016,-008,-013, ! dy1(1:10,27,15)
     . 000,-001, 000,-004, 000, 002, 000, 026,-001,-006, ! dy1(1:10,28,15)
     .-008,-002,-003,-009,-006, 002, 003, 022,-013,-011, ! dy1(1:10,29,15)
     . 004,-008,-007,-018,-013,-014,-003, 011,-015, 000, ! dy1(1:10,30,15)
     .-020,-017,-008,-013,-022,-008,-008, 007,-022,-021, ! dy1(1:10,31,15)
     .-011,-015,-002,-009,-010,-008,-003, 013,-012,-009, ! dy1(1:10,32,15)
     . 003,-008, 004, 000,-005,-002, 009, 014,-006, 003, ! dy1(1:10,33,15)
     . 000, 000, 005,-005,-008, 002, 009, 010,-003, 009, ! dy1(1:10,34,15)
     . 015, 006, 012, 008, 001, 012, 014, 009, 010, 009, ! dy1(1:10,35,15)
     . 021, 003, 007, 012, 003, 019, 009, 013, 014, 022, ! dy1(1:10,36,15)
     . 036, 036, 021, 017, 018, 039, 025, 009, 035, 033, ! dy1(1:10,37,15)
     . 015,-005,-021, 005,-001, 021,-008,-022,-005, 030, ! dy1(1:10,38,15)
     . 056, 076, 005, 019, 040, 055, 048, 045, 029, 061, ! dy1(1:10, 1,16)
     . 029, 021,-006, 000, 011, 020, 006, 005, 002, 028, ! dy1(1:10, 2,16)
     .-001, 019,-008, 009, 007, 022, 009,-001,-012, 016, ! dy1(1:10, 3,16)
     . 011, 009,-009, 004, 004, 007, 006,-005,-009, 003, ! dy1(1:10, 4,16)
     .-010, 000,-002,-004,-005,-003,-001,-020,-022, 000, ! dy1(1:10, 5,16)
     .-018,-003,-008, 002,-004,-003,-003,-022,-014,-005, ! dy1(1:10, 6,16)
     . 000,-011,-005,-005,-009,-011, 003,-038,-008,-007, ! dy1(1:10, 7,16)
     .-019,-010,-019,-015,-019,-017,-012,-051,-018,-012, ! dy1(1:10, 8,16)
     .-020,-024,-010,-006,-019,-021,-012,-043,-016,-019, ! dy1(1:10, 9,16)
     .-024,-009,-017,-012,-017,-025,-011,-042,-009,-015, ! dy1(1:10,10,16)
     .-026,-014, 002,-010,-014,-023,-005,-033,-009,-012, ! dy1(1:10,11,16)
     .-014, 000, 001, 002, 000,-013, 001,-019, 000,-001, ! dy1(1:10,12,16)
     .-016, 001, 014, 008, 005,-005, 005, 001, 002,-002, ! dy1(1:10,13,16)
     . 000, 005, 017, 010, 013,-002, 013, 005, 008, 003, ! dy1(1:10,14,16)
     . 000, 014, 014, 022, 021, 006, 011, 015, 004, 013, ! dy1(1:10,15,16)
     .-001, 013, 013, 018, 015, 000, 011, 015, 001, 008, ! dy1(1:10,16,16)
     . 013, 015, 021, 025, 027, 014, 018, 022, 010, 027, ! dy1(1:10,17,16)
     . 011, 003, 011, 023, 024, 017, 015, 021, 006, 010, ! dy1(1:10,18,16)
     . 006, 013, 021, 022, 024, 013, 012, 030, 005, 019, ! dy1(1:10,19,16)
     .-008,-020,-011,-006,-004,-012,-017, 004,-023,-020, ! dy1(1:10,20,16)
     .-022,-016,-011,-005,-007,-016,-016, 011,-024,-016, ! dy1(1:10,21,16)
     .-009,-014, 002,-001,-006,-015,-022, 013,-024,-025, ! dy1(1:10,22,16)
     .-011,-009,-010,-008,-010,-009,-017, 011,-019,-020, ! dy1(1:10,23,16)
     .-025,-029,-013,-017,-019,-029,-034, 001,-023,-029, ! dy1(1:10,24,16)
     .-004,-007,-005,-014,-010,-018,-016, 006,-009,-023, ! dy1(1:10,25,16)
     .-001,-015,-007,-011,-008,-011,-016, 009,-010,-022, ! dy1(1:10,26,16)
     .-011,-006,-001,-010,-006,-006,-007, 011,-006,-017, ! dy1(1:10,27,16)
     .-013,-015,-004,-011,-013,-007,-008, 005,-008,-013, ! dy1(1:10,28,16)
     . 005, 004, 007, 001, 000,-001, 003, 018,-002,-003, ! dy1(1:10,29,16)
     . 000,-004,-002,-005,-009, 001, 007, 019,-003,-008, ! dy1(1:10,30,16)
     . 011,-004, 000,-011,-011,-003, 008, 013,-017, 004, ! dy1(1:10,31,16)
     . 012, 008, 000,-010,-013,-001, 007, 018,-011, 004, ! dy1(1:10,32,16)
     . 008, 002, 004, 002,-006,-001, 012, 022,-001, 005, ! dy1(1:10,33,16)
     . 016, 007, 007, 008,-001, 005, 014, 016, 001, 008, ! dy1(1:10,34,16)
     . 020, 019, 004, 002,-002, 014, 012, 008, 005, 021, ! dy1(1:10,35,16)
     . 039, 022, 015, 024, 012, 023, 023, 016, 023, 022, ! dy1(1:10,36,16)
     . 058, 057, 031, 027, 028, 050, 028, 021, 043, 062, ! dy1(1:10,37,16)
     . 038, 026, 010, 010,-004, 025, 010,-002, 017, 031, ! dy1(1:10,38,16)
     . 086, 066, 041, 050, 056, 087, 049, 048, 033, 060, ! dy1(1:10, 1,17)
     . 028, 041, 001, 023, 026, 045, 018, 020, 001, 041, ! dy1(1:10, 2,17)
     . 024, 032, 007, 013, 017, 034, 018, 012, 001, 028, ! dy1(1:10, 3,17)
     . 018, 025, 007, 012, 020, 025, 011, 001,-004, 029, ! dy1(1:10, 4,17)
     . 014, 023,-002, 017, 012, 015, 014,-007, 001, 020, ! dy1(1:10, 5,17)
     . 000, 003,-005, 000,-004, 000, 001,-027,-003, 003, ! dy1(1:10, 6,17)
     .-018,-006,-009, 000,-008,-008,-005,-038,-004, 000, ! dy1(1:10, 7,17)
     .-026,-022,-022,-008,-015,-022,-011,-045,-009,-012, ! dy1(1:10, 8,17)
     .-021,-007,-010,-011,-012,-020,-007,-043,-001,-004, ! dy1(1:10, 9,17)
     .-017,-017,-010,-005,-010,-015,-013,-037, 000,-007, ! dy1(1:10,10,17)
     .-035,-011,-009,-010,-011,-023,-012,-034, 002,-015, ! dy1(1:10,11,17)
     .-001,-011,-004,-006,-005,-019,-008,-027, 002,-001, ! dy1(1:10,12,17)
     .-028,-006, 003, 001, 001,-012, 000,-004, 006,-002, ! dy1(1:10,13,17)
     .-019, 002, 006, 005, 005,-011, 001, 000, 009,-003, ! dy1(1:10,14,17)
     . 008, 012, 017, 012, 017, 002, 012, 013, 014, 006, ! dy1(1:10,15,17)
     .-006, 007, 015, 016, 018, 006, 009, 015, 008, 005, ! dy1(1:10,16,17)
     . 000, 009, 019, 017, 020, 010, 009, 020, 012, 010, ! dy1(1:10,17,17)
     . 005, 009, 020, 017, 020, 014, 015, 020, 014, 005, ! dy1(1:10,18,17)
     . 012, 008, 020, 020, 024, 014, 014, 032, 010, 007, ! dy1(1:10,19,17)
     .-016,-015,-005,-003,-006,-017,-017, 003,-022,-020, ! dy1(1:10,20,17)
     .-001,-007,-002,-008,-004,-013,-012, 003,-019,-024, ! dy1(1:10,21,17)
     .-007,-008, 000, 004,-002,-007,-010, 003,-017,-016, ! dy1(1:10,22,17)
     .-006,-001, 005,-003,-005,-015,-012, 000,-016,-014, ! dy1(1:10,23,17)
     .-002,-001,-004,-011,-008,-015,-013,-007,-016,-029, ! dy1(1:10,24,17)
     .-011,-003, 001, 003,-007,-015,-014,-006,-013,-010, ! dy1(1:10,25,17)
     .-012,-015,-011,-011,-016,-024,-019,-014,-020,-018, ! dy1(1:10,26,17)
     .-007,-013,-005,-005,-015,-014,-013,-015,-011,-023, ! dy1(1:10,27,17)
     .-005,-007,-003,-015,-010,-004,-007,-015,-009,-014, ! dy1(1:10,28,17)
     . 003, 000, 000,-002,-009, 000, 000,-013,-004,-014, ! dy1(1:10,29,17)
     . 009, 006, 006,-004,-008, 000, 007,-004, 001, 004, ! dy1(1:10,30,17)
     .-001, 006,-007,-002,-018,-009, 003, 002,-013,-006, ! dy1(1:10,31,17)
     . 002,-002, 000, 003,-014,-007, 010, 003,-011,-007, ! dy1(1:10,32,17)
     . 031, 021, 005,-003,-009, 000, 015, 015, 000, 014, ! dy1(1:10,33,17)
     . 033, 025, 004, 005,-008, 011, 014, 022, 000, 013, ! dy1(1:10,34,17)
     . 055, 032, 015, 020, 005, 025, 024, 023, 016, 041, ! dy1(1:10,35,17)
     . 064, 047, 023, 027, 012, 041, 028, 017, 030, 045, ! dy1(1:10,36,17)
     . 083, 071, 037, 035, 034, 067, 042, 027, 051, 080, ! dy1(1:10,37,17)
     . 034, 049, 010, 006, 010, 046, 016,-008, 023, 052, ! dy1(1:10,38,17)
     . 105, 093, 052, 069, 083, 112, 085, 089, 061, 116, ! dy1(1:10, 1,18)
     . 100, 084, 022, 044, 060, 102, 062, 061, 037, 082, ! dy1(1:10, 2,18)
     . 067, 062, 011, 031, 041, 070, 033, 039, 020, 061, ! dy1(1:10, 3,18)
     . 039, 047, 000, 025, 021, 052, 031, 026, 027, 047, ! dy1(1:10, 4,18)
     . 016, 022, 004, 005, 015, 035, 016, 013, 018, 043, ! dy1(1:10, 5,18)
     . 014, 022,-001, 011, 014, 019, 017,-002, 022, 029, ! dy1(1:10, 6,18)
     . 006, 006,-007, 004, 002, 001, 009,-024, 026, 010, ! dy1(1:10, 7,18)
     .-012, 001,-013,-010,-008,-006, 000,-028, 014, 001, ! dy1(1:10, 8,18)
     .-003,-013,-010, 005,-009,-006,-003,-035, 025, 000, ! dy1(1:10, 9,18)
     .-015,-002,-003,-004,-001,-013,-003,-027, 028,-002, ! dy1(1:10,10,18)
     .-010,-006,-010,-001,-003,-016,-004,-021, 018, 005, ! dy1(1:10,11,18)
     .-026,-001,-004, 000, 000,-009,-002,-008, 024,-015, ! dy1(1:10,12,18)
     .-030,-011,-015,-015,-007,-025,-018,-023, 008,-023, ! dy1(1:10,13,18)
     .-047,-016,-021,-016,-012,-026,-017,-018, 005,-024, ! dy1(1:10,14,18)
     .-013, 005, 009, 006, 012, 005, 005, 019, 027,-004, ! dy1(1:10,15,18)
     .-003, 016, 022, 015, 023, 019, 014, 029, 033, 004, ! dy1(1:10,16,18)
     . 010, 015, 023, 018, 019, 015, 013, 031, 026, 002, ! dy1(1:10,17,18)
     . 014, 022, 026, 019, 022, 018, 014, 030, 021, 010, ! dy1(1:10,18,18)
     . 011, 023, 036, 015, 020, 015, 018, 012, 021,-006, ! dy1(1:10,19,18)
     .-016, 000, 000,-015,-011,-016,-013,-018,-006,-032, ! dy1(1:10,20,18)
     .-003, 000, 002, 001,-003,-010,-009,-017,-009,-017, ! dy1(1:10,21,18)
     . 004, 001, 004,-013,-008,-015,-009,-028,-006,-025, ! dy1(1:10,22,18)
     . 000, 000, 000, 000,-007,-009,-009,-038,-016,-021, ! dy1(1:10,23,18)
     . 000, 007, 007, 013,-002,-003, 000,-037,-009,-005, ! dy1(1:10,24,18)
     . 015, 001, 004, 009,-008,-002, 000,-047,-008,-011, ! dy1(1:10,25,18)
     .-020,-008,-009,-004,-018,-011,-013,-062,-022,-026, ! dy1(1:10,26,18)
     . 001,-001,-010,-011,-021,-011,-011,-072,-015,-010, ! dy1(1:10,27,18)
     .-017,-007,-009,-009,-028,-012,-013,-066,-009,-016, ! dy1(1:10,28,18)
     . 018, 011, 004,-008,-017,-003, 001,-062, 003,-001, ! dy1(1:10,29,18)
     . 016, 015, 001, 014,-012, 008, 001,-042, 000, 007, ! dy1(1:10,30,18)
     . 066, 025, 010, 006,-017, 015, 017,-022, 012, 019, ! dy1(1:10,31,18)
     . 038, 037, 006, 008,-016, 010, 011,-011, 015, 027, ! dy1(1:10,32,18)
     . 048, 035, 002, 026,-002, 024, 019, 001, 019, 027, ! dy1(1:10,33,18)
     . 066, 057, 020, 030, 005, 033, 034, 017, 035, 050, ! dy1(1:10,34,18)
     . 066, 065, 020, 035, 007, 040, 029, 029, 041, 062, ! dy1(1:10,35,18)
     . 093, 075, 028, 049, 032, 066, 047, 035, 061, 096, ! dy1(1:10,36,18)
     . 135, 117, 060, 065, 049, 096, 066, 040, 091, 118, ! dy1(1:10,37,18)
     . 118, 085, 032, 044, 027, 075, 046, 016, 065, 100/ ! dy1(1:10,38,18)

      data dy2 / 
     . 051, 024, 098, 090, 014, 028, 083,-018, 050, 011, ! dy2(1:10, 1, 1)
     . 049, 018, 099, 073, 011, 014, 063,-020, 048, 005, ! dy2(1:10, 2, 1)
     . 019, 000, 068, 044,-018,-013, 032,-066, 009,-030, ! dy2(1:10, 3, 1)
     . 000,-043, 022, 007,-038,-027,-003,-083,-011,-036, ! dy2(1:10, 4, 1)
     .-019,-036, 005, 000,-047,-041,-002,-102,-018,-050, ! dy2(1:10, 5, 1)
     .-012,-042,-012,-012,-058,-050,-010,-109,-033,-048, ! dy2(1:10, 6, 1)
     .-050,-059,-035,-044,-067,-071,-031,-128,-029,-065, ! dy2(1:10, 7, 1)
     .-055,-056,-033,-032,-063,-055,-031,-130,-031,-056, ! dy2(1:10, 8, 1)
     .-029,-054,-031,-030,-061,-054,-034,-118,-032,-057, ! dy2(1:10, 9, 1)
     .-022,-033,-033,-028,-048,-047,-025,-103,-019,-043, ! dy2(1:10,10, 1)
     .-014,-032,-026,-020,-042,-032,-023,-086,-019,-029, ! dy2(1:10,11, 1)
     .-009,-032,-033,-019,-035,-027,-021,-064,-014,-034, ! dy2(1:10,12, 1)
     .-002,-009,-028,-015,-016,-011,-015,-027,-008,-005, ! dy2(1:10,13, 1)
     . 005,-002,-010, 001,-003, 003,-008,-020, 002,-008, ! dy2(1:10,14, 1)
     . 040, 021, 000, 010, 009, 023, 011, 022, 020, 013, ! dy2(1:10,15, 1)
     . 032, 033, 015, 019, 028, 030, 013, 041, 022, 029, ! dy2(1:10,16, 1)
     . 041, 040, 028, 037, 034, 042, 032, 056, 031, 037, ! dy2(1:10,17, 1)
     . 061, 057, 033, 030, 040, 040, 025, 062, 030, 034, ! dy2(1:10,18, 1)
     . 062, 053, 035, 041, 035, 042, 034, 068, 032, 049, ! dy2(1:10,19, 1)
     . 047, 032, 021, 010, 030, 014, 019, 039, 020, 025, ! dy2(1:10,20, 1)
     . 027, 024, 003, 010, 011,-001, 004, 026,-005, 021, ! dy2(1:10,21, 1)
     . 018, 007, 000, 006, 007,-013,-005, 008,-008, 003, ! dy2(1:10,22, 1)
     . 003, 003,-021, 004,-011,-027,-013,-006,-025, 000, ! dy2(1:10,23, 1)
     .-027,-027,-032,-020,-026,-044,-033,-030,-053,-018, ! dy2(1:10,24, 1)
     .-038,-039,-043,-034,-040,-051,-048,-046,-060,-035, ! dy2(1:10,25, 1)
     .-033,-032,-022,-018,-019,-037,-025,-028,-039,-013, ! dy2(1:10,26, 1)
     .-021,-028,-027,-021,-025,-035,-025,-033,-029,-025, ! dy2(1:10,27, 1)
     .-031,-024,-020,-015,-010,-031,-015,-023,-016,-014, ! dy2(1:10,28, 1)
     .-025,-017,-013,-014,-012,-031,-013,-025,-013,-024, ! dy2(1:10,29, 1)
     .-011,-017, 002, 001, 000,-008,-001,-003, 002,-009, ! dy2(1:10,30, 1)
     .-007,-008, 003, 007, 002,-019, 010, 005, 011, 002, ! dy2(1:10,31, 1)
     .-001, 009, 017, 016, 014,-006, 012, 022, 021, 000, ! dy2(1:10,32, 1)
     . 000, 010, 026, 018, 020, 002, 031, 038, 029, 003, ! dy2(1:10,33, 1)
     .-019,-003, 029, 015, 019,-002, 028, 034, 035, 017, ! dy2(1:10,34, 1)
     . 017, 017, 025, 027, 019,-003, 028, 031, 037, 003, ! dy2(1:10,35, 1)
     .-024,-002, 020, 013, 023, 005, 021, 014, 035, 015, ! dy2(1:10,36, 1)
     .-021,-015, 004,-003, 005, 000, 016,-040, 016,-027, ! dy2(1:10,37, 1)
     .-029,-019,-022,-029,-010,-007, 000,-056, 023,-011, ! dy2(1:10,38, 1)
     . 063, 013, 092, 054, 014, 029, 043, 026, 039, 003, ! dy2(1:10, 1, 2)
     . 036, 028, 090, 068, 019, 027, 053, 023, 034, 004, ! dy2(1:10, 2, 2)
     . 026,-002, 042, 034,-007, 007, 017,-006, 007,-010, ! dy2(1:10, 3, 2)
     .-005,-014, 030, 016,-019,-006, 014,-018, 003,-014, ! dy2(1:10, 4, 2)
     . 005,-023, 006,-001,-027,-017, 000,-033,-011,-017, ! dy2(1:10, 5, 2)
     .-022,-025,-015,-019,-041,-024,-010,-054,-019,-028, ! dy2(1:10, 6, 2)
     .-015,-038,-017,-014,-037,-021,-011,-064,-011,-015, ! dy2(1:10, 7, 2)
     .-014,-040,-025,-019,-047,-028,-018,-068,-025,-033, ! dy2(1:10, 8, 2)
     .-025,-042,-031,-034,-040,-025,-026,-068,-018,-028, ! dy2(1:10, 9, 2)
     .-025,-041,-038,-031,-039,-026,-034,-064,-023,-031, ! dy2(1:10,10, 2)
     .-011,-034,-037,-028,-026,-025,-026,-046,-014,-035, ! dy2(1:10,11, 2)
     . 000,-010,-029,-021,-016,-012,-025,-025,-013,-016, ! dy2(1:10,12, 2)
     . 010,-003,-011,-011, 002, 007,-009, 006, 004,-016, ! dy2(1:10,13, 2)
     . 020, 006,-009,-007, 009, 010,-008, 019, 013, 002, ! dy2(1:10,14, 2)
     . 024, 029, 011, 009, 025, 021, 004, 044, 016, 007, ! dy2(1:10,15, 2)
     . 043, 036, 024, 016, 034, 033, 014, 062, 025, 014, ! dy2(1:10,16, 2)
     . 037, 048, 031, 021, 037, 036, 016, 070, 030, 022, ! dy2(1:10,17, 2)
     . 047, 059, 042, 033, 047, 047, 028, 079, 035, 035, ! dy2(1:10,18, 2)
     . 047, 058, 036, 036, 043, 032, 028, 065, 035, 031, ! dy2(1:10,19, 2)
     . 031, 049, 027, 034, 032, 028, 015, 058, 021, 018, ! dy2(1:10,20, 2)
     . 010, 027, 006, 018, 015, 006, 001, 022, 001, 003, ! dy2(1:10,21, 2)
     .-006, 002,-017, 001,-008,-008,-015,-003,-021, 000, ! dy2(1:10,22, 2)
     .-016,-014,-022,-014,-015,-023,-024,-021,-040,-025, ! dy2(1:10,23, 2)
     .-030,-032,-040,-031,-033,-037,-043,-045,-055,-033, ! dy2(1:10,24, 2)
     .-051,-046,-045,-042,-038,-045,-044,-064,-067,-037, ! dy2(1:10,25, 2)
     .-015,-018,-024,-017,-017,-013,-024,-041,-041,-020, ! dy2(1:10,26, 2)
     .-023,-018,-026,-029,-026,-023,-028,-048,-048,-012, ! dy2(1:10,27, 2)
     .-006,-025,-021,-018,-020,-019,-025,-040,-036,-024, ! dy2(1:10,28, 2)
     .-011,-010,-006,-010,-009,-013,-016,-024,-018,-011, ! dy2(1:10,29, 2)
     .-003,-008,-002,-003,-011,-013,-009,-014,-008,-005, ! dy2(1:10,30, 2)
     .-011,-005, 005,-001, 000,-001, 000, 006, 006,-005, ! dy2(1:10,31, 2)
     . 001, 001, 017, 011, 008, 011, 010, 017, 018, 005, ! dy2(1:10,32, 2)
     . 000, 006, 020, 015, 011, 004, 011, 030, 020, 007, ! dy2(1:10,33, 2)
     . 007, 019, 027, 031, 020, 014, 023, 040, 029, 007, ! dy2(1:10,34, 2)
     .-014, 001, 020, 015, 015, 012, 015, 031, 020, 015, ! dy2(1:10,35, 2)
     . 011, 003, 021, 020, 013, 013, 020, 026, 022, 012, ! dy2(1:10,36, 2)
     .-037,-002, 006, 003, 010, 013, 005, 004, 011, 005, ! dy2(1:10,37, 2)
     .-046,-047, 001,-004,-006, 003,-014,-024,-007,-011, ! dy2(1:10,38, 2)
     . 055, 042, 079, 075, 040, 060, 076, 054, 053, 036, ! dy2(1:10, 1, 3)
     . 073, 032, 079, 071, 035, 051, 058, 058, 045, 036, ! dy2(1:10, 2, 3)
     . 017, 009, 043, 028, 003, 017, 028, 031, 012, 011, ! dy2(1:10, 3, 3)
     . 009,-002, 021, 017,-002, 015, 017, 017, 010, 017, ! dy2(1:10, 4, 3)
     . 014, 000, 010, 007,-010, 017, 016, 008, 010, 010, ! dy2(1:10, 5, 3)
     .-003,-025,-006,-004,-025,-005,-002,-017,-007, 000, ! dy2(1:10, 6, 3)
     .-006,-028,-021,-016,-026,-007,-008,-035,-012,-019, ! dy2(1:10, 7, 3)
     .-015,-024,-022,-015,-019,-002,-009,-035,-002,-006, ! dy2(1:10, 8, 3)
     .-006,-039,-029,-025,-031,-018,-025,-047,-017,-022, ! dy2(1:10, 9, 3)
     .-009,-024,-027,-020,-020,-014,-022,-028,-009,-015, ! dy2(1:10,10, 3)
     .-004,-008,-018,-017,-009,-006,-012,-013, 003,-008, ! dy2(1:10,11, 3)
     .-010,-012,-021,-023,-003,-007,-022,-002,-006,-018, ! dy2(1:10,12, 3)
     . 027, 011,-006, 002, 016, 005,-001, 031, 010, 008, ! dy2(1:10,13, 3)
     . 022, 025, 012, 011, 029, 022, 013, 042, 020, 010, ! dy2(1:10,14, 3)
     . 017, 014, 005, 002, 025, 009,-002, 044, 015, 010, ! dy2(1:10,15, 3)
     . 041, 041, 028, 017, 041, 032, 016, 064, 026, 020, ! dy2(1:10,16, 3)
     . 057, 057, 051, 037, 052, 041, 034, 080, 043, 033, ! dy2(1:10,17, 3)
     . 054, 059, 039, 039, 049, 036, 026, 070, 038, 034, ! dy2(1:10,18, 3)
     . 054, 058, 045, 050, 047, 040, 034, 069, 041, 034, ! dy2(1:10,19, 3)
     . 038, 048, 031, 036, 039, 029, 029, 056, 034, 034, ! dy2(1:10,20, 3)
     . 021, 019, 008, 012, 012, 003, 000, 017, 001, 001, ! dy2(1:10,21, 3)
     .-006,-006,-008,-005, 001,-006,-007,-005,-016,-005, ! dy2(1:10,22, 3)
     .-021,-017,-022,-020,-017,-019,-024,-030,-032,-013, ! dy2(1:10,23, 3)
     .-037,-033,-044,-036,-036,-042,-043,-064,-060,-040, ! dy2(1:10,24, 3)
     .-029,-041,-053,-045,-040,-038,-044,-074,-070,-035, ! dy2(1:10,25, 3)
     .-017,-021,-027,-021,-019,-020,-023,-055,-048,-014, ! dy2(1:10,26, 3)
     .-023,-018,-021,-022,-026,-020,-028,-053,-050,-027, ! dy2(1:10,27, 3)
     .-027,-015,-021,-021,-018,-014,-021,-049,-034,-007, ! dy2(1:10,28, 3)
     .-003,-011,-005,-007,-011,-006,-012,-033,-019,-009, ! dy2(1:10,29, 3)
     .-010,-002,-004,-006,-010,-002,-013,-025,-013,-013, ! dy2(1:10,30, 3)
     .-009, 001, 009, 000, 005, 005,-002, 002, 014, 014, ! dy2(1:10,31, 3)
     . 004, 018, 018, 009, 014, 016, 013, 012, 021, 006, ! dy2(1:10,32, 3)
     . 012, 011, 015, 016, 012, 008, 008, 024, 024, 013, ! dy2(1:10,33, 3)
     . 029, 014, 034, 028, 023, 024, 022, 043, 031, 023, ! dy2(1:10,34, 3)
     . 010, 029, 033, 030, 027, 020, 023, 044, 032, 026, ! dy2(1:10,35, 3)
     .-018, 004, 016, 019, 016, 016, 009, 038, 014, 015, ! dy2(1:10,36, 3)
     .-003, 016, 020, 019, 019, 027, 016, 020, 023, 021, ! dy2(1:10,37, 3)
     .-033, 000, 003, 005, 001, 003,-005, 001, 007, 000, ! dy2(1:10,38, 3)
     . 066, 064, 058, 052, 026, 042, 042, 057, 035, 019, ! dy2(1:10, 1, 4)
     . 047, 036, 060, 050, 031, 037, 044, 057, 035, 037, ! dy2(1:10, 2, 4)
     . 008, 006, 024, 022, 001, 018, 024, 039, 010, 017, ! dy2(1:10, 3, 4)
     . 022, 006, 013, 009,-006, 017, 016, 030, 009, 015, ! dy2(1:10, 4, 4)
     .-002, 000, 003,-001,-005, 015, 018, 020, 003, 012, ! dy2(1:10, 5, 4)
     .-001,-008,-007,-005,-010, 006, 002, 003,-001, 008, ! dy2(1:10, 6, 4)
     . 009,-004, 000, 003,-004, 009, 000,-001, 008, 004, ! dy2(1:10, 7, 4)
     . 026,-007, 001, 004,-008, 013, 010,-002, 005,-004, ! dy2(1:10, 8, 4)
     .-021,-025,-032,-023,-024,-015,-028,-023,-015,-015, ! dy2(1:10, 9, 4)
     .-005,-023,-025,-020,-010,-010,-018,-011,-006,-010, ! dy2(1:10,10, 4)
     . 012, 004, 003,-001, 009, 006,-003, 004, 011, 005, ! dy2(1:10,11, 4)
     . 005, 001,-015,-013, 006,-002,-014, 004, 000,-010, ! dy2(1:10,12, 4)
     . 022, 012, 016, 002, 026, 010, 001, 039, 015, 016, ! dy2(1:10,13, 4)
     . 025, 035, 020, 010, 045, 025, 011, 050, 030, 019, ! dy2(1:10,14, 4)
     . 037, 028, 018, 005, 029, 015, 009, 048, 018, 011, ! dy2(1:10,15, 4)
     . 036, 036, 034, 026, 043, 025, 018, 059, 031, 026, ! dy2(1:10,16, 4)
     . 048, 044, 041, 043, 053, 040, 029, 068, 039, 031, ! dy2(1:10,17, 4)
     . 032, 042, 038, 034, 041, 031, 025, 063, 034, 024, ! dy2(1:10,18, 4)
     . 039, 037, 037, 032, 039, 030, 027, 048, 038, 028, ! dy2(1:10,19, 4)
     . 033, 031, 037, 027, 035, 027, 021, 047, 032, 025, ! dy2(1:10,20, 4)
     .-009, 006, 000,-001, 007, 001,-002, 002,-002,-003, ! dy2(1:10,21, 4)
     .-008,-003,-009,-007,-004,-005,-007,-015,-011,-010, ! dy2(1:10,22, 4)
     .-012,-014,-024,-021,-015,-012,-018,-038,-030,-015, ! dy2(1:10,23, 4)
     .-037,-046,-038,-038,-039,-036,-037,-065,-054,-031, ! dy2(1:10,24, 4)
     .-035,-030,-039,-040,-041,-030,-036,-075,-055,-033, ! dy2(1:10,25, 4)
     .-015,-008,-013,-013,-014,-008,-014,-053,-031,-011, ! dy2(1:10,26, 4)
     .-015,-012,-022,-015,-025,-018,-021,-059,-041,-017, ! dy2(1:10,27, 4)
     . 000,-012,-001,-005,-014,-005,-018,-048,-024,-017, ! dy2(1:10,28, 4)
     . 001, 004, 000, 000,-001, 005,-007,-034,-016, 002, ! dy2(1:10,29, 4)
     .-015,-002,-003,-006,-009, 009,-015,-029,-008,-005, ! dy2(1:10,30, 4)
     . 007, 013, 011, 006, 005, 002,-003,-010, 012, 000, ! dy2(1:10,31, 4)
     . 026, 029, 028, 024, 015, 027, 008, 010, 029, 026, ! dy2(1:10,32, 4)
     .-001, 008, 014, 014, 013, 015, 010, 018, 022, 018, ! dy2(1:10,33, 4)
     . 017, 035, 021, 020, 021, 019, 017, 035, 036, 019, ! dy2(1:10,34, 4)
     . 005, 018, 026, 027, 019, 023, 016, 041, 023, 032, ! dy2(1:10,35, 4)
     .-007, 012, 015, 013, 010, 013, 004, 032, 017, 017, ! dy2(1:10,36, 4)
     .-006,-001, 008, 008, 009, 017, 000, 025, 007, 010, ! dy2(1:10,37, 4)
     .-031,-020,-008,-008,-007, 005,-019, 003,-002, 005, ! dy2(1:10,38, 4)
     . 022, 022, 029, 025, 003, 014, 036, 061, 010, 039, ! dy2(1:10, 1, 5)
     . 029, 026, 039, 041, 021, 038, 042, 069, 020, 031, ! dy2(1:10, 2, 5)
     . 025, 010, 009, 009, 000, 019, 016, 048, 003, 025, ! dy2(1:10, 3, 5)
     . 011,-001, 003, 005, 002, 022, 014, 041, 004, 021, ! dy2(1:10, 4, 5)
     . 020, 001, 006,-001,-001, 014, 007, 025, 004, 012, ! dy2(1:10, 5, 5)
     . 017,-012, 004, 002, 000, 015, 011, 019, 006, 016, ! dy2(1:10, 6, 5)
     . 011, 005,-005,-002, 005, 015, 007, 008, 009, 018, ! dy2(1:10, 7, 5)
     . 008, 008, 000, 009, 005, 014, 008, 014, 011, 024, ! dy2(1:10, 8, 5)
     . 009,-016,-013,-010,-007,-001,-009,-006,-005,-001, ! dy2(1:10, 9, 5)
     . 001,-008,-017,-011,-010,-009,-018,-005,-004,-009, ! dy2(1:10,10, 5)
     . 009, 005,-004,-002, 013, 001,-002, 012, 013, 008, ! dy2(1:10,11, 5)
     . 009, 007,-001,-008, 009, 000,-006, 015, 009, 002, ! dy2(1:10,12, 5)
     . 038, 034, 014, 012, 033, 015, 013, 039, 024, 018, ! dy2(1:10,13, 5)
     . 043, 032, 032, 024, 045, 030, 026, 053, 037, 023, ! dy2(1:10,14, 5)
     . 033, 029, 027, 021, 036, 013, 013, 050, 025, 018, ! dy2(1:10,15, 5)
     . 046, 035, 030, 036, 044, 023, 023, 057, 036, 019, ! dy2(1:10,16, 5)
     . 040, 035, 040, 036, 045, 029, 030, 060, 038, 028, ! dy2(1:10,17, 5)
     . 021, 021, 021, 019, 026, 015, 013, 035, 024, 015, ! dy2(1:10,18, 5)
     . 014, 019, 021, 012, 029, 020, 019, 031, 029, 010, ! dy2(1:10,19, 5)
     . 016, 017, 010, 005, 018, 011, 011, 024, 019, 009, ! dy2(1:10,20, 5)
     .-009,-011,-003,-007,-006,-007,-003,-008,-005,-008, ! dy2(1:10,21, 5)
     .-022,-015,-012,-017,-015,-013,-011,-028,-019,-019, ! dy2(1:10,22, 5)
     .-020,-022,-017,-021,-015,-014,-014,-040,-024,-016, ! dy2(1:10,23, 5)
     .-033,-032,-027,-033,-032,-024,-026,-059,-043,-029, ! dy2(1:10,24, 5)
     .-025,-024,-030,-026,-037,-025,-029,-070,-043,-031, ! dy2(1:10,25, 5)
     .-004,-009,-002,-002,-015, 000,-003,-038,-016,-002, ! dy2(1:10,26, 5)
     .-020,-020,-018,-006,-021,-013,-024,-055,-033,-010, ! dy2(1:10,27, 5)
     .-002,-003,-001, 002,-010, 004,-009,-043,-019,-002, ! dy2(1:10,28, 5)
     .-009,-001,-002,-003,-003, 008,-018,-040,-005, 010, ! dy2(1:10,29, 5)
     .-014,-009,-009,-012,-008, 001,-020,-037,-012,-008, ! dy2(1:10,30, 5)
     . 020, 023, 016, 019, 013, 028, 003,-002, 019, 017, ! dy2(1:10,31, 5)
     . 035, 032, 021, 026, 020, 039, 012, 008, 032, 033, ! dy2(1:10,32, 5)
     . 004, 024, 016, 016, 016, 017, 008, 016, 019, 022, ! dy2(1:10,33, 5)
     . 014, 017, 017, 020, 016, 021, 014, 025, 022, 030, ! dy2(1:10,34, 5)
     .-004, 011, 009, 010, 009, 011, 007, 034, 018, 015, ! dy2(1:10,35, 5)
     . 012, 001, 000, 008, 008, 007, 002, 032, 007, 019, ! dy2(1:10,36, 5)
     .-034,-005, 000,-002,-002, 006,-002, 025, 000, 004, ! dy2(1:10,37, 5)
     .-042,-022,-013,-026,-014, 000,-025, 009,-027,-013, ! dy2(1:10,38, 5)
     . 021,-003, 000, 008,-001, 016, 010, 049,-022, 022, ! dy2(1:10, 1, 6)
     . 026, 014, 008, 004,-003, 009, 018, 055,-001, 022, ! dy2(1:10, 2, 6)
     . 019, 000, 008, 003, 001, 020, 008, 043, 003, 019, ! dy2(1:10, 3, 6)
     . 015, 001,-001, 000,-003, 009, 009, 033,-004, 016, ! dy2(1:10, 4, 6)
     . 012, 002, 001, 004, 002, 020, 014, 029, 009, 018, ! dy2(1:10, 5, 6)
     . 006, 000,-007,-009,-009, 004,-004, 012,-006, 010, ! dy2(1:10, 6, 6)
     . 029, 006, 004, 011, 007, 019, 013, 019, 016, 025, ! dy2(1:10, 7, 6)
     . 018, 004, 003, 008, 010, 022, 010, 023, 019, 024, ! dy2(1:10, 8, 6)
     . 001,-008,-013,-007,-001, 000,-008, 005,-002, 005, ! dy2(1:10, 9, 6)
     . 008,-007,-007, 001, 003, 000,-003, 008, 001, 001, ! dy2(1:10,10, 6)
     . 028, 011, 010, 006, 017, 011, 003, 022, 015, 012, ! dy2(1:10,11, 6)
     . 004, 001,-001,-002, 009, 000,-002, 020, 012, 005, ! dy2(1:10,12, 6)
     . 029, 020, 022, 011, 026, 008, 012, 035, 021, 009, ! dy2(1:10,13, 6)
     . 037, 029, 027, 032, 033, 022, 020, 045, 033, 029, ! dy2(1:10,14, 6)
     . 025, 030, 022, 023, 026, 016, 016, 040, 023, 004, ! dy2(1:10,15, 6)
     . 031, 022, 025, 030, 031, 025, 019, 045, 029, 013, ! dy2(1:10,16, 6)
     . 016, 025, 029, 024, 029, 016, 025, 047, 034, 016, ! dy2(1:10,17, 6)
     . 022, 015, 022, 019, 027, 019, 024, 035, 027, 012, ! dy2(1:10,18, 6)
     . 006, 012, 017, 008, 018, 010, 017, 019, 025, 017, ! dy2(1:10,19, 6)
     . 000, 018, 009, 006, 004, 006, 018, 012, 017, 000, ! dy2(1:10,20, 6)
     .-019,-021,-010,-015,-011,-014,-003,-022,-012,-016, ! dy2(1:10,21, 6)
     .-010,-010,-006,-011,-015,-010,-002,-026,-011,-016, ! dy2(1:10,22, 6)
     .-019,-019,-015,-017,-023,-017,-007,-043,-022,-022, ! dy2(1:10,23, 6)
     .-026,-026,-025,-023,-030,-026,-020,-054,-033,-030, ! dy2(1:10,24, 6)
     .-025,-029,-023,-018,-030,-017,-023,-052,-032,-025, ! dy2(1:10,25, 6)
     .-008,-005, 003, 002,-007, 003, 000,-033,-009,-003, ! dy2(1:10,26, 6)
     .-012,-010,-017,-017,-020,-010,-017,-046,-026,-022, ! dy2(1:10,27, 6)
     .-014,-010,-007,-013,-011, 000,-019,-039,-015,-003, ! dy2(1:10,28, 6)
     . 011, 014, 010, 000, 000, 016,-003,-026, 001, 004, ! dy2(1:10,29, 6)
     . 011, 004,-003, 000,-005, 011,-017,-027,-008, 008, ! dy2(1:10,30, 6)
     .-003, 028, 015, 018, 018, 020, 000,-005, 015, 017, ! dy2(1:10,31, 6)
     . 029, 032, 029, 029, 032, 036, 020, 014, 025, 032, ! dy2(1:10,32, 6)
     . 011, 020, 008, 013, 012, 009, 000, 006, 010, 022, ! dy2(1:10,33, 6)
     . 011, 019, 004, 008, 013, 009, 004, 023, 014, 016, ! dy2(1:10,34, 6)
     .-004, 017, 009, 007, 012, 004, 008, 029, 012, 029, ! dy2(1:10,35, 6)
     .-016, 008,-002,-004,-001,-003,-005, 025, 004, 009, ! dy2(1:10,36, 6)
     .-024,-018,-015,-012,-011,-003,-011, 021,-017, 000, ! dy2(1:10,37, 6)
     . 008,-022,-017,-027,-026,-034,-034, 003,-032,-020, ! dy2(1:10,38, 6)
     . 019,-010,-009,-015,-024,-007, 005, 022,-021, 002, ! dy2(1:10, 1, 7)
     .-005,-009,-013,-005,-025,-001,-004, 030,-021, 000, ! dy2(1:10, 2, 7)
     . 013,-003,-004,-010,-004, 012, 008, 023,-009, 016, ! dy2(1:10, 3, 7)
     . 004, 000,-003,-006,-001, 011, 002, 022, 004, 014, ! dy2(1:10, 4, 7)
     . 020, 010, 005, 000, 002, 018, 007, 022, 003, 012, ! dy2(1:10, 5, 7)
     . 001,-009,-012,-006,-006, 007,-004, 009, 000, 005, ! dy2(1:10, 6, 7)
     . 016, 007,-001, 000, 000, 012,-002, 012, 002, 015, ! dy2(1:10, 7, 7)
     . 028, 017, 000, 004, 005, 017, 002, 018, 008, 016, ! dy2(1:10, 8, 7)
     . 002,-016,-012, 004, 000, 003,-006, 007, 001, 005, ! dy2(1:10, 9, 7)
     . 009,-004,-003,-003, 002, 000,-007, 009, 006, 002, ! dy2(1:10,10, 7)
     . 022, 016, 011, 012, 015, 012, 009, 029, 020, 015, ! dy2(1:10,11, 7)
     . 019,-002, 000, 012, 004, 000,-003, 024, 009, 004, ! dy2(1:10,12, 7)
     . 001, 019, 014, 025, 017, 009, 003, 031, 020, 012, ! dy2(1:10,13, 7)
     . 022, 024, 021, 037, 026, 017, 015, 039, 032, 016, ! dy2(1:10,14, 7)
     . 022, 015, 012, 021, 021, 009, 010, 032, 023, 007, ! dy2(1:10,15, 7)
     . 020, 027, 021, 017, 021, 011, 013, 027, 026, 013, ! dy2(1:10,16, 7)
     . 018, 029, 020, 021, 027, 012, 021, 035, 032, 017, ! dy2(1:10,17, 7)
     . 007, 005, 018, 007, 016, 002, 017, 017, 020, 002, ! dy2(1:10,18, 7)
     . 000, 017, 014, 002, 005, 007, 022, 006, 018,-005, ! dy2(1:10,19, 7)
     .-018,-007, 004,-007, 000,-009, 002,-005, 006,-009, ! dy2(1:10,20, 7)
     .-025,-016,-011,-010,-020,-015,-005,-022,-009,-019, ! dy2(1:10,21, 7)
     .-022,-024,-018,-011,-024,-023,-006,-033,-016,-021, ! dy2(1:10,22, 7)
     .-022,-019,-020,-018,-027,-020,-008,-040,-021,-029, ! dy2(1:10,23, 7)
     .-037,-025,-024,-028,-035,-027,-020,-048,-026,-031, ! dy2(1:10,24, 7)
     .-018,-031,-021,-036,-031,-026,-025,-053,-029,-028, ! dy2(1:10,25, 7)
     .-013,-001,-003,-009,-011,-008,-005,-029,-009,-007, ! dy2(1:10,26, 7)
     .-019,-018,-006,-012,-018,-006,-015,-038,-014,-017, ! dy2(1:10,27, 7)
     . 000, 000,-005,-003,-009,-001,-013,-029,-010,-008, ! dy2(1:10,28, 7)
     . 006, 008, 002, 008, 005, 017,-005,-019,-001, 015, ! dy2(1:10,29, 7)
     .-006, 009, 000, 001, 005, 015,-011,-024,-008, 000, ! dy2(1:10,30, 7)
     . 011, 016, 004, 013, 018, 020, 000,-007, 001, 017, ! dy2(1:10,31, 7)
     . 030, 044, 022, 030, 033, 036, 016, 009, 019, 040, ! dy2(1:10,32, 7)
     .-009, 012, 000, 001, 011, 007,-006,-001,-006, 009, ! dy2(1:10,33, 7)
     . 007, 012, 004, 009, 015, 010, 003, 015, 000, 020, ! dy2(1:10,34, 7)
     .-006, 014,-001, 004, 008, 002, 006, 024, 008, 020, ! dy2(1:10,35, 7)
     .-013,-012,-020,-009,-004,-014,-006, 018,-013, 009, ! dy2(1:10,36, 7)
     .-069,-017,-018,-029,-027,-022,-025, 013,-018,-006, ! dy2(1:10,37, 7)
     .-054,-033,-050,-044,-052,-036,-037, 000,-045,-009, ! dy2(1:10,38, 7)
     .-046,-040,-016,-051,-025,-025,-024, 002,-038,-001, ! dy2(1:10, 1, 8)
     . 005,-018,-023,-023,-022,-011,-011, 005,-026,-006, ! dy2(1:10, 2, 8)
     . 005,-013,-011,-012,-017, 002,-005, 004,-020, 006, ! dy2(1:10, 3, 8)
     . 013,-004,-005,-008,-004, 007, 001, 004,-007, 020, ! dy2(1:10, 4, 8)
     . 015, 004,-008,-005,-007, 011, 000, 004,-003, 023, ! dy2(1:10, 5, 8)
     . 012, 000,-006,-004, 003, 004,-005, 000, 000, 010, ! dy2(1:10, 6, 8)
     . 017, 006,-005, 007, 000, 018, 003, 011, 006, 007, ! dy2(1:10, 7, 8)
     . 012, 003,-008, 003,-005, 009,-004, 011, 008, 001, ! dy2(1:10, 8, 8)
     . 000,-005,-011,-005,-002,-002,-008, 003, 000, 005, ! dy2(1:10, 9, 8)
     . 013,-004,-010,-001,-006,-008,-009, 006, 000,-004, ! dy2(1:10,10, 8)
     . 005, 000, 006, 020, 005, 001, 002, 022, 010, 006, ! dy2(1:10,11, 8)
     . 024, 003, 013, 027, 006, 004, 004, 027, 015, 008, ! dy2(1:10,12, 8)
     . 006,-004, 011, 029, 010,-006, 010, 028, 019, 008, ! dy2(1:10,13, 8)
     . 032, 018, 024, 022, 019, 012, 018, 039, 029, 019, ! dy2(1:10,14, 8)
     . 002, 005, 010, 008, 008,-005, 007, 018, 016, 000, ! dy2(1:10,15, 8)
     . 017, 017, 016, 014, 018, 005, 016, 021, 026, 005, ! dy2(1:10,16, 8)
     . 011, 017, 020, 016, 017, 008, 022, 020, 024, 000, ! dy2(1:10,17, 8)
     . 000, 016, 013, 007, 002, 000, 015, 004, 018,-001, ! dy2(1:10,18, 8)
     .-003, 001, 011, 008, 009,-007, 017, 004, 013, 003, ! dy2(1:10,19, 8)
     .-007,-002, 000, 003,-005,-009, 011,-008,-001,-015, ! dy2(1:10,20, 8)
     .-037,-026,-017,-018,-028,-021,-010,-030,-015,-034, ! dy2(1:10,21, 8)
     .-028,-023,-019,-024,-026,-020,-010,-033,-017,-030, ! dy2(1:10,22, 8)
     .-027,-024,-018,-026,-031,-026,-016,-031,-014,-028, ! dy2(1:10,23, 8)
     .-037,-038,-018,-035,-033,-027,-023,-041,-027,-039, ! dy2(1:10,24, 8)
     .-028,-029,-031,-028,-036,-028,-026,-044,-025,-031, ! dy2(1:10,25, 8)
     .-002,-017,-006,-009,-010, 001,-003,-020,-003,-009, ! dy2(1:10,26, 8)
     . 000, 002,-004, 001,-009,-004,-008,-019,-012,-006, ! dy2(1:10,27, 8)
     .-010, 003, 005, 005, 004, 009,-004,-015,-004, 003, ! dy2(1:10,28, 8)
     . 019, 020, 008, 013, 011, 014,-001,-007,-005, 005, ! dy2(1:10,29, 8)
     . 001, 007,-002, 001, 015, 015,-006,-011,-014, 007, ! dy2(1:10,30, 8)
     .-007, 026, 008, 012, 016, 024, 001, 002,-005, 017, ! dy2(1:10,31, 8)
     . 014, 029, 017, 019, 033, 036, 005, 012, 005, 027, ! dy2(1:10,32, 8)
     . 014, 018,-003, 005, 016, 009,-001, 008,-003, 021, ! dy2(1:10,33, 8)
     .-001, 018,-005, 006, 015, 009, 003, 013, 000, 020, ! dy2(1:10,34, 8)
     .-005, 006,-008,-002, 009, 005, 001, 021,-004, 018, ! dy2(1:10,35, 8)
     .-037,-018,-019,-022,-009,-009,-009, 010,-018, 003, ! dy2(1:10,36, 8)
     .-034,-037,-035,-044,-035,-046,-028,-005,-038,-026, ! dy2(1:10,37, 8)
     .-044,-073,-045,-057,-054,-051,-042,-022,-064,-053, ! dy2(1:10,38, 8)
     .-001,-029,-031,-030,-036,-025,-015,-019,-032,-017, ! dy2(1:10, 1, 9)
     .-011,-034,-020,-032,-030,-015,-020,-015,-026,-015, ! dy2(1:10, 2, 9)
     . 000,-017,-027,-033,-022,-010,-016,-023,-030, 007, ! dy2(1:10, 3, 9)
     . 008,-011,-016,-016,-012, 003,-002,-014,-013, 004, ! dy2(1:10, 4, 9)
     . 022, 002,-010,-006,-005, 007, 000,-005,-006, 015, ! dy2(1:10, 5, 9)
     . 003,-007,-017,-009,-014, 003,-014,-015,-011,-002, ! dy2(1:10, 6, 9)
     . 020, 010,-008,-007, 003, 000,-003, 000, 006, 024, ! dy2(1:10, 7, 9)
     . 017, 006,-005,-007, 005, 017,-003, 003, 003, 016, ! dy2(1:10, 8, 9)
     .-005,-015,-016,-008,-019,-009,-018,-004,-009,-007, ! dy2(1:10, 9, 9)
     .-011,-026,-013, 005,-015,-010,-016, 004,-006,-012, ! dy2(1:10,10, 9)
     . 012, 002, 002, 028, 001, 004,-003, 021, 006, 002, ! dy2(1:10,11, 9)
     .-009,-009,-003, 022,-007,-014,-002, 015, 003, 002, ! dy2(1:10,12, 9)
     . 010, 000, 007, 012, 003,-004, 005, 020, 013, 005, ! dy2(1:10,13, 9)
     . 007, 005, 015, 026, 008, 000, 012, 025, 020, 018, ! dy2(1:10,14, 9)
     . 003, 004, 010, 004, 003,-008, 013, 015, 017,-006, ! dy2(1:10,15, 9)
     .-002, 009, 011, 006, 001,-011, 010, 009, 013, 001, ! dy2(1:10,16, 9)
     . 000, 008, 013, 014, 004,-002, 019, 009, 021,-002, ! dy2(1:10,17, 9)
     .-001, 006, 012, 018, 001,-002, 016, 000, 010,-003, ! dy2(1:10,18, 9)
     .-014, 004, 011, 018,-003, 003, 017,-001, 010,-008, ! dy2(1:10,19, 9)
     .-026,-003, 000,-006,-017,-011, 005,-012, 000,-019, ! dy2(1:10,20, 9)
     .-012,-016,-005,-015,-018,-013,-005,-015,-006,-020, ! dy2(1:10,21, 9)
     .-029,-028,-007,-032,-030,-020,-013,-024,-013,-031, ! dy2(1:10,22, 9)
     .-030,-032,-024,-031,-030,-026,-022,-029,-020,-031, ! dy2(1:10,23, 9)
     .-032,-030,-022,-027,-039,-029,-027,-030,-023,-041, ! dy2(1:10,24, 9)
     .-032,-025,-021,-020,-022,-018,-022,-020,-017,-026, ! dy2(1:10,25, 9)
     .-008, 005, 001, 004,-005, 002,-004,-001, 000,-010, ! dy2(1:10,26, 9)
     .-006,-003,-010,-001,-001, 000,-012,-007,-007,-004, ! dy2(1:10,27, 9)
     . 009, 014, 001, 006, 007, 003,-002,-004,-006, 000, ! dy2(1:10,28, 9)
     . 022, 017, 007, 013, 024, 025, 005, 006,-002, 016, ! dy2(1:10,29, 9)
     . 002, 017, 000, 010, 015, 021,-004, 000,-013, 006, ! dy2(1:10,30, 9)
     . 011, 021, 003, 002, 026, 022, 000, 004,-005, 021, ! dy2(1:10,31, 9)
     . 024, 035, 013, 018, 033, 038, 010, 019, 006, 028, ! dy2(1:10,32, 9)
     . 000, 012,-008, 000, 019, 017,-001, 008,-007, 017, ! dy2(1:10,33, 9)
     . 014, 011, 004, 002, 021, 019, 009, 022, 005, 026, ! dy2(1:10,34, 9)
     .-009, 000,-007, 003, 007, 000, 011, 024, 003, 024, ! dy2(1:10,35, 9)
     .-025,-023,-027,-020,-010,-021,-010, 002,-021, 003, ! dy2(1:10,36, 9)
     .-069,-028,-033,-038,-031,-025,-013,-001,-027,-016, ! dy2(1:10,37, 9)
     .-044,-031,-048,-064,-052,-053,-041,-027,-061,-036, ! dy2(1:10,38, 9)
     .-023,-030,-023,-041,-039,-035,-035,-045,-033,-006, ! dy2(1:10, 1,10)
     . 006,-021,-023,-023,-028,-033,-018,-048,-026,-004, ! dy2(1:10, 2,10)
     . 000,-033,-035,-026,-038,-027,-024,-047,-034,-011, ! dy2(1:10, 3,10)
     . 008,-014,-025,-019,-026,-015,-025,-041,-036, 006, ! dy2(1:10, 4,10)
     . 007, 001,-009,-003, 000,-005,-006,-020,-008, 013, ! dy2(1:10, 5,10)
     . 000, 007,-018,-020,-010,-001,-010,-012,-013, 011, ! dy2(1:10, 6,10)
     . 000,-009,-010, 000,-010, 008,-008,-010,-002, 001, ! dy2(1:10, 7,10)
     . 004, 005,-011,-001,-004,-001,-007,-008,-005,-002, ! dy2(1:10, 8,10)
     .-018,-022,-028, 000,-024,-020,-032,-025,-028,-015, ! dy2(1:10, 9,10)
     .-006,-013,-024, 007,-027,-018,-028,-024,-024,-014, ! dy2(1:10,10,10)
     .-015,-025,-016, 007,-019,-019,-010,-008,-009,-009, ! dy2(1:10,11,10)
     .-012,-029,-016,-004,-020,-021,-014,-010,-014,-018, ! dy2(1:10,12,10)
     .-018,-026,-007,-005,-014,-021,-008,-007,-003,-003, ! dy2(1:10,13,10)
     . 003,-010,-002,-007,-007,-014, 002, 001, 003,-005, ! dy2(1:10,14,10)
     .-015,-020, 002, 009,-014,-031,-001,-008,-002,-004, ! dy2(1:10,15,10)
     . 000,-008, 015, 022,-002,-007, 016, 004, 012,-001, ! dy2(1:10,16,10)
     .-006, 008, 022, 038, 002,-009, 021, 008, 011, 003, ! dy2(1:10,17,10)
     . 005, 002, 014, 015,-011, 000, 011, 000, 011,-009, ! dy2(1:10,18,10)
     . 000, 007, 022,-006, 002,-003, 021, 006, 013,-001, ! dy2(1:10,19,10)
     .-015,-011, 013,-021,-011,-012, 005,-009,-004,-023, ! dy2(1:10,20,10)
     .-025,-019,-006,-025,-025,-024,-010,-011,-011,-026, ! dy2(1:10,21,10)
     .-028,-013,-004,-018,-023,-018,-012,-012,-012,-023, ! dy2(1:10,22,10)
     .-007,-009,-002,-009,-027,-018,-008,-007,-007,-032, ! dy2(1:10,23,10)
     .-024,-019,-013,-011,-019,-019,-021,-008,-018,-020, ! dy2(1:10,24,10)
     .-023,-018,-020,-016,-028,-026,-026, 001,-016,-028, ! dy2(1:10,25,10)
     .-001,-009, 002, 004,-002,-005,-002, 020, 004,-004, ! dy2(1:10,26,10)
     . 000,-005,-017,-002,-001,-004,-016, 004,-021,-006, ! dy2(1:10,27,10)
     . 006, 005,-012,-003, 013, 000,-012, 009,-012,-001, ! dy2(1:10,28,10)
     . 017, 015,-002, 007, 019, 002,-008, 010,-008, 007, ! dy2(1:10,29,10)
     . 010, 011,-005, 001, 022, 021,-005, 009,-006, 024, ! dy2(1:10,30,10)
     . 018, 020, 003, 007, 030, 035, 004, 015, 004, 021, ! dy2(1:10,31,10)
     . 044, 023, 017, 016, 046, 054, 025, 027, 025, 039, ! dy2(1:10,32,10)
     . 004, 026,-001, 001, 025, 026, 005, 014, 011, 028, ! dy2(1:10,33,10)
     . 006, 006,-006,-003, 016, 012, 013, 009, 011, 032, ! dy2(1:10,34,10)
     . 000, 009,-007,-006, 009, 011, 008, 009, 001, 008, ! dy2(1:10,35,10)
     .-006,-013,-014,-014,-003,-009, 002, 002,-014, 005, ! dy2(1:10,36,10)
     .-049,-043,-031,-036,-032,-033,-016,-036,-043,-025, ! dy2(1:10,37,10)
     .-033,-085,-067,-053,-052,-044,-026,-057,-054,-045, ! dy2(1:10,38,10)
     .-023,-008, 026, 023, 006, 013,-006,-006, 000, 014, ! dy2(1:10, 1,11)
     . 043, 016, 025, 009, 014, 007, 018,-010, 016, 023, ! dy2(1:10, 2,11)
     . 035, 009, 000, 009, 005, 009, 016,-009, 008, 018, ! dy2(1:10, 3,11)
     . 051, 021, 016, 022, 021, 034, 017, 000, 008, 047, ! dy2(1:10, 4,11)
     . 046, 036, 017, 022, 025, 025, 026, 013, 024, 039, ! dy2(1:10, 5,11)
     . 048, 030, 019, 022, 026, 035, 026, 019, 022, 041, ! dy2(1:10, 6,11)
     . 049, 049, 023, 028, 024, 035, 028, 031, 037, 039, ! dy2(1:10, 7,11)
     . 013, 024, 007, 019, 013, 017, 005, 016, 014, 023, ! dy2(1:10, 8,11)
     . 021,-010,-002, 026,-001, 006,-002,-007,-001, 006, ! dy2(1:10, 9,11)
     .-005, 001,-005, 020,-003, 001,-007,-002,-004, 000, ! dy2(1:10,10,11)
     . 027, 001, 011, 025, 006, 012, 007, 010, 009, 005, ! dy2(1:10,11,11)
     . 003,-006,-001, 016, 000, 003, 006, 006, 004, 008, ! dy2(1:10,12,11)
     . 010, 002, 012, 011, 005,-003, 015, 015, 018, 009, ! dy2(1:10,13,11)
     . 003, 000, 018, 016, 013, 006, 019, 017, 021, 003, ! dy2(1:10,14,11)
     . 019, 007, 025, 036, 010, 001, 026, 013, 024, 011, ! dy2(1:10,15,11)
     . 025, 034, 054, 068, 031, 019, 046, 031, 046, 023, ! dy2(1:10,16,11)
     . 028, 029, 041, 053, 022, 019, 038, 027, 030, 018, ! dy2(1:10,17,11)
     . 010, 020, 032, 030, 011, 011, 026, 017, 030, 008, ! dy2(1:10,18,11)
     . 016, 017, 029,-004, 008, 018, 034, 015, 019, 002, ! dy2(1:10,19,11)
     . 011, 015, 025,-003, 001, 003, 019, 004, 019,-005, ! dy2(1:10,20,11)
     . 021, 010, 025, 011, 002, 012, 015, 011, 016, 000, ! dy2(1:10,21,11)
     . 007, 019, 020, 013, 000, 012, 015, 013, 014, 008, ! dy2(1:10,22,11)
     . 000, 005, 018, 015, 002, 009, 010, 018, 016,-005, ! dy2(1:10,23,11)
     .-004, 010, 020, 012,-001, 007, 004, 021, 011,-002, ! dy2(1:10,24,11)
     .-022,-008,-006, 000, 000,-006,-006, 019,-008,-011, ! dy2(1:10,25,11)
     .-009, 005, 003,-004, 006, 001,-001, 018, 000,-003, ! dy2(1:10,26,11)
     . 000, 005,-012, 003, 003,-001,-016, 015,-009,-002, ! dy2(1:10,27,11)
     . 019, 021, 003, 014, 022, 014, 004, 021, 003, 012, ! dy2(1:10,28,11)
     . 017, 022, 011, 022, 037, 024, 009, 030, 009, 032, ! dy2(1:10,29,11)
     . 033, 021, 003, 011, 035, 038, 005, 019, 003, 026, ! dy2(1:10,30,11)
     . 024, 038, 018, 017, 037, 045, 022, 032, 015, 048, ! dy2(1:10,31,11)
     . 021, 023, 007, 015, 033, 033, 014, 023, 008, 032, ! dy2(1:10,32,11)
     . 001, 006,-007,-004, 024, 035, 007, 010, 010, 021, ! dy2(1:10,33,11)
     .-001, 012,-004,-001, 019, 014, 006, 009, 008, 020, ! dy2(1:10,34,11)
     .-003, 003,-013,-002, 005, 016, 009, 010,-004, 012, ! dy2(1:10,35,11)
     .-015,-036,-023,-031,-011,-017,-013,-014,-015, 002, ! dy2(1:10,36,11)
     .-058,-036,-036,-036,-042,-039,-018,-030,-048,-028, ! dy2(1:10,37,11)
     .-079,-061,-038,-045,-041,-048,-023,-037,-056,-035, ! dy2(1:10,38,11)
     . 004,-006, 011,-003,-008,-004, 014,-006, 005, 006, ! dy2(1:10, 1,12)
     . 016,-013, 000, 003,-011,-025, 008,-012, 006, 012, ! dy2(1:10, 2,12)
     . 045, 000, 008, 003, 000, 005, 010,-018, 007, 022, ! dy2(1:10, 3,12)
     . 020, 022, 007, 004, 011, 003, 015,-008, 018, 022, ! dy2(1:10, 4,12)
     . 047, 025, 012, 006, 020, 019, 018, 003, 017, 029, ! dy2(1:10, 5,12)
     . 032, 025, 019, 015, 016, 019, 013, 004, 015, 028, ! dy2(1:10, 6,12)
     . 056, 033, 020, 043, 029, 028, 025, 018, 025, 030, ! dy2(1:10, 7,12)
     . 036, 011, 003, 029, 006, 014, 008, 003, 007, 019, ! dy2(1:10, 8,12)
     .-011, 004,-008, 006,-009,-002,-010,-014,-011,-003, ! dy2(1:10, 9,12)
     . 000,-011,-006,-003,-006, 000,-006,-016,-009,-001, ! dy2(1:10,10,12)
     .-011,-002,-001, 005, 001, 006, 007,-008,-005, 004, ! dy2(1:10,11,12)
     .-006,-013,-007,-005,-009,-004,-005,-017,-009,-004, ! dy2(1:10,12,12)
     .-020,-015,-004, 013,-007, 001, 005,-004,-002,-005, ! dy2(1:10,13,12)
     .-010,-014, 001, 020, 000,-004, 008,-006, 001, 002, ! dy2(1:10,14,12)
     .-011,-005, 004, 025,-011,-002, 000,-003, 000,-009, ! dy2(1:10,15,12)
     .-006,-013, 005, 004,-012,-005, 006,-003,-003,-011, ! dy2(1:10,16,12)
     . 004, 002, 021, 002,-002, 001, 015, 005, 014, 000, ! dy2(1:10,17,12)
     . 010, 008, 024,-004, 000, 004, 013, 008, 006, 000, ! dy2(1:10,18,12)
     . 000, 002, 023, 015, 003,-001, 011, 009, 016,-005, ! dy2(1:10,19,12)
     .-011, 001, 018, 004,-006, 000, 006, 003, 002,-011, ! dy2(1:10,20,12)
     . 000,-006, 007, 001,-013,-016,-001, 000,-008,-011, ! dy2(1:10,21,12)
     . 006, 004, 012, 008,-004,-006,-001, 006, 001,-004, ! dy2(1:10,22,12)
     .-002, 003, 008, 007,-005,-008,-001, 012,-001,-007, ! dy2(1:10,23,12)
     . 004, 002, 003, 004, 000,-012,-007, 011,-005,-009, ! dy2(1:10,24,12)
     . 014, 011, 008, 003, 015, 006,-002, 029, 000, 006, ! dy2(1:10,25,12)
     . 021, 012, 004, 010, 013, 003,-003, 025, 001,-001, ! dy2(1:10,26,12)
     .-010, 004,-004, 001, 016, 006,-002, 024,-001, 002, ! dy2(1:10,27,12)
     . 010, 014,-005, 000, 024, 011,-004, 024, 003, 010, ! dy2(1:10,28,12)
     . 034, 027, 002, 013, 036, 030, 010, 036, 020, 032, ! dy2(1:10,29,12)
     . 007, 019, 003, 002, 028, 023, 005, 027, 012, 024, ! dy2(1:10,30,12)
     . 029, 026, 002, 015, 042, 044, 024, 031, 030, 034, ! dy2(1:10,31,12)
     . 014, 017,-005, 000, 030, 044, 016, 020, 027, 025, ! dy2(1:10,32,12)
     .-009,-004,-013,-016, 010, 013, 002, 005, 005, 019, ! dy2(1:10,33,12)
     .-006,-003,-012,-007, 012, 020, 016, 010, 014, 023, ! dy2(1:10,34,12)
     .-024,-015,-015,-018, 003, 002, 005,-006, 007, 011, ! dy2(1:10,35,12)
     .-032,-007,-017,-019,-018,-011,-002,-018,-019,-013, ! dy2(1:10,36,12)
     .-068,-050,-029,-041,-030,-026,-013,-049,-035,-021, ! dy2(1:10,37,12)
     .-030,-023,-023,-035,-035,-020,-018,-059,-036,-035, ! dy2(1:10,38,12)
     .-006,-006, 015,-019,-010,-024, 013,-008, 033, 035, ! dy2(1:10, 1,13)
     .-022,-010, 016,-001,-008, 005, 013, 003, 028, 015, ! dy2(1:10, 2,13)
     . 041,-001,-005,-013,-003,-028, 003,-018, 006, 011, ! dy2(1:10, 3,13)
     . 021,-001, 002,-002, 003, 007, 009,-002, 010, 027, ! dy2(1:10, 4,13)
     . 027, 022, 016, 005, 010, 007, 019, 003, 018, 027, ! dy2(1:10, 5,13)
     . 030, 013,-002, 014, 007, 006, 009,-001, 002, 017, ! dy2(1:10, 6,13)
     . 006, 002,-002, 010, 004, 014, 004, 005, 003, 013, ! dy2(1:10, 7,13)
     . 005,-008,-014, 004,-008, 000,-008,-012,-021,-007, ! dy2(1:10, 8,13)
     . 000,-005,-017,-020,-011,-001,-012,-026,-022,-010, ! dy2(1:10, 9,13)
     .-016,-010,-020,-015,-013,-005,-013,-028,-022,-005, ! dy2(1:10,10,13)
     . 007, 009, 000,-002, 002, 006, 004,-017,-008, 005, ! dy2(1:10,11,13)
     .-017,-004,-008, 003,-005, 003,-001,-020,-012,-005, ! dy2(1:10,12,13)
     .-001,-014,-007, 019,-005,-007, 001,-012,-006,-003, ! dy2(1:10,13,13)
     .-012,-007,-002, 013,-009,-001, 008,-012, 000,-006, ! dy2(1:10,14,13)
     .-018,-007, 000, 000,-010,-012,-003,-012,-010,-013, ! dy2(1:10,15,13)
     .-012,-001, 006,-012,-008,-003, 002,-008,-005,-004, ! dy2(1:10,16,13)
     .-003,-005, 006,-015,-010,-009,-002,-006,-003,-018, ! dy2(1:10,17,13)
     .-004,-003, 011, 000,-009,-012,-002,-005,-003,-010, ! dy2(1:10,18,13)
     . 001, 004, 012, 000,-006,-008, 003, 005,-003,-008, ! dy2(1:10,19,13)
     . 003,-002, 009, 004,-012,-018,-006,-004,-008,-015, ! dy2(1:10,20,13)
     .-009,-005, 005, 004,-008,-016,-009, 004,-002,-015, ! dy2(1:10,21,13)
     .-002, 000, 006, 009,-004,-010,-008, 010,-007,-011, ! dy2(1:10,22,13)
     . 004, 010, 005,-002, 001,-009,-009, 009,-008,-006, ! dy2(1:10,23,13)
     .-011, 000, 002,-002, 005,-007,-012, 014,-003,-008, ! dy2(1:10,24,13)
     . 014, 011, 000, 002, 013,-003,-004, 015, 000, 000, ! dy2(1:10,25,13)
     . 011, 015, 002, 000, 016, 001,-001, 018, 005, 002, ! dy2(1:10,26,13)
     . 003, 009,-005, 000, 019, 007,-005, 021, 009, 006, ! dy2(1:10,27,13)
     . 017, 013, 005, 003, 028, 016, 003, 031, 012, 012, ! dy2(1:10,28,13)
     . 016, 018,-002, 010, 037, 024, 009, 035, 018, 020, ! dy2(1:10,29,13)
     . 013, 004,-015,-005, 020, 021, 012, 027, 018, 019, ! dy2(1:10,30,13)
     . 005, 010,-001,-003, 026, 030, 015, 031, 027, 031, ! dy2(1:10,31,13)
     . 015, 010,-015,-004, 027, 025, 014, 025, 022, 031, ! dy2(1:10,32,13)
     .-016,-009,-017,-014, 004, 023, 006, 012, 011, 014, ! dy2(1:10,33,13)
     .-002, 002,-006,-004, 008, 018, 010, 008, 011, 007, ! dy2(1:10,34,13)
     .-017, 005,-003,-008, 004, 007, 011, 001, 006, 008, ! dy2(1:10,35,13)
     . 008,-028,-020,-022,-019,-014,-009,-028,-011,-004, ! dy2(1:10,36,13)
     .-052,-019,-006,-026,-031,-025,-012,-055,-037,-035, ! dy2(1:10,37,13)
     .-067,-067,-023,-044,-042,-019,-034,-055,-047,-033, ! dy2(1:10,38,13)
     . 015, 011, 031,-021, 001, 000, 012, 015, 037,-010, ! dy2(1:10, 1,14)
     . 027,-029, 013,-012,-019,-034, 016, 011, 028, 000, ! dy2(1:10, 2,14)
     . 014, 003, 005, 002,-007,-004, 008, 011, 017, 023, ! dy2(1:10, 3,14)
     . 054, 003,-001,-010, 001,-007, 011, 002, 014, 009, ! dy2(1:10, 4,14)
     . 011, 006, 000, 019, 006, 004, 016, 013, 018, 019, ! dy2(1:10, 5,14)
     . 019, 015, 011, 016, 008, 006, 012, 008, 010, 017, ! dy2(1:10, 6,14)
     . 043, 007,-002, 000, 004, 011, 009, 004, 000, 022, ! dy2(1:10, 7,14)
     . 012,-006,-010,-010,-010, 000,-002,-006,-006, 003, ! dy2(1:10, 8,14)
     .-008,-030,-034,-027,-025,-021,-022,-037,-038,-020, ! dy2(1:10, 9,14)
     . 001,-021,-022,-024,-012,-006,-014,-033,-027,-009, ! dy2(1:10,10,14)
     . 008, 005,-008, 002,-001, 008, 003,-020,-010, 003, ! dy2(1:10,11,14)
     . 002, 000,-010, 010,-005,-002,-002,-021,-012,-004, ! dy2(1:10,12,14)
     .-010, 007,-005, 000,-005, 001, 005,-018,-006,-001, ! dy2(1:10,13,14)
     .-006,-005,-012,-004,-006,-004, 000,-019,-008,-009, ! dy2(1:10,14,14)
     .-006, 001, 000,-020,-006, 000, 003,-011,-010,-008, ! dy2(1:10,15,14)
     .-002, 001, 002,-019,-009,-006, 000,-010,-010,-020, ! dy2(1:10,16,14)
     . 003, 007, 007,-002,-002,-002, 004,-001, 001,-004, ! dy2(1:10,17,14)
     .-005,-003, 004, 000,-010,-012,-004,-006,-007,-013, ! dy2(1:10,18,14)
     . 003, 015, 011, 015,-002,-017,-003, 003,-002,-017, ! dy2(1:10,19,14)
     . 000, 000, 000, 005,-007,-016,-015,-003,-012,-018, ! dy2(1:10,20,14)
     .-015, 000, 000,-002,-009,-028,-013, 000,-018,-022, ! dy2(1:10,21,14)
     . 012, 002, 000,-002,-002,-018,-014, 007,-012,-015, ! dy2(1:10,22,14)
     .-008, 007, 000, 004, 008,-015,-009, 019,-004,-003, ! dy2(1:10,23,14)
     . 002, 010, 005,-002, 006,-010,-010, 018,-001,-012, ! dy2(1:10,24,14)
     . 005, 012, 000, 014, 021, 006,-003, 023, 014, 007, ! dy2(1:10,25,14)
     . 009, 028, 017, 009, 028, 011, 009, 033, 019, 017, ! dy2(1:10,26,14)
     . 005, 003, 003, 007, 021, 009,-002, 023, 017, 012, ! dy2(1:10,27,14)
     . 003, 004,-013,-003, 021, 011, 000, 026, 018, 015, ! dy2(1:10,28,14)
     . 007, 007,-005, 000, 019, 017, 009, 032, 032, 027, ! dy2(1:10,29,14)
     . 000,-002,-018,-005, 016, 015, 001, 021, 024, 026, ! dy2(1:10,30,14)
     . 017, 013,-003, 002, 023, 035, 025, 032, 040, 033, ! dy2(1:10,31,14)
     . 016, 016,-005,-002, 022, 036, 017, 028, 032, 026, ! dy2(1:10,32,14)
     .-010,-008,-010,-010, 002, 018, 007, 007, 019, 014, ! dy2(1:10,33,14)
     .-018,-002,-003,-005, 003, 020, 007, 006, 011, 023, ! dy2(1:10,34,14)
     . 010,-007, 002,-012,-004, 008, 012,-009, 005, 003, ! dy2(1:10,35,14)
     .-049,-019,-012,-012,-017,-012,-003,-032,-016,-007, ! dy2(1:10,36,14)
     .-039,-043,-019,-035,-039,-033,-021,-070,-039,-048, ! dy2(1:10,37,14)
     .-048,-037,-026,-031,-041,-035,-019,-076,-040,-042, ! dy2(1:10,38,14)
     .-021,-035, 013, 018,-030,-021, 018, 035, 024, 015, ! dy2(1:10, 1,15)
     . 022,-030, 003, 003,-022,-032, 003, 031, 027, 010, ! dy2(1:10, 2,15)
     . 000,-023, 007,-003,-010,-026, 003, 028, 019, 003, ! dy2(1:10, 3,15)
     . 039, 022, 011, 029, 011, 001, 021, 039, 028, 035, ! dy2(1:10, 4,15)
     . 024, 004, 010, 013, 000, 002, 019, 024, 014, 019, ! dy2(1:10, 5,15)
     . 017,-012,-013, 000,-008,-011, 000, 002,-005,-004, ! dy2(1:10, 6,15)
     . 026, 007, 003, 002, 002, 009, 010, 013, 007, 016, ! dy2(1:10, 7,15)
     . 003,-016,-017,-014,-012,-010,-007,-008,-012,-003, ! dy2(1:10, 8,15)
     .-026,-022,-034,-037,-018,-014,-018,-029,-032,-017, ! dy2(1:10, 9,15)
     .-020,-023,-030,-018,-022,-019,-019,-031,-023,-019, ! dy2(1:10,10,15)
     . 007,-006,-025, 003,-009, 000,-009,-023,-019, 000, ! dy2(1:10,11,15)
     .-011,-005,-020,-003,-009,-006,-007,-027,-014,-013, ! dy2(1:10,12,15)
     . 008, 004,-004,-011, 000, 004, 003,-011,-008, 005, ! dy2(1:10,13,15)
     . 000, 010,-005,-013,-003, 000, 009,-014,-002,-008, ! dy2(1:10,14,15)
     .-004,-004,-011,-033,-008,-006,-006,-023,-013,-009, ! dy2(1:10,15,15)
     .-004, 008, 004,-001,-004,-001, 000,-012,-009,-006, ! dy2(1:10,16,15)
     . 010, 011, 003, 000, 000,-001, 003,-001,-012,-007, ! dy2(1:10,17,15)
     .-009, 003, 004, 008,-001,-013,-005,-001,-008,-016, ! dy2(1:10,18,15)
     . 007, 005, 002, 006, 004,-008,-002, 006,-007, 004, ! dy2(1:10,19,15)
     .-008,-001,-003, 000,-003,-016,-007,-001,-015,-016, ! dy2(1:10,20,15)
     . 008, 014,-002,-001, 004,-010,-006, 015,-007,-002, ! dy2(1:10,21,15)
     . 023, 018, 006, 000, 015,-005,-005, 023, 000, 006, ! dy2(1:10,22,15)
     . 003, 009, 007,-001, 011,-003,-009, 023, 004,-004, ! dy2(1:10,23,15)
     .-010, 010,-003, 001, 014,-006,-005, 022, 005, 002, ! dy2(1:10,24,15)
     . 011, 009, 010,-002, 022, 008,-001, 027, 022, 007, ! dy2(1:10,25,15)
     . 023, 016, 011, 001, 024, 008, 003, 029, 025, 013, ! dy2(1:10,26,15)
     .-013,-006,-014,-007, 006,-001,-003, 021, 016,-002, ! dy2(1:10,27,15)
     .-006,-006,-007,-007, 015, 011,-001, 023, 026, 009, ! dy2(1:10,28,15)
     . 011, 011,-004, 003, 017, 011, 012, 032, 041, 021, ! dy2(1:10,29,15)
     .-009,-007,-018,-007, 010, 018, 009, 027, 036, 013, ! dy2(1:10,30,15)
     . 011, 017, 001, 007, 018, 025, 020, 028, 039, 024, ! dy2(1:10,31,15)
     . 000, 008,-012,-003, 010, 021, 016, 022, 040, 015, ! dy2(1:10,32,15)
     .-009,-012, 000,-011,-004, 011, 009, 001, 007, 008, ! dy2(1:10,33,15)
     . 006, 012, 008, 006, 007, 022, 019, 001, 021, 003, ! dy2(1:10,34,15)
     .-009, 007, 013, 000, 000, 007, 016,-012, 002, 014, ! dy2(1:10,35,15)
     .-020,-010,-002,-012,-015,-012, 002,-041,-015,-018, ! dy2(1:10,36,15)
     .-057,-067,-012,-023,-035,-037,-017,-082,-040,-033, ! dy2(1:10,37,15)
     .-048,-042, 003,-020,-030,-024,-019,-075,-047,-050, ! dy2(1:10,38,15)
     . 012,-014, 046, 025,-014,-023, 028, 071, 065, 003, ! dy2(1:10, 1,16)
     .-020,-024, 040, 021,-020,-018, 013, 065, 047, 012, ! dy2(1:10, 2,16)
     . 002,-022, 007, 000,-021,-020,-002, 033, 023, 004, ! dy2(1:10, 3,16)
     . 013,-012, 004, 021,-010,-012, 008, 048, 020, 009, ! dy2(1:10, 4,16)
     . 006, 008, 018, 026, 017,-004, 022, 042, 034, 010, ! dy2(1:10, 5,16)
     . 016, 002, 001,-003,-002, 001, 007, 021, 010, 018, ! dy2(1:10, 6,16)
     . 012, 006, 000,-004, 006, 007, 009, 003, 006, 009, ! dy2(1:10, 7,16)
     . 000,-010,-011,-013,-004, 000, 003,-001,-006, 015, ! dy2(1:10, 8,16)
     . 001,-028,-032,-016,-025,-025,-018,-034,-026,-023, ! dy2(1:10, 9,16)
     .-004,-020,-019,-001,-012,-007,-010,-027,-023,-003, ! dy2(1:10,10,16)
     .-005, 008,-007, 007,-005,-003,-001,-021,-010,-013, ! dy2(1:10,11,16)
     . 016, 004,-020,-010,-005,-004,-006,-025,-016, 003, ! dy2(1:10,12,16)
     . 008, 006,-002,-026, 000,-003, 001,-029,-012,-017, ! dy2(1:10,13,16)
     . 004, 010,-004,-024,-002, 002,-001,-026,-019,-006, ! dy2(1:10,14,16)
     .-002,-010,-021,-020,-010,-013,-008,-037,-020,-027, ! dy2(1:10,15,16)
     .-009,-001,-017,-014,-012,-011,-009,-027,-021,-015, ! dy2(1:10,16,16)
     .-001,-004,-006, 000,-008,-019,-010,-019,-019,-021, ! dy2(1:10,17,16)
     .-016, 000,-009, 002,-006,-018,-002,-009,-021,-004, ! dy2(1:10,18,16)
     .-003, 005, 000,-012,-004,-024,-008,-001,-010,-020, ! dy2(1:10,19,16)
     . 001, 011,-005,-014,-003,-015,-010,-001,-013,-012, ! dy2(1:10,20,16)
     . 007, 000,-007,-014, 000,-018,-012, 008,-015,-016, ! dy2(1:10,21,16)
     .-007,-003,-008,-002, 005,-015,-009, 015,-012,-015, ! dy2(1:10,22,16)
     . 003, 006,-003,-004, 009,-011,-011, 016,-005,-006, ! dy2(1:10,23,16)
     . 001,-005,-007,-016, 007,-007,-015, 015,-002,-009, ! dy2(1:10,24,16)
     .-020,-005,-008,-011, 005,-006,-006, 017, 012,-007, ! dy2(1:10,25,16)
     .-010,-003,-007,-009, 009, 000,-008, 021, 011,-004, ! dy2(1:10,26,16)
     .-009,-008,-001,-012, 010, 004, 000, 024, 026, 002, ! dy2(1:10,27,16)
     .-009,-002,-014,-011, 011, 007, 005, 025, 033, 007, ! dy2(1:10,28,16)
     . 004, 000, 001,-005, 012, 022, 014, 028, 043, 024, ! dy2(1:10,29,16)
     .-008,-009,-007,-014, 003, 004, 005, 016, 031, 010, ! dy2(1:10,30,16)
     .-006, 022, 008, 012, 018, 035, 029, 028, 044, 027, ! dy2(1:10,31,16)
     . 000, 008, 013, 006, 015, 029, 023, 027, 043, 025, ! dy2(1:10,32,16)
     .-023, 006, 002, 000, 000, 018, 016,-001, 015, 008, ! dy2(1:10,33,16)
     .-004, 000, 011, 004,-004, 014, 019,-013, 006, 004, ! dy2(1:10,34,16)
     .-020,-017, 010,-001,-004, 002, 015,-034,-009,-015, ! dy2(1:10,35,16)
     .-025,-020, 005,-007,-014,-017, 005,-063,-028,-021, ! dy2(1:10,36,16)
     .-031,-024, 013,-004,-023,-022,-002,-079,-040,-036, ! dy2(1:10,37,16)
     .-048,-047, 001,-016,-035,-025,-015,-098,-040,-035, ! dy2(1:10,38,16)
     . 014,-023, 047, 034,-011,-030, 030, 093, 061, 000, ! dy2(1:10, 1,17)
     . 022,-006, 053, 034,-014,-012, 020, 097, 067,-016, ! dy2(1:10, 2,17)
     . 006,-019, 027, 033,-015,-024, 004, 063, 034,-002, ! dy2(1:10, 3,17)
     .-002,-017, 023, 029,-010,-011, 012, 058, 030,-011, ! dy2(1:10, 4,17)
     . 021, 001, 016, 022,-003, 007, 011, 056, 022, 023, ! dy2(1:10, 5,17)
     . 016,-013, 013, 010, 005, 002, 016, 033, 020, 013, ! dy2(1:10, 6,17)
     . 010, 044, 011, 021, 015, 022, 028, 028, 021, 023, ! dy2(1:10, 7,17)
     . 015,-018,-003, 002,-006,-003, 005, 002, 000,-005, ! dy2(1:10, 8,17)
     .-020,-010,-021,-001,-011,-012,-006,-026,-017,-006, ! dy2(1:10, 9,17)
     . 005,-007,-017,-012,-011,-017,-010,-034,-017,-021, ! dy2(1:10,10,17)
     . 025, 006,-006,-010, 000, 005, 004,-021,-011, 008, ! dy2(1:10,11,17)
     . 012, 002, 000,-019, 002, 000, 002,-028,-007,-010, ! dy2(1:10,12,17)
     . 002, 014,-006,-013, 001, 000, 007,-032,-012, 007, ! dy2(1:10,13,17)
     . 016,-001,-011,-013,-006,-006, 002,-036,-016,-002, ! dy2(1:10,14,17)
     .-014, 002,-008,-009,-004,-009,-006,-041,-024,-012, ! dy2(1:10,15,17)
     .-007,-007,-015,-002,-008,-020,-006,-038,-029,-014, ! dy2(1:10,16,17)
     .-030,-011,-019,-013,-009,-020,-009,-034,-026,-025, ! dy2(1:10,17,17)
     . 005,-003,-007,-009,-001,-014,-006,-021,-019,-014, ! dy2(1:10,18,17)
     . 004, 005,-008,-018, 005,-012,-004,-001,-015,-014, ! dy2(1:10,19,17)
     .-018, 000,-011,-020,-001,-022,-010,-011,-032,-023, ! dy2(1:10,20,17)
     . 004, 000,-014,-012,-001,-021,-012, 001,-021,-009, ! dy2(1:10,21,17)
     .-002, 005,-013,-007, 008,-013,-007, 004,-011,-011, ! dy2(1:10,22,17)
     . 001, 002,-007,-010, 006,-010,-007, 015,-007,-007, ! dy2(1:10,23,17)
     .-009,-002,-004,-010, 007,-001,-005, 017, 004,-005, ! dy2(1:10,24,17)
     . 004, 001,-003,-011, 011, 009,-001, 022, 012, 004, ! dy2(1:10,25,17)
     .-005,-004, 001,-009, 009, 011, 000, 021, 016,-011, ! dy2(1:10,26,17)
     .-009, 000,-011,-014, 006, 011, 000, 020, 019, 003, ! dy2(1:10,27,17)
     .-010,-002, 000,-006, 005, 009, 007, 020, 028,-003, ! dy2(1:10,28,17)
     . 000, 009,-001, 004, 013, 017, 015, 028, 036, 017, ! dy2(1:10,29,17)
     .-005, 009, 003, 000, 000, 021, 016, 027, 035, 013, ! dy2(1:10,30,17)
     . 016, 012, 022, 011, 017, 028, 026, 021, 034, 019, ! dy2(1:10,31,17)
     . 019, 003, 008, 012, 008, 018, 019, 013, 021, 016, ! dy2(1:10,32,17)
     .-014, 000, 007,-006,-016, 008, 006,-032,-002,-012, ! dy2(1:10,33,17)
     .-001, 001, 015,-002,-003, 011, 015,-036,-009, 001, ! dy2(1:10,34,17)
     .-020, 003, 022, 003,-003, 006, 017,-051,-005,-007, ! dy2(1:10,35,17)
     .-036,-013, 019,-005,-012,-009, 005,-079,-031,-020, ! dy2(1:10,36,17)
     .-057,-048, 012,-015,-033,-034,-004,-114,-058,-044, ! dy2(1:10,37,17)
     .-048,-028, 018,-010,-024,-031,-005,-112,-071,-047, ! dy2(1:10,38,17)
     . 034, 010, 126, 063,-001, 007, 028, 147, 094, 015, ! dy2(1:10, 1,18)
     . 000, 012, 131, 111, 012,-017, 043, 147, 094, 005, ! dy2(1:10, 2,18)
     .-003,-012, 067, 080,-002,-028, 019, 108, 062, 003, ! dy2(1:10, 3,18)
     . 025, 008, 072, 054, 008,-002, 023, 091, 053, 005, ! dy2(1:10, 4,18)
     .-007,-008, 042, 031, 010,-007, 020, 067, 034, 010, ! dy2(1:10, 5,18)
     . 004, 007, 030, 005, 000,-013, 008, 041, 014,-009, ! dy2(1:10, 6,18)
     . 043, 000, 026, 043, 022, 017, 028, 023, 021, 012, ! dy2(1:10, 7,18)
     . 007, 008, 019, 034, 008,-003, 016, 019, 012, 009, ! dy2(1:10, 8,18)
     . 026, 003, 009, 004,-002,-004, 001,-019,-003,-017, ! dy2(1:10, 9,18)
     . 003,-001, 009, 002, 001, 001, 010,-025,-004, 003, ! dy2(1:10,10,18)
     . 000, 022, 016,-018, 004,-001, 013,-040,-008, 000, ! dy2(1:10,11,18)
     . 022, 013, 004,-002, 002, 005, 011,-039,-011, 014, ! dy2(1:10,12,18)
     . 043, 038, 023, 022, 016, 016, 018,-040, 000, 012, ! dy2(1:10,13,18)
     . 004, 039, 011, 013, 012, 010, 016,-043,-013, 006, ! dy2(1:10,14,18)
     . 028, 014, 010, 005, 012, 000, 017,-045,-005, 000, ! dy2(1:10,15,18)
     . 020, 028, 010, 003, 019, 005, 014,-033,-011, 010, ! dy2(1:10,16,18)
     . 027, 011, 003,-004, 012, 000, 009,-035,-013, 005, ! dy2(1:10,17,18)
     .-002, 019,-002,-015, 015,-003, 007,-035,-016, 000, ! dy2(1:10,18,18)
     . 004, 003, 011, 011, 008,-008, 007,-013,-017, 006, ! dy2(1:10,19,18)
     . 009,-006,-008,-009,-002,-030,-012,-026,-028,-007, ! dy2(1:10,20,18)
     . 018, 011, 000, 000, 012, 000, 004,-001,-020,-013, ! dy2(1:10,21,18)
     . 004, 007, 006,-013, 005,-003,-008, 001,-021,-001, ! dy2(1:10,22,18)
     . 004, 003,-006,-005, 005, 004,-005, 005,-010, 005, ! dy2(1:10,23,18)
     . 012, 026, 004,-007, 015, 017, 007, 013, 000, 004, ! dy2(1:10,24,18)
     .-006, 014, 000, 000, 004, 019, 006, 015, 013, 002, ! dy2(1:10,25,18)
     . 006, 024, 003,-005, 012, 021, 010, 012, 009, 015, ! dy2(1:10,26,18)
     .-006, 010, 008, 001,-002, 024, 006, 013, 016, 000, ! dy2(1:10,27,18)
     .-010, 005, 005, 000, 007, 024, 013, 017, 017, 011, ! dy2(1:10,28,18)
     . 003, 013, 015,-002, 000, 032, 024, 010, 017,-003, ! dy2(1:10,29,18)
     . 000,-008, 004, 000,-004, 017, 008,-008,-002,-005, ! dy2(1:10,30,18)
     .-003, 038, 036, 028, 000, 024, 032,-008,-001, 017, ! dy2(1:10,31,18)
     . 003, 027, 027, 006,-002, 019, 025,-023,-001,-008, ! dy2(1:10,32,18)
     .-032,-008, 002, 002,-008, 002, 002,-057,-035,-022, ! dy2(1:10,33,18)
     .-032,-006, 023,-001,-010, 000, 013,-077,-024,-023, ! dy2(1:10,34,18)
     .-032,-023, 029, 007,-014, 008, 012,-095,-042,-027, ! dy2(1:10,35,18)
     .-034,-022, 029, 001,-021,-010,-004,-125,-055,-039, ! dy2(1:10,36,18)
     .-057,-041, 041,-010,-034,-004,-013,-156,-083,-062, ! dy2(1:10,37,18)
     .-110,-040, 041,-031,-052, 001,-004,-167,-067,-105/ ! dy2(1:10,38,18)

c----------------------------------------------------------------------


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     POLYNOMIAL COEFFICIENTS
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c
c J to preserve Jena's structure as much as possible, I copy
c   these over ; but I could have just used them directly below....
c
      do i = 1, 9
         a(i) = apoly(i,chip,1+filter)
         b(i) = bpoly(i,chip,1+filter)
         enddo


CCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     CORRECTIONS
C     
CCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCC


c
c J compute where the input point is in terms of the x array 
c   indices ; the goal is to just find the one below and the 
c   one above; the challenge is that there are zones of 6, 
c   and there's a repeat at the end of one and the beginning
c   of the next ;  "nx" is the lower index and "fx" tells you
c   where you are between "nx" and "nx+1"
c
      zx    = int((xr-23.5)/675+100)-100
      xr_zx = 23.5 + zx*675
      rzx   = (xr-xr_zx)/112.5 - 0.5
      if (rzx.lt.0.000) rzx = 0.000
      if (rzx.gt.4.999) rzx = 4.999

      rx    = 2 + (zx)*6 + rzx


      if (rx.lt. 2) rx =  2.0001
      if (rx.gt.37) rx = 36.9999
      if (zx.eq.-1) rx =  1.0001
      if (zx.eq. 6) rx = 37.9999


      nx = int(rx)
      fx = (xr         -x_tab(nx)) /
     .     (x_tab(nx+1)-x_tab(nx))

      if (zx.eq.-1) then
          nx =  1
          fx =  0.0001
          endif
      if (zx.eq. 6) then
          nx = 37
          fx =  0.9999
          endif


c
c the y-index stuff is a bit easier, but only
c a bit.  Again, "ny" is the lower index and "fy"
c tells you where you are between "ny" and "ny+1"
c

      ry = 2 + (yr-171)/114
      if (ry.lt.01.0001) ry = 01.0001
      if (ry.gt.17.9999) ry = 17.9999
  
      if (chip.eq.1) then
         if (yr.le. 911.and.ry.gt. 7.999) ry =  7.999 
         if (yr.gt. 911.and.ry.lt. 9.001) ry =  9.001
         ny = int(ry)
         fy = (yr          -y1_tab(ny)) /
     .        (y1_tab(ny+1)-y1_tab(ny))
         endif

      if (chip.eq.2) then
         if (yr.le.1140.and.ry.gt. 9.999) ry =  9.999 
         if (yr.gt.1140.and.ry.lt.11.001) ry = 11.001
         ny = int(ry)
         fy = (yr          -y2_tab(ny)) /
     .        (y2_tab(ny+1)-y2_tab(ny))
         endif



      if (chip.eq.1) then
         dxt = (1-fx)*(1-fy)*dx1(01+filter,nx,  ny  )*0.001
     .       + (1-fx)*( fy )*dx1(01+filter,nx,  ny+1)*0.001
     .       + ( fx )*(1-fy)*dx1(01+filter,nx+1,ny  )*0.001
     .       + ( fx )*( fy )*dx1(01+filter,nx+1,ny+1)*0.001
         dyt = (1-fx)*(1-fy)*dy1(01+filter,nx  ,ny  )*0.001
     .       + (1-fx)*( fy )*dy1(01+filter,nx  ,ny+1)*0.001
     .       + ( fx )*(1-fy)*dy1(01+filter,nx+1,ny  )*0.001
     .       + ( fx )*( fy )*dy1(01+filter,nx+1,ny+1)*0.001
  199    format(a6,3x,f10.4,1x,i2,1x,i2,
     .             5x,4f10.4,5x,4f8.4)
         endif

      if (chip.eq.2) then
         dxt = (1-fx)*(1-fy)*dx2(01+filter,nx  ,ny  )*0.001
     .       + (1-fx)*( fy )*dx2(01+filter,nx  ,ny+1)*0.001
     .       + ( fx )*(1-fy)*dx2(01+filter,nx+1,ny  )*0.001
     .       + ( fx )*( fy )*dx2(01+filter,nx+1,ny+1)*0.001
         dyt = (1-fx)*(1-fy)*dy2(01+filter,nx  ,ny  )*0.001
     .       + (1-fx)*( fy )*dy2(01+filter,nx  ,ny+1)*0.001
     .       + ( fx )*(1-fy)*dy2(01+filter,nx+1,ny  )*0.001
     .       + ( fx )*( fy )*dy2(01+filter,nx+1,ny+1)*0.001
          endif


C     MEMO
C     filter:  
C     0 -> F225W
C     1 -> F275W
C     2 -> F336W
C     3 -> F390W
C     4 -> F438W
C     5 -> F555W
C     6 -> F606W
C     7 -> F775W
C     8 -> F814W
C     9 -> F850L


CCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     POLYNOMIAL
C    
CCCCCCCCCCCCCCCCCCCCCCCCCC


c     ----> input xt, yt ----> output xc, yc (polynomial correction)
c
c J I did find it strange that you used (xt,yt) here and not
c   (xr,yr)...
c

                                ! normalized positions

      x_cen=2048
      y_cen=1025
      xt = xr + dxt
      yt = yr + dyt
      rx=(xt-x_cen)/x_cen
      ry=(yt-y_cen)/y_cen

                                !corrections

      dxp = (a(3)*rx*rx)
     .     +(a(4)*rx*ry)
     .     +(a(5)*ry*ry)
     .     +(a(6)*rx*rx*rx)
     .     +(a(7)*rx*rx*ry)
     .     +(a(8)*rx*ry*ry)
     .     +(a(9)*ry*ry*ry)

      dyp = (b(1)*rx)
     .     +(b(2)*ry)
     .     +(b(3)*rx*rx)
     .     +(b(4)*rx*ry)
     .     +(b(5)*ry*ry)
     .     +(b(6)*rx*rx*rx)
     .     +(b(7)*rx*rx*ry)
     .     +(b(8)*rx*ry*ry)
     .     +(b(9)*ry*ry*ry)

      xc = xr + dxp + dxt
      yc = yr + dyp + dyt

      return
      end ! subroutine WFC3_UVIS_gcj










c--------------------------------------------------------
c
c This is a wrapper that takes Jay's more compact 
c    formulation of Jena's correction and puts it into
c    a meta coordinate system.  
c
c The properties of this system are that 
c    1) the bottom chip is in the same frame as Jena's, 
c       but has y --> y+200 to prevent negative coordinates 
c    2) the top chip is transformed using linear transformations 
c       into the frame implied by the bottom chip
c
      subroutine wfc3uv_gcj_meta_k(xr,yr,chip,filter,xc,yc)
      implicit none

      real*8 xr, yr
      integer chip
      integer filter
      real*8 xc, yc

      real*8 xru, yru
      real*8 xcu, ycu
      real*8 fx, fy
      real*8 skew1
      real*8 skew2
      real*8 a(9), b(9)

      real    rx,  ry
      integer nx,  ny
      real    ddx, ddy

      real filt_shift_x(10)
      data filt_shift_x / 0.200,  0.265,  0.440,  0.280, -0.075,
     .                   -0.085,  0.000,  0.260, -0.070,  0.155/

      real filt_shift_y(10)
      data filt_shift_y / 0.265,  0.120, -0.315, -0.330, -0.520,
     .                   -0.685,  0.000, -1.480, -0.455, -0.920/

      real filt_scl(10)
      data filt_scl / 1.000110, 0.999997, 0.999851, 0.999894, 1.000007,
     .                1.000121, 0.999997, 1.000004, 1.000048, 1.000158/

      real skewxx_(10), skewxy_(10), offxxx_(10)
      real skewyx_(10), skewyy_(10), offyyy_(10)
      real*8 xch, ych


      data offxxx_ /   -2.024,   -2.001,   -2.023,   -2.001,   -2.003,
     .                 -1.996,   -2.020,   -1.975,   -2.048,   -1.969/

      data offyyy_ / 2073.362, 2073.354, 2073.383, 2073.411, 2073.402,
     .               2073.434, 2073.410, 2073.409, 2073.463, 2073.389/

      data skewxx_ /  -12.152,  -12.123,  -12.068,  -12.098,  -12.097,
     .                -12.063,  -12.035,  -11.995,  -12.038,  -12.027/

      data skewxy_ /   -2.355,   -2.318,   -2.384,   -2.391,   -2.368,
     .                 -2.356,   -2.367,   -2.346,   -2.384,   -2.299/

      data skewyx_ /    2.306,    2.310,    2.367,    2.397,    2.367,
     .                  2.360,    2.363,    2.339,    2.387,    2.300/

      data skewyy_ /  -12.146,  -12.116,  -12.068,  -12.110,  -12.096,
     .                -12.046,  -12.031,  -11.987,  -12.033,  -12.017/


      if (chip.ne.1.and.chip.ne.2) stop 'meta_k: invalid chip'

      if (chip.eq.1) then
           xru = xr
           yru = yr
           call wfc3uv_gcj(xru,yru,1,filter,xcu,ycu)
           xc  = xcu
           yc  = ycu
           endif

      if (chip.eq.2) then
           xru = xr
           yru = yr
           call wfc3uv_gcj(xru,yru,2,filter,xcu,ycu)

           fx  = (xcu-2048.00)/2048.0
           fy  = (ycu-1024.00)/2048.0

           xc  = xcu + offxxx_(1+filter) 
     .               + skewxx_(1+filter)*fx
     .               + skewxy_(1+filter)*fy

           yc  = ycu + offyyy_(1+filter) 
     .               + skewyx_(1+filter)*fx
     .               + skewyy_(1+filter)*fy


         endif

      xc = xc + filt_shift_x(1+filter)
      yc = yc + filt_shift_y(1+filter)



c-------------------------------------------------------
c
c include the field expansion
c
      xc = 2048 + (xc-2048)/filt_scl(1+filter)
      yc = 2048 + (yc-2048)/filt_scl(1+filter)


c-----------------------------------------------------
c
c J this is just to keep things from going negative
c   coordinate-wise.
c
      yc = yc + 200 

      return
      end




c--------------------------------------------------------
c
c This particular wrapper takes a position in the 
c    4096x4096 abutted frame and finds the meta-frame 
c    position.
c
      subroutine wfc3uv_gcj_meta(xr,yr,filter,xc,yc)
      implicit none
      real*8 xr, yr
      integer filter
      real*8 xc, yc

      real*8 yru
      integer chip


      if (yr.lt.2048) then
         yru = yr
         call wfc3uv_gcj_meta_k(xr ,yr ,1,filter,xc,yc)
         return
         endif 

      if (yr.gt.2048) then
         yru = yr-2048
         call wfc3uv_gcj_meta_k(xr ,yru,2,filter,xc,yc)
         return
         endif 

      return          
      end

                                              
c*********************************************
c****                                         
c**** #include "/home/jayander/FORTRAN/MAC_OSX/FITS/readfits_i2r.f"
c****                                         
c*********************************************
                                              
 
 
      subroutine readfits_i2r(FILE,pix,NX,NY)
      implicit none
 
      character*80 FILE
      integer NX, NY
      real*4 pix(NX,NY)
 
      character*70 INFO(10)
      common / fitsinfo / INFO
 
      integer naxes
      integer laxis(3)
 
      character*8  field
      character*20 stream
 
      real*4 pixu
      integer nbyte0
      integer nbyteE
      integer nbyte1
      integer nbyte2
      integer nbper
      integer i,ios, k
      integer ii, jj
      integer n
      integer NXU, NYU
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffc,buffb)
 
      integer*2 ibuff(1440)
 
      integer ifirst, i1, i2
      integer j
 
      integer np1, np2, npt
      integer nextend
      integer nread
      real bscale, bzero
      integer bitpix
 
      integer r2i
 
      logical LINUX
      data LINUX/.true./
      logical DIAG
      data DIAG /.false./
 
 
      character*70 HDR(25)
      common/HDR/HDR
      character*80 FILEU
 
      FILEU = FILE
      do i = 75,2,-1
         if (FILE(i:i+4).eq.'.fits') FILEU = FILE(1:i+4)
         enddo
 
      open(10,file=FILEU,status='old',
     .     err=900,recl=2880,form='UNFORMATTED',
     .     access='DIRECT')
 
      if (DIAG) print*,'...opened'
 
      bscale = 1
      bzero  = 0
 
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
      HDR(24) = '  1.000 '
 
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
         if (field.eq.'FILTNAM1') INFO(02) = stream
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
 
      if (BITPIX.ne.16) then
         print*,'readfits_i2...: '
         print*,'                '
         print*,'   you called a routine to read in an'
         print*,'   unsigned i2 image, and the image you'
         print*,'   gave it has BITPIX = ',BITPIX
         print*,'   '
         print*,'   FILEU: ',FILEU
         print*,'   '
         stop
         endif
 
      nbper = 2*laxis(1)*laxis(2)
      npt   =   laxis(1)*laxis(2)
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
      if (laxis(1).gt.NX.and.laxis(2).gt.NY) then
         print*,' not enough image space! '
         print*,' '
         print*,'  laxis1: ',laxis(1)
         print*,'      NX: ',NX
         print*,' '
         print*,'  laxis2: ',laxis(2)
         print*,'      NY: ',NY
         print*,' '
         stop
         endif
 
      NXU = laxis(1)
      NYU = laxis(2)
 
      if (DIAG) then
         print*,'     NX: ',NX
         print*,'     NY: ',NY
         print*,'  NBPER: ',nbper
         print*,'  NBYT1: ',nbyte1
         print*,'  NBYT2: ',nbyte2
         print*,' IFIRST: ',ifirst
         print*,'     I1: ',i1
         print*,'     I2: ',i2
         print*,'    NPT: ',NPT
         endif
      do i = i1, i2, 1
         read(10,rec=i,iostat=ios) buffc
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/2 + 1
         np2 = (nbyteE-nbyte1)/2 + 1
         np2 = min(np2,npt)
         call buff2pix_i2r(buffb,ibuff,0001,1440)
         do n = np1, np2, 1
            jj = n/NXU   + 1
            ii = n-NXU*(jj-1)
            pixu = ibuff(n-np1+1)*bscale+bzero
            pix(ii,jj) = pixu
            enddo
         enddo
 
      if (DIAG) then
         print*,'  NBPER: ',nbper
         print*,'  NBYT1: ',nbyte1
         print*,'  NBYT2: ',nbyte2
         print*,' IFIRST: ',ifirst
         print*,'     I1: ',i1
         print*,'     I2: ',i2
         print*,'    NPT: ',NPT
         endif
 
 
      return
 
 900  continue
      print*,'READFITS_I2 ERROR'
      print*,'   FILEU: ',FILEU
      stop
 
      end
 
 
 
      subroutine buff2pix_i2r(buff,pix,n1,nt)
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
 
 
                                              
c*********************************************
c****                                         
c**** #include "/home/jayander/FORTRAN/MAC_OSX/BUFF/buff_work.f"
c****                                         
c*********************************************
                                              
 
 
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
 
 
c-------------------------------------------------------
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
 
 
 
 
      subroutine buff2pix_i4(buff,pix,n1,nt)
      implicit none
      byte buff(2880)
      integer*4 pix(*)
      integer n1,nt
 
      byte b(4)
      integer ii
      equivalence(ii,b)
 
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
         if (npu.ge.1.and.npu.le.nt) pix(npu) = ii
         enddo
 
      return
      end
 
 
 
      subroutine buff2pix_i2(buff,pix,n1,nt)
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
 
 
c-------------------------------------------------------
c
c
      subroutine buff2pix_i2_edge(buff,pix,n1,nt,
     .                            NXP,NYP,NXF,NYF)
      implicit none
      byte buff(2880)
      integer NXP,NYP
      integer*2 pix(NXP,NYP)
      integer n1,nt
      integer NXF,NYF
 
      real pbuff(720)
      common /sneaky/pbuff
      integer i
      integer npu, nbu
      integer NX, NY
 
      byte b(4)
      integer ii
      equivalence(ii,b)
 
      do i = 1, 1440
         npu = n1+i-1
         nbu = (i-1)*2
         if (npu.ge.1.and.npu.le.nt) then
            NX = npu - (npu-1)/NXF*NXF
            NY = 1 + (npu-1)/NXF
            if (NX.le.NXP.and.NY.le.NYP) then
               if (.not.(_LINUX_)) then
                  b(1) = buff(nbu+1)
                  b(2) = buff(nbu+2)
                  endif
               if ((_LINUX_)) then
                  b(2) = buff(nbu+1)
                  b(1) = buff(nbu+2)
                  endif
               if (npu.ge.1.and.npu.le.nt) pix(NX,NY) = ii
               endif
            endif
         enddo
 
      return
      end
 
 
      subroutine pix2buff_i2(buff,pix,n1,nt)
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
 
 
      subroutine pix2buff_i4(buff,pix,n1,nt)
      implicit none
      byte buff(2880)
      integer*4 pix(*)
      integer n1,nt
 
      byte b(4)
      integer ii
      equivalence(ii,b)
 
      integer i, npu, nbu
 
      do i = 1, 720
         npu = n1+i-1
         nbu = (i-1)*4
         if (npu.ge.1.and.npu.le.nt) ii = pix(npu)
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
 
 
 
 
 
