 
c---------------------------------------------------------
c
c **********************************************
c **********************************************
c **********************************************
c *****                                    *****
c *****                                    *****
c *****           img2xym_WFC.09x10.F      *****
c *****                                    *****
c *****                                    *****
c **********************************************
c **********************************************
c **********************************************
c
c
c This is a FORTRAN program that will read in an image and find all
c the stars in that image that meet some simple criteria in terms
c of brightness and isolation.  It will then use a 9x10 array of PSFs
c to measure positions and fluxes for the stars.  Note the routine is designed
c to measure stars that are relatively isolated.  This means that if there
c are nearby neighbors within (say) 2 pixels, the finding and measurement
c of both the stars will be compromised.
c
c a typical run will look like:
c
c ./img2xym_WFC.09x10.e 5 50 999999 "PSFIMAGE.F606W.fits"
c
c             [PERT] [SUBT] [QSEL] [APSTUDY] [CONST] j92301veq_flt.fits
c
c WARNING: be sure to search on "LINUX" in this code and to set every flag
c appropriately.
c
c
c
c ***************
c *             *
c * PARAMETERS  *
c *             *
c ***************
c
c   (1) "5" refers to the isolation index; it will only find stars
c           in pixels that have no brighter neighbors within 5 pixels
c           (if this were "-5", then it would be even more strict; not
c           only could the peak not have anything brighter within
c           5 pixels, but it couldn't even have anything *reasonably*
c           bright within 5 pixels.  More stringent.
c
c   (2) "50" refers to the minimum flux over sky a star must have to
c           be included in the find;  for each local max, it adds up
c           the 4 brightest pixels about it and sees if this sum is
c           least "50" DN above sky.
c
c   (3) "999999" This  means that a star can have an essentially infinite
c           flux and still be included in the sample; it measures
c           saturated stars reasonably well.  If you don't want to
c           measure saturated stars, set this at 54999.
c
c   (4) "PSF..." This is the name of the PSF to be used.  You have a
c           couple options here.  The PSF is designed for use on _flt images,
c           but you can just have it find centroid positions and use crude
c           aperture photometry.  Run the program without arguments to see
c           how to call it with aperture photometry.
c
c   (5) "j..." list of images (can include asterisks to do multiple files,
c           very very handy!).  These images can be in the form of _flt's
c           or what I usually use "_WI4.fits"... it's an integer*4 version
c           of the _flts.  It is much more easily compressible than real*4's.
c
c   (+) If any of the arguments #5 onwards is "PERT" then the routine will
c       construct a PSF perturbation to add to the library PSF.
c
c   (+) If any of the arguments #5 onwards is "SUB" then the routine will
c       generate a file "subtract.fits" that removes the identified sources
c       from the image.  This file will be overwritten if there are multiple
c       image arguments given.
c
c
c **************
c *            *
c * DISCUSSION *
c *            *
c **************
c
c There are 3 basic parameters that a point-source can have: the 2-d position
c and the flux.  In order to measure these parameters, you also have to have
c an idea of the sky value.
c
c In this routine I measure:
c
c a) SKY: By taking an annulus as close as possible to the star's central
c         pixel.  I go out to where the star's flux is less than 10% of the
c         sky value.  I then use an annulus of 3 pixels in width with the
c         inner radius at the 10% location.  For faint stars, this means
c         measuing sky between 2 and 5 pixels.  For bright stars, it means
c         measuring sky between 8 and 11 pixels.  For saturated stars, sky
c         is measured farther out.  In all cases, though, the PSF is used
c         to estimate the contribution of the star to the sky annulus, so
c         that we can get an accurate idea of what the background is.
c
c b) POSITION: To find the best position, I use only the inner 5x5 pixels
c         centered on the star's brightest pixel.  To find the optimum
c         position, I do a simple grid-search for the position that minimizes
c         the sum of the residuals.  In this sum, each pixel is weighted by
c         1/P(i,j), which is 1/sigma**2 from the expected poisson noise in
c         that pixel.
c
c c) FLUX: Once we have a position and a sky value, the last thing to solve
c         for is the flux.  I again use proper noise-based weighting.  This
c         reduces to PSF-based "aperture" photometry for the inner pixels
c         of bright stars, but reduces to PSF-weighting for the lower S/N
c         pixels, such as the central pixels of faint stars.  In choosing
c         the optimal aperture, I start with the inner 3x3 pixels.  Then I
c         include successive annuli in the aperture, so long as the successive
c         annuli do not have pixels in them that are inconsistent with the
c         flux (as determined by the central pixels).  I should add that
c         the fluxes I get out are quite good, but it's hard to say for
c         sure that they are "optimal" for every situation.  Often non-random
c         concerns can be more important than the random issues.  Such
c         non-random issues are:  errors in the PSF, crowding by brighter
c         and fainter neighbors, slopes in the background, fluctuations
c         in the background, CRs, warm pixels, etc.  I do not do any
c         correction for warm pixels; I figure that stars that are affected
c         by such will be terribly measured, and that will come out when
c         you compare them against other (unaffected) observations.
c
c The PSF has been constructed by img2psf_WFC.09x10, a program that would
c be very very hard to export.  I have found that it has systematically
c good fluxes to within 1% or so.  The PSF is just stored in a simple
c fits file, and is supersampled by a factor of 4.
c
c This routine also allows for a "perturbation" of the PSF.  I've noticed
c that the PSF does change from exposure to exposure. The pertpsf will start
c with the spatially varying 09x10 PSF, which is a function of (dx,dy,ix,iy)
c then will add to it a PERTPSF, which doesn't vary spatially, but is just
c a function of (dx,dy).  This will allow the PSF to be slightly different
c from exposure to exposure.   To have the routine derive a PERTPSF from
c each image, then specify "+" on the command line before specifying the PSF
c file.
c
c OUTPUTS FOR THE PROGRAM:
c
c a)  To screen:            I've tried to make it as terse as possible.  It
c                           outputs every 100th star, but also outputs all
c                           the saturated stars.  Unfortunately, I have not
c                           tried to make all the output intelligible.  That
c                           would take too long, and it's already taken
c                           an eternity to get this to market.
c
c b)  File#1:  imname.xym:  this contains lots of the background info
c                           on how the program was called and what PSFs
c                           it used.  After that, it gives the raw position
c                           and psf-fitted instrumental magnitude
c                           (it also outputs .xymu which contains only the
c                           unsaturated stars).  The 4th column is QFIT,
c                           the quality of fit parameter, which records
c                           the fractional disagreement between the model
c                           and the image pixels.  Closer to zero is better.
c
c c)  File#2:  imname.xymc: gives the distortion-corrected position and
c                           the pixel-area corrected instrumental magnitude,
c                           followed by the raw values and the quality-of-fit
c                           estimate
c
c d)  File#3:   LOG.PERTS.fits:  will report for all images in the same
c                           fun what the perturbation PSF was.  This can
c                           be informative to look at.
c
c e)  File#4:   subtract.fits: if the 'SUBT' flag is set on the image name,
c                           then it will create an image that has the
c                           measured sources subtracted.
c
c
c Don't hesitate to contact me with questoins...
c
c      Jay Anderson
c
c      Nov 18, 2005
c      jay@eeyore.rice.edu
c      (713)348-3594
c
c---------------------------------------------------------
 
c------------------------------------------------------
c
c COMPILED PARAMETERS
c
c------------------------------------------------------
 
c------------------------------------------------------
c
c this is the image size; you should probably change it
c if your read-in image is different
c
#define _PXDIMX_     4096
#define _PXDIMY_     4096
 
c------------------------------------------------
c
c max number of stars to be used in perturbing the
c PSF
c
#define _NLMAX_      19999
 
c------------------------------------------------
c
c max # of perturbed PSFs at a time; should be set
c to "1" in this program, since we deal with only
c one image at a time
c
#define _NIMMAX_     1
#define _NIMMAXX_    1
 
 
c------------------------------------------------
c
#define _LINUX_ .true.
 
c------------------------------------------------
c
c I simply don't have time to document this program
c as much as I might like.  Rather than anticipate
c questions people may have, I'd rather just discuss the
c obvious things and let people come to me with questions.
c
c
      program img2xym_WFC
      implicit none
 
      character*200 STRING
      character*80  PROGNAME
      character*80  PSFFILE
      character*80  FILEFITS
      character*80  INFILE
      character*80  OUTPUT
      character*80  OUTPUTU
      character*80  OUTPUTC
      character*80  OUTPUTA
      character*80  DIRECT
      character*80  PREFIX
      character*7   SUFFIX
      character*80  OUTPERT
 
      integer m, oi, oj
      integer i0, j0
      integer ii, jj, ir
      integer iii,jjj
      integer  i,  j
      integer iu
      integer NSAT, NUNSAT
      logical SATD, trip
 
      real qfit
      real dx,dy
      integer hobs, fnd_hloc
      integer pobs, fnd_ploc
      integer HMIN
      real    PMAX
      real    FMIN
 
      real crat
      common /crat/ crat
      integer mmeth
      common /mmeth/ mmeth
 
      integer NARGX, NARGY, NARGU, NARGSKIP
      integer iargc, NARG, NARGs
      integer NOUT
 
      integer BITPIX
      common /BITPIX_/BITPIX
 
      real*8 xr, yr, fr, sr, mr
      real*8 dr
      real*8 xro,yro,fro,sro
      real*8 dclip
      real   z_xyoptXX
      real   apl(10)
 
      real*8 xgc, ygc
      real*8 mgc
      real   wfc3uv_zpa
 
      real*8 fro3x3, sro3x3
      real*8 fro5x5, sro5x5
 
      real ss
      real rpsf_photijk
      real HMB, histmode
      real mbar_sky
      real mbar_skyopt
      real minsky
 
c
c the inner and outer radii used to find
c sky
c
      integer irmin, irmax
      common /SKYINFO_/irmin,irmax
 
      integer hh, hhist(10)
      character*20 hhwhy(10)
      data hhwhy/
     .     '..1..not a local max',
     .     '..2..near saturation',
     .     '..3..low bad........',
     .     '..4..not iso enough.',
     .     '..5..not enough flux',
     .     '..6..near low bad...',
     .     '..7..near the edge..',
     .     '..8..near neighbor..',
     .     '..9..qfit rejected..',
     .     '.10..good!..........'/
 
 
      real apphot
 
 
      real pixo(_PXDIMX_,_PXDIMY_)
      real pixc(_PXDIMX_,_PXDIMY_)
      real pixp(_PXDIMX_,_PXDIMY_)
      real pixx(_PXDIMX_,_PXDIMY_)
      real pixh(_PXDIMX_,_PXDIMY_)
      real pixz(8412,2070)
      real pixe
 
      real plist(25)
 
      real fout, fxout, fyout
      logical inimage
      real fest, festr
 
      real zzz
      common /zzz/zzz
      real reff
      common /reff/reff
 
      logical DOQUIK
      data    DOQUIK / .false. /
 
      logical DOGSN
      data    DOGSN / .false. /
 
      logical DOAPPHOT
      data    DOAPPHOT / .false. /
      real    RAP
      integer SKI,SKO
      integer NST
 
      logical IMSUB
 
c-----------------------------------------
c
c contains info from the fits image header
c
      character*70 HDR(25)
      common/HDR/HDR
      character*70 INFO(10)
      common / fitsinfo / INFO
      integer GAIN
      real   rGAIN
      integer IFILT
 
      integer HIFLAG
      common / HIFLAG_ / HIFLAG
      data HIFLAG / 65000 /
 
      integer LOFLAG
      common / LOFLAG_ / LOFLAG
      data LOFLAG / -0250 /
 
      common/NARG_/NARG
 
c----------------------------------
c
c here are the 9x9 fiducial PSFs
c
 
c     real      psfbar(101,101)
c     real      psf09x10(101,101,09,10)
c     common /pss09x10_/psf09x10
 
 
      real psf78u(101,101,7,8)
      common /psf78_/psf78u
      real psfnim(101,101,1)
      common /psfnim_/psfnim
 
      real psfbar(101,101)
 
      integer NIMu
      common /NIMu_/NIMu
      data NIMu/0/
 
      logical DOCONST
      logical DOPERT
      logical DOQSEL
      logical DONOCR
      logical DOAPSTUDY
 
      logical APFIRST
      common /APFIRST_/APFIRST
      data    APFIRST/.true./
 
      real pertimg(501,1501)  ! can do up to 150 images...
 
 
c------------------------------------------------
c
c these lines are crucial; they tell me where the
c fiducial PSFs are specified.
c
      integer ixcl(09)
      integer iycl(10)
      common /ixycl/ixcl,iycl
      data ixcl /0000,0512,1024,1536,2168,2800,3192,3584,4096/
      data iycl /0000,0512,1024,1536,2048,2048,2560,3072,3584,4096/
 
      integer iit, jjt
      integer ninput
      data ninput/0/
 
      integer ipsf,jpsf
      real    xpsf,ypsf
      real*8  psum
 
 
      real  rr
      real  nlc(20)
      data  nlc/ 0.000, 0.000, 0.000, 0.000, 0.000,
     .           0.000, 0.000, 0.000, 0.000, 0.000,
     .           0.000, 0.000, 0.000, 0.000, 0.025,
     .           0.050, 0.075, 0.100, 0.125, 0.150/
 
      do i = 001, 0501
      do j = 001, 1501
         pertimg(i,j) = 0.
         enddo
         enddo
 
      NARGs = iargc()
 
      if (iargc().lt.5) then
         print*,'                                        '
         print*,'This routine takes 5 args:              '
         print*,'                                        '
         print*,'img2xym HMIN FMIN PMAX PSFFILE '
         print*,'                    [PERT] [SUBT] [QSEL]'
         print*,'                    [APSTUDY] [CONST] IMG.fits'
         print*,'                                        '
         print*,'1) HMIN: dominance of peak              '
         print*,'2) FMIN: min peak excess                '
         print*,'3) PMAX: max pix included               '
         print*,'                                        '
         print*,'4) PSFFILE: PSFIMAGE.FILT.fits          '
         print*,'                                        '
         print*,'               -or-                     '
         print*,'                                        '
         print*,'   "APPHOT 2.99 6 9"          '
         print*,'   ---> for aperture photy    '
         print*,'        w/ap (r<2.99) and     '
         print*,'        sky bet 6 & 9         '
         print*,'                              '
         print*,'      -or-                    '
         print*,'                              '
         print*,'   "QUICK"                    '
         print*,'   ---> for the fastest, dirtiest      '
         print*,'        possible reduction             '
         print*,'                                       '
         print*,'PERT -- before using the PSF, this flag tells the  '
         print*,'        program to tweak the PSF to better fit the '
         print*,'        PSF in the image                           '
         print*,'                                             '
         print*,'SUBT -- turning this flag on will generate a '
         print*,'        subtracted image after measuring     '
         print*,'        all the stars in each image          '
         print*,'                                             '
         print*,'QSEL -- turning this flag on will require a  '
         print*,'        peak to have a qfit of <=0.5 to make '
         print*,'        it into the final list               '
         print*,'                                             '
         print*,'APSTUDY turning this flag on will generate an '
         print*,'        additional file named "-.aps" that will'
         print*,'        report the photometry through a variety'
         print*,'        of apertures.  Useful for demonstrating'
         print*,'        that the PSF-fitting does not suffer from'
         print*,'        errors in the PSF core.               '
         print*,'                                             '
         print*,'CONST - this flag will use an "average" PSF across'
         print*,'        the chip; it is mostly usful just to show '
         print*,'        in the residual image how important it is '
         print*,'        to use the 9x10 PSF.                 '
         print*,'                                             '
         print*,'5+ IMG.fits:                         '
         print*,'   ---> sequence of .fits files to   '
         print*,'        analyze (images can be in    '
         print*,'        many different formats)      '
         print*,'                                     '
         print*,'   COMPILED PARAMETERS:              '
         print*,'                                     '
         print*,'      PXDIMX: ',_PXDIMX_
         print*,'      PXDIMY: ',_PXDIMY_
         print*,' '
         stop
         endif
 
 
      call getarg(0,PROGNAME)
 
      call getarg(1,STRING)
      read(STRING,*) HMIN
 
      call getarg(2,STRING)
      read(STRING,*) FMIN
 
      call getarg(3,STRING)
      read(STRING,*) PMAX
 
      call getarg(4,PSFFILE)
 
      print*,' '
      print*,' '
      write(*,'(''#          '')')
      write(*,'(''#          '')')
      write(*,'(''# OUTPUT FROM PROGRAM img2xym_WFC  '')')
      write(*,'(''#          '')')
      write(*,'(''#          '',a80)')  PROGNAME
      write(*,'(''#          '')')
      write(*,'(''#    HMIN: '',i3  )') HMIN
      write(*,'(''#    FMIN: '',f9.1)') FMIN
      write(*,'(''#    PMAX: '',f9.1)') PMAX
      write(*,'(''# PSFFILE: '',a80)')  PSFFILE
      write(*,'(''#  DOQUIK: '',l1)')   DOQUIK
      write(*,'(''#DOAPPHOT: '',l1)')   DOAPPHOT
      write(*,'(''#          '')')
 
 
      IMSUB  = .false.
      DOPERT = .false.
      DOQSEL = .false.
      DONOCR = .false.
      DOCONST = .false.
      DOAPSTUDY = .false.
      do NARG = 5, NARGs
         call getarg(NARG,INFILE)
         if (INFILE(1:4).eq.'SUBT')    IMSUB     = .true.
         if (INFILE(1:4).eq.'PERT')    DOPERT    = .true.
         if (INFILE(1:4).eq.'QSEL')    DOQSEL    = .true.
         if (INFILE(1:4).eq.'NOCR')    DONOCR    = .true.
         if (INFILE(1:5).eq.'CONST')   DOCONST   = .true.
         if (INFILE(1:7).eq.'APSTUDY') DOAPSTUDY = .true.
         if (INFILE(1:4).eq.'NOHI')    HIFLAG = 9999999
         enddo
 
      print*,'             '
      print*,'      IMSUB: ',IMSUB
      print*,'     DOPERT: ',DOPERT
      print*,'     DOQSEL: ',DOQSEL
      print*,'    DOCONST: ',DOCONST
      print*,'  DOAPSTUDY: ',DOAPSTUDY
      print*,'             '
 
      if (PSFFILE(1:5).eq.'QUICK') then
          NIMu = 0
          print*,'---> NIMu: ',NIMu
          DOQUIK = .true.
          goto 777
          endif
      if (PSFFILE(1:5).eq.'GSN') then
          NIMu = 0
          print*,'---> NIMu: ',NIMu
          DOGSN = .true.
          goto 777
          endif
      if (PSFFILE(1:6).eq.'APPHOT') then
          write(*,'(''# PSFFILE: '',80a)') PSFFILE
          read(PSFFILE(7:80),*) RAP,SKI,SKO
          write(*,'(''# ---> APERTURE RAP: '',f8.4)') RAP
          write(*,'(''# --->    INNER SKI: '',i3)') SKI
          write(*,'(''# --->    OUTER SKO: '',i3)') SKO
          DOAPPHOT = .true.
          NIMu = 0
          print*,'---> NIMu: ',NIMu
          goto 777
          endif
 
 
 
       NIMu = 0
       print*,' '
       write(*,'(''READ IN PSF: '',a80)') PSFFILE
       call readin_psf(PSFFILE,psf78u)
 
c
c compute an average PSF
c
       call avg_psf07x08(psf78u,psfbar)
 
c
c USE A CONSTANT PSF
c
       if (DOCONST) then
          call copy_01x01_into_07x08(psfbar,psf78u)
          print*,' '
          print*,'USING A SPATIALLY CONSTANT PSF...'
          print*,' '
          endif
 
 
       do i = 1, 101
       do j = 1, 101
          psfnim(i,j,1) = 0.
          enddo
          enddo
 
 
 777  continue
 
 
 
      NUNSAT = 0
 
      NST  = 0
 
 
c
c go through each of the command line images
c
      NARGSKIP = 0
      do NARG = 5, NARGs
 
         call getarg(NARG,INFILE)
         if (INFILE(1:4).eq.'SUBT'.or.
     .       INFILE(1:4).eq.'PERT'.or.
     .       INFILE(1:4).eq.'QSEL'.or.
     .       INFILE(1:4).eq.'NOCR'.or.
     .       INFILE(1:5).eq.'CONST'.or.
     .       INFILE(1:7).eq.'APSTUDY') then
               NARGSKIP = NARGSKIP + 1
               goto 999
               endif
 
         IFILT = 0 ! start out not knowing which filter
 
         call getarg(NARG,INFILE)
c         if (INFILE(1:1).eq.'=') then
c            IMSUB = .true.
c            INFILE = INFILE(2:80)
c            endif
         call dirstrip(INFILE,direct,prefix,suffix)
         write(*,'(''#          '')')
         write(*,'(''#  INFILE: '',80a)') INFILE
         write(*,'(''#     DIR: '',80a)') DIRECT
         write(*,'(''#     PRE: '',80a)') PREFIX
         write(*,'(''#    SUFF: '',80a)') SUFFIX
         write(*,'(''#          '')')
 
         print*,'  '
         write(*,'(''ENTER WFCREAD: '',80a)') INFILE(1:80)
         print*,'  '
         if (PREFIX(11:13).eq.'drz') then
            call readfits_r4e(INFILE,pixc,_PXDIMX_,_PXDIMY_,1)
            goto 334
            endif
         if (PREFIX(11:13).eq.'raz'.or.
     .       PREFIX(11:13).eq.'rac') then
            call readfits_r4(INFILE,pixz,8412,2070)
            do i = 0001, 2048
            do j = 0001, 2048
               pixc(i     ,j     ) = pixz(25+i+0*2103,j)
               pixc(4096-i,j     ) = pixz(25+i+1*2103,j)
               pixc(i     ,4097-j) = pixz(25+i+2*2103,j)
               pixc(4097-i,4097-j) = pixz(25+i+3*2103,j)
               enddo
               enddo
            goto 334
            endif
         call WFCREAD(INFILE,pixc)
         if (PREFIX(11:13).eq.'WJ2'.or.PREFIX(11:13).eq.'WJ3') then
            print*,' '
            print*,'Stored in clever space-saving'
            print*,'format that requires a little'
            print*,'bitty correction. '
            print*,' '
            iu = 0
            do i = 1, 4096
            do j = 1, 4096
               pixe = pixc(i,j)
               if (pixe.gt.55000) then
                  iu = iu + 1
                  pixe = 55000 + (pixe-55000)*5
                  pixc(i,j) = pixe
                  endif
               pixo(i,j) = pixc(i,j)
               enddo
               enddo
            print*,'  '
            print*,'  NUMPIX FIXED: ',iu
            print*,'  '
            endif
c        call sat_fix(pixc)
         call find_ifilt(IFILT,HDR)
 
c
c look at the gain flag to determine what the
c highest good pixel is...
c
         GAIN = 1
         print*,'INFO(10): ',INFO(10)
         do i = 1, 70
            if (INFO(10)(i:i).eq.'1') goto 330
            if (INFO(10)(i:i).eq.'2') goto 330
            enddo
         goto 333
 330     continue
         read(INFO(10),*,end=333) rGAIN
         GAIN = rGAIN + 0.1
 333     continue
         HIFLAG = 55000
c        if (GAIN.eq.2.and.BITPIX.ne.16) HIFLAG = 70000
         if (GAIN.eq.2) HIFLAG = 70000
         write(*,'(''#  BITPIX: '',i3)') BITPIX
         write(*,'(''#  HIFLAG: '',i6)') HIFLAG
         write(*,'(''#    GAIN: '',i3)') GAIN
 334     continue
 
         iu = -1
         do i = 1,76
            if (prefix(i:i).ne.' ') iu = i
            enddo
         output  = prefix(1:iu) // '.xym'
c        outputc = prefix(1:iu) // '.xymc'
         if (DOQUIK) output = prefix(1:iu) // '.xymq'
         if (DOGSN)  output = prefix(1:iu) // '.xymg'
         if (DOAPSTUDY) outputa = prefix(1:iu) // '.aps'
 
         write(*,'(''#  OUTPUT: '',80a)') OUTPUT
 
         HMB = histmode(50,_PXDIMX_-50,50,_PXDIMY_-50,pixc)
         write(*,'(''#  HISTMODE...  pixc (W/ SKYSUB): '',f8.2)') HMB
 
c        if (PMAX.gt.HIFLAG) call procsatn(pixc,HIFLAG*1.)
 
         if (.false.) then
            call writfits_r4('pixc.fits',pixc,_PXDIMX_,_PXDIMY_)
            if (.true.) stop
            endif
 
         open(77,file=output ,status='unknown')
         if ((.not.DOQUIK).and.(.not.DOGSN)) then
c            open(78,file=outputc,status='unknown') ! dist corrected
c            write(78,478)
c            write(78,477)
c            write(78,478)
  477        format('# xcor  ',1x,'  ycor  ',1x,'  mcor  ',1x,
     .              '  xraw  ',1x,'  yraw  ',1x,'  mraw  ',1x,
     .              '  qfit  ')
  478        format('--------',1x,'--------',1x,'--------',1x,
     .              '--------',1x,'--------',1x,'--------',1x,
     .              '--------')
             endif
         if (DOAPSTUDY) then
             APFIRST = .true.
             open(33,file=outputa,status='unknown') ! output aperture info
             endif
 
 
         write(77,'(''#          '')')
         write(77,'(''# OUTPUT FROM PROGRAM img2xym_WFC  '')')
         write(77,'(''#          '')')
         write(77,'(''#          '',a80)')  PROGNAME
         write(77,'(''#          '')')
         write(77,'(''#    HMIN: '',i3  )') HMIN
         write(77,'(''#    FMIN: '',f9.1)') FMIN
         write(77,'(''#    PMAX: '',f9.1)') PMAX
         write(77,'(''# PSFFILE: '',a80)')  PSFFILE
         write(77,'(''#  INFILE: '',a80)')  INFILE
         write(77,'(''#          '')')
         write(77,'(''#  BITPIX: '',i3)') BITPIX
         write(77,'(''#  HIFLAG: '',i6)') HIFLAG
         write(77,'(''#    GAIN: '',i3)') GAIN
         write(77,'(''#  DOQUIK: '',l1)') DOQUIK
         write(77,'(''#   DOGSN: '',l1)') DOGSN
         write(77,'(''#DOAPPHOT: '',l1)') DOAPPHOT
         write(77,'(''#             '')')
         write(77,'(''#        HMB: '',f9.1)') HMB
         write(77,'(''#      IFILT: '',i3)') IFILT
         write(77,'(''#             '')')
         write(77,'(''#             '')')
         write(77,'(''#      IMSUB: '',l1)') IMSUB
         write(77,'(''#     DOPERT: '',l1)') DOPERT
         write(77,'(''#     DOQSEL: '',l1)') DOQSEL
         write(77,'(''#    DOCONST: '',l1)') DOCONST
         write(77,'(''#  DOAPSTUDY: '',l1)') DOAPSTUDY
         write(77,'(''#             '')')
 
         ss = rpsf_photijk(0.,0.,2048,2048)
 
         write(77,'(''# '')')
         write(77,'(''# PSF CENTERS (x 1e5): psf78u(51,51,*,*) '')')
         write(77,'(''# '')')
         write(77,'(''#   '',9i7)') (i,i=1,7)
         write(77,'(''# '')')
         do j = 08,01,-1
            write(77,'(''# '',i3,1x,9i7)')
     .             j,(int(1e5*psf78u(51,51,i,j)),i=1,7)
            enddo
         write(77,'(''# '')')
         write(77,'(''#   '',9i7)') (i,i=1,7)
         write(77,'(''# '')')
 
         write(77,110)
         write(77,'(''# INPUT CENTRAL PSF: psf78u(*,*,5,5) '')')
         write(77,110)
         write(77,112) (i,i=36,66)
         write(77,110)
         do j = 61, 41, -1
            write(77,111) j,(int(1e5*psf78u(i,j,5,5)),i=36,66)
            enddo
         write(77,110)
         write(77,112) (i,i=36,66)
         write(77,110)
 
         write(77,110)
         write(77,'(''# INPUT CENTRAL PSF: psf78u(*,*,5,5) '')')
         write(77,110)
         write(77,122) (i,i=51-28,51+28,4)
         write(77,110)
         do j = 51+28, 51-28, -4
            write(77,121) j,(int(1e5*psf78u(i,j,5,5)),i=51-28,51+28,4)
            enddo
         write(77,110)
         write(77,122) (i,i=51-28,51+28,4)
         write(77,110)
 
 
c----------------------------------------------------------
c
c below, find the perturbation PSF, and output
c it graphically into the LOG file...
c
         if (DOPERT) then
            print*,'---> NIMu: ',NIMu
            do i = 1, 101
            do j = 1, 101
               psfnim(i,j,1) = 0.
               enddo
               enddo
            call find_psfpert(pixc,psfnim,HMB,psfbar)
            NARGU = NARG-4 - NARGSKIP
            NARGY = 1 + (NARGU-1)/10
            NARGX = NARGU - 10*(NARGY-1)
            i0 = 26 + (NARGX-1)*50
            j0 = 26 + (NARGY-1)*50
            do ii = -20,20,4
               pertimg(i0-25,j0+ii) = -0.01
               pertimg(i0+25,j0+ii) = -0.01
               pertimg(i0+ii,j0-25) = -0.01
               pertimg(i0+ii,j0+25) = -0.01
               enddo
            do i = -24, 24
            do j = -24, 24
               ii = i0 + i
               jj = j0 + j
               pertimg(i0+i,j0+j) = psfnim(51+i,51+j,1)
               xr = xr + psfnim(51+i,51+j,1)
               enddo
               enddo
            OUTPERT = 'LOG.perts.fits'
            call writfits_r4(OUTPERT,pertimg,501,j0+25)
            print*,'  psfnim(51,51,1): ',psfnim(51,51,1)
            print*,'          '
            endif
 
         write(77,110)
         write(77,'(''# PERT-PSF: psfnim(*,*,1) '')')
         write(77,110)
         write(77,112) (i,i=41,61)
         write(77,110)
         do j = 61, 41, -1
            write(77,111) j,(int(1e5*psfnim(i,j,1)),i=41,61)
            enddo
         write(77,110)
         write(77,112) (i,i=41,61)
         write(77,110)
 
 110     format('# ')
 111     format('# ',i3,1x,5i4,8i5,6i6,8i5,4i4)
 112     format('# ',3x,1x,5i4,8i5,6i6,8i5,4i4)
 
 121     format('# ',i3,1x,15i6)
 122     format('# ',3x,1x,15i6)
 
         write(*,'(''# '')')
         write(*,'(''# PSFCEN  PRE: '',f8.5)') psfbar(51,51)
         write(*,'(''# PSFCEN PERT: '',f8.5)') psfnim(51,51,1)
         write(*,'(''# PSFCEN POST: '',f8.5)') psfbar(51,51)+
     .                                         psfnim(51,51,1)
         write(*,'(''# '')')
 
         write(77,'(''# '')')
         write(77,'(''# PSFCEN  PRE: '',f8.5)') psfbar(51,51)
         write(77,'(''# PSFCEN PERT: '',f8.5)') psfnim(51,51,1)
         write(77,'(''# PSFCEN POST: '',f8.5)') psfbar(51,51)+
     .                                          psfnim(51,51,1)
         write(77,'(''# '')')
 
         write(77,'(''# PSF ENCLOSED ENERGY... '')')
         write(77,'(''# '')')
         do ir = 1, 12
            psum = 0.
            do ipsf = 001, 101
            do jpsf = 001, 101
               xpsf = 0.25*(ipsf-51)
               ypsf = 0.25*(jpsf-51)
               if (xpsf**2+ypsf**2.lt.ir**2) then
                  psum = psum + psfbar(ipsf,jpsf)
                  endif
               enddo
               enddo
            write(77,'(''# '',i2.2,1x,f10.7)') ir,psum/16.0
            enddo
         write(77,'(''# '')')
         write(77,'(''# '')')
 
 
c-------------------------------------------------------
c
c now, go through and measure the saturated stars; this new
c version of it takes the PSF into account, and is most accurately
c cone when PERT is used, since it allows account to be taken of
c how flux sloshes from the core to the outskirts, which without
c correction would result in a different zeropoint for saturated
c and unsaturated stars, since unsaturated stars use 5x5-pixel
c apertures, and saturated one generally use larger ones.
c
         do i = 0001, 4096
         do j = 0001, 4096
            pixo(i,j) = pixc(i,j)
            enddo
            enddo
         call peak_sat(pixc,pixp,pixx,HIFLAG*1.0)
         call sat_phot(pixc,pixp,pixx)
         call find_ifilt(IFILT,HDR)
 
c        call writfits_r4('pixo.fits',pixo,4096,4096)
c        call writfits_r4('pixp.fits',pixp,4096,4096)
c        call writfits_r4('pixx.fits',pixx,4096,4096)
c        call writfits_r4('pixc.fits',pixc,4096,4096)
 
 
 
 
 
 
 
 
 
 
 
 
c-------------------------------------------------------
c
c pixh will hold the subtracted image, in case the user
c wants it output
c
         do i = 001, 4096
         do j = 001, 4096
            pixh(i,j) = min(pixc(i,j),HIFLAG*1.)
            enddo
            enddo
 
         print*,'  '
         print*,'  '
         print*,'*************************************'
         print*,'*************************************'
         print*,'**   '
         print*,'**  GO THRU IMAGE AND FIND PEAKS... '
         print*,'**   '
         print*,'**   PMAX: ',PMAX
         print*,'**   NIMu: ',NIMu
         print*,'**   '
         print*,'*************************************'
         print*,'*************************************'
         print*,'  '
         print*,'  '
 
         do hh = 1, 10
            hhist(hh) = 0
            enddo
 
         print*,'                      '
         print*,'  KEY TO PIXEL INFO:  '
         print*,'                      '
         print*,'    P08   P04    P06  '
         print*,'    P03   PCEN   P02  '
         print*,'    P07   P05    P09  '
         print*,'                      '
 
         write( *,378)
         write( *,278)
         write( *,378)
 
c--------------------------------------
c go through the image pixel by pixel to
c see if each one qualifies as a source worth
c "finding" and measuring
c
         NOUT = 0
         do jjt = 5, _PXDIMY_-5
         do iit = 5, _PXDIMX_-5
 
 
            hh = 1
            ii = iit
            jj = jjt
            if (pixc(ii,jj).lt.1) goto 444
 
            xr = 0.
            yr = 0.
            mr = 0.
            fr = 0.
            sr = 0.
            if (jj.ge.2046.and.jj.le.2050) goto 444
            qfit  = 0.00
            mmeth = 0
            irmin = 0
            irmax = 0
 
c-------------------------
c skip if not a local max
c
            hh = 1
            do i = ii-1,ii+1
            do j = jj-1,jj+1
               if (pixc(i,j).gt.pixc(ii,jj)) goto 444
               enddo
               enddo
 
c-------------------------
c skip if too close to bigger saturated pixels
c
            hh = 2
            if (pixc(ii,jj).ge.HIFLAG) then
               do i = ii-2,ii+2
               do j = jj-3,jj+3
                  if ((i.ne.ii.or.j.ne.jj).and.
     .                pixc(i,j).gt.pixc(ii,jj)) goto 444
                  enddo
                  enddo
               endif
 
c-------------------------
c skip if too close to LOBAD pixel
c
            hh = 3
            if (pixc(ii,jj).lt.LOFLAG) goto 444
 
c-------------------------
c skip if too close to brighter neighbor
c
            hh = 4
            if (pixc(ii,jj).gt.PMAX) goto 444
            if (pixc(ii,jj).gt.PMAX.and.PMAX.lt.HIFLAG) goto 444
            hobs = fnd_hloc(ii,jj,pixc)
            if (hobs.lt.abs(HMIN)) goto 444
 
c-------------------------
c skip if not enough flux in this peak
c
            hh = 5
            FEST = pixc(ii,jj) +
     .             max(pixc(ii+1,jj)+pixc(ii+1,jj+1)+pixc(ii,jj+1),
     .                 pixc(ii+1,jj)+pixc(ii+1,jj-1)+pixc(ii,jj-1),
     .                 pixc(ii-1,jj)+pixc(ii-1,jj+1)+pixc(ii,jj+1),
     .                 pixc(ii-1,jj)+pixc(ii-1,jj-1)+pixc(ii,jj-1))
            minsky = min(
     .          (pixc(ii+1,jj+1)+pixc(ii  ,jj+1)+pixc(ii-1,jj+1))/3,
     .          (pixc(ii+1,jj-1)+pixc(ii  ,jj-1)+pixc(ii-1,jj-1))/3,
     .          (pixc(ii+1,jj+1)+pixc(ii+1,jj  )+pixc(ii+1,jj-1))/3,
     .          (pixc(ii-1,jj+1)+pixc(ii-1,jj  )+pixc(ii-1,jj-1))/3)
            FESTR = FEST
            if (pixc(ii,jj).gt.HIFLAG) minsky = HMB
            FEST  = FEST - 4*minsky
            if (pixc(ii,jj).lt.HIFLAG.and.FEST.lt.FMIN) goto 444
 
c-------------------------
c skip if not enough flux in this peak
c
            hh = 6
            NSAT = 0
            do i = ii-3,ii+3
            do j = jj-3,jj+3
               if (pixc(ii,jj).lt.HIFLAG.and.
     .             pixc(i,j).le.LOFLAG) goto 444
               if (pixc(i,j).gt.HIFLAG) NSAT = NSAT + 1
               enddo
               enddo
 
c-------------------------
c skip if too close to image edge
c
            hh = 7
            sr = min(HMB,minsky)
            if (pixc(ii,jj).gt.HIFLAG) sr = HMB
            do m = 1, 25
               if (.not.inimage(ii+oi(m),jj+oj(m))) goto 444
               plist(m) = pixc(ii+oi(m),jj+oj(m))
               enddo
 
c-------------------------
c found peak worth measuring!
c
 
 
c-------------------------
c
c centroid positions, 3x3 aperture flux (quick-n-easy)
c
 
            fxout = (pixc(ii+1,jj)-pixc(ii-1,jj))/2/
     .              (pixc(ii  ,jj)-min(pixc(ii+1,jj),pixc(ii-1,jj)))
            fyout = (pixc(ii,jj+1)-pixc(ii,jj-1))/2/
     .              (pixc(ii,jj  )-min(pixc(ii,jj+1),pixc(ii,jj-1)))
            fout  = pixc(ii+1,jj+1)+pixc(ii  ,jj+1)+pixc(ii-1,jj+1)+
     .              pixc(ii+1,jj  )+pixc(ii  ,jj  )+pixc(ii-1,jj  )+
     .              pixc(ii+1,jj-1)+pixc(ii  ,jj-1)+pixc(ii-1,jj-1)-
     .              9*sr
            SATD = .false.
 
            if (DOQUIK) then
               xr = ii + fxout
               yr = jj + fyout
               fr = FEST
               sr = minsky
               goto 137 ! if quick reduction, we're done!
               endif
 
            if (DOGSN) then
               xr = ii
               yr = jj
               sr = mbar_sky(ii,jj,10,15,pixc)
               call gsn2Dfit_5x5(xr,yr,fr,dr,sr,pixc)
               qfit = 0.
               goto 137 ! if quick reduction, we're done!
               endif
 
            if (fout.lt.1) fout = FEST
            if (.not.(fxout.lt. 1.0)) fxout = 0.
            if (.not.(fxout.gt.-1.0)) fxout = 0.
            if (.not.(fyout.lt. 1.0)) fyout = 0.
            if (.not.(fyout.gt.-1.0)) fyout = 0.
            if (.not.( fout.gt.-1.0))  fout = 1.
 
            if (pixc(ii,jj).gt.HIFLAG) then
               if (pixc(ii,jj).lt.199999) goto 444
               SATD = .true.
               xr = ii
               yr = jj
               fr = pixc(ii,jj)
               mr = -2.5*log10(fr)
               if (jj.gt.2048) then ! top chip nonlinear
                  rr = -mr
                  if (rr.gt.19.99) rr = 19.99
                  if (rr.lt.01.01) rr =  1.01
                  ir = int(-mr)
                  if (ir.gt.14) then
                     mr = mr - (nlc(ir)+(rr-ir)*(nlc(ir+1)-nlc(ir)))
                     fr = 10**(-mr/2.5)
                     endif
                  endif
               mmeth = 9
               goto 137
               endif
 
            if (DOAPPHOT) then
               if (RAP.le.0) then  ! backdoor... easiest possible thing to do,
                  fr = fout        ! measure the brightest 4/9 pixels...
                  xr = ii + fxout
                  yr = jj + fyout
                  goto 137
                  endif
               xr = ii + fxout
               yr = jj + fyout
               fr = apphot(ii,jj,RAP,SKI,SKO,pixc,ss)
               sr = ss
               if (DONOCR) then
                  do iii = ii-1, ii+1
                  do jjj = jj-1, jj+1
                     if (pixc(iii,jjj)-sr.lt.
     .                   0.08*(pixc(ii,jj)-sr))goto 444
                     enddo
                     enddo
                  endif
               goto 137
               endif
 
            xro = ii + fxout
            yro = jj + fyout
 
c-------------------------------------
c
c crude fit
c
            call find_xyz5x5(xro,yro,fro,sro,qfit,pixc)
            fro5x5 = fro
            sro5x5 = sro
 
            fro3x3 = 0.
            sro3x3 = 0.
 
            xr = xro
            yr = yro
            fr = fro5x5
            sr = sro5x5
            if (fr.lt.1) then
               call find_xyz3x3(xro,yro,fro,sro,qfit,pixc)
               fro3x3 = fro
               sro3x3 = sro
               fr = fro3x3
               sr = sro3x3
               endif
 
c-----------------------------------------------------
c find a good sky, using the closest decent pixels
c (this does subtract off a PSF first)
c
c           sr = mbar_sky(ii,jj,10,15,pixc)
 
            sr = mbar_skyopt(xro,yro,pixc)
 
            if (DONOCR) then
               do iii = ii-1, ii+1
               do jjj = jj-1, jj+1
                  if (pixc(iii,jjj)-sr.lt.0.08*(pixc(ii,jj)-sr))goto 444
                  enddo
                  enddo
               endif
 
            hh = 8
            if (HMIN.lt.0.and.pixc(ii,jj).lt.HIFLAG) then
                ss = sr
                pobs = fnd_ploc(ii,jj,pixc,ss)
                if (pobs.lt.abs(HMIN)) goto 444
                endif
 
 
c-----------------------------------------------------
c find my best postion, then flux
c
            xr = xro
            yr = yro
            call find_xyzXX(xr,yr,fr,sr,qfit,pixc)
            fr = z_xyoptXX(xr,yr,sr,pixc,2.5)
 
            hh = 9
            if (DOQSEL.and.(.not.SATD)) then
               if (qfit.lt.0.0) goto 444
               if (qfit.gt.0.5) goto 444
               endif
 
            hh = 10
 
            if (fr.lt.1) then  ! if the "best" flux is found
               xr = xro        ! to be negative, then go with something
               yr = yro        ! cruder; we know it's a peak, so it *should*
               fr = fro        ! have a positve flux
               sr = sro
               mmeth = 0
               endif
 
  137       continue
            if (.not.(fr.gt.1.00)) fr = 1.0
            if (.not.(fr.lt.1e19)) fr = 1.0
            mr = -2.5*log10(fr)
            NOUT = NOUT + 1
            if (NOUT.eq.NOUT/2500*2500) then
               write( *,378)
               write( *,278)
               write( *,378)
               endif
            if (NOUT.eq.NOUT/100*100.or.NOUT.lt.100.or.SATD)
     .          write( *,178) NOUT,xr,yr,mr,sr,mmeth,irmin,irmax,
     .                        (int(plist(m)),m=1,9),
     .                        int(minsky),int(FESTR),int(FEST)
 
 
            if (DOQUIK)
     .         write(77,'(f8.3,1x,f8.3,1x,f8.3,1x,i8)')
     .            dclip(xr,0000.0d0,9999.0d0),
     .            dclip(yr,0000.0d0,9999.0d0),
     .            dclip(mr, -99.9d0,   0.0d0)
            if (DOGSN)
     .         write(77,'(f8.3,1x,f8.3,1x,f8.3,1x,f8.1,1x,i8)')
     .            dclip(xr,0000.0d0,9999.0d0),
     .            dclip(yr,0000.0d0,9999.0d0),
     .            dclip(mr, -99.9d0,   0.0d0),
     .            dclip(sr, -99.9d0,9999.0d0),
     .            int(pixc(ii,jj)+0.5)
            if ((.not.DOQUIK).and.(.not.DOGSN)) then
               write(77,377)
     .               dclip(xr,0000.0d0,9999.0d0),
     .               dclip(yr,0000.0d0,9999.0d0),
     .               dclip(mr, -99.9d0,   0.0d0),
     .               dclip(qfit* 1.0d0,-0.9d0,   9.9d0),
     .               pixc(ii,jj)+0.5
               call wfc3uv_gc(xr,yr,xgc,ygc,0)
               mgc = mr + wfc3uv_zpa(xr,yr)
c              write(78,177)
c    .               dclip(xgc,-999.0d0,9999.0d0),
c    .               dclip(ygc,-999.0d0,9999.0d0),
c    .               dclip(mgc, -99.9d0,   0.0d0),
c    .               dclip(xr ,0000.0d0,9999.0d0),
c    .               dclip(yr ,0000.0d0,9999.0d0),
c    .               dclip(mr , -99.9d0,   0.0d0),qfit
               endif
 377        format(f8.3,1x,f8.3,1x,f8.4,1x,f8.5,1x,f10.2,1x,f10.2)
 177        format(f8.3,1x,f8.3,1x,f8.4,1x,
     .             f8.3,1x,f8.3,1x,f8.4,1x,f8.5)
 178        format(i7,1x,f8.3,1x,f8.3,1x,f7.3,1x,f6.1,'|',i1,i2,i2,
     .                      '|',i7,'|',4i6,'|',4i6,'|',3i8,1x,1l1)
 278        format('  NOUT ',1x,    ! i7
     .             '  XCEN  ',1x,   ! f8.3
     .             '  YCEN  ',1x,   ! f8.3
     .             '  MAG  ',1x,    ! f7.3
     .             '  SKY ' ,'|', ! f6.1
     .             'MiSoS'  ,'|', ! i1,1x,i1,1x,i2
     .             '   PCEN','|',   ! i7
     .             '   P02',        ! i6
     .             '   P03',
     .             '   P04',
     .             '   P05','|',
     .             '   P06',
     .             '   P07',
     .             '   P08',
     .             '   P09','|',
     .             '    mSKY',
     .             '    F2x2',
     .             ' F2x2-SK')
 378        format('-------',1x,'--------',1x,'--------',1x,
     .             '-------',1x,'------','|','- - -','|',
     .             '-------','|',
     .             ' ----- ----- ----- -----','|',
     .             ' ----- ----- ----- -----','|',
     .             ' ------- ------- -------')
 
 
            if (.not.SATD.and.DOAPSTUDY) then
               call apstudy(ii,jj,xr,yr,fr,sr,pixc,apl)
               endif
 
            if (IMSUB.and.(.not.SATD)) then
               do i = max(1,ii-10),min(_PXDIMX_,ii+10)
               do j = max(1,jj-10),min(_PXDIMX_,jj+10)
                  dx = i-xr
                  dy = j-yr
                  if (sqrt(dx**2+dy**2).lt.10) then
                     pixe = fr*rpsf_photijk(dx,dy,ii,jj)
                     if (pixh(i,j).lt.HIFLAG)
     .                   pixh(i,j) = pixh(i,j) - pixe
                     endif
                  enddo
                  enddo
               endif
 
c            do i = max(1,ii-5),min(_PXDIMX_,ii+5)
c            do j = max(1,jj-5),min(_PXDIMX_,jj+5)
c               write(99,199) NARG-4,i,j,int(pixc(i,j)+0.5),xr,yr,fr,sr
c 199           format(i3,1x,i4,1x,i4,1x,i10,1x,
c     .                f8.2,1x,f8.2,1x,f12.0,1x,f9.3)
c               enddo
c               enddo
 
            pixc(ii,jj) = pixc(ii,jj)+1
 444        continue
            hhist(hh) = hhist(hh) + 1
 123        format(f9.3,1x,f9.3,1x,f8.4,1x,f9.1,1x,f8.2,5x,
     .             f9.3,1x,f9.3,1x,f8.4,1x,f9.1,1x,f8.2,1x,i6,1x,f6.2)
            enddo
            enddo
 
c----------------------------------------
c
c report on why some pixels weren't identified
c as stars
c
         do hh = 1, 10
            write( *,277) hh,hhist(hh),hhist(hh)/40.96/4096,hhwhy(hh)
            write(77,277) hh,hhist(hh),hhist(hh)/40.96/4096,hhwhy(hh)
 277        format('#  ',i2,1x,i8,1x,f8.3,'%',3x,a20)
            enddo
 
         close(77)
c        if (.not.DOQUIK) close(78)
 
 
c-----------------------------------------
c
c generate a subtracted image, if asked to
c
         if (IMSUB) then
            FILEFITS = prefix(1:iu) // '_SUB.fits'
            call writfits_r4(FILEFITS,pixh,_PXDIMX_,_PXDIMY_)
            endif
 
 999     continue
         enddo! NARG
 
       end
 
 
 
 
 
 
 
c---------------------------------------------
c
c given a pixel (i,j), how far out in radius do
c you have to go to find a brighter pixel?
c
      integer function fnd_hloc(i,j,pix)
      implicit none
 
      integer i, j
      real pix(_PXDIMX_,_PXDIMY_)
 
      integer m, h
      integer ooi, ooj
      integer  oi,  oj
      real    pcen
      integer ii, jj
      logical inimage
 
      pcen = pix(i,j)
      do m = 2, 251
         ooi = oi(m)
         ooj = oj(m)
         h = sqrt((ooi**2+ooj**2)*1.0)
         ii = i + ooi
         jj = j + ooj
         if (.not.inimage(ii,jj)) goto 1
         if (pix(ii,jj).gt.pix(i,j)) goto 1
         enddo
 
 1    continue
      fnd_hloc = h
      return
      end
 
 
 
c------------------------------------------------
c
c This is similar to fnd_hloc, but it is more stringent.
c It not only requires this pixel to be the brightest within
c a certain radius, but it also requires it to be MUCH
c brighter than its surrouding pixels.  If there are any
c "suspicious" pixels nearby, it reports them as crowding
c neighbors.
c
      integer function fnd_ploc(i,j,pix,SKY)
      implicit none
 
      integer i, j
      real pix(_PXDIMX_,_PXDIMY_)
      real    SKY
 
      integer m, h
      integer ooi, ooj
      integer  oi,  oj
      real    pcen
      integer ii, jj
      logical inimage
 
 
 
      pcen = pix(i,j)
      do m = 2, 251
         ooi = oi(m)
         ooj = oj(m)
         h = sqrt((ooi**2+ooj**2)*1.0)
         ii = i + ooi
         jj = j + ooj
         if (.not.inimage(ii,jj)) goto 1
         if (pix(ii,jj).gt.pix(i,j)) goto 1
         if (m.gt.09.and.pix(ii,jj)-SKY.gt.0.50*(pcen-SKY)) goto 1
         if (m.gt.21.and.pix(ii,jj)-SKY.gt.0.25*(pcen-SKY)) goto 1
         if (m.gt.37.and.pix(ii,jj)-SKY.gt.0.10*(pcen-SKY)) goto 1
         if (m.gt.45.and.pix(ii,jj)-SKY.gt.0.05*(pcen-SKY)) goto 1
         enddo
 
 1    continue
      fnd_ploc = h
      return
      end
 
 
 
c----------------------------------------
c
c Returns "true" if a given pixels is within
c the confines of the image.
c
      logical function inimage(i,j)
      implicit none
 
      integer i, j
 
      inimage = .true.
 
      if (i.lt.       1) inimage = .false.
      if (i.gt._PXDIMX_) inimage = .false.
      if (j.lt.       1) inimage = .false.
      if (j.gt._PXDIMY_) inimage = .false.
 
      return
      end
 
 
c--------------------------------------
c
c bubble-sorts a real*4 list into ascending order
c
      subroutine rbubble(r1,NTOT)
      implicit none
      real r1(1)
      real temp
      integer NTOT
      integer n
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
 
c--------------------------------------
c
c bubble-sorts a real*8 list into ascending order
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
 
 
c-------------------------------------------------
c
c this routine uses a histogram-method to find the mode
c of a large number of pixels, those within (i1:i2,j1:j2)
c of the array pixarr (only use pixels with values that
c are between -10 and 4990)
c
      real function histmode(i1,i2,j1,j2,pixarr)
      implicit none
 
      integer i1,i2,j1,j2
      real pixarr(_PXDIMX_,_PXDIMY_)
 
      integer i, j, h
      integer mhi, mlo
      real hist_summ
 
      integer hist(5000)
      integer hcum(5000)
      integer NTOT
      real root
 
      real    pcum(101)
      integer ipc
      integer imin
      real    dmin
      real    ptot
 
      integer HIFLAG
      common /HIFLAG_/HIFLAG
      integer LOFLAG
      common /LOFLAG_/LOFLAG
 
      do h = 1,5000
         hist(h) = 0
         hcum(h) = 0
         enddo
 
      NTOT = 0
      PTOT = 0
      do i = i1,i2
      do j = j1,j2
         if (pixarr(i,j).ge.HIFLAG) goto 2
         if (pixarr(i,j).le.LOFLAG) goto 2
         h = pixarr(i,j)+10.5
         if (h.lt.   1) goto 2
         if (h.gt.5000) goto 2
         hist(h) = hist(h) + 1
         NTOT = NTOT + 1
         PTOT = PTOT + pixarr(i,j)
 2       continue
         enddo
         enddo
 
 
      hcum(1) = hist(1)/2
      do h = 2, 5000
         hcum(h) = hcum(h-1) + (hist(h)+hist(h-1))/2
         ipc = 100.*hcum(h)/NTOT + 1.5
         do i = ipc, 101
            pcum(i) = h
            enddo
         enddo
      do ipc = 1, 101  ! goes from 0% to 100%
         pcum(ipc) = root(NTOT*(ipc-1)/100.,hcum,5000)
         enddo
 
      dmin = pcum(40+1)-pcum(1)
      imin = 1
      do i = 2, 101-40
         if (pcum(40+i)-pcum(i).lt.dmin) then
            dmin = pcum(40+i)-pcum(i)
            imin = i
            endif
         enddo
 
      mhi = pcum(40+imin)+3
      mlo = pcum(imin)   -3
      print*,'--->  MLO ',mlo
      print*,'--->  MHI ',mhi
 
      histmode = 0.
      hist_summ = 0.
      do h = mlo, mhi
         hist_summ = hist_summ + hist(h)
         histmode = histmode + hist(h)*(h-10)
 4       format(1x,i4,1x,i4,1x,i8,1x,f8.4,' %  ',
     .          1x,i8,1x,f8.4,' %  ')
         enddo
 
      histmode = histmode / hist_summ
      print*,'---> HISTMODE: ',HISTMODE
 
      if (hist_summ.le.0) histmode = PTOT / NTOT
      print*,'--->    MLO-10: ',MLO-10
      print*,'--->    MHI-10: ',MHI-10
 
      NTOT = 0
      PTOT = 0
      do i = i1,i2
      do j = j1,j2
         if (pixarr(i,j).ge.(mlo-10).and.
     .       pixarr(i,j).le.(mhi-10)) then
             NTOT = NTOT + 1
             PTOT = PTOT + pixarr(i,j)
             endif
         enddo
         enddo
      histmode = PTOT/NTOT
      print*,'         PTOT: ',PTOT
      print*,'         NTOT: ',NTOT
      print*,'---> HISTMODE: ',HISTMODE
 
      return
 
      end
 
 
c--------------------------------------
c
c return the value of n where list(n) = v ;
c list must be monotonic increasing
c
      real function root(v,list,NTOT)
      implicit none
      real    v
      integer list(1)
      integer NTOT
      integer nl
 
      do nl = 1, NTOT-2
         if (v.lt.list(nl+1)) goto 3
         enddo
 3    continue
 
      root = nl
      if (list(nl).ne.list(nl+1))
     .    root = nl + 1.*(v-list(nl))/(list(nl+1)-list(nl))
 
      return
      end
 
 
 
c---------------------------------------------------
c
c takes a filename and splits it off into directory,
c prefix, and suffix
c
      subroutine dirstrip(FULL,DIRECT,PREFIX,SUFFIX)
      implicit none
 
 
      character*80 FULL
      character*80 DIRECT
      character*80 PREFIX
      character*7  SUFFIX
 
      integer i
 
      integer ndir
      integer ndot
      integer nbrk
      integer nend
 
      ndir = 0
      ndot = 0
      nbrk = 0
      nend = 0
      do i = 1,80
         if (FULL(i:i).eq.'/')  ndir = i
         if (FULL(i:i).eq.'.')  ndot = i
         if (FULL(i:i).eq.'[')  nbrk = i
         if (FULL(i:i).eq.' '.and.
     .       nend.eq.0)  nend = i-1
         enddo
 
      if (ndot.eq.0) then
         print*,' '
         print*,'DIRSTRIP: '
         print*,' '
         print*,'THE IMAGE MUST HAVE A DOT IN ITS NAME.'
         print*,'NAME: ',FULL
         print*,' '
         stop
         endif
 
 
      DIRECT = FULL(1:ndir)
      PREFIX = FULL(ndir+1:ndot-1)
      SUFFIX = FULL(ndot+1:nend  )
 
      return
      end
 
c---------------------------------------
c
c Returns rval if rval is between rlo and rhi,
c otherwise it returns the closer endpoint.  It
c also make sure rval is not a NaN
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
 
 
 
c-------------------------------------------
c
c same as rclip, but for double precision numbers
c
      real*8 function dclip(dval,dlo,dhi)
      implicit none
 
      real*8 dval
      real*8 dlo
      real*8 dhi
 
      dclip = dval
      if (dclip.gt.dhi) dclip = dhi
      if (dclip.lt.dlo) dclip = dlo
      if (.not.(dclip.lt.dhi).and.
     .    .not.(dclip.gt.dlo)) dclip = dhi
 
      return
      end
 
 
 
c---------------------------------------------------
c
c This routine will find an iteratively clipped mean sky
c value from an annulis about a point.  The point is (ixc,iyc).
c It will take the pixels between r=inr and r=ior and will
c analyze them for a robust average.
c
      real function mbar_sky(ixc,iyc,inr,ior,pixarr)
      implicit none
 
      integer ixc, iyc
      integer inr, ior
      real pixarr(_PXDIMX_,_PXDIMY_)
 
      integer i  , j
      integer ixh, iyh
      real    rij
      integer nuse
 
      logical inimage
      real sklist(99999)
      integer nsk
      real bar, sig
 
      integer HIFLAG
      common / HIFLAG_ / HIFLAG
 
      integer LOFLAG
      common / LOFLAG_ / LOFLAG
 
 
      nsk = 0
      do i = -ior+1, ior-1
      do j = -ior+1, ior-1
         rij = i**2+j**2
         if (rij.ge.ior**2) goto 333
         if (rij.lt.inr**2) goto 333
         ixh = ixc + i
         iyh = iyc + j
         if (.not.inimage(ixh,iyh)) goto 333
         if (pixarr(ixh,iyh).le.LOFLAG) goto 333
         if (pixarr(ixh,iyh).ge.HIFLAG) goto 333
         nsk = nsk + 1
         if (nsk.gt.99999) goto 333
         sklist(nsk) = pixarr(ixh,iyh)
 333     continue
         enddo
         enddo
 
      call rbarsigs(sklist,nsk,bar,sig,nuse,1.75)
 
      mbar_sky = bar
 
      return
      end
 
 
 
c----------------------------------------
c
c this routine will do simple aperture photometry
c about a pixel (i,j)
c
 
      real function apphot(i,j,rap,ri,ro,pix,sky)
      implicit none
 
      integer i, j
      real    rap
      integer ri, ro
      real pix(_PXDIMX_,_PXDIMY_)
      real    sky
 
 
      integer ii, jj
      real mbar_sky
      integer ntot
      real    ftot
 
      sky = mbar_sky(i,j,ri,ro,pix)
 
      ntot = 0
      ftot = 0
      do ii = i-int(rap+0.5),i+int(rap+0.5)
      do jj = j-int(rap+0.5),j+int(rap+0.5)
         if (ii.lt.      01) goto 1
         if (jj.lt.      01) goto 1
         if (ii.gt._PXDIMX_) goto 1
         if (jj.gt._PXDIMY_) goto 1
         if ((ii-i)**2+(jj-j)**2.gt.rap**2) goto 1
         ftot = ftot + pix(ii,jj)-sky
         ntot = ntot + 1
 1       continue
         enddo
         enddo
 
      apphot = ftot
 
      return
      end
 
 
c
c this routine will assume a sky and will fit
c a gaussian to the inner 5x5 pixels of a star
c
 
      subroutine gsn2Dfit_5x5(x,y,z,d,s,pix)
      implicit none
 
      real*8 x,y,z,d,s
      real pix(_PXDIMX_,_PXDIMY_)
 
      integer i,ii,icen,i0
      integer j,jj,jcen,j0
      real ptot, psum
      real xtot, ytot
 
      i0 = int(x+0.5)
      j0 = int(y+0.5)
 
      z = 0
      d = 0
 
      if (i0.lt.6) return
      if (j0.lt.6) return
      if (i0.gt._PXDIMX_-5) return
      if (j0.gt._PXDIMY_-5) return
 
c      print*,' '
c      print*,' x: ',x
c      print*,' y: ',y
c      print*,' s: ',s
c      print*,' p: ',pix(i0,j0)
c      print*,' '
c      print*,' '
c
c      print*,' '
c      write(*,110) (i,i=i0-5,i0+5)
c      print*,' '
c      do j = j0+5,j0-5,-1
c         write(*,111) j,(int(pix(i,j)+0.5),i=i0-5,i0+5)
c  110    format(5x,1x,11i10)
c  111    format(i5,1x,11i10)
c         enddo
c      print*,' '
c      write(*,110) (i,i=i0-5,i0+5)
c      print*,' '
 
      ptot = 0
      icen = i0
      jcen = j0
      do i = i0-1,i0+1
      do j = j0-1,j0+1
         psum = 0
         do ii = i-2, i+2
         do jj = j-2, j+2
            psum = psum + pix(ii,jj)
            enddo
            enddo
         if (psum.gt.ptot) then
            ptot = psum
            icen = i
            jcen = j
            endif
         enddo
         enddo
 
      z = ptot - 25*s
 
      xtot = 0.
      ytot = 0.
      do i = icen-2, icen+2
      do j = jcen-2, jcen+2
         xtot = xtot + (pix(i,j)-s)*(i-icen)
         ytot = ytot + (pix(i,j)-s)*(j-jcen)
         enddo
         enddo
      x = icen + xtot/ptot
      y = jcen + ytot/ptot
 
c      print*,' x: ',x
c      print*,' y: ',y
c      print*,' z: ',z
c      print*,' s: ',s
c      print*,' d: ',d
 
      return
      end
 
 
 
 
 
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
 
 
      write(*,'(''ENTER readfits_i2: '',80a)') FILEU
 
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
         call buff2pix_i2(buffb,ibuff,0001,1440)
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
      stop
 
      end
 
 
 
c--------------------------------------------------------
c
c This routine reads in images stored in integer*4 format...
c
c to use:
c
c   input parameters:
c
c         FILE     is the filename
c         NX,NY    are the dimensions of the image to be
c                  read in
c
c
c  output parameters:
c
c         pix( , ) is an integer*4 image of dimensions NX,NY
c
c
c
      subroutine readfits_i4r(FILE,pix,NX,NY)
      implicit none
 
      integer NX, NY
      character*80 FILE
      real*4 pix(NX,NY)
 
      character*70 INFO(10)
      common / fitsinfo / INFO
 
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
      integer ii, jj
      integer n
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffb,buffc)
 
      integer*4 ibuff(720)
 
      integer ifirst, i1, i2
      integer j
 
      integer np1, np2, npt
      integer nextend
      integer nread
      real bscale, bzero
      integer bitpix
 
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
 
      if (DIAG) then
         print*,'enter readfits...'
         print*,'FILE: ',FILE(1:60)
         endif
 
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
      do i = 1, 23
         HDR(i) = ' '
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
         if (field.eq.'FILTNAM1'.and.stream(1:1).eq.'F')
     .                            INFO(02) = stream
         if (field.eq.'FILTNAM2'.and.stream(1:1).eq.'F')
     .                            INFO(02) = stream
         if (field.eq.'FILTER1 '.and.stream(1:1).eq.'F')
     .                            INFO(02) = stream
         if (field.eq.'FILTER2 '.and.stream(1:1).eq.'F')
     .                            INFO(02) = stream
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
 
 
      if (laxis(1).ne.NX.or.laxis(2).ne.NY) then
         print*,'FITS image not the expected dimensions: '
         print*,'    input NX: ',nx
         print*,'    input NY: ',ny
         print*,'              '
         print*,'    laxis(1): ',laxis(1)
         print*,'    laxis(2): ',laxis(2)
         print*,'  '
         stop
         endif
 
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
 
      if (BITPIX.ne.32) then
         print*,'                '
         print*,'readfits_i4...: '
         print*,'                '
         print*,'   you called a routine to read in a'
         print*,'   long i*4 image, and the image you'
         print*,'   gave it has BITPIX = ',BITPIX
         print*,'   '
         print*,'   FILE: ',FILE
         print*,'   '
         stop
         endif
 
      nbper = 4*laxis(1)*laxis(2)
      npt   =   laxis(1)*laxis(2)
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
      do i = i1, i2, 1
         read(10,rec=i,iostat=ios) buffc
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/4 + 1
         np2 = (nbyteE-nbyte1)/4 + 1
         np2 = min(np2,npt)
         call buff2pix_i4(buffb,ibuff,0001,0720)
         do n = np1, np2, 1
            jj = n/NX   + 1
            ii = n-NX*(jj-1)
            pix(ii,jj) = ibuff(n-np1+1)*bscale+bzero
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
      print*,'READFITS_I4R ERROR'
      print*,' FILEU: ',FILEU
      stop
 
      end
 
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
 
      character*8  field
      character*20 stream
 
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
         stream = buffc(k*80+11:k*80+31)
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
         print*,'prob'
         print*,'---> readfits_r4 called and BITPIX.ne.-32...'
         print*,'---> BITPIX = ',BITPIX
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
      print*,'READFITS ERROR: '
      print*,'  FILE : ',FILE
      print*,'  FILEU: ',FILEU
      stop
 
      end
 
 
 
c
c
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
 
      FILEU = FILE
      do i = 75,2,-1
         if (FILE(i:i+4).eq.'.fits') FILEU = FILE(1:i+4)
         enddo
 
      if (DIAG) then
         print*,'enter readfits_r4...'
         print*,'FILE: ',FILE(1:60)
         endif
 
      open(10,file=FILE,status='old',
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
 
 
 
 
 
 
      subroutine writfits_r4(FILE,pix,PXDIMX,PXDIMY)
      implicit none
 
      character*80 FILE
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
 
      character*80 FILEU
      character*70 HDR(25)
      common/HDR/HDR
 
      FILEU = FILE
      do i = 75,2,-1
         if (FILE(i:i+4).eq.'.fits') FILEU = FILE(1:i+4)
         enddo
 
c      write(*,'(''writfits_r4:  '',i5,''x'',i5,1x,80a)')
c     .            PXDIMX,PXDIMY,FILEU
 
      open(10,file=FILEU,status='unknown',
     .     err=900,recl=2880,form='UNFORMATTED',
     .     access='DIRECT')
 
 
      i = 1
      write(buffc( 0*80+1: 1*80),'(''SIMPLE  = T      '')')
      write(buffc( 1*80+1: 2*80),'(''BITPIX  = -32    '')')
      write(buffc( 2*80+1: 3*80),'(''NAXIS   ='',i12)') 2
      write(buffc( 3*80+1: 4*80),'(''NAXIS1  ='',i12)') PXDIMX
      write(buffc( 4*80+1: 5*80),'(''NAXIS2  ='',i12)') PXDIMY
      write(buffc( 5*80+1: 6*80),'(''DATATYPE='',9a)')
     .                          ' ''REAL*4'''
      write(buffc( 6*80+1: 7*80),'(''DATE    ='',11a)')
     .                          ' ''28/01/00'''
      write(buffc(07*80+1:08*80),'(''COMMENT  '',a05)') '     '
      write(buffc(08*80+1:09*80),'(''COMMENT  '',a05)') '     '
      write(buffc(09*80+1:10*80),'(''CRPIX1  ='',a20)') HDR(01)
      write(buffc(10*80+1:11*80),'(''CRPIX2  ='',a20)') HDR(02)
      write(buffc(11*80+1:12*80),'(''CRVAL1  ='',a20)') HDR(03)
      write(buffc(12*80+1:13*80),'(''CRVAL2  ='',a20)') HDR(04)
      write(buffc(13*80+1:14*80),'(''CTYPE1  ='',a20)') HDR(05)
      write(buffc(14*80+1:15*80),'(''CTYPE2  ='',a20)') HDR(06)
      write(buffc(15*80+1:16*80),'(''CD1_1   ='',a20)') HDR(07)
      write(buffc(16*80+1:17*80),'(''CD1_2   ='',a20)') HDR(08)
      write(buffc(17*80+1:18*80),'(''CD2_1   ='',a20)') HDR(09)
      write(buffc(18*80+1:19*80),'(''CD2_2   ='',a20)') HDR(10)
      write(buffc(19*80+1:20*80),'(''ORIENTAT='',a20)') HDR(11)
      write(buffc(20*80+1:21*80),'(''PA_APER ='',a20)') HDR(12)
      write(buffc(21*80+1:22*80),'(''PA_V3   ='',a20)') HDR(13)
      write(buffc(22*80+1:23*80),'(''DATE-OBS='',a20)') HDR(14)
      write(buffc(23*80+1:24*80),'(''TIME-OBS='',a20)') HDR(15)
      write(buffc(24*80+1:25*80),'(''EXPTIME ='',a20)') HDR(16)
      write(buffc(25*80+1:26*80),'(''ROOTNAME='',a20)') HDR(17)
      write(buffc(26*80+1:27*80),'(''TARGNAME='',a20)') HDR(18)
      write(buffc(27*80+1:28*80),'(''RA_TARG ='',a20)') HDR(19)
      write(buffc(28*80+1:29*80),'(''DEC_TARG='',a20)') HDR(20)
      write(buffc(29*80+1:30*80),'(''PROPOSID='',a20)') HDR(21)
      write(buffc(30*80+1:31*80),'(''FILTER1 ='',a20)') HDR(22)
      write(buffc(31*80+1:32*80),'(''FILTER2 ='',a20)') HDR(23)
      write(buffc(33*80+1:34*80),'(''VAFACTOR='',a20)') HDR(24)
      write(buffc(32*80+1:33*80),'(''CCDGAIN ='',a20)') HDR(25)
      write(buffc(33*80+1:34*80),'(''COMMENT  '',a05)') '     '
      write(buffc(34*80+1:35*80),'(''COMMENT  '',a05)') '     '
      write(buffc(35*80+1:36*80),'(''END      '')')
 
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
 
 
      subroutine WFCFLTREAD(FILE,pix)
      implicit none
      real pix(4096,4096)
      character*80 FILE
 
      call readfits_WFC(FILE,pix(0001,0001),1)
      call readfits_WFC(FILE,pix(0001,2049),4)
 
      return
      end
 
 
 
      subroutine readfits_WFC(FILE,pix,nimg)
      implicit none
 
      character*80 FILE
      real pix(4096,2048)
 
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
 
      if (laxis(2).eq.2051) laxis(2) = 2048
      if (laxis(1).ne.4096.or.
     .    laxis(2).ne.2048) then
          print*,'  laxis1: ',laxis(1)
          print*,'  laxis2: ',laxis(2)
          print*,'    4096: ',4096
          print*,'    2048: ',2048
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
         call buff2pix_r4(buffb,pix,np1,npt)
         enddo
 
      return
 
 900  continue
      print*,'readfits_WFC: READFITS ERROR'
      stop
 
      end
 
 
c-----------------------------------------------------------------
c
c NOTE that if the image size is larger than 4096x4096, only
c a 4096x4096 image will be read in.  I did want to make it
c able to do aperture photometry on the drz images.
c
c
      subroutine WFCDRZREAD(FILENAME,pix)
      implicit none
 
      character*80 FILENAME
      real pix(4096,4096)
 
      real pixlist(720)
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
      integer j
 
      character*70 INFO(10)
      common / fitsinfo / INFO
 
      character*2880 buffc
      byte buffb(2880)
      equivalence (buffc,buffb)
 
      integer ifirst, i1, i2
 
      integer np1, np2, npt
      integer nextend
      integer nread
      real bscale, bzero
      integer bitpix
 
      real    i2r
 
      integer ii, nn, nx, ny
      real EXPTIME
 
      logical LINUX
      data    LINUX/.true./
 
 
      print*,'WFCDRZREAD -- FILENAME: ',FILENAME
 
      open(10,file=FILENAME,status='old',
     .     err=900,recl=2880,form='UNFORMATTED',
     .     access='DIRECT')
 
      bscale = 1.0
      bzero  = 0.0
 
      naxes = -1
      laxis(1) = 1
      laxis(2) = 1
      laxis(3) = 1
      nextend = 0
 
      EXPTIME = 1.0
 
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
         if (field.eq.'EXPTIME ') read(stream,*) EXPTIME
         if (field.eq.'FILTER1 '.and.
     .       stream(1:1).ne.'C'.and.stream(2:2).ne.'C')
     .       read(stream,*) INFO(2)
         if (field.eq.'FILTER2 '.and.
     .       stream(1:1).ne.'C'.and.stream(2:2).ne.'C')
     .       read(stream,*) INFO(2)
         if (field.eq.'END     ') goto 101
 
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
 
         enddo
      goto 100
 101  continue
 
      nread = nread + 1
      if (abs(BITPIX).ne.32) goto 100
 
      ifirst = i+1
      i1 = i
      i2 = i
 
      if (laxis(1).ne.4096.or.
     .    laxis(2).ne.4096) then
          print*,'                       '
          print*,'***************        '
          print*,'**  WARNING  **        '
          print*,'***************        '
          print*,'                       '
          print*,'The image:             '
          print*,'                       '
          write(*,'(10x,80a)') FILENAME
          print*,'                       '
          print*,'has dimensions:        '
          print*,'                       '
          print*,'     laxis1: ',laxis(1)
          print*,'     laxis2: ',laxis(2)
          print*,'                       '
          print*,'This routine is designed to deal with '
          print*,'  4096x4096 images.  But it will      '
          print*,'  read the image into as much of      '
          print*,'  the 4096x4096 array as possible, and'
          print*,'  will operate on that.               '
          print*,'                       '
          print*,'                       '
          print*,'***************        '
          print*,'                       '
          endif
 
      nbper  = 4*laxis(1)*laxis(2)
      npt    = laxis(1)*laxis(2)
      nbyte1 = 1
      nbyte2 = nbper
      i1 = i+1 + nbyte1/2880
      i2 = i+1 + nbyte2/2880
 
      do i = i1, i2, 1
         read(10,rec=i,iostat=ios) buffc
         nbyte0 = (i-ifirst)*2880+   1
         nbyteE = (i-ifirst)*2880+2880
         np1 = (nbyte0-nbyte1)/4 + 1
         np2 = (nbyteE-nbyte1)/4 + 1
         call buff2pix_r4(buffb,pixlist,001,720)
         do ii = 001, 720
            nn = np1 + (ii-1)
            ny = 1 + nn/laxis(1)
            nx = nn-(ny-1)*laxis(1)
            if (nx.ge.000.and.nx.le.4096.and.
     .          ny.ge.000.and.ny.le.4096) then
                pix(nx,ny) = pixlist(ii)
                endif
            enddo
         enddo
 
c      print*,' '
c      print*,'MULTIPLY BY EXPTIME...'
c      print*,' '
c      print*,'   EXPTIME = ',EXPTIME
c      print*,' '
c      do i = 0001, 4096
c      do j = 0001, 4096
c         pix(i,j) = pix(i,j)*EXPTIME
c         enddo
c         enddo
 
      return
 
 900  continue
      write(*,*)
      write(*,'(''FILE NOT FOUND: '')')
      write(*,'(80a)') FILENAME
      write(*,*)
      stop
 
      end
 
 
 
c---------------------------------------
c
c This routine will read in a generic WFC image
c (either a _W, a _Wi4, a _flt, or a _drz image)
c into a 4096x4096 array.
c
c
      subroutine WFCREAD(INFILE,pix)
      implicit none
 
      character*80 INFILE
      real pix(4096,4096)
 
      character*20 stream
      integer BITPIX
      integer NEXTEND
      logical interior
      integer i,j
 
      common /BITPIX_/BITPIX
 
      print*,'                 '
      print*,'ENTER WFCREAD... ',INFILE(1:40)
      print*,'                 '
 
      call query_hdr(INFILE,'BITPIX  ',stream)
      read(stream,*) BITPIX
      stream = '0'
      call query_hdr(INFILE,'NEXTEND ',stream)
      read(stream,*) NEXTEND
 
      print*,' '
      print*,'    BITPIX : ',BITPIX
      print*,'    NEXTEND: ',NEXTEND
      print*,' '
c
c _W.fits image (int*2)
c
      if (BITPIX.eq.+16.and.NEXTEND.eq.0) then
         write(*,'(''INFILE_W:   readfits_i2r:  '',80a)') INFILE
         call readfits_i2r(infile,pix,4096,4096)
         return
         endif
 
c
c _Wi4.fits image (int*4)
c
      if (BITPIX.eq.+32) then
         write(*,'(''INFILE_Wi4: readfits_i2r:  '',80a)') INFILE
         call readfits_i4r(infile,pix,4096,4096)
         return
         endif
 
c
c _drz.fits image (real*4)
c
      if (BITPIX.eq.16.and.NEXTEND.eq.3) then
         write(*,'(''INFILE_drz:   WFCDRZREAD:  '',80a)') INFILE
         call WFCDRZREAD(INFILE,pix)
         print*,' '
         print*,'FIX THE ZEROS IN THE DRZ...'
         do i = 0006, 4090
         if (i.eq.i/512*512) write(*,'(4x,''i: '',i4.4)') i
         do j = 0006, 4090
            if (pix(i,j).eq.0.00) then
                pix(i,j) = -999.9
                if (interior(i,j,pix)) pix(i,j) = 99999.0
                endif
            enddo
            enddo
         return
         endif
 
c
c _flt.fits image (real*4 in extensions 1 and 4...)
c
      if (BITPIX.eq.+16.and.NEXTEND.eq.6) then
         write(*,'(''INFILE_flt:   WFCFLTREAD:  '',80a)') INFILE
         call WFCFLTREAD(INFILE,pix)
         return
         endif
 
      if (BITPIX.eq.-32) then
         write(*,'(''INFILE:      readfits_r4:  '',80a)') INFILE
         call readfits_r4(INFILE,pix,4096,4096)
         return
         endif
 
 
 
      print*,'UNABLE TO (WFC)READ IN IMAGE: ',INFILE
      stop
      end
 
 
 
 
c----------------------------------------------------------
c
c the drizzled images have "0.000" values for (1) saturated
c   pixels and (2) pixels at the edges of the field that aren't
c   mapped into the actual image.  I want the #1 pixels to
c   have a HIBAD falue, but not the #2 pixels.  So this routine
c   will determine if given pixel is interior, or not.
c
      logical function interior(i,j,pix)
      implicit none
 
      integer i,j
      real pix(4096,4096)
      integer jc
      integer k
 
      interior = .false.
 
      jc = 2055 + 137*(i-0110)/(4050-0110)
      if (abs(j-jc).lt.30) return
 
      do k = i-1,0001,-1
         if (pix(k,j).gt.0) goto 11
         enddo
      return
 
 11   continue
      do k = i+1,4096,+1
         if (pix(k,j).gt.0) goto 12
         enddo
      return
 
 12   continue
      do k = j+1,0001,-1
         if (pix(i,k).gt.0) goto 13
         enddo
      return
 
 13   continue
      do k = j+1,4096,+1
         if (pix(i,k).gt.0) goto 14
         enddo
      return
 
 14   continue
      interior = .true.
      return
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
 
 
 
c#include "/home/jayander/FORTRAN/CODE/ROUTINES/WFC.09x10/iclip.f"
 
c--------------------------------
c
c take the psf in an image format (as it is
c read in from the fits file) and convert it into
c an array format (which is how I use it)
c
      subroutine psfimg_to_psf09x10(psfimg,psf09x10)
      implicit none
      real psf09x10(101,101,09,10)
      real psfimg(901,1001)
 
      integer ir, jr
      integer ipsf
      integer jpsf
 
 
      do ir = 01, 09
      do jr = 01, 10
         do ipsf = 001, 101
         do jpsf = 001, 101
            psf09x10(ipsf,jpsf,ir,jr) =
     .         psfimg(ipsf+(ir-1)*100,
     .                jpsf+(jr-1)*100)
            enddo
            enddo
         enddo
         enddo
 
      return
      end
 
 
 
c--------------------------------
c
c reads in a PSF image file
c
      subroutine readin_psfimage(psfimage,FITSFILE)
      implicit none
      real psfimage(901,1001)
      character*80 FITSFILE
      call readfits_r4(FITSFILE,psfimage,901,1001)
      return
      end
 
 
c------------------------------------------------
c
c find the average PSF across the chip
c
      subroutine avg_psf07x08(psf07x08,psfbar)
      implicit none
      real psf07x08(101,101,07,08)
      real psfbar(101,101)
      integer ipsf, jpsf
      integer ir, jr
 
      do ipsf = 001, 101
      do jpsf = 001, 101
         psfbar(ipsf,jpsf) = 0.
         do ir = 1, 07
         do jr = 1, 08
            psfbar(ipsf,jpsf) = psfbar(ipsf,jpsf) +
     .                          psf07x08(ipsf,jpsf,ir,jr)
            enddo
            enddo
         psfbar(ipsf,jpsf) = psfbar(ipsf,jpsf)/7/8
         enddo
         enddo
 
      return
      end
 
c------------------------------
c
c
      subroutine copy_01x01_into_07x08(psf,psf07x08)
      implicit none
      real psf(101,101)
      real psf07x08(101,101,07,08)
      integer ir  , jr
      integer ipsf, jpsf
 
      do ir = 01, 07
      do jr = 01, 08
         do ipsf = 001, 101
         do jpsf = 001, 101
            psf07x08(ipsf,jpsf,ir,jr) = psf(ipsf,jpsf)
            enddo
            enddo
         enddo
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
 
 
 
c-------------------------------------------
c
c this will just call rpsf_photij; it used to
c be that this routine added the perturbation, but
c now I add the perturbation within rpsf_photij
c
      real function rpsf_photijk(dx,dy,i,j)
      implicit none
      real dx, dy
      integer i, j
      real pertpsf(51,51)
      common /pertpsf_/pertpsf
 
      real rpsf_phot
      real rpsf_photij
      real supp
 
c
c this is the perturbative part
c
      real psfnim(101,101,_NIMMAXX_)
      common /psfnim_/psfnim
 
      integer NIMu
      common /NIMu_/NIMu
 
 
      rpsf_photijk = rpsf_photij(dx,dy,i,j)
 
 
c------------------------------------------------------------
c
c if we're using a perturbed PSF, add the perturbation here;
c
      if (NIMu.ne.0.and.sqrt(dx**2+dy**2).lt.6.0) then
         supp = rpsf_phot(dx,dy,psfnim(1,1,NIMu))
         rpsf_photijk = rpsf_photijk + supp
         endif
 
      return
      end
 
 
c----------------------------------------------
c
c uses a chisq minimization technique to find the optimal
c values for (x,y,z) using the inner 5x5 pixels; it
c takes the sky value (SKIN) as a given.  It does a very
c crude grid search to find the chisq minimim.
c
      subroutine find_xyzXX(x,y,z,SKIN,qfit,pix)
      implicit none
 
      real*8 x,y,z
 
      real*8 SKIN
 
      real pix(_PXDIMX_,_PXDIMY_)
      real qfit
 
      real psfx(5,5)
      real pixx(5,5)
      real sig2(5,5)
      real psfx_(5,5)
 
      integer mmeth
      common / mmeth / mmeth
 
      integer NIT
      integer ii, jj, i, j
 
      real rpsf_photijk
      real ftot
      real ptot
      real ftotw
      real ptotw
 
      real fxu, fx0, fxout
      real fyu, fy0, fyout
      real  zu,  zout
      real  dx, dy
 
      real  dd(9)
      data dd /0.100,0.050,0.030,0.020,
     .         0.010,0.005,0.003,0.002,
     .         0.001/
 
      real dxn(9)
      real dyn(9)
      data dxn/0.000,+1.000,-1.000,+0.000,+0.000,
     .               +0.707,+0.707,-0.707,-0.707/
      data dyn/0.000,+0.000,+0.000,+1.000,-1.000,
     .               +0.707,-0.707,-0.707,+0.707/
 
      integer nnn
 
      real err, chimin, chitot
      real abstot, absmin
      real pixtot
      real pixval
 
      integer NDONE
      data NDONE/0/
 
      mmeth = 5
 
      ii = int(x+0.5)
      jj = int(y+0.5)
 
      pixtot = 0.
      do i = 1, 5
      do j = 1, 5
         pixval    = pix(ii+i-3,jj+j-3)
         pixx(i,j) = pixval - SKIN
         pixtot = pixtot + pixx(i,j)
         sig2(i,j) = max(pixval,0.00) + 25.00
         enddo
         enddo
 
 
 
      fxout = x-ii
      fyout = y-jj
      do NIT = 01, 09
         fx0 = fxout
         fy0 = fyout
         chimin = 9e30
         do nnn = 01, 09
            fxu = fx0 + dd(NIT)*dxn(nnn)
            fyu = fy0 + dd(NIT)*dyn(nnn)
            ftot  = 0.
            ptot  = 0.
            ftotw = 0.
            ptotw = 0.
            do i = 1, 5
            do j = 1, 5
               dx = (i-3) - fxu
               dy = (j-3) - fyu
               psfx(i,j) = rpsf_photijk(dx,dy,ii,jj)
               ftot  = ftot  + psfx(i,j)
               ptot  = ptot  + pixx(i,j)
               ftotw = ftotw + psfx(i,j)*psfx(i,j)/sig2(i,j)
               ptotw = ptotw + pixx(i,j)*psfx(i,j)/sig2(i,j)
               enddo
               enddo
            zu = ptot/ftot
            chitot = 0.
            abstot = 0.
            do i = 1, 5
            do j = 1, 5
               err  = abs(pixx(i,j)-zu*psfx(i,j))
               chitot = chitot + err**2/sig2(i,j)
               abstot = abstot + err
               enddo
               enddo
            if (chitot.lt.chimin) then
               chimin = chitot
               fxout = fxu
               fyout = fyu
               zout  =  zu
               absmin = abstot/pixtot
               do i = 1, 5
               do j = 1, 5
                  psfx_(i,j) = psfx(i,j)
                  enddo
                  enddo
               endif
            enddo!nnn
         if (NIT.ge.2.and.zout.le.  100) goto 5
         if (NIT.ge.4.and.zout.le. 1000) goto 5
         if (NIT.ge.6.and.zout.le.10000) goto 5
         enddo!NIT
 5    continue
 
      x = ii + fxout
      y = jj + fyout
      z = zout
 
      qfit = sqrt(chitot/(25-3))
      qfit = absmin
      if (qfit.gt.9.99) qfit = 9.99
 
      return
      end
 
 
c----------------------------------------------------
c
c this routine will fit to the center using only
c stars within the fitting radius, and will taper
c out to the annulus; this is useful in the psf-pert
c routine, since you'd like to be as centered as possible
c on the very core.
c
c    find_xyzXX_chi_CENFOCUS
c
      subroutine find_xyzXX_CENFOCUS(x,y,z,SKIN,qfit,pix,RINR,ROTR)
      implicit none
 
      real*8 x,y,z
 
      real*8 SKIN
 
      real pix(_PXDIMX_,_PXDIMY_)
      real qfit
      real RINR,ROTR
 
      integer mmeth
      common / mmeth / mmeth
 
      integer NIT
      integer ii, jj, i, j
 
      real rpsf_photijk
      real ftot
      real ptot
      real fxu, fx0, fxout
      real fyu, fy0, fyout
      real  zu,       zout
      real  dx, dy, dr
 
      real  dd(9)
      data dd /0.100,0.050,0.030,0.020,
     .         0.010,0.005,0.003,0.002,
     .         0.001/
 
      real dxn(9)
      real dyn(9)
      data dxn/0.000,+1.000,-1.000,+0.000,+0.000,
     .               +0.707,+0.707,-0.707,-0.707/
      data dyn/0.000,+0.000,+0.000,+1.000,-1.000,
     .               +0.707,-0.707,-0.707,+0.707/
 
      integer nnn
 
      real err, chimin, chitot
      real abstot, absmin
 
      real pixval
      real psfx(5,5)
      real pixx(5,5)
      real sig2(5,5)
      real wwww(5,5),wtu
      real wx, wy
 
      real sumpw
      real sumfw
      real sumw
      real p1,p2,p3,p4,pmx
      real pval, fval
 
      integer NDONE
      data NDONE/0/
 
      mmeth = 5
 
      ii = int(x+0.5)
      jj = int(y+0.5)
 
      do i = 1, 5
      do j = 1, 5
         pixval    = pix(ii+i-3,jj+j-3)
         pixx(i,j) = pixval - SKIN
         sig2(i,j) = max(pixval,0.00) + 25.00
         wwww(i,j) = 0.
         enddo
         enddo
 
      p1 = pixx(4,3)+pixx(4,4)+pixx(3,4)
      p2 = pixx(2,3)+pixx(2,4)+pixx(3,4)
      p3 = pixx(2,3)+pixx(2,2)+pixx(3,2)
      p4 = pixx(4,3)+pixx(4,2)+pixx(3,2)
      pmx = max(p1,p2,p3,p4)
 
      wwww(3,3) = 1
      if (p1.eq.pmx) then
         wwww(4,3) = 1.00
         wwww(4,4) = 1.00
         wwww(3,4) = 1.00
         endif
      if (p2.eq.pmx) then
         wwww(2,3) = 1.00
         wwww(2,4) = 1.00
         wwww(3,4) = 1.00
         endif
      if (p3.eq.pmx) then
         wwww(2,3) = 1.00
         wwww(2,2) = 1.00
         wwww(3,2) = 1.00
         endif
      if (p4.eq.pmx) then
         wwww(4,3) = 1.00
         wwww(4,2) = 1.00
         wwww(3,2) = 1.00
         endif
 
 
      fxout = x-ii
      fyout = y-jj
      do NIT = 01, 09
         fx0 = fxout
         fy0 = fyout
         chimin = 9e30
         do nnn = 01, 09
            fxu = fx0 + dd(NIT)*dxn(nnn)
            fyu = fy0 + dd(NIT)*dyn(nnn)
            sumpw = 0.
            sumfw = 0.
            sumw  = 0.
            do i = 1, 5
            do j = 1, 5
               dx = (i-3) - fxu
               dy = (j-3) - fyu
               dr = sqrt(dx**2+dy**2)
               wtu = (ROTR-dr)/(ROTR-RINR)
               if (wtu.lt.0.0) wtu = 0.0
               if (wtu.gt.1.0) wtu = 1.0
               wx = 1.5-abs(dx)
               if (wx.lt.0.0) wx = 0.0
               if (wx.gt.1.0) wx = 1.0
               wy = 1.5-abs(dy)
               if (wy.lt.0.0) wy = 0.0
               if (wy.gt.1.0) wy = 1.0
               psfx(i,j) = rpsf_photijk(dx,dy,ii,jj)
               enddo
               enddo
            sumpw = 0.
            sumfw = 0.
            sumw  = 0.
            do i = 1, 5
            do j = 1, 5
               sumpw = sumpw + wwww(i,j)*pixx(i,j)
               sumfw = sumfw + wwww(i,j)*psfx(i,j)
               sumw  = sumw  + wwww(i,j)
               enddo
               enddo
            zu = sumpw/sumfw
            chitot = 0.
            abstot = 0.
            do i = 1, 5
            do j = 1, 5
               err  = abs(pixx(i,j)-zu*psfx(i,j))
               chitot = chitot + wwww(i,j)*err**2
               abstot = abstot + wwww(i,j)*err
               enddo
               enddo
            chitot = chitot/sumw
            if (chitot.lt.chimin) then
               chimin = chitot
               fxout = fxu
               fyout = fyu
               zout  =  zu
               absmin = abstot/sumpw
               endif
            enddo!nnn
         if (NIT.ge.2.and.zout.le.  100) goto 5
         if (NIT.ge.4.and.zout.le. 1000) goto 5
         if (NIT.ge.6.and.zout.le.10000) goto 5
         enddo!NIT
 5    continue
 
      x = ii + fxout
      y = jj + fyout
      z = zout
 
      qfit = sqrt(chitot/(25-3))
      qfit = absmin
      if (qfit.gt.9.99) qfit = 9.99
 
c
c find an aperture flux
c
      ptot = 0.
      ftot = 0.
      do i = -4, 4
      do j = -4, 4
         if (i**2+j**2.le.4.5**2) then
            dx = ii+i-x
            dy = jj+j-y
            ptot = ptot + pix(ii+i,jj+j)-SKIN
            ftot = ftot + rpsf_photijk(dx,dy,ii,jj)
            endif
         enddo
         enddo
      z = ptot/ftot
 
 
c
c an estimate of the quality of fit
c
      qfit = 0.
      do i = -4, 4
      do j = -4, 4
         if (i**2+j**2.le.4.5**2) then
            dx = ii+i-x
            dy = jj+j-y
            pval = pix(ii+i,jj+j)-SKIN
            fval = rpsf_photijk(dx,dy,ii,jj)
            err = (pval-z*fval)
            qfit = qfit + abs(err)/z
            endif
         enddo
         enddo
 
      return
      end
 
 
c
c this has been modified to return the flux
c within the inner 3x3 pixels (unless they
c agree to better than 2%)
c
 
      subroutine find_xyz3x3(x,y,z,s,qfit,pix)
      implicit none
 
      real*8 x,y,z,s
      real   qfit
 
      real pix(_PXDIMX_,_PXDIMY_)
 
      integer mmeth
      common / mmeth / mmeth
 
      integer irmin, irmax
      common / SKYINFO_ / irmin,irmax
 
      integer i, ii
      integer j, jj
      real dx, dy
 
      real flxout, skyout
      real rpsf_photijk
      real r
 
      real psf_list(9)
      real pix_list(9)
 
      integer n
 
      ii = int(x+0.5)
      jj = int(y+0.5)
 
      n = 0
      do i = ii-1, ii+1
      do j = jj-1, jj+1
         n = n + 1
         dx = x-i
         dy = y-j
         psf_list(n) = rpsf_photijk(dx,dy,ii,jj)
         pix_list(n) = pix(i,j)
         enddo
         enddo
 
      call lsq_star_fit_3x3(psf_list,pix_list,flxout,skyout,r)
 
      irmin = 0
      irmax = 0
      mmeth = 0
 
      z = flxout
      s = skyout
 
c      print*,'  3x3:  x: ',x
c      print*,'  3x3:  y: ',y
c      print*,'  3x3:  z: ',z
c      print*,'  3x3:  s: ',s
 
      qfit = 0.999
      return
      end
 
 
      subroutine lsq_star_fit_3x3(d,p,flx,sky,r)
      implicit none
 
      real d(9)
      real p(9)
      real flx, sky
      real r
 
      call LINFIT_3x3(d,p,9,sky,flx,r)
 
      return
      end
 
 
 
 
      SUBROUTINE LINFIT_3x3(X,Y,NPTS,A,B,R)
      INTEGER NPTS
      REAL SUM, SUMX, SUMY, SUMX2, SUMXY, SUMY2
      REAL DELTA
      REAL X(1), Y(1)
      real r
 
   11 SUM = 0
      SUMX = 0
      SUMY = 0
      SUMX2 = 0
      SUMXY = 0
      SUMY2 = 0
   21 DO 59  I = 1, NPTS
           XI = X(I)
           YI = Y(I)
           SUM   = SUM   + 1
           SUMX  = SUMX  + XI
           SUMY  = SUMY  + YI
           SUMX2 = SUMX2 + XI*XI
           SUMXY = SUMXY + XI*YI
           SUMY2 = SUMY2 + YI*YI
   59 CONTINUE
C
C                CALCULATE COEFFICIENTS AND STANDARD DEVIATIONS
C
 
   51 DELTA = SUM*SUMX2 - SUMX*SUMX
      A = (SUMX2*SUMY - SUMX*SUMXY) / DELTA
   53 B = (SUMXY*SUM  - SUMX*SUMY ) / DELTA
      R = 0.00
   71 IF (DELTA*(SUM*SUMY2 - SUMY*SUMY).GE.0) THEN
          R = (SUM*SUMXY - SUMX*SUMY) /
     C        SQRT(DELTA*(SUM*SUMY2 - SUMY*SUMY))
          ENDIF
 
 
      RETURN
      END
 
 
 
c
c this has been modified to return the flux
c within the inner 5x5 pixels (unless they
c agree to better than 2%)
c
 
      subroutine find_xyz5x5(x,y,z,s,qfit,pix)
      implicit none
 
      real*8 x,y,z,s
      real   qfit
 
      real pix(_PXDIMX_,_PXDIMY_)
 
 
      integer mmeth
      common / mmeth / mmeth
 
      integer irmin, irmax
      common / SKYINFO_ / irmin,irmax
 
      integer i, ii
      integer j, jj
      real dx, dy
 
      real flxout, skyout
      real rpsf_photijk
      real r
 
      real psf_list(25)
      real pix_list(25)
 
      integer n, NTOT
 
      ii = int(x+0.5)
      jj = int(y+0.5)
 
      n = 0
      do i = ii-2, ii+2
      do j = jj-2, jj+2
         if (pix(i,j).le.pix(ii,jj)) then
            n = n + 1
            dx = x-i
            dy = y-j
            psf_list(n) = rpsf_photijk(dx,dy,ii,jj)
            pix_list(n) = pix(i,j)
            endif
         enddo
         enddo
      NTOT = n
 
      call LINFIT(psf_list,pix_list,NTOT,skyout,flxout,r)
 
      irmin = 0
      irmax = 0
      mmeth = 0
 
      z = flxout
      s = skyout
 
      qfit = 0.999
      return
      end
 
 
 
      SUBROUTINE LINFIT(X,Y,NPTS,A,B,R)
      INTEGER NPTS
      REAL SUM, SUMX, SUMY, SUMX2, SUMXY, SUMY2
      REAL DELTA
      REAL X(1), Y(1)
      real r
 
   11 SUM = 0
      SUMX = 0
      SUMY = 0
      SUMX2 = 0
      SUMXY = 0
      SUMY2 = 0
   21 DO 59  I = 1, NPTS
           XI = X(I)
           YI = Y(I)
           SUM   = SUM   + 1
           SUMX  = SUMX  + XI
           SUMY  = SUMY  + YI
           SUMX2 = SUMX2 + XI*XI
           SUMXY = SUMXY + XI*YI
           SUMY2 = SUMY2 + YI*YI
   59 CONTINUE
C
C                CALCULATE COEFFICIENTS AND STANDARD DEVIATIONS
C
 
   51 DELTA = SUM*SUMX2 - SUMX*SUMX
      A = (SUMX2*SUMY - SUMX*SUMXY) / DELTA
   53 B = (SUMXY*SUM  - SUMX*SUMY ) / DELTA
      R = 0.00
   71 IF (DELTA*(SUM*SUMY2 - SUMY*SUMY).GE.0) THEN
          R = (SUM*SUMXY - SUMX*SUMY) /
     C        SQRT(DELTA*(SUM*SUMY2 - SUMY*SUMY))
          ENDIF
 
 
      RETURN
      END
 
 
 
 
 
c-----------------------------------------
c
c this will measure the sky *and* remove the
c star profile before it does; this way it
c can measure a sky near a star which still
c may be very slightly affected by the star
c itself
c
      real function mbar_skyopt(xc,yc,pix)
      implicit none
 
      real*8  xc, yc
      real pix(_PXDIMX_,_PXDIMY_)
 
      real bar1, sig1
 
      integer ixc, iyc
      integer i, j
      real dx, dy
 
      real FCEN, PCEN
      real ZZZ
      common/zzz/zzz
      real SKY
      real PIXSUM
      real PSFSUM
      integer NNNSUM
      integer n1, n1u
      integer n2
      real sl1(1000)
      logical inimage
      real rpsf_photij
      real reff
      common /reff/reff
      integer irmin, irmax, irnxt
      common /SKYINFO_/irmin, irmax
 
 
      integer HIFLAG
      common / HIFLAG_ / HIFLAG
 
      integer LOFLAG
      common / LOFLAG_ / LOFLAG
 
      ixc = int(xc+0.5)
      iyc = int(yc+0.5)
 
      PCEN = 0.
      FCEN = 0.
      do i = -2, 2
      do j = -2, 2
         dx = ixc+i-xc
         dy = iyc+j-yc
         PCEN = PCEN + rpsf_photij(dx,dy,ixc,iyc)
         FCEN = FCEN + pix(ixc+i,iyc+j)
         enddo
         enddo
 
c      print*,' '
c      print*,' IXC: ',IXC
c      print*,' IYC: ',IYC
c      print*,' '
c      print*,'pCEN: ',pix(ixc,iyc)
c      print*,'FCEN: ',FCEN
c      print*,'PCEN: ',PCEN
c      print*,'ZCEN: ',FCEN/PCEN
c      print*,' '
c      print*,' '
 
      irmin = 2
 5    continue
      irmax = irmin+1
      PIXSUM = 0.
      PSFSUM = 0.
      NNNSUM = 0
      do i = -irmax, irmax
      do j = -irmax, irmax
         dx = ixc+i-xc
         dy = iyc+j-yc
         if (i**2+j**2.ge.(irmin-0.5)**2.and.
     .       i**2+j**2.lt.(irmax+0.5)**2.and.
     .       inimage(ixc+i,iyc+j)) then
                PIXSUM = PIXSUM + pix(ixc+i,iyc+j)
                PSFSUM = PSFSUM + rpsf_photij(dx,dy,ixc,iyc)
                NNNSUM = NNNSUM + 1
                endif
         enddo
         enddo
 
      ZZZ = 0.0
      SKY = (PIXSUM-ZZZ*PSFSUM)/NNNSUM
      ZZZ = (FCEN-25*SKY)/PCEN
      if (ZZZ.lt.0) ZZZ = 0.
      SKY = (PIXSUM-ZZZ*PSFSUM)/NNNSUM
      ZZZ = (FCEN-25*SKY)/PCEN
      if (ZZZ.lt.0) ZZZ = 0.
      SKY = (PIXSUM-ZZZ*PSFSUM)/NNNSUM
      ZZZ = (FCEN-25*SKY)/PCEN
      if (ZZZ.lt.0) ZZZ = 0.
      SKY = (PIXSUM-ZZZ*PSFSUM)/NNNSUM
      ZZZ = (FCEN-25*SKY)/PCEN
      if (ZZZ.lt.0) ZZZ = 0.
 
      if (PSFSUM*ZZZ/NNNSUM.gt.0.03*(PIXSUM/NNNSUM+15).and.
     .   irmin.lt.16) then
         irmin = irmin + 1
         goto 5
         endif
 
 
 
c-------------------------------------------------
c
c ok, now we have the inner radius; where
c the star contribution is less than 10% of the
c sky brightness
c
      irmin = irmin
      irmax = irmin+3
  6   continue
      irnxt = irmax+1
      n1 = 0
      n2 = 0
      reff = 0.
      do i = -irmax-1, irmax+1
      do j = -irmax-1, irmax+1
         dx = ixc+i-xc
         dy = iyc+j-yc
         if (i**2+j**2.ge.(irmin-0.5)**2.and.
     .       i**2+j**2.lt.(irmax+0.5)**2.and.
     .       inimage(ixc+i,iyc+j)) then
                if (n1.lt.1000) n1 = n1 + 1
                sl1(n1) = pix(ixc+i,iyc+j)
     .                   - ZZZ*rpsf_photij(dx,dy,ixc,iyc)
                reff = reff + sqrt(dx**2+dy**2)
                endif
         enddo
         enddo
 
      reff = reff/n1
 
      call barsiggg(sl1,n1,bar1,sig1,n1u)
 
      mbar_skyopt = bar1
 
      return
      end
 
 
c------------------------------------------------
c
c
      subroutine barsiggg(xlist,NTOT,bar,sig,NUSE)
      implicit none
 
      integer NTOT
      real xlist(NTOT)
      real bar
      real sig
      integer NUSE
 
      integer n
      real*8    bsum, ssum
      integer nsum
      integer NIT
 
 
      bar = 0.e0
      sig = 9e9
      do NIT = 1, 20
         bsum = 0.
         ssum = 0.
         nsum = 0.
         do n = 1, NTOT
            if (abs(xlist(n)-bar).le.1.65*sig) then
               bsum = bsum + xlist(n)
               ssum = ssum + abs(xlist(n)-bar)
               nsum = nsum + 1
               endif
            enddo
         if (nsum.gt.0) bar = bsum / nsum
         if (nsum.gt.1) sig = ssum/(nsum-1)*1.1
         enddo
      NUSE = nsum
      if (nsum.le.1) sig = 0.999
 
      return
      end
 
 
c-----------------------------------------------
c
 
      real function z_xyoptXX(xin,yin,skin,pix,RADFIT)
      implicit none
      real*8 xin
      real*8 yin
      real*8 skin
      real pix(4096,4096)
      real    RADFIT
 
      real   rpsf_photijk
 
      real    dx, dy
      integer ix, iy
      integer i,  j
 
      integer NL, NLs
      real dl(500)
      real pl(500)
      real fl(500)
      real sl(500)
      real wl(500)
 
      real*8 sumpf
      real*8 sumff
      real*8 sumpfw
      real*8 sumffw
      real zuse
 
      integer NARG
      common /NARG_/NARG
      integer mmeth
      common /mmeth/ mmeth
 
 
      ix = int(xin+0.5)
      iy = int(yin+0.5)
 
      NLs = 0
      do i = ix-int(RADFIT+1),ix+int(RADFIT+1)
      do j = iy-int(RADFIT+1),iy+int(RADFIT+1)
         dx = i-xin
         dy = j-yin
         if (sqrt(dx**2+dy**2).lt.RADFIT) then
            NLs = NLs + 1
            if (NLs.gt.500) then
               print*,'z_opt: NLs.gt.500'
               stop
               endif
            pl(NLs) = pix(i,j)-skin
            fl(NLs) = rpsf_photijk(dx,dy,ix,iy)
            dl(NLs) = sqrt(dx**2+dy**2)
            endif
         enddo
         enddo
 
      sumpf = 0.
      sumff = 0.
      do NL = 01, NLs
         sumpf = sumpf + fl(NL)*pl(NL)
         sumff = sumff + fl(NL)*fl(NL)
         enddo
      zuse = sumpf/sumff
      if (zuse.lt.0) zuse = 0.
 
      do NL = 1, NLs
         sl(NL) = sqrt(zuse*fl(NL)+skin+15)
         wl(NL) = 1.00/sl(NL)**2
         enddo
 
      sumpfw = 0.
      sumffw = 0.
      do NL = 01, NLs
         sumpfw = sumpfw + fl(NL)*pl(NL)*wl(NL)
         sumffw = sumffw + fl(NL)*fl(NL)*wl(NL)
         enddo
      zuse = sumpfw/sumffw
      if (zuse.lt.0) zuse = 0.
 
      z_xyoptXX = zuse
 
      return
      end
 
 
      subroutine procsatn(pixarr,UPR_LIM)
      implicit none
 
      real pixarr(_PXDIMX_,_PXDIMY_)
      real UPR_LIM
 
      integer i, j
      real rpsf_photijk
 
      real dx, dy
      real ri
      real pixval,pixtot
      real psfval,psftot
      real zlist(9999)
      integer nl, nt
      integer ii, jj
      integer ir
      real reff
      integer iru
      real    fru
      integer nsat
 
c
c this is the radial profile of the WFC F675W PSF...
c
c   from r = 0 to r = 50 _flt pixels...
c
      real rprof(51)
      data  rprof /
     .   0.27000000,0.07940824,0.02189085,0.00431118,0.00219986,
     .   0.00184045,0.00074173,0.00039569,0.00032109,0.00023005,
     .   0.00017415,0.00016298,0.00013543,0.00010287,0.00008278,
     .   0.00007402,0.00005836,0.00004262,0.00003057,0.00002460,
     .   0.00002157,0.00001737,0.00001614,0.00001478,0.00001383,
     .   0.00001180,0.00001164,0.00001120,0.00000997,0.00000976,
     .   0.00000882,0.00000834,0.00000775,0.00000671,0.00000619,
     .   0.00000547,0.00000536,0.00000503,0.00000473,0.00000469,
     .   0.00000431,0.00000448,0.00000425,0.00000403,0.00000381,
     .   0.00000336,0.00000298,0.00000237,0.00000193,0.00000173,
     .   0.00000040 /
 
      write(*,'(''#  '')')
      write(*,'(''# ENTER PROCSATN: '')')
      write(*,'(''#        UPR_LIM: '',f8.2)') UPR_LIM
      write(*,'(''#  '')')
 
      nsat = 0
 
      do j = _PXDIMY_-2, 3, -1
      if (j.eq.j/128*128)
     .write(*,'(''#  '',i4.4,'' out of 4096.  NSAT: '',i6)') j,nsat
      do i = 3, _PXDIMX_-2
         if (pixarr(i,j).gt.UPR_LIM) then
            nsat = nsat + 1
            ir = 1
  5         continue
            if (ir.lt.10) ir = ir + 1
            nl = 0
            nt = 0
            psftot = 0.
            pixtot = 0.
            do ii = i-ir, i+ir
            do jj = j-ir, j+ir
               ri = sqrt(1.*(ii-i)**2+(jj-j)**2)
               if (ii.gt.1.and.ii.lt.(_PXDIMX_).and.
     .             jj.gt.1.and.jj.lt.(_PXDIMY_).and.
     .             ri.lt.ir) then
                   nt = nt + 1
                   if (pixarr(ii,jj).gt.UPR_LIM) goto 6
                   if (pixarr(ii+1,jj).gt.UPR_LIM.and.
     .                 pixarr(ii+2,jj).gt.UPR_LIM.and.
     .                 pixarr(ii+3,jj).gt.UPR_LIM) goto 6
                   if (pixarr(ii-1,jj).gt.UPR_LIM.and.
     .                 pixarr(ii-2,jj).gt.UPR_LIM.and.
     .                 pixarr(ii-3,jj).gt.UPR_LIM) goto 6
                   if (pixarr(ii,jj+1).gt.UPR_LIM.and.
     .                 pixarr(ii,jj+2).gt.UPR_LIM.and.
     .                 pixarr(ii,jj+3).gt.UPR_LIM) goto 6
                   if (pixarr(ii,jj-1).gt.UPR_LIM.and.
     .                 pixarr(ii,jj-2).gt.UPR_LIM.and.
     .                 pixarr(ii,jj-3).gt.UPR_LIM) goto 6
                   nl = nl + 1
                   pixval = pixarr(ii,jj)
                   reff = ri
                   if (reff.gt.3.5) reff = 3.5 + (reff-3.5)*0.5
                   iru = int(ri)
                   fru = ri-int(ri)
                   psfval = rprof(iru+1)+
     .               (fru)*(rprof(iru+2)-rprof(iru+1))
                   if (ri.lt.3) then
                      dx = ii-i
                      dy = jj-j
                      psfval = max(rpsf_photijk(dx,dy,ii,jj),
     .                             0.85*psfval)
                      endif
                   pixtot = pixtot + pixval
                   psftot = psftot + psfval
                   zlist(nl) = pixval/psfval
c                   if (i.eq.2060.and.j.eq.2095)
c     .                 write(*,339) ii,jj,ri,pixval,psfval,zlist(nl)
c 339               format(i4,1x,i4,1x,f8.2,1x,f9.1,1x,f8.6,1x,f10.1)
  6                continue
                   endif
                enddo
                enddo
            if (nl.lt.nt/2.and.ir.lt.10) goto 5
           !call rbubble(zlist,nl)
            pixarr(i,j) = pixtot/psftot*rprof(1)
            if (rpsf_photijk(0.,0.,ii,jj).gt.0) then
               pixarr(i,j) = pixtot/psftot*rpsf_photijk(0.,0.,ii,jj)
               endif
            if (psftot.le.0) pixarr(i,j) = UPR_LIM
            if (.not.(pixarr(i,j).gt.UPR_LIM)) pixarr(i,j) = UPR_LIM
            endif
c         if (i.eq.2060.and.j.eq.2095)  then
c              print*,'---> pixarr(i,j): ',pixarr(i,j)
c              stop
c              endif
         enddo!i
         enddo!j
 
 
      return
      end
 
 
      subroutine init_psf(psf)
      real psf(21,21)
      integer i, j
      real dx, dy
      real r
 
      real rprof(51)
      data  rprof /
     .   0.27000000,0.07940824,0.02189085,0.00431118,0.00219986,
     .   0.00184045,0.00074173,0.00039569,0.00032109,0.00023005,
     .   0.00017415,0.00016298,0.00013543,0.00010287,0.00008278,
     .   0.00007402,0.00005836,0.00004262,0.00003057,0.00002460,
     .   0.00002157,0.00001737,0.00001614,0.00001478,0.00001383,
     .   0.00001180,0.00001164,0.00001120,0.00000997,0.00000976,
     .   0.00000882,0.00000834,0.00000775,0.00000671,0.00000619,
     .   0.00000547,0.00000536,0.00000503,0.00000473,0.00000469,
     .   0.00000431,0.00000448,0.00000425,0.00000403,0.00000381,
     .   0.00000336,0.00000298,0.00000237,0.00000193,0.00000173,
     .   0.00000040 /
 
 
      do i = 01, 21
      do j = 01, 21
         dx = (i-11)/4.0
         dy = (j-11)/4.0
         r = sqrt(dx**2+dy**2)
         iru = int(r)
         fru = r-int(r)
         psf(i,j) = rprof(iru+1)+(fru)*(rprof(iru+2)-rprof(iru+1))
         enddo
         enddo
 
      return
      end
 
c-------------------------------------------------------
c
c this block of routines is designed to smooth a derived
c PSF perturbation, for use in the img2xym_WFC.09x10
c routine; the idea is that the PSF changes from exposure
c to exposure, but by and large it doesn't change much,
c and the change (in absolute PSF level) is about the same
c across the chip, so we can treat the PSF as:
c
c PSF_image(dx,dy,ix,iy) = PSF_library(dx,dy,ix,iy) + PERT_image(dx,dy)
c
c this routine ensures that the PERT part is smooth
c
c
      subroutine smoo_psfpert(mi,mo)
      implicit none
 
      real mi(101,101)
      real mo(101,101)
 
      real m0(101,101)
      real ma(101,101)
      real mb(101,101)
      real mc(101,101)
 
      integer i, j
      real    r
 
 
      do i = 1, 101
      do j = 1, 101
         m0(i,j) = mi(i,j)
         enddo
         enddo
 
      call sm5quad(mi,ma,101,101)
      call sm7quad(mi,mb,101,101)
      call sm7plan(mi,mc,101,101)
 
      do i = 1, 101
      do j = 1, 101
         r = sqrt((i-51.)**2+(j-51.)**2)/4.0
         mo(i,j) = mc(i,j)
         if (r.lt.07.5) mo(i,j) = (mb(i,j)+mc(i,j))/2
         if (r.lt.06.5) mo(i,j) =  mb(i,j)
         if (r.lt.05.5) mo(i,j) = (ma(i,j)+mb(i,j))/2
         if (r.lt.04.5) mo(i,j) =  ma(i,j)
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
      subroutine sm7quad(r,s,NX,NY)
      implicit none
 
      integer NX, NY
      real r(NX,NY)
      real s(NX,NY)
 
      integer i, j
      integer im, jm
 
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
 
 
 
 
c-----------------------------------------------------
c
c This routine will find the position and flux for a
c saturated star by identifying the closest unsaturated
c pixels and fitting a PSF to them, solving for flux
c and position.  The sky value (s) is input.
c
      subroutine find_xyz_SAT(x,y,z,s,qfit,pix)
      implicit none
 
      real*8 x,y,z,s
      real   qfit
      real   pix(_PXDIMX_,_PXDIMY_)
 
      integer i0, j0
      integer i, j
      real    plist(9999)
      real    flist(9999)
      real    qlist(9999)
      integer ilist(9999)
      integer jlist(9999)
      real    zlist(9999)
 
      real*8 x0, y0
      integer NT, NL, NLs
      real mbar_sky
      real rpsf_photijk
 
      real dx, dy, r, RMAX
      real QOFMIN, QOFDX, QOFDY, QOFFLX, QOF, dQOF
      real*8 PTOT, FTOT, FLX
      integer NNN
      common /nqq/nnn
      data NNN/0/
 
      integer HIFLAG
      common / HIFLAG_ / HIFLAG
 
      integer LOFLAG
      common / LOFLAG_ / LOFLAG
 
      real*8 RTOT, REFF, PEFF
      real bar, sig
      integer NIT
      real    DD
      real    ddx, ddy
      integer idx, idy
 
      i0 = int(x+0.5)
      j0 = int(y+0.5)
 
      s = mbar_sky(i0,j0,25,35,pix)
 
c      print*,'FIND_XYZ_SAT: ',i0,j0,s,pix(i0,j0)
 
      RMAX = 2.5
  1   continue
      NT = 0
      NL = 0
      RMAX = RMAX + 0.25
      PTOT = 0.
      RTOT = 0.
      do i = i0-int(RMAX+0.5),i0+int(RMAX+0.5)
      do j = j0-int(RMAX+0.5),j0+int(RMAX+0.5)
         r = sqrt(((i-i0)**2+(j-j0)**2)*1.0)
         if (r.ge.RMAX) goto 2
         if (i.lt.2.or.i.gt.(_PXDIMX_-1)) goto 2
         if (j.lt.2.or.j.gt.(_PXDIMY_-1)) goto 2
         NT = NT + 1
         if (pix(i,j  ).gt.HIFLAG) goto 2
         if (pix(i,j+1).gt.HIFLAG) goto 2
         if (pix(i,j-1).gt.HIFLAG) goto 2
         if (pix(i,j  ).lt.LOFLAG) goto 2
         if (pix(i,j+1).lt.LOFLAG) goto 2
         if (pix(i,j-1).lt.LOFLAG) goto 2
         if (pix(i+1,j).gt.HIFLAG.and.pix(i+2,j).gt.HIFLAG) goto 2
         if (pix(i-1,j).gt.HIFLAG.and.pix(i-2,j).gt.HIFLAG) goto 2
         NL = NL + 1
         plist(NL) = pix(i,j)-s
         ilist(NL) = i
         jlist(NL) = j
         PTOT = PTOT + plist(NL)
         RTOT = RTOT + r
 2       continue
         enddo
         enddo
c      print*,'---> ',RMAX,NT,NL
      if (NL.lt.0.75*NT.and.RMAX.lt.9) goto 1
 
      NLs = NL
      REFF = RTOT/NLs
 
      x0 = i0
      y0 = j0
      do NIT = 1, 4
         if (NIT.eq.2.and.QOFFLX.gt.10**(15.0/2.5)) then
            s = mbar_sky(i0,j0,25,35,pix)
            endif
         if (NIT.eq.2.and.QOFFLX.gt.10**(17.5/2.5)) then
            s = mbar_sky(i0,j0,30,40,pix)
            endif
         if (NIT.eq.1) DD = 0.30
         if (NIT.eq.2) DD = 0.10
         if (NIT.eq.3) DD = 0.03
         if (NIT.eq.4) DD = 0.01
         QOFMIN = 9.9e9
         QOFDX  = 0.0
         QOFDY  = 0.0
         QOFFLX = 0.0
         do idx = -4,4,1
         do idy = -4,4,1
            dx = idx*dd
            dy = idx*dd
            FTOT = 0.
            do NL = 1, NLs
               ddx = ilist(NL)-x0-dx
               ddy = jlist(NL)-y0-dy
               flist(NL) = rpsf_photijk(ddx,ddy,i0,j0)
               FTOT = FTOT + plist(NL)/flist(NL)
               zlist(NL) = plist(NL)/flist(NL)
               enddo
            FLX = FTOT/NL
            call barsig_FSAT(zlist,NLs,bar,sig)
            FLX = bar
            QOF = 0.
            do NL = 1, NLs
               dQOF = abs(plist(NL)-FLX*flist(NL))
               if (dQOF.gt.plist(NL)*0.20) dQOF = plist(NL)*0.20
               QOF = QOF + dQOF
              enddo
            if (QOF.lt.QOFMIN) then
               QOFMIN = QOF
               QOFDX  = DX
               QOFDY  = DY
               QOFFLX = FLX
               PEFF = 0.
               do NL = 1, NLs
                  qlist(NL) = flist(NL)
                  PEFF = PEFF + flist(NL)
                  enddo
               PEFF = PEFF/NLs
               endif
            enddo
            enddo
         x0 = x0 + QOFDX
         y0 = y0 + QOFDY
         enddo!NIT
 
      x = x0
      y = y0
      z = QOFFLX
      qfit = QOFMIN/PTOT
 
      if (z.lt.1) z = 1
 
      write(44,144) x,y,-2.5*log10(z),s,z,REFF,PEFF,qfit
 144  format(f9.3,1x,f9.3,1x,f9.4,1x,f8.2,1x,f14.1,1x,
     .       f8.3,1x,f10.7,1x,f8.4)
 
      return
      end
 
c------------------------------------------------
c
      subroutine barsig_FSAT(xlist,NTOT,bar,sig)
      implicit none
 
      integer NTOT
      real xlist(NTOT)
      real bar
      real sig
 
      integer n
      real*8  bsum, ssum
      integer NIT
      integer NSUM
      real sig0
 
      bar = 0.e0
      sig = 9e9
      NIT = 0
  1   continue
      NIT = NIT + 1
      sig0 = sig
      bsum = 0.
      ssum = 0.
      nsum = 0
      do n = 1, NTOT
         if (abs(xlist(n)-bar).le.3.50*sig) then
            bsum = bsum + xlist(n)
            ssum = ssum + abs(xlist(n)-bar)
            nsum = nsum + 1
            endif
         enddo
      if (nsum.gt.0) bar = bsum/nsum
      if (nsum.gt.1) sig = ssum/nsum*1.10
      if (NIT.lt.20) goto 1
      if (nsum.le.1) sig = 0.999
 
      return
      end
 
 
 
 
 
      subroutine zero_fix(pix)
      implicit none
      real    pix(4096,4096)
 
      integer hist1(25), hist2(25), hist3(25)
      integer i, j
      integer PM2, PM1, P00, PP1, PP2
      integer EM2, EM1, E00, EP1, EP2
      integer ERR
      real    f
      real    A, B, C
      real    fmin, emin
      real    rand
      integer fi
 
      print*,'          '
      print*,'          '
      print*,'ZERO_FIX: (for antiquated way of storing images)'
      print*,'          '
      print*,'          '
 
      call find_hist(pix,hist1)
 
      do i = 0001, 4096
      do j = 0001, 4096
        !write(*,'(2i5,1x,f10.3)') i,j,pix(i,j)
         if (pix(i,j).lt.-0.9) pix(i,j) = pix(i,j)-1
         enddo
         enddo
 
      call find_hist(pix,hist2)
 
c      print*,' '
c      do h = 01, 25
c         hz = h-13
c         write(*,'(i3,1x,i3,1x,2i10)') h,hz,hist1(h),hist2(h)
c         enddo
c      print*,' '
 
      fmin = 0
      emin = 9e9
c      print*,' '
      do fi = 0, 100, 1
         f = fi*0.01
         PM2 = hist2(11)
         PM1 = hist2(13)*( f )
         P00 = hist2(13)*(1-f)
         PP1 = hist2(14)
         PP2 = hist2(15)
         A = P00
         B = ((PP1-PM1)/2.0 + (PP2-PM2)/4.0)/2.0
         C = (PP2+PM2-2*A)/8.0
 
c-----------------------------------------------------
c         A = -0.08571*PM2
c     .       +0.34286*PM1
c     .       +0.48572*P00
c     .       +0.34286*PP1
c     .       -0.08571*PP2
c         B = -0.20000*PM2
c     .       -0.10000*PM1
c     .       +0.00000*P00
c     .       +0.10000*PP1
c     .       +0.20000*PP2
c         C = +0.71430*PM2
c     .       -0.35710*PM1
c     .       -0.71430*P00
c     .       -0.35710*PP1
c     .       +0.71430*PP2
c         print*,' '
c         write(*,'(30x,3f15.2)') A,B,C
c         A = P00
c         B = ((PP1-PM1)/2.0 + (PP2-PM2)/4.0)/2.0
c         C = (PP2+PM2-2*A)/8.0
c         write(*,'(30x,3f15.2)') A,B,C
c         print*,' '
c-----------------------------------------------------
 
         EM2 = A - 2*B + 4*C
         EM1 = A -   B +   C
         E00 = A
         EP1 = A +   B +   C
         EP2 = A + 2*B + 4*C
         ERR = abs(PM2-EM2)
     .       + abs(PM1-EM1)
     .       + abs(P00-E00)
     .       + abs(PP1-EP1)
     .       + abs(PP2-EP2)
c         write(*,'(f8.4,10x,5i8,10x,5i8,10x,i12)') f,
c     .             PM2,PM1,P00,PP1,PP2,
c     .             EM2,EM1,E00,EP1,EP2,
c     .             ERR
         if (ERR.lt.emin) then
            emin = ERR
            fmin = f
            endif
         enddo
c      print*,' '
c      print*,'fmin: ',fmin
c      print*,' '
 
      do i = 0001, 4096
      do j = 0001, 4096
         if (int(pix(i,j)).eq.0) then
            if (rand().lt.fmin) pix(i,j) = -1
            endif
         enddo
         enddo
 
      call find_hist(pix,hist3)
 
c     print*,' '
c      do h = 01, 25
c        hz = h-13
c        write(*,'(i3,1x,i3,1x,3i10)') h,hz,hist1(h),hist2(h),hist3(h)
c        enddo
c      print*,' '
 
      return
 
      end
 
 
c-----------------------------------
c
c
      subroutine find_hist(pix,hist)
      implicit none
 
      real pix(4096,4096)
      integer hist(25)
 
      integer h
      integer i, j
 
      do h = 01, 25
         hist(h) = 0
         enddo
 
      do i = 0001, 4096
      do j = 0001, 4096
         h = 13 + pix(i,j)
         if (h.ge.01.and.h.le.25) hist(h) = hist(h)+1
         enddo
         enddo
 
      return
      end
 
 
 
c------------------------------------------------
c
c this routine will use the header info to figure
c out which filter number the observations were
c taken with, so that the correct distortion
c correction can be used...
c
c  0) use no table correction, only the polynomial
c  1) F435W
c  2) F475W
c  3) F555W
c  4) F606W
c  5) F814W
c  6) F850L
c  7) F625W
c  8) F658N
c  9) F775W
c 10) F660N
c 11) F502N
c 12) F550M
c
c
      subroutine find_ifilt(IFILT,HDR)
      implicit none
 
      integer IFILT
      character*70 HDR(25)
 
      character*5 FILTNAME(12)
      data FILTNAME/'F435W','F475W','F555W','F606W','F814W',
     .              'F850L','F625W','F658N','F775W','F660N',
     .              'F502N','F550M'/
      integer i
 
      IFILT = 0
c      do h = 1, 25
c        write(*,'(i3,1x,a20,1x,a5)') h,HDR(H),HDR(H)(2:6)
c         enddo
 
      do i = 1, 12
         if (HDR(22)(1:5).eq.FILTNAME(i)) IFILT = i
         if (HDR(22)(2:6).eq.FILTNAME(i)) IFILT = i
         if (HDR(23)(1:5).eq.FILTNAME(i)) IFILT = i
         if (HDR(23)(2:6).eq.FILTNAME(i)) IFILT = i
         enddo
 
      print*,'--->              '
      print*,'---> find_ifilt:  '
      if (IFILT.ne.0) print*,'--->      IFILT: ',IFILT,FILTNAME(IFILT)
      if (IFILT.eq.0) print*,'--->      IFILT: ',IFILT,'NOFILT'
      print*,'--->              '
      return
      end
 
 
 
      subroutine rbarsigs(xlist,NTOT,bar,sig,NUSE,SIGCLIP)
      implicit none
 
      integer NTOT
      real xlist(NTOT)
      real bar
      real sig
      integer NUSE
      real SIGCLIP
 
      integer n
      real*8    bsum, ssum
      integer nsum, nsumo
      integer NIT
 
 
      nsum = 0
      bar = 0.e0
      sig = 9e9
      do NIT = 1, 20
         nsumo = nsum
         bsum = 0.
         ssum = 0.
         nsum = 0.
         do n = 1, NTOT
            if (abs(xlist(n)-bar).le.SIGCLIP*sig) then
               bsum = bsum + xlist(n)
               ssum = ssum + abs(xlist(n)-bar)
               nsum = nsum + 1
               endif
            enddo
         if (nsum.gt.0) bar = bsum / nsum
         if (nsum.gt.1) sig = ssum/(nsum-1)
         if (nsum.lt.0.35*NTOT.and.NIT.ge.3) return
         if (nsum.eq.nsumo.and.NIT.ge.3) goto 1
         enddo
  1   continue
      NUSE = nsum
      if (nsum.le.1) sig = 0.999
 
      return
      end
 
 
c------------------------------------------------------------
c
c this routine will derive a perturbation for a library
c PSF so that it will better fit an image
c
      subroutine find_psfpert(pixc,psfnim,hmb,psfbar)
      implicit none
 
      real pixc(_PXDIMX_,_PXDIMY_)
      real psfnim(101,101)
      real hmb
      real psfbar(101,101)
 
      integer i, ii, iii
      integer j, jj, jjj
 
      real*8 xr, yr, fr, sr
      real   qfit
 
      real   dx, dy
      real   rpsf_photijk
 
      real    FEST
      integer hobs, fnd_hloc
      integer pobs, fnd_ploc
      real    mbar_sky
      integer NL, NLs, NLu
      integer il(_NLMAX_)
      integer jl(_NLMAX_)
      real*8  xl(_NLMAX_)
      real*8  yl(_NLMAX_)
      real    ml(_NLMAX_)
      real    zl(_NLMAX_)
      real    sl(_NLMAX_)
      real    ql(_NLMAX_), qbar, qsig
      integer ul(_NLMAX_)
      real    pixl(9,9,_NLMAX_)
      real    psfl(9,9,_NLMAX_)
      real    psfq
 
      integer NSs, NSu
      real    dpsf(_NLMAX_)
 
      integer ipsf, jpsf
      real    xpsf, ypsf, rpsf
      real    bar(101,101)
      real    sig(101,101)
      real    new(101,101)
      real    smu(101,101)
      real    image(401,401)
      integer NIT
      real    fx, fy
      integer ir, nr
      real    presum,psfsum,prtsum
 
      integer NIMu
      common /NIMu_/NIMu
      integer NPERT
      common /NPERT_/NPERT
      data NPERT/0/
      character*80 FILENAME
      real sum
      real*8 psfnimtot
      real*8 psfbartot
      real*8 spsfpsf
      real*8 snimpsf
      real   FMIN
 
      NPERT = NPERT + 1
 
      print*,'        '
      print*,'        '
      print*,'ENTER FIND_PSFPERT: '
      print*,'        '
      print*,'   FIND THE AVERAGE DIFFERENCE BETWEEN THE      '
      print*,'   INPUT LIBRARY PSF AND THE PSF IN THIS IMAGE  '
      print*,'        '
 
 
      do i = 001, 101
      do j = 001, 101
         psfnim(i,j) = 0.
         enddo
         enddo
 
 
      open(85,file='LOG.PSFPERT.XYM',status='unknown')
      FMIN = 25000
  3   continue
      print*,'             '
      print*,' FIND STARS TO USE: '
      print*,'          FMIN: ',FMIN
      print*,'           HMB: ',HMB
      print*,'             '
      NLs = 0
      do jj = 30, _PXDIMX_-30
      do ii = 30, _PXDIMX_-30
         if (pixc(ii,jj).lt.FMIN/4) goto 222
         if (pixc(ii,jj).gt.54999) goto 222
 
         FEST = pixc(ii,jj) +
     .          max(pixc(ii+1,jj)+pixc(ii+1,jj+1)+pixc(ii,jj+1),
     .              pixc(ii+1,jj)+pixc(ii+1,jj-1)+pixc(ii,jj-1),
     .              pixc(ii-1,jj)+pixc(ii-1,jj+1)+pixc(ii,jj+1),
     .              pixc(ii-1,jj)+pixc(ii-1,jj-1)+pixc(ii,jj-1))
     .        - hmb*4
 
         if (FEST.gt.75000) goto 222
         if (FEST.lt. FMIN) goto 222
         hobs = fnd_hloc(ii,jj,pixc)
         pobs = fnd_ploc(ii,jj,pixc,hmb)
         if (hobs.lt.9) goto 222
         if (pobs.lt.7) goto 222
         fx = (pixc(ii+1,jj)-pixc(ii-1,jj))/2/
     .        (pixc(ii  ,jj)-min(pixc(ii+1,jj),pixc(ii-1,jj)))
         fy = (pixc(ii,jj+1)-pixc(ii,jj-1))/2/
     .        (pixc(ii,jj  )-min(pixc(ii,jj+1),pixc(ii,jj-1)))
         xr = ii + fx
         yr = jj + fy
         sr = mbar_sky(ii,jj,15,25,pixc)
c
c using the following routine, the fit for the star's center will
c use only the central pixels, so that the center of the PSFPERT
c will be flat (it will only account for the height of the core, but
c won't try to impose an offset).  The flux is set using the inner
c r=3.5 pixels
c
         call find_xyzXX_CENFOCUS(xr,yr,fr,sr,qfit,pixc,1.0,1.5)
c        write(73,173) xr,yr,fr,sr,qfit,pixc(ii,jj)
c        write( *,173) xr,yr,fr,sr,qfit,pixc(ii,jj)
  173    format(10x,f8.2,1x,f8.2,1x,f9.1,1x,f8.3,1x,f8.4,1x,f9.2)
         if (qfit.gt.0.5) goto 222  ! this may help prevent CRs from
         NLs = NLs + 1              ! being included
         if (NLs.gt._NLMAX_) then
            print*,'  NLs    : ', NLs
            print*,'       ii: ', ii
            print*,'       jj: ', jj
            print*,'  _NLMAX_: ',_NLMAX_
            stop
            endif
         il(NLs) = ii
         jl(NLs) = jj
         xl(NLs) = xr
         yl(NLs) = yr
         ml(NLs) = -2.5*log10(max(fr,1.))
         zl(NLs) = fr
         sl(NLs) = sr
         ql(NLs) = qfit
         ul(NLs) = 1
         do i = 1, 9
         do j = 1, 9
            iii = ii+(i-5)
            jjj = jj+(j-5)
            dx  = iii-xl(NLs)
            dy  = jjj-yl(NLs)
            pixl(i,j,NLs) = pixc(iii,jjj) - sl(NLs)
            psfl(i,j,NLs) = rpsf_photijk(dx,dy,ii,jj)
            enddo
            enddo
         if ((NLs.lt.00010).or.
     .       (NLs.lt.00100.and.NLs.eq.NLs/0010*0010).or.
     .       (NLs.lt.01000.and.NLs.eq.NLs/0100*0100).or.
     .       (NLs.lt.10000.and.NLs.eq.NLs/1000*1000))
     .   write( *,133) xl(NLs),yl(NLs),ml(NLs),zl(NLs),sl(NLs),
     .                 NLs,il(NLs),jl(NLs),
     .                 int(pixc(il(NLs),jl(NLs))),
     .                 psfl(5,5,NLs),ql(NLs)
c         write(85,133) xl(NLs),yl(NLs),ml(NLs),zl(NLs),sl(NLs),
c     .                 NLs,il(NLs),jl(NLs),
c     .                 int(pixc(il(NLs),jl(NLs))),
c     .                 psfl(5,5,NLs),ql(NLs),
c     .                 xr0,yr0,-2.5*log10(max(fr0,1.0))
 133     format(f9.3,1x,f9.3,1x,f8.4,1x,f10.1,1x,f7.2,3x,
     .                 i6,3x,i4,1x,i4,1x,i5,1x,f8.6,1x,f8.6,1x,
     .                 f9.3,1x,f9.3,1x,f8.4)
 222     continue
         enddo
         enddo
c      close(85)
 
 
      call rbarsigs(ql,NLs,qbar,qsig,NLu,3.)
 
 
c
c only use the ones that are similar (not CRs or whatever)
c
      do NL = 1, NLs
         ul(NL) = 1
         if (ql(NL).gt.1.5*qbar+2*qsig) ul(NL) = 0
         enddo
 
 
      print*,'     '
      print*,'    NLs: ',NLs
      print*,'   qbar: ',qbar
      print*,'   qsig: ',qsig
      print*,'   qlim: ',1.5*qbar+2*qsig
      print*,'    NLu: ',NLu
      print*,'   NIMu: ',NIMu
      print*,'     '
 
      if (NLs.lt.035.and.FMIN.gt.2000) then
         print*,'     '
         print*,'NOT ENOUGH STARS FOR PERT: ',NLs
         print*,'     FOR A FLUX CUTOFF OF: ',FMIN
         print*,'     '
         FMIN = FMIN/2
         print*,'WILL TRY AGAIN WITH FMIN: ',FMIN
         print*,'     '
         print*,'     '
         goto 3
         endif
 
 
      print*,'                        '
      print*,'                        '
      print*,'Now, iterate to find the best psf-adjustment...'
      print*,'                        '
      do NIT = 1,4
      print*,'                        '
      print*,'************************'
      print*,'*****                   '
      print*,'***** FIND PSF_PERT     '
      print*,'*****                   '
      print*,'*****    NIT: ',NIT,rpsf_photijk(0.,0.,2048,2048)
      print*,'*****                   '
      print*,'************************'
      print*,'                        '
 
      print*,'                 '
      print*,'  psfnim(51,51): ',psfnim(51,51)
      print*,'                 '
 
      do NL = 1, NLs
         do i = 1, 9
         do j = 1, 9
            iii = il(NL)+(i-5)
            jjj = jl(NL)+(j-5)
            dx  = iii-xl(NL)
            dy  = jjj-yl(NL)
            NIMu = 0
            psfq         = rpsf_photijk(dx,dy,il(NL),jl(NL))
            NIMu = 1
            psfl(i,j,NL) = rpsf_photijk(dx,dy,il(NL),jl(NL))
            enddo
            enddo
         enddo
 
      do jpsf = 001, 101
      do ipsf = 001, 101
         bar(ipsf,jpsf) = 0.
         sig(ipsf,jpsf) = 0.
         xpsf = (ipsf-51)/4.0
         ypsf = (jpsf-51)/4.0
         rpsf = sqrt(xpsf**2+ypsf**2)
         if (rpsf.lt.5.0) then
            NSs = 0
            do NL = 1, NLs
               xr = xl(NL) + xpsf
               yr = yl(NL) + ypsf
               ii = int(xr+0.5)
               jj = int(yr+0.5)
               i  = 5 + ii - il(NL)
               j  = 5 + jj - jl(NL)
               dx = ii-xl(NL)
               dy = jj-yl(NL)
               if (i.ge.1.and.i.le.9.and.
     .             j.ge.1.and.j.le.9.and.
c    .             abs(dx-xpsf).le.0.30.and.
c    .             abs(dy-ypsf).le.0.30.and.
     .             ul(NL).eq.1) then
                   NSs = NSs + 1
                   dpsf(NSs) = pixl(i,j,NL)/zl(NL) - psfl(i,j,NL)
                   endif
               enddo
            call rbarsigs(dpsf,NSs,bar(ipsf,jpsf),sig(ipsf,jpsf),NSu,2.)
            if (ipsf.eq.51.and.(jpsf-51).eq.(jpsf-51)/2*2)
     .      write( *,174) ipsf,jpsf,rpsf,
     .                    bar(ipsf,jpsf),sig(ipsf,jpsf),NSs,
     .                    psfnim(ipsf,jpsf)
 174        format(i2,1x,i2,1x,f8.3,1x,f8.5,1x,f8.5,1x,i5,1x,f8.5)
            endif
         enddo
         enddo
      print*,'                 '
      print*,' CENTRAL VALUE:   psfnim(51,51): ',psfnim(51,51)
      print*,'                 '
 
      sum = 0.
      do ipsf = 001, 101
      do jpsf = 001, 101
         xpsf = (ipsf-51)/4.0
         ypsf = (jpsf-51)/4.0
         rpsf = sqrt(xpsf**2+ypsf**2)
         new(ipsf,jpsf) = psfnim(ipsf,jpsf) + 0.90*bar(ipsf,jpsf)
         sum = sum + new(ipsf,jpsf)
         enddo
         enddo
      print*,'                 '
      print*,'  psfnim(51,51): ',psfnim(51,51)
      print*,'            sum: ',sum/16.0
      print*,'                 '
 
      call smoo_psfpert(new,smu)
 
      do ipsf = 001, 101
      do jpsf = 001, 101
         image(ipsf+000,jpsf+100*(NIT-1)) = psfnim(ipsf,jpsf)
         image(ipsf+100,jpsf+100*(NIT-1)) =    bar(ipsf,jpsf)
         image(ipsf+200,jpsf+100*(NIT-1)) =    smu(ipsf,jpsf)-
     .                                      psfnim(ipsf,jpsf)
         image(ipsf+300,jpsf+100*(NIT-1)) =    smu(ipsf,jpsf)
         enddo
         enddo
      FILENAME = 'LOG.pert.fits'
      call writfits_r4(FILENAME,image,401,401)
 
      psfnimtot = 0.
      psfbartot = 0.
      spsfpsf = 0.
      snimpsf = 0.
      do ipsf = 01, 101
      do jpsf = 01, 101
         psfnim(ipsf,jpsf) = smu(ipsf,jpsf)
         psfnimtot = psfnimtot + psfnim(ipsf,jpsf)
         psfbartot = psfbartot + psfbar(ipsf,jpsf)
         spsfpsf = spsfpsf + psfbar(ipsf,jpsf)*psfbar(ipsf,jpsf)
         snimpsf = snimpsf + psfnim(ipsf,jpsf)*psfbar(ipsf,jpsf)
         enddo
         enddo
      print*,'                 '
      print*,'  psfnim(51,51): ',psfnim(51,51)
      print*,'     new(51,51): ',new(51,51)
      print*,'     smu(51,51): ',smu(51,51)
      print*,'     bar(51,51): ',bar(51,51)
      print*,'                 '
      print*,'----> psfnimtot: ',psfnimtot/16
      print*,'----> psfbartot: ',psfbartot/16
      print*,'---->       rat: ',snimpsf/spsfpsf
      print*,'                 '
      psfnimtot = 0.
      do ipsf = 01, 101
      do jpsf = 01, 101
         psfnimtot = psfnimtot + psfnim(ipsf,jpsf)
         enddo
         enddo
      print*,'----> psfnimtot: ',psfnimtot/16
      print*,'  psfnim(51,51): ',psfnim(51,51)
      print*,'                 '
      print*,' AVG DESIRED ADJUSTMENT x10000  '
      print*,'                 '
      do jpsf = 61, 41, -1
         write(*,'(i2,1x,21i6)') jpsf,(int(10000*bar(ipsf,jpsf)+0.5),
     .      ipsf=41,61,1)
         enddo
      print*,'                 '
      print*,'                 '
      print*,' AVG TOTAL PSF x10000  '
      print*,'                 '
      do jpsf = 61, 41, -1
         write(*,'(i2,1x,21i6)') jpsf,(int(10000*psfbar(ipsf,jpsf)+0.5),
     .      ipsf=41,61,1)
         enddo
      print*,'                 '
      print*,'                 '
      print*,' NEW TOTAL ADJUSTMENT x10000  '
      print*,'                 '
      write(*,'(i2,1x,21i6)') jpsf,(ipsf,ipsf=41,61,1)
      print*,'                 '
      do jpsf = 61, 41, -1
         write(*,'(i2,1x,21i6)') jpsf,(int(10000*psfnim(ipsf,jpsf)+0.5),
     .      ipsf=41,61,1)
         enddo
 
 
      enddo!NIT
 
      close(67)
 
c     print*,'         '
c     print*,'     4    bar    (new-smu)   (smu-old)   smu  '
c     print*,'     3    bar    (new-smu)   (smu-old)   smu  '
c     print*,'     2    bar    (new-smu)   (smu-old)   smu  '
c     print*,'     1    bar    (new-smu)   (smu-old)   smu  '
 
 
c
c show the radial behavior
c
      print*,' '
      print*,'PERTPSF RADIAL PROFILE: '
      print*,' '
  136 format('.. .....    ..........    ..........',
     .       '    ..........    ..........')
  137 format('IR  NPIX     UNPERTPSF      NEWPSF  ',
     .       '     JUSTPERT      TOTPERT  ')
      write(*,136)
      write(*,137)
      write(*,136)
      do ir = 0, 12, 1
         nr = 0
         presum = 0
         psfsum = 0
         prtsum = 0
         do ipsf = 001, 101
         do jpsf = 001, 101
            xpsf = (ipsf-51)/4.0
            ypsf = (jpsf-51)/4.0
            rpsf = sqrt(xpsf**2+ypsf**2)
            if (ir-rpsf.gt.-0.5.and.ir-rpsf.le.0.5) then
               nr = nr + 1
               NIMu = 0
               presum = presum + rpsf_photijk(xpsf,ypsf,2048,2048)
               NIMu = 1
               psfsum = psfsum + rpsf_photijk(xpsf,ypsf,2048,2048)
               prtsum = prtsum + psfnim(ipsf,jpsf)
               endif
            enddo
            enddo
         write(*,'(i2,1x,i5,4x,f10.7,4x,f10.7,4x,f10.7,4x,f10.7)')
     .       ir,nr,presum/nr,psfsum/nr,prtsum/nr,prtsum/16
         enddo
      print*,' '
 
      open(83,file='LOG.PSFPERT.OUT',status='unknown')
      do ipsf = 001, 101
      do jpsf = 001, 101
         xpsf = (ipsf-51)/4.0
         ypsf = (jpsf-51)/4.0
         rpsf = sqrt(xpsf**2+ypsf**2)
         if (rpsf.lt.8) write(83,183) ipsf,jpsf,psfnim(ipsf,jpsf)
         enddo
         enddo
 183  format(i3,1x,i3,1x,f9.7)
 
      do ipsf = 001, 101
         psfnim(ipsf,001) = 0.
         psfnim(ipsf,101) = 0.
         psfnim(001,ipsf) = 0.
         psfnim(101,ipsf) = 0.
         enddo
 
 
      return
      end
 
 
 
 
 
c----------------------------------------------------
c
c This routine will generate
c
c
c
c
c apl(01):   1-pixel  r< 1 (central pixel)
c apl(02):   9-pixel, r< 2 (inner 3x3 pixels)
c apl(03):  25-pixel, r< 3 (inner 5x5 pixels)
c apl(04):  45-pixel, r< 4 (inner 7x7 pixels, minus corners)
c apl(05):  69-pixel, r< 5
c apl(06): 109-pixel, r< 6
c apl(07): 145-pixel, r< 7
c apl(08): 193-pixel, r< 8
c apl(09): 249-pixel, r< 9
c apl(10): 305-pixel, r<10
c
c The output of this routine  will detail how the measured flux for a star may
c vary depending on the aperture used to measure it.
c
c Here, we measure stars with apertures with outer radii between 1 pixel and
c 10 pixels.  The 1-pixel aperture has just the central pixel for the star.
c The 10-pixel aperture includes the 305 pixels within r=10.
c
c The aperure includes all pixels that are reasonably consistent with the
c stellar profile.  Pixels that disagree slightly from the PSF are included
c so that the result will *not* be very sensitive to the PSF model.
c
c The model is used to normalize all these aperture-magnitudes to the same
c system.  Essentially the PSF provides an "aperture correction", but does not
c contribute to any harsher pixel-weighting.  You add up the sum of the flux
c over sky, the sum of the PSF over the aperture (I always only use whole pixels!)
c and divide the two to get the properly scaled flux for the star, which
c should correspond to an aperture of 10 pixels.
c
c This can be useful, becuase it will tell us how the PSF-fit flux (the 3rd column)
c is or is not consistent with more rigid, aperture-based fluxes.  For instance,
c you can plot mfit-vs-map_05 to see if there is any flux out beyond the PSF-fit
c radius that the PSF does not account for.  The problem with using larger apertures
c is that there is a larger chance that it will be corrupted by either a CR or a
c neighbor.
c
c In such an examination, it is good to look at the lower bound of the mag residual.
c
c
c
 
 
      subroutine apstudy(ic,jc,x,y,z,s0,pixc,apl)
      implicit none
 
      integer ic,jc
      real*8  x,y,z,s0
      real pixc(_PXDIMX_,_PXDIMY_)
      real apl(10)
 
      real*8  ss
      real    mbar_sky
      real    pixarr(21,21)
      real    psfarr(21,21)
      real    sigarr(21,21)
      integer usearr(21,21)
      integer gudarr(21,21)
      integer nrej1, nrej2
 
      real rpsf_photijk
      integer i,ii
      integer j,jj
      real*8 sumpix
      real*8 sumpsf
      integer np, nu
      integer ipix, jpix
 
      integer rap
      real dx, dy
      real zr
      real sig1, sig2
      logical APFIRST
      common /APFIRST_/APFIRST
 
      if (APFIRST) then
         APFIRST = .false.
         write(33,'(''#                                           '')')
         write(33,'(''# SEARCH ON APNOTE in the code...           '')')
         write(33,'(''#                                           '')')
         write(33,135)
         write(33,134)
         endif
 
      if (ic.lt.15.or.ic.gt.4080.or.
     .    jc.lt.15.or.jc.gt.4080) then
          return
          endif
 
      do rap = 1, 10
         apl(rap) = 0
         enddo
 
      ss = mbar_sky(ic,jc,15,25,pixc)
 
      do i = -10, 10
      do j = -10, 10
         dx = ic+i-x
         dy = jc+j-y
         ipix = ic+i
         jpix = jc+j
         pixarr(11+i,11+j) = pixc(ipix,jpix)
         psfarr(11+i,11+j) = rpsf_photijk(dx,dy,ic,jc)
         sig1 = sqrt(max(z*psfarr(11+i,11+j)+ss,0.0)+14.0)
         sig2 = 0.05*z*psfarr(11+i,11+j)
         sigarr(11+i,11+j) = sqrt(sig1**2+sig2**2)
         enddo
         enddo
 
c
c find the bad pixels...
c
      nrej1 = 0
      do i = -10, 10
      do j = -10, 10
         gudarr(11+i,11+j) = 1
         if (i**2+j**2.ge.3**2) then ! don't reject anything in inner 5x5
            if (abs(pixarr(11+i,11+j)-ss-z*psfarr(11+i,11+j)).gt.
     .            7.5*sigarr(11+i,11+j)) then
               gudarr(11+i,11+j) = 0
               endif
            endif
         nrej1 = nrej1 + (gudarr(11+i,11+j)-1)
         enddo
         enddo
 
c
c reject things within 1 pixel of a bad pixel
c
      nrej2 = 0
      do i = -10, 10
      do j = -10, 10
         usearr(11+i,11+j) = 1
         if (i**2+j**2.ge.3**2) then ! don't reject anything in inner 5x5
            do ii = max(-10,i-1), min(10,i+1)
            do jj = max(-10,j-1), min(10,j+1)
               if (gudarr(11+ii,11+jj).eq.0) usearr(11+i,11+j) = 0
               enddo
               enddo
            endif
         nrej2 = nrej2 + (usearr(11+i,11+j)-1)
         enddo
         enddo
 
      do rap = 1, 10
         sumpix = 0.
         sumpsf = 0.
         np = 0
         nu = 0
         do i = -10, 10
         do j = -10, 10
            if (i**2+j**2.lt.rap**2) then
               np = np + 1
               nu = nu + usearr(11+i,11+j)
               sumpix =sumpix + (pixarr(11+i,11+j)-ss)*usearr(11+i,11+j)
               sumpsf =sumpsf +  psfarr(11+i,11+j)   *usearr(11+i,11+j)
               endif
            enddo
            enddo
         zr = sumpix/sumpsf
         apl(rap) = -2.5*log10(max(zr,1.0))
         enddo
 
      write(33,133) x,y,-2.5*log10(max(z,1.0)),ss,(apl(rap),rap=1,10)
  133 format(f8.2,1x,f8.2,1x,f7.3,1x,f8.2,4x,10(f7.3,1x))
  134 format('#.......',1x,'........',1x,'.......',1x,'........',4x,
     .       '.......',1x,'.......',1x,'.......',1x,'.......',1x,
     .       '.......',1x,'.......',1x,'.......',1x,'.......',1x,
     .       '.......',1x,'.......')
  135 format('# xraw  ',1x,'  yraw  ',1x,' mfit  ',1x,'sky15-25',4x,
     .       'map_r01',1x,'map_r02',1x,'map_r03',1x,'map_r04',1x,
     .       'map_r05',1x,'map_r06',1x,'map_r07',1x,'map_r08',1x,
     .       'map_r09',1x,'map_r10')
 
 
      return
      end
 
 
 
 
      subroutine find_qfit(xr,yr,sr,qfit,pixc)
      implicit none
      real*8 xr,yr,sr
      real qfit
      real pixc(4096,4096)
 
      integer i,j
 
      real plist(25)
      real flist(25)
      integer L
      real dx, dy
      real rpsf_photijk
      real ptot,ftot
      real flx, perr
 
      L = 0
      ptot = 0.
      ftot = 0.
      do i = int(xr+0.5)-2,int(xr+0.5)+2
      do j = int(yr+0.5)-2,int(yr+0.5)+2
         L = L + 1
         dx = i-xr
         dy = j-yr
         plist(L) = pixc(i,j)-sr
         flist(L) = rpsf_photijk(dx,dy,i,j)
         ptot = ptot + plist(L)
         ftot = ftot + flist(L)
         enddo
         enddo
      flx = ptot/ftot
      perr = 0.
      do L = 1, 25
         perr = perr + abs(plist(L)-flx*flist(L))
         enddo
      qfit = perr/flx
      if (qfit.gt.9.99) qfit = 9.99
      if (qfit.lt.0.00) qfit = 0.00
 
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
 
 
 
 
 
 
 
      subroutine wfc3uv_gc(xr,yr,xc,yc,kin)
      implicit none
 
      real*8 xr,yr
      real*8 xc,yc
      integer kin
 
      real xru,yru
      real rx, ry
      real fx, fy
 
      real x, y
      integer kinu
 
 
 
c---------------------------
c
c bottom chip actuall chip#2
c
      real poly1x(10)
      data poly1x / 2048.001, 2052.406,   64.737,  12.024, -5.873,
     .                -0.147,    0.205,   -0.038,   0.051,  0.012/
      real poly1y(10)
      data poly1y / 1026.000,    0.000, 1024.000,  -0.134,  5.906,
     .                -3.184,    0.032,    0.058,  -0.021,  0.005/
 
c---------------------------
c
c bottom chip actuall chip#2
c
      real poly2x(10)
      data poly2x / 2176.839, 2041.115,   62.891,  11.830, -5.883,
     .                -0.294,    0.006,   -0.072,  -0.057, -0.068/
      real poly2y(10)
      data poly2y / 3095.615,   13.118, 1011.904,   0.255,  6.046,
     .                -2.779,    0.270,    0.210,   0.100,  0.125/
 
      real*8 xloc,  yloc
      real*8 xlocu, ylocu
      real*8 dxadj, dyadj
 
      real xarr(65,66)
      common /wfc3xarr_/xarr
 
      real yarr(65,66)
      common /wfc3yarr_/yarr
 
      real    ri, rj
      integer i,  j
      real    fi, fj
 
      xc = xr
      yc = yr
 
      xru = xr
      yru = yr
 
      kinu = abs(kin)
      if (kinu.eq.0) then
         if (yr.lt.2048) then
            yru  = yr
            kinu = 2
            endif
         if (yr.ge.2048) then
            yru = yr-2048.00
            kinu = 1
            endif
         endif
 
      xloc = xr
      yloc = yr
 
 
      ylocu = yloc-int((yloc-1)/2048)*2048
      ri = 1 + (xloc/64)
      rj = 1 + (ylocu/64)
      if (yloc.gt.2048) rj = rj + 33
 
      if (ri.lt. 1.01) ri =  1.01
      if (rj.lt. 1.01) rj =  1.01
      if (ri.gt.64.99) ri = 64.99
      if (rj.gt.65.99) rj = 65.99
 
      i = int(ri)
      j = int(rj)
      fi = ri-i
      fj = rj-j
 
      if (i.gt.64) i = 64
      if (j.gt.65) j = 65
 
      dxadj = (1-fi)*(1-fj)*xarr(i  ,j  )
     .      + (1-fi)*( fj )*xarr(i  ,j+1)
     .      + ( fi )*(1-fj)*xarr(i+1,j  )
     .      + ( fi )*( fj )*xarr(i+1,j+1)
 
      dyadj = (1-fi)*(1-fj)*yarr(i  ,j  )
     .      + (1-fi)*( fj )*yarr(i  ,j+1)
     .      + ( fi )*(1-fj)*yarr(i+1,j  )
     .      + ( fi )*( fj )*yarr(i+1,j+1)
 
 
      if (kinu.eq.1) then  ! top chip
          x = (xru-2048.00)/2048.00
          y = (yru-1024.00)/1024.00
          xc = poly2x(01)
     .       + poly2x(02)*x
     .       + poly2x(03)*y
     .       + poly2x(04)*x*x
     .       + poly2x(05)*x*y
     .       + poly2x(06)*y*y
     .       + poly2x(07)*x*x*x
     .       + poly2x(08)*x*x*y
     .       + poly2x(09)*x*y*y
     .       + poly2x(10)*y*y*y
          yc = poly2y(01)
     .       + poly2y(02)*x
     .       + poly2y(03)*y
     .       + poly2y(04)*x*x
     .       + poly2y(05)*x*y
     .       + poly2y(06)*y*y
     .       + poly2y(07)*x*x*x
     .       + poly2y(08)*x*x*y
     .       + poly2y(09)*x*y*y
     .       + poly2y(10)*y*y*y
 
          xc = xc - dxadj
          yc = yc - dyadj
          return
          endif
 
      if (kinu.eq.2) then
         x = (xru-2048.00)/2048.00
         y = (yru-1024.00)/1024.00
          xc = poly1x(01)
     .       + poly1x(02)*x
     .       + poly1x(03)*y
     .       + poly1x(04)*x*x
     .       + poly1x(05)*x*y
     .       + poly1x(06)*y*y
     .       + poly1x(07)*x*x*x
     .       + poly1x(08)*x*x*y
     .       + poly1x(09)*x*y*y
     .       + poly1x(10)*y*y*y
          yc = poly1y(01)
     .       + poly1y(02)*x
     .       + poly1y(03)*y
     .       + poly1y(04)*x*x
     .       + poly1y(05)*x*y
     .       + poly1y(06)*y*y
     .       + poly1y(07)*x*x*x
     .       + poly1y(08)*x*x*y
     .       + poly1y(09)*x*y*y
     .       + poly1y(10)*y*y*y
 
          xc = xc - dxadj
          yc = yc - dyadj
          return
          endif
 
      stop 'wfc3uv_gc'
      end
 
 
 
      block data wfc3_xarr
      real xarr(65,66)
      common /wfc3xarr_/xarr
      data xarr /
c XARR; ROW j = 01
     .-0.027,-0.023,-0.018,-0.014,-0.005,-0.007,-0.007,-0.006,
     .-0.010,-0.009,-0.002, 0.003, 0.011, 0.012, 0.009, 0.007,
     . 0.008, 0.011, 0.011, 0.011, 0.011, 0.006, 0.001,-0.003,
     .-0.007,-0.008,-0.013,-0.021,-0.028,-0.034,-0.033,-0.029,
     .-0.020,-0.016,-0.014,-0.019,-0.024,-0.029,-0.030,-0.028,
     .-0.028,-0.025,-0.018,-0.011,-0.005, 0.002, 0.006, 0.011,
     . 0.015, 0.017, 0.021, 0.019, 0.020, 0.015, 0.019, 0.021,
     . 0.020, 0.016, 0.012, 0.008, 0.010, 0.016, 0.019, 0.023,
     . 0.028,
c XARR; ROW j = 02
     .-0.022,-0.018,-0.014,-0.010,-0.002,-0.002,-0.001,-0.001,
     .-0.004,-0.002, 0.003, 0.012, 0.019, 0.021, 0.020, 0.019,
     . 0.020, 0.023, 0.022, 0.023, 0.022, 0.016, 0.009, 0.005,
     . 0.002, 0.001,-0.005,-0.013,-0.021,-0.025,-0.025,-0.022,
     .-0.016,-0.013,-0.011,-0.017,-0.023,-0.030,-0.032,-0.030,
     .-0.029,-0.026,-0.021,-0.014,-0.009,-0.002,-0.001, 0.004,
     . 0.008, 0.011, 0.014, 0.013, 0.013, 0.012, 0.014, 0.014,
     . 0.013, 0.009, 0.005, 0.001, 0.001, 0.006, 0.009, 0.013,
     . 0.017,
c XARR; ROW j = 03
     .-0.019,-0.015,-0.011,-0.008,-0.001, 0.001, 0.002, 0.002,
     .-0.002, 0.001, 0.007, 0.018, 0.024, 0.027, 0.026, 0.026,
     . 0.027, 0.030, 0.029, 0.030, 0.030, 0.023, 0.015, 0.011,
     . 0.008, 0.008, 0.000,-0.008,-0.016,-0.019,-0.019,-0.017,
     .-0.013,-0.010,-0.008,-0.016,-0.023,-0.030,-0.034,-0.032,
     .-0.031,-0.029,-0.024,-0.018,-0.013,-0.006,-0.008,-0.003,
     . 0.001, 0.005, 0.007, 0.007, 0.006, 0.008, 0.008, 0.008,
     . 0.006, 0.002,-0.002,-0.006,-0.007,-0.003, 0.000, 0.004,
     . 0.009,
c XARR; ROW j = 04
     .-0.018,-0.015,-0.011,-0.007,-0.001, 0.002, 0.002, 0.003,
     .-0.001, 0.002, 0.007, 0.020, 0.025, 0.029, 0.029, 0.029,
     . 0.031, 0.034, 0.033, 0.035, 0.034, 0.027, 0.019, 0.015,
     . 0.012, 0.011, 0.004,-0.004,-0.012,-0.014,-0.016,-0.014,
     .-0.011,-0.008,-0.007,-0.016,-0.024,-0.032,-0.036,-0.034,
     .-0.034,-0.032,-0.028,-0.023,-0.019,-0.011,-0.014,-0.010,
     .-0.006,-0.002,-0.001, 0.001,-0.001, 0.003, 0.003, 0.002,
     . 0.001,-0.003,-0.008,-0.011,-0.013,-0.009,-0.006,-0.003,
     . 0.002,
c XARR; ROW j = 05
     .-0.020,-0.015,-0.010,-0.006,-0.003, 0.000, 0.000, 0.000,
     .-0.003,-0.001, 0.005, 0.016, 0.024, 0.027, 0.025, 0.024,
     . 0.027, 0.030, 0.032, 0.033, 0.032, 0.027, 0.021, 0.014,
     . 0.010, 0.008, 0.005,-0.003,-0.011,-0.016,-0.016,-0.014,
     .-0.009,-0.009,-0.011,-0.019,-0.028,-0.034,-0.038,-0.040,
     .-0.042,-0.040,-0.035,-0.030,-0.026,-0.021,-0.020,-0.018,
     .-0.013,-0.012,-0.009,-0.008,-0.007,-0.006,-0.003,-0.004,
     .-0.005,-0.009,-0.012,-0.016,-0.020,-0.014,-0.012,-0.010,
     .-0.009,
c XARR; ROW j = 06
     .-0.029,-0.024,-0.018,-0.013,-0.008,-0.006,-0.005,-0.005,
     .-0.008,-0.008, 0.002, 0.012, 0.019, 0.021, 0.019, 0.018,
     . 0.022, 0.025, 0.027, 0.030, 0.030, 0.026, 0.020, 0.014,
     . 0.009, 0.008, 0.004,-0.003,-0.008,-0.013,-0.014,-0.014,
     .-0.010,-0.009,-0.011,-0.019,-0.028,-0.037,-0.041,-0.041,
     .-0.045,-0.044,-0.042,-0.038,-0.032,-0.027,-0.028,-0.024,
     .-0.022,-0.019,-0.016,-0.015,-0.013,-0.011,-0.009,-0.005,
     .-0.005,-0.011,-0.016,-0.022,-0.022,-0.014,-0.013,-0.014,
     .-0.015,
c XARR; ROW j = 07
     .-0.037,-0.032,-0.028,-0.023,-0.017,-0.014,-0.012,-0.011,
     .-0.015,-0.014,-0.004, 0.006, 0.013, 0.016, 0.014, 0.014,
     . 0.019, 0.022, 0.025, 0.030, 0.031, 0.025, 0.019, 0.015,
     . 0.010, 0.010, 0.005,-0.001,-0.006,-0.009,-0.012,-0.011,
     .-0.010,-0.010,-0.012,-0.018,-0.028,-0.036,-0.042,-0.043,
     .-0.046,-0.046,-0.044,-0.041,-0.037,-0.031,-0.030,-0.029,
     .-0.026,-0.025,-0.021,-0.020,-0.017,-0.015,-0.011,-0.009,
     .-0.007,-0.012,-0.015,-0.021,-0.022,-0.013,-0.013,-0.015,
     .-0.018,
c XARR; ROW j = 08
     .-0.045,-0.041,-0.036,-0.031,-0.023,-0.021,-0.019,-0.018,
     .-0.019,-0.017,-0.011,-0.001, 0.008, 0.013, 0.013, 0.015,
     . 0.016, 0.022, 0.025, 0.032, 0.032, 0.026, 0.019, 0.016,
     . 0.012, 0.014, 0.007, 0.002,-0.002,-0.004,-0.006,-0.008,
     .-0.011,-0.010,-0.013,-0.019,-0.025,-0.035,-0.041,-0.044,
     .-0.048,-0.048,-0.047,-0.043,-0.038,-0.035,-0.034,-0.033,
     .-0.032,-0.029,-0.027,-0.025,-0.021,-0.019,-0.014,-0.012,
     .-0.009,-0.012,-0.016,-0.020,-0.021,-0.016,-0.016,-0.016,
     .-0.017,
c XARR; ROW j = 09
     .-0.050,-0.046,-0.042,-0.037,-0.028,-0.026,-0.022,-0.020,
     .-0.022,-0.019,-0.012,-0.002, 0.009, 0.016, 0.017, 0.018,
     . 0.020, 0.025, 0.029, 0.034, 0.033, 0.028, 0.021, 0.017,
     . 0.015, 0.016, 0.011, 0.008, 0.003,-0.002,-0.004,-0.006,
     .-0.010,-0.012,-0.014,-0.021,-0.026,-0.035,-0.042,-0.045,
     .-0.048,-0.049,-0.047,-0.045,-0.042,-0.039,-0.040,-0.037,
     .-0.036,-0.034,-0.031,-0.031,-0.027,-0.021,-0.018,-0.012,
     .-0.009,-0.011,-0.014,-0.018,-0.020,-0.014,-0.013,-0.014,
     .-0.015,
c XARR; ROW j = 10
     .-0.051,-0.049,-0.046,-0.041,-0.032,-0.028,-0.024,-0.020,
     .-0.022,-0.020,-0.013,-0.001, 0.011, 0.019, 0.020, 0.023,
     . 0.025, 0.029, 0.031, 0.035, 0.033, 0.029, 0.022, 0.016,
     . 0.016, 0.018, 0.015, 0.010, 0.004,-0.001,-0.003,-0.009,
     .-0.012,-0.016,-0.016,-0.024,-0.030,-0.037,-0.043,-0.045,
     .-0.047,-0.049,-0.048,-0.047,-0.044,-0.040,-0.040,-0.040,
     .-0.039,-0.038,-0.035,-0.033,-0.028,-0.025,-0.019,-0.013,
     .-0.010,-0.010,-0.012,-0.015,-0.016,-0.011,-0.011,-0.012,
     .-0.015,
c XARR; ROW j = 11
     .-0.054,-0.050,-0.046,-0.040,-0.036,-0.029,-0.023,-0.021,
     .-0.020,-0.019,-0.010, 0.003, 0.015, 0.024, 0.026, 0.027,
     . 0.032, 0.033, 0.034, 0.037, 0.034, 0.030, 0.024, 0.020,
     . 0.020, 0.021, 0.018, 0.013, 0.006, 0.002,-0.002,-0.007,
     .-0.012,-0.015,-0.020,-0.023,-0.029,-0.037,-0.044,-0.046,
     .-0.050,-0.048,-0.048,-0.048,-0.044,-0.041,-0.040,-0.040,
     .-0.040,-0.040,-0.037,-0.034,-0.030,-0.026,-0.020,-0.015,
     .-0.010,-0.010,-0.011,-0.013,-0.015,-0.006,-0.007,-0.010,
     .-0.015,
c XARR; ROW j = 12
     .-0.052,-0.052,-0.051,-0.047,-0.041,-0.031,-0.023,-0.019,
     .-0.019,-0.017,-0.006, 0.008, 0.021, 0.031, 0.032, 0.035,
     . 0.038, 0.038, 0.036, 0.041, 0.038, 0.033, 0.026, 0.023,
     . 0.021, 0.024, 0.019, 0.013, 0.007, 0.002,-0.001,-0.006,
     .-0.011,-0.014,-0.018,-0.023,-0.029,-0.039,-0.045,-0.048,
     .-0.051,-0.049,-0.050,-0.049,-0.045,-0.042,-0.041,-0.039,
     .-0.041,-0.039,-0.037,-0.035,-0.031,-0.026,-0.021,-0.015,
     .-0.010,-0.008,-0.010,-0.012,-0.011,-0.003,-0.002,-0.004,
     .-0.007,
c XARR; ROW j = 13
     .-0.057,-0.057,-0.054,-0.050,-0.043,-0.035,-0.027,-0.020,
     .-0.020,-0.015,-0.004, 0.010, 0.024, 0.036, 0.038, 0.039,
     . 0.041, 0.041, 0.041, 0.042, 0.037, 0.034, 0.030, 0.027,
     . 0.025, 0.025, 0.018, 0.015, 0.008, 0.003,-0.001,-0.004,
     .-0.010,-0.012,-0.016,-0.022,-0.030,-0.037,-0.045,-0.047,
     .-0.051,-0.050,-0.051,-0.048,-0.044,-0.039,-0.040,-0.039,
     .-0.039,-0.039,-0.037,-0.034,-0.031,-0.026,-0.021,-0.017,
     .-0.011,-0.008,-0.007,-0.008,-0.007,-0.002, 0.001, 0.004,
     . 0.007,
c XARR; ROW j = 14
     .-0.052,-0.051,-0.048,-0.043,-0.037,-0.029,-0.020,-0.014,
     .-0.011,-0.007, 0.002, 0.016, 0.031, 0.041, 0.045, 0.044,
     . 0.043, 0.044, 0.042, 0.042, 0.040, 0.038, 0.035, 0.033,
     . 0.033, 0.033, 0.025, 0.018, 0.012, 0.007, 0.006, 0.002,
     .-0.004,-0.005,-0.008,-0.016,-0.024,-0.033,-0.040,-0.042,
     .-0.047,-0.047,-0.046,-0.041,-0.035,-0.032,-0.033,-0.030,
     .-0.029,-0.029,-0.027,-0.027,-0.024,-0.021,-0.019,-0.016,
     .-0.011,-0.009,-0.007,-0.006,-0.006,-0.002, 0.003, 0.008,
     . 0.015,
c XARR; ROW j = 15
     .-0.044,-0.041,-0.037,-0.032,-0.025,-0.017,-0.008,-0.002,
     . 0.000, 0.004, 0.014, 0.026, 0.039, 0.048, 0.050, 0.048,
     . 0.047, 0.046, 0.046, 0.046, 0.043, 0.044, 0.043, 0.042,
     . 0.042, 0.041, 0.033, 0.027, 0.020, 0.013, 0.012, 0.010,
     . 0.008, 0.007, 0.003,-0.007,-0.016,-0.025,-0.031,-0.032,
     .-0.038,-0.039,-0.035,-0.029,-0.023,-0.017,-0.018,-0.016,
     .-0.014,-0.015,-0.014,-0.013,-0.012,-0.012,-0.012,-0.010,
     .-0.008,-0.007,-0.006,-0.004,-0.003, 0.002, 0.007, 0.014,
     . 0.023,
c XARR; ROW j = 16
     .-0.030,-0.028,-0.024,-0.019,-0.009,-0.002, 0.007, 0.013,
     . 0.015, 0.019, 0.026, 0.037, 0.047, 0.054, 0.055, 0.053,
     . 0.052, 0.052, 0.049, 0.049, 0.049, 0.051, 0.052, 0.052,
     . 0.052, 0.050, 0.042, 0.035, 0.028, 0.022, 0.020, 0.020,
     . 0.018, 0.017, 0.013, 0.004,-0.006,-0.013,-0.021,-0.021,
     .-0.028,-0.026,-0.022,-0.016,-0.009, 0.000,-0.002, 0.000,
     .-0.001,-0.002,-0.001, 0.002, 0.000,-0.005,-0.005,-0.007,
     .-0.005,-0.005,-0.004,-0.004,-0.001, 0.000, 0.008, 0.018,
     . 0.032,
c XARR; ROW j = 17
     .-0.015,-0.012,-0.007,-0.002, 0.005, 0.015, 0.023, 0.030,
     . 0.030, 0.032, 0.040, 0.049, 0.057, 0.061, 0.062, 0.060,
     . 0.060, 0.058, 0.057, 0.057, 0.058, 0.058, 0.060, 0.062,
     . 0.063, 0.059, 0.051, 0.044, 0.036, 0.031, 0.029, 0.029,
     . 0.028, 0.027, 0.021, 0.015, 0.006,-0.004,-0.010,-0.012,
     .-0.019,-0.016,-0.014,-0.007, 0.002, 0.009, 0.007, 0.009,
     . 0.008, 0.008, 0.007, 0.008, 0.006, 0.002,-0.002,-0.006,
     .-0.006,-0.006,-0.006,-0.006,-0.005,-0.003, 0.006, 0.020,
     . 0.038,
c XARR; ROW j = 18
     .-0.003, 0.000, 0.004, 0.009, 0.017, 0.027, 0.033, 0.041,
     . 0.043, 0.044, 0.050, 0.057, 0.064, 0.070, 0.069, 0.067,
     . 0.069, 0.068, 0.065, 0.064, 0.063, 0.063, 0.065, 0.066,
     . 0.067, 0.063, 0.055, 0.048, 0.039, 0.033, 0.032, 0.030,
     . 0.030, 0.028, 0.024, 0.015, 0.007,-0.003,-0.008,-0.011,
     .-0.018,-0.016,-0.013,-0.006, 0.003, 0.009, 0.005, 0.006,
     . 0.006, 0.005, 0.005, 0.005, 0.004,-0.002,-0.006,-0.011,
     .-0.011,-0.010,-0.010,-0.012,-0.009,-0.009, 0.002, 0.019,
     . 0.041,
c XARR; ROW j = 19
     . 0.000, 0.007, 0.014, 0.020, 0.022, 0.029, 0.035, 0.043,
     . 0.047, 0.049, 0.052, 0.059, 0.065, 0.070, 0.071, 0.070,
     . 0.072, 0.072, 0.072, 0.066, 0.064, 0.063, 0.064, 0.066,
     . 0.066, 0.060, 0.052, 0.046, 0.038, 0.030, 0.027, 0.026,
     . 0.024, 0.022, 0.019, 0.009, 0.001,-0.006,-0.013,-0.014,
     .-0.024,-0.025,-0.019,-0.013,-0.005, 0.001,-0.006,-0.008,
     .-0.006,-0.007,-0.005,-0.003,-0.006,-0.010,-0.014,-0.019,
     .-0.020,-0.015,-0.016,-0.015,-0.015,-0.014,-0.002, 0.018,
     . 0.043,
c XARR; ROW j = 20
     . 0.010, 0.015, 0.020, 0.024, 0.026, 0.032, 0.040, 0.046,
     . 0.050, 0.052, 0.055, 0.062, 0.068, 0.070, 0.071, 0.071,
     . 0.073, 0.073, 0.072, 0.065, 0.063, 0.059, 0.059, 0.060,
     . 0.060, 0.056, 0.050, 0.041, 0.033, 0.025, 0.021, 0.017,
     . 0.015, 0.014, 0.011, 0.002,-0.006,-0.013,-0.019,-0.021,
     .-0.033,-0.033,-0.027,-0.021,-0.015,-0.010,-0.018,-0.018,
     .-0.016,-0.015,-0.015,-0.010,-0.013,-0.017,-0.021,-0.025,
     .-0.026,-0.022,-0.021,-0.019,-0.018,-0.017,-0.004, 0.016,
     . 0.043,
c XARR; ROW j = 21
     . 0.009, 0.014, 0.019, 0.024, 0.027, 0.032, 0.037, 0.044,
     . 0.048, 0.049, 0.053, 0.059, 0.065, 0.067, 0.066, 0.066,
     . 0.068, 0.070, 0.070, 0.064, 0.059, 0.055, 0.054, 0.054,
     . 0.054, 0.050, 0.042, 0.035, 0.027, 0.018, 0.014, 0.010,
     . 0.008, 0.006, 0.003,-0.006,-0.012,-0.019,-0.026,-0.027,
     .-0.039,-0.039,-0.034,-0.029,-0.022,-0.017,-0.026,-0.026,
     .-0.025,-0.023,-0.019,-0.015,-0.018,-0.023,-0.025,-0.028,
     .-0.027,-0.024,-0.024,-0.023,-0.018,-0.018,-0.004, 0.017,
     . 0.045,
c XARR; ROW j = 22
     . 0.005, 0.010, 0.015, 0.020, 0.022, 0.027, 0.031, 0.038,
     . 0.043, 0.045, 0.050, 0.055, 0.059, 0.061, 0.061, 0.060,
     . 0.063, 0.065, 0.067, 0.062, 0.058, 0.053, 0.051, 0.050,
     . 0.050, 0.044, 0.037, 0.029, 0.022, 0.013, 0.008, 0.005,
     . 0.002, 0.000,-0.003,-0.010,-0.015,-0.023,-0.031,-0.034,
     .-0.043,-0.043,-0.038,-0.033,-0.027,-0.023,-0.030,-0.030,
     .-0.029,-0.027,-0.023,-0.018,-0.020,-0.025,-0.027,-0.029,
     .-0.027,-0.024,-0.023,-0.022,-0.016,-0.015,-0.001, 0.021,
     . 0.049,
c XARR; ROW j = 23
     . 0.002, 0.008, 0.013, 0.017, 0.019, 0.022, 0.026, 0.032,
     . 0.037, 0.041, 0.046, 0.051, 0.056, 0.057, 0.058, 0.058,
     . 0.060, 0.062, 0.066, 0.061, 0.059, 0.052, 0.048, 0.048,
     . 0.047, 0.041, 0.033, 0.025, 0.018, 0.010, 0.006, 0.001,
     .-0.003,-0.004,-0.006,-0.011,-0.017,-0.029,-0.034,-0.037,
     .-0.045,-0.044,-0.038,-0.034,-0.028,-0.025,-0.031,-0.031,
     .-0.030,-0.029,-0.024,-0.020,-0.021,-0.025,-0.027,-0.028,
     .-0.027,-0.023,-0.020,-0.018,-0.013,-0.012, 0.003, 0.025,
     . 0.053,
c XARR; ROW j = 24
     .-0.001, 0.005, 0.010, 0.013, 0.016, 0.018, 0.020, 0.026,
     . 0.033, 0.037, 0.042, 0.048, 0.053, 0.057, 0.057, 0.056,
     . 0.058, 0.062, 0.065, 0.063, 0.058, 0.054, 0.049, 0.046,
     . 0.044, 0.039, 0.030, 0.024, 0.016, 0.009, 0.006, 0.000,
     .-0.005,-0.006,-0.008,-0.013,-0.020,-0.031,-0.038,-0.040,
     .-0.046,-0.044,-0.036,-0.032,-0.026,-0.025,-0.031,-0.029,
     .-0.027,-0.029,-0.024,-0.022,-0.023,-0.025,-0.028,-0.029,
     .-0.025,-0.019,-0.017,-0.016,-0.011,-0.010, 0.006, 0.030,
     . 0.061,
c XARR; ROW j = 25
     .-0.002, 0.005, 0.010, 0.014, 0.012, 0.014, 0.017, 0.024,
     . 0.029, 0.035, 0.038, 0.048, 0.054, 0.058, 0.058, 0.057,
     . 0.057, 0.062, 0.066, 0.063, 0.060, 0.055, 0.051, 0.047,
     . 0.044, 0.035, 0.027, 0.021, 0.014, 0.005, 0.002,-0.004,
     .-0.006,-0.009,-0.013,-0.019,-0.024,-0.035,-0.043,-0.046,
     .-0.048,-0.042,-0.034,-0.030,-0.024,-0.024,-0.027,-0.027,
     .-0.026,-0.026,-0.023,-0.021,-0.023,-0.027,-0.029,-0.028,
     .-0.022,-0.015,-0.014,-0.011,-0.005, 0.000, 0.017, 0.040,
     . 0.070,
c XARR; ROW j = 26
     .-0.002, 0.004, 0.009, 0.012, 0.011, 0.013, 0.015, 0.021,
     . 0.028, 0.031, 0.037, 0.045, 0.053, 0.056, 0.056, 0.055,
     . 0.057, 0.059, 0.063, 0.062, 0.060, 0.056, 0.050, 0.044,
     . 0.042, 0.034, 0.025, 0.016, 0.009, 0.002,-0.002,-0.008,
     .-0.010,-0.014,-0.017,-0.024,-0.031,-0.041,-0.048,-0.050,
     .-0.051,-0.043,-0.034,-0.028,-0.024,-0.022,-0.025,-0.023,
     .-0.022,-0.021,-0.021,-0.017,-0.021,-0.026,-0.026,-0.026,
     .-0.020,-0.012,-0.011,-0.006, 0.001, 0.009, 0.025, 0.048,
     . 0.076,
c XARR; ROW j = 27
     . 0.001, 0.005, 0.009, 0.011, 0.009, 0.014, 0.016, 0.020,
     . 0.027, 0.029, 0.035, 0.043, 0.049, 0.053, 0.055, 0.053,
     . 0.054, 0.057, 0.061, 0.061, 0.058, 0.054, 0.049, 0.043,
     . 0.039, 0.032, 0.021, 0.012, 0.006,-0.001,-0.005,-0.011,
     .-0.016,-0.019,-0.023,-0.028,-0.034,-0.046,-0.053,-0.054,
     .-0.053,-0.045,-0.035,-0.031,-0.024,-0.021,-0.023,-0.021,
     .-0.018,-0.018,-0.015,-0.014,-0.018,-0.024,-0.025,-0.024,
     .-0.018,-0.010,-0.007,-0.001, 0.008, 0.015, 0.033, 0.056,
     . 0.084,
c XARR; ROW j = 28
     . 0.005, 0.008, 0.011, 0.012, 0.012, 0.015, 0.019, 0.022,
     . 0.027, 0.031, 0.034, 0.041, 0.046, 0.051, 0.052, 0.052,
     . 0.051, 0.053, 0.057, 0.059, 0.057, 0.050, 0.045, 0.040,
     . 0.035, 0.027, 0.017, 0.007, 0.000,-0.007,-0.013,-0.018,
     .-0.023,-0.026,-0.030,-0.035,-0.043,-0.052,-0.057,-0.058,
     .-0.057,-0.049,-0.040,-0.032,-0.025,-0.021,-0.022,-0.018,
     .-0.016,-0.014,-0.010,-0.009,-0.013,-0.018,-0.024,-0.022,
     .-0.015,-0.006,-0.003, 0.004, 0.013, 0.024, 0.042, 0.066,
     . 0.095,
c XARR; ROW j = 29
     . 0.012, 0.014, 0.016, 0.018, 0.020, 0.020, 0.022, 0.025,
     . 0.029, 0.033, 0.035, 0.039, 0.045, 0.049, 0.050, 0.050,
     . 0.049, 0.050, 0.051, 0.053, 0.050, 0.044, 0.039, 0.035,
     . 0.029, 0.025, 0.012, 0.003,-0.005,-0.013,-0.019,-0.024,
     .-0.031,-0.034,-0.035,-0.040,-0.047,-0.054,-0.059,-0.059,
     .-0.058,-0.051,-0.042,-0.035,-0.027,-0.022,-0.021,-0.015,
     .-0.012,-0.009,-0.002,-0.002,-0.006,-0.012,-0.016,-0.016,
     .-0.010,-0.002, 0.002, 0.012, 0.022, 0.032, 0.051, 0.076,
     . 0.106,
c XARR; ROW j = 30
     . 0.016, 0.020, 0.023, 0.025, 0.022, 0.026, 0.027, 0.031,
     . 0.037, 0.035, 0.036, 0.043, 0.045, 0.051, 0.051, 0.048,
     . 0.047, 0.047, 0.045, 0.046, 0.042, 0.037, 0.033, 0.028,
     . 0.025, 0.020, 0.010, 0.001,-0.008,-0.015,-0.021,-0.029,
     .-0.035,-0.037,-0.037,-0.041,-0.046,-0.053,-0.056,-0.056,
     .-0.058,-0.052,-0.043,-0.034,-0.024,-0.019,-0.019,-0.010,
     .-0.003,-0.001, 0.005, 0.005, 0.002,-0.002,-0.009,-0.009,
     .-0.003, 0.007, 0.012, 0.020, 0.030, 0.045, 0.064, 0.088,
     . 0.117,
c XARR; ROW j = 31
     . 0.023, 0.027, 0.030, 0.033, 0.031, 0.035, 0.036, 0.040,
     . 0.045, 0.044, 0.044, 0.049, 0.051, 0.054, 0.054, 0.049,
     . 0.047, 0.046, 0.044, 0.044, 0.039, 0.035, 0.031, 0.027,
     . 0.024, 0.018, 0.009, 0.000,-0.009,-0.017,-0.022,-0.030,
     .-0.035,-0.037,-0.037,-0.039,-0.043,-0.049,-0.052,-0.052,
     .-0.054,-0.047,-0.038,-0.029,-0.019,-0.014,-0.011,-0.003,
     . 0.004, 0.007, 0.013, 0.014, 0.012, 0.007, 0.002, 0.002,
     . 0.008, 0.017, 0.022, 0.030, 0.040, 0.054, 0.073, 0.097,
     . 0.126,
c XARR; ROW j = 32
     . 0.032, 0.036, 0.040, 0.043, 0.042, 0.046, 0.048, 0.052,
     . 0.056, 0.055, 0.054, 0.058, 0.059, 0.060, 0.058, 0.053,
     . 0.048, 0.046, 0.045, 0.044, 0.040, 0.035, 0.032, 0.028,
     . 0.025, 0.018, 0.009, 0.000,-0.009,-0.017,-0.022,-0.029,
     .-0.032,-0.034,-0.034,-0.034,-0.036,-0.042,-0.045,-0.046,
     .-0.046,-0.039,-0.030,-0.021,-0.011,-0.005,-0.001, 0.006,
     . 0.012, 0.016, 0.023, 0.026, 0.023, 0.018, 0.014, 0.015,
     . 0.022, 0.029, 0.033, 0.040, 0.051, 0.063, 0.082, 0.106,
     . 0.135,
c XARR; ROW j = 33
     . 0.043, 0.048, 0.052, 0.055, 0.057, 0.059, 0.062, 0.067,
     . 0.069, 0.071, 0.069, 0.069, 0.070, 0.069, 0.065, 0.058,
     . 0.049, 0.047, 0.047, 0.046, 0.043, 0.037, 0.035, 0.031,
     . 0.027, 0.020, 0.010, 0.002,-0.007,-0.015,-0.021,-0.027,
     .-0.027,-0.028,-0.028,-0.027,-0.027,-0.031,-0.035,-0.038,
     .-0.036,-0.028,-0.019,-0.010,-0.001, 0.007, 0.013, 0.018,
     . 0.022, 0.026, 0.034, 0.039, 0.037, 0.031, 0.030, 0.032,
     . 0.038, 0.044, 0.046, 0.052, 0.062, 0.072, 0.091, 0.115,
     . 0.144,
c XARR; ROW j = 34
     .-0.169,-0.167,-0.166,-0.164,-0.159,-0.155,-0.150,-0.142,
     .-0.141,-0.140,-0.129,-0.116,-0.103,-0.092,-0.088,-0.085,
     .-0.078,-0.077,-0.075,-0.080,-0.070,-0.059,-0.043,-0.029,
     .-0.017,-0.023,-0.023,-0.025,-0.031,-0.033,-0.037,-0.041,
     .-0.045,-0.048,-0.054,-0.056,-0.060,-0.061,-0.063,-0.056,
     .-0.062,-0.057,-0.050,-0.034,-0.020,-0.005,-0.004, 0.007,
     . 0.016, 0.028, 0.038, 0.046, 0.058, 0.067, 0.075, 0.084,
     . 0.089, 0.086, 0.085, 0.081, 0.083, 0.077, 0.081, 0.088,
     . 0.097,
c XARR; ROW j = 35
     .-0.153,-0.152,-0.150,-0.147,-0.144,-0.140,-0.136,-0.130,
     .-0.129,-0.127,-0.117,-0.107,-0.094,-0.085,-0.083,-0.080,
     .-0.075,-0.073,-0.070,-0.077,-0.068,-0.056,-0.042,-0.028,
     .-0.018,-0.024,-0.025,-0.030,-0.037,-0.041,-0.045,-0.051,
     .-0.056,-0.058,-0.063,-0.065,-0.068,-0.068,-0.068,-0.062,
     .-0.067,-0.062,-0.053,-0.038,-0.024,-0.011,-0.008, 0.001,
     . 0.011, 0.022, 0.031, 0.039, 0.050, 0.059, 0.067, 0.075,
     . 0.079, 0.079, 0.078, 0.077, 0.079, 0.079, 0.084, 0.091,
     . 0.101,
c XARR; ROW j = 36
     .-0.141,-0.139,-0.136,-0.134,-0.131,-0.127,-0.124,-0.119,
     .-0.118,-0.117,-0.108,-0.099,-0.087,-0.079,-0.078,-0.077,
     .-0.072,-0.070,-0.066,-0.074,-0.065,-0.053,-0.040,-0.027,
     .-0.018,-0.025,-0.027,-0.034,-0.042,-0.047,-0.051,-0.058,
     .-0.063,-0.065,-0.069,-0.071,-0.073,-0.073,-0.070,-0.066,
     .-0.070,-0.064,-0.054,-0.040,-0.026,-0.015,-0.011,-0.003,
     . 0.006, 0.017, 0.026, 0.034, 0.044, 0.053, 0.061, 0.069,
     . 0.073, 0.075, 0.073, 0.075, 0.078, 0.080, 0.086, 0.094,
     . 0.104,
c XARR; ROW j = 37
     .-0.130,-0.128,-0.125,-0.122,-0.120,-0.116,-0.114,-0.109,
     .-0.110,-0.108,-0.100,-0.094,-0.083,-0.075,-0.075,-0.074,
     .-0.070,-0.068,-0.063,-0.070,-0.062,-0.051,-0.038,-0.025,
     .-0.018,-0.026,-0.028,-0.037,-0.046,-0.052,-0.056,-0.063,
     .-0.068,-0.068,-0.071,-0.074,-0.074,-0.074,-0.070,-0.066,
     .-0.070,-0.064,-0.053,-0.040,-0.027,-0.016,-0.012,-0.005,
     . 0.004, 0.015, 0.024, 0.032, 0.041, 0.050, 0.059, 0.067,
     . 0.069, 0.072, 0.071, 0.074, 0.077, 0.081, 0.087, 0.096,
     . 0.107,
c XARR; ROW j = 38
     .-0.124,-0.119,-0.115,-0.112,-0.111,-0.108,-0.106,-0.101,
     .-0.102,-0.103,-0.097,-0.089,-0.081,-0.076,-0.076,-0.075,
     .-0.070,-0.066,-0.062,-0.068,-0.059,-0.047,-0.035,-0.024,
     .-0.019,-0.028,-0.032,-0.041,-0.050,-0.058,-0.062,-0.066,
     .-0.068,-0.069,-0.070,-0.072,-0.071,-0.070,-0.067,-0.064,
     .-0.067,-0.061,-0.052,-0.039,-0.026,-0.015,-0.012,-0.006,
     . 0.002, 0.015, 0.025, 0.034, 0.044, 0.051, 0.061, 0.068,
     . 0.072, 0.072, 0.072, 0.076, 0.081, 0.080, 0.088, 0.098,
     . 0.112,
c XARR; ROW j = 39
     .-0.121,-0.115,-0.110,-0.107,-0.105,-0.102,-0.101,-0.096,
     .-0.099,-0.100,-0.093,-0.087,-0.079,-0.074,-0.075,-0.075,
     .-0.068,-0.064,-0.058,-0.062,-0.052,-0.044,-0.031,-0.022,
     .-0.018,-0.029,-0.034,-0.044,-0.052,-0.059,-0.063,-0.066,
     .-0.066,-0.065,-0.064,-0.063,-0.064,-0.063,-0.060,-0.056,
     .-0.060,-0.055,-0.048,-0.034,-0.020,-0.010,-0.009,-0.002,
     . 0.004, 0.017, 0.027, 0.036, 0.046, 0.053, 0.061, 0.069,
     . 0.071, 0.071, 0.073, 0.077, 0.082, 0.082, 0.089, 0.098,
     . 0.110,
c XARR; ROW j = 40
     .-0.116,-0.110,-0.104,-0.100,-0.097,-0.094,-0.095,-0.091,
     .-0.094,-0.095,-0.087,-0.081,-0.075,-0.070,-0.072,-0.070,
     .-0.064,-0.059,-0.052,-0.055,-0.047,-0.037,-0.025,-0.017,
     .-0.015,-0.026,-0.033,-0.042,-0.049,-0.055,-0.057,-0.059,
     .-0.061,-0.057,-0.054,-0.054,-0.053,-0.055,-0.053,-0.048,
     .-0.053,-0.050,-0.042,-0.030,-0.016,-0.006,-0.005, 0.004,
     . 0.011, 0.021, 0.031, 0.039, 0.049, 0.057, 0.062, 0.068,
     . 0.071, 0.071, 0.074, 0.075, 0.081, 0.080, 0.086, 0.096,
     . 0.109,
c XARR; ROW j = 41
     .-0.113,-0.105,-0.098,-0.093,-0.091,-0.089,-0.087,-0.086,
     .-0.088,-0.087,-0.082,-0.072,-0.066,-0.062,-0.063,-0.062,
     .-0.057,-0.050,-0.045,-0.048,-0.039,-0.029,-0.017,-0.011,
     .-0.010,-0.023,-0.029,-0.036,-0.043,-0.050,-0.050,-0.051,
     .-0.053,-0.048,-0.047,-0.045,-0.044,-0.043,-0.045,-0.039,
     .-0.047,-0.043,-0.036,-0.024,-0.012, 0.000, 0.001, 0.011,
     . 0.020, 0.028, 0.037, 0.043, 0.052, 0.060, 0.063, 0.069,
     . 0.072, 0.073, 0.073, 0.074, 0.078, 0.075, 0.082, 0.092,
     . 0.105,
c XARR; ROW j = 42
     .-0.112,-0.100,-0.090,-0.084,-0.082,-0.081,-0.080,-0.080,
     .-0.081,-0.081,-0.074,-0.063,-0.055,-0.053,-0.055,-0.053,
     .-0.047,-0.041,-0.036,-0.042,-0.035,-0.023,-0.014,-0.006,
     .-0.006,-0.019,-0.023,-0.029,-0.036,-0.039,-0.042,-0.042,
     .-0.040,-0.040,-0.039,-0.037,-0.037,-0.035,-0.035,-0.036,
     .-0.040,-0.036,-0.030,-0.016,-0.004, 0.007, 0.012, 0.019,
     . 0.028, 0.036, 0.043, 0.049, 0.056, 0.060, 0.067, 0.070,
     . 0.072, 0.073, 0.074, 0.075, 0.077, 0.073, 0.078, 0.086,
     . 0.099,
c XARR; ROW j = 43
     .-0.109,-0.095,-0.084,-0.077,-0.076,-0.074,-0.074,-0.074,
     .-0.076,-0.076,-0.067,-0.057,-0.046,-0.044,-0.045,-0.044,
     .-0.036,-0.032,-0.029,-0.035,-0.029,-0.020,-0.010,-0.003,
     . 0.000,-0.010,-0.013,-0.019,-0.025,-0.027,-0.028,-0.029,
     .-0.030,-0.030,-0.030,-0.028,-0.028,-0.028,-0.028,-0.028,
     .-0.034,-0.030,-0.022,-0.008, 0.005, 0.015, 0.020, 0.028,
     . 0.036, 0.044, 0.047, 0.051, 0.058, 0.062, 0.068, 0.071,
     . 0.073, 0.072, 0.071, 0.072, 0.072, 0.067, 0.071, 0.080,
     . 0.092,
c XARR; ROW j = 44
     .-0.103,-0.090,-0.079,-0.072,-0.070,-0.067,-0.066,-0.067,
     .-0.071,-0.070,-0.063,-0.049,-0.040,-0.035,-0.034,-0.033,
     .-0.028,-0.023,-0.023,-0.027,-0.022,-0.016,-0.005, 0.002,
     . 0.005,-0.001,-0.003,-0.007,-0.010,-0.014,-0.017,-0.019,
     .-0.021,-0.022,-0.022,-0.023,-0.022,-0.023,-0.022,-0.023,
     .-0.028,-0.023,-0.014,-0.001, 0.013, 0.023, 0.027, 0.034,
     . 0.039, 0.047, 0.050, 0.054, 0.059, 0.063, 0.069, 0.072,
     . 0.075, 0.074, 0.072, 0.068, 0.068, 0.060, 0.065, 0.075,
     . 0.089,
c XARR; ROW j = 45
     .-0.101,-0.085,-0.073,-0.065,-0.061,-0.060,-0.062,-0.062,
     .-0.067,-0.065,-0.057,-0.046,-0.034,-0.026,-0.026,-0.023,
     .-0.019,-0.016,-0.016,-0.023,-0.019,-0.011, 0.000, 0.008,
     . 0.013, 0.009, 0.008, 0.005, 0.004,-0.001,-0.003,-0.009,
     .-0.014,-0.015,-0.016,-0.017,-0.016,-0.018,-0.017,-0.015,
     .-0.020,-0.014,-0.005, 0.008, 0.022, 0.031, 0.032, 0.040,
     . 0.043, 0.050, 0.052, 0.054, 0.061, 0.068, 0.071, 0.075,
     . 0.076, 0.075, 0.073, 0.069, 0.064, 0.056, 0.060, 0.069,
     . 0.084,
c XARR; ROW j = 46
     .-0.097,-0.082,-0.070,-0.061,-0.057,-0.054,-0.054,-0.057,
     .-0.062,-0.062,-0.053,-0.039,-0.027,-0.019,-0.020,-0.017,
     .-0.014,-0.012,-0.010,-0.018,-0.014,-0.003, 0.009, 0.018,
     . 0.025, 0.020, 0.021, 0.019, 0.016, 0.011, 0.007, 0.000,
     .-0.006,-0.009,-0.008,-0.009,-0.009,-0.009,-0.009,-0.006,
     .-0.010,-0.005, 0.005, 0.018, 0.030, 0.040, 0.039, 0.046,
     . 0.050, 0.054, 0.057, 0.060, 0.066, 0.073, 0.077, 0.079,
     . 0.081, 0.078, 0.075, 0.067, 0.063, 0.055, 0.057, 0.064,
     . 0.078,
c XARR; ROW j = 47
     .-0.095,-0.079,-0.067,-0.058,-0.054,-0.049,-0.049,-0.051,
     .-0.055,-0.054,-0.045,-0.031,-0.018,-0.011,-0.011,-0.012,
     .-0.009,-0.007,-0.003,-0.009,-0.004, 0.008, 0.020, 0.032,
     . 0.040, 0.035, 0.035, 0.031, 0.026, 0.021, 0.015, 0.010,
     . 0.006, 0.002, 0.002, 0.002, 0.002, 0.004, 0.004, 0.004,
     . 0.000, 0.004, 0.015, 0.028, 0.038, 0.047, 0.047, 0.051,
     . 0.057, 0.060, 0.063, 0.066, 0.074, 0.078, 0.083, 0.085,
     . 0.085, 0.080, 0.074, 0.066, 0.062, 0.051, 0.052, 0.059,
     . 0.072,
c XARR; ROW j = 48
     .-0.088,-0.074,-0.062,-0.053,-0.049,-0.044,-0.043,-0.042,
     .-0.046,-0.046,-0.034,-0.021,-0.008,-0.002,-0.005,-0.005,
     .-0.002, 0.001, 0.006, 0.001, 0.009, 0.023, 0.037, 0.049,
     . 0.055, 0.049, 0.048, 0.044, 0.036, 0.031, 0.026, 0.022,
     . 0.020, 0.018, 0.018, 0.018, 0.019, 0.019, 0.019, 0.019,
     . 0.014, 0.016, 0.025, 0.036, 0.047, 0.055, 0.055, 0.058,
     . 0.063, 0.067, 0.069, 0.073, 0.079, 0.082, 0.087, 0.089,
     . 0.087, 0.081, 0.073, 0.064, 0.058, 0.044, 0.045, 0.051,
     . 0.064,
c XARR; ROW j = 49
     .-0.077,-0.065,-0.055,-0.047,-0.042,-0.036,-0.033,-0.033,
     .-0.036,-0.033,-0.023,-0.009, 0.002, 0.006, 0.004, 0.004,
     . 0.008, 0.013, 0.018, 0.016, 0.027, 0.040, 0.055, 0.065,
     . 0.070, 0.063, 0.060, 0.052, 0.046, 0.040, 0.036, 0.035,
     . 0.032, 0.032, 0.033, 0.033, 0.032, 0.029, 0.026, 0.027,
     . 0.020, 0.022, 0.027, 0.037, 0.047, 0.054, 0.054, 0.058,
     . 0.063, 0.069, 0.069, 0.070, 0.076, 0.078, 0.082, 0.084,
     . 0.081, 0.071, 0.061, 0.048, 0.041, 0.028, 0.027, 0.032,
     . 0.043,
c XARR; ROW j = 50
     .-0.063,-0.052,-0.043,-0.035,-0.030,-0.024,-0.020,-0.017,
     .-0.021,-0.019,-0.009, 0.003, 0.012, 0.017, 0.015, 0.015,
     . 0.021, 0.027, 0.033, 0.035, 0.045, 0.057, 0.070, 0.078,
     . 0.082, 0.073, 0.065, 0.058, 0.053, 0.046, 0.043, 0.043,
     . 0.045, 0.048, 0.048, 0.045, 0.042, 0.038, 0.035, 0.034,
     . 0.026, 0.025, 0.028, 0.035, 0.044, 0.050, 0.049, 0.053,
     . 0.058, 0.064, 0.064, 0.064, 0.070, 0.073, 0.077, 0.077,
     . 0.073, 0.062, 0.052, 0.036, 0.027, 0.018, 0.014, 0.015,
     . 0.021,
c XARR; ROW j = 51
     .-0.047,-0.034,-0.024,-0.016,-0.012,-0.008,-0.007,-0.002,
     .-0.004,-0.003, 0.005, 0.013, 0.021, 0.026, 0.025, 0.027,
     . 0.034, 0.042, 0.049, 0.048, 0.056, 0.067, 0.078, 0.084,
     . 0.086, 0.075, 0.067, 0.060, 0.056, 0.049, 0.050, 0.050,
     . 0.053, 0.057, 0.059, 0.055, 0.051, 0.046, 0.040, 0.039,
     . 0.031, 0.030, 0.031, 0.035, 0.039, 0.046, 0.044, 0.048,
     . 0.052, 0.057, 0.058, 0.059, 0.065, 0.071, 0.074, 0.075,
     . 0.069, 0.058, 0.046, 0.031, 0.018, 0.008, 0.003, 0.002,
     . 0.005,
c XARR; ROW j = 52
     .-0.026,-0.017,-0.010,-0.003, 0.001, 0.006, 0.009, 0.010,
     . 0.008, 0.007, 0.012, 0.021, 0.028, 0.031, 0.031, 0.035,
     . 0.044, 0.053, 0.060, 0.056, 0.060, 0.070, 0.077, 0.081,
     . 0.083, 0.071, 0.067, 0.061, 0.056, 0.053, 0.053, 0.055,
     . 0.059, 0.064, 0.066, 0.063, 0.056, 0.051, 0.046, 0.045,
     . 0.036, 0.034, 0.035, 0.035, 0.038, 0.041, 0.039, 0.043,
     . 0.048, 0.051, 0.054, 0.056, 0.063, 0.068, 0.073, 0.075,
     . 0.069, 0.058, 0.045, 0.029, 0.015, 0.004,-0.002,-0.005,
     .-0.004,
c XARR; ROW j = 53
     .-0.019,-0.009,-0.001, 0.006, 0.011, 0.016, 0.017, 0.017,
     . 0.014, 0.011, 0.017, 0.023, 0.030, 0.036, 0.038, 0.041,
     . 0.051, 0.059, 0.066, 0.059, 0.060, 0.065, 0.068, 0.073,
     . 0.075, 0.066, 0.064, 0.060, 0.056, 0.053, 0.053, 0.058,
     . 0.063, 0.067, 0.070, 0.066, 0.060, 0.056, 0.050, 0.049,
     . 0.042, 0.039, 0.038, 0.038, 0.038, 0.039, 0.037, 0.040,
     . 0.043, 0.048, 0.050, 0.051, 0.060, 0.065, 0.072, 0.075,
     . 0.069, 0.054, 0.041, 0.025, 0.014,-0.002,-0.010,-0.015,
     .-0.017,
c XARR; ROW j = 54
     .-0.013,-0.003, 0.005, 0.011, 0.016, 0.018, 0.020, 0.020,
     . 0.013, 0.008, 0.015, 0.022, 0.031, 0.038, 0.041, 0.048,
     . 0.058, 0.066, 0.068, 0.057, 0.055, 0.058, 0.062, 0.067,
     . 0.070, 0.065, 0.065, 0.062, 0.059, 0.053, 0.053, 0.058,
     . 0.062, 0.067, 0.071, 0.067, 0.062, 0.056, 0.053, 0.053,
     . 0.047, 0.042, 0.040, 0.038, 0.037, 0.037, 0.034, 0.037,
     . 0.038, 0.044, 0.045, 0.045, 0.053, 0.060, 0.068, 0.071,
     . 0.063, 0.048, 0.035, 0.018, 0.005,-0.011,-0.021,-0.029,
     .-0.035,
c XARR; ROW j = 55
     .-0.012, 0.000, 0.009, 0.016, 0.018, 0.019, 0.017, 0.015,
     . 0.007, 0.004, 0.010, 0.020, 0.032, 0.040, 0.044, 0.054,
     . 0.063, 0.067, 0.067, 0.054, 0.051, 0.052, 0.056, 0.063,
     . 0.070, 0.066, 0.067, 0.063, 0.059, 0.051, 0.049, 0.053,
     . 0.056, 0.061, 0.065, 0.060, 0.055, 0.050, 0.046, 0.047,
     . 0.041, 0.037, 0.032, 0.031, 0.030, 0.031, 0.028, 0.028,
     . 0.029, 0.034, 0.034, 0.032, 0.040, 0.044, 0.052, 0.054,
     . 0.045, 0.028, 0.014,-0.002,-0.015,-0.031,-0.042,-0.053,
     .-0.062,
c XARR; ROW j = 56
     .-0.008, 0.004, 0.013, 0.018, 0.020, 0.020, 0.014, 0.011,
     . 0.006, 0.002, 0.010, 0.021, 0.035, 0.046, 0.052, 0.057,
     . 0.064, 0.065, 0.063, 0.049, 0.046, 0.047, 0.056, 0.065,
     . 0.074, 0.071, 0.069, 0.065, 0.059, 0.049, 0.045, 0.044,
     . 0.046, 0.053, 0.057, 0.053, 0.048, 0.043, 0.039, 0.041,
     . 0.034, 0.030, 0.026, 0.024, 0.023, 0.023, 0.018, 0.020,
     . 0.021, 0.023, 0.021, 0.018, 0.024, 0.031, 0.035, 0.034,
     . 0.024, 0.008,-0.006,-0.024,-0.040,-0.055,-0.068,-0.078,
     .-0.087,
c XARR; ROW j = 57
     .-0.004, 0.008, 0.017, 0.022, 0.023, 0.019, 0.014, 0.009,
     . 0.003, 0.002, 0.011, 0.025, 0.042, 0.054, 0.058, 0.059,
     . 0.063, 0.062, 0.058, 0.045, 0.043, 0.047, 0.059, 0.069,
     . 0.078, 0.073, 0.071, 0.065, 0.056, 0.042, 0.035, 0.035,
     . 0.038, 0.046, 0.049, 0.045, 0.040, 0.037, 0.033, 0.032,
     . 0.025, 0.021, 0.018, 0.016, 0.015, 0.013, 0.009, 0.010,
     . 0.011, 0.010, 0.008, 0.004, 0.008, 0.013, 0.015, 0.012,
     . 0.004,-0.012,-0.029,-0.048,-0.065,-0.080,-0.093,-0.104,
     .-0.113,
c XARR; ROW j = 58
     . 0.000, 0.012, 0.020, 0.025, 0.023, 0.019, 0.015, 0.008,
     . 0.005, 0.004, 0.016, 0.033, 0.048, 0.058, 0.059, 0.057,
     . 0.058, 0.056, 0.051, 0.039, 0.039, 0.048, 0.061, 0.072,
     . 0.079, 0.073, 0.069, 0.059, 0.047, 0.032, 0.024, 0.024,
     . 0.031, 0.038, 0.043, 0.037, 0.032, 0.028, 0.026, 0.022,
     . 0.017, 0.014, 0.009, 0.008, 0.008, 0.005, 0.003, 0.003,
     . 0.000, 0.000,-0.006,-0.012,-0.008,-0.007,-0.005,-0.011,
     .-0.023,-0.038,-0.054,-0.072,-0.090,-0.106,-0.120,-0.133,
     .-0.144,
c XARR; ROW j = 59
     . 0.004, 0.014, 0.020, 0.024, 0.021, 0.018, 0.014, 0.011,
     . 0.008, 0.007, 0.023, 0.038, 0.053, 0.059, 0.057, 0.051,
     . 0.052, 0.048, 0.046, 0.036, 0.040, 0.050, 0.064, 0.072,
     . 0.079, 0.070, 0.062, 0.049, 0.036, 0.020, 0.014, 0.016,
     . 0.022, 0.031, 0.035, 0.031, 0.025, 0.019, 0.016, 0.014,
     . 0.007, 0.003,-0.001,-0.002,-0.002,-0.004,-0.009,-0.011,
     .-0.014,-0.016,-0.023,-0.031,-0.030,-0.030,-0.030,-0.037,
     .-0.048,-0.067,-0.081,-0.101,-0.117,-0.133,-0.148,-0.162,
     .-0.176,
c XARR; ROW j = 60
     . 0.010, 0.016, 0.020, 0.022, 0.018, 0.015, 0.016, 0.013,
     . 0.012, 0.016, 0.028, 0.043, 0.056, 0.057, 0.052, 0.047,
     . 0.045, 0.044, 0.044, 0.036, 0.044, 0.052, 0.064, 0.071,
     . 0.072, 0.059, 0.049, 0.035, 0.022, 0.009, 0.004, 0.006,
     . 0.013, 0.021, 0.025, 0.020, 0.014, 0.008, 0.005, 0.002,
     .-0.005,-0.010,-0.013,-0.014,-0.013,-0.017,-0.022,-0.026,
     .-0.032,-0.035,-0.043,-0.051,-0.053,-0.056,-0.058,-0.063,
     .-0.076,-0.092,-0.111,-0.130,-0.146,-0.165,-0.178,-0.190,
     .-0.199,
c XARR; ROW j = 61
     . 0.016, 0.018, 0.020, 0.020, 0.016, 0.016, 0.017, 0.018,
     . 0.017, 0.022, 0.033, 0.044, 0.053, 0.053, 0.045, 0.039,
     . 0.039, 0.039, 0.042, 0.038, 0.044, 0.052, 0.061, 0.064,
     . 0.061, 0.045, 0.034, 0.021, 0.010,-0.002,-0.007,-0.005,
     . 0.003, 0.011, 0.015, 0.009, 0.001,-0.004,-0.009,-0.011,
     .-0.019,-0.021,-0.025,-0.026,-0.027,-0.031,-0.039,-0.044,
     .-0.050,-0.056,-0.065,-0.076,-0.079,-0.082,-0.083,-0.092,
     .-0.104,-0.124,-0.140,-0.158,-0.174,-0.195,-0.209,-0.220,
     .-0.228,
c XARR; ROW j = 62
     . 0.020, 0.017, 0.015, 0.014, 0.012, 0.015, 0.018, 0.021,
     . 0.023, 0.024, 0.035, 0.045, 0.050, 0.047, 0.038, 0.030,
     . 0.033, 0.036, 0.041, 0.037, 0.043, 0.048, 0.054, 0.055,
     . 0.049, 0.032, 0.020, 0.007,-0.002,-0.012,-0.016,-0.014,
     .-0.009,-0.001, 0.004,-0.003,-0.012,-0.017,-0.023,-0.025,
     .-0.031,-0.036,-0.039,-0.040,-0.042,-0.049,-0.057,-0.063,
     .-0.070,-0.077,-0.087,-0.102,-0.104,-0.106,-0.111,-0.119,
     .-0.133,-0.152,-0.168,-0.188,-0.204,-0.222,-0.237,-0.249,
     .-0.260,
c XARR; ROW j = 63
     . 0.026, 0.018, 0.013, 0.011, 0.010, 0.017, 0.020, 0.025,
     . 0.025, 0.025, 0.033, 0.040, 0.044, 0.042, 0.032, 0.023,
     . 0.028, 0.035, 0.042, 0.037, 0.041, 0.044, 0.046, 0.044,
     . 0.036, 0.017, 0.008,-0.003,-0.011,-0.021,-0.026,-0.024,
     .-0.017,-0.012,-0.010,-0.018,-0.026,-0.029,-0.033,-0.037,
     .-0.045,-0.048,-0.052,-0.053,-0.059,-0.065,-0.075,-0.081,
     .-0.091,-0.100,-0.115,-0.126,-0.130,-0.134,-0.139,-0.147,
     .-0.163,-0.182,-0.201,-0.223,-0.238,-0.252,-0.269,-0.284,
     .-0.298,
c XARR; ROW j = 64
     . 0.030, 0.022, 0.017, 0.014, 0.012, 0.018, 0.022, 0.026,
     . 0.024, 0.022, 0.030, 0.036, 0.040, 0.038, 0.029, 0.021,
     . 0.028, 0.034, 0.041, 0.035, 0.038, 0.039, 0.039, 0.035,
     . 0.027, 0.009,-0.001,-0.011,-0.020,-0.030,-0.036,-0.034,
     .-0.028,-0.023,-0.022,-0.030,-0.039,-0.042,-0.046,-0.051,
     .-0.059,-0.064,-0.070,-0.073,-0.079,-0.088,-0.099,-0.108,
     .-0.119,-0.130,-0.146,-0.158,-0.162,-0.167,-0.172,-0.182,
     .-0.198,-0.219,-0.239,-0.262,-0.279,-0.294,-0.311,-0.327,
     .-0.342,
c XARR; ROW j = 65
     . 0.037, 0.028, 0.022, 0.019, 0.015, 0.019, 0.026, 0.027,
     . 0.022, 0.018, 0.024, 0.032, 0.036, 0.034, 0.027, 0.022,
     . 0.029, 0.034, 0.039, 0.032, 0.033, 0.032, 0.032, 0.027,
     . 0.020, 0.002,-0.009,-0.018,-0.027,-0.038,-0.045,-0.045,
     .-0.040,-0.035,-0.034,-0.042,-0.050,-0.055,-0.060,-0.065,
     .-0.075,-0.082,-0.090,-0.097,-0.102,-0.114,-0.127,-0.138,
     .-0.150,-0.164,-0.180,-0.194,-0.199,-0.205,-0.210,-0.222,
     .-0.237,-0.260,-0.283,-0.306,-0.325,-0.343,-0.361,-0.377,
     .-0.392,
c XARR; ROW j = 66
     . 0.045, 0.036, 0.030, 0.027, 0.019, 0.021, 0.032, 0.027,
     . 0.018, 0.011, 0.017, 0.026, 0.032, 0.030, 0.027, 0.025,
     . 0.032, 0.035, 0.037, 0.029, 0.028, 0.025, 0.025, 0.018,
     . 0.013,-0.003,-0.014,-0.023,-0.033,-0.047,-0.053,-0.055,
     .-0.052,-0.047,-0.047,-0.053,-0.061,-0.069,-0.073,-0.080,
     .-0.091,-0.102,-0.113,-0.124,-0.129,-0.144,-0.158,-0.173,
     .-0.186,-0.203,-0.218,-0.234,-0.240,-0.246,-0.251,-0.265,
     .-0.280,-0.306,-0.331,-0.354,-0.377,-0.400,-0.417,-0.434,
     .-0.449/
      end
 
      block data wfc3_yarr
      real yarr(65,66)
      common /wfc3yarr_/yarr
      data yarr /
c YARR; ROW j = 01
     .-0.285,-0.267,-0.252,-0.238,-0.223,-0.217,-0.209,-0.204,
     .-0.202,-0.197,-0.188,-0.175,-0.158,-0.145,-0.139,-0.132,
     .-0.124,-0.115,-0.105,-0.102,-0.090,-0.076,-0.060,-0.048,
     .-0.041,-0.040,-0.036,-0.035,-0.031,-0.027,-0.018,-0.008,
     . 0.000, 0.008, 0.012, 0.016, 0.016, 0.015, 0.010, 0.007,
     . 0.000,-0.005,-0.009,-0.011,-0.012,-0.016,-0.024,-0.029,
     .-0.035,-0.041,-0.049,-0.056,-0.059,-0.060,-0.059,-0.059,
     .-0.061,-0.064,-0.067,-0.070,-0.076,-0.080,-0.079,-0.076,
     .-0.072,
c YARR; ROW j = 02
     .-0.288,-0.272,-0.257,-0.244,-0.233,-0.225,-0.217,-0.212,
     .-0.209,-0.203,-0.192,-0.177,-0.160,-0.147,-0.139,-0.131,
     .-0.124,-0.116,-0.106,-0.102,-0.090,-0.077,-0.063,-0.050,
     .-0.041,-0.041,-0.037,-0.035,-0.030,-0.026,-0.016,-0.007,
     . 0.003, 0.011, 0.015, 0.018, 0.018, 0.017, 0.015, 0.011,
     . 0.004, 0.001,-0.001,-0.002,-0.003,-0.007,-0.015,-0.020,
     .-0.026,-0.033,-0.041,-0.048,-0.052,-0.051,-0.050,-0.050,
     .-0.052,-0.056,-0.059,-0.062,-0.068,-0.071,-0.071,-0.068,
     .-0.064,
c YARR; ROW j = 03
     .-0.286,-0.271,-0.257,-0.244,-0.236,-0.227,-0.219,-0.215,
     .-0.210,-0.203,-0.190,-0.174,-0.157,-0.143,-0.134,-0.126,
     .-0.120,-0.113,-0.103,-0.098,-0.087,-0.075,-0.061,-0.048,
     .-0.037,-0.037,-0.033,-0.031,-0.026,-0.020,-0.011,-0.002,
     . 0.008, 0.016, 0.020, 0.022, 0.022, 0.022, 0.021, 0.018,
     . 0.010, 0.009, 0.008, 0.008, 0.006, 0.002,-0.006,-0.011,
     .-0.017,-0.024,-0.032,-0.041,-0.044,-0.043,-0.041,-0.041,
     .-0.043,-0.048,-0.051,-0.055,-0.059,-0.062,-0.062,-0.060,
     .-0.056,
c YARR; ROW j = 04
     .-0.279,-0.264,-0.251,-0.240,-0.232,-0.223,-0.216,-0.211,
     .-0.205,-0.197,-0.184,-0.167,-0.149,-0.134,-0.125,-0.117,
     .-0.111,-0.104,-0.095,-0.091,-0.079,-0.068,-0.054,-0.041,
     .-0.030,-0.030,-0.025,-0.022,-0.016,-0.011,-0.001, 0.007,
     . 0.017, 0.025, 0.028, 0.029, 0.029, 0.029, 0.029, 0.026,
     . 0.018, 0.018, 0.018, 0.018, 0.017, 0.012, 0.003,-0.001,
     .-0.008,-0.016,-0.024,-0.033,-0.036,-0.034,-0.032,-0.032,
     .-0.034,-0.039,-0.043,-0.046,-0.049,-0.053,-0.054,-0.052,
     .-0.048,
c YARR; ROW j = 05
     .-0.263,-0.251,-0.240,-0.230,-0.221,-0.214,-0.208,-0.199,
     .-0.194,-0.184,-0.171,-0.154,-0.136,-0.120,-0.112,-0.105,
     .-0.098,-0.089,-0.082,-0.079,-0.068,-0.056,-0.041,-0.027,
     .-0.018,-0.017,-0.010,-0.005, 0.002, 0.008, 0.016, 0.025,
     . 0.033, 0.039, 0.041, 0.042, 0.040, 0.041, 0.041, 0.039,
     . 0.033, 0.032, 0.030, 0.032, 0.032, 0.027, 0.016, 0.009,
     . 0.002,-0.005,-0.014,-0.021,-0.025,-0.025,-0.022,-0.022,
     .-0.023,-0.028,-0.032,-0.034,-0.037,-0.043,-0.042,-0.039,
     .-0.033,
c YARR; ROW j = 06
     .-0.247,-0.235,-0.224,-0.215,-0.206,-0.199,-0.193,-0.183,
     .-0.178,-0.168,-0.152,-0.134,-0.115,-0.100,-0.093,-0.086,
     .-0.078,-0.071,-0.063,-0.063,-0.052,-0.039,-0.025,-0.011,
     .-0.001,-0.001, 0.006, 0.011, 0.019, 0.027, 0.034, 0.041,
     . 0.047, 0.051, 0.052, 0.052, 0.052, 0.052, 0.052, 0.049,
     . 0.044, 0.044, 0.044, 0.043, 0.041, 0.035, 0.026, 0.018,
     . 0.008, 0.001,-0.007,-0.015,-0.019,-0.019,-0.018,-0.016,
     .-0.019,-0.023,-0.025,-0.027,-0.030,-0.034,-0.033,-0.030,
     .-0.024,
c YARR; ROW j = 07
     .-0.229,-0.217,-0.206,-0.196,-0.189,-0.182,-0.174,-0.167,
     .-0.160,-0.150,-0.134,-0.116,-0.098,-0.082,-0.073,-0.068,
     .-0.061,-0.053,-0.046,-0.044,-0.033,-0.021,-0.006, 0.005,
     . 0.014, 0.016, 0.021, 0.026, 0.032, 0.038, 0.046, 0.051,
     . 0.055, 0.057, 0.058, 0.057, 0.058, 0.058, 0.059, 0.057,
     . 0.053, 0.052, 0.052, 0.051, 0.048, 0.042, 0.033, 0.024,
     . 0.014, 0.007,-0.002,-0.009,-0.010,-0.012,-0.012,-0.010,
     .-0.011,-0.015,-0.017,-0.020,-0.022,-0.025,-0.023,-0.019,
     .-0.013,
c YARR; ROW j = 08
     .-0.211,-0.201,-0.191,-0.182,-0.174,-0.169,-0.159,-0.152,
     .-0.146,-0.135,-0.118,-0.100,-0.081,-0.067,-0.059,-0.051,
     .-0.044,-0.037,-0.029,-0.027,-0.018,-0.005, 0.007, 0.018,
     . 0.027, 0.027, 0.033, 0.038, 0.042, 0.046, 0.052, 0.058,
     . 0.061, 0.063, 0.061, 0.061, 0.062, 0.063, 0.063, 0.063,
     . 0.059, 0.058, 0.058, 0.056, 0.052, 0.049, 0.038, 0.031,
     . 0.023, 0.014, 0.005,-0.001,-0.004,-0.004,-0.003,-0.002,
     .-0.005,-0.008,-0.011,-0.013,-0.014,-0.015,-0.014,-0.010,
     .-0.006,
c YARR; ROW j = 09
     .-0.201,-0.189,-0.179,-0.170,-0.165,-0.156,-0.149,-0.140,
     .-0.133,-0.123,-0.106,-0.087,-0.068,-0.053,-0.043,-0.037,
     .-0.030,-0.022,-0.014,-0.012,-0.001, 0.008, 0.021, 0.033,
     . 0.041, 0.041, 0.045, 0.048, 0.051, 0.056, 0.061, 0.065,
     . 0.067, 0.069, 0.067, 0.068, 0.067, 0.067, 0.066, 0.066,
     . 0.062, 0.061, 0.059, 0.059, 0.056, 0.053, 0.044, 0.036,
     . 0.029, 0.022, 0.014, 0.006, 0.004, 0.002, 0.006, 0.006,
     . 0.002,-0.003,-0.004,-0.007,-0.009,-0.010,-0.008,-0.004,
     . 0.002,
c YARR; ROW j = 10
     .-0.186,-0.177,-0.168,-0.160,-0.152,-0.145,-0.137,-0.128,
     .-0.121,-0.108,-0.092,-0.073,-0.055,-0.038,-0.028,-0.020,
     .-0.013,-0.006, 0.000, 0.003, 0.013, 0.023, 0.035, 0.045,
     . 0.053, 0.052, 0.056, 0.058, 0.061, 0.062, 0.066, 0.068,
     . 0.072, 0.075, 0.073, 0.071, 0.071, 0.070, 0.070, 0.068,
     . 0.062, 0.060, 0.061, 0.060, 0.058, 0.055, 0.047, 0.042,
     . 0.038, 0.030, 0.021, 0.015, 0.012, 0.012, 0.013, 0.012,
     . 0.008, 0.002,-0.002,-0.005,-0.009,-0.007,-0.005,-0.002,
     . 0.004,
c YARR; ROW j = 11
     .-0.176,-0.168,-0.159,-0.151,-0.141,-0.132,-0.125,-0.114,
     .-0.108,-0.096,-0.078,-0.057,-0.039,-0.023,-0.014,-0.007,
     . 0.001, 0.008, 0.014, 0.013, 0.023, 0.032, 0.045, 0.054,
     . 0.061, 0.059, 0.063, 0.065, 0.067, 0.068, 0.070, 0.073,
     . 0.077, 0.080, 0.078, 0.077, 0.076, 0.074, 0.073, 0.070,
     . 0.064, 0.061, 0.061, 0.061, 0.059, 0.056, 0.049, 0.044,
     . 0.041, 0.035, 0.026, 0.020, 0.017, 0.016, 0.019, 0.018,
     . 0.012, 0.003,-0.003,-0.006,-0.008,-0.009,-0.008,-0.004,
     . 0.001,
c YARR; ROW j = 12
     .-0.165,-0.157,-0.150,-0.141,-0.131,-0.123,-0.113,-0.102,
     .-0.094,-0.083,-0.066,-0.046,-0.027,-0.012,-0.003, 0.004,
     . 0.012, 0.019, 0.023, 0.022, 0.031, 0.040, 0.051, 0.062,
     . 0.068, 0.066, 0.068, 0.069, 0.070, 0.070, 0.072, 0.074,
     . 0.078, 0.082, 0.082, 0.079, 0.077, 0.076, 0.074, 0.072,
     . 0.065, 0.062, 0.060, 0.061, 0.057, 0.056, 0.049, 0.044,
     . 0.040, 0.037, 0.029, 0.023, 0.020, 0.019, 0.021, 0.021,
     . 0.015, 0.006,-0.002,-0.007,-0.010,-0.009,-0.008,-0.006,
     .-0.003,
c YARR; ROW j = 13
     .-0.153,-0.146,-0.139,-0.131,-0.121,-0.112,-0.102,-0.093,
     .-0.084,-0.072,-0.055,-0.036,-0.020,-0.005, 0.005, 0.012,
     . 0.019, 0.025, 0.029, 0.031, 0.038, 0.047, 0.058, 0.067,
     . 0.072, 0.072, 0.073, 0.072, 0.072, 0.070, 0.070, 0.073,
     . 0.076, 0.081, 0.083, 0.079, 0.077, 0.075, 0.074, 0.073,
     . 0.068, 0.064, 0.063, 0.063, 0.061, 0.059, 0.052, 0.046,
     . 0.042, 0.039, 0.031, 0.026, 0.024, 0.023, 0.024, 0.024,
     . 0.019, 0.010, 0.003,-0.004,-0.007,-0.009,-0.008,-0.004,
     . 0.001,
c YARR; ROW j = 14
     .-0.148,-0.140,-0.131,-0.122,-0.114,-0.104,-0.096,-0.086,
     .-0.077,-0.066,-0.051,-0.033,-0.018,-0.006, 0.004, 0.011,
     . 0.018, 0.024, 0.029, 0.030, 0.039, 0.045, 0.053, 0.061,
     . 0.066, 0.064, 0.065, 0.065, 0.064, 0.061, 0.061, 0.064,
     . 0.070, 0.077, 0.077, 0.075, 0.074, 0.073, 0.072, 0.071,
     . 0.065, 0.064, 0.065, 0.063, 0.060, 0.057, 0.048, 0.044,
     . 0.040, 0.034, 0.027, 0.022, 0.019, 0.018, 0.018, 0.016,
     . 0.011, 0.003,-0.005,-0.011,-0.015,-0.018,-0.017,-0.015,
     .-0.011,
c YARR; ROW j = 15
     .-0.147,-0.136,-0.126,-0.116,-0.109,-0.098,-0.091,-0.082,
     .-0.073,-0.064,-0.050,-0.034,-0.020,-0.009, 0.001, 0.007,
     . 0.015, 0.022, 0.027, 0.028, 0.036, 0.041, 0.047, 0.051,
     . 0.053, 0.052, 0.053, 0.053, 0.053, 0.051, 0.050, 0.055,
     . 0.063, 0.070, 0.073, 0.071, 0.067, 0.066, 0.064, 0.063,
     . 0.059, 0.057, 0.057, 0.057, 0.054, 0.050, 0.041, 0.034,
     . 0.030, 0.025, 0.018, 0.011, 0.006, 0.002, 0.002, 0.002,
     .-0.003,-0.011,-0.019,-0.026,-0.031,-0.034,-0.032,-0.028,
     .-0.022,
c YARR; ROW j = 16
     .-0.146,-0.135,-0.124,-0.113,-0.103,-0.094,-0.087,-0.077,
     .-0.071,-0.062,-0.047,-0.032,-0.020,-0.010,-0.003, 0.004,
     . 0.014, 0.020, 0.024, 0.026, 0.030, 0.034, 0.037, 0.039,
     . 0.039, 0.038, 0.039, 0.039, 0.039, 0.037, 0.038, 0.043,
     . 0.053, 0.063, 0.067, 0.064, 0.062, 0.058, 0.057, 0.055,
     . 0.050, 0.049, 0.050, 0.048, 0.045, 0.039, 0.031, 0.025,
     . 0.021, 0.015, 0.008, 0.001,-0.005,-0.009,-0.013,-0.015,
     .-0.021,-0.029,-0.036,-0.043,-0.048,-0.048,-0.048,-0.046,
     .-0.043,
c YARR; ROW j = 17
     .-0.151,-0.137,-0.123,-0.111,-0.100,-0.091,-0.084,-0.076,
     .-0.069,-0.060,-0.047,-0.034,-0.022,-0.012,-0.005, 0.000,
     . 0.008, 0.016, 0.022, 0.022, 0.025, 0.025, 0.027, 0.027,
     . 0.026, 0.025, 0.026, 0.026, 0.027, 0.025, 0.027, 0.035,
     . 0.045, 0.055, 0.061, 0.058, 0.055, 0.052, 0.050, 0.048,
     . 0.043, 0.041, 0.041, 0.040, 0.036, 0.030, 0.021, 0.017,
     . 0.013, 0.009, 0.003,-0.004,-0.012,-0.017,-0.021,-0.024,
     .-0.028,-0.036,-0.044,-0.051,-0.056,-0.057,-0.057,-0.055,
     .-0.051,
c YARR; ROW j = 18
     .-0.152,-0.136,-0.123,-0.111,-0.098,-0.090,-0.083,-0.075,
     .-0.066,-0.058,-0.046,-0.034,-0.023,-0.014,-0.007, 0.000,
     . 0.006, 0.014, 0.018, 0.021, 0.024, 0.024, 0.023, 0.021,
     . 0.021, 0.021, 0.022, 0.022, 0.023, 0.021, 0.022, 0.031,
     . 0.042, 0.052, 0.059, 0.058, 0.054, 0.052, 0.050, 0.048,
     . 0.043, 0.042, 0.039, 0.038, 0.034, 0.030, 0.022, 0.018,
     . 0.013, 0.011, 0.007, 0.002,-0.007,-0.013,-0.020,-0.021,
     .-0.023,-0.030,-0.038,-0.045,-0.050,-0.053,-0.052,-0.049,
     .-0.044,
c YARR; ROW j = 19
     .-0.157,-0.138,-0.122,-0.108,-0.098,-0.088,-0.080,-0.073,
     .-0.066,-0.055,-0.043,-0.030,-0.020,-0.012,-0.005, 0.003,
     . 0.009, 0.015, 0.021, 0.025, 0.029, 0.029, 0.028, 0.026,
     . 0.027, 0.028, 0.030, 0.029, 0.030, 0.028, 0.028, 0.037,
     . 0.046, 0.056, 0.063, 0.063, 0.059, 0.058, 0.055, 0.055,
     . 0.051, 0.049, 0.048, 0.046, 0.042, 0.039, 0.034, 0.030,
     . 0.026, 0.023, 0.018, 0.014, 0.007, 0.001,-0.004,-0.007,
     .-0.010,-0.014,-0.021,-0.028,-0.032,-0.038,-0.038,-0.036,
     .-0.031,
c YARR; ROW j = 20
     .-0.158,-0.137,-0.119,-0.104,-0.096,-0.085,-0.077,-0.070,
     .-0.061,-0.051,-0.038,-0.026,-0.016,-0.008, 0.000, 0.007,
     . 0.012, 0.018, 0.024, 0.029, 0.033, 0.034, 0.033, 0.034,
     . 0.035, 0.036, 0.036, 0.037, 0.038, 0.036, 0.037, 0.043,
     . 0.053, 0.063, 0.069, 0.069, 0.066, 0.065, 0.064, 0.063,
     . 0.058, 0.058, 0.058, 0.057, 0.054, 0.051, 0.045, 0.042,
     . 0.039, 0.035, 0.029, 0.025, 0.018, 0.013, 0.008, 0.004,
     . 0.001,-0.003,-0.010,-0.017,-0.021,-0.027,-0.029,-0.029,
     .-0.027,
c YARR; ROW j = 21
     .-0.151,-0.129,-0.112,-0.097,-0.089,-0.078,-0.071,-0.062,
     .-0.053,-0.042,-0.030,-0.019,-0.009, 0.001, 0.007, 0.014,
     . 0.019, 0.025, 0.029, 0.035, 0.039, 0.043, 0.043, 0.044,
     . 0.044, 0.046, 0.046, 0.047, 0.048, 0.046, 0.046, 0.053,
     . 0.063, 0.072, 0.077, 0.076, 0.074, 0.073, 0.071, 0.070,
     . 0.066, 0.065, 0.064, 0.063, 0.062, 0.060, 0.054, 0.050,
     . 0.045, 0.041, 0.037, 0.033, 0.026, 0.021, 0.015, 0.012,
     . 0.008, 0.004,-0.003,-0.010,-0.016,-0.023,-0.025,-0.025,
     .-0.023,
c YARR; ROW j = 22
     .-0.142,-0.121,-0.102,-0.087,-0.079,-0.069,-0.060,-0.051,
     .-0.042,-0.033,-0.019,-0.008, 0.002, 0.010, 0.017, 0.021,
     . 0.028, 0.032, 0.034, 0.039, 0.045, 0.050, 0.052, 0.052,
     . 0.052, 0.053, 0.054, 0.053, 0.054, 0.053, 0.053, 0.060,
     . 0.068, 0.076, 0.082, 0.080, 0.079, 0.077, 0.075, 0.072,
     . 0.067, 0.066, 0.066, 0.067, 0.066, 0.064, 0.059, 0.055,
     . 0.050, 0.047, 0.042, 0.038, 0.033, 0.027, 0.020, 0.015,
     . 0.012, 0.006,-0.002,-0.009,-0.015,-0.021,-0.025,-0.027,
     .-0.028,
c YARR; ROW j = 23
     .-0.131,-0.110,-0.092,-0.077,-0.068,-0.058,-0.048,-0.040,
     .-0.032,-0.023,-0.012,-0.001, 0.007, 0.016, 0.023, 0.027,
     . 0.032, 0.036, 0.039, 0.045, 0.050, 0.055, 0.057, 0.057,
     . 0.056, 0.056, 0.056, 0.055, 0.054, 0.052, 0.053, 0.061,
     . 0.069, 0.078, 0.082, 0.081, 0.077, 0.076, 0.071, 0.068,
     . 0.062, 0.062, 0.061, 0.064, 0.063, 0.064, 0.060, 0.056,
     . 0.051, 0.048, 0.044, 0.042, 0.038, 0.032, 0.023, 0.018,
     . 0.015, 0.009, 0.001,-0.008,-0.016,-0.021,-0.026,-0.029,
     .-0.030,
c YARR; ROW j = 24
     .-0.121,-0.100,-0.082,-0.068,-0.059,-0.048,-0.039,-0.032,
     .-0.023,-0.016,-0.006, 0.003, 0.011, 0.016, 0.023, 0.029,
     . 0.033, 0.037, 0.040, 0.045, 0.052, 0.056, 0.058, 0.057,
     . 0.057, 0.056, 0.055, 0.053, 0.052, 0.051, 0.052, 0.059,
     . 0.067, 0.075, 0.080, 0.080, 0.076, 0.073, 0.067, 0.063,
     . 0.058, 0.057, 0.058, 0.060, 0.060, 0.062, 0.060, 0.055,
     . 0.052, 0.049, 0.046, 0.045, 0.041, 0.034, 0.028, 0.022,
     . 0.016, 0.009, 0.001,-0.007,-0.015,-0.023,-0.028,-0.031,
     .-0.032,
c YARR; ROW j = 25
     .-0.116,-0.095,-0.077,-0.062,-0.054,-0.043,-0.034,-0.027,
     .-0.019,-0.011,-0.003, 0.005, 0.011, 0.016, 0.024, 0.030,
     . 0.034, 0.037, 0.041, 0.046, 0.052, 0.057, 0.057, 0.057,
     . 0.057, 0.056, 0.054, 0.052, 0.049, 0.047, 0.049, 0.056,
     . 0.064, 0.073, 0.078, 0.076, 0.072, 0.069, 0.064, 0.060,
     . 0.055, 0.053, 0.054, 0.056, 0.056, 0.058, 0.056, 0.052,
     . 0.051, 0.048, 0.045, 0.044, 0.041, 0.035, 0.029, 0.020,
     . 0.012, 0.006,-0.002,-0.011,-0.020,-0.030,-0.036,-0.040,
     .-0.042,
c YARR; ROW j = 26
     .-0.112,-0.092,-0.075,-0.060,-0.051,-0.041,-0.031,-0.023,
     .-0.014,-0.007, 0.001, 0.008, 0.013, 0.019, 0.027, 0.032,
     . 0.036, 0.040, 0.041, 0.045, 0.050, 0.053, 0.054, 0.054,
     . 0.053, 0.052, 0.050, 0.049, 0.047, 0.044, 0.046, 0.051,
     . 0.061, 0.070, 0.074, 0.073, 0.069, 0.066, 0.063, 0.061,
     . 0.056, 0.055, 0.055, 0.056, 0.055, 0.055, 0.052, 0.050,
     . 0.048, 0.045, 0.042, 0.041, 0.037, 0.032, 0.024, 0.016,
     . 0.007, 0.000,-0.008,-0.018,-0.029,-0.039,-0.045,-0.050,
     .-0.053,
c YARR; ROW j = 27
     .-0.111,-0.091,-0.075,-0.060,-0.053,-0.041,-0.031,-0.021,
     .-0.013,-0.005, 0.002, 0.011, 0.017, 0.022, 0.028, 0.032,
     . 0.038, 0.041, 0.042, 0.043, 0.047, 0.048, 0.049, 0.049,
     . 0.048, 0.047, 0.048, 0.046, 0.044, 0.042, 0.043, 0.052,
     . 0.061, 0.070, 0.074, 0.073, 0.069, 0.068, 0.065, 0.064,
     . 0.060, 0.059, 0.056, 0.057, 0.054, 0.053, 0.048, 0.045,
     . 0.042, 0.040, 0.038, 0.036, 0.031, 0.025, 0.018, 0.009,
     . 0.001,-0.007,-0.018,-0.027,-0.038,-0.049,-0.057,-0.062,
     .-0.067,
c YARR; ROW j = 28
     .-0.112,-0.095,-0.080,-0.066,-0.056,-0.045,-0.033,-0.023,
     .-0.014,-0.008, 0.003, 0.011, 0.018, 0.023, 0.028, 0.031,
     . 0.037, 0.041, 0.041, 0.041, 0.043, 0.043, 0.045, 0.043,
     . 0.042, 0.044, 0.044, 0.043, 0.042, 0.038, 0.040, 0.048,
     . 0.058, 0.067, 0.072, 0.070, 0.067, 0.066, 0.066, 0.066,
     . 0.063, 0.060, 0.058, 0.057, 0.052, 0.050, 0.043, 0.040,
     . 0.038, 0.035, 0.031, 0.029, 0.024, 0.019, 0.012, 0.001,
     .-0.007,-0.015,-0.025,-0.036,-0.048,-0.060,-0.068,-0.073,
     .-0.076,
c YARR; ROW j = 29
     .-0.118,-0.100,-0.084,-0.069,-0.061,-0.049,-0.037,-0.027,
     .-0.018,-0.009, 0.001, 0.011, 0.017, 0.022, 0.027, 0.031,
     . 0.037, 0.040, 0.040, 0.039, 0.040, 0.042, 0.042, 0.040,
     . 0.038, 0.037, 0.038, 0.039, 0.038, 0.035, 0.036, 0.043,
     . 0.055, 0.065, 0.069, 0.068, 0.064, 0.063, 0.064, 0.064,
     . 0.060, 0.059, 0.056, 0.054, 0.049, 0.045, 0.040, 0.037,
     . 0.034, 0.030, 0.025, 0.022, 0.018, 0.012, 0.005,-0.005,
     .-0.015,-0.024,-0.033,-0.045,-0.057,-0.069,-0.078,-0.084,
     .-0.087,
c YARR; ROW j = 30
     .-0.131,-0.108,-0.089,-0.074,-0.068,-0.055,-0.043,-0.031,
     .-0.021,-0.013,-0.004, 0.006, 0.014, 0.018, 0.024, 0.027,
     . 0.032, 0.036, 0.039, 0.039, 0.042, 0.042, 0.041, 0.039,
     . 0.039, 0.038, 0.038, 0.037, 0.035, 0.032, 0.035, 0.042,
     . 0.054, 0.064, 0.068, 0.065, 0.061, 0.060, 0.059, 0.058,
     . 0.054, 0.052, 0.050, 0.047, 0.043, 0.041, 0.036, 0.033,
     . 0.031, 0.026, 0.019, 0.017, 0.010, 0.004,-0.003,-0.015,
     .-0.023,-0.032,-0.044,-0.055,-0.068,-0.083,-0.091,-0.096,
     .-0.098,
c YARR; ROW j = 31
     .-0.142,-0.119,-0.099,-0.081,-0.073,-0.060,-0.048,-0.036,
     .-0.026,-0.019,-0.009, 0.001, 0.009, 0.013, 0.019, 0.023,
     . 0.028, 0.033, 0.036, 0.037, 0.040, 0.040, 0.039, 0.037,
     . 0.036, 0.034, 0.034, 0.032, 0.030, 0.027, 0.029, 0.037,
     . 0.049, 0.058, 0.062, 0.059, 0.055, 0.052, 0.050, 0.048,
     . 0.044, 0.041, 0.040, 0.038, 0.034, 0.033, 0.028, 0.026,
     . 0.023, 0.017, 0.010, 0.006,-0.001,-0.008,-0.016,-0.027,
     .-0.036,-0.045,-0.056,-0.067,-0.080,-0.095,-0.104,-0.110,
     .-0.112,
c YARR; ROW j = 32
     .-0.155,-0.130,-0.109,-0.091,-0.079,-0.065,-0.053,-0.042,
     .-0.033,-0.026,-0.016,-0.006, 0.002, 0.006, 0.013, 0.018,
     . 0.025, 0.030, 0.033, 0.034, 0.038, 0.038, 0.037, 0.034,
     . 0.032, 0.031, 0.029, 0.026, 0.023, 0.020, 0.021, 0.029,
     . 0.041, 0.050, 0.055, 0.051, 0.046, 0.042, 0.038, 0.035,
     . 0.031, 0.028, 0.026, 0.025, 0.023, 0.022, 0.019, 0.016,
     . 0.013, 0.007, 0.000,-0.007,-0.015,-0.023,-0.031,-0.042,
     .-0.052,-0.060,-0.070,-0.081,-0.094,-0.109,-0.118,-0.124,
     .-0.127,
c YARR; ROW j = 33
     .-0.169,-0.143,-0.121,-0.101,-0.084,-0.071,-0.059,-0.049,
     .-0.042,-0.035,-0.025,-0.015,-0.007,-0.002, 0.006, 0.012,
     . 0.021, 0.027, 0.029, 0.031, 0.036, 0.036, 0.035, 0.031,
     . 0.028, 0.026, 0.024, 0.019, 0.015, 0.011, 0.011, 0.019,
     . 0.032, 0.040, 0.045, 0.041, 0.034, 0.028, 0.023, 0.018,
     . 0.015, 0.011, 0.009, 0.009, 0.009, 0.010, 0.007, 0.004,
     . 0.001,-0.005,-0.013,-0.023,-0.031,-0.041,-0.050,-0.058,
     .-0.070,-0.078,-0.085,-0.097,-0.110,-0.124,-0.134,-0.140,
     .-0.144,
c YARR; ROW j = 34
     . 0.110, 0.137, 0.162, 0.183, 0.201, 0.214, 0.227, 0.238,
     . 0.249, 0.256, 0.261, 0.263, 0.263, 0.261, 0.261, 0.258,
     . 0.255, 0.249, 0.245, 0.240, 0.232, 0.218, 0.203, 0.189,
     . 0.178, 0.167, 0.152, 0.137, 0.127, 0.119, 0.115, 0.114,
     . 0.115, 0.116, 0.118, 0.117, 0.115, 0.110, 0.107, 0.105,
     . 0.108, 0.102, 0.092, 0.079, 0.068, 0.060, 0.058, 0.048,
     . 0.036, 0.025, 0.015, 0.009,-0.002,-0.013,-0.022,-0.031,
     .-0.038,-0.044,-0.052,-0.057,-0.060,-0.059,-0.053,-0.043,
     .-0.029,
c YARR; ROW j = 35
     . 0.101, 0.126, 0.149, 0.169, 0.184, 0.197, 0.207, 0.217,
     . 0.227, 0.234, 0.239, 0.240, 0.240, 0.240, 0.240, 0.238,
     . 0.234, 0.230, 0.225, 0.219, 0.209, 0.197, 0.183, 0.170,
     . 0.161, 0.150, 0.137, 0.125, 0.114, 0.105, 0.099, 0.096,
     . 0.096, 0.097, 0.097, 0.097, 0.096, 0.094, 0.095, 0.096,
     . 0.099, 0.094, 0.085, 0.073, 0.063, 0.054, 0.049, 0.040,
     . 0.028, 0.018, 0.008, 0.000,-0.012,-0.022,-0.031,-0.040,
     .-0.047,-0.053,-0.059,-0.063,-0.064,-0.065,-0.059,-0.048,
     .-0.034,
c YARR; ROW j = 36
     . 0.094, 0.117, 0.138, 0.155, 0.169, 0.182, 0.189, 0.198,
     . 0.208, 0.215, 0.218, 0.220, 0.220, 0.220, 0.221, 0.218,
     . 0.214, 0.210, 0.204, 0.198, 0.187, 0.176, 0.162, 0.151,
     . 0.142, 0.132, 0.121, 0.111, 0.099, 0.089, 0.083, 0.078,
     . 0.077, 0.078, 0.078, 0.078, 0.077, 0.078, 0.081, 0.084,
     . 0.088, 0.084, 0.074, 0.064, 0.053, 0.044, 0.038, 0.029,
     . 0.017, 0.008,-0.002,-0.011,-0.022,-0.032,-0.042,-0.050,
     .-0.057,-0.063,-0.067,-0.071,-0.070,-0.072,-0.065,-0.054,
     .-0.039,
c YARR; ROW j = 37
     . 0.089, 0.110, 0.128, 0.144, 0.156, 0.168, 0.173, 0.181,
     . 0.190, 0.196, 0.199, 0.201, 0.201, 0.202, 0.202, 0.199,
     . 0.195, 0.191, 0.184, 0.177, 0.166, 0.155, 0.143, 0.132,
     . 0.124, 0.114, 0.104, 0.094, 0.083, 0.072, 0.065, 0.060,
     . 0.058, 0.060, 0.060, 0.060, 0.060, 0.063, 0.067, 0.071,
     . 0.074, 0.070, 0.061, 0.050, 0.039, 0.030, 0.023, 0.014,
     . 0.003,-0.005,-0.014,-0.024,-0.034,-0.043,-0.053,-0.060,
     .-0.068,-0.073,-0.076,-0.079,-0.078,-0.079,-0.072,-0.060,
     .-0.045,
c YARR; ROW j = 38
     . 0.084, 0.103, 0.118, 0.132, 0.144, 0.152, 0.161, 0.169,
     . 0.176, 0.179, 0.183, 0.186, 0.188, 0.188, 0.189, 0.185,
     . 0.180, 0.173, 0.165, 0.158, 0.149, 0.136, 0.126, 0.116,
     . 0.107, 0.099, 0.087, 0.075, 0.066, 0.056, 0.047, 0.042,
     . 0.040, 0.043, 0.044, 0.046, 0.047, 0.050, 0.054, 0.057,
     . 0.061, 0.058, 0.047, 0.033, 0.021, 0.010, 0.006,-0.002,
     .-0.013,-0.023,-0.031,-0.038,-0.047,-0.057,-0.066,-0.074,
     .-0.080,-0.086,-0.090,-0.090,-0.089,-0.087,-0.080,-0.071,
     .-0.058,
c YARR; ROW j = 39
     . 0.083, 0.099, 0.113, 0.125, 0.135, 0.143, 0.150, 0.157,
     . 0.161, 0.164, 0.165, 0.167, 0.169, 0.171, 0.171, 0.166,
     . 0.160, 0.152, 0.145, 0.138, 0.127, 0.117, 0.106, 0.097,
     . 0.088, 0.078, 0.067, 0.054, 0.045, 0.035, 0.028, 0.025,
     . 0.024, 0.025, 0.027, 0.029, 0.031, 0.034, 0.036, 0.038,
     . 0.042, 0.036, 0.024, 0.011,-0.002,-0.012,-0.015,-0.025,
     .-0.033,-0.041,-0.048,-0.054,-0.061,-0.072,-0.080,-0.087,
     .-0.091,-0.096,-0.098,-0.100,-0.098,-0.094,-0.087,-0.077,
     .-0.064,
c YARR; ROW j = 40
     . 0.083, 0.097, 0.110, 0.120, 0.127, 0.133, 0.139, 0.142,
     . 0.146, 0.146, 0.145, 0.146, 0.147, 0.147, 0.149, 0.145,
     . 0.138, 0.130, 0.121, 0.113, 0.103, 0.093, 0.080, 0.070,
     . 0.062, 0.052, 0.040, 0.028, 0.019, 0.013, 0.009, 0.007,
     . 0.006, 0.008, 0.011, 0.013, 0.015, 0.016, 0.018, 0.018,
     . 0.022, 0.015, 0.001,-0.013,-0.026,-0.034,-0.036,-0.046,
     .-0.054,-0.060,-0.065,-0.068,-0.076,-0.085,-0.092,-0.097,
     .-0.100,-0.103,-0.105,-0.105,-0.101,-0.098,-0.090,-0.079,
     .-0.065,
c YARR; ROW j = 41
     . 0.084, 0.095, 0.105, 0.112, 0.118, 0.123, 0.124, 0.127,
     . 0.128, 0.127, 0.123, 0.122, 0.122, 0.123, 0.124, 0.121,
     . 0.113, 0.105, 0.097, 0.090, 0.080, 0.067, 0.054, 0.043,
     . 0.036, 0.028, 0.016, 0.006,-0.002,-0.006,-0.006,-0.007,
     .-0.006,-0.004,-0.002,-0.001,-0.001,-0.001,-0.001,-0.001,
     .-0.001,-0.009,-0.022,-0.036,-0.048,-0.056,-0.058,-0.065,
     .-0.071,-0.076,-0.080,-0.084,-0.092,-0.100,-0.105,-0.109,
     .-0.110,-0.110,-0.110,-0.109,-0.104,-0.101,-0.092,-0.080,
     .-0.065,
c YARR; ROW j = 42
     . 0.085, 0.093, 0.099, 0.104, 0.109, 0.112, 0.112, 0.111,
     . 0.112, 0.109, 0.103, 0.100, 0.097, 0.098, 0.101, 0.096,
     . 0.088, 0.081, 0.073, 0.066, 0.056, 0.042, 0.029, 0.019,
     . 0.012, 0.004,-0.006,-0.014,-0.019,-0.021,-0.021,-0.022,
     .-0.020,-0.017,-0.017,-0.018,-0.019,-0.022,-0.022,-0.024,
     .-0.024,-0.033,-0.046,-0.058,-0.068,-0.075,-0.076,-0.081,
     .-0.086,-0.090,-0.093,-0.097,-0.105,-0.113,-0.118,-0.120,
     .-0.118,-0.117,-0.116,-0.114,-0.108,-0.104,-0.093,-0.079,
     .-0.062,
c YARR; ROW j = 43
     . 0.085, 0.090, 0.093, 0.096, 0.098, 0.099, 0.097, 0.098,
     . 0.096, 0.091, 0.085, 0.081, 0.077, 0.077, 0.078, 0.072,
     . 0.065, 0.057, 0.049, 0.041, 0.030, 0.016, 0.004,-0.005,
     .-0.011,-0.018,-0.026,-0.031,-0.035,-0.037,-0.037,-0.037,
     .-0.035,-0.032,-0.033,-0.035,-0.039,-0.041,-0.043,-0.045,
     .-0.046,-0.054,-0.066,-0.076,-0.086,-0.091,-0.092,-0.096,
     .-0.102,-0.105,-0.108,-0.111,-0.120,-0.127,-0.130,-0.130,
     .-0.127,-0.125,-0.122,-0.117,-0.111,-0.106,-0.093,-0.077,
     .-0.057,
c YARR; ROW j = 44
     . 0.086, 0.086, 0.086, 0.086, 0.086, 0.085, 0.084, 0.083,
     . 0.081, 0.074, 0.067, 0.062, 0.058, 0.056, 0.056, 0.049,
     . 0.041, 0.034, 0.024, 0.016, 0.006,-0.008,-0.018,-0.026,
     .-0.032,-0.036,-0.043,-0.048,-0.050,-0.052,-0.053,-0.053,
     .-0.053,-0.051,-0.051,-0.054,-0.058,-0.061,-0.063,-0.065,
     .-0.066,-0.072,-0.084,-0.093,-0.101,-0.106,-0.107,-0.112,
     .-0.116,-0.119,-0.120,-0.124,-0.131,-0.138,-0.139,-0.137,
     .-0.132,-0.129,-0.124,-0.119,-0.112,-0.100,-0.086,-0.069,
     .-0.050,
c YARR; ROW j = 45
     . 0.086, 0.083, 0.080, 0.078, 0.075, 0.071, 0.069, 0.067,
     . 0.063, 0.057, 0.049, 0.044, 0.039, 0.036, 0.034, 0.027,
     . 0.019, 0.009, 0.000,-0.008,-0.020,-0.029,-0.039,-0.044,
     .-0.049,-0.054,-0.059,-0.063,-0.066,-0.069,-0.072,-0.073,
     .-0.072,-0.070,-0.070,-0.072,-0.075,-0.079,-0.080,-0.082,
     .-0.081,-0.087,-0.097,-0.106,-0.114,-0.119,-0.121,-0.126,
     .-0.129,-0.130,-0.132,-0.135,-0.141,-0.146,-0.146,-0.144,
     .-0.137,-0.131,-0.126,-0.120,-0.111,-0.098,-0.082,-0.062,
     .-0.038,
c YARR; ROW j = 46
     . 0.081, 0.077, 0.073, 0.069, 0.061, 0.058, 0.054, 0.050,
     . 0.045, 0.038, 0.031, 0.026, 0.021, 0.017, 0.014, 0.006,
     .-0.005,-0.014,-0.025,-0.032,-0.040,-0.050,-0.058,-0.062,
     .-0.065,-0.068,-0.073,-0.078,-0.083,-0.088,-0.090,-0.092,
     .-0.090,-0.087,-0.087,-0.088,-0.091,-0.092,-0.093,-0.095,
     .-0.094,-0.100,-0.110,-0.120,-0.129,-0.134,-0.135,-0.139,
     .-0.143,-0.144,-0.144,-0.143,-0.148,-0.153,-0.151,-0.147,
     .-0.140,-0.133,-0.128,-0.120,-0.107,-0.094,-0.075,-0.053,
     .-0.027,
c YARR; ROW j = 47
     . 0.077, 0.071, 0.066, 0.060, 0.053, 0.047, 0.039, 0.034,
     . 0.029, 0.023, 0.015, 0.008, 0.003, 0.000,-0.003,-0.013,
     .-0.025,-0.035,-0.045,-0.051,-0.058,-0.066,-0.072,-0.076,
     .-0.079,-0.081,-0.086,-0.092,-0.099,-0.107,-0.110,-0.111,
     .-0.108,-0.104,-0.102,-0.102,-0.103,-0.104,-0.105,-0.107,
     .-0.107,-0.114,-0.124,-0.135,-0.145,-0.151,-0.153,-0.157,
     .-0.158,-0.157,-0.156,-0.155,-0.158,-0.161,-0.159,-0.152,
     .-0.145,-0.137,-0.130,-0.121,-0.107,-0.093,-0.071,-0.046,
     .-0.016,
c YARR; ROW j = 48
     . 0.071, 0.064, 0.057, 0.050, 0.045, 0.036, 0.028, 0.020,
     . 0.015, 0.007, 0.001,-0.006,-0.013,-0.017,-0.021,-0.032,
     .-0.043,-0.053,-0.063,-0.068,-0.073,-0.081,-0.086,-0.089,
     .-0.092,-0.097,-0.105,-0.113,-0.120,-0.128,-0.133,-0.134,
     .-0.129,-0.122,-0.117,-0.117,-0.118,-0.120,-0.121,-0.123,
     .-0.124,-0.131,-0.143,-0.154,-0.164,-0.169,-0.170,-0.172,
     .-0.172,-0.169,-0.168,-0.164,-0.167,-0.166,-0.162,-0.156,
     .-0.148,-0.139,-0.130,-0.118,-0.103,-0.089,-0.066,-0.038,
     .-0.007,
c YARR; ROW j = 49
     . 0.067, 0.060, 0.053, 0.045, 0.039, 0.029, 0.019, 0.011,
     . 0.002,-0.005,-0.012,-0.020,-0.025,-0.031,-0.039,-0.049,
     .-0.058,-0.066,-0.074,-0.080,-0.084,-0.091,-0.096,-0.101,
     .-0.105,-0.112,-0.123,-0.133,-0.141,-0.148,-0.152,-0.150,
     .-0.144,-0.135,-0.131,-0.129,-0.130,-0.131,-0.133,-0.137,
     .-0.139,-0.147,-0.157,-0.166,-0.175,-0.179,-0.179,-0.180,
     .-0.179,-0.174,-0.171,-0.166,-0.166,-0.164,-0.160,-0.151,
     .-0.140,-0.130,-0.120,-0.108,-0.092,-0.072,-0.048,-0.021,
     . 0.010,
c YARR; ROW j = 50
     . 0.055, 0.049, 0.042, 0.034, 0.025, 0.014, 0.005,-0.005,
     .-0.013,-0.024,-0.032,-0.039,-0.046,-0.052,-0.059,-0.068,
     .-0.076,-0.081,-0.087,-0.091,-0.096,-0.102,-0.109,-0.114,
     .-0.122,-0.132,-0.144,-0.154,-0.160,-0.166,-0.168,-0.164,
     .-0.156,-0.149,-0.142,-0.141,-0.142,-0.145,-0.147,-0.152,
     .-0.153,-0.158,-0.165,-0.172,-0.178,-0.181,-0.179,-0.178,
     .-0.177,-0.172,-0.166,-0.161,-0.158,-0.155,-0.148,-0.136,
     .-0.124,-0.112,-0.100,-0.087,-0.070,-0.047,-0.022, 0.006,
     . 0.038,
c YARR; ROW j = 51
     . 0.033, 0.028, 0.021, 0.014, 0.003,-0.008,-0.019,-0.027,
     .-0.037,-0.048,-0.056,-0.063,-0.068,-0.073,-0.079,-0.085,
     .-0.090,-0.095,-0.097,-0.100,-0.105,-0.112,-0.121,-0.131,
     .-0.140,-0.151,-0.162,-0.170,-0.176,-0.179,-0.177,-0.171,
     .-0.165,-0.156,-0.152,-0.151,-0.155,-0.157,-0.160,-0.162,
     .-0.161,-0.163,-0.166,-0.170,-0.173,-0.173,-0.171,-0.170,
     .-0.167,-0.163,-0.158,-0.150,-0.146,-0.139,-0.129,-0.116,
     .-0.103,-0.089,-0.076,-0.059,-0.040,-0.019, 0.007, 0.035,
     . 0.067,
c YARR; ROW j = 52
     . 0.011, 0.007, 0.001,-0.006,-0.018,-0.027,-0.038,-0.048,
     .-0.058,-0.068,-0.076,-0.083,-0.088,-0.090,-0.092,-0.096,
     .-0.100,-0.101,-0.103,-0.105,-0.110,-0.120,-0.132,-0.142,
     .-0.151,-0.160,-0.168,-0.174,-0.178,-0.179,-0.174,-0.169,
     .-0.163,-0.158,-0.157,-0.159,-0.163,-0.166,-0.167,-0.167,
     .-0.164,-0.164,-0.164,-0.163,-0.163,-0.161,-0.157,-0.157,
     .-0.154,-0.149,-0.142,-0.134,-0.127,-0.120,-0.107,-0.092,
     .-0.077,-0.063,-0.049,-0.031,-0.011, 0.012, 0.037, 0.065,
     . 0.095,
c YARR; ROW j = 53
     . 0.002,-0.003,-0.010,-0.018,-0.027,-0.038,-0.050,-0.058,
     .-0.068,-0.078,-0.086,-0.092,-0.095,-0.094,-0.092,-0.094,
     .-0.097,-0.098,-0.101,-0.106,-0.113,-0.124,-0.135,-0.145,
     .-0.152,-0.159,-0.163,-0.165,-0.165,-0.166,-0.164,-0.161,
     .-0.158,-0.156,-0.157,-0.161,-0.165,-0.169,-0.168,-0.167,
     .-0.162,-0.159,-0.157,-0.155,-0.152,-0.149,-0.146,-0.143,
     .-0.141,-0.134,-0.127,-0.119,-0.110,-0.099,-0.085,-0.069,
     .-0.054,-0.039,-0.022,-0.004, 0.017, 0.038, 0.064, 0.093,
     . 0.125,
c YARR; ROW j = 54
     .-0.001,-0.006,-0.012,-0.020,-0.028,-0.039,-0.050,-0.059,
     .-0.069,-0.080,-0.086,-0.089,-0.088,-0.086,-0.083,-0.085,
     .-0.086,-0.089,-0.096,-0.103,-0.112,-0.124,-0.133,-0.140,
     .-0.145,-0.147,-0.150,-0.150,-0.150,-0.149,-0.150,-0.151,
     .-0.152,-0.153,-0.156,-0.161,-0.164,-0.166,-0.164,-0.161,
     .-0.157,-0.154,-0.149,-0.146,-0.142,-0.138,-0.133,-0.129,
     .-0.125,-0.118,-0.110,-0.100,-0.090,-0.076,-0.061,-0.044,
     .-0.028,-0.014, 0.003, 0.022, 0.044, 0.067, 0.092, 0.120,
     . 0.149,
c YARR; ROW j = 55
     . 0.013, 0.007, 0.000,-0.008,-0.015,-0.027,-0.037,-0.048,
     .-0.059,-0.069,-0.075,-0.076,-0.072,-0.068,-0.067,-0.068,
     .-0.073,-0.080,-0.087,-0.098,-0.107,-0.116,-0.124,-0.128,
     .-0.128,-0.129,-0.128,-0.128,-0.129,-0.132,-0.135,-0.138,
     .-0.142,-0.145,-0.149,-0.153,-0.155,-0.154,-0.152,-0.150,
     .-0.145,-0.142,-0.139,-0.135,-0.131,-0.125,-0.117,-0.112,
     .-0.107,-0.098,-0.087,-0.075,-0.062,-0.050,-0.033,-0.016,
     . 0.002, 0.017, 0.032, 0.050, 0.073, 0.099, 0.125, 0.152,
     . 0.181,
c YARR; ROW j = 56
     . 0.027, 0.021, 0.015, 0.006,-0.006,-0.016,-0.027,-0.037,
     .-0.048,-0.058,-0.061,-0.060,-0.056,-0.053,-0.051,-0.055,
     .-0.063,-0.073,-0.085,-0.096,-0.104,-0.110,-0.114,-0.113,
     .-0.111,-0.110,-0.111,-0.112,-0.114,-0.120,-0.125,-0.130,
     .-0.134,-0.138,-0.139,-0.141,-0.140,-0.138,-0.137,-0.135,
     .-0.130,-0.126,-0.123,-0.119,-0.116,-0.108,-0.099,-0.092,
     .-0.082,-0.073,-0.062,-0.049,-0.035,-0.020,-0.001, 0.016,
     . 0.033, 0.048, 0.065, 0.084, 0.107, 0.132, 0.159, 0.188,
     . 0.218,
c YARR; ROW j = 57
     . 0.037, 0.031, 0.023, 0.013, 0.002,-0.011,-0.023,-0.031,
     .-0.041,-0.048,-0.050,-0.049,-0.043,-0.037,-0.037,-0.044,
     .-0.056,-0.070,-0.080,-0.091,-0.097,-0.100,-0.101,-0.098,
     .-0.093,-0.092,-0.094,-0.098,-0.106,-0.112,-0.117,-0.123,
     .-0.126,-0.127,-0.127,-0.127,-0.126,-0.124,-0.120,-0.117,
     .-0.113,-0.111,-0.107,-0.101,-0.095,-0.086,-0.076,-0.069,
     .-0.057,-0.045,-0.033,-0.021,-0.007, 0.009, 0.030, 0.049,
     . 0.065, 0.080, 0.097, 0.117, 0.141, 0.168, 0.196, 0.225,
     . 0.256,
c YARR; ROW j = 58
     . 0.044, 0.040, 0.034, 0.025, 0.010,-0.004,-0.016,-0.024,
     .-0.029,-0.034,-0.035,-0.032,-0.028,-0.024,-0.025,-0.035,
     .-0.050,-0.063,-0.075,-0.083,-0.085,-0.085,-0.083,-0.078,
     .-0.075,-0.076,-0.081,-0.088,-0.096,-0.104,-0.109,-0.115,
     .-0.114,-0.113,-0.112,-0.111,-0.110,-0.108,-0.102,-0.098,
     .-0.094,-0.091,-0.087,-0.079,-0.071,-0.061,-0.050,-0.040,
     .-0.028,-0.015,-0.003, 0.011, 0.026, 0.043, 0.063, 0.082,
     . 0.099, 0.114, 0.132, 0.153, 0.177, 0.204, 0.233, 0.264,
     . 0.297,
c YARR; ROW j = 59
     . 0.053, 0.050, 0.043, 0.035, 0.021, 0.006,-0.006,-0.011,
     .-0.015,-0.017,-0.016,-0.014,-0.010,-0.008,-0.013,-0.025,
     .-0.040,-0.053,-0.063,-0.069,-0.069,-0.066,-0.061,-0.058,
     .-0.057,-0.062,-0.069,-0.077,-0.087,-0.095,-0.099,-0.100,
     .-0.098,-0.094,-0.093,-0.092,-0.089,-0.086,-0.080,-0.075,
     .-0.071,-0.067,-0.060,-0.053,-0.043,-0.033,-0.021,-0.010,
     . 0.004, 0.017, 0.030, 0.043, 0.058, 0.075, 0.096, 0.116,
     . 0.133, 0.151, 0.170, 0.192, 0.217, 0.244, 0.274, 0.306,
     . 0.341,
c YARR; ROW j = 60
     . 0.064, 0.060, 0.055, 0.046, 0.033, 0.019, 0.009, 0.005,
     . 0.002, 0.001, 0.001, 0.005, 0.007, 0.005,-0.002,-0.016,
     .-0.029,-0.038,-0.047,-0.050,-0.047,-0.043,-0.039,-0.037,
     .-0.039,-0.046,-0.055,-0.064,-0.073,-0.080,-0.082,-0.081,
     .-0.076,-0.072,-0.069,-0.068,-0.065,-0.061,-0.055,-0.050,
     .-0.045,-0.040,-0.031,-0.020,-0.011, 0.002, 0.014, 0.026,
     . 0.039, 0.054, 0.067, 0.081, 0.096, 0.111, 0.132, 0.153,
     . 0.172, 0.191, 0.211, 0.234, 0.261, 0.290, 0.322, 0.355,
     . 0.390,
c YARR; ROW j = 61
     . 0.078, 0.073, 0.066, 0.058, 0.046, 0.033, 0.027, 0.024,
     . 0.022, 0.020, 0.018, 0.020, 0.020, 0.016, 0.009,-0.004,
     .-0.016,-0.023,-0.027,-0.029,-0.025,-0.020,-0.018,-0.018,
     .-0.023,-0.032,-0.039,-0.046,-0.053,-0.058,-0.059,-0.056,
     .-0.049,-0.044,-0.041,-0.039,-0.036,-0.033,-0.027,-0.021,
     .-0.014,-0.006, 0.003, 0.013, 0.025, 0.038, 0.052, 0.063,
     . 0.078, 0.093, 0.107, 0.121, 0.135, 0.151, 0.173, 0.193,
     . 0.214, 0.235, 0.256, 0.281, 0.308, 0.339, 0.371, 0.407,
     . 0.445,
c YARR; ROW j = 62
     . 0.085, 0.083, 0.079, 0.072, 0.059, 0.050, 0.044, 0.044,
     . 0.042, 0.038, 0.036, 0.035, 0.032, 0.028, 0.022, 0.011,
     . 0.003,-0.001,-0.004,-0.004, 0.000, 0.001, 0.003, 0.001,
     .-0.005,-0.012,-0.021,-0.025,-0.027,-0.029,-0.028,-0.024,
     .-0.017,-0.011,-0.007,-0.006,-0.003, 0.002, 0.007, 0.013,
     . 0.020, 0.029, 0.040, 0.053, 0.064, 0.077, 0.091, 0.104,
     . 0.120, 0.136, 0.149, 0.164, 0.178, 0.194, 0.218, 0.240,
     . 0.262, 0.284, 0.306, 0.331, 0.361, 0.392, 0.427, 0.465,
     . 0.505,
c YARR; ROW j = 63
     . 0.095, 0.095, 0.092, 0.088, 0.079, 0.069, 0.065, 0.066,
     . 0.062, 0.058, 0.052, 0.047, 0.047, 0.043, 0.038, 0.031,
     . 0.027, 0.025, 0.028, 0.026, 0.026, 0.025, 0.022, 0.019,
     . 0.014, 0.008, 0.005, 0.005, 0.003, 0.004, 0.005, 0.011,
     . 0.017, 0.023, 0.025, 0.027, 0.030, 0.034, 0.042, 0.048,
     . 0.055, 0.064, 0.076, 0.089, 0.103, 0.116, 0.129, 0.143,
     . 0.161, 0.179, 0.193, 0.207, 0.221, 0.239, 0.263, 0.288,
     . 0.312, 0.336, 0.360, 0.386, 0.418, 0.453, 0.488, 0.525,
     . 0.564,
c YARR; ROW j = 64
     . 0.106, 0.107, 0.107, 0.104, 0.097, 0.089, 0.086, 0.087,
     . 0.083, 0.079, 0.072, 0.066, 0.065, 0.063, 0.060, 0.055,
     . 0.051, 0.051, 0.053, 0.051, 0.049, 0.046, 0.042, 0.039,
     . 0.035, 0.030, 0.028, 0.030, 0.030, 0.032, 0.035, 0.040,
     . 0.048, 0.055, 0.059, 0.061, 0.065, 0.070, 0.079, 0.085,
     . 0.093, 0.103, 0.115, 0.129, 0.143, 0.158, 0.174, 0.189,
     . 0.209, 0.228, 0.244, 0.259, 0.274, 0.294, 0.320, 0.346,
     . 0.372, 0.398, 0.424, 0.452, 0.485, 0.521, 0.558, 0.596,
     . 0.636,
c YARR; ROW j = 65
     . 0.117, 0.121, 0.122, 0.121, 0.116, 0.111, 0.109, 0.110,
     . 0.107, 0.101, 0.094, 0.089, 0.087, 0.086, 0.085, 0.081,
     . 0.079, 0.080, 0.079, 0.076, 0.073, 0.067, 0.061, 0.058,
     . 0.055, 0.052, 0.052, 0.055, 0.057, 0.060, 0.063, 0.069,
     . 0.078, 0.087, 0.092, 0.097, 0.102, 0.109, 0.118, 0.126,
     . 0.134, 0.143, 0.156, 0.170, 0.185, 0.202, 0.220, 0.238,
     . 0.260, 0.282, 0.299, 0.316, 0.333, 0.355, 0.383, 0.412,
     . 0.440, 0.468, 0.495, 0.527, 0.561, 0.597, 0.635, 0.674,
     . 0.715,
c YARR; ROW j = 66
     . 0.130, 0.136, 0.139, 0.140, 0.136, 0.135, 0.134, 0.135,
     . 0.132, 0.126, 0.118, 0.114, 0.111, 0.112, 0.113, 0.112,
     . 0.109, 0.111, 0.106, 0.101, 0.096, 0.087, 0.081, 0.076,
     . 0.076, 0.076, 0.077, 0.081, 0.085, 0.086, 0.091, 0.095,
     . 0.107, 0.118, 0.126, 0.133, 0.141, 0.150, 0.160, 0.168,
     . 0.176, 0.185, 0.198, 0.212, 0.229, 0.248, 0.270, 0.291,
     . 0.316, 0.339, 0.360, 0.379, 0.398, 0.422, 0.452, 0.484,
     . 0.515, 0.545, 0.575, 0.609, 0.645, 0.680, 0.719, 0.759,
     . 0.802/
      end
 
 
c----------------------------------------
c
c
c
      real function wfc3uv_zpa(xr,yr)
      implicit none
 
      real*8 xr, yr
 
      real    ri,rj
      integer ii,jj
      real    fi,fj
 
      real wfc3uv_zp_arr(64,64)
      common /wfc3uv_zp_arr_/wfc3uv_zp_arr
 
      ri = 1 + (xr-32)/64.0
      rj = 1 + (yr-32)/64.0
      if (ri.le. 1.0) ri =  1.00
      if (rj.le. 1.0) rj =  1.00
      if (ri.ge.64.0) ri = 63.99
      if (rj.ge.64.0) rj = 63.99
      ii = int(ri)
      jj = int(rj)
      fi = ri-ii
      fj = rj-jj
 
      wfc3uv_zpa = (1-fi)*(1-fj)*wfc3uv_zp_arr(ii  ,jj  )
     .           + (1-fi)*( fj )*wfc3uv_zp_arr(ii  ,jj+1)
     .           + ( fi )*(1-fj)*wfc3uv_zp_arr(ii+1,jj  )
     .           + ( fi )*( fj )*wfc3uv_zp_arr(ii+1,jj+1)
 
      return
      end
 
 
 
c--------------------------------
c
c
      block data wfc3uv_zp_arr_BD
      real wfc3uv_zp_arr(64,64)
      common /wfc3uv_zp_arr_/wfc3uv_zp_arr
      data wfc3uv_zp_arr/
     .-0.0008,-0.0014,-0.0020,-0.0026,-0.0033,-0.0038,-0.0044,-0.0051,
     .-0.0055,-0.0060,-0.0066,-0.0071,-0.0078,-0.0084,-0.0090,-0.0095,
     .-0.0101,-0.0107,-0.0113,-0.0119,-0.0126,-0.0132,-0.0137,-0.0143,
     .-0.0148,-0.0155,-0.0161,-0.0167,-0.0172,-0.0177,-0.0183,-0.0188,
     .-0.0194,-0.0200,-0.0208,-0.0214,-0.0220,-0.0225,-0.0230,-0.0236,
     .-0.0242,-0.0247,-0.0252,-0.0258,-0.0264,-0.0271,-0.0277,-0.0283,
     .-0.0289,-0.0295,-0.0302,-0.0308,-0.0314,-0.0319,-0.0326,-0.0332,
     .-0.0339,-0.0345,-0.0351,-0.0356,-0.0361,-0.0368,-0.0374,-0.0380,
     .-0.0001,-0.0007,-0.0013,-0.0018,-0.0025,-0.0031,-0.0037,-0.0043,
     .-0.0048,-0.0053,-0.0058,-0.0064,-0.0071,-0.0077,-0.0082,-0.0088,
     .-0.0094,-0.0100,-0.0106,-0.0112,-0.0119,-0.0125,-0.0130,-0.0136,
     .-0.0141,-0.0148,-0.0154,-0.0160,-0.0165,-0.0170,-0.0176,-0.0181,
     .-0.0188,-0.0194,-0.0201,-0.0207,-0.0213,-0.0218,-0.0223,-0.0229,
     .-0.0235,-0.0240,-0.0246,-0.0252,-0.0258,-0.0265,-0.0270,-0.0276,
     .-0.0283,-0.0289,-0.0295,-0.0302,-0.0307,-0.0313,-0.0319,-0.0326,
     .-0.0332,-0.0338,-0.0345,-0.0350,-0.0355,-0.0361,-0.0368,-0.0374,
     . 0.0006, 0.0000,-0.0006,-0.0011,-0.0018,-0.0024,-0.0030,-0.0036,
     .-0.0041,-0.0046,-0.0051,-0.0057,-0.0063,-0.0070,-0.0075,-0.0081,
     .-0.0086,-0.0093,-0.0099,-0.0105,-0.0112,-0.0118,-0.0123,-0.0129,
     .-0.0134,-0.0141,-0.0147,-0.0153,-0.0158,-0.0164,-0.0169,-0.0175,
     .-0.0181,-0.0187,-0.0195,-0.0201,-0.0207,-0.0212,-0.0217,-0.0223,
     .-0.0229,-0.0234,-0.0240,-0.0246,-0.0251,-0.0259,-0.0264,-0.0270,
     .-0.0276,-0.0283,-0.0289,-0.0295,-0.0301,-0.0307,-0.0313,-0.0319,
     .-0.0326,-0.0332,-0.0338,-0.0344,-0.0349,-0.0355,-0.0361,-0.0367,
     . 0.0014, 0.0008, 0.0001,-0.0004,-0.0010,-0.0017,-0.0022,-0.0028,
     .-0.0033,-0.0038,-0.0043,-0.0050,-0.0056,-0.0063,-0.0069,-0.0074,
     .-0.0079,-0.0086,-0.0092,-0.0098,-0.0105,-0.0111,-0.0116,-0.0122,
     .-0.0127,-0.0134,-0.0140,-0.0145,-0.0150,-0.0156,-0.0161,-0.0167,
     .-0.0174,-0.0180,-0.0188,-0.0194,-0.0199,-0.0205,-0.0210,-0.0216,
     .-0.0221,-0.0227,-0.0233,-0.0239,-0.0244,-0.0252,-0.0258,-0.0264,
     .-0.0270,-0.0276,-0.0282,-0.0288,-0.0294,-0.0300,-0.0307,-0.0313,
     .-0.0319,-0.0325,-0.0331,-0.0337,-0.0342,-0.0349,-0.0354,-0.0360,
     . 0.0020, 0.0015, 0.0009, 0.0003,-0.0003,-0.0009,-0.0015,-0.0021,
     .-0.0027,-0.0031,-0.0036,-0.0042,-0.0049,-0.0055,-0.0061,-0.0066,
     .-0.0072,-0.0078,-0.0084,-0.0091,-0.0097,-0.0104,-0.0110,-0.0115,
     .-0.0120,-0.0127,-0.0133,-0.0139,-0.0144,-0.0149,-0.0155,-0.0161,
     .-0.0168,-0.0174,-0.0182,-0.0188,-0.0193,-0.0199,-0.0204,-0.0210,
     .-0.0216,-0.0221,-0.0227,-0.0233,-0.0239,-0.0246,-0.0252,-0.0258,
     .-0.0264,-0.0270,-0.0277,-0.0283,-0.0289,-0.0295,-0.0301,-0.0307,
     .-0.0314,-0.0320,-0.0326,-0.0332,-0.0336,-0.0343,-0.0349,-0.0355,
     . 0.0027, 0.0021, 0.0016, 0.0010, 0.0003,-0.0002,-0.0008,-0.0015,
     .-0.0020,-0.0024,-0.0030,-0.0036,-0.0043,-0.0049,-0.0055,-0.0060,
     .-0.0066,-0.0072,-0.0078,-0.0084,-0.0091,-0.0097,-0.0103,-0.0109,
     .-0.0114,-0.0121,-0.0127,-0.0133,-0.0139,-0.0144,-0.0150,-0.0156,
     .-0.0163,-0.0169,-0.0176,-0.0182,-0.0188,-0.0193,-0.0198,-0.0205,
     .-0.0210,-0.0216,-0.0221,-0.0227,-0.0233,-0.0240,-0.0246,-0.0252,
     .-0.0258,-0.0264,-0.0270,-0.0276,-0.0282,-0.0288,-0.0294,-0.0300,
     .-0.0307,-0.0313,-0.0320,-0.0325,-0.0329,-0.0337,-0.0343,-0.0349,
     . 0.0033, 0.0027, 0.0021, 0.0016, 0.0009, 0.0003,-0.0002,-0.0009,
     .-0.0014,-0.0019,-0.0024,-0.0030,-0.0037,-0.0044,-0.0049,-0.0054,
     .-0.0060,-0.0066,-0.0071,-0.0078,-0.0085,-0.0091,-0.0097,-0.0103,
     .-0.0108,-0.0115,-0.0121,-0.0127,-0.0133,-0.0139,-0.0145,-0.0151,
     .-0.0156,-0.0163,-0.0170,-0.0176,-0.0182,-0.0188,-0.0193,-0.0199,
     .-0.0204,-0.0210,-0.0216,-0.0222,-0.0227,-0.0234,-0.0240,-0.0245,
     .-0.0251,-0.0257,-0.0264,-0.0270,-0.0276,-0.0281,-0.0288,-0.0294,
     .-0.0301,-0.0307,-0.0313,-0.0319,-0.0323,-0.0330,-0.0337,-0.0343,
     . 0.0038, 0.0033, 0.0027, 0.0022, 0.0015, 0.0009, 0.0003,-0.0003,
     .-0.0008,-0.0013,-0.0018,-0.0024,-0.0030,-0.0037,-0.0043,-0.0049,
     .-0.0054,-0.0060,-0.0065,-0.0072,-0.0079,-0.0085,-0.0090,-0.0096,
     .-0.0101,-0.0109,-0.0115,-0.0121,-0.0127,-0.0132,-0.0138,-0.0145,
     .-0.0150,-0.0157,-0.0163,-0.0169,-0.0176,-0.0181,-0.0187,-0.0193,
     .-0.0198,-0.0204,-0.0210,-0.0216,-0.0222,-0.0228,-0.0234,-0.0240,
     .-0.0245,-0.0251,-0.0257,-0.0263,-0.0269,-0.0275,-0.0281,-0.0287,
     .-0.0295,-0.0301,-0.0307,-0.0313,-0.0318,-0.0325,-0.0331,-0.0337,
     . 0.0045, 0.0039, 0.0033, 0.0028, 0.0022, 0.0016, 0.0010, 0.0003,
     .-0.0001,-0.0006,-0.0012,-0.0017,-0.0024,-0.0031,-0.0036,-0.0042,
     .-0.0047,-0.0054,-0.0059,-0.0066,-0.0072,-0.0079,-0.0085,-0.0090,
     .-0.0095,-0.0102,-0.0108,-0.0115,-0.0121,-0.0127,-0.0133,-0.0139,
     .-0.0145,-0.0150,-0.0157,-0.0163,-0.0169,-0.0175,-0.0181,-0.0187,
     .-0.0193,-0.0198,-0.0204,-0.0210,-0.0216,-0.0222,-0.0227,-0.0233,
     .-0.0239,-0.0245,-0.0251,-0.0257,-0.0262,-0.0269,-0.0275,-0.0281,
     .-0.0288,-0.0295,-0.0301,-0.0307,-0.0312,-0.0319,-0.0325,-0.0332,
     . 0.0050, 0.0045, 0.0039, 0.0034, 0.0028, 0.0023, 0.0016, 0.0010,
     . 0.0004, 0.0000,-0.0005,-0.0011,-0.0017,-0.0024,-0.0030,-0.0036,
     .-0.0042,-0.0048,-0.0053,-0.0061,-0.0067,-0.0073,-0.0079,-0.0084,
     .-0.0090,-0.0097,-0.0103,-0.0109,-0.0115,-0.0121,-0.0127,-0.0133,
     .-0.0138,-0.0144,-0.0150,-0.0157,-0.0163,-0.0169,-0.0174,-0.0180,
     .-0.0186,-0.0192,-0.0198,-0.0203,-0.0209,-0.0216,-0.0222,-0.0228,
     .-0.0233,-0.0239,-0.0245,-0.0251,-0.0257,-0.0263,-0.0269,-0.0275,
     .-0.0282,-0.0289,-0.0295,-0.0301,-0.0306,-0.0314,-0.0320,-0.0326,
     . 0.0057, 0.0051, 0.0046, 0.0040, 0.0035, 0.0029, 0.0023, 0.0017,
     . 0.0011, 0.0006, 0.0001,-0.0005,-0.0011,-0.0019,-0.0024,-0.0030,
     .-0.0036,-0.0042,-0.0048,-0.0055,-0.0061,-0.0067,-0.0073,-0.0078,
     .-0.0083,-0.0091,-0.0097,-0.0103,-0.0109,-0.0115,-0.0121,-0.0127,
     .-0.0133,-0.0138,-0.0144,-0.0151,-0.0157,-0.0163,-0.0168,-0.0174,
     .-0.0179,-0.0186,-0.0192,-0.0197,-0.0204,-0.0210,-0.0216,-0.0222,
     .-0.0228,-0.0233,-0.0239,-0.0245,-0.0251,-0.0257,-0.0263,-0.0269,
     .-0.0276,-0.0282,-0.0289,-0.0295,-0.0299,-0.0307,-0.0314,-0.0320,
     . 0.0063, 0.0057, 0.0052, 0.0046, 0.0041, 0.0035, 0.0029, 0.0022,
     . 0.0017, 0.0013, 0.0007, 0.0001,-0.0006,-0.0013,-0.0018,-0.0024,
     .-0.0031,-0.0037,-0.0042,-0.0049,-0.0055,-0.0061,-0.0067,-0.0072,
     .-0.0078,-0.0085,-0.0091,-0.0097,-0.0103,-0.0109,-0.0115,-0.0121,
     .-0.0127,-0.0133,-0.0139,-0.0145,-0.0151,-0.0157,-0.0162,-0.0168,
     .-0.0173,-0.0179,-0.0185,-0.0190,-0.0196,-0.0203,-0.0209,-0.0216,
     .-0.0221,-0.0227,-0.0233,-0.0239,-0.0244,-0.0251,-0.0257,-0.0262,
     .-0.0269,-0.0275,-0.0282,-0.0288,-0.0293,-0.0300,-0.0306,-0.0312,
     . 0.0068, 0.0063, 0.0058, 0.0052, 0.0047, 0.0041, 0.0035, 0.0028,
     . 0.0023, 0.0018, 0.0012, 0.0006, 0.0000,-0.0008,-0.0014,-0.0020,
     .-0.0026,-0.0032,-0.0037,-0.0044,-0.0049,-0.0056,-0.0062,-0.0068,
     .-0.0074,-0.0081,-0.0086,-0.0092,-0.0098,-0.0104,-0.0110,-0.0116,
     .-0.0121,-0.0127,-0.0134,-0.0140,-0.0145,-0.0151,-0.0156,-0.0163,
     .-0.0168,-0.0173,-0.0179,-0.0185,-0.0191,-0.0198,-0.0204,-0.0210,
     .-0.0216,-0.0222,-0.0228,-0.0234,-0.0240,-0.0246,-0.0252,-0.0258,
     .-0.0265,-0.0271,-0.0277,-0.0283,-0.0289,-0.0295,-0.0301,-0.0308,
     . 0.0074, 0.0069, 0.0064, 0.0058, 0.0053, 0.0047, 0.0041, 0.0034,
     . 0.0028, 0.0023, 0.0018, 0.0012, 0.0005,-0.0002,-0.0008,-0.0014,
     .-0.0020,-0.0026,-0.0031,-0.0038,-0.0043,-0.0050,-0.0056,-0.0062,
     .-0.0068,-0.0075,-0.0081,-0.0087,-0.0093,-0.0098,-0.0104,-0.0110,
     .-0.0115,-0.0121,-0.0128,-0.0134,-0.0140,-0.0146,-0.0151,-0.0157,
     .-0.0162,-0.0168,-0.0173,-0.0179,-0.0185,-0.0192,-0.0198,-0.0205,
     .-0.0211,-0.0216,-0.0223,-0.0229,-0.0235,-0.0242,-0.0248,-0.0253,
     .-0.0260,-0.0266,-0.0272,-0.0278,-0.0284,-0.0290,-0.0295,-0.0301,
     . 0.0080, 0.0075, 0.0070, 0.0065, 0.0059, 0.0053, 0.0047, 0.0040,
     . 0.0035, 0.0030, 0.0024, 0.0018, 0.0012, 0.0004,-0.0002,-0.0008,
     .-0.0013,-0.0019,-0.0025,-0.0031,-0.0037,-0.0044,-0.0050,-0.0056,
     .-0.0062,-0.0069,-0.0075,-0.0081,-0.0087,-0.0092,-0.0098,-0.0104,
     .-0.0109,-0.0115,-0.0122,-0.0128,-0.0133,-0.0139,-0.0144,-0.0151,
     .-0.0156,-0.0161,-0.0167,-0.0173,-0.0179,-0.0187,-0.0192,-0.0198,
     .-0.0205,-0.0210,-0.0216,-0.0223,-0.0229,-0.0235,-0.0242,-0.0248,
     .-0.0254,-0.0260,-0.0266,-0.0272,-0.0278,-0.0283,-0.0289,-0.0295,
     . 0.0086, 0.0081, 0.0076, 0.0071, 0.0065, 0.0059, 0.0053, 0.0046,
     . 0.0041, 0.0035, 0.0030, 0.0024, 0.0017, 0.0011, 0.0004,-0.0002,
     .-0.0008,-0.0013,-0.0019,-0.0025,-0.0031,-0.0037,-0.0043,-0.0050,
     .-0.0056,-0.0063,-0.0069,-0.0075,-0.0080,-0.0085,-0.0091,-0.0097,
     .-0.0103,-0.0109,-0.0115,-0.0121,-0.0127,-0.0133,-0.0138,-0.0145,
     .-0.0149,-0.0155,-0.0161,-0.0167,-0.0173,-0.0180,-0.0186,-0.0192,
     .-0.0198,-0.0203,-0.0209,-0.0216,-0.0223,-0.0228,-0.0235,-0.0240,
     .-0.0246,-0.0252,-0.0258,-0.0264,-0.0271,-0.0276,-0.0281,-0.0286,
     . 0.0093, 0.0087, 0.0082, 0.0077, 0.0071, 0.0065, 0.0059, 0.0052,
     . 0.0047, 0.0042, 0.0036, 0.0030, 0.0024, 0.0017, 0.0011, 0.0005,
     .-0.0001,-0.0007,-0.0013,-0.0018,-0.0024,-0.0030,-0.0036,-0.0042,
     .-0.0049,-0.0055,-0.0061,-0.0067,-0.0072,-0.0078,-0.0084,-0.0089,
     .-0.0095,-0.0102,-0.0108,-0.0114,-0.0120,-0.0125,-0.0131,-0.0137,
     .-0.0142,-0.0148,-0.0153,-0.0159,-0.0165,-0.0173,-0.0178,-0.0184,
     .-0.0190,-0.0196,-0.0201,-0.0208,-0.0214,-0.0221,-0.0227,-0.0232,
     .-0.0238,-0.0244,-0.0250,-0.0256,-0.0262,-0.0267,-0.0272,-0.0277,
     . 0.0099, 0.0094, 0.0088, 0.0082, 0.0077, 0.0071, 0.0066, 0.0059,
     . 0.0053, 0.0048, 0.0043, 0.0037, 0.0031, 0.0024, 0.0018, 0.0013,
     . 0.0006, 0.0000,-0.0006,-0.0011,-0.0017,-0.0022,-0.0028,-0.0034,
     .-0.0041,-0.0047,-0.0053,-0.0059,-0.0065,-0.0070,-0.0076,-0.0082,
     .-0.0088,-0.0094,-0.0101,-0.0107,-0.0113,-0.0118,-0.0123,-0.0130,
     .-0.0135,-0.0140,-0.0145,-0.0151,-0.0157,-0.0165,-0.0170,-0.0176,
     .-0.0182,-0.0188,-0.0194,-0.0200,-0.0207,-0.0212,-0.0219,-0.0224,
     .-0.0229,-0.0236,-0.0242,-0.0247,-0.0254,-0.0259,-0.0264,-0.0269,
     . 0.0106, 0.0100, 0.0095, 0.0089, 0.0083, 0.0078, 0.0072, 0.0066,
     . 0.0060, 0.0054, 0.0049, 0.0043, 0.0037, 0.0031, 0.0025, 0.0019,
     . 0.0013, 0.0007, 0.0000,-0.0005,-0.0011,-0.0016,-0.0022,-0.0027,
     .-0.0034,-0.0041,-0.0047,-0.0052,-0.0058,-0.0063,-0.0069,-0.0075,
     .-0.0081,-0.0087,-0.0094,-0.0100,-0.0106,-0.0112,-0.0117,-0.0124,
     .-0.0128,-0.0133,-0.0139,-0.0145,-0.0151,-0.0159,-0.0164,-0.0169,
     .-0.0176,-0.0182,-0.0187,-0.0194,-0.0201,-0.0207,-0.0213,-0.0218,
     .-0.0224,-0.0230,-0.0236,-0.0242,-0.0248,-0.0253,-0.0258,-0.0264,
     . 0.0113, 0.0108, 0.0102, 0.0096, 0.0090, 0.0085, 0.0079, 0.0073,
     . 0.0067, 0.0061, 0.0056, 0.0050, 0.0043, 0.0037, 0.0031, 0.0026,
     . 0.0020, 0.0014, 0.0007, 0.0001,-0.0004,-0.0009,-0.0015,-0.0021,
     .-0.0027,-0.0034,-0.0040,-0.0046,-0.0052,-0.0057,-0.0063,-0.0069,
     .-0.0074,-0.0081,-0.0088,-0.0093,-0.0099,-0.0105,-0.0111,-0.0118,
     .-0.0122,-0.0127,-0.0133,-0.0139,-0.0145,-0.0153,-0.0158,-0.0164,
     .-0.0170,-0.0176,-0.0181,-0.0189,-0.0195,-0.0201,-0.0207,-0.0213,
     .-0.0218,-0.0225,-0.0231,-0.0237,-0.0243,-0.0247,-0.0252,-0.0257,
     . 0.0120, 0.0114, 0.0109, 0.0102, 0.0097, 0.0091, 0.0086, 0.0080,
     . 0.0073, 0.0068, 0.0063, 0.0057, 0.0050, 0.0044, 0.0038, 0.0032,
     . 0.0026, 0.0020, 0.0013, 0.0007, 0.0002,-0.0003,-0.0009,-0.0015,
     .-0.0022,-0.0028,-0.0034,-0.0040,-0.0046,-0.0051,-0.0057,-0.0063,
     .-0.0069,-0.0075,-0.0082,-0.0087,-0.0094,-0.0100,-0.0105,-0.0113,
     .-0.0117,-0.0122,-0.0128,-0.0133,-0.0140,-0.0148,-0.0152,-0.0158,
     .-0.0164,-0.0170,-0.0176,-0.0183,-0.0189,-0.0195,-0.0201,-0.0207,
     .-0.0213,-0.0219,-0.0225,-0.0231,-0.0238,-0.0241,-0.0247,-0.0252,
     . 0.0127, 0.0121, 0.0115, 0.0109, 0.0103, 0.0098, 0.0092, 0.0086,
     . 0.0080, 0.0074, 0.0068, 0.0062, 0.0056, 0.0050, 0.0044, 0.0038,
     . 0.0032, 0.0027, 0.0019, 0.0014, 0.0007, 0.0002,-0.0003,-0.0009,
     .-0.0016,-0.0023,-0.0029,-0.0035,-0.0041,-0.0046,-0.0052,-0.0058,
     .-0.0063,-0.0069,-0.0076,-0.0082,-0.0089,-0.0094,-0.0100,-0.0107,
     .-0.0112,-0.0117,-0.0123,-0.0128,-0.0135,-0.0142,-0.0147,-0.0153,
     .-0.0158,-0.0164,-0.0170,-0.0176,-0.0183,-0.0189,-0.0195,-0.0200,
     .-0.0206,-0.0212,-0.0219,-0.0225,-0.0232,-0.0235,-0.0240,-0.0246,
     . 0.0133, 0.0127, 0.0121, 0.0115, 0.0109, 0.0103, 0.0098, 0.0092,
     . 0.0086, 0.0080, 0.0074, 0.0068, 0.0061, 0.0055, 0.0049, 0.0044,
     . 0.0038, 0.0032, 0.0025, 0.0019, 0.0013, 0.0007, 0.0002,-0.0004,
     .-0.0010,-0.0017,-0.0023,-0.0029,-0.0035,-0.0040,-0.0046,-0.0052,
     .-0.0058,-0.0064,-0.0070,-0.0076,-0.0083,-0.0088,-0.0094,-0.0101,
     .-0.0105,-0.0110,-0.0117,-0.0122,-0.0129,-0.0136,-0.0140,-0.0146,
     .-0.0153,-0.0158,-0.0164,-0.0170,-0.0177,-0.0183,-0.0188,-0.0194,
     .-0.0200,-0.0207,-0.0213,-0.0218,-0.0225,-0.0229,-0.0234,-0.0239,
     . 0.0138, 0.0132, 0.0126, 0.0120, 0.0114, 0.0109, 0.0104, 0.0098,
     . 0.0092, 0.0086, 0.0080, 0.0074, 0.0068, 0.0061, 0.0055, 0.0050,
     . 0.0044, 0.0038, 0.0032, 0.0025, 0.0020, 0.0014, 0.0008, 0.0002,
     .-0.0004,-0.0011,-0.0016,-0.0023,-0.0029,-0.0034,-0.0040,-0.0046,
     .-0.0051,-0.0057,-0.0064,-0.0070,-0.0077,-0.0082,-0.0087,-0.0093,
     .-0.0098,-0.0104,-0.0110,-0.0116,-0.0123,-0.0130,-0.0135,-0.0140,
     .-0.0147,-0.0152,-0.0158,-0.0165,-0.0171,-0.0177,-0.0183,-0.0188,
     .-0.0194,-0.0201,-0.0207,-0.0213,-0.0220,-0.0224,-0.0229,-0.0234,
     . 0.0144, 0.0138, 0.0132, 0.0126, 0.0120, 0.0115, 0.0110, 0.0104,
     . 0.0098, 0.0092, 0.0087, 0.0081, 0.0074, 0.0068, 0.0062, 0.0056,
     . 0.0051, 0.0045, 0.0038, 0.0032, 0.0025, 0.0019, 0.0013, 0.0008,
     . 0.0001,-0.0005,-0.0011,-0.0016,-0.0022,-0.0028,-0.0034,-0.0040,
     .-0.0046,-0.0052,-0.0058,-0.0064,-0.0070,-0.0076,-0.0080,-0.0086,
     .-0.0091,-0.0096,-0.0103,-0.0109,-0.0116,-0.0123,-0.0128,-0.0134,
     .-0.0141,-0.0146,-0.0152,-0.0159,-0.0166,-0.0171,-0.0177,-0.0182,
     .-0.0188,-0.0195,-0.0201,-0.0207,-0.0213,-0.0218,-0.0223,-0.0228,
     . 0.0150, 0.0144, 0.0138, 0.0131, 0.0126, 0.0120, 0.0115, 0.0110,
     . 0.0103, 0.0098, 0.0093, 0.0087, 0.0081, 0.0074, 0.0068, 0.0062,
     . 0.0057, 0.0051, 0.0044, 0.0038, 0.0031, 0.0025, 0.0019, 0.0014,
     . 0.0007, 0.0001,-0.0004,-0.0010,-0.0016,-0.0021,-0.0027,-0.0033,
     .-0.0039,-0.0045,-0.0051,-0.0057,-0.0064,-0.0069,-0.0074,-0.0079,
     .-0.0084,-0.0090,-0.0097,-0.0103,-0.0110,-0.0117,-0.0122,-0.0128,
     .-0.0135,-0.0140,-0.0146,-0.0154,-0.0160,-0.0165,-0.0171,-0.0176,
     .-0.0182,-0.0189,-0.0195,-0.0201,-0.0207,-0.0212,-0.0217,-0.0222,
     . 0.0155, 0.0149, 0.0143, 0.0137, 0.0132, 0.0126, 0.0121, 0.0116,
     . 0.0109, 0.0104, 0.0099, 0.0093, 0.0087, 0.0080, 0.0074, 0.0068,
     . 0.0063, 0.0057, 0.0051, 0.0044, 0.0037, 0.0032, 0.0026, 0.0020,
     . 0.0013, 0.0007, 0.0002,-0.0004,-0.0010,-0.0015,-0.0022,-0.0027,
     .-0.0033,-0.0039,-0.0045,-0.0051,-0.0058,-0.0063,-0.0067,-0.0073,
     .-0.0078,-0.0084,-0.0090,-0.0096,-0.0103,-0.0110,-0.0116,-0.0122,
     .-0.0128,-0.0134,-0.0141,-0.0147,-0.0154,-0.0159,-0.0165,-0.0170,
     .-0.0176,-0.0183,-0.0188,-0.0194,-0.0200,-0.0205,-0.0210,-0.0216,
     . 0.0161, 0.0155, 0.0150, 0.0144, 0.0138, 0.0132, 0.0127, 0.0121,
     . 0.0116, 0.0110, 0.0104, 0.0099, 0.0093, 0.0086, 0.0080, 0.0074,
     . 0.0069, 0.0063, 0.0057, 0.0050, 0.0044, 0.0038, 0.0032, 0.0026,
     . 0.0020, 0.0013, 0.0008, 0.0002,-0.0004,-0.0009,-0.0015,-0.0021,
     .-0.0027,-0.0032,-0.0039,-0.0045,-0.0051,-0.0057,-0.0062,-0.0067,
     .-0.0072,-0.0078,-0.0084,-0.0090,-0.0097,-0.0104,-0.0109,-0.0115,
     .-0.0122,-0.0127,-0.0134,-0.0141,-0.0147,-0.0153,-0.0158,-0.0164,
     .-0.0170,-0.0177,-0.0182,-0.0188,-0.0194,-0.0199,-0.0204,-0.0209,
     . 0.0166, 0.0161, 0.0156, 0.0149, 0.0144, 0.0138, 0.0133, 0.0128,
     . 0.0121, 0.0115, 0.0110, 0.0104, 0.0099, 0.0092, 0.0086, 0.0080,
     . 0.0074, 0.0069, 0.0064, 0.0057, 0.0051, 0.0045, 0.0039, 0.0033,
     . 0.0027, 0.0020, 0.0015, 0.0009, 0.0003,-0.0003,-0.0009,-0.0015,
     .-0.0020,-0.0026,-0.0032,-0.0039,-0.0045,-0.0050,-0.0056,-0.0062,
     .-0.0067,-0.0073,-0.0079,-0.0084,-0.0091,-0.0098,-0.0102,-0.0109,
     .-0.0115,-0.0121,-0.0128,-0.0135,-0.0141,-0.0147,-0.0153,-0.0158,
     .-0.0163,-0.0170,-0.0176,-0.0182,-0.0188,-0.0193,-0.0198,-0.0203,
     . 0.0173, 0.0167, 0.0161, 0.0155, 0.0151, 0.0144, 0.0139, 0.0134,
     . 0.0127, 0.0121, 0.0116, 0.0110, 0.0105, 0.0098, 0.0091, 0.0086,
     . 0.0081, 0.0075, 0.0069, 0.0063, 0.0057, 0.0051, 0.0045, 0.0039,
     . 0.0033, 0.0026, 0.0020, 0.0014, 0.0009, 0.0003,-0.0004,-0.0009,
     .-0.0014,-0.0020,-0.0026,-0.0033,-0.0039,-0.0045,-0.0050,-0.0057,
     .-0.0061,-0.0067,-0.0073,-0.0078,-0.0085,-0.0092,-0.0097,-0.0103,
     .-0.0110,-0.0115,-0.0122,-0.0129,-0.0136,-0.0142,-0.0147,-0.0152,
     .-0.0158,-0.0164,-0.0170,-0.0176,-0.0181,-0.0186,-0.0192,-0.0197,
     . 0.0179, 0.0173, 0.0167, 0.0161, 0.0157, 0.0151, 0.0145, 0.0140,
     . 0.0133, 0.0127, 0.0122, 0.0116, 0.0110, 0.0104, 0.0097, 0.0092,
     . 0.0087, 0.0081, 0.0075, 0.0069, 0.0063, 0.0057, 0.0051, 0.0045,
     . 0.0039, 0.0033, 0.0026, 0.0020, 0.0014, 0.0009, 0.0002,-0.0003,
     .-0.0009,-0.0014,-0.0020,-0.0026,-0.0033,-0.0039,-0.0045,-0.0051,
     .-0.0055,-0.0061,-0.0067,-0.0073,-0.0079,-0.0086,-0.0091,-0.0097,
     .-0.0104,-0.0109,-0.0116,-0.0123,-0.0130,-0.0136,-0.0141,-0.0146,
     .-0.0152,-0.0159,-0.0164,-0.0170,-0.0175,-0.0180,-0.0186,-0.0191,
     . 0.0185, 0.0179, 0.0174, 0.0168, 0.0163, 0.0157, 0.0152, 0.0145,
     . 0.0139, 0.0133, 0.0128, 0.0122, 0.0116, 0.0109, 0.0103, 0.0098,
     . 0.0093, 0.0087, 0.0081, 0.0075, 0.0069, 0.0064, 0.0058, 0.0052,
     . 0.0045, 0.0039, 0.0033, 0.0026, 0.0020, 0.0015, 0.0009, 0.0003,
     .-0.0002,-0.0008,-0.0014,-0.0020,-0.0027,-0.0033,-0.0039,-0.0045,
     .-0.0050,-0.0056,-0.0061,-0.0067,-0.0073,-0.0079,-0.0085,-0.0091,
     .-0.0098,-0.0103,-0.0110,-0.0117,-0.0124,-0.0130,-0.0135,-0.0140,
     .-0.0147,-0.0153,-0.0158,-0.0164,-0.0170,-0.0174,-0.0180,-0.0185,
     . 0.0188, 0.0182, 0.0176, 0.0170, 0.0164, 0.0158, 0.0153, 0.0146,
     . 0.0140, 0.0135, 0.0130, 0.0124, 0.0118, 0.0111, 0.0105, 0.0100,
     . 0.0093, 0.0088, 0.0080, 0.0077, 0.0071, 0.0066, 0.0060, 0.0054,
     . 0.0046, 0.0041, 0.0035, 0.0029, 0.0023, 0.0017, 0.0010, 0.0004,
     .-0.0001,-0.0008,-0.0013,-0.0019,-0.0024,-0.0030,-0.0034,-0.0042,
     .-0.0046,-0.0051,-0.0055,-0.0061,-0.0067,-0.0076,-0.0080,-0.0086,
     .-0.0092,-0.0098,-0.0104,-0.0110,-0.0116,-0.0122,-0.0128,-0.0135,
     .-0.0142,-0.0148,-0.0154,-0.0158,-0.0165,-0.0170,-0.0176,-0.0181,
     . 0.0195, 0.0189, 0.0183, 0.0177, 0.0171, 0.0165, 0.0159, 0.0152,
     . 0.0147, 0.0142, 0.0136, 0.0131, 0.0124, 0.0117, 0.0111, 0.0106,
     . 0.0100, 0.0094, 0.0086, 0.0083, 0.0078, 0.0072, 0.0066, 0.0060,
     . 0.0052, 0.0047, 0.0040, 0.0034, 0.0029, 0.0023, 0.0016, 0.0010,
     . 0.0005,-0.0001,-0.0007,-0.0013,-0.0018,-0.0023,-0.0028,-0.0036,
     .-0.0040,-0.0045,-0.0050,-0.0056,-0.0062,-0.0070,-0.0075,-0.0081,
     .-0.0086,-0.0092,-0.0099,-0.0104,-0.0110,-0.0116,-0.0122,-0.0129,
     .-0.0135,-0.0142,-0.0147,-0.0153,-0.0159,-0.0164,-0.0170,-0.0175,
     . 0.0202, 0.0196, 0.0189, 0.0183, 0.0178, 0.0171, 0.0166, 0.0159,
     . 0.0153, 0.0148, 0.0142, 0.0137, 0.0131, 0.0124, 0.0118, 0.0112,
     . 0.0106, 0.0101, 0.0093, 0.0089, 0.0084, 0.0079, 0.0073, 0.0066,
     . 0.0058, 0.0053, 0.0046, 0.0040, 0.0035, 0.0029, 0.0022, 0.0017,
     . 0.0011, 0.0005,-0.0001,-0.0006,-0.0012,-0.0017,-0.0023,-0.0030,
     .-0.0034,-0.0039,-0.0045,-0.0051,-0.0057,-0.0064,-0.0070,-0.0075,
     .-0.0081,-0.0087,-0.0093,-0.0099,-0.0104,-0.0110,-0.0116,-0.0123,
     .-0.0129,-0.0136,-0.0141,-0.0147,-0.0153,-0.0158,-0.0164,-0.0169,
     . 0.0208, 0.0202, 0.0196, 0.0190, 0.0184, 0.0178, 0.0173, 0.0166,
     . 0.0160, 0.0155, 0.0149, 0.0144, 0.0138, 0.0131, 0.0125, 0.0119,
     . 0.0113, 0.0107, 0.0099, 0.0096, 0.0091, 0.0085, 0.0079, 0.0072,
     . 0.0064, 0.0059, 0.0052, 0.0046, 0.0041, 0.0035, 0.0029, 0.0023,
     . 0.0018, 0.0012, 0.0006, 0.0001,-0.0005,-0.0010,-0.0016,-0.0024,
     .-0.0028,-0.0033,-0.0039,-0.0045,-0.0052,-0.0059,-0.0064,-0.0069,
     .-0.0075,-0.0081,-0.0087,-0.0093,-0.0099,-0.0105,-0.0111,-0.0118,
     .-0.0124,-0.0130,-0.0135,-0.0141,-0.0147,-0.0152,-0.0158,-0.0164,
     . 0.0216, 0.0209, 0.0203, 0.0196, 0.0191, 0.0184, 0.0179, 0.0172,
     . 0.0166, 0.0161, 0.0155, 0.0149, 0.0143, 0.0136, 0.0130, 0.0125,
     . 0.0119, 0.0113, 0.0105, 0.0102, 0.0096, 0.0091, 0.0085, 0.0078,
     . 0.0069, 0.0064, 0.0057, 0.0052, 0.0046, 0.0041, 0.0035, 0.0030,
     . 0.0024, 0.0018, 0.0012, 0.0006, 0.0001,-0.0005,-0.0011,-0.0018,
     .-0.0023,-0.0029,-0.0034,-0.0040,-0.0046,-0.0054,-0.0059,-0.0065,
     .-0.0069,-0.0075,-0.0081,-0.0087,-0.0093,-0.0099,-0.0105,-0.0111,
     .-0.0118,-0.0123,-0.0129,-0.0135,-0.0141,-0.0146,-0.0151,-0.0157,
     . 0.0222, 0.0216, 0.0210, 0.0203, 0.0197, 0.0190, 0.0185, 0.0177,
     . 0.0172, 0.0167, 0.0160, 0.0154, 0.0148, 0.0141, 0.0136, 0.0131,
     . 0.0125, 0.0119, 0.0111, 0.0107, 0.0102, 0.0096, 0.0089, 0.0082,
     . 0.0074, 0.0069, 0.0063, 0.0057, 0.0052, 0.0047, 0.0041, 0.0036,
     . 0.0030, 0.0025, 0.0019, 0.0012, 0.0006, 0.0001,-0.0005,-0.0013,
     .-0.0017,-0.0023,-0.0028,-0.0034,-0.0040,-0.0048,-0.0052,-0.0059,
     .-0.0063,-0.0069,-0.0075,-0.0081,-0.0087,-0.0093,-0.0099,-0.0105,
     .-0.0111,-0.0117,-0.0122,-0.0128,-0.0135,-0.0139,-0.0145,-0.0150,
     . 0.0229, 0.0222, 0.0216, 0.0209, 0.0203, 0.0196, 0.0190, 0.0183,
     . 0.0178, 0.0173, 0.0166, 0.0160, 0.0154, 0.0147, 0.0142, 0.0136,
     . 0.0131, 0.0125, 0.0117, 0.0114, 0.0108, 0.0102, 0.0095, 0.0088,
     . 0.0080, 0.0075, 0.0069, 0.0064, 0.0058, 0.0054, 0.0048, 0.0042,
     . 0.0038, 0.0031, 0.0025, 0.0019, 0.0012, 0.0006, 0.0001,-0.0007,
     .-0.0012,-0.0017,-0.0022,-0.0028,-0.0034,-0.0042,-0.0046,-0.0052,
     .-0.0057,-0.0063,-0.0069,-0.0075,-0.0081,-0.0088,-0.0093,-0.0099,
     .-0.0105,-0.0111,-0.0116,-0.0122,-0.0129,-0.0133,-0.0139,-0.0144,
     . 0.0235, 0.0229, 0.0222, 0.0215, 0.0209, 0.0203, 0.0196, 0.0190,
     . 0.0184, 0.0179, 0.0173, 0.0166, 0.0160, 0.0153, 0.0148, 0.0142,
     . 0.0137, 0.0131, 0.0123, 0.0119, 0.0114, 0.0108, 0.0102, 0.0095,
     . 0.0086, 0.0082, 0.0076, 0.0071, 0.0065, 0.0060, 0.0054, 0.0049,
     . 0.0043, 0.0037, 0.0031, 0.0024, 0.0018, 0.0012, 0.0006,-0.0001,
     .-0.0006,-0.0011,-0.0016,-0.0022,-0.0028,-0.0035,-0.0040,-0.0045,
     .-0.0051,-0.0057,-0.0063,-0.0069,-0.0075,-0.0081,-0.0087,-0.0093,
     .-0.0099,-0.0105,-0.0111,-0.0116,-0.0123,-0.0127,-0.0133,-0.0138,
     . 0.0242, 0.0235, 0.0228, 0.0221, 0.0215, 0.0208, 0.0202, 0.0196,
     . 0.0190, 0.0186, 0.0180, 0.0174, 0.0166, 0.0160, 0.0154, 0.0149,
     . 0.0143, 0.0137, 0.0129, 0.0125, 0.0120, 0.0114, 0.0108, 0.0101,
     . 0.0093, 0.0089, 0.0083, 0.0077, 0.0072, 0.0066, 0.0060, 0.0055,
     . 0.0049, 0.0043, 0.0037, 0.0030, 0.0024, 0.0018, 0.0012, 0.0005,
     . 0.0001,-0.0005,-0.0009,-0.0015,-0.0021,-0.0028,-0.0034,-0.0039,
     .-0.0045,-0.0052,-0.0058,-0.0063,-0.0070,-0.0075,-0.0081,-0.0087,
     .-0.0093,-0.0099,-0.0105,-0.0110,-0.0117,-0.0122,-0.0127,-0.0131,
     . 0.0248, 0.0241, 0.0234, 0.0227, 0.0221, 0.0215, 0.0208, 0.0202,
     . 0.0197, 0.0192, 0.0186, 0.0180, 0.0173, 0.0166, 0.0160, 0.0155,
     . 0.0149, 0.0142, 0.0135, 0.0131, 0.0126, 0.0120, 0.0114, 0.0108,
     . 0.0100, 0.0095, 0.0089, 0.0083, 0.0078, 0.0072, 0.0066, 0.0060,
     . 0.0054, 0.0048, 0.0042, 0.0036, 0.0030, 0.0024, 0.0018, 0.0011,
     . 0.0007, 0.0002,-0.0003,-0.0009,-0.0015,-0.0022,-0.0028,-0.0034,
     .-0.0039,-0.0046,-0.0052,-0.0057,-0.0063,-0.0069,-0.0075,-0.0081,
     .-0.0087,-0.0093,-0.0099,-0.0104,-0.0111,-0.0114,-0.0119,-0.0125,
     . 0.0255, 0.0248, 0.0241, 0.0234, 0.0227, 0.0220, 0.0214, 0.0207,
     . 0.0203, 0.0198, 0.0193, 0.0186, 0.0179, 0.0172, 0.0166, 0.0161,
     . 0.0155, 0.0148, 0.0141, 0.0137, 0.0132, 0.0127, 0.0121, 0.0114,
     . 0.0107, 0.0102, 0.0096, 0.0090, 0.0083, 0.0077, 0.0071, 0.0065,
     . 0.0060, 0.0054, 0.0048, 0.0043, 0.0036, 0.0031, 0.0025, 0.0018,
     . 0.0014, 0.0009, 0.0004,-0.0002,-0.0009,-0.0016,-0.0021,-0.0028,
     .-0.0033,-0.0040,-0.0046,-0.0051,-0.0057,-0.0063,-0.0069,-0.0075,
     .-0.0081,-0.0087,-0.0093,-0.0099,-0.0105,-0.0109,-0.0114,-0.0118,
     . 0.0260, 0.0254, 0.0247, 0.0240, 0.0233, 0.0226, 0.0220, 0.0213,
     . 0.0208, 0.0204, 0.0199, 0.0193, 0.0186, 0.0178, 0.0173, 0.0167,
     . 0.0161, 0.0154, 0.0147, 0.0143, 0.0139, 0.0133, 0.0127, 0.0121,
     . 0.0113, 0.0108, 0.0102, 0.0096, 0.0089, 0.0083, 0.0077, 0.0071,
     . 0.0066, 0.0060, 0.0054, 0.0049, 0.0043, 0.0037, 0.0032, 0.0024,
     . 0.0020, 0.0015, 0.0009, 0.0003,-0.0003,-0.0011,-0.0015,-0.0022,
     .-0.0028,-0.0034,-0.0040,-0.0045,-0.0050,-0.0057,-0.0063,-0.0069,
     .-0.0075,-0.0081,-0.0088,-0.0093,-0.0099,-0.0103,-0.0108,-0.0113,
     . 0.0266, 0.0260, 0.0253, 0.0246, 0.0240, 0.0233, 0.0226, 0.0220,
     . 0.0215, 0.0210, 0.0205, 0.0199, 0.0192, 0.0185, 0.0179, 0.0173,
     . 0.0167, 0.0161, 0.0154, 0.0150, 0.0145, 0.0140, 0.0134, 0.0127,
     . 0.0119, 0.0115, 0.0108, 0.0102, 0.0095, 0.0089, 0.0083, 0.0077,
     . 0.0071, 0.0066, 0.0061, 0.0055, 0.0049, 0.0043, 0.0037, 0.0030,
     . 0.0026, 0.0021, 0.0015, 0.0009, 0.0002,-0.0006,-0.0011,-0.0017,
     .-0.0022,-0.0028,-0.0034,-0.0039,-0.0045,-0.0051,-0.0057,-0.0063,
     .-0.0070,-0.0076,-0.0082,-0.0088,-0.0094,-0.0098,-0.0103,-0.0107,
     . 0.0272, 0.0265, 0.0259, 0.0252, 0.0246, 0.0239, 0.0233, 0.0226,
     . 0.0221, 0.0217, 0.0211, 0.0205, 0.0198, 0.0191, 0.0185, 0.0179,
     . 0.0173, 0.0168, 0.0160, 0.0157, 0.0152, 0.0146, 0.0140, 0.0133,
     . 0.0125, 0.0120, 0.0113, 0.0106, 0.0100, 0.0094, 0.0088, 0.0083,
     . 0.0077, 0.0072, 0.0066, 0.0061, 0.0055, 0.0048, 0.0042, 0.0036,
     . 0.0031, 0.0026, 0.0020, 0.0014, 0.0008, 0.0000,-0.0005,-0.0010,
     .-0.0016,-0.0022,-0.0028,-0.0033,-0.0039,-0.0045,-0.0051,-0.0058,
     .-0.0064,-0.0070,-0.0076,-0.0081,-0.0088,-0.0092,-0.0097,-0.0102,
     . 0.0278, 0.0272, 0.0265, 0.0259, 0.0253, 0.0246, 0.0240, 0.0233,
     . 0.0227, 0.0223, 0.0218, 0.0211, 0.0204, 0.0197, 0.0191, 0.0186,
     . 0.0180, 0.0175, 0.0168, 0.0164, 0.0159, 0.0153, 0.0146, 0.0139,
     . 0.0131, 0.0125, 0.0118, 0.0112, 0.0106, 0.0101, 0.0095, 0.0090,
     . 0.0084, 0.0079, 0.0073, 0.0067, 0.0061, 0.0055, 0.0049, 0.0041,
     . 0.0037, 0.0032, 0.0027, 0.0021, 0.0015, 0.0007, 0.0002,-0.0003,
     .-0.0009,-0.0016,-0.0021,-0.0026,-0.0032,-0.0038,-0.0044,-0.0050,
     .-0.0057,-0.0063,-0.0069,-0.0075,-0.0081,-0.0084,-0.0090,-0.0095,
     . 0.0283, 0.0276, 0.0270, 0.0264, 0.0257, 0.0251, 0.0245, 0.0238,
     . 0.0233, 0.0228, 0.0222, 0.0216, 0.0209, 0.0202, 0.0196, 0.0191,
     . 0.0186, 0.0180, 0.0174, 0.0170, 0.0164, 0.0158, 0.0152, 0.0144,
     . 0.0136, 0.0130, 0.0124, 0.0118, 0.0112, 0.0107, 0.0102, 0.0096,
     . 0.0091, 0.0084, 0.0078, 0.0072, 0.0066, 0.0060, 0.0054, 0.0047,
     . 0.0043, 0.0038, 0.0033, 0.0028, 0.0022, 0.0015, 0.0010, 0.0004,
     .-0.0002,-0.0009,-0.0014,-0.0019,-0.0025,-0.0031,-0.0037,-0.0043,
     .-0.0050,-0.0056,-0.0062,-0.0067,-0.0073,-0.0078,-0.0083,-0.0088,
     . 0.0287, 0.0281, 0.0274, 0.0268, 0.0262, 0.0255, 0.0250, 0.0243,
     . 0.0237, 0.0233, 0.0227, 0.0221, 0.0215, 0.0208, 0.0203, 0.0198,
     . 0.0193, 0.0187, 0.0180, 0.0176, 0.0170, 0.0164, 0.0157, 0.0150,
     . 0.0141, 0.0136, 0.0130, 0.0125, 0.0119, 0.0114, 0.0109, 0.0103,
     . 0.0098, 0.0091, 0.0084, 0.0078, 0.0072, 0.0066, 0.0060, 0.0054,
     . 0.0049, 0.0044, 0.0040, 0.0034, 0.0029, 0.0022, 0.0017, 0.0011,
     . 0.0005,-0.0002,-0.0008,-0.0013,-0.0018,-0.0024,-0.0030,-0.0037,
     .-0.0044,-0.0050,-0.0056,-0.0061,-0.0067,-0.0072,-0.0077,-0.0082,
     . 0.0293, 0.0286, 0.0280, 0.0274, 0.0268, 0.0262, 0.0256, 0.0249,
     . 0.0244, 0.0239, 0.0233, 0.0227, 0.0221, 0.0215, 0.0210, 0.0205,
     . 0.0200, 0.0194, 0.0186, 0.0182, 0.0177, 0.0170, 0.0163, 0.0156,
     . 0.0148, 0.0144, 0.0138, 0.0133, 0.0127, 0.0122, 0.0117, 0.0111,
     . 0.0105, 0.0098, 0.0091, 0.0084, 0.0078, 0.0072, 0.0067, 0.0060,
     . 0.0056, 0.0051, 0.0045, 0.0040, 0.0035, 0.0028, 0.0023, 0.0017,
     . 0.0011, 0.0005,-0.0001,-0.0006,-0.0012,-0.0018,-0.0024,-0.0031,
     .-0.0037,-0.0044,-0.0050,-0.0056,-0.0061,-0.0066,-0.0072,-0.0077,
     . 0.0300, 0.0294, 0.0288, 0.0282, 0.0276, 0.0269, 0.0263, 0.0257,
     . 0.0251, 0.0246, 0.0241, 0.0235, 0.0229, 0.0223, 0.0218, 0.0213,
     . 0.0207, 0.0201, 0.0193, 0.0188, 0.0182, 0.0176, 0.0170, 0.0164,
     . 0.0156, 0.0152, 0.0146, 0.0141, 0.0136, 0.0130, 0.0124, 0.0118,
     . 0.0111, 0.0105, 0.0097, 0.0091, 0.0085, 0.0079, 0.0074, 0.0067,
     . 0.0062, 0.0057, 0.0051, 0.0046, 0.0040, 0.0034, 0.0029, 0.0023,
     . 0.0017, 0.0011, 0.0005, 0.0000,-0.0006,-0.0012,-0.0018,-0.0025,
     .-0.0033,-0.0038,-0.0045,-0.0050,-0.0056,-0.0061,-0.0067,-0.0072,
     . 0.0307, 0.0301, 0.0295, 0.0289, 0.0283, 0.0277, 0.0271, 0.0264,
     . 0.0258, 0.0254, 0.0248, 0.0243, 0.0237, 0.0231, 0.0225, 0.0221,
     . 0.0214, 0.0207, 0.0198, 0.0194, 0.0188, 0.0182, 0.0177, 0.0171,
     . 0.0164, 0.0160, 0.0153, 0.0148, 0.0142, 0.0136, 0.0131, 0.0124,
     . 0.0117, 0.0111, 0.0103, 0.0097, 0.0091, 0.0086, 0.0080, 0.0073,
     . 0.0068, 0.0063, 0.0057, 0.0051, 0.0046, 0.0039, 0.0034, 0.0028,
     . 0.0023, 0.0017, 0.0011, 0.0006, 0.0000,-0.0005,-0.0012,-0.0020,
     .-0.0027,-0.0033,-0.0039,-0.0045,-0.0051,-0.0056,-0.0062,-0.0068,
     . 0.0316, 0.0310, 0.0303, 0.0297, 0.0291, 0.0285, 0.0278, 0.0271,
     . 0.0266, 0.0262, 0.0256, 0.0251, 0.0245, 0.0238, 0.0233, 0.0227,
     . 0.0220, 0.0213, 0.0204, 0.0200, 0.0195, 0.0189, 0.0184, 0.0179,
     . 0.0171, 0.0167, 0.0160, 0.0154, 0.0147, 0.0142, 0.0137, 0.0130,
     . 0.0124, 0.0117, 0.0110, 0.0104, 0.0098, 0.0093, 0.0087, 0.0080,
     . 0.0075, 0.0069, 0.0063, 0.0057, 0.0051, 0.0045, 0.0040, 0.0034,
     . 0.0029, 0.0023, 0.0017, 0.0013, 0.0007, 0.0001,-0.0006,-0.0014,
     .-0.0021,-0.0026,-0.0033,-0.0039,-0.0045,-0.0050,-0.0056,-0.0062,
     . 0.0322, 0.0316, 0.0309, 0.0302, 0.0296, 0.0289, 0.0284, 0.0277,
     . 0.0272, 0.0268, 0.0263, 0.0257, 0.0251, 0.0244, 0.0238, 0.0232,
     . 0.0225, 0.0217, 0.0209, 0.0205, 0.0200, 0.0195, 0.0190, 0.0185,
     . 0.0177, 0.0172, 0.0165, 0.0159, 0.0152, 0.0146, 0.0141, 0.0135,
     . 0.0129, 0.0123, 0.0116, 0.0111, 0.0105, 0.0099, 0.0094, 0.0086,
     . 0.0081, 0.0075, 0.0069, 0.0063, 0.0058, 0.0051, 0.0046, 0.0041,
     . 0.0036, 0.0029, 0.0023, 0.0018, 0.0012, 0.0007, 0.0000,-0.0008,
     .-0.0015,-0.0021,-0.0027,-0.0033,-0.0039,-0.0044,-0.0050,-0.0055,
     . 0.0327, 0.0321, 0.0314, 0.0307, 0.0301, 0.0294, 0.0288, 0.0282,
     . 0.0277, 0.0273, 0.0268, 0.0263, 0.0257, 0.0250, 0.0243, 0.0237,
     . 0.0229, 0.0223, 0.0215, 0.0211, 0.0206, 0.0202, 0.0197, 0.0191,
     . 0.0183, 0.0177, 0.0170, 0.0163, 0.0156, 0.0151, 0.0146, 0.0140,
     . 0.0136, 0.0129, 0.0122, 0.0116, 0.0110, 0.0105, 0.0100, 0.0092,
     . 0.0087, 0.0081, 0.0075, 0.0070, 0.0064, 0.0057, 0.0053, 0.0047,
     . 0.0041, 0.0035, 0.0028, 0.0024, 0.0018, 0.0012, 0.0005,-0.0002,
     .-0.0009,-0.0015,-0.0022,-0.0027,-0.0033,-0.0038,-0.0044,-0.0050,
     . 0.0333, 0.0327, 0.0320, 0.0313, 0.0306, 0.0300, 0.0294, 0.0289,
     . 0.0284, 0.0280, 0.0275, 0.0269, 0.0262, 0.0255, 0.0248, 0.0242,
     . 0.0235, 0.0228, 0.0221, 0.0218, 0.0214, 0.0209, 0.0203, 0.0197,
     . 0.0188, 0.0182, 0.0175, 0.0168, 0.0161, 0.0156, 0.0152, 0.0147,
     . 0.0142, 0.0135, 0.0128, 0.0122, 0.0116, 0.0111, 0.0105, 0.0098,
     . 0.0093, 0.0087, 0.0081, 0.0076, 0.0070, 0.0064, 0.0059, 0.0053,
     . 0.0047, 0.0040, 0.0034, 0.0030, 0.0024, 0.0018, 0.0010, 0.0004,
     .-0.0003,-0.0010,-0.0016,-0.0022,-0.0027,-0.0033,-0.0039,-0.0044,
     . 0.0339, 0.0332, 0.0326, 0.0319, 0.0313, 0.0306, 0.0300, 0.0295,
     . 0.0290, 0.0287, 0.0281, 0.0275, 0.0268, 0.0260, 0.0253, 0.0248,
     . 0.0241, 0.0235, 0.0228, 0.0225, 0.0220, 0.0216, 0.0209, 0.0202,
     . 0.0193, 0.0187, 0.0180, 0.0174, 0.0167, 0.0162, 0.0158, 0.0154,
     . 0.0148, 0.0142, 0.0134, 0.0128, 0.0123, 0.0117, 0.0111, 0.0105,
     . 0.0099, 0.0093, 0.0088, 0.0082, 0.0076, 0.0070, 0.0064, 0.0058,
     . 0.0053, 0.0046, 0.0040, 0.0035, 0.0029, 0.0023, 0.0016, 0.0009,
     . 0.0002,-0.0003,-0.0010,-0.0015,-0.0021,-0.0027,-0.0033,-0.0039,
     . 0.0344, 0.0338, 0.0331, 0.0324, 0.0319, 0.0313, 0.0307, 0.0302,
     . 0.0296, 0.0293, 0.0287, 0.0281, 0.0272, 0.0265, 0.0258, 0.0253,
     . 0.0247, 0.0242, 0.0235, 0.0232, 0.0227, 0.0221, 0.0214, 0.0207,
     . 0.0199, 0.0193, 0.0186, 0.0180, 0.0174, 0.0169, 0.0165, 0.0160,
     . 0.0155, 0.0148, 0.0141, 0.0134, 0.0129, 0.0123, 0.0117, 0.0110,
     . 0.0105, 0.0099, 0.0094, 0.0089, 0.0082, 0.0076, 0.0071, 0.0064,
     . 0.0059, 0.0052, 0.0046, 0.0041, 0.0035, 0.0028, 0.0022, 0.0015,
     . 0.0008, 0.0002,-0.0004,-0.0009,-0.0015,-0.0020,-0.0026,-0.0032,
     . 0.0349, 0.0343, 0.0337, 0.0330, 0.0325, 0.0320, 0.0314, 0.0308,
     . 0.0303, 0.0298, 0.0292, 0.0285, 0.0277, 0.0270, 0.0264, 0.0259,
     . 0.0253, 0.0248, 0.0241, 0.0238, 0.0232, 0.0227, 0.0219, 0.0212,
     . 0.0203, 0.0198, 0.0192, 0.0186, 0.0181, 0.0176, 0.0171, 0.0167,
     . 0.0161, 0.0154, 0.0147, 0.0141, 0.0135, 0.0129, 0.0123, 0.0117,
     . 0.0112, 0.0106, 0.0100, 0.0094, 0.0088, 0.0082, 0.0076, 0.0070,
     . 0.0064, 0.0058, 0.0051, 0.0046, 0.0040, 0.0035, 0.0028, 0.0021,
     . 0.0014, 0.0008, 0.0002,-0.0003,-0.0010,-0.0015,-0.0020,-0.0025,
     . 0.0353, 0.0348, 0.0342, 0.0336, 0.0331, 0.0326, 0.0320, 0.0314,
     . 0.0308, 0.0303, 0.0297, 0.0290, 0.0282, 0.0275, 0.0270, 0.0266,
     . 0.0260, 0.0255, 0.0248, 0.0244, 0.0238, 0.0231, 0.0224, 0.0217,
     . 0.0209, 0.0204, 0.0198, 0.0193, 0.0187, 0.0183, 0.0178, 0.0173,
     . 0.0167, 0.0161, 0.0153, 0.0147, 0.0141, 0.0135, 0.0130, 0.0123,
     . 0.0117, 0.0112, 0.0106, 0.0101, 0.0094, 0.0087, 0.0082, 0.0076,
     . 0.0070, 0.0063, 0.0057, 0.0053, 0.0046, 0.0040, 0.0034, 0.0027,
     . 0.0020, 0.0015, 0.0008, 0.0003,-0.0003,-0.0009,-0.0014,-0.0019,
     . 0.0359, 0.0353, 0.0348, 0.0343, 0.0338, 0.0332, 0.0326, 0.0320,
     . 0.0313, 0.0309, 0.0302, 0.0295, 0.0288, 0.0281, 0.0276, 0.0273,
     . 0.0267, 0.0262, 0.0255, 0.0250, 0.0243, 0.0236, 0.0229, 0.0222,
     . 0.0214, 0.0210, 0.0205, 0.0200, 0.0194, 0.0189, 0.0184, 0.0179,
     . 0.0173, 0.0166, 0.0158, 0.0152, 0.0147, 0.0141, 0.0135, 0.0128,
     . 0.0123, 0.0117, 0.0111, 0.0105, 0.0099, 0.0092, 0.0087, 0.0081,
     . 0.0075, 0.0068, 0.0062, 0.0058, 0.0052, 0.0046, 0.0039, 0.0032,
     . 0.0026, 0.0021, 0.0014, 0.0009, 0.0003,-0.0002,-0.0008,-0.0014,
     . 0.0364, 0.0359, 0.0354, 0.0348, 0.0344, 0.0338, 0.0332, 0.0325,
     . 0.0319, 0.0314, 0.0308, 0.0301, 0.0294, 0.0288, 0.0282, 0.0279,
     . 0.0273, 0.0267, 0.0259, 0.0254, 0.0248, 0.0241, 0.0235, 0.0228,
     . 0.0220, 0.0216, 0.0210, 0.0205, 0.0199, 0.0194, 0.0189, 0.0184,
     . 0.0178, 0.0171, 0.0164, 0.0158, 0.0153, 0.0147, 0.0141, 0.0134,
     . 0.0129, 0.0123, 0.0117, 0.0111, 0.0105, 0.0098, 0.0093, 0.0087,
     . 0.0081, 0.0074, 0.0069, 0.0065, 0.0059, 0.0053, 0.0046, 0.0039,
     . 0.0033, 0.0027, 0.0021, 0.0016, 0.0010, 0.0004,-0.0001,-0.0007,
     . 0.0369, 0.0364, 0.0359, 0.0353, 0.0349, 0.0344, 0.0337, 0.0330,
     . 0.0324, 0.0320, 0.0314, 0.0308, 0.0301, 0.0294, 0.0288, 0.0285,
     . 0.0279, 0.0273, 0.0264, 0.0260, 0.0253, 0.0247, 0.0240, 0.0233,
     . 0.0226, 0.0221, 0.0216, 0.0210, 0.0204, 0.0199, 0.0194, 0.0189,
     . 0.0183, 0.0177, 0.0169, 0.0164, 0.0158, 0.0152, 0.0147, 0.0140,
     . 0.0134, 0.0128, 0.0122, 0.0116, 0.0110, 0.0104, 0.0098, 0.0093,
     . 0.0087, 0.0080, 0.0075, 0.0071, 0.0065, 0.0059, 0.0053, 0.0046,
     . 0.0039, 0.0034, 0.0028, 0.0023, 0.0017, 0.0011, 0.0005,-0.0001,
     . 0.0375, 0.0370, 0.0365, 0.0359, 0.0355, 0.0350, 0.0342, 0.0336,
     . 0.0330, 0.0326, 0.0320, 0.0314, 0.0306, 0.0300, 0.0295, 0.0291,
     . 0.0284, 0.0278, 0.0270, 0.0265, 0.0258, 0.0252, 0.0245, 0.0239,
     . 0.0231, 0.0227, 0.0221, 0.0215, 0.0209, 0.0204, 0.0199, 0.0194,
     . 0.0188, 0.0182, 0.0175, 0.0169, 0.0164, 0.0158, 0.0152, 0.0145,
     . 0.0139, 0.0133, 0.0127, 0.0122, 0.0115, 0.0109, 0.0104, 0.0098,
     . 0.0092, 0.0087, 0.0081, 0.0077, 0.0071, 0.0066, 0.0058, 0.0053,
     . 0.0045, 0.0040, 0.0034, 0.0029, 0.0023, 0.0017, 0.0011, 0.0006/
      end
 
 
 
 
      subroutine sat_fix(pixr)
      implicit none
 
      real      pixr(4096,4096)
      real      pixs(4096,4096)
      real      pixx(4096,4096)
 
      integer*2 pixi(4096,4096)
      integer*4 pixii
      byte      pixu(4096,4096)
 
      integer i
      integer j
      integer NARGs, NARG, iargc
 
      integer L, LL, Ls
      integer iil(9999)
      integer jjl(9999)
      real    frl(9999)
      real    rrl(9999)
      integer nrl(9999)
      real    srl(9999)
 
      logical trip
      character*80 FILEI
      character*80 FILEO
 
      real UPR_LIM
 
      integer iu
      real pixe
      integer il(99999), jl(99999)
      integer NSATs
      integer NCENs
      real sr, mbar_sky
      real fr, fru, fr0
      real apphot_sat
      real apphot_satX
      integer nr, nr0
      integer ir
 
      real    reff
      integer ireff
      integer ii, jj
      integer ncirc
 
      NARGs = iargc()
 
 144  format(1x,i4,1x,i4,1x,i7,1x,i7)
 
      do i = 0001, 4096
         pixr(i,2046) = -750.0
         pixr(i,2047) = -750.0
         pixr(i,2048) = -750.0
         pixr(i,2049) = -750.0
         pixr(i,2050) = -750.0
         pixr(i,2051) = -750.0
         enddo
 
 
      print*,' '
      print*,'ASSUME GAIN OF TWO.'
      print*,' '
      UPR_LIM = 62500.0
 
      do i = 0001, 4096
      do j = 0001, 4096
         pixu(i,j) = 0
         pixs(i,j) = pixr(i,j)
         pixx(i,j) = pixr(i,j)
         if (abs(pixr(i,j)-99999).le.100) then
            if (pixr(i+1,j  ).eq.pixr(i-1,j  ).and.
     .          pixr(i  ,j+1).eq.pixr(i  ,j-1).and.
     .          pixr(i+1,j+1).eq.pixr(i-1,j-1).and.
     .          pixr(i+1,j-1).eq.pixr(i-1,j+1)) then
                print*,'SAT_FIX:  ALREADY FIXED!'
                return
                endif
            endif
         enddo
         enddo
 
      NSATs = 0
      NCENs = 0
      do i = 1, 4096
         if (i.eq.i/128*128) write(*,'(''  i: '',i4.4,1x,2i9)')
     .                        i,NSATs,NCENs
         do j = 1, 4096
            if (pixr(i,j).gt.UPR_LIM.and.pixu(i,j).eq.0) then
               call find_sat(i,j,pixr,pixu,UPR_LIM,il,jl,Ls)
               endif
            if (pixr(i,j).ge.99998) NSATs = NSATs + 1
            if (pixr(i,j).eq.99999) NCENs = NCENs + 1
            enddo
         enddo
 
      do i = 0001, 4096
      do j = 0001, 4096
         if (pixr(i,j).ge.99998) pixx(i,j) = -750
         enddo
         enddo
 
      Ls = 0
      do i = 0001+10, 4096-10
      do j = 0001+10, 4096-10
         if (pixr(i,j).eq.99999) then
            sr = mbar_sky(i,j,15,20,pixs)
            fr  = apphot_sat(i,j,pixs,sr,nr)
            Ls = Ls + 1
            iil(Ls) = i
            jjl(Ls) = j
            frl(Ls) = fr
            nrl(Ls) = nr
            srl(Ls) = sr
            reff = sqrt(fr/90000/3.14159)
            rrl(Ls) = reff
            endif
         enddo
         enddo
 
      call cosort_riiirr(frl,iil,jjl,nrl,srl,rrl,Ls)
 
      do L = 1, Ls
         write( *,888) iil(L),jjl(L),
     .                 -2.5*log10(max(frl(L),1.0)),
     .                 srl(L),frl(L),nrl(L),L
c         write(98,888) iil(L),jjl(L),
c     .                 -2.5*log10(max(frl(L),1.0)),
c     .                 srl(L),frl(L),nrl(L),L
  888    format(i4.4,1x,i4.4,1x,f10.4,1x,f8.2,
     .               1x,f12.1,1x,i8,1x,i5.5)
         enddo
 
      do L = 1, Ls
         i  = iil(L)
         j  = jjl(L)
         fr = frl(L)
         reff = sqrt(fr/90000/3.14159)
         ireff = reff + 0.90
         if (ireff.lt.1) ireff = 1
         ncirc = 0
         do ii = -ireff, +ireff
         do jj = -ireff, +ireff
            if (ii**2+jj**2.le.(ireff+0.5)**2)
     .          ncirc = ncirc + 1
            enddo
            enddo
         fru = fr-100000
         do ii = -ireff-1, +ireff+1
         do jj = -ireff-1, +ireff+1
            if (ii**2+jj**2.le.(ireff+1.5)**2)
     .          pixx(i+ii,j+jj) = -750
            enddo
            enddo
         do LL = L+1, Ls
            if ((iil(L)-iil(LL))**2+
     .          (jjl(L)-jjl(LL))**2.lt.
     .          (rrl(L)+rrl(LL)+3)**2) goto 777
            enddo
         if (abs(j-2048.5).lt.reff+5) goto 777
         if (j+0001.lt.reff+5) goto 777
         if (4095-j.lt.reff+5) goto 777
         do ii = -ireff, +ireff
         do jj = -ireff, +ireff
            if (ii**2+jj**2.le.(ireff+0.5)**2)
     .          pixx(i+ii,j+jj) = fru/(ncirc-1)
            enddo
            enddo
         pixx(i,j) = 100000
  777    continue
         enddo
 
 
      if (.false.) then
         open(55,file='LOG.peaksat_CEN.reg',status='unknown')
         write(55,'(''# Region file format: DS9 version 3.0'')')
         write(55,'(''global color=yellow'')')
 
         open(56,file='LOG.peaksat_ALL.reg',status='unknown')
         write(56,'(''# Region file format: DS9 version 3.0'')')
         write(56,'(''global color=yellow'')')
         do i = 1, 4096
         do j = 1, 4096
            if (pixr(i,j).eq.99999) write(55,112) i,j,1.0,'red'
            if (pixr(i,j).ge.99998) write(56,112) i,j,0.5,'yellow'
            enddo
            enddo
         close(55)
         close(56)
  112    format('image;circle( ',i4,',',i4,',',f6.3,') # color=',a7)
         endif
 
      do i = 1, 4096
      do j = 1, 4096
         pixr(i,j) = pixx(i,j)
         enddo
         enddo
 
 
c      if (.true.) then
c         print*,'---> re-anal: '
c         do i = 0001+10, 4096-10
c         do j = 0001+10, 4096-10
c            if (pixx(i,j).gt.95000) then
c               ir = 1
c    3          ir = ir + 1
c               fr = 0
c               trip = .false.
c               do ii = -ir, ir
c               do jj = -ir, ir
c                  if (ii**2+jj**2.le.(ir+0.5)**2) then
c                     if (pixx(i+ii,j+jj).gt.0)
c     .                   fr = fr + pixx(i+ii,j+jj)
c                     if (pixx(i+ii,j+jj).lt.0) trip = .true.
c                     endif
c                  enddo
c                  enddo
c               print*,i,j,ir,fr,trip
c               if (.not.trip) goto 3
c               write(99,120) i,j,-2.5*log10(max(fr,1))
c  120          format(i4.4,1x,i4.4,1x,f10.4)
c               endif
c            enddo
c            enddo
c         endif
 
 
      return
      end
 
 
 
 
      subroutine find_sat(i0,j0,pixr,pixu,UPR_LIM,il,jl,Ls)
      implicit none
 
      integer i0
      integer j0
      real pixr(4096,4096)
      byte pixu(4096,4096)
      real UPR_LIM
      integer il(99999)
      integer jl(99999)
      integer ml(99999)
      integer Ls
 
      integer Lo, l
      integer jmin, jmax
      integer ii,jj
      integer iii,jjj
      integer oi, oj, m
 
      integer imax, mmax
 
      Ls = 0
      if (pixr(i0,j0).lt.UPR_LIM) return
 
      Ls = 1
      il(1) = i0
      jl(1) = j0
      pixu(i0,j0) = 1
      jmin = 0001
      jmax = 2048
      if (j0.ge.2049) then
         jmin = 2049
         jmax = 4096
         endif
    1 continue
      Lo = Ls
      do l = 1, Lo
         do ii = max(il(L)-1,0001), min(il(L)+1,4096)
         do jj = max(jl(L)-1,jmin), min(jl(L)+1,jmax)
            if (pixu(ii,jj).eq.0) then
               if (pixr(ii,jj).gt.55000) then
                  pixu(ii,jj) = 1
                  Ls = Ls + 1
                  if (Ls.gt.99999) then
                     print*,'  Ls wants to be > 999999... '
                     endif
                  il(Ls) = ii
                  jl(Ls) = jj
                  m = 0
    3             continue
                  m = m + 1
                  iii = ii + oi(m)
                  jjj = jj + oj(m)
                  if (iii.ge.0001.and.iii.le.4096.and.
     .                jjj.ge.jmin.and.jjj.le.jmax) then
                      if (pixr(iii,jjj).gt.55000) goto 3
                      endif
                  ml(Ls) = m
                  endif
               endif
            enddo
            enddo
         enddo
      if (Ls.gt.Lo) goto 1
 
      imax = il(1)
      jmax = jl(1)
      mmax = ml(1)
      do L = 1, Ls
         if (ml(L).gt.mmax) then
            imax = il(L)
            jmax = jl(L)
            mmax = ml(L)
            endif
         enddo
      do L = 1, Ls
         pixr(il(L),jl(L)) = 99998
         enddo
      pixr(imax,jmax) = 99999
 
      return
      end
 
 
 
c--------------------------------------------------
c
c This routine orders the pixels by successive distance
c from each other.  If you step through the array
c like this:
c
c
c     do m = 1, 99
c        i = ixc + oi(m)
c        j = iyc + oj(m)
c        ...
c        enddo
c
c
c Then, the pixels (i,j) will be in order of increasing
c distance from the central pixel (ixc,iyc).
c
c This can be useful for many applications.
c
c
      integer function oi(m)
      implicit none
      integer m
 
      integer qm, fm
      integer im
      integer mm
      integer oim, ojm
      integer nm
 
      integer oil(91)
      data oil / 0,1,1,2,2,1,2,3,3,1,3,2,4,4,
     .               1,3,4,2,3,5,4,5,1,5,2,4,
     .           3,5,6,1,6,6,2,5,4,3,6,7,5,1,
     .               7,6,4,2,7,7,3,6,5,8,1,4,
     .           8,7,8,2,6,8,3,7,5,4,8,9,1,9,
     .                 9,7,2,6,5,8,3,9,4,9,7,
     .           9,5,8,6,8,7,8,9,6,9,7,9,8,9 /
 
      integer ojl(91)
      data ojl / 0,0,1,0,1,2,2,0,1,3,2,3,0,1,4,
     .                   3,2,4,4,0,3,1,5,2,5,4,
     .           5,3,0,6,1,2,6,4,5,6,3,0,5,7,1,
     .                   4,6,7,2,3,7,5,6,0,8,7,
     .           1,4,2,8,6,3,8,5,7,8,4,0,9,1,2,
     .                     6,9,7,8,5,9,3,9,4,7,
     .             5,9,6,8,7,8,8,6,9,7,9,8,9,9 /
 
 
      qm = m+2 - 4*int((m+2)/4)
      fm = int((m+6)/4)
      if (fm.le.91) then
         if (qm.eq.0) oi =  oil(fm)
         if (qm.eq.1) oi = -oil(fm)
         if (qm.eq.2) oi = -ojl(fm)
         if (qm.eq.3) oi =  ojl(fm)
         return
         endif
 
      im = (sqrt(1.*m-1)-1)/2
      mm = (2*im+1)**2+1
      nm =  (m-mm)/4
      oim = nm-im
      ojm = im
 
      if (qm.eq.0) oi =  oim
      if (qm.eq.1) oi = -oim
      if (qm.eq.2) oi = -ojm
      if (qm.eq.3) oi =  ojm
 
      return
      end
 
c----------------------------------------------------------------
c
c the y-analog to oi(m)
c
      integer function oj(m)
      implicit none
      integer m
 
      integer qm, fm
      integer im, mm, nm, oim, ojm
 
      integer oil(91)
      data oil / 0,1,1,2,2,1,2,3,3,1,3,2,4,4,
     .               1,3,4,2,3,5,4,5,1,5,2,4,
     .           3,5,6,1,6,6,2,5,4,3,6,7,5,1,
     .               7,6,4,2,7,7,3,6,5,8,1,4,
     .           8,7,8,2,6,8,3,7,5,4,8,9,1,9,
     .                 9,7,2,6,5,8,3,9,4,9,7,
     .           9,5,8,6,8,7,8,9,6,9,7,9,8,9 /
 
      integer ojl(91)
      data ojl / 0,0,1,0,1,2,2,0,1,3,2,3,0,1,4,
     .                   3,2,4,4,0,3,1,5,2,5,4,
     .           5,3,0,6,1,2,6,4,5,6,3,0,5,7,1,
     .                   4,6,7,2,3,7,5,6,0,8,7,
     .           1,4,2,8,6,3,8,5,7,8,4,0,9,1,2,
     .                     6,9,7,8,5,9,3,9,4,7,
     .             5,9,6,8,7,8,8,6,9,7,9,8,9,9 /
 
      qm = m+2 - 4*int((m+2)/4)
      fm = int((m+6)/4)
      if (fm.le.91) then
         if (qm.eq.0) oj =  ojl(fm)
         if (qm.eq.1) oj = -ojl(fm)
         if (qm.eq.2) oj =  oil(fm)
         if (qm.eq.3) oj = -oil(fm)
         return
         endif
 
      im = (sqrt(1.*m-1)-1)/2
      mm = (2*im+1)**2+1
      nm =  (m-mm)/4
      oim = nm-im
      ojm = im
 
      if (qm.eq.0) oj =  ojm
      if (qm.eq.1) oj = -ojm
      if (qm.eq.2) oj =  oim
      if (qm.eq.3) oj = -oim
 
      return
      end
 
c----------------------------------------------------------------
c
c should be obvious
c
      real function orm(m)
      implicit none
      integer m
 
      integer oi, oj
 
      orm = sqrt(1.*(oi(m)**2 + oj(m)**2))
 
      return
      end
 
 
 
c--------------------------------------------------------
c
c sorts in ASCENDING ORDER!!!
c
      subroutine cosort_riiirr(r1,i2,i3,i4,r5,r6,NTOT)
      implicit none
      real r1(1)
      integer i2(1), i3(1), i4(1)
      real r5(1), r6(1)
      integer NTOT
 
      integer n
      logical change
 777  continue
      change = .false.
      do n = 1, NTOT-1
         if (r1(n).gt.r1(n+1)) then
            call swap(r1(n),r1(n+1))
            call swip(i2(n),i2(n+1))
            call swip(i3(n),i3(n+1))
            call swip(i4(n),i4(n+1))
            call swap(r5(n),r5(n+1))
            call swap(r6(n),r6(n+1))
            change = .true.
            endif
         enddo
      if (change) goto 777
      end
 
 
      subroutine swap(x,y)
      implicit none
      real x, y, t
      t = x
      x = y
      y = t
      return
      end
 
 
      subroutine swip(x,y)
      implicit none
      integer x, y, t
      t = x
      x = y
      y = t
      return
      end
 
 
 
 
      real function apphot_sat(ii,jj,pixc,sr,nsat)
      implicit none
      integer ii,jj
      real pixc(4096,4096)
      real sr
      integer nsat
 
      integer u1(101,101)
      integer u2(101,101)
      integer nadd, naddo
      integer iu,ju
      integer iii,jjj
      integer nadd2
      real ftot, f5
      integer i,j
 
 
      nsat = 0
      ftot = 0
      do i = -5, 5
      do j = -5, 5
         if (i**2+j**2.lt.5.5**2) then
            ftot = ftot + (pixc(ii+i,jj+j)-sr)
            if (pixc(ii+i,jj+j).gt.70000) nsat = nsat + 1
            endif
         enddo
         enddo
      apphot_sat = ftot
      if (nsat.le.3) return
 
      do i = 001, 101
      do j = 001, 101
         u1(i,j) = 0
         u2(i,j) = 0
         enddo
         enddo
 
      do i = -2, 2
      do j = -2, 2
         if (pixc(ii+i,jj+j).gt.50000) u1(51+i,51+j) = 1
         enddo
         enddo
 
      naddo = 0
 
      nadd = 0
   1  continue
      naddo =  nadd
      do i = -49, 49
      do j = -49, 49
         if (i+ii.le.0002) goto 3
         if (j+jj.le.0002) goto 3
         if (i+ii.ge.4095) goto 3
         if (j+jj.ge.4095) goto 3
         iu = max(min(i+ii,4095),2)
         ju = max(min(j+jj,4095),2)
         if (u1(51+i,51+j).eq.1) then
            do iii = -1,+1
            do jjj = -1,+1
               if (pixc(iu+iii,ju+jjj).gt.50000) u1(51+i+iii,51+j+jjj)=1
               enddo
               enddo
            endif
    3    continue
         enddo
         enddo
 
      nadd = 0
      do i = 001, 101
      do j = 001, 101
         nadd = nadd + u1(i,j)
         enddo
         enddo
 
      if (nadd.gt.naddo) goto 1
 
 
      nadd2 = 0
      do i = 002, 100
      do j = 002, 100
         u2(i,j) = 0
         if (u1(i  ,j  ).eq.1.or.
     .       u1(i+1,j  ).eq.1.or.
     .       u1(i-1,j  ).eq.1.or.
     .       u1(i  ,j+1).eq.1.or.
     .       u1(i  ,j-1).eq.1.or.
     .       u1(i+1,j+1).eq.1.or.
     .       u1(i-1,j+1).eq.1.or.
     .       u1(i+1,j-1).eq.1.or.
     .       u1(i-1,j-1).eq.1) u2(i,j) = 1
         if ((i-51)**2+(j-51)**2.le.5.5**2) u2(i,j) = 1
         nadd2 = nadd2 + u2(i,j)
         enddo
         enddo
 
      f5   = 0
      ftot = 0
      nsat = 0
      do i = 001, 101
      do j = 001, 101
         if ((i-51)**2+(j-51)**2.lt.5.5**2) then
            f5 = f5 + (pixc(ii+(i-51),jj+(j-51))-sr)
            endif
         if (u2(i,j).eq.1.and.
     .      ii+(i-51).gt.0002.and.
     .      jj+(j-51).gt.0002.and.
     .      ii+(i-51).lt.4095.and.
     .      jj+(j-51).lt.4095) then
            ftot = ftot + (pixc(ii+(i-51),jj+(j-51))-sr)
            if (pixc(ii+(i-51),jj+(j-51)).gt.70000) nsat = nsat+1
            endif
         enddo
         enddo
 
      apphot_sat = ftot
      nsat = nadd2
 
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "SATAP/apphot_sat2.f"
c****                                         
c*********************************************
                                              
 
      real function apphot_sat2(ii,jj,pixc,sr,nsat,rtot,pixq)
      implicit none
      integer ii,jj
      real pixc(4096,4096)
      real sr
      integer nsat
      real    rtot
      real pixq(4096,4096)
 
      integer u1(101,1001)
      integer u2(101,1001)
      integer nadd, naddo
      integer iu,ju
      integer iii,jjj
      integer nadd2
      real ftot
      integer i,j
      real rpsf_photij
      integer hh(1000), hbin
      integer HIFLAG
      common / HIFLAG_ / HIFLAG
 
      if (.false.) then
         do hbin = 1, 1000
            hh(hbin) = 0
            enddo
         do i = 0001, 4096
         do j = 0001, 4096
            hbin = pixc(i,j)/100 + 0.5
            if (hbin.lt.0001) hbin = 0001
            if (hbin.gt.1000) hbin = 1000
            hh(hbin) = hh(hbin) + 1
            enddo
            enddo
         do hbin = 1, 1000
            write(91,'(i7,1x,i10)') hbin*100,hh(hbin)
            enddo
         stop
         endif
 
      nsat = 0
      ftot = 0
      rtot = 0
      do i = -2, 2
      do j = -2, 2
         if (i**2+j**2.lt.5.5**2) then
            ftot = ftot + (pixc(ii+i,jj+j)-sr)
            rtot = rtot + rpsf_photij(1.*i,1.*j,ii,jj)
            if (pixc(ii+i,jj+j).gt.HIFLAG) nsat = nsat + 1
            endif
         enddo
         enddo
      apphot_sat2 = ftot/rtot
      if (nsat.le.6) then
         write(99,199) 1,ii,jj,pixc(ii,jj),apphot_sat2,ftot,rtot,nsat
  199    format(1x,i1,1x,i4.4,1x,i4.4,1x,f8.1,
     .                1x,f10.1,1x,f10.1,1x,f8.5,1x,i5)
         pixq(ii,jj) = apphot_sat2
         return
         endif
 
 
      nsat = 0
      ftot = 0
      rtot = 0
      do i = -5, 5
      do j = -5, 5
         if (i**2+j**2.lt.5.5**2) then
            ftot = ftot + (pixc(ii+i,jj+j)-sr)
            rtot = rtot + rpsf_photij(1.*i,1.*j,ii,jj)
            if (pixc(ii+i,jj+j).gt.HIFLAG) nsat = nsat + 1
            endif
         enddo
         enddo
      apphot_sat2 = ftot/rtot
      if (nsat.le.10) then
         write(99,199) 2,ii,jj,pixc(ii,jj),apphot_sat2,ftot,rtot,nsat
         pixq(ii,jj) = apphot_sat2
         return
         endif
 
      do i = 001, 0101
      do j = 001, 1001
         u1(i,j) = 0
         u2(i,j) = 0
         enddo
         enddo
 
      do i = -2, 2
      do j = -2, 2
         if (pixc(ii+i,jj+j).gt.HIFLAG) u1(051+i,501+j) = 1
         enddo
         enddo
 
      naddo = 0
 
      nadd = 0
   1  continue
      naddo =  nadd
      do i = -049, 049
      do j = -499, 499
         if (i+ii.le.0002) goto 3
         if (j+jj.le.0002) goto 3
         if (jj.lt.2049.and.j+jj.gt.2046) goto 3
         if (jj.gt.2048.and.j+jj.lt.2051) goto 3
         if (i+ii.ge.4095) goto 3
         if (j+jj.ge.4095) goto 3
         iu = max(min(i+ii,4095),2)
         ju = max(min(j+jj,4095),2)
         if (u1(051+i,501+j).eq.1) then
            do iii = -1,+1
            do jjj = -1,+1
               if (pixc(iu+iii,ju+jjj).gt.HIFLAG)
     .             u1(051+i+iii,501+j+jjj)=1
               enddo
               enddo
            endif
    3    continue
         enddo
         enddo
 
      nadd = 0
      do i = 001, 0101
      do j = 001, 1001
         nadd = nadd + u1(i,j)
         enddo
         enddo
 
      if (nadd.gt.naddo) goto 1
 
 
      nadd2 = 0
      do i = 003, 0098
      do j = 003, 0998
         u2(i,j) = 0
         if (ii+(i-051).ge.0001.and.
     .       ii+(i-051).le.4096.and.
     .       jj+(j-501).ge.0001.and.
     .       jj+(j-501).le.4096) then
             if (pixc(ii+(i-051),jj+(j-501)).lt.HIFLAG) then
                if (u1(i  ,j  ).eq.1.or.
     .              u1(i+1,j  ).eq.1.or.
     .              u1(i-1,j  ).eq.1.or.
     .              u1(i  ,j+1).eq.1.or.
     .              u1(i  ,j-1).eq.1.or.
     .              u1(i+1,j+1).eq.1.or.
     .              u1(i-1,j+1).eq.1.or.
     .              u1(i+1,j-1).eq.1.or.
     .              u1(i-1,j-1).eq.1.or.
     .              u1(i-2,j  ).eq.1.or.
     .              u1(i+2,j  ).eq.1.or.
     .              u1(i  ,j-2).eq.1.or.
     .              u1(i  ,j+2).eq.1) u2(i,j) = 1
                    endif
                endif
         if ((i-051)**2+(j-501)**2.le.5.5**2) u2(i,j) = 1
         nadd2 = nadd2 + u2(i,j)
         enddo
         enddo
 
      ftot = 0
      rtot = 0
      nsat = 0
      do i = 001, 0101
      do j = 001, 1001
         if (u2(i,j).eq.1.and.
     .      ii+(i-051).gt.0002.and.
     .      jj+(j-501).gt.0002.and.
     .      ii+(i-051).lt.4095.and.
     .      jj+(j-501).lt.4095) then
            ftot = ftot + (pixc(ii+(i-51),jj+(j-501))-sr)
            rtot = rtot + rpsf_photij(1.*(i-51),1.*(j-501),ii,jj)
            if (pixc(ii+(i-51),jj+(j-501)).gt.HIFLAG) nsat = nsat+1
            endif
         enddo
         enddo
 
      do i = 001, 0101
      do j = 001, 1001
         if (u2(i,j).eq.1) pixq(ii+(i-51),jj+(j-501)) = 60000.0
         if (u1(i,j).eq.1) pixq(ii+(i-51),jj+(j-501)) = 65000.0
         enddo
         enddo
 
      apphot_sat2 = ftot/rtot
      pixq(ii,jj) = apphot_sat2
      nsat = nadd2
      write(99,199) 2,ii,jj,pixc(ii,jj),apphot_sat2,ftot,rtot,nsat
 
      return
      end
 
 
                                              
c*********************************************
c****                                         
c**** #include "SATAP/max_contig.f"
c****                                         
c*********************************************
                                              
 
 
      subroutine max_contig(pixp,pixx,UPR_LIM)
      implicit none
 
      real    pixp(4096,4096)
      real    pixx(4096,4096)
      integer i ,j
      real    UPR_LIM
 
      integer jmin
      integer jmax
 
      integer L, Ls, Lo, LL
      integer il(99999)
      integer jl(99999)
      real    pmax
      integer NIT, NCHG
 
      do i = 0001, 4096
      do j = 0001, 4096
         pixx(i,j) = pixp(i,j)
         enddo
         enddo
 
      NIT  = 1
    1 NIT  = NIT + 1
      NCHG = 0
      do j = 0002, 2047,  1
      do i = 0002, 4095,  1
         if (pixx(i,j).gt.UPR_LIM) then
            if (pixx(i,j).lt.pixx(i+1,j  )) then
               pixx(i,j) = pixx(i+1,j  )
               NCHG = NCHG + 1
               endif
            if (pixx(i,j).lt.pixx(i-1,j  )) then
               pixx(i,j) = pixx(i-1,j  )
               NCHG = NCHG + 1
               endif
            if (pixx(i,j).lt.pixx(i  ,j+1)) then
               pixx(i,j) = pixx(i  ,j+1)
               NCHG = NCHG + 1
               endif
            if (pixx(i,j).lt.pixx(i  ,j-1)) then
               pixx(i,j) = pixx(i  ,j-1)
               NCHG = NCHG + 1
               endif
            endif
         enddo
         enddo
      do j = 2047, 0002, -1
      do i = 0002, 4095,  1
         if (pixx(i,j).gt.UPR_LIM) then
            if (pixx(i,j).lt.pixx(i+1,j  )) then
               pixx(i,j) = pixx(i+1,j  )
               NCHG = NCHG + 1
               endif
            if (pixx(i,j).lt.pixx(i-1,j  )) then
               pixx(i,j) = pixx(i-1,j  )
               NCHG = NCHG + 1
               endif
            if (pixx(i,j).lt.pixx(i  ,j+1)) then
               pixx(i,j) = pixx(i  ,j+1)
               NCHG = NCHG + 1
               endif
            if (pixx(i,j).lt.pixx(i  ,j-1)) then
               pixx(i,j) = pixx(i  ,j-1)
               NCHG = NCHG + 1
               endif
            endif
         enddo
         enddo
      do j = 2050, 4095,  1
      do i = 0002, 4095,  1
         if (pixx(i,j).gt.UPR_LIM) then
            if (pixx(i,j).lt.pixx(i+1,j  )) then
               pixx(i,j) = pixx(i+1,j  )
               NCHG = NCHG + 1
               endif
            if (pixx(i,j).lt.pixx(i-1,j  )) then
               pixx(i,j) = pixx(i-1,j  )
               NCHG = NCHG + 1
               endif
            if (pixx(i,j).lt.pixx(i  ,j+1)) then
               pixx(i,j) = pixx(i  ,j+1)
               NCHG = NCHG + 1
               endif
            if (pixx(i,j).lt.pixx(i  ,j-1)) then
               pixx(i,j) = pixx(i  ,j-1)
               NCHG = NCHG + 1
               endif
            endif
         enddo
         enddo
      do j = 4095, 2050, -1
      do i = 0002, 4095,  1
         if (pixx(i,j).gt.UPR_LIM) then
            if (pixx(i,j).lt.pixx(i+1,j  )) then
               pixx(i,j) = pixx(i+1,j  )
               NCHG = NCHG + 1
               endif
            if (pixx(i,j).lt.pixx(i-1,j  )) then
               pixx(i,j) = pixx(i-1,j  )
               NCHG = NCHG + 1
               endif
            if (pixx(i,j).lt.pixx(i  ,j+1)) then
               pixx(i,j) = pixx(i  ,j+1)
               NCHG = NCHG + 1
               endif
            if (pixx(i,j).lt.pixx(i  ,j-1)) then
               pixx(i,j) = pixx(i  ,j-1)
               NCHG = NCHG + 1
               endif
            endif
         enddo
         enddo
      write(*,'(''  NIT: '',i6,1x,i8)') NIT,NCHG
      if (NCHG.gt.0) goto 1
 
      return
      end
 
                                              
c*********************************************
c****                                         
c**** #include "SATAP/peak_sat.f"
c****                                         
c*********************************************
                                              
 
c---------------------------------------------------
c
c this routine will go through the image pixr() and
c generate two supplemental images:  pixp() and pixx().
c pixp() gives a sense of how saturated each saturated
c pixel is, and pixx() tells you what the maximum value
c is for a touching saturated pixel ; if pixx.eq.pixp,
c then the pixel is the maximum of its distribution
c
      subroutine peak_sat(pixr,pixp,pixx,UPR_LIM)
      implicit none
 
      real pixr(4096,4096)
      real pixp(4096,4096)
      real pixx(4096,4096)
      real UPR_LIM
 
      integer i,j
      integer ii,jj
      integer imin, imax
      integer jmin, jmax
      integer r
 
      integer nsat, ntot
      real    ptota
      real    ptotb
      real    ptot
 
      print*,'ENTER PEAK_SAT...'
 
      do i = 0005, 4090
      if (i.eq.i/100*100) print*,'  i: ',i
      do j = 0005, 4090
         pixp(i,j) = 0
         if (j.le.2048) then
            jmin = 0001
            jmax = 2048
            endif
         if (j.ge.2049) then
            jmin = 2049
            jmax = 4096
            endif
         if (pixr(i,j).gt.UPR_LIM) then
            r = 1
            ptota = 0
    3       ptota = ptota + pixr(i,j-r)
            if (pixr(i,j-r).gt.0.75*UPR_LIM.and.j-r-1.gt.jmin) then
                r = r + 1
                goto 3
                endif
            r = 1
            ptotb = 0
    2       ptotb = ptotb + pixr(i,j+r)
            if (pixr(i,j+r).gt.0.75*UPR_LIM.and.j+r+1.lt.jmax) then
                r = r + 1
                goto 2
                endif
            r = 2
    1       continue
            nsat = 0
            ntot = 0
            ptot = 0
            do ii = -r, +r
            do jj = -r, +r
               if (i+r.lt.0001.or.
     .             i+r.gt.4096.or.
     .             j+r.lt.jmin.or.
     .             j+r.gt.jmax) then
                   pixp(i,j) = -750
                   goto 666
                   endif
               if (ii**2+jj**2.le.(r+0.5)**2) then
                  ntot = ntot + 1
                  ptot = ptot + max(pixr(i+ii,j+jj),0.0)
                  if (pixr(i+ii,j+jj).gt.UPR_LIM) nsat = nsat + 1
                  endif
               enddo
               enddo
            if (nsat.gt.0.85*ntot) then
               r = r + 1
               goto 1
               endif
            pixp(i,j) = ptot + ptota + ptotb
            endif
  666    continue
         enddo
         enddo
 
      print*,'               '
      print*,'FIND MAX_CONTIG...'
      print*,'---> go through the image pixel by pixel and determine'
      print*,'     the brightest pixp pixel that is saturatedly '
      print*,'     contiguous'
      print*,'               '
      call max_contig(pixp,pixx,UPR_LIM)
 
      print*,'               '
      print*,'               '
      end
 
 
 
                                              
c*********************************************
c****                                         
c**** #include "SATAP/sat_phot.f"
c****                                         
c*********************************************
                                              
 
c---------------------------------------------
c
c
c
      subroutine sat_phot(pixr,pixp,pixx)
      implicit none
 
      real      pixr(4096,4096)
      real      pixp(4096,4096)
      real      pixx(4096,4096)
      real      pixs(4096,4096)
 
      integer i
      integer j
      integer NARGs, NARG, iargc
 
      integer Ls
 
      logical trip
      character*80 FILEI
      character*80 FILEO
 
      real UPR_LIM
 
      integer iu
      real pixe
      integer il(99999), jl(99999)
      integer NSATs
      integer NCENs
      real sr, mbar_sky
      real zr
      real apphot_sat2, rtot
      integer nr, nr0
      integer ir
 
      real    reff
      integer ireff
      integer ii, jj
      integer ncirc
      real    rpsf_photij
      integer HIFLAG
      common / HIFLAG_ / HIFLAG
 
      integer ntot, nap
      real*8 ptot
      real*8 ftot
      real pixpu
      integer iii, jjj
 
      UPR_LIM = HIFLAG
 
      do i = 0001, 4096
         pixr(i,2048) = -750
         pixr(i,2049) = -750
         enddo
 
      do i = 0001, 4096
      do j = 0001, 4096
         pixs(i,j) = pixr(i,j)
         enddo
         enddo
 
      Ls = 0
      do i = 0001, 4096
      do j = 0001, 4096
         if (pixp(i,j).gt.HIFLAG.and.
     .       pixp(i,j).eq.pixx(i,j)) then
             Ls = Ls + 1
             write(92,*) i,j,pixp(i,j),Ls
             endif
         enddo
         enddo
      close(92)
 
      open(55,file='LOG.peaksat_CEN.reg',status='unknown')
      write(55,'(''# Region file format: DS9 version 3.0'')')
      write(55,'(''global color=yellow'')')
 
      print*,'Do apphot...'
      Ls = 0
      do i = 0001+10, 4096-10
      do j = 0001+10, 4096-10
         if (pixp(i,j).gt.UPR_LIM) then
             if (pixp(i,j).lt.pixx(i,j)) goto 47
             do ii = -7, 7
             do jj = -7, 7
                if (ii**2+jj**2.le.7.5**2.and.
     .              pixp(i,j).lt.pixp(i+ii,j+jj)) goto 47
                enddo
                enddo
             sr = mbar_sky(i,j,25,30,pixs)
             pixpu = pixp(i,j)
             nap  = 0.
             ntot = 0.
             ptot = 0.
             ftot = 0.
             do iii = max(i-50,0001+2), min(i+50,4096-2)
             do jjj = 0001+2, 4096-2
                if (pixx(iii  ,jjj  ).eq.pixpu.or.
     .              pixx(iii+1,jjj  ).eq.pixpu.or.
     .              pixx(iii-1,jjj  ).eq.pixpu.or.
     .              pixx(iii  ,jjj+1).eq.pixpu.or.
     .              pixx(iii  ,jjj-1).eq.pixpu.or.
     .              pixx(iii+1,jjj+1).eq.pixpu.or.
     .              pixx(iii+1,jjj-1).eq.pixpu.or.
     .              pixx(iii-1,jjj+1).eq.pixpu.or.
     .              pixx(iii-1,jjj-1).eq.pixpu.or.
     .              pixx(iii-2,jjj-2).eq.pixpu.or.
     .              pixx(iii-2,jjj-1).eq.pixpu.or.
     .              pixx(iii-2,jjj  ).eq.pixpu.or.
     .              pixx(iii-2,jjj+1).eq.pixpu.or.
     .              pixx(iii-2,jjj+2).eq.pixpu.or.
     .              pixx(iii-1,jjj+2).eq.pixpu.or.
     .              pixx(iii  ,jjj+2).eq.pixpu.or.
     .              pixx(iii+1,jjj+2).eq.pixpu.or.
     .              pixx(iii+2,jjj+2).eq.pixpu.or.
     .              pixx(iii+2,jjj+1).eq.pixpu.or.
     .              pixx(iii+2,jjj  ).eq.pixpu.or.
     .              pixx(iii+2,jjj-1).eq.pixpu.or.
     .              pixx(iii+2,jjj-2).eq.pixpu.or.
     .              pixx(iii+1,jjj-2).eq.pixpu.or.
     .              pixx(iii  ,jjj-2).eq.pixpu.or.
     .              pixx(iii-1,jjj-2).eq.pixpu) then
                    ptot = ptot + (pixr(iii,jjj)-sr)
                    ftot = ftot + rpsf_photij(1.*(iii-i),
     .                                        1.*(jjj-j),i,j)
                    nap  = nap  + 1
                    if (pixr(iii,jjj).gt.UPR_LIM) ntot = ntot + 1
                    endif
                enddo
                enddo
             rtot = ftot
             zr = ptot/ftot
             Ls = Ls + 1
             write( *,888) i,j,-2.5*log10(max(zr,1.0)),
     .                     sr,zr,ntot,nap,Ls,ftot,pixp(i,j), pixx(i,j)
             write(98,888) i,j,-2.5*log10(max(zr,1.0)),
     .                     sr,zr,ntot,nap,Ls,ftot,pixp(i,j), pixx(i,j)
  888        format(i4.4,1x,i4.4,1x,f10.4,1x,f8.2,
     .                   3x,f12.1,1x,i8,1x,i8,1x,i5.5,1x,f10.6,
     .                   3x,f12.1,1x,f12.1)
             pixr(i,j) = zr
             write(55,112) i,j,0.5,'blue'
  112        format('image;circle( ',i4,',',i4,',',f6.3,') # color=',a7)
             endif
   47    continue
         enddo
         enddo
      close(98)
      close(55)
 
      return
      end
 
 
 
