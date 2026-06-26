       program VI_HST_ogle_man_match4

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       parameter(maxstars=30000,nbmax=100,Nfs=100)
       real*8 x(maxstars),xe(maxstars),y(maxstars),ye(maxstars)
       real*8 x0(maxstars),y0(maxstars)
       real*8 xe0(maxstars),ye0(maxstars)
       real*8 xec(maxstars),yec(maxstars)
       real*8 amag(maxstars),ain(maxstars)
       real*8 amag0(maxstars),ain0(maxstars)
       real*8 bmag0(maxstars),bmag(maxstars)
       real*8 duma0(maxstars),duma(maxstars)
       real*8 dumb0(maxstars),dumb(maxstars)
       real*8 tr1(2),tr2(2),delxfit(2)
       integer iatype(maxstars),iseq(maxstars),iti(maxstars)
       integer iatype0(maxstars),ibtype0(maxstars),ibtype(maxstars)
       integer ityp_fid_I(maxstars),ityp_fid_V(maxstars),it(maxstars)
       real*8 fid_x(maxstars),fid_y(maxstars),fid_in(maxstars)
       real*8 fida_x(maxstars),fida_y(maxstars)
       real*8 fid_V(maxstars),fid_I(maxstars)
       integer ityp_fid_I0(maxstars),ityp_fid_V0(maxstars)
       real*8 fid_x0(maxstars),fid_y0(maxstars),fid_in0(maxstars)
       real*8 fid_V0(maxstars),fid_I0(maxstars)
       real*8 xnew(maxstars),ynew(maxstars)
       real*8 x_cmd(20),y_cmd(20),x_H(20),y_H(20)
       real*8 fI_cmd(20),HH(20)
       integer ibi(nbmax),jbi(nbmax),imat(nbmax),ibmat(nbmax)
       character*160 inline,inline2
       character*80 mapfile,psfstar_I_file,psfstar_V_file
       character*80 matchup_file_I,matchup_file_V
       double precision ra_arcsec,dec_arcsec,ra_t_arcsec,dec_t_arcsec
       double precision dec_deg,cos_dec,dec_t_deg,cos_t_dec
       double precision ra_deg
       double precision rah,ram,ras,decd,decm,decs,rad
       double precision Io_t
       double precision A1I(Nfs),Imag_fsky(Nfs),A1V(Nfs),Vmag_fsky(Nfs)
       double precision chi2_I(Nfs),chi2_V(Nfs)
       double precision lg_chi2_Imx(Nfs),lg_chi2_Vmx(Nfs)
       integer icalnum(maxstars), inumcal(maxstars), icals(maxstars)
       integer icalnum0(maxstars)

       data icals / maxstars*0 /
       data icalnum0 / maxstars*0 /

c      hardwired parameters
c      --------------------
       scale_rat = 0.26d0/0.04d0   ! OGLE I-band/HST I-band
       as_pix_hst = 0.0d04
       as_pix_ogle = 0.26d0
       mx_nmatch = 200
       mn_match = 20
       del_x_mx = 100.d0
       del_y_mx = 100.d0
       ratin_mx = 9.5d0
       del_x_tol = 20.0d0
       del_y_tol = 20.0d0
       theta_tol = 0.1d0
       r_tol = 3.0d0
       err_min = 0.20d0
       errc_min = 0.3d0
       tolout = 3.0d0
       chi_thresh = 25.0d0
       nbtotmin = 4
       pi = atan(-1.d0)

       errc_min2 = errc_min**2
       err_min2 = err_min**2

c   mb08379 I-band case:
c   log10(chi2) < lg_chi2max = a*Imag_fsky + b + 5.15
c                 3.4  = a*(-13) + b
c                 2.65 = a*(-10) + b
c                 0.75 = -3*a   so a = -0.25
c                 2.65 = 2.5 + b  so b = 0.15
c      log10(chi2) < -0.25*Imag_fsky + 0.15
c    use a = A_lg_chi2_I, A_lg_chi2_V in this code, but impliment B in the next code
c      log10(chi2) + 0.25*Imag_fsky = lg_chi2_Imx < 0.15
c      

       write(6,*) 'enter match radius in arc sec:'
       read(5,*) rad_m
       rad_m2 = rad_m**2
       rad_m2hst = rad_m2/as_pix_hst**2
       rad_m2ogle = rad_m2/as_pix_ogle**2
       write(6,*) 'enter target coordinates RA: h m s DEC: deg m s'
       read(5,*) rah,ram,ras,decd,decm,decs
       if(decd.lt.0.d0.or.inline(20:20).eq.'-'.or.
     &                    inline(21:21).eq.'-') then
         decd = decd - decm/60.d0 - decs/3600.d0
       else
         decd = decd + decm/60.d0 + decs/3600.d0
       endif
       dec_t_arcsec = decd*3600.d0
       cos_dec = cos(decd*pi/180.d0)
       rah = rah + ram/60.d0 + ras/3600.d0
       rad = rah*15.d0
       ra_t_arcsec = rad*cos_dec*3600.d0

       write(6,*) 'enter hwhma, hwhmb (arc sec):'
       read(5,*) hwhma, hwhmb
       hwhm_a = hwhma/as_pix_ogle
       hwhm_a2 = hwhm_a**2
       hwhm_b = hwhmb/as_pix_ogle
       hwhm_b2 = hwhm_b**2

       write(6,*) 'enter A_lg_chi2_I, A_lg_chi2_V:'
       read(5,*) A_lg_chi2_I, A_lg_chi2_V

c      enter I and V PSF fit mag file names and read files
c      ---------------------------------------------------
       write(6,*) 'enter PSF fit I mag file name:'
       read(5,'(a)') psfstar_I_file
       write(6,*) 'enter PSF fit V mag file name:'
       read(5,'(a)') psfstar_V_file

       open(unit=7,file=psfstar_I_file,status='old')
       kstar = 0
       do i = 1,999999
         read(7,'(a)',end=3,err=3) inline
         lin = lenc(inline)
         if(inline(1:16).eq.'# BEST-FIT MODEL') then
           read(inline(27:lin),*) kstar
           if(kstar.gt.Nfs) stop 'kstar_I > Nfs'
         elseif(inline(1:11).eq.'# chi2MIN =') then
           if (index(inline(12:lin), '*') == 0) then
            read(inline(12:lin), *) chi2_I(kstar)
           else
            print *, "Overflow detected, setting chi2_I(", kstar, 
     &        ") to -1.0"
            chi2_I(kstar) = -1.0
           endif
ccc           read(inline(12:lin),*) chi2_I(kstar)
         elseif(inline(1:11).eq.'#  A1_min =') then
           read(inline(12:lin),*) A1I(kstar)
           Imag_fsky(kstar) = -2.5d0*log10(A1I(kstar))
ccc           lg_chi2_Imx = log10(chi2_I(kstar)) + 0.25*Imag_fsky(kstar)
           lg_chi2_Imx(kstar) = log10(chi2_I(kstar))
     &                        + 0.25d0*Imag_fsky(kstar)
         endif
       enddo
 3     continue
       close(7)
       open(unit=8,file=psfstar_V_file,status='old')
       kstar = 0
       do i = 1,999999
         read(8,'(a)',end=4,err=4) inline
         lin = lenc(inline)
         if(inline(1:16).eq.'# BEST-FIT MODEL') then
           read(inline(27:lin),*) kstar
           if(kstar.gt.Nfs) stop 'kstar_V > Nfs'
         elseif(inline(1:11).eq.'# chi2MIN =') then
           if (index(inline(12:lin), '*') == 0) then
            read(inline(12:lin), *) chi2_I(kstar)
           else
            print *, "Overflow detected, setting chi2_I(", kstar, 
     &        ") to -1.0"
            chi2_I(kstar) = -1.0
           endif
ccc           read(inline(12:lin),*) chi2_V(kstar)
         elseif(inline(1:11).eq.'#  A1_min =') then
           read(inline(12:lin),*) A1V(kstar)
           Vmag_fsky(kstar) = -2.5d0*log10(A1V(kstar))
           lg_chi2_Vmx(kstar) = log10(chi2_V(kstar))
     &                        + 0.25d0*Vmag_fsky(kstar)
         endif
       enddo
 4     continue
       close(8)

c      open and read CMD file
c      ----------------------
       write(6,*) 'enter OGLE-III map file name:'
       read(5,'(a)') mapfile
       open(unit=1,file=mapfile,status='old')
       nf = 0
       r2min_arcsec = 1.e10
       do i=1,999999
         read(1,'(a)',err=11,end=11) inline
         len=lenc(inline)
         read(inline(1:7),*) id_in
         read(inline(8:9),*) rah
         read(inline(11:12),*) ram
         read(inline(14:19),*) ras
         read(inline(20:22),*) decd
         read(inline(24:25),*) decm
         read(inline(27:31),*) decs
         read(inline(31:len),*) xo,yo,fV,fVmI,fI,n1,n2,sV,n3,n4,sI
         if(decd.lt.0.d0.or.inline(20:20).eq.'-'.or.
     &                      inline(21:21).eq.'-') then
           decd = decd - decm/60.d0 - decs/3600.d0
         else
           decd = decd + decm/60.d0 + decs/3600.d0
         endif
         dec_arcsec = decd*3600.d0
         cos_dec = cos(decd*pi/180.d0)
         rah = rah + ram/60.d0 + ras/3600.d0
         rad = rah*15.d0
         ra_arcsec = rad*cos_dec*3600.d0
         r2 = (ra_arcsec-ra_t_arcsec)**2 + (dec_arcsec-dec_t_arcsec)**2
         if(r2.lt.r2min_arcsec) then
           r2min_arcsec = r2
c          x = N,  y = E
           xo_t = xo
           yo_t = yo
         endif
       enddo
 11    continue
       rewind(1)
       nstot = 0
       do i=1,999999
         read(1,'(a)',err=21,end=21) inline
         len=lenc(inline)
         read(inline(31:len),*) xo,yo,fV,fVmI,fI,n1,n2,sV,n3,n4,sI
         if(fV.lt.30.d0.and.fI.lt.30.d0.and.sV.lt.0.3d0.and.sI.lt.0.3d0)
     &       then
           dx = xo - xo_t
           dy = yo - yo_t
           if(dx**2+dy**2.le.rad_m2ogle) then
             nstot = nstot + 1
             x0(nstot) = xo
             y0(nstot) = yo
             bmag0(nstot) = fV
             amag0(nstot) = fI
             ain0(nstot) = 10**(0.4d0*(21.d0-amag0(nstot)))
             iatype0(nstot) = 1
             ibtype0(nstot) = nstot
             xe0(nstot) = err_min
             ye0(nstot) = err_min
             xe(nstot) = err_min
             ye(nstot) = err_min
             xec(nstot) = err_min
             yec(nstot) = err_min
           endif
         endif
       enddo
 21    continue
       close(1)

       write(6,*) nstot,' stars read in from OGLE file'

       call sort_by_intensity(x0,y0,amag0,bmag0,duma0,ain0,dumb0,
     &                    iatype0,ibtype0,nstot,x,y,amag,bmag,duma,
     &                    ain,dumb,iatype,ibtype,it)


c      enter the matched coordinates
c      -----------------------------
       write(6,*) 'enter matched coordinates'
       write(6,*) 'format x_hst y_hst x_ogle y_ogle I_hst I_ogle,',
     &            ' 0 0 0 0 0 0 to end:'
       nman = 1
       do
         read(5,*) x_cmd(nman),y_cmd(nman),x_H(nman),y_H(nman),
     &           fI_cmd(nman),HH(nman)
         if(x_cmd(nman)**2+x_H(nman)**2.le.0.d0) exit
         nman = nman + 1
       enddo
       nman = nman - 1
       ratin_av = 0.d0
       do i = 1,nman
         fI_man = 10**(0.4d0*(21.d0-fI_cmd(i)))
         fH_man = 10**(0.4d0*(21.d0-HH(i)))
         rat = fH_man/fI_man
         ratin_av = ratin_av + rat
         xe(i) = errc_min
         ye(i) = errc_min
         ibi(i) = i
         jbi(i) = i
       enddo
       ratin_av = ratin_av/nman

c      open and read Jay Anderson format MATCHUP.XYMEEE files (F814W & F555W)
c      ----------------------------------------------------------------------
       write(6,*) 'enter matchup_file_I and matchup_file_V names:'
       read(5,'(a)') matchup_file_I
       read(5,'(a)') matchup_file_V
       open(unit=2,file=matchup_file_I,status='old')
       open(unit=3,file=matchup_file_V,status='old')

c      MATCHUP file lists for each star
c      x, y, mag, x_sig, y_sig, mag_sig, NIMf, NIMg, NIMm, NIMMINu,NPK, pki(NPK),pkj(NPK), 
c        calibration star # (only for calibration stars)
c        each line of the matchup_file_I file covers the same star as the matchup_file_V file
c      SEAN NOTE: The comment about columns names above is out-dated, inconsistent with newer
c      MATCHUP.XYM formatting. I need to make a modification to account for older format.
       nf = 0
       ncal = 0
       do i=1,99999
         read(2,'(a)',err=31,end=31) inline
         read(3,'(a)',err=31,end=31) inline2
         if(inline(1:1).ne.'#') then
           read(inline,*) xxx,yyy,amagin,xxxerr,yyyerr,amagerr
           lin = lenc(inline)
           if(lin.gt.133.and.inline(140:140).ne.'t') then
             read(inline(135:140),*) icalin
             if(icalin.gt.0) then
               ncal = ncal + 1
               icalnum0(i) = icalin
             else
               icalnum0(i) = 0
             endif
           endif
         endif
         if(inline2(1:1).ne.'#') then
           read(inline2,*) xx,yy,bmagin,xxerr,yyerr,bmagerr
         endif
         nf = nf + 1
         fid_x0(i) = xxx
         fid_y0(i) = yyy
         fid_I0(i) = amagin
         fid_V0(i) = bmagin
         fid_in0(i) = 10**(0.4d0*(21.d0-fid_I0(i)))
         x0(nf) = xxx
         y0(nf) = yyy
       enddo
 31    continue
       close(2)
       close(3)
       n_fid_tot = nf

       write(6,*) n_fid_tot,' stars read in from MATCHUP files'

       call sort_by_intensity(fid_x0,fid_y0,fid_I0,fid_V0,duma0,
     &                    fid_in0,dumb0,icalnum0,ityp_fid_V0,
     &                    n_fid_tot,fid_x,fid_y,fid_I,fid_V,duma,
     &                    fid_in,dumb,icalnum,ityp_fid_V,it)

       do i = 1,n_fid_tot
          ic = icalnum(i)
          if(ic.gt.0) inumcal(ic) = i
       enddo

c      fit the precise transformation for the manually entered coords
c      --------------------------------------------------------------
       CALL fldtrans(nbmax,nman,ibi,jbi,maxstars,x_H,xe,y_H,ye,
     &               x_cmd,y_cmd,tr1,tr2,delxfit,tolout,nbtotmin,ier)

c   write out transformation parameters
c   -----------------------------------
       write(6,*) 'initial transformation parameters:'
       write(6,*)' TR11=',tr1(1),' TR12=',tr2(1),' TR21=',tr1(2),
     &      ' TR22=',tr2(2),' delX=',delxfit(1),' delY=',delxfit(2)

c  make the translations from the template coordinates to the field coords.
c  ------------------------------------------------------------------------
       do 300 i=1,n_fid_tot
          fida_x(i)=tr1(1)*fid_x(i)+tr2(1)*fid_y(i)+delxfit(1)
          fida_y(i)=tr1(2)*fid_x(i)+tr2(2)*fid_y(i)+delxfit(2)
 300   continue

c   match Fiducial HST stars to OGLE stars based on crude transformation
c   --------------------------------------------------------------------
       nbtot = 0
       j0 = 1
       rat_mx = 4.0d0*ratin_av
       rat_mn = 0.25d0*ratin_av
       do 400 i=1,n_fid_tot
          do 380 j=j0,nstot
             rat = ain(j)/fid_in(i)
             if(rat.gt.rat_mx) then
               j0=min(j+1,nstot)
             elseif(rat.ge.rat_mn) then
ccc               if(iatype(j).eq.1.and.ityp_fid_I(i).eq.1.and.
ccc     &                               ityp_fid_V(i).eq.1) then
                 dx=x(j)-fida_x(i)
                 dy=y(j)-fida_y(i)
                 chi=(dx/xec(j))**2+(dy/yec(j))**2
                 if(chi.le.chi_thresh) then
                   nbtot=nbtot+1
                   ibi(nbtot)=i
                   jbi(nbtot)=j
                   if(nbtot.ge.nbmax) go to 401
ccc                 endif
               endif
             else
               go to 381
             endif
 380      continue
 381      continue
 400   continue
 401   continue

c      fit the precise transformation
c      ------------------------------
       CALL fldtrans(nbmax,nbtot,ibi,jbi,maxstars,x,xe,y,ye,
     &               fida_x,fida_y,tr1,tr2,delxfit,tolout,nbtotmin,ier)

c   write out transformation parameters
c   -----------------------------------
       write(6,*) '2nd transformation:'
       write(6,*)' TR11=',tr1(1),' TR12=',tr2(1),' TR21=',tr1(2),
     &      ' TR22=',tr2(2),' delX=',delxfit(1),' delY=',delxfit(2)

       open(unit=3,file='VI_HST_ogle_all_matches4.dat',status='new')
       open(unit=4,file='VI_HST_ogle_Cal_matches4.dat',status='new')
       lIfile = lenc(psfstar_I_file)
       lVfile = lenc(psfstar_V_file)
       write(3,'(a)') '# I_fsky file =',psfstar_I_file(1:lIfile)
       write(4,'(a)') '# I_fsky file =',psfstar_I_file(1:lIfile)
       write(3,'(a)') '# V_fsky file =',psfstar_V_file(1:lVfile)
       write(4,'(a)') '# V_fsky file =',psfstar_V_file(1:lVfile)
       write(3,994) rad_m, hwhma, hwhmb, A_lg_chi2_I, A_lg_chi2_V
       write(4,994) rad_m, hwhma, hwhmb, A_lg_chi2_I, A_lg_chi2_V
       write(3,996)
       write(4,997)

c   do the final match with the precise transformation with
c   all stars within 1.0 HWHM of the OGLE star
c   -------------------------------------------------------
       i0=1
       do 500 j=1,nstot
          nmatch = 0
          nmatchb = 0
          do 450 i=i0,n_fid_tot
             rat = ain(j)/fid_in(i)
             if(rat.lt.0.1d0*rat_mn) then
               i0=min(i+1,n_fid_tot)
ccc             elseif(rat.gt.5.d0*rat_mx) then
ccc               go to 451
             else
               xnew(i)=tr1(1)*fida_x(i)+tr2(1)*fida_y(i)+delxfit(1)
               ynew(i)=tr1(2)*fida_x(i)+tr2(2)*fida_y(i)+delxfit(2)
               dx=x(j)-xnew(i)
               dy=y(j)-ynew(i)
               dr2=dx**2+dy**2
               if(dr2.le.hwhm_a2) then
                 nmatch=nmatch+1
                 imat(nmatch)=i
                 if(icalnum(i).gt.0) then
                   if(icals(j).gt.10000) then
                   elseif(icals(j).gt.100) then
                     icals(j) = 10000*icalnum(i) + icals(j)
                   elseif(icals(j).gt.0) then
                     icals(j) = 100*icalnum(i) + icals(j)
                   else
                     icals(j) = icalnum(i)
                   endif
                 endif
               elseif(dr2.le.hwhm_b2) then
                 nmatchb = nmatchb + 1
                 ibmat(nmatchb)=i
               endif
             endif
 450      continue
 451      continue
          if(nmatch.eq.0) then
            sImag=99.999d0
            sVmag=99.999d0
            sImag1=99.999d0
            sVmag1=99.999d0
          elseif(nmatch.eq.1) then
            sImag=fid_I(imat(1))
            sVmag=fid_V(imat(1))
            sImag1=fid_I(imat(1))
            sVmag1=fid_V(imat(1))
          else
            sImag1=fid_I(imat(1))
            sVmag1=fid_V(imat(1))
            ampI=0.d0
            ampV=0.d0
            do im=1,nmatch
              i=imat(im)
              ampI = ampI+10**(0.4d0*(21.d0-fid_I(i)))
              ampV = ampV+10**(0.4d0*(21.d0-fid_V(i)))
            enddo
            sImag=21.d0-2.5d0*log10(ampI)
            sVmag=21.d0-2.5d0*log10(ampV)
          endif
          if(nmatchb.eq.0) then
            bImag=99.999d0
            bVmag=99.999d0
          elseif(nmatchb.eq.1) then
            bImag=fid_I(ibmat(1))
            bVmag=fid_V(ibmat(1))
          else
            bampI=0.d0
            bampV=0.d0
            do im=1,nmatchb
              i=ibmat(im)
              bampI = bampI+10**(0.4d0*(21.d0-fid_I(i)))
              bampV = bampV+10**(0.4d0*(21.d0-fid_V(i)))
            enddo
            bImag=21.d0-2.5d0*log10(bampI)
            bVmag=21.d0-2.5d0*log10(bampV)
          endif
          sVmIo=sVmag-sImag
          sVomVc = sVmag-bmag(j)
          sIomIc = sImag-amag(j)
          sV1omVc = sVmag1-bmag(j)
          sI1omIc = sImag1-amag(j)
          write(3,998) icals(j),x(j),y(j),nmatch,nmatchb,
     &       amag(j),bmag(j),sImag1,sImag,bImag,sVmag1,sVmag,bVmag,
     &       sI1omIc,sIomIc,sV1omVc,sVomVc
          if(icals(j).gt.0) then
            ii = icals(j) - 100*(icals(j)/100)
            sIomIfs = amag(j) - Imag_fsky(ii)
            sVomVfs = bmag(j) - Vmag_fsky(ii)
            write(4,999) icals(j),x(j),y(j),nmatch,nmatchb,amag(j),
     &        bmag(j),sImag1,sImag,Imag_fsky(ii),sIomIfs,
     &        lg_chi2_Imx(ii),bImag,sVmag1,sVmag,Vmag_fsky(ii),sVomVfs,
     &        lg_chi2_Vmx(ii),bVmag,sI1omIc,sIomIc,sV1omVc,sVomVc
          endif
 500   continue
       close(3)
 994   format('# match radius in arc sec =',f8.2,' HWHMa =',f7.2,
     &        ' HWHMb =',f7.2,'  A_lg_chi2_I, A_lg_chi2_V =',2f8.3)
 996   format('# Cal#   x        y     nmat nbmat I_ogle  V_ogle ',
     &        ' I_hst1  I_hst   I_hstB  V_hst1  V_hst   V_hstB ',
     &        ' Io-Ih1  Io-Ih   Vo-Vh1  Vo-Vh')
 997   format('# Cal#   x        y     nmat nbmat I_ogle  V_ogle ',
     &        ' I_hst1  I_hst   I_hfs   Io-Ihfs lg_c2Imx I_hstB ',
     &        ' V_hst1  V_hst   V_hfs   Vo-Vhfs lg_c2Vmx V_hstB ',
     &        ' Io-Ih1  Io-Ih   Vo-Vh1  Vo-Vh')
 998   format(i6,2F9.3,2i5,12f8.3)
 999   format(i6,2F9.3,2i5,5f8.3,f9.4,5f8.3,f9.4,6f8.3)

       stop
       END

c==============================================================================
C $Id: sort_by_intensity.F,v 1.1 1993/03/20 17:54:52 bennett Exp $ $Log: sort_by
c Revision 1.5  1992/10/21  12:39:59  bennett
c
      SUBROUTINE sort_by_intensity(x,y,fI,fV,fH,in,fK,itI,itV,
     &             nstars,tx,ty,tfI,tfV,tfH,tin,tfK,ittI,ittV,it)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real*8 in(nstars),tin(nstars)
      integer it(nstars)
      real*8 x(nstars),y(nstars),tx(nstars),ty(nstars)
      real*8 fI(nstars),fV(nstars),tfI(nstars),tfV(nstars)
      real*8 fH(nstars),tfH(nstars),fK(nstars),tfK(nstars)
      integer itI(nstars),itV(nstars),ittI(nstars),ittV(nstars)


      call indexx(nstars,in,it)

      do 10 i=1,nstars
         tin(i)=in(it(1+nstars-i))
         tx(i)=x(it(1+nstars-i))
         ty(i)=y(it(1+nstars-i))
         tfI(i)=fI(it(1+nstars-i))
         tfV(i)=fV(it(1+nstars-i))
         tfH(i)=fH(it(1+nstars-i))
         tfK(i)=fK(it(1+nstars-i))
         ittI(i)=itI(it(1+nstars-i))
         ittV(i)=itV(it(1+nstars-i))
10    continue

      return
      end
c==============================================================================
C $Id: fldtrans.F,v 1.2 1993/04/23 10:19:41 bennett Exp $ $Log: fldtrans.F,v $
c Revision 1.7  1992/10/21  12:39:14  bennett
c
       SUBROUTINE fldtrans(nbmax,nbtot,ibi,jbi,maxstars,x,xe,y,ye,
     &               fida_x,fida_y,tr1,tr2,delx,tolout,nbtotmin,ier)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (MMAX = 11)
      PARAMETER(NMAX=200,toler=1.d-16,MA=3)
      DIMENSION Aa(MA),V(MA,MA),U(NMAX,MA),W(MA),B(NMAX)
      real*8 tr1(2),tr2(2),delx(2)
       dimension ibi(nbmax),ibi0(nmax),ibout(nmax)
       dimension jbi(nbmax),jbi0(nmax)
       dimension fida_x(maxstars),fida_y(maxstars)
       dimension x(maxstars),xe(maxstars),y(maxstars),ye(maxstars)
C:
C:  Find field transformation via SVD decomposition
C:             (guts are from Num. Rec. SVDFIT.FOR).
C:

      ier=0

      tolout2=tolout*tolout

 1    continue

c   x-position fit
c   --------------
      DO 12 ib=1,nbtot
        i=ibi(ib)
        j=jbi(ib)
        u(ib,1)=1./xe(j)
        u(ib,2)=fida_x(i)*u(ib,1)
        u(ib,3)=fida_y(i)*u(ib,1)
        B(ib)=x(j)*u(ib,1)
12    CONTINUE
      CALL SVDCMP(U,nbtot,MA,nmax,MA,W,V,ier)
      if(ier.ne.0) RETURN
      WMAX=0.
      DO 13 J=1,MA
        IF(W(J).GT.WMAX) WMAX=W(J)
13    CONTINUE
      wTHRESH=toler*WMAX
      DO 14 J=1,MA
        IF(W(J).LT.wTHRESH) W(J)=0.
14    CONTINUE
      CALL SVBKSB(U,W,V,nbtot,MA,nmax,MA,B,Aa)

c     find the RMS error and throw out "outliers"
c     ------------------------------------------
      nout=0
      fitCHISQ=0.
      DO 15 ib=1,nbtot
        i=ibi(ib)
        j=jbi(ib)
        ibi0(ib)=i
        jbi0(ib)=j
        SUM=Aa(1)+Aa(2)*fida_x(i)+Aa(3)*fida_y(i)
        delxsig2=(x(j)-SUM)**2
        fitCHISQ=fitCHISQ+delxsig2
15    CONTINUE
      delxsig2av=fitchisq/nbtot

      DO 16 ib=1,nbtot
        i=ibi(ib)
        j=jbi(ib)
        ibi0(ib)=i
        jbi0(ib)=j
        SUM=Aa(1)+Aa(2)*fida_x(i)+Aa(3)*fida_y(i)
        delsig2=(x(j)-SUM)**2
ccc           write(6,*) ' star #',i,' x fit error ',x(j)-SUM
        if(delsig2.gt.tolout2*delxsig2av) then
           nout=nout+1
           ibout(ib)=1
        else
           ibout(ib)=0
        endif
16    CONTINUE

      if(nout.gt.0) then
c        throw out "outliers"
c        -------------------
         ib=0
         do 17 ib0=1,nbtot
            i=ibi0(ib0)
            if(ibout(ib0).eq.0) then
               ib=ib+1
               ibi(ib)=ibi0(ib0)
               jbi(ib)=jbi0(ib0)
            endif
 17      continue
         nbtot=ib
         if(nbtot.lt.nbtotmin) then
            write(6,*) ' TOO MANY OUTLIERS IN FIELD TRANS FIT'
            ier=1
            RETURN
         endif
         write(6,*) ' retrying FIELD TRANS FIT due to'
     &              ,nout,' x outliers'
         go to 1
      endif

c   update x-transformation parameters
c   ----------------------------------
      delx(1)=Aa(1)
      tr1(1)=Aa(2)
      tr2(1)=Aa(3)

 101  continue

c   y-position fit
c   --------------
      DO 112 ib=1,nbtot
        i=ibi(ib)
        j=jbi(ib)
        u(ib,1)=1./ye(j)
        u(ib,2)=fida_x(i)*u(ib,1)
        u(ib,3)=fida_y(i)*u(ib,1)
        B(ib)=y(j)*u(ib,1)
112   CONTINUE
      CALL SVDCMP(U,nbtot,MA,nmax,MA,W,V,ier)
      if(ier.ne.0) RETURN
      WMAX=0.
      DO 113 J=1,MA
        IF(W(J).GT.WMAX) WMAX=W(J)
113   CONTINUE
      wTHRESH=toler*WMAX
      DO 114 J=1,MA
        IF(W(J).LT.wTHRESH) W(J)=0.
114   CONTINUE
      CALL SVBKSB(U,W,V,nbtot,MA,nmax,MA,B,Aa)

c     find the average error and throw out "outliers"
c     -----------------------------------------------
      nout=0
      fitCHISQ=0.
      DO 115 ib=1,nbtot
        i=ibi(ib)
        j=jbi(ib)
        ibi0(ib)=i
        jbi0(ib)=j
        SUM=Aa(1)+Aa(2)*fida_x(i)+Aa(3)*fida_y(i)
        delysig2=(y(j)-SUM)**2
        fitCHISQ=fitCHISQ+delysig2
115   CONTINUE
      delysig2av=fitchisq/nbtot

      DO 116 ib=1,nbtot
        i=ibi(ib)
        j=jbi(ib)
        ibi0(ib)=i
        jbi0(ib)=j
        SUM=Aa(1)+Aa(2)*fida_x(i)+Aa(3)*fida_y(i)
        delsig2=(y(j)-SUM)**2
ccc           write(6,*) ' star #',i,' y fit error ',y(j)-SUM
        if(delsig2.gt.tolout2*delysig2av) then
           nout=nout+1
           ibout(ib)=1
        else
           ibout(ib)=0
        endif
116   CONTINUE

      if(nout.gt.0) then
c        throw out "outliers"
c        -------------------
         ib=0
         do 117 ib0=1,nbtot
            i=ibi0(ib0)
            if(ibout(ib0).eq.0) then
               ib=ib+1
               ibi(ib)=ibi0(ib0)
               jbi(ib)=jbi0(ib0)
            endif
 117     continue
         nbtot=ib
         if(nbtot.lt.nbtotmin) then
            write(6,*) ' TOO MANY OUTLIERS IN FIELD TRANS FIT'
            ier=1
            RETURN
         endif
         write(6,*) ' retrying FIELD TRANS FIT due to'
     &              ,nout,' y outliers'
         go to 1
      endif

c   update y-transformation parameters
c   ----------------------------------
      delx(2)=Aa(1)
      tr1(2)=Aa(2)
      tr2(2)=Aa(3)

       avxerr=sqrt(delxsig2av)
       avyerr=sqrt(delysig2av)
       write(6,*) nbtot,' stars used for coord. trans. fit:'
       write(6,*) 'mean x error =',avxerr,' mean y error =',avyerr

       RETURN
       END
c==============================================================================
C $Id: svdcmp.F,v 1.4 1992/10/21 12:40:05 bennett Exp $ Initial revision
c Initial revision
c
      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V,ier)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=100)
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

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=100)
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
c==============================================================================
C $Id: indexx.F,v 1.1 1993/03/20 17:47:28 bennett Exp $ $Log: indexx.F,v $
c
      SUBROUTINE INDEXX(N,ARRIN,INDX)
c
c   from Numerical Recipes p. 233
c
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
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
