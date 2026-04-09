       program fit_HST_Iogle_col

       parameter(nmax=50,nparmax=6)
       real a(nparmax),covar(nparmax,nparmax)
       real E_mag(nmax),I_mag(nmax),V_mag(nmax),E_err(nmax)
       real E_mag_mod(nmax)
       real E_mag_targ,I_mag_targ,V_mag_targ,E_err_targ
       integer num(nmax),nmat(nmax),nbmat(nmax)
       real I_ogle(nmax),V_ogle(nmax),I_hst1(nmax),I_hst(nmax),
     &      I_hfs(nmax),IomIhfs(nmax),lg_c2Imx(nmax),I_hstB(nmax),
     &      V_hst1(nmax),V_hst(nmax),V_hfs(nmax),
     &      VomVhfs(nmax),lg_c2Vmx(nmax),V_hstB(nmax),IomIh1(nmax),
     &      IomIh(nmax),VomVh1(nmax),VomVh(nmax)
       real Ihstm1(nmax),IhstmB(nmax),Vhstm1(nmax),VhstmB(nmax)
       real I_1mhst_mx, V_1mhst_mx, I_Bmhst1_mx, V_Bmhst1_mx
       real lg_c2Imax, lg_c2Vmax
       real IomIhfs_av, IomIhfs_rms, IomIhfs_sig
       real VomVhfs_av, VomVhfs_rms, VomVhfs_sig
       integer igood(nmax),igoodI(nmax),igoodV(nmax)
       integer id(nmax),lista(nparmax)
       character*80 infile1,infile2,outfile
       character*200 inline

       external colmagfit

       ma = nparmax
       ncvm = nparmax

       write(6,*) 'enter input file name:'
       read(5,'(a)') infile1
ccc       read(5,'(a)') infile2
ccc       read(5,'(a)') outfile
ccc       lin1 = lenc(infile1)
ccc       lin2 = lenc(infile2)
ccc       outfile = '# '//infile1(1:lin1)//' '//infile2(1:lin2)
       write(6,*) 'enter iorder, errorI, errorV, VmImin, VmImax,',
     &            ' Imax, chi2thresh:'
       read(5,*) iorder,errorI,errorV,VmImin,VmImax,fImax,chi2thresh
       write(6,*)
     &    'enter I_1mhst_mx, V_1mhst_mx, I_Bmhst1_mx, V_Bmhst1_mx:'
       read(5,*) I_1mhst_mx, V_1mhst_mx, I_Bmhst1_mx, V_Bmhst1_mx
       write(6,*) 'enter lg_c2Imax, lg_c2Vmax'
       read(5,*) lg_c2Imax, lg_c2Vmax

       open(unit=1,file=infile1,status='old')
       do i = 1,6
         read(1,*)
       enddo
       nr = 0
       do i=1,nmax
         read(1,'(a)',end=11,err=11) inline
         lin = lenc(inline)
         read(inline(5:lin),*) num(i),x,y,nmat(i),nbmat(i),I_ogle(i),
     &     V_ogle(i),I_hst1(i),I_hst(i),I_hfs(i),IomIhfs(i),lg_c2Imx(i),
     &     I_hstB(i),V_hst1(i),V_hst(i),V_hfs(i),VomVhfs(i),
     &     lg_c2Vmx(i),V_hstB(i),IomIh1(i),IomIh(i),VomVh1(i),VomVh(i)
         VmI = V_ogle(i) - I_ogle(i)
         if(I_ogle(i).le.fImax.and.VmI.ge.VmImin.and.VmI.le.VmImax) then
           nr = nr + 1
           id(nr) = num(i)
           Ihstm1(nr) = I_hst1(nr)-I_hst(nr)
           IhstmB(nr) = I_hst1(nr)-I_hstB(nr)
           Vhstm1(nr) = V_hst1(nr)-V_hst(nr)
           VhstmB(nr) = V_hst1(nr)-V_hstB(nr)
         endif
       enddo
 11    continue
       close(1)

       ngood = 0
       ngoodI = 0
       ngoodV = 0
       nbadIhm1 = 0
       nbadVhm1 = 0
       nbadIh1mB = 0
       nbadVh1mB = 0
       nbadc2Imx = 0
       nbadc2Vmx = 0

       do i = 1,nr
         igI = 0
         igV = 0
         if(I_hst1(i)-I_hst(i).ge.I_1mhst_mx) then
           nbadIhm1 = nbadIhm1 + 1
         elseif(I_hstB(i)-I_hst1(i).le.I_Bmhst1_mx) then
           nbadIh1mB = nbadIh1mB + 1
         elseif(lg_c2Imx(i).ge.lg_c2Imax) then
           nbadc2Imx = nbadc2Imx + 1
         else
           ngoodI = ngoodI + 1
           igoodI(ngoodI) = i
           igI = 1
         endif
         if(V_hst1(i)-V_hst(i).ge.V_1mhst_mx) then
           nbadVhm1 = nbadVhm1 + 1
         elseif(V_hstB(i)-V_hst1(i).le.V_Bmhst1_mx) then
           nbadVh1mB = nbadVh1mB + 1
         elseif(lg_c2Vmx(i).ge.lg_c2Vmax) then
           nbadc2Vmx = nbadc2Vmx + 1
         else
           ngoodV = ngoodV + 1
           igoodV(ngoodV) = i
           igV = 1
         endif
         if(igI*igV.gt.0) then
           ngood = ngood + 1
           igood(ngood) = i
         endif
       enddo

       write(6,941) nbadIhm1,nbadIh1mB,nbadc2Imx
       write(6,942) nbadVhm1,nbadVh1mB,nbadc2Vmx
 941   format('I:',i4,' fail HWHMa blend test;',i4,
     &        ' fail HWHMb blend test;',i4,' fail chi^2 test')
 942   format('V:',i4,' fail HWHMa blend test;',i4,
     &        ' fail HWHMb blend test;',i4,' fail chi^2 test')

c      report selected calibration stars
c      ---------------------------------
       write(6,*) ngoodI,' I-band calibration stars:'
       write(6,951) (num(igoodI(j)), j = 1,ngoodI)
       write(6,*) ngoodV,' V-band calibration stars:'
       write(6,951) (num(igoodV(j)), j = 1,ngoodV)
       write(6,*) ngood,' dual-band calibration stars:'
       write(6,951) (num(igood(j)), j = 1,ngood)

 951   format(20i4)

       open(unit=2,file='fit_HST_IV_ogle_col.log',position='append')
       write(2,'(a)') 
     & '==============================================================='
       write(2,'(a)') infile1
       write(2,981) iorder,errorI,errorV,VmImin,VmImax
       write(2,982) fImax,chi2thresh
       write(2,983) I_1mhst_mx, V_1mhst_mx, I_Bmhst1_mx, V_Bmhst1_mx,
     &              lg_c2Imax, lg_c2Vmax

 981   format('# iorder =',i3,' errorI, errorV =',2f7.3,'  ',f7.3,
     &   ' <= V-I <=',f7.3)
 982   format('# I <=',f7.3,'      chi2thresh =',f8.3)
 983   format('#     Max           Max           Max           Max',/,
     &   '# I_hst1-I_hst  V_hst1-V_hst  I_hstB-I_hst1 V_hstB-V_hst1',
     &   '   lg_c2Imax,    lg_c2Vmax',/,6f14.4)

       write(2,*) ngoodI,' I-band calibration stars:'
       write(2,951) (num(igoodI(j)), j = 1,ngoodI)
       write(2,*) ngoodV,' V-band calibration stars:'
       write(2,951) (num(igoodV(j)), j = 1,ngoodV)
       write(2,*) ngood,' dual-band calibration stars:'
       write(2,951) (num(igood(j)), j = 1,ngood)

c      determine single color offsets
c      ------------------------------
       IomIhfs_av = 0.
       do j = 1,ngoodI
         IomIhfs_av = IomIhfs_av + IomIhfs(igoodI(j))
       enddo
       IomIhfs_av = IomIhfs_av/ngoodI
       IomIhfs_rms = 0.
       do j = 1,ngoodI
         IomIhfs_rms = IomIhfs_rms + (IomIhfs(igoodI(j))-IomIhfs_av)**2
       enddo
       IomIhfs_rms = sqrt(IomIhfs_rms/ngoodI)
       if(ngoodI.gt.1) then
         fgoodm1 = ngood - 1
         IomIhfs_sig = IomIhfs_rms/sqrt(fgoodm1)
       else
         IomIhfs_sig = 0.
       endif
       write(2,952) IomIhfs_av, IomIhfs_sig, IomIhfs_rms
       write(6,952) IomIhfs_av, IomIhfs_sig, IomIhfs_rms
 952   format('I_ogle-I_hst =',f9.4,' +-',f8.4,'   RMS =',f8.4)

       VomVhfs_av = 0.
       do j = 1,ngoodV
         VomVhfs_av = VomVhfs_av + VomVhfs(igoodV(j))
       enddo
       VomVhfs_av = VomVhfs_av/ngoodV
       VomVhfs_rms = 0.
       do j = 1,ngoodV
         VomVhfs_rms = VomVhfs_rms + (VomVhfs(igoodV(j))-VomVhfs_av)**2
       enddo
       VomVhfs_rms = sqrt(VomVhfs_rms/ngoodV)
       if(ngoodV.gt.1) then
         fgoodVm1 = ngoodV - 1
         VomVhfs_sig = VomVhfs_rms/sqrt(fgoodVm1)
       else
         VomVhfs_sig = 0.
       endif
       write(2,953) VomVhfs_av, VomVhfs_sig, VomVhfs_rms
       write(6,953) VomVhfs_av, VomVhfs_sig, VomVhfs_rms
 953   format('V_ogle-V_hst =',f9.4,' +-',f8.4,'   RMS =',f8.4)

c      Fit for I_ogle as a function of I, V hst - fit_sky photometry
c      -------------------------------------------------------------
       n = ngood
       do j = 1,ngood
         i = igood(j)
         E_mag(j) = I_ogle(i)
         I_mag(j) = I_hfs(i)
         V_mag(j) = V_hfs(i)
         E_err(j) = errorI
         id(j) = num(i)
       enddo
       
 15    continue

       do i=1,3
         lista(i) = i
         a(i) = 1.
       enddo
       lista(1) = 1
       lista(2) = 3
       mfit = 2
       do i=3,5
         lista(i) = 0
         a(i) = 0.
       enddo

       call LFIT_VI(I_mag,V_mag,E_mag,E_err,n,a,ma,lista,mfit,covar,
     &             ncvm,chi2,colmagfit)

       if(iorder.gt.1) then
         mfit = iorder + 1
         do i=3,mfit
           lista(i) = i+1
         enddo
         do i=mfit+1,6
           lista(i) = 0
           a(i) = 0.
         enddo
         call LFIT_VI(I_mag,V_mag,E_mag,E_err,n,a,ma,lista,mfit,covar,
     &             ncvm,chi2,colmagfit)

       endif

       VmI_targ = V_mag_targ-I_mag_targ
       E_mag_modt = a(1)+a(2)*I_mag_targ+a(3)*VmI_targ
     &             +a(4)*VmI_targ**2+a(5)*VmI_targ**3
       chi2max = 0.
       chi2sum = 0.
       varsum = 0.
       do i=1,n
         VmI = V_mag(i)-I_mag(i)
         E_mag_mod(i) = a(1)+a(2)*I_mag(i)+a(3)*VmI
     &                 +a(4)*VmI**2+a(5)*VmI**3
         var = (E_mag_mod(i)-E_mag(i))**2
         chi2 = ((E_mag_mod(i)-E_mag(i))/E_err(i))**2
         if(chi2.gt.chi2max) then
           chi2max = chi2
           imax = i
         endif
         chi2sum = chi2sum + chi2
         varsum = varsum + var
       enddo
       if(chi2max.gt.chi2thresh) then
         write(6,*) 'removing star',id(imax),'  chi^2 =',chi2max
         if(imax.ne.n) then
           do i = imax,n-1
             id(i) = id(i+1)
             E_mag(i) = E_mag(i+1)
             I_mag(i) = I_mag(i+1)
             V_mag(i) = V_mag(i+1)
             E_err(i) = E_err(i+1)
           enddo
         endif
         n = n - 1
         go to 15
       endif
       rms = sqrt(varsum/n)

       write(6,990) n,chi2sum,rms,(a(i), i=1,mfit+1)
       write(2,990) n,chi2sum,rms,(a(i), i=1,mfit+1)
       write(2,992)
       write(2,994) (id(i),E_mag(i),E_mag_mod(i),E_err(i),
     &             ((E_mag(i)-E_mag_mod(i))/E_err(i))**2,
     &               I_mag(i),V_mag(i),V_mag(i)-I_mag(i), i=1,n)

 990   format(
     &    '#  N     chi2      RMS   Io = const    *Ihst    *(V-I)hst ',
     &        '*(V-I)I**2 *(V-I)I**3',/,i5,f10.3,f9.5,5f11.6)
 991   format(
     &    '#  N     chi2      RMS   Vo = const    *Ihst    *(V-I)hst ',
     &        '*(V-I)I**2 *(V-I)I**3',/,i5,f10.3,f9.5,5f11.6)
 992   format('#id      I_ogle   I_oglem  I_og_err chi^2   ',
     &        'I_hst   V_hst   V-I')
 993   format('#id      V_ogle   V_oglem  V_og_err chi^2   ',
     &        'I_hst   V_hst   V-I')
 994   format(i7,2f9.4,f8.4,f9.3,3f8.3)

c      Fit for V_ogle as a function of I, V hst - fit_sky photometry
c      -------------------------------------------------------------
       n = ngood
       do j = 1,ngood
         i = igood(j)
         E_mag(j) = V_ogle(i)
         I_mag(j) = I_hfs(i)
         V_mag(j) = V_hfs(i)
         E_err(j) = errorV
       enddo
       
 25    continue

       do i=1,3
         lista(i) = i
         a(i) = 1.
       enddo
       lista(1) = 1
       lista(2) = 3
       mfit = 2
       do i=3,5
         lista(i) = 0
         a(i) = 0.
       enddo

       call LFIT_VI(I_mag,V_mag,E_mag,E_err,n,a,ma,lista,mfit,covar,
     &             ncvm,chi2,colmagfit)

       if(iorder.gt.1) then
         mfit = iorder + 1
         do i=3,mfit
           lista(i) = i+1
         enddo
         do i=mfit+1,6
           lista(i) = 0
           a(i) = 0.
         enddo
         call LFIT_VI(I_mag,V_mag,E_mag,E_err,n,a,ma,lista,mfit,covar,
     &             ncvm,chi2,colmagfit)

       endif

       VmI_targ = V_mag_targ-I_mag_targ
       E_mag_modt = a(1)+a(2)*I_mag_targ+a(3)*VmI_targ
     &             +a(4)*VmI_targ**2+a(5)*VmI_targ**3
       chi2max = 0.
       chi2sum = 0.
       varsum = 0.
       do i=1,n
         VmI = V_mag(i)-I_mag(i)
         E_mag_mod(i) = a(1)+a(2)*I_mag(i)+a(3)*VmI
     &                 +a(4)*VmI**2+a(5)*VmI**3
         var = (E_mag_mod(i)-E_mag(i))**2
         chi2 = ((E_mag_mod(i)-E_mag(i))/E_err(i))**2
         if(chi2.gt.chi2max) then
           chi2max = chi2
           imax = i
         endif
         chi2sum = chi2sum + chi2
         varsum = varsum + var
       enddo
       if(chi2max.gt.chi2thresh) then
         write(6,*) 'removing star',id(imax),'  chi^2 =',chi2max
         if(imax.ne.n) then
           do i = imax,n-1
             id(i) = id(i+1)
             E_mag(i) = E_mag(i+1)
             I_mag(i) = I_mag(i+1)
             V_mag(i) = V_mag(i+1)
             E_err(i) = E_err(i+1)
           enddo
         endif
         n = n - 1
         go to 25
       endif
       rms = sqrt(varsum/n)

       write(6,991) n,chi2sum,rms,(a(i), i=1,mfit+1)
       write(2,991) n,chi2sum,rms,(a(i), i=1,mfit+1)
       write(2,993)

       write(2,994) (id(i),E_mag(i),E_mag_mod(i),E_err(i),
     &             ((E_mag(i)-E_mag_mod(i))/E_err(i))**2,
     &               I_mag(i),V_mag(i),V_mag(i)-I_mag(i), i=1,n)
       write(2,'(a)') '#'

       stop
       end

       subroutine colmagfit(I_mag,V_mag,afunc,ma)

       real afunc(ma)
       real I_mag,V_mag

       VmI = V_mag - I_mag
       afunc(1) = 1.
       afunc(2) = I_mag
       afunc(3) = VmI
       afunc(4) = VmI**2
       afunc(5) = VmI**3

       return
       end

      SUBROUTINE LFIT_VI(I_mag,V_mag,Y,SIG,NDATA,A,MA,LISTA,MFIT,COVAR,
     * NCVM,CHISQ,FUNCS)
c     replaced X with I_mag and V_mag
      PARAMETER (MMAX=50)
      real I_mag(NDATA),V_mag(NDATA)
      DIMENSION Y(NDATA),SIG(NDATA),A(MA),LISTA(MA),
     *    COVAR(NCVM,NCVM),BETA(MMAX),AFUNC(MMAX)
      external FUNCS

      KK=MFIT+1
      DO 12 J=1,MA 
        IHIT=0
        DO 11 K=1,MFIT
          IF (LISTA(K).EQ.J) IHIT=IHIT+1
11      CONTINUE
        IF (IHIT.EQ.0) THEN
          LISTA(KK)=J
          KK=KK+1
        ELSE IF (IHIT.GT.1) THEN
          write(6,*) 'Improper set in LISTA:'
          write(6,*) lista
          stop
        ENDIF
12    CONTINUE
      IF (KK.NE.(MA+1)) then
        write(6,*) 'Improper set in LISTA (2):'
        write(6,*) lista
        stop
      endif
      DO 14 J=1,MFIT
        DO 13 K=1,MFIT
          COVAR(J,K)=0.
13      CONTINUE
        BETA(J)=0.
14    CONTINUE
      DO 18 I=1,NDATA
        CALL FUNCS(I_mag(I),V_mag(I),AFUNC,MA)
        YM=Y(I)
        IF(MFIT.LT.MA) THEN
          DO 15 J=MFIT+1,MA
            YM=YM-A(LISTA(J))*AFUNC(LISTA(J))
15        CONTINUE
        ENDIF
        SIG2I=1./SIG(I)**2
        DO 17 J=1,MFIT
          WT=AFUNC(LISTA(J))*SIG2I
          DO 16 K=1,J
            COVAR(J,K)=COVAR(J,K)+WT*AFUNC(LISTA(K))
16        CONTINUE
          BETA(J)=BETA(J)+YM*WT
17      CONTINUE
18    CONTINUE
      IF (MFIT.GT.1) THEN
        DO 21 J=2,MFIT
          DO 19 K=1,J-1
            COVAR(K,J)=COVAR(J,K)
19        CONTINUE
21      CONTINUE
      ENDIF
      CALL GAUSSJ(COVAR,MFIT,NCVM,BETA,1,1)
      DO 22 J=1,MFIT
        A(LISTA(J))=BETA(J)
22    CONTINUE
      CHISQ=0.
      DO 24 I=1,NDATA
        CALL FUNCS(I_mag(I),V_mag(I),AFUNC,MA)
        SUM=0.
        DO 23 J=1,MA
          SUM=SUM+A(J)*AFUNC(J)
23      CONTINUE
        CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
24    CONTINUE
      CALL COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      RETURN
      END

      SUBROUTINE COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      DIMENSION COVAR(NCVM,NCVM),LISTA(MFIT)
      DO 12 J=1,MA-1
        DO 11 I=J+1,MA
          COVAR(I,J)=0.
11      CONTINUE
12    CONTINUE
      DO 14 I=1,MFIT-1
        DO 13 J=I+1,MFIT
          IF(LISTA(J).GT.LISTA(I)) THEN
            COVAR(LISTA(J),LISTA(I))=COVAR(I,J)
          ELSE
            COVAR(LISTA(I),LISTA(J))=COVAR(I,J)
          ENDIF
13      CONTINUE
14    CONTINUE
      SWAP=COVAR(1,1)
      DO 15 J=1,MA
        COVAR(1,J)=COVAR(J,J)
        COVAR(J,J)=0.
15    CONTINUE
      COVAR(LISTA(1),LISTA(1))=SWAP
      DO 16 J=2,MFIT
        COVAR(LISTA(J),LISTA(J))=COVAR(1,J)
16    CONTINUE
      DO 18 J=2,MA
        DO 17 I=1,J-1
          COVAR(I,J)=COVAR(J,I)
17      CONTINUE
18    CONTINUE
      RETURN
      END

      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      PARAMETER (NMAX=50)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                stop 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) stop 'Singular matrix.'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
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

