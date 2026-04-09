       program HST_Keck_coord_trans

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       parameter(nbmax=30,nsmax=30)
       real*8 x(nsmax),xe(nsmax),y(nsmax),ye(nsmax)
       real*8 x0(nsmax),y0(nsmax)
       real*8 xe0(nsmax),ye0(nsmax)
       real*8 xenm(nsmax),yenm(nsmax)
       real*8 amag(nsmax),ain(nsmax)
       real*8 amag0(nsmax),ain0(nsmax)
       real*8 bmag0(nsmax),bmag(nsmax)
       real*8 duma0(nsmax),duma(nsmax)
       real*8 dumb0(nsmax),dumb(nsmax)
       real*8 tr1(2),tr2(2),delxfit(2)
       integer iatype(nsmax),iseq(nsmax),iti(nsmax)
       integer iatype0(nsmax),ibtype0(nsmax),ibtype(nsmax)
       integer ityp_fid_I(nsmax),ityp_fid_V(nsmax),it(nsmax)
       real*8 fid_x(nsmax),fid_y(nsmax),fid_in(nsmax)
       real*8 fida_x(nsmax),fida_y(nsmax)
       real*8 fid_V(nsmax),fid_I(nsmax)
       integer ityp_fid_I0(nsmax),ityp_fid_V0(nsmax)
       real*8 fid_x0(nsmax),fid_y0(nsmax),fid_in0(nsmax)
       real*8 fid_V0(nsmax),fid_I0(nsmax)
       real*8 hst_mag(nsmax),hst_msig(nsmax)
       real*8 xnew(nsmax),ynew(nsmax)
       real*8 x_keck(nsmax),y_keck(nsmax),x_hst(nsmax),y_hst(nsmax)
       real*8 fK_keck(nsmax),fK_sig(nsmax)
       real*8 x_hst0(nsmax),y_hst0(nsmax),hst_mag0(nsmax)
       real*8 hst_msig0(nsmax),x_keck0(nsmax),y_keck0(nsmax)
       real*8 fK_keck0(nsmax),fK_sig0(nsmax)
       integer kd0(nsmax)
       integer ibi(nbmax),jbi(nbmax),imat(nbmax),i0(nsmax),kd(nsmax)
       character*120 inline
       character*80 infile
       double precision ra_arcsec,dec_arcsec,ra_t_arcsec,dec_t_arcsec
       double precision dec_deg,cos_dec,dec_t_deg,cos_t_dec
       double precision ra_deg
       double precision rah,ram,ras,decd,decm,decs,rad
       double precision Io_t

c      hardwired parameters
c      --------------------
       scale_rat = 0.371/0.26   ! CTIO I-band/OGLE I-band
       as_pix_o = 0.26
       as_pix_m = 0.371
       mx_nmatch = 200
       mn_match = 20
       del_x_mx = 100.
       del_y_mx = 100.
       ratin_mx = 9.5
       del_x_tol = 20.0
       del_y_tol = 20.0
       theta_tol = 0.1
       r_tol = 3.0
       nfast = 901
       nslow = 901
       err_min = 0.02
       tolout = 3.0
       chi_thresh = 25.0
       nbtotmin = 4

       err_min2 = err_min**2

c      open and read matched star file
c      -------------------------------
       write(6,*) 'enter input file name:'
       read(5,'(a)') infile
       write(6,*) 'enter chi2cut, err_min:'
       read(5,*) chi2cut, err_min
       open(unit=1,file=infile,status='old')
       nf = 0
       read(1,*)
       read(1,*)
       itarg = 0
       nt = 0
       nm = 0
       do j = 1,99999
         read(1,*,end=11,err=11) xh,yh,hmag,xherr,yherr,hmag_sig,kid,
     &             xk,yk,fmagk,fmsigk,itarg
         if(itarg.ne.1) then
           nm = nm + 1
           x_hst(nm) = xh
           y_hst(nm) = yh
           hst_mag(nm) = hmag
           xenm(nm) = xherr
           yenm(nm) = yherr
           hst_msig(nm) = hmag_sig
           kd(nm) = kid
           x_keck(nm) = xk
           y_keck(nm) = yk
           fK_keck(nm) = fmagk
           fK_sig(nm) = fmsigk
          x_hst0(nm) = x_hst(nm)
          y_hst0(nm) = y_hst(nm)
          hst_mag0(nm) = hst_mag(nm)
          hst_msig0(nm) = hst_msig(nm)
          kd0(nm) = kd(nm)
          x_keck0(nm) = x_keck(nm)
          y_keck0(nm) = y_keck(nm)
          fK_keck0(nm) = fK_keck(nm)
          fK_sig0(nm) = fK_sig(nm) 
         endif
       enddo
 11    continue
       nmt = nm
       nm = nm - nt
       nm0 = nm
       write(6,*) 'read in',nm,' matched star coordinates'
       write(6,*) '    and',nt,' target star coordinates'
       if(nm.gt.nsmax) then
         stop 'nsmax too small'
       endif
       do i = 1,nm
         xe(i) = sqrt(xenm(i)**2+err_min**2)
         ye(i) = sqrt(yenm(i)**2+err_min**2)
         ibi(i) = i
         jbi(i) = i
         i0(i) = i
       enddo

 22    continue

c      fit the precise transformation
c      ------------------------------
       CALL fldtrans(nbmax,nm,ibi,jbi,nsmax,x_hst,xe,y_hst,ye,
     &               x_keck,y_keck,tr1,tr2,delxfit,tolout,nbtotmin,ier)

c   write out transformation parameters
c   -----------------------------------
       write(6,*) 'initial transformation parameters:'
       write(6,*)' TR11=',tr1(1),' TR12=',tr2(1),' TR21=',tr1(2),
     &      ' TR22=',tr2(2),' delX=',delxfit(1),' delY=',delxfit(2)

c      determine chi^2 and worst data point
c      ------------------------------------
       chi2max = 0.d0
       chi2sum = 0.d0
       rms_x = 0.d0
       rms_y = 0.d0
       ichi2 = 0
       do i = 1,nm
         xnew(i)=tr1(1)*x_keck(i)+tr2(1)*y_keck(i)+delxfit(1)
         ynew(i)=tr1(2)*x_keck(i)+tr2(2)*y_keck(i)+delxfit(2)
         chi2 = ((xnew(i)-x_hst(i))/xe(i))**2
     &        + ((ynew(i)-y_hst(i))/ye(i))**2
         write(6,*) 'i, chi^2 =',i0(i),chi2
         chi2sum = chi2sum + chi2
         if(chi2.gt.chi2max.and.chi2.gt.chi2cut) then
           chi2max = chi2
           ichi2 = i
         endif
         rms_x = rms_x + (xnew(i)-x_hst(i))**2
         rms_y = rms_y + (ynew(i)-y_hst(i))**2
       enddo
       rms_x = sqrt(rms_x/nm)
       rms_y = sqrt(rms_y/nm)

       write(6,991) chi2sum,rms_x,rms_y,nm
c      remove worst data point
c      -----------------------
       if(ichi2.gt.0) then
         write(6,992) i0(ichi2),x_hst(ichi2),y_hst(ichi2),chi2max
         do i = ichi2,nm-1
           x_hst(i) = x_hst(i+1)
           y_hst(i) = y_hst(i+1)
           xe(i) = xe(i+1)
           ye(i) = ye(i+1)
           x_keck(i) = x_keck(i+1)
           y_keck(i) = y_keck(i+1)
           fK_keck(i) = fK_keck(i+1)
           fK_sig(i) = fK_sig(i+1)
           hst_mag(i) = hst_mag(i+1)
           hst_msig(i) = hst_msig(i+1)
           kd(i) = kd(i+1)
           i0(i) = i0(i+1)
         enddo
         nm = nm - 1
         go to 22
       endif

c      write out final transformed coordinates
c      ---------------------------------------
       write(6,993)
       do i = 1,nmt
         xnew(i)=tr1(1)*x_keck0(i)+tr2(1)*y_keck0(i)+delxfit(1)
         ynew(i)=tr1(2)*x_keck0(i)+tr2(2)*y_keck0(i)+delxfit(2)
         dx = xnew(i) - x_hst0(i)
         dy = ynew(i) - y_hst0(i)
         if(i.le.nm.or.i.ge.nm0) then
           write(6,994) i,xnew(i),x_hst0(i),ynew(i),y_hst0(i),dx,dy,
     &       xenm(i),yenm(i),hst_mag0(i),hst_msig0(i),kd0(i),
     &       fK_keck0(i),fK_sig0(i)
         endif
       enddo

 991   format('chi^2 = ',f11.3,'  RMS_x =',f10.6,'  RMS_y =',f10.6,
     &        '   for',i4,' stars')
 992   format(' removing point',i4,' at',2f9.4,'  with chi^2 =',f11.3)
 993   format('#    x_tran    x_hst     y_tran    y_hst     dx      ',
     &     'dy      x_sig   y_sig    hst_mag hst_sig K#  K-mag  K-sig')
 994   format(i3,4f10.4,4f8.4,f9.4,f8.4,i4,2f7.3)

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
