      implicit none
     
      call readem
      call interpem

      stop
      end
 
**************************************************************************************************

      subroutine readem
c
c Read the data file allem.d containing emissivities on a grid of densities and temperatures
c log10(Te)=3.0(0.1)4.5,    log10(Ne)=2.0(0.5)6.0
c
      implicit none

      integer mxd,mxt,maxtr
      parameter (mxd=9, mxt=16,maxtr=10000)
      integer mult,nlo,nup,ind,nd,nt,id,k,Llo,Lup,is,dell,ntr,ik
      integer lineid(maxtr),point(maxtr),status
      character*1 splo,spup,spec(6),kard
      real*8 em(maxtr,mxd,mxt),tl(mxt),dl(mxd),emx(mxt),
     *       alam(maxtr),alamx

      common/emdata/em,tl,dl,nd,nt
      common/linedata/alam,lineid,point,ntr

      data tl/3.0,3.1, 3.2,3.3,3.4,3.5,
     *     3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5/
      data dl/2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0/
      data spec/'S','P','D','F','G','H'/

      nd=9
      nt=16
      ntr=1032

      status=access('allem.d','')

      if(status.ne.0) then
        write(6,*) 'File allem.d is required and does not exist'
        stop
      endif

      open(unit=10,status='old',file='allem.d')

      do id=1,nd

        do ik=1,9
          read(10,*) kard
        enddo           

        do ik=1,ntr

          read(10,10) mult,nlo,splo,nup,spup,alamx,
     *                        (emx(k),k=1,nt)
          is=mult/2
          do k=1,6
            if(splo.eq.spec(k)) Llo=k-1
            if(spup.eq.spec(k)) Lup=k-1
          enddo
          dell=(Llo-Lup+1)/2
          if(id.eq.1) then
            ind=300*nup+50*nlo+10*Llo+2*dell+is
            alam(ik)=alamx
            lineid(ik)=ind
          endif
          do k=1,nt
            em(ik,id,k)=log10(emx(k))
          enddo

        enddo
      enddo

      return

   10 format(i2,i3,1x,a1,i4,1x,a1,11x,f11.2,1x,20e10.3)
   11 format(i5,i2,i3,1x,a1,i4,1x,a1,f11.2,11x,1x,1p20e10.3)

      end

*****************************************************************************************

      subroutine linelists

      implicit none
      integer maxtr
      parameter(maxtr=10000)
      character*1 spec(6)
      integer lineid(maxtr),point(maxtr),used(maxtr),lineidk,
     *        ntr,itr,jtr,nup,nlo,Lup,Llo,is,k,dell
      real*8 alam(maxtr),alamin

      common/linedata/alam,lineid,point,ntr

      data spec/'S','P','D','F','G','H'/

c
c  Sort lines by increasing wavelength
c
      do itr=1,ntr
         used(itr)=0
      enddo
      do itr=1,ntr
         alamin=1.d10
         do jtr=1,ntr
            if(alam(jtr).lt.alamin.and.used(jtr).eq.0) then
              alamin=alam(jtr)
              point(itr)=jtr
            endif
         enddo
         used(point(itr))=1
      enddo

c
c  Unpack quantum numbers and create wavelength ordered line list
c
      open(unit=10,status='unknown',file='HeIwlist')

      write(10,101) 
      do itr=1,ntr

         k=lineid(point(itr))
         nup=k/300
         k=k-300*nup
         nlo=k/50
         k=k-50*nlo
         Llo=k/10
         k=k-10*Llo
         dell=k/2
         Lup=Llo+1-2*dell
         is=k-2*dell
         is=2*is+1
         write(10,100) itr,alam(point(itr)),nup,is,spec(Lup+1),
     *       nlo,is,spec(Llo+1)

      enddo

 100  format(1x,i5,f15.2,5x,i3,1x,i1,a1,1x,'--',1x,i3,1x,i1,a1,i15)
 101  format(1x,'Wavelength table'/1x,'----------------'/
     *    1x,'Index',5x,'Wavelength',7x,'n SL',6x,'n SL'/
     *   12x,'Angstrom'/)
      return
      end

*******************************************************************************************
 
      subroutine interpem

c
c  Interpolate in emissivities for transition number itr at user density and temperature, ud, ut
c  Result is uem
c

      implicit none

      integer mxd,mxt,maxtr,iopt,maxout
      parameter (mxd=9, mxt=16,maxtr=10000,maxout=100)

      integer lineid(maxtr),point(maxtr),uind(maxout),ind,nd,nt,id,k,
     *        ntr,itr,status,nutr,iutr,iuem3
      real*8 em(maxtr,mxd,mxt),tl(mxt),dl(mxd),ut,ud,uem,TwoDnt
      real*8 alam(maxtr),uwave(maxout),ulam,udl,utl,uem1,uem2,uem3,
     *       refem

      common/emdata/em,tl,dl,nd,nt
      common/linedata/alam,lineid,point,ntr

      write(6,10)
 10   format(//'Program to provide emissivities for helium '
     *   ' recombination lines'/
     *   ' in Case B (Del Zanna and Storey, 2022, MNRAS)' /
     *   'Emissivities  are available for all transitions nu->nl,'
     *   ' where nu=2-25, nl=2-5',/
     *   'Allowed density (Ne[/cm**3]) range, log10(Ne)=2.0-6.0',/
     *   'Allowed temperature (Te [K]) range, log10(Te)=3.0-4.5',//)     

 1    write(6,11)  
 11   format('Choose one of three options' //
     *  '1: Generate a list of available transitions and  wavelengths.',
     *  / '   The file HeIwlist will be created. The listed wavelength,', 
     *  / '   to 2 d.p., is used to identify a particular transition' //
     *  '2: Specify a  transition by wavelength as listed (to 2 d.p.)',
     *  / '   in the file Hewlist. e.g. 4471.50', 
     *  / '   You will be prompted for a density and temperature'//
     *  '3: Provide a file containing a list of transition',
     *  / '   wavelengths as listed (to 2 d.p.) in the file',
     *  / '   Hewdata, one wavelength per line.',
     *  / '   You will be prompted for a density and temperature'//  
     *  '4: Exit' //
     *  'Enter an integer')

      read(5,*) iopt
      if(iopt.lt.1.or.iopt.gt.4) go to 1

      if(iopt.eq.4) stop

c
c Option 1
c
      if(iopt.eq.1) then

        call linelists

        go to 1
c
c Option 2
c
      elseif(iopt.eq.2) then

 20      write(6,210)
 210     format(1x,'Enter the transition wavelength in'
     *              ' Angstrom to 2 d.p.')
        read(5,*) ulam
        do itr=1,ntr
           if(abs(alam(itr)-ulam).lt.0.01) go to 21
        enddo
        write(6,*) 'Wavelength not found'
        go to 1
 21     write(6,211) 
 211    format(1x,'Enter a density (per cm**3) and temperature (K)')
        read(5,*) ud,ut

        udl=log10(ud)
        utl=log10(ut)
        refem=TwoDnt(610,udl,utl,3)             
        uem1=TwoDnt(itr,udl,utl,1)
        uem2=TwoDnt(itr,udl,utl,2)
        uem3=TwoDnt(itr,udl,utl,3)
c Correct "round-to-even" feature in gftn
        iuem3=log10(uem3)-6
        uem3=uem3+10.d0**iuem3
        refem=uem3/refem
        write(6,212) ud,ut,uem3,refem,refem
        
 212    format(//'Supplied density and temperature ',2f8.0,//
     *  'Emissivity (erg.cm**3/s) with 6 point interpolation',/ 
     *  'in both density and temperature:',10x,1pe10.2,/
     *  'Relative to 4471: ',0pf10.3,' = ',1pe10.3//)

        write(6,213) 
 213    format('1: New density and temperature'//
     *         '2: New transition '//
     *         '3: Return to main menu'//
     *         'Enter an integer'//)
        read(5,*) iopt
        if(iopt.ne.1.and.iopt.ne.2.and.iopt.ne.3) go to 1
        if(iopt.eq.1) go to 21
        if(iopt.eq.2) go to 20
        if(iopt.eq.3) go to 1
        go to 1
c
c  Option 3
c
      elseif(iopt.eq.3) then

        write(6,310)
 310    format('List of HeI wavelengths is required in file',
     :     ' HeIwdata'/
     :     'First record contains number of lines (integer)'/
     :     'One wavelength per line (to two d.p.)'/
     :     'Output in HeIemdata'//)

        status=access('HeIwdata','')

        if(status.ne.0) then
          write(6,311)
 311      format('File HeIwdata does not exist'//)
          go to 1
        endif

        open(unit=15,status='unknown',file='HeIwdata')
        open(unit=20,status='unknown',file='HeIemdata')

        rewind(15)

* read line wavelengths
* note that lines are identified by their wavelength as 
* printed in the Hewlist. file and given to 2 d.p.

        nutr=1
 30     read(15,*,end=32) uwave(nutr)
        do itr=1,ntr
          if(abs(alam(itr)-uwave(nutr)).lt.0.01) then
            uind(nutr)=itr
            go to 31
          endif
        enddo
        write(6,313) nutr,uwave(nutr)         
        go to 1

 31     write(6,312) nutr,uwave(nutr)

 312    format(1x'line number',i5,' wavelength [A] ',f12.2,'  found')
 313    format(1x'line number',i5,' wavelength [A] ',f12.2,' not found')

        nutr=nutr+1

        if(nutr.gt.maxout) then
          write(6,314) nutr,maxout
          go to 1
        endif
 314    format(1x,' Number of transitions ',i5,' exceeds dimension',i5)
        go to 30

 32     nutr=nutr-1
        write(6,*) 'number of lines = ',nutr

        write(6,211) 
        read(5,*) ud,ut

        udl=log10(ud)
        utl=log10(ut)
        
        write(6,315) ud,ut
        write(20,315) ud,ut
 315    format(1x,'Supplied density and temperature ',2f8.0//
     *  'Emissivities (erg.cm**3/s) with 6 point interpolation',
     *  /'in both density and temperature:'/
     *   'Final column is emissivity relative to 4471A line'//)

        do iutr=1,nutr

           itr=uind(iutr)
           refem=TwoDnt(610,udl,utl,3)             
           uem1=TwoDnt(itr,udl,utl,1)
           uem2=TwoDnt(itr,udl,utl,2)
           uem3=TwoDnt(itr,udl,utl,3)
c Correct "round-to-even" feature in gftn
           iuem3=log10(uem3)-6
           uem3=uem3+10.d0**iuem3
           refem=uem3/refem
           write(6,316) iutr,uwave(iutr),uem3,refem
           write(20,316) iutr,uwave(iutr),uem3,refem
 316       format(1x'line number',i5,' wavelength [A] ',f12.2,1p2e12.2)

        enddo

        write(6,317) 
 317    format(//'1: New density and temperature'//
     *           '2: Return to main menu'//
     *           'Enter an integer'//)
        read(5,*) iopt
        if(iopt.ne.1.and.iopt.ne.2) go to 1
        if(iopt.eq.1) go to 32
        if(iopt.eq.2) go to 1

      endif

      return

      end

*******************************************************************************


      double precision function OneDt(itr,in,ut,np)

      implicit none

      integer mxt,mxd,mxout,maxtr
      parameter (maxtr=10000,mxout=60,mxt=16,mxd=9)

      real*8 tl(mxt),dl(mxd),em(maxtr,mxd,mxt),alam(maxtr)
      real*8 ft(mxt),xt(mxt),fvt,ut
      integer lineid(maxtr),point(maxtr),itr,in,it,nt,nd,ier,np,ntr

      common/emdata/em,tl,dl,nd,nt
      common/linedata/alam,lineid,point,ntr

      if(ut.lt.2.6.or.ut.gt.4.5) then
         write(6,*) 'log10(temperature)=',ut,' outside allowable',
     *             ' range of 2.6-4.5'
         stop
      endif 

      do it=1,nt
         ft(it)=em(itr,in,it)
         xt(it)=tl(it)
      enddo

      call lagr(ft,xt,fvt,ut,np,nt,mxt,ier)
      OneDt=fvt

      return

      end

*********************************************************************************

      double precision function OneDn(itr,ud,ut,np)

      implicit none

      integer mxt,mxd,mxout,maxtr
      parameter (maxtr=10000,mxout=60,mxt=16,mxd=9)

      real*8 tl(mxt),dl(mxd),fn(mxd),an(mxd),
     :       em(maxtr,mxd,mxt)
      real*8 fvd,ud,ut,OneDt
      integer itr,np,nt,nd,in,ier

      common/emdata/em,tl,dl,nd,nt

      if(ud.lt.2.0.or.ud.gt.6.0) then
         write(6,*) 'log10(density)=',ud,' outside allowable',
     *              ' range of 2.0-6.0'
         stop
      endif 

      do in=1,nd
         fn(in)=OneDt(itr,in,ut,np)
         an(in)=dl(in)
      enddo


      call lagr(fn,an,fvd,ud,np,nd,mxd,ier)
      OneDn=fvd

      return

      end

********************************************************************************************

      double precision function TwoDnt(itr,ud,ut,np)

      implicit none

      real*8 ud,ut,OneDn
      integer itr,np
      
      TwoDnt=10.d0**OneDn(itr,ud,ut,np)

      return
      end

***********************************************************************************************

      SUBROUTINE LAGR(F,X,FV,XV,NP,NI,ND,IER)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  SUBROUTINE TO PERFORM LAGRANGIAN INTERPOLATION OF A SET OF
C  TABULATED FUNCTION VALUES F(X(I))
C
C  INPUTS
C  F(I)  = ARRAY OF TABULATED FUNCTION VALUES
C  X(I)  = ARRAY OF CORRESPONDING X VALUES
C  NI    = NUMBER OF TABULATED VALUES
C  2*NP  = NUMBER OF POINTS TO BE USED IN INTERPOLATION
C  XV    = VALUE OF X AT WHICH FUNCTION VALUE IS REQUIRED
C
C  OUTPUTS
C  FV    = INTERPOLATED VALUE
C  IER   = ERROR CODE (.NE.0 IF XV OUTSIDE RANGE OF X)
C
      DIMENSION F(ND),X(ND),XIQ(20)

      LG=2*NP
      IER=0
      IF(XV.LT.X(1).OR.XV.GT.X(NI)) THEN
        IER=1
        RETURN
      ENDIF
C
C  LOCATE XV IN X
      DO 5999 I=2,NI
        IF(X(I).GE.XV) THEN
C
C  SELECT SUBSET OF X FOR INTERPOLATION, TAKING CARE AT ENDS
          IN=NI-LG+1
          ILOW=I-NP
          IF(ILOW.LT.1) ILOW=1
          IF(ILOW.GT.IN) ILOW=IN
          DO 5998 ITAU=1,LG
            XIQ(ITAU)=X(ILOW+ITAU-1)
 5998     CONTINUE
C
C  CONSTRUCT INTERPOLATION COEFFICIENTS
          FV=0.D0
          DO 5997 ITAU=1,LG
            FL=PHI(XIQ,LG,ITAU,XV)/PHI(XIQ,LG,ITAU,XIQ(ITAU))
            FV=FV+FL*F(ILOW+ITAU-1)
 5997     CONTINUE
          RETURN
        ENDIF
 5999 CONTINUE
      END

*******************************************************************************

      DOUBLE PRECISION FUNCTION PHI(XIQ,LG,ITAU,X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XIQ(20)
      FK=1.D0
      DO 5999 L=1,LG
        IF(L.NE.ITAU) THEN
          FK=FK*(X-XIQ(L))
        ENDIF
 5999 CONTINUE
      PHI=FK
      RETURN
      END
      
******************************************************************************

