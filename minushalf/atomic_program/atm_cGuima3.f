c
      subroutine run()
cGuima Nesta versao se le o numero 'maxit' maximo de interaccoes, logo depois
cGuima dos orbitais de valencia. Ver os input INP.pg e INP.pt. Estes arquivos
cGuima devem ser renomeados par INP.
c
cGuima Adicionada leitura de potencial que se adiciona ao pseudopotencial
cGuima gerado. O potencial a adicionar est� no arquivo 'adiciona'.
cGuima A instruccao      inquire(file='adiciona',exist=lexist)
cGuima                   if(lexist)  open(unit=21,file='adiciona')
cGuima verifica se o arquivo existe. Se existir o arquivo � aberto e lido.
c
cGuima Cria arquivo 'VTOTAL' com o potencial da Eq. de Sch. ou Dirac.
cGuima Cria arquivo 'psfun.Guima' com as funccoes de onda ae, pg e pt.

cGuima Calcula tambem os valores das medias <r**2>, <r**4> e auto-energia
cGuima eletrost�tica para os orbitais de valencia.
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'charge.h'
      include 'elecpot.h'
      include 'energy.h'
c
C     .. Parameters ..
      double precision tol
      parameter (tol=1.D-8)
      double precision zero, one
      parameter (zero=0.D0,one=1.D0)
C     ..
C     .. Local Scalars ..
      double precision aa, dv, dvmax, t1, t2, xmixo, zsold
      integer i, icon2, iconv, iiter, iter, 
     &        itsm, maxit, nconf, jobold
      character icold*2, naold*2
      character pot_id*40, headline*79
C     ..
C     .. Local Arrays ..
      double precision econf(100)
      integer nn(norbmx)
C     ..
C     .. Arrays in Common ..
      double precision vn1d(nrmax), 
     &                 vn11d(nrmax), vn2d(nrmax), vn22d(nrmax),
     &                 vn1u(nrmax), vn11u(nrmax), vn2u(nrmax),
     &                 vn22u(nrmax), wk1(nrmax), wk2(nrmax)
C     ..
      common /mixer/ vn1d, vn11d, vn2d, vn22d,
     &               vn1u, vn11u, vn2u, vn22u,
     &               wk1, wk2
c
C     .. External Subroutines ..
      external hsc, ker, tm2
      double precision second
      external dmixp, dsolv1, dsolv2, etotal, ext, input,
     &         prdiff, pseudo, velect, vionic, second, denplot
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, exp, log, sqrt
C     ..
c
c  Initialize the timer.
c
      t1 = second()
cGuima inserccao
      open(unit=27,file='psfun.guima')
cGuima fim da inserccao
c
c  Startup values for doing multiple input data sets.
c
      naold = '  '
      icold = '  '
      zsold = zero
      jobold = -1
      nconf = 0
      nr = nrmax
c
c      open files
c
      open(unit=5,file='INP',status='old',form='formatted')
      open(unit=6,file='OUT',status='unknown',form='formatted')
c
c   Print version
c
      call prversion
c
c
c   Start of loop over configurations or different jobs.
c   Note that it is possible to mix different kinds of
c   jobs (for example, an all-electron calculation followed
c   by a pseudopotential test).
c
   10 continue
      norb = norbmx
c
c   Read the input data.
c   New: allow directives before each job
c
      call check_directives(5)
cGuima modificado
      Call Input(maxit)
      close(5)
cGuima termina modificaccao
c
c     Prevent computing energy differences between different
c     kinds of jobs. Print those differences only at
c     the end of a series.
c
      if (nconf .gt. 0) then
         if ((job .ne. jobold) .or.
     $       (naold .ne. nameat) .or.
     $       (icold .ne. icorr) .or.
     $       (zsold .ne. zsh))    then

            if ((job.lt.1.or.job.gt.3) .and.
     $           (nconf .ge. 2)) then
               call prdiff(nconf,econf,jobold)
            endif
            nconf = 0
         endif
      endif

      if (job .lt. 0) then
c  Stop - no more data in input file,
         t2 = second()
         write(6,9000) t2 - t1
 9000    format(//' The total time for the calculation is ',f12.5,
     &         ' seconds')
         return
      end if
c
c     Print out info gathered by Input
c
      call Header
c
c
c  Jump charge density 'set up' and ionic data input if
c  configuration test.
c
c  All-electron calculations must have a properly
c  setup starting charge.
c  Make sure it is done even if we change jobs
c  to "ae".
c
      itsm = znuc/9 + 3
c
      if (zsold .eq. zsh .and.
     $    naold .eq. nameat  .and.
     $    jobold .eq. job    .and.
     &    (job.lt.1.or.job.gt.4)) then

c        do nothing

      else

         if (job .lt. 4) then

            write(6,'(a)') 'Setting up initial charge'
c
c           Set up the initial charge density.
c           cdd and cdu  =  (4 pi r**2 rho(r))/2

c           The charge density setup (aa) is scaled with
c           an empirical formula that reduces the
c           number of iterations needed for the screening
c           potential convergence.

            aa = sqrt(sqrt(znuc))/2 + one
            do 20 i = 1, nr
               cdd(i) = zel*aa**3*exp(-aa*r(i))*r(i)**2/4
               cdu(i) = cdd(i)
   20       continue

         end if
c
c        set up ionic potentials
c
         write(6,'(a)') 'Setting up ionic potential'
         call Vionic
c
      end if
c
c   Set up the electronic potential.
c
      call velect(0,0,ispp,zel)
c
      do 30 i = 1, nr
         vid(i) = vod(i)
         viu(i) = vou(i)
   30 continue
c
c   Start the iteration loop for electronic convergence.
c
      iconv = 0
      icon2 = 0
cGuima modificado
c      maxit = 100
cGuima termina modificaccao
c
c    The screening potential mixing parameter is
c    an empirical function of the nuclear charge.
c    Larger atoms require a slower convergence
c    then smaller atoms.
c
      xmixo = one/log(znuc+7*one)
c
      do 50 iter = 1, maxit
c
         if (iter .eq. maxit) iconv = 1
c
c  compute orbitals
c     
cag      Use another array to avoid passing no directly (it is
cag      in a common block shared by atm and dsolvX)
cag
         do i=1,norb
            nn(i) = no(i)
         enddo
cag
         if (icon2 .lt. 2) then
            call dsolv1(1,norb,nn)
         else
            call dsolv2(iter,iconv,ispp,1,norb,ncore,nn)
         end if
c
c  set up output electronic potential from charge density
c
         call velect(iter,iconv,ispp,zel)
c
c  check for convergence
c
         if (iconv .gt. 0) go to 60
         dvmax = zero
         do 40 i = 1, nr
            dv = (vod(i)-vid(i))/(1.D0+vod(i)+vou(i))
            if (abs(dv) .gt. dvmax) dvmax = abs(dv)
            dv = (vou(i)-viu(i))/(1.D0+vou(i)+vod(i))
            if (abs(dv) .gt. dvmax) dvmax = abs(dv)
   40    continue
         icon2 = icon2 + 1
         if (dvmax .le. tol) iconv = 1
c
c  Mix the input and output electronic potentials.
c
c    The screening potential is initially mixed with a
c    percentage of old and new for itsm iterations.
c    This brings the potential to a stable region
c    after which an Anderson's extrapolation scheme
c    is used.
c
         if (iter .lt. itsm) then
            iiter = 2
         else
            iiter = iter - itsm + 3
         end if
c
        call dmixp(vod,vid,xmixo,iiter,3,nr,wk1,wk2,
     1             vn1d,vn11d,vn2d,vn22d)
        call dmixp(vou,viu,xmixo,iiter,3,nr,wk1,wk2,
     1             vn1u,vn11u,vn2u,vn22u)
c
   50 continue
c
c   End of iteration of electronic convergence loop.
c
      write(6,9010) dvmax, xmixo
 9010 format(/' potential not converged - dvmax =',d10.4,' xmixo =',
     &      f5.3)
      call ext(1)
c
   60 continue
      write(6,9020) icon2
 9020 format(/' Total number of iterations needed for',
     &      ' electron screening potential is ',i2,/)
c
c  Find the total energy.
c
      call etotal(1,norb)
c
c  Plot the charge density.
c
      call denplot
c
c   Compute the logarithmic derivative as a function of energy 
c
      if (logder_radius .gt. 0.d0) call logder(ncore+1,norb,'AE')
c
c   Replace the valence charge density. Kludgeish...
c
c      if (job .eq. 5) then
c         job = 6
c         call Vionic
c         job = 5
c      endif
c
c     Better way...
c
      if (job .eq. 5) then
         call change_valence
      endif
c
c  Pseudopotential generation.
c
      if (job .ge. 1 .and. job .le. 3) then
c
ctwb
        ifcore = job - 1
ctwb
        if (scheme .eq. 0) then 
c HSC
          pot_id = 'Hamann, Schluter, and Chiang'
          write(headline,8000)
 8000     format(' nl    s    eigenvalue',6x,'rc',10x,'cl',9x,
     &           'gamma',7x,'delta')
c
          call Pseudo(pot_id,headline,hsc)
c
        else if (scheme .eq. 1) then
c KER
          pot_id = 'Kerker'                        
          write(headline,8020)
 8020     format(' nl    s    eigenvalue',6x,'rc',
     &             4x,6x,'cdrc',7x,'delta')
c
          call Pseudo(pot_id,headline,ker)
c
        else if (scheme .eq. 6) then
c TM2
          pot_id = 'Troullier-Martins'                        
          write(headline,8040)
 8040     format(' nl    s    eigenvalue',6x,'rc',
     &             10x,'cdrc',7x,'delta')
c
          call Pseudo(pot_id,headline,tm2)
c
        else
              write(6,*) ' Only HSC, KER, TM2 and TM4 supported.'
              stop 'NOFLAV'
c
        end if
c

      end if
c
      nconf = nconf + 1
      econf(nconf) = etot(10)
c
c   End loop of configuration.
c
      naold = nameat
      icold = icorr
      zsold = zsh
      jobold = job

      go to 10
c
      end








c
c $Id: nucl_z.f,v 1.2 1997/05/22 17:32:21 wdpgaara Exp $
c
c $Log: nucl_z.f,v $
c Revision 1.2  1997/05/22 17:32:21  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:55  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      double precision function nucl_z(name)
c
c      function determines the nuclear charge of an element
c
c  njtj  ***  modifications  ***
c    All elements from H to Lr are included
c  njtj  ***  modifications  ***
c
C     .. Parameters ..
      double precision one
      parameter (one=1.D0)
C     ..
C     .. Scalar Arguments ..
      character name*2
C     ..
      double precision charge
C     .. External Subroutines ..
      external ext
C     ..
      if (name .eq. 'H ' .or. name .eq. ' H') then
         charge = 1*one
      else if (name .eq. 'He') then
         charge = 2*one
      else if (name .eq. 'Li') then
         charge = 3*one
      else if (name .eq. 'Be') then
         charge = 4*one
      else if (name .eq. 'B ' .or. name .eq. ' B') then
         charge = 5*one
      else if (name .eq. 'C ' .or. name .eq. ' C') then
         charge = 6*one
      else if (name .eq. 'N ' .or. name .eq. ' N') then
         charge = 7*one
      else if (name .eq. 'O ' .or. name .eq. ' O') then
         charge = 8*one
      else if (name .eq. 'F ' .or. name .eq. ' F') then
         charge = 9*one
      else if (name .eq. 'Ne') then
         charge = 10*one
      else if (name .eq. 'Na') then
         charge = 11*one
      else if (name .eq. 'Mg') then
         charge = 12*one
      else if (name .eq. 'Al') then
         charge = 13*one
      else if (name .eq. 'Si') then
         charge = 14*one
      else if (name .eq. 'P ' .or. name .eq. ' P') then
         charge = 15*one
      else if (name .eq. 'S ' .or. name .eq. ' S') then
         charge = 16*one
      else if (name .eq. 'Cl') then
         charge = 17*one
      else if (name .eq. 'Ar') then
         charge = 18*one
      else if (name .eq. 'K ' .or. name .eq. ' K') then
         charge = 19*one
      else if (name .eq. 'Ca') then
         charge = 20*one
      else if (name .eq. 'Sc') then
         charge = 21*one
      else if (name .eq. 'Ti') then
         charge = 22*one
      else if (name .eq. 'V ' .or. name .eq. ' V') then
         charge = 23*one
      else if (name .eq. 'Cr') then
         charge = 24*one
      else if (name .eq. 'Mn') then
         charge = 25*one
      else if (name .eq. 'Fe') then
         charge = 26*one
      else if (name .eq. 'Co') then
         charge = 27*one
      else if (name .eq. 'Ni') then
         charge = 28*one
      else if (name .eq. 'Cu') then
         charge = 29*one
      else if (name .eq. 'Zn') then
         charge = 30*one
      else if (name .eq. 'Ga') then
         charge = 31*one
      else if (name .eq. 'Ge') then
         charge = 32*one
      else if (name .eq. 'As') then
         charge = 33*one
      else if (name .eq. 'Se') then
         charge = 34*one
      else if (name .eq. 'Br') then
         charge = 35*one
      else if (name .eq. 'Kr') then
         charge = 36*one
      else if (name .eq. 'Rb') then
         charge = 37*one
      else if (name .eq. 'Sr') then
         charge = 38*one
      else if (name .eq. 'Y ' .or. name .eq. ' Y') then
         charge = 39*one
      else if (name .eq. 'Zr') then
         charge = 40*one
      else if (name .eq. 'Nb') then
         charge = 41*one
      else if (name .eq. 'Mo') then
         charge = 42*one
      else if (name .eq. 'Tc') then
         charge = 43*one
      else if (name .eq. 'Ru') then
         charge = 44*one
      else if (name .eq. 'Rh') then
         charge = 45*one
      else if (name .eq. 'Pd') then
         charge = 46*one
      else if (name .eq. 'Ag') then
         charge = 47*one
      else if (name .eq. 'Cd') then
         charge = 48*one
      else if (name .eq. 'In') then
         charge = 49*one
      else if (name .eq. 'Sn') then
         charge = 50*one
      else if (name .eq. 'Sb') then
         charge = 51*one
      else if (name .eq. 'Te') then
         charge = 52*one
      else if (name .eq. 'I ' .or. name .eq. ' I') then
         charge = 53*one
      else if (name .eq. 'Xe') then
         charge = 54*one
      else if (name .eq. 'Cs') then
         charge = 55*one
      else if (name .eq. 'Ba') then
         charge = 56*one
      else if (name .eq. 'La') then
         charge = 57*one
      else if (name .eq. 'Ce') then
         charge = 58*one
      else if (name .eq. 'Pr') then
         charge = 59*one
      else if (name .eq. 'Nd') then
         charge = 60*one
      else if (name .eq. 'Pm') then
         charge = 61*one
      else if (name .eq. 'Sm') then
         charge = 62*one
      else if (name .eq. 'Eu') then
         charge = 63*one
      else if (name .eq. 'Gd') then
         charge = 64*one
      else if (name .eq. 'Tb') then
         charge = 65*one
      else if (name .eq. 'Dy') then
         charge = 66*one
      else if (name .eq. 'Ho') then
         charge = 67*one
      else if (name .eq. 'Er') then
         charge = 68*one
      else if (name .eq. 'Tm') then
         charge = 69*one
      else if (name .eq. 'Yb') then
         charge = 70*one
      else if (name .eq. 'Lu') then
         charge = 71*one
      else if (name .eq. 'Hf') then
         charge = 72*one
      else if (name .eq. 'Ta') then
         charge = 73*one
      else if (name .eq. 'W ' .or. name .eq. ' W') then
         charge = 74*one
      else if (name .eq. 'Re') then
         charge = 75*one
      else if (name .eq. 'Os') then
         charge = 76*one
      else if (name .eq. 'Ir') then
         charge = 77*one
      else if (name .eq. 'Pt') then
         charge = 78*one
      else if (name .eq. 'Au') then
         charge = 79*one
      else if (name .eq. 'Hg') then
         charge = 80*one
      else if (name .eq. 'Tl') then
         charge = 81*one
      else if (name .eq. 'Pb') then
         charge = 82*one
      else if (name .eq. 'Bi') then
         charge = 83*one
      else if (name .eq. 'Po') then
         charge = 84*one
      else if (name .eq. 'At') then
         charge = 85*one
      else if (name .eq. 'Rn') then
         charge = 86*one
      else if (name .eq. 'Fr') then
         charge = 87*one
      else if (name .eq. 'Ra') then
         charge = 88*one
      else if (name .eq. 'Ac') then
         charge = 89*one
      else if (name .eq. 'Th') then
         charge = 90*one
      else if (name .eq. 'Pa') then
         charge = 91*one
      else if (name .eq. ' U' .or. name .eq. 'U ') then
         charge = 92*one
      else if (name .eq. 'Np') then
         charge = 93*one
      else if (name .eq. 'Pu') then
         charge = 94*one
      else if (name .eq. 'Am') then
         charge = 95*one
      else if (name .eq. 'Cm') then
         charge = 96*one
      else if (name .eq. 'Bk') then
         charge = 97*one
      else if (name .eq. 'Cf') then
         charge = 98*one
      else if (name .eq. 'Es') then
         charge = 99*one
      else if (name .eq. 'Fm') then
         charge = 100*one
      else if (name .eq. 'Md') then
         charge = 101*one
      else if (name .eq. 'No') then
         charge = 102*one
      else if (name .eq. 'Lr') then
         charge = 103*one
      else
         write(6,9000) name
 9000    format(//'element ',a2,' unknown')
         call ext(200)
      end if
c
      nucl_z = charge
c
      return
c
      end
c
c $Id: difnrl.f,v 1.3 1997/05/22 17:32:06 wdpgaara Exp $
c
c $Log: difnrl.f,v $
c Revision 1.3  1997/05/22 17:32:06  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.2  1994/02/18  01:26:11  garcia
c *** empty log message ***
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine difnrl(iter,iorb,v,ar,br,n,l,spin,eigv,iflag)
c
      implicit none
c
      include 'radial.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
c
c    difnrl integrates the Schroedinger equation
c    if finds the eigenvalue eigv, the wavefunction ar
c    and the derivative br = d(ar)/dr
c
c    iorb:       orbital number
c    n, l, spin: orbital's quantum numbers
c    eigv:       eigenvalue
c
c  njtj  ***  modifications  ***
c    This routine has had major modifications.  Some
c    of the data used inside the main loop has been
c    calculated outside the main loop to reduce the number
c    of operations(uses extra array space to gain speed)
c    The predictor-corrector functions have been put
c    into a array.
c    The iflag variable was added to indicate nonconvergence
c    for other programs.  It has no use in the atom program
c    and can be removed by the user.
c    All output from the routine is compatible to
c    the Berkeley/Sverre Froyen version.
c  njtj  ***  modifications  ***
c
c  twb   ***  modifications  ***
c    Underflow trap fixed
c  twb   ***  modifications  ***
c
c  njtj
c  &&&  Machine dependent Parameter
c  &&&    The value of expzer is machine dependent.
c  &&&    The user must switch in the correct value for
c  &&&    the machine in use from the list, or find
c  &&&    it for their machine.
c  &&&  Machine dependent Parameter
c  njtj
c
c  Integration coefficients
c
c
c  njtj  *** start modification  ***
c    Arrays added to gain speed.
c
c
c  njtj  ***  end modification  ***
c
C     .. Parameters ..
      double precision zero, pnine, two, etol
      parameter (zero=0.D0,pnine=0.9D0,two=2.D0,etol=-1.D-7)
      double precision tol
      parameter (tol=1.D-10)
      double precision abc1, abc2, abc3, abc4, abc5, amc0, amc1, amc2,
     &                 amc3, amc4
      parameter (abc1=190.1D0/72,abc2=-138.7D0/36,abc3=10.9D0/3,
     &          abc4=-63.7D0/36,abc5=25.1D0/72,amc0=25.1D0/72,
     &          amc1=32.3D0/36,amc2=-1.1D0/3,amc3=5.3D0/36,
     &          amc4=-1.9D0/72)
C     ..
C     .. Scalar Arguments ..
      double precision eigv, spin
      integer iflag, iorb, iter, n, l
C     ..
C     .. Array Arguments ..
      double precision ar(*), br(*), v(*)
C     ..
C     .. Local Scalars ..
      double precision aa, alf, arc, arctp, arp, bb, brc, brctp, brp,
     &                 dev, emax, emin, evold, expzer, factor, fb0, fb1,
     &                 temp, var0, vev, vzero, zeff
      integer icount, istop, itmax, j, j1, j2, j3, j4, j5, juflow, ll,
     &        lp, nctp, ninf, ninf1, ninf2, ninf3, ninf4, nodes
C     ..
C     .. Local Arrays ..
      double precision rabrlo(5), rlp(5)
C     ..
C     .. External Subroutines ..
      external ext
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, exp, log, sign, sqrt
C     ..
C     .. Arrays in Common ..
      double precision fa(nrmax), fb(nrmax), rab2(nrmax)
C     ..
C     .. Common blocks ..
      common  rab2, fa, fb
C     ..
c------Machine dependent parameter-
c------Require exp(-2*expzer) to be within the range of the machine
c
c      AG: Let TINY be the smallest normalized number:
c
c      Then: expzer = log(1.d0/sqrt(TINY))
c
c      We can obtain TINY by a call to smach (LINPACK/SCILIB):
c
c      TINY = smach(2)
c
c      Beware of machines with IEEE arithmetic (Sun, etc). 
c
      expzer = 3.7D2
cApollo      expzer = 3.7D2
cSun      expzer = 3.7D2
cVax      expzer = 44.D0
Cray      expzer =  2.8E3
c
c  njtj  *** major modification start  ***
c
c    Loop data calculated outside loop to gain speed.
c
      itmax = 100
      iflag = 0
      lp = l + 1
      ar(1) = zero
      if (l .eq. 0) then
         br(1) = b*a
      else
         br(1) = zero
      end if
      do 10 j = 2, nr
         ar(j) = zero
         br(j) = zero
   10 continue
c
c     Startup for predictor-corrector (rR goes as r^l near r=0)
c
      do 30 j = 2, 5
         rlp(j) = r(j)**lp
         rabrlo(j) = rab(j)*r(j)**l
   30 continue
c
      do 50 j = 1, nr
         rab2(j) = rab(j)*rab(j)
   50 continue
c
c   set underflow trap
c   twb *** begin modification ***
c
      juflow = 1
      do 60 j = 2, nr
         if (lp*abs(log(r(j))) .lt. sqrt(expzer)) go to 70
         juflow = j
   60 continue
   70 continue
c
c  twb *** end modification ***
c  njtj  *** end major modification  ***
c
c   determine effective charge and vzero for startup of
c   outward integration
c
c   ar = r**(l+1) * (1 + aa r + bb r**2 + ... )
c   aa = -znuc / lp     bb = (-2 znuc aa + v(0) - e)/(4 l + 6)
c
      zeff = zero
      if (spin .lt. 0.1D0 .and. viod(lp,2) .lt. -0.1D0) zeff = znuc
      if (spin .gt. 0.1D0 .and. viou(lp,2) .lt. -0.1D0) zeff = znuc
      aa = -zeff/lp
      vzero = -2*zeff*aa
      if (zeff .eq. zero) then
         if (spin .lt. 0.1D0) then
            vzero = vzero + viod(lp,2)/r(2)
         else
            vzero = vzero + viou(lp,2)/r(2)
         end if
      end if
      if (spin .lt. 0.1D0) then
         vzero = vzero + vid(2)
      else
         vzero = vzero + viu(2)
      end if
      var0 = zero
      if (l .eq. 0) var0 = -2*zeff
      if (l .eq. 1) var0 = two
c
      emax = zero
      emin = -two*100000
      if (eigv .gt. emax) eigv = emax
   80 continue
      if (itmax .lt. 2) write(6,9000) iorb, iter, eigv, nodes
 9000 format(' iorb =',i3,' iter =',i3,' ev =',d18.10,' nodes =',i2)
      if (itmax .eq. 0) then
         iflag = 1
c
         return
c
      end if
      if (eigv .gt. zero) then
         write(6,9010) iorb
         call ext(620+iorb)
      end if
 9010 format(//' error in difnrl - ev(',i2,') greater then v(infinty)')
c
c   find practical infinity ninf and classical turning
c   point nctp for orbital
c
      icount = 0
   90 continue
      icount = icount + 1
      do 100 j = nr, 2, -1
         temp = v(j) - eigv
         if (temp .lt. zero) temp = zero
         if (r(j)*sqrt(temp) .lt. expzer) go to 110
  100 continue
  110 continue
      ninf = j
      nctp = ninf - 5
      do 120 j = 2, ninf - 5
         if (v(j) .lt. eigv) nctp = j
  120 continue
      if (eigv .ge. etol*10) nctp = ninf - 5
      if (eigv .ge. etol) eigv = zero
      if (nctp .le. 6) then
         eigv = pnine*eigv
         if (icount .gt. 100) then
            write(6,9020) iorb
            call ext(650+iorb)
         end if
c
         go to 90
c
      end if
 9020 format(//'error in difnrl - cannot find the classical ',
     &      /' turning point for orbital ',i2)
c
c   outward integration from 1 to nctp
c   startup
c
      bb = (vzero-eigv)/(4*lp+2)
      do 130 j = 2, 5
         ar(j) = rlp(j)*(1+(aa+bb*r(j))*r(j))
         br(j) = rabrlo(j)*(lp+(aa*(lp+1)+bb*(lp+2)*r(j))*r(j))
  130 continue
c
c  njtj  ***  start major modification  ***
c    Predictor-corrector array added.
c
      fa(1) = br(1)
      fb(1) = b*br(1) + rab2(1)*var0
      fa(2) = br(2)
      fb(2) = b*br(2) + rab2(2)*(v(2)-eigv)*ar(2)
      fa(3) = br(3)
      fb(3) = b*br(3) + rab2(3)*(v(3)-eigv)*ar(3)
      fa(4) = br(4)
      fb(4) = b*br(4) + rab2(4)*(v(4)-eigv)*ar(4)
      fa(5) = br(5)
      fb(5) = b*br(5) + rab2(5)*(v(5)-eigv)*ar(5)
c
c   integration loop
c
      nodes = 0
      do 140 j = 6, nctp
c
c   predictor (Adams-Bashforth)
c
         j1 = j - 1
         j2 = j - 2
         j3 = j - 3
         j4 = j - 4
         j5 = j - 5
         vev = v(j) - eigv
         arp = ar(j1) + abc1*fa(j1) + abc2*fa(j2) + abc3*fa(j3) +
     &         abc4*fa(j4) + abc5*fa(j5)
         brp = br(j1) + abc1*fb(j1) + abc2*fb(j2) + abc3*fb(j3) +
     &         abc4*fb(j4) + abc5*fb(j5)
         fb1 = b*brp + rab2(j)*vev*arp
c
c   corrector (Adams-Moulton)
c
         arc = ar(j1) + amc0*brp + amc1*fa(j1) + amc2*fa(j2) +
     &         amc3*fa(j3) + amc4*fa(j4)
         brc = br(j1) + amc0*fb1 + amc1*fb(j1) + amc2*fb(j2) +
     &         amc3*fb(j3) + amc4*fb(j4)
         fb0 = b*brc + rab2(j)*vev*arc
c
c   error reduction step
c
         ar(j) = arc + amc0*(brc-brp)
         br(j) = brc + amc0*(fb0-fb1)
         fa(j) = br(j)
         fb(j) = b*br(j) + rab2(j)*vev*ar(j)
c
c   count nodes - if no underflow
c
         if (j .gt. juflow .and. ar(j)*ar(j-1) .lt.
     &       zero) nodes = nodes + 1
  140 continue
c
c  njtj  ***  end major modification  ***
c
      arctp = ar(nctp)
      brctp = br(nctp)
c
c   end outward integration
c
c   if number of nodes correct, start inward integration
c   else modify energy stepwise and try again
c
      if (nodes .ne. n-l-1) then
         if (nodes .lt. n-l-1) then
c
c  too few nodes; increase ev
c
            if (eigv .gt. emin) emin = eigv
            eigv = eigv - eigv/10
         else
c
c  too many nodes; decrease ev
c
            if (eigv .lt. emax) emax = eigv
            eigv = eigv + eigv/10
         end if
         itmax = itmax - 1
c
         go to 80
c
      end if
c
c   inward integration from ninf to nctp
c   startup
c
      do 150 j = ninf, ninf - 4, -1
         alf = v(j) - eigv
         if (alf .lt. zero) alf = zero
         alf = sqrt(alf)
         ar(j) = exp(-alf*r(j))
         br(j) = -rab(j)*alf*ar(j)
  150 continue
c
c  njtj  ***  start major modification  ***
c    Array for predictor-corrector added.
c
      fa(ninf) = br(ninf)
      fb(ninf) = b*br(ninf) + rab2(ninf)*(v(ninf)-eigv)*ar(ninf)
      ninf1 = ninf - 1
      fa(ninf1) = br(ninf1)
      fb(ninf1) = b*br(ninf1) + rab2(ninf1)*(v(ninf1)-eigv)*
     &            ar(ninf1)
      ninf2 = ninf - 2
      fa(ninf2) = br(ninf2)
      fb(ninf2) = b*br(ninf2) + rab2(ninf2)*(v(ninf2)-eigv)*
     &            ar(ninf2)
      ninf3 = ninf - 3
      fa(ninf3) = br(ninf3)
      fb(ninf3) = b*br(ninf3) + rab2(ninf3)*(v(ninf3)-eigv)*
     &            ar(ninf3)
      ninf4 = ninf - 4
      fa(ninf4) = br(ninf4)
      fb(ninf4) = b*br(ninf4) + rab2(ninf4)*(v(ninf4)-eigv)*
     &            ar(ninf4)
c
c   integration loop
c
      istop = ninf - nctp
      if (istop .lt. 5) go to 170
      do 160 j = ninf - 5, nctp, -1
c
c   predictor (Adams-Bashforth)
c
         j1 = j + 1
         j2 = j + 2
         j3 = j + 3
         j4 = j + 4
         j5 = j + 5
         vev = v(j) - eigv
         arp = ar(j1) - (abc1*fa(j1)+abc2*fa(j2)+abc3*fa(j3)+
     &         abc4*fa(j4)+abc5*fa(j5))
         brp = br(j1) - (abc1*fb(j1)+abc2*fb(j2)+abc3*fb(j3)+
     &         abc4*fb(j4)+abc5*fb(j5))
         fb0 = b*brp + rab2(j)*vev*arp
c
c   corrector (Adams-Moulton)
c
         arc = ar(j1) - (amc0*brp+amc1*fa(j1)+amc2*fa(j2)+amc3*fa(j3)+
     &         amc4*fa(j4))
         brc = br(j1) - (amc0*fb0+amc1*fb(j1)+amc2*fb(j2)+amc3*fb(j3)+
     &         amc4*fb(j4))
c
         fb1 = b*brc + rab2(j)*vev*arc
c
c   error reduction step
c
         ar(j) = arc - amc0*(brc-brp)
         br(j) = brc - amc0*(fb1-fb0)
         fa(j) = br(j)
         fb(j) = b*br(j) + rab2(j)*vev*ar(j)
  160 continue
c
c   end inward integration
c
c  njtj  *** end major modification  ***
c
c   rescale ar and br outside nctp to match ar(nctp) from
c   outward integration
c
  170 continue
      factor = arctp/ar(nctp)
      do 180 j = nctp, ninf
         ar(j) = factor*ar(j)
         br(j) = factor*br(j)
  180 continue
c
c   find normalizing factor
c
      factor = zero
      ll = 4
      do 190 j = 2, ninf
         factor = factor + ll*ar(j)*ar(j)*rab(j)
         ll = 6 - ll
  190 continue
      factor = factor/3
c
c   modify eigenvalue ev
c
      dev = arctp*(brctp-br(nctp))/(factor*rab(nctp))
      if (5*abs(dev) .gt. -eigv) dev = sign(eigv,dev)/5
      itmax = itmax - 1
      evold = eigv
      eigv = eigv + dev
      if (eigv .gt. emax) eigv = (evold+emax)/2
      if (eigv .lt. emin) eigv = (evold+emin)/2
      if (abs(dev) .gt. tol*(1-eigv)) go to 80
c
c   normalize wavefunction and change br from d(ar)/dj to d(ar)/dr
c
      factor = 1/sqrt(factor)
      do 200 j = 1, ninf
         ar(j) = factor*ar(j)
         br(j) = factor*br(j)/rab(j)
  200 continue
c
      return
c
      end
C
c $Id: difrel.f,v 1.3 1997/05/22 17:32:07 wdpgaara Exp $
c
c $Log: difrel.f,v $
c Revision 1.3  1997/05/22 17:32:07  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.2  1994/02/18  01:26:18  garcia
c *** empty log message ***
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine difrel(iter,iorb,v,ar,br,n,l,spin,eigv)
c
      implicit none
c
      include 'radial.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
c     
c  difrel integrates the relativistic Dirac equation
c  it finds the eigenvalue ev, the major and minor component
c  of the wavefunction, ar and br.  It uses an intial guess
c  for the eigenvalues from dsolv1
c
c  njtj  ***  modifications  ***
c    This routine has major modifications.
c    1)The data needed inside the loops has been calculated
c    outside the main loop(increases speed for non-opt
c    compilers, i.e. dumb compilers).
c    2)The predict/correct values are placed in an array.
c    Output is unchanged
c  njtj  ***  modifications  ***
c
c  twb   ***  modifications  ***
c    Underflow trap fixed
c  twb   ***  modifications  ***
c  njtj
c  &&&  Machine dependent Parameter
c  &&&    The value of expzer is machine dependent.
c  &&&    The user must switch in the correct value for
c  &&&    the machine in use from the list, or find
c  &&&    it for their machine.
c  &&&  Machine dependent Parameter
c  njtj
c
c  Tolerance
c
c
c  Integration coefficients
c
c------Machine dependent parameter-
c------Require exp(-2*expzer) to be within the range of the machine
c
C     .. Parameters ..
      double precision zero, pnine, one, ai
      parameter (zero=0.D0,pnine=0.9D0,one=1.D0,ai=2*137.0360411D0)
      double precision etol
      parameter (etol=-1.D-7)
      double precision tol
      parameter (tol=1.D-10)
      double precision abc1, abc2, abc3, abc4, abc5, amc0, amc1, amc2,
     &                 amc3, amc4
      parameter (abc1=190.1D0/72,abc2=-138.7D0/36,abc3=10.9D0/3,
     &          abc4=-63.7D0/36,abc5=25.1D0/72,amc0=25.1D0/72,
     &          amc1=32.3D0/36,amc2=-1.1D0/3,amc3=5.3D0/36,
     &          amc4=-1.9D0/72)
C     ..
C     .. Scalar Arguments ..
      double precision eigv, spin
      integer iorb, iter, n, l
C     ..
C     .. Array Arguments ..
      double precision ar(*), br(*), v(*)
C     ..
C     .. Local Scalars ..
      double precision a1, a2, ai2, alf, arc, arin, arout, arp, arpin,
     &                 arpout, az, b0, b1, b2, brc, brp, dev, emax,
     &                 emin, evold, evv, evvai2, expzer, factor, faj,
     &                 fbj, s, temp, vzero
      integer icount, istop, itmax, j, juflow, ka, ll, nctp, ninf, nodes
C     ..
C     .. Local Arrays ..
      double precision rs(5)
C     ..
C     .. External Subroutines ..
      external ext
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, exp, log, sign, sqrt
C     ..
C     .. Arrays in Common ..
      double precision fa(nrmax), fb(nrmax), rabai(nrmax), rabkar(nrmax)
C     ..
C     .. Common blocks ..
      common  rabkar, rabai, fa, fb
C     ..
      expzer = 3.7D2
cApollo      expzer = 3.7D2
cSun      expzer = 3.7D2
cVax      expzer = 44.D0
Cray      expzer = 2.8E3
c
      itmax = 100
      ai2 = ai*ai
      az = znuc/(2*ai)
      ka = l + 1
      if (spin .lt. 0.1D0 .and. l .ne. 0) ka = -l
c
c  determine effective charge and vzero for startup of
c  outward integration
c  ar = r**s * (1  + a1 r + a2 r**2 + ... )
c  br = r**s * (b0 + b1 r + b2 r**2 + ... )
c  s = sqrt (ka**2 - az**2)    b0 = - az / (s + ka)
c  an = (az (v0 - e) a(n-1) - (s + n + ka) (v0 - e - ai**2) b(n-1))
c        / (n ai (2 s + n))
c  bn = ((v0 - e) a(n-1) - 2 znuc an ) / ( ai (s + n + ka))
c
      s = sqrt(ka*ka-az*az)
      if (ka .gt. 0) then
         b0 = -az/(s+ka)
      else
         b0 = (s-ka)/az
      end if
      if (spin .lt. 0.1D0) then
         vzero = vid(2)
      else
         vzero = viu(2)
      end if
c
c  njtj  ***  start major modification  ***
c    Loop data calculated only once.
c    Set ar() and br() to zero.
c
      do 10 j = 1, nr
         ar(j) = zero
         br(j) = zero
   10 continue
      do 20 j = 2, nr
         rabkar(j) = rab(j)*ka/r(j)
   20 continue
      do 30 j = 2, nr
         rabai(j) = rab(j)/ai
   30 continue
      do 40 j = 2, 5
         rs(j) = r(j)**s
   40 continue
c
c  set the underflow trap.
c  twb *** begin modification ***
c
      juflow = 1
      do 50 j = 2, nr
         if (s*abs(log(r(j))) .lt. sqrt(expzer)) go to 60
         juflow = j
   50 continue
   60 continue
c
c  twb *** end modification ***
c  njtj *** end major modification  ***
c
      emax = zero
      emin = -one*100000
      if (eigv .gt. emax) eigv = emax
   70 continue
      if (itmax .lt. 2) write(6,9000) iorb, iter, eigv, nodes
 9000 format(' iorb =',i3,' iter =',i3,' ev =',d18.10,' nodes =',i2)
      if (itmax .eq. 0) return
      if (eigv .gt. zero) then
         write(6,9010) iorb
         call ext(620+iorb)
      end if
 9010 format(//' error in difrel - ev(',i2,') greater then v(infinty)')
c
c  Find practical infinity ninf and classical turning
c  point nctp for orbital.
c
      icount = 0
   80 continue
      icount = icount + 1
      do 90 j = nr, 2, -1
         temp = v(j) - eigv
         if (temp .lt. zero) temp = zero
         if (r(j)*sqrt(temp) .lt. expzer) go to 100
   90 continue
  100 continue
      ninf = j
      nctp = ninf - 5
      do 110 j = 2, ninf - 5
         if (v(j) .lt. eigv) nctp = j
  110 continue
      if (eigv .ge. etol*100) nctp = ninf - 5
      if (eigv .ge. etol) eigv = zero
      if (nctp .le. 6) then
         eigv = pnine*eigv
         if (icount .gt. 100) then
            write(6,9020) iorb
            call ext(650+iorb)
         end if
c
         go to 80
c
      end if
 9020 format(//'error in difrel - cannot find classical',
     &      /'turning point in orbital ',i2)
c
c  Outward integration from 1 to nctp, startup.
c
      a1 = (az*(vzero-eigv)-(s+1+ka)*(vzero-eigv-ai2)*b0)/
     &     (ai*(2*s+1))
      b1 = ((vzero-eigv)-2*znuc*a1)/(ai*(s+1+ka))
      a2 = (az*(vzero-eigv)*a1-(s+2+ka)*(vzero-eigv-ai2)*b1)/
     &     (2*ai*(2*s+2))
      b2 = ((vzero-eigv)*a1-2*znuc*a2)/(ai*(s+2+ka))
      do 120 j = 2, 5
         ar(j) = rs(j)*(1+(a1+a2*r(j))*r(j))
         br(j) = rs(j)*(b0+(b1+b2*r(j))*r(j))
  120 continue
      fa(1) = zero
      fb(1) = zero
      fa(2) = rabkar(2)*ar(2) + (eigv-v(2)+ai2)*br(2)*rabai(2)
      fb(2) = -rabkar(2)*br(2) - (eigv-v(2))*ar(2)*rabai(2)
      fa(3) = rabkar(3)*ar(3) + (eigv-v(3)+ai2)*br(3)*rabai(3)
      fb(3) = -rabkar(3)*br(3) - (eigv-v(3))*ar(3)*rabai(3)
      fa(4) = rabkar(4)*ar(4) + (eigv-v(4)+ai2)*br(4)*rabai(4)
      fb(4) = -rabkar(4)*br(4) - (eigv-v(4))*ar(4)*rabai(4)
      fa(5) = rabkar(5)*ar(5) + (eigv-v(5)+ai2)*br(5)*rabai(5)
      fb(5) = -rabkar(5)*br(5) - (eigv-v(5))*ar(5)*rabai(5)
c
c  Intergration loop.
c
      nodes = 0
      do 130 j = 6, nctp
c
c  Predictor (Adams-Bashforth).
c
         evvai2 = eigv - v(j) + ai2
         evv = eigv - v(j)
         arp = ar(j-1) + abc1*fa(j-1) + abc2*fa(j-2) + abc3*fa(j-3) +
     &         abc4*fa(j-4) + abc5*fa(j-5)
         brp = br(j-1) + abc1*fb(j-1) + abc2*fb(j-2) + abc3*fb(j-3) +
     &         abc4*fb(j-4) + abc5*fb(j-5)
         fa(j) = rabkar(j)*arp + evvai2*brp*rabai(j)
         fb(j) = -rabkar(j)*brp - evv*arp*rabai(j)
c
c  Corrector (Adams-Moulton).
c
         arc = ar(j-1) + amc0*fa(j) + amc1*fa(j-1) + amc2*fa(j-2) +
     &         amc3*fa(j-3) + amc4*fa(j-4)
         brc = br(j-1) + amc0*fb(j) + amc1*fb(j-1) + amc2*fb(j-2) +
     &         amc3*fb(j-3) + amc4*fb(j-4)
         faj = rabkar(j)*arc + evvai2*brc*rabai(j)
         fbj = -rabkar(j)*brc - evv*arc*rabai(j)
c
c  Error reduction step.
c
         ar(j) = arc + amc0*(faj-fa(j))
         br(j) = brc + amc0*(fbj-fb(j))
         fa(j) = rabkar(j)*ar(j) + evvai2*br(j)*rabai(j)
         fb(j) = -rabkar(j)*br(j) - evv*ar(j)*rabai(j)
c
c  Count nodes - if no underflow.
c
         if (j .gt. juflow .and. ar(j)*ar(j-1) .lt.
     &       zero) nodes = nodes + 1
  130 continue
      arout = ar(nctp)
      arpout = fa(nctp)
c
c  End outward integration.
c  If number of nodes correct, start inward integration
c  else modify energy stepwise and try again.
c
      if (nodes .ne. n-l-1) then
c
c  too many nodes decrease ev
c
         if (nodes .gt. n-l-1) then
            if (eigv .lt. emax) emax = eigv
            eigv = eigv + eigv/10
c
c  too few nodes increase ev
c
         else
            if (eigv .gt. emin) emin = eigv
            eigv = eigv - eigv/10
         end if
         itmax = itmax - 1
c
         go to 70
c
      end if
c
c  Inward integration from ninf to nctp startup.
c
      do 140 j = ninf, ninf - 4, -1
         alf = v(j) - eigv
         if (alf .lt. zero) alf = zero
         alf = sqrt(alf)
         ar(j) = exp(-alf*r(j))
         br(j) = ai*(alf+ka/r(j))*ar(j)/(v(j)-eigv-ai2)
  140 continue
      fa(ninf) = rabkar(ninf)*ar(ninf) +
     &           (eigv-v(ninf)+ai2)*br(ninf)*rabai(ninf)
      fb(ninf) = -rabkar(ninf)*br(ninf) -
     &           (eigv-v(ninf))*ar(ninf)*rabai(ninf)
      fa(ninf-1) = rabkar(ninf-1)*ar(ninf-1) +
     &             (eigv-v(ninf-1)+ai2)*br(ninf-1)*rabai(ninf-1)
      fb(ninf-1) = -rabkar(ninf-1)*br(ninf-1) -
     &             (eigv-v(ninf-1))*ar(ninf-1)*rabai(ninf-1)
      fa(ninf-2) = rabkar(ninf-2)*ar(ninf-2) +
     &             (eigv-v(ninf-2)+ai2)*br(ninf-2)*rabai(ninf-2)
      fb(ninf-2) = -rabkar(ninf-2)*br(ninf-2) -
     &             (eigv-v(ninf-2))*ar(ninf-2)*rabai(ninf-2)
      fa(ninf-3) = rabkar(ninf-3)*ar(ninf-3) +
     &             (eigv-v(ninf-3)+ai2)*br(ninf-3)*rabai(ninf-3)
      fb(ninf-3) = -rabkar(ninf-3)*br(ninf-3) -
     &             (eigv-v(ninf-3))*ar(ninf-3)*rabai(ninf-3)
      fa(ninf-4) = rabkar(ninf-4)*ar(ninf-4) +
     &             (eigv-v(ninf-4)+ai2)*br(ninf-4)*rabai(ninf-4)
      fb(ninf-4) = -rabkar(ninf-4)*br(ninf-4) -
     &             (eigv-v(ninf-4))*ar(ninf-4)*rabai(ninf-4)
c
c  Integration loop.
c
      istop = ninf - nctp
      if (istop .lt. 5) go to 160
      do 150 j = ninf - 5, nctp, -1
c
c  Predictor (Adams-Bashforth).
c
         evvai2 = eigv - v(j) + ai2
         evv = eigv - v(j)
         arp = ar(j+1) - (abc1*fa(j+1)+abc2*fa(j+2)+abc3*fa(j+3)+
     &         abc4*fa(j+4)+abc5*fa(j+5))
         brp = br(j+1) - (abc1*fb(j+1)+abc2*fb(j+2)+abc3*fb(j+3)+
     &         abc4*fb(j+4)+abc5*fb(j+5))
         fa(j) = rabkar(j)*arp + evvai2*brp*rabai(j)
         fb(j) = -rabkar(j)*brp - evv*arp*rabai(j)
c
c  Corrector (Adams-Moulton).
c
         arc = ar(j+1) - (amc0*fa(j)+amc1*fa(j+1)+amc2*fa(j+2)+
     &         amc3*fa(j+3)+amc4*fa(j+4))
         brc = br(j+1) - (amc0*fb(j)+amc1*fb(j+1)+amc2*fb(j+2)+
     &         amc3*fb(j+3)+amc4*fb(j+4))
         faj = rabkar(j)*arc + evvai2*brc*rabai(j)
         fbj = -rabkar(j)*brc - evv*arc*rabai(j)
c
c  Error reduction step.
c
         ar(j) = arc + amc0*(faj-fa(j))
         br(j) = brc + amc0*(fbj-fb(j))
         fa(j) = rabkar(j)*ar(j) + evvai2*br(j)*rabai(j)
         fb(j) = -rabkar(j)*br(j) - evv*ar(j)*rabai(j)
  150 continue
  160 continue
      arin = ar(nctp)
      arpin = fa(nctp)
c
c  End inward integration
c  Rescale ar and br outside nctp to match ar(nctp) from
c  outward integration.
c
      factor = arout/arin
      do 170 j = nctp, ninf
         ar(j) = factor*ar(j)
         br(j) = factor*br(j)
  170 continue
      arpin = factor*arpin
c
c  Find the normalizing factor.
c
      factor = zero
      ll = 4
      do 180 j = 2, ninf
         factor = factor + ll*(ar(j)*ar(j)+br(j)*br(j))*rab(j)
         ll = 6 - ll
  180 continue
      factor = factor/3
c
c  Modify the eigenvalue ev.
c
      dev = arout*(arpout-arpin)/(factor*rab(nctp))
      if (5*abs(dev) .gt. -eigv) dev = sign(eigv,dev)/5
      itmax = itmax - 1
      evold = eigv
      eigv = eigv + dev
      if (eigv .gt. emax) then
         eigv = (evold+emax)/2
      else if (eigv .lt. emin) then
         eigv = (evold+emin)/2
      end if
      if (abs(dev) .gt. tol*(1-eigv)) go to 70
c
c  Normalize the wavefunction.
c
      factor = 1/sqrt(factor)
      do 190 j = 1, ninf
         ar(j) = factor*ar(j)
         br(j) = factor*br(j)
  190 continue
c
      return
c
      end
C
c $Id: dmixp.f,v 1.4 1999/02/25 13:48:42 wdpgaara Exp $
c
      subroutine dmixp(a,b,beta,icy,id,nmsh,c,d,vn1,vn12,vn2,vn22)
c
      implicit none
c
C*    ADAPTED FROM K.C.PANDEY
C*    USING ANDERSON'S EXTRAPOLATION SCHEME
C*    EQS 4.1-4.9,4.15-4.18 OF
C*    D.G.ANDERSON J.ASSOC.COMPUTING MACHINERY,12,547(1965)
c
C*    COMPUTES A NEW VECTOR IN A ITERATIVE SCHEME
c
C*    INPUT A=NEWPOT B=OLDPOT
C*    OUTPUT A=A-B B=NEWPOT
C*    BETA=MIXING,IN=ITER. NUMBER
C*    ID=1,2 OR 3 DIFF CONV METH.
C*    ICY CYCLE NUMBER ,ICY=1 ON FIRST/ZEROTH CALL
C*    C,D WORK ARRAYS OF SIZE NMSH
C*    VN1,VN12,VN2,VN22 STORAGE ARRAYS OF SIZE NMSH
C
c
c     Modified by Alberto Garcia to make use of Level 1 BLAS
c
C     .. Parameters ..
      double precision uze, um, detol
      parameter (uze=0.0D0,um=1.0D0,detol=1.D-9)
C     ..
C     .. Scalar Arguments ..
      double precision beta
      integer icy, id, nmsh
C     ..
C     .. Array Arguments ..
      double precision a(*), b(*), c(*), d(*)
C     ..
C     .. Local Scalars ..
      double precision a2, bt1, bt2, d11, d12, d22, det, dett, r2, rd1m,
     &                 rd2m, t1, t2, x
      integer in
C     ..
C     .. External Subroutines ..
c     BLAS level 1 (SCILIB)
      double precision sdot
      external sdot, saxpy, sscal, scopy
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs
C     ..
      double precision vn1(*), vn12(*), vn2(*), vn22(*)
C     ..
      in = icy - 1
      if (in .eq. 0) then
         call saxpy(nmsh,um,a,1,b,1)
c
         return
c
      end if
      call saxpy(nmsh,-um,b,1,a,1)
      r2 = sdot(nmsh,a,1,a,1)
      if (id .eq. 1) then
         call saxpy(nmsh,beta,a,1,b,1)
c
         return
c
      end if
      if (in .eq. 1) then
         call scopy(nmsh,a,1,vn1,1)
         call scopy(nmsh,b,1,vn2,1)
         call saxpy(nmsh,beta,a,1,b,1)
c
         return
c
      end if
      call scopy(nmsh,vn1,1,c,1)
      if (id .eq. 3 .and. in .gt. 2) then
         call scopy(nmsh,vn12,1,d,1)
      end if
      call scopy(nmsh,a,1,vn1,1)
      if (id .gt. 2 .and. in .gt. 1) then
         call scopy(nmsh,c,1,vn12,1)
      end if
      call saxpy(nmsh,-um,a,1,c,1)
      d11 = sdot(nmsh,c,1,c,1)
      rd1m = sdot(nmsh,a,1,c,1)
      if (in .le. 2 .or. id .le. 2) then
         t1 = -rd1m/d11
         x = um - t1
         bt1 = beta*t1
         call sscal(nmsh,beta,a,1)
         call saxpy(nmsh,bt1,c,1,a,1)
         call scopy(nmsh,vn2,1,d,1)
         call saxpy(nmsh,t1,d,1,a,1)
         call scopy(nmsh,b,1,vn2,1)
         if (id .gt. 2 .and. in .eq. 2) then
            call scopy(nmsh,d,1,vn22,1)
         end if
         call sscal(nmsh,x,b,1)
         call saxpy(nmsh,1.d0,a,1,b,1)
c
         return
c
      end if
      call saxpy(nmsh,-um,a,1,d,1)
      d22 = sdot(nmsh,d,1,d,1)
      d12 = sdot(nmsh,c,1,d,1)
      rd2m = sdot(nmsh,a,1,d,1)
      a2 = d11*d22
      det = a2 - d12*d12
      dett = det/a2
      if (abs(dett) .ge. detol) then
         t1 = (-rd1m*d22+rd2m*d12)/det
         t2 = (rd1m*d12-rd2m*d11)/det
      else
         t1 = -rd1m/d11
         t2 = uze
      end if
      x = um - t1 - t2
      bt1 = beta*t1
      bt2 = beta*t2
      call sscal(nmsh,beta,a,1)
      call saxpy(nmsh,bt1,c,1,a,1)
      call saxpy(nmsh,bt2,d,1,a,1)
      call saxpy(nmsh,t1,vn2,1,a,1)
      call saxpy(nmsh,t2,vn22,1,a,1)
      call scopy(nmsh,vn2,1,vn22,1)
      call scopy(nmsh,b,1,vn2,1)
      call sscal(nmsh,x,b,1)
      call saxpy(nmsh,1.d0,a,1,b,1)
c
      return
c
      end







C
c $Id: dsolv1.f,v 1.3 1999/02/25 13:49:13 wdpgaara Exp $
c
      subroutine dsolv1(nfirst,nlast,nn)
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'ion.h'
      include 'charge.h'
      include 'elecpot.h'
      include 'energy.h'
c
c   dsolv1 finds the (non)-relativistic wave function
c   using finite differences and matrix diagonalization.
c   An initial guess for the eigenvalues need not be supplied.
c
C     .. Parameters ..
      double precision zero, one, pone, opf
      parameter (zero=0.D0,one=1.D0,pone=0.1D0,opf=1.5D0)
C     ..
C     .. Scalar Arguments ..
      integer nfirst, nlast
C     ..
c     nn will generally be 'no'
c
      integer nn(*)
C     ..
C     .. Local Scalars ..
      double precision bl, bu, c1, c2, denr, eps
      integer i, ierr, j, k, ki, kn, l, llp, nrm
C     ..
C     .. Local Arrays ..
      double precision e(10)
      integer ind(10), nmax(2,5)
C     ..
C     .. External Subroutines ..
      external tinvit, tridib
C     ..
      double precision d(nrmax), dk(nrmax), rv1(nrmax), rv2(nrmax),
     &                 rv3(nrmax), rv4(nrmax), rv5(nrmax), sd(nrmax),
     &                 sd2(nrmax), z(6*nrmax)
C     ..
c
c   Initialize the charge density arrays.
c
      do 10 i = 1, nr
         cdd(i) = zero
         cdu(i) = zero
   10 continue
c
c   Find the max n given l and s.
c   Zero spin is treated as down.
c
      do 40 i = 1, 2
         do 30 j = 1, lmax
            nmax(i,j) = 0
            do 20 k = nfirst, nlast
               if (nn(k) .le. 0) go to 20
               if (lo(k) .ne. j-1) go to 20
               if ((so(k)-pone)*(i-opf) .lt. zero) go to 20
               nmax(i,j) = nn(k)
   20       continue
   30    continue
   40 continue
c
c   Set up hamiltonian matrix for kinetic energy.
c   Only the diagonal depends on the potential.
c
      c2 = -one/b**2
      c1 = -2*one*c2 + one/4
      dk(1) = c1/(r(2)+a)**2
      sd(1) = zero
      sd2(1) = zero
      do 50 i = 3, nr
         dk(i-1) = c1/(r(i)+a)**2
         sd(i-1) = c2/((r(i)+a)*(r(i-1)+a))
         sd2(i-1) = sd(i-1)**2
   50 continue
c
c   Start loop over spin down=1 and spin up=2.
c
      nrm = nr - 1
      do 100 i = 1, 2
c
c   Start loop over s p d... states.
c
         do 90 j = 1, lmax
            if (nmax(i,j) .eq. 0) go to 90
            llp = j*(j-1)
            do 60 k = 2, nr
               if (i .eq. 1) then
                  d(k-1) = dk(k-1) + (viod(j,k)+llp/r(k))/r(k) + vid(k)
               else
                  d(k-1) = dk(k-1) + (viou(j,k)+llp/r(k))/r(k) + viu(k)
               end if
   60       continue
c
c   Diagonalize the matrix.
c
            eps = -one
            call tridib(nrm,eps,d,sd,sd2,bl,bu,1,nmax(i,j),e,ind,ierr,
     &                  rv4,rv5)
            if (ierr .ne. 0) write(6,9000) ierr
 9000       format(/' error in tridib ****** ierr =',i3,/)
            call tinvit(nrm,nrm,d,sd,sd2,nmax(i,j),e,ind,z,ierr,rv1,rv2,
     &                  rv3,rv4,rv5)
            if (ierr .ne. 0) write(6,9010) ierr
 9010       format(/' error in tinvit ****** ierr =',i3,/)
c
c   Save the energy levels and add to charge density.
c
            ki = 1
            kn = 0
            do 80 k = nfirst, nlast
               if (nn(k) .le. 0) go to 80
               if (lo(k) .ne. j-1) go to 80
               if ((so(k)-pone)*(i-opf) .lt. zero) go to 80
               ev(k) = e(ki)
               do 70 l = 2, nr
                  denr = zo(k)*z(kn+l-1)**2/rab(l)
                  if (i .eq. 1) then
                     cdd(l) = cdd(l) + denr
                  else
                     cdu(l) = cdu(l) + denr
                  end if
   70          continue
               ki = ki + 1
               kn = kn + nrm
   80       continue
   90    continue
  100 continue
c
c   End loop over s p and d states.
c
      return
c
      end
c
c $Id: pseudo.f,v 1.4 1999/02/26 14:26:46 wdpgaara Exp $
c
      subroutine Pseudo(pot_id,headline,ps_generator)
c
c     Generates the pseudopotential. The particular flavor is
c     determined by ps_generator.
c
      implicit none
c
      character pot_id*40, headline*79
      external ps_generator
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'elecpot.h'
c
C     .. Parameters ..
      double precision zero, one
      parameter (zero=0.D0,one=1.D0)
C     ..
C     .. Local Scalars ..
      integer i, llp, lp
C     ..
C     .. Local Arrays ..
      double precision ar(nrmax), br(nrmax)
C     ..
C     .. External Subroutines ..
      external ext, wf, wrapup
C     ..
c
      do 10 i = 1, 5
         indd(i) = 0
         indu(i) = 0
   10 continue
c
      if (ncore .eq. norb) return
      if (job .ne. 1 .and. job .ne. 2 .and. job .ne. 3) return
c
c
c  read rc(s),rc(p),rc(d),rc(f),cfac,rcfac
c
c  cfac is used for the pseudocore - the pseudocore stops where
c  the core charge density equals cfac times the renormalized
c  valence charge density (renormalized to make the atom neutral).
c
c  If cfac is input as negative, the full core charge is used,
c  if cfac is input as zero, it is set equal to one.
c
c  rcfac is used for the pseudocore cut off radius.  If set
c  to less than or equal to zero cfac is used.  cfac must be
c  set to greater than zero.
c
      read(5,9000) (rc_input(i),i=1,4), cfac, rcfac
 9000 format(6f10.5)
      if (cfac .eq. zero) cfac = one
c
c   Reset vod and vou to zero.  They are here used to store
c   the pseudo valence charge density.
c
      do 20 i = 1, nr
         vod(i) = zero
         vou(i) = zero
   20 continue
c
c  Print the heading.
c
      write(6,9010) nameat, pot_id, headline
 9010 format(//1x,a2,' pseudopotential generation: ',a,//,a,/)
c
c      start loop over valence orbitals
c
      ncp = ncore + 1
c
      do 220 i = ncp, norb
c
         lp = lo(i) + 1
         llp = lo(i)*lp
         if (down(i)) then
            if (indd(lp) .ne. 0) then
               write(6,9020) 'down', lp - 1
               call ext(800+lp)
            else
               indd(lp) = i
            end if
         else
            if (indu(lp) .ne. 0) then
               write(6,9020) 'up',lp - 1
               call ext(810+lp)
            else
               indu(lp) = i
            end if
         end if
 9020    format(//' error in pseudo - two ',a4,
     &            ' spin orbitals of the same ',
     &         /'angular momentum (',i1,') exist')
c
c
c      Find all electron wave function and its nodes and
c      extrema.
c
         call Wf(i,ar,br)
c
c      Find the pseudopotential
c
         call ps_generator(i,ar,br)
c
  220 continue
c
c  End loop over valence orbitals.
c
      call wrapup(pot_id)
c
      return
c
      end
c
      subroutine wf(i,ar,br)
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
      include 'energy.h'
c
c     Solves the wave equation for orbital i and
c     finds the extrema and nodes.
c
      double precision zero, ai, pnine
      parameter (zero=0.d0,ai=2*137.0360411D0,pnine=0.9d0)
c
      integer i, nextr
      double precision rextr, rzero
c
      double precision ar(*), br(*)
c
      integer iflag, ist, j, lp, llp, ka
      double precision arp, arpm
      double precision v(nrmax)
c
      lp = lo(i) + 1
      llp = lo(i)*lp
c
      do 10 j = 1, nr
         ar(j) = 0.d0
   10 continue
      if (down(i)) then
         do 20 j = 2, nr
            v(j) = viod(lp,j)/r(j) + vid(j)
   20    continue
      else
         do 30 j = 2, nr
            v(j) = viou(lp,j)/r(j) + viu(j)
   30    continue
      end if
c
      if ( .not. relativistic) then
c
c           Add 'centrifugal term'
c
         do 40 j = 2, nr
            v(j) = v(j) + llp/r(j)**2
   40    continue
      end if
c
      if ( .not. relativistic) then
         call difnrl(0,i,v,ar,br,no(i),lo(i),so(i),ev(i),iflag)
      else
         call difrel(0,i,v,ar,br,no(i),lo(i),so(i),ev(i))
      end if
c
c     Plot and make the wavefunction 'upright'
c
      ist = nint(sign(1.d0,ar(nr-85)))
cGuima modifica
      call potrw(ar,r,nr,lo(i),1,ist,rc(lp))
cGuima fim da modificaccao
c
      do 50 j = 1, nr
         ar(j) = ar(j)*ist
         br(j) = br(j)*ist
   50 continue
c
c  Find the last zero and extremum.
c
      ka = lo(i) + 1
      if (down(i) .and. lo(i) .ne. 0) ka = -lo(i)
      nextr = no(i) - lo(i)
      rzero = zero
      arp = br(2)
c
      if (relativistic) then
         if (down(i)) then
            arp = ka*ar(2)/r(2) + (ev(i)-viod(lp,2)/r(2)-vid(2)+ai*ai)*
     &            br(2)/ai
         else
            arp = ka*ar(2)/r(2) + (ev(i)-viou(lp,2)/r(2)-viu(2)+ai*ai)*
     &            br(2)/ai
         end if
      end if
c
      do 60 j = 3, nr
c
         if (nextr .eq. 0) go to 70
c
         if (ar(j-1)*ar(j) .le. zero) rzero = (ar(j)*r(j-1)-
     &       ar(j-1)*r(j))/(ar(j)-ar(j-1))
         arpm = arp
         arp = br(j)
c
         if (relativistic) then
            if (down(i)) then
               arp = ka*ar(j)/r(j) + (ev(i)-viod(lp,j)/r(j)-vid(j)+
     &               ai*ai)*br(j)/ai
            else
               arp = ka*ar(j)/r(j) + (ev(i)-viou(lp,j)/r(j)-viu(j)+
     &               ai*ai)*br(j)/ai
            end if
         end if
c
         if (arp*arpm .le. zero) then
            rextr = (arp*r(j-1)-arpm*r(j))/(arp-arpm)
            nextr = nextr - 1
         end if
c
   60 continue
   70 continue
c
c  Check rc, if outside bounds reset.
c
      if (rzero .lt. r(2)) rzero = r(2)
c
c  Check rc if inside rzero,
c  reset to .9 between rmax and rzero if inside
c  if rc(lp) is negative, rc(lp) is percent of way
c  betweeen rzero and rmax.
c
      if (rc_input(lp) .gt. rzero) then
c
c           do nothing
c
      else if (rc_input(lp) .ge. zero) then
c
c        rc is inside the node...
c
         write(6,'(a,3f8.3)')
     $     ' Requested rc inside node ! ** rzero, rextr, rc_input:',
     $                          rzero, rextr, rc_input(lp)
         rc_input(lp) = rzero + pnine*(rextr-rzero)
         write(6,'(a,f6.3)') ' rc changed to ', rc_input(lp)
c
      else
c
c        compute rc as a fraction of rextr-rzero
c
         rc_input(lp) = rzero - rc_input(lp)*(rextr-rzero)
         write(6,'(a,f5.2)') 'rc set to ', rc_input(lp)
c
      end if
c
c        Warn the user if rc > rextr
c
      if (rc_input(lp) .ge. rextr) write(6,9000) rc_input(lp), rextr
 9000 format(' Core radius (',f5.2,
     &      ') outside wfn extremum (',f5.2,')')
c
      return
c
      end


c
c $Id: hsc.f,v 1.2 1997/05/22 17:32:15 wdpgaara Exp $
c
c $Log: hsc.f,v $
c Revision 1.2  1997/05/22 17:32:15  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine hsc(i,ar,br)
c
c     Implements the HSC method.
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
      include 'energy.h'
c
      double precision zero, one
      parameter (zero=0.D0,one=1.D0)
      double precision small, small2, small3, pzfive
      parameter (small=1.D-13,small2=1.D-10,small3=1.D-18,pzfive=.05D0)
      double precision pfive, small4
      parameter (pfive=0.5D0,small4=1.D-6)
c
      integer i
      double precision ar(*), br(*)
c
      double precision aa, ag, cl, dcl, delta, dev, devold, eviae, fjm1,
     &                 gamma, gg, gpp, rra, rrc, rrp
c
      integer iflag, ist, j, j3rc, k, ll, llp, lp
      character id*1
      integer nops(norbmx)
c
      double precision f(nrmax), g(nrmax), v(nrmax), arps(nrmax)
c
      external difnrl, dsolv2, potrw, wtrans
      intrinsic nint, sign
c
c     Do everything non-relativistically...
c
      if (polarized) then
         id = 's'
      else
         id = ' '
      endif
c
c  Reset the n quantum number to give the proper number of
c  nodes (0) for the pseudowavefunction. Zero out the rest.
c
      do 10 j = 1, norb
         nops(j) = 0
   10 continue
      lp = lo(i) + 1
      llp = lo(i) * lp
      rc(lp) = rc_input(lp)
      nops(i) = lp
c
c  njtj  ***  modification start  ***
c  Set up the functions f(r/rc) and g(r/rc) and
c  modify the ionic potential.
c
      aa = 4*one
      dcl = -6*one*lp
      cl = dcl
c
      do 20 j = 1, nr
         rrc = r(j)/rc(lp)
         rra = rrc**aa
         f(j) = zero
         if (rra .lt. 88*one) f(j) = exp(-rra)
         g(j) = rrc**lp*f(j)
         fjm1 = one - f(j)
         if (fjm1 .lt. small4) fjm1 = (one-pfive*rra)*rra
         if (down(i)) then
            viod(lp,j) = fjm1*viod(lp,j) - f(j)*r(j)*vid(j) +
     &                   dcl*r(j)*f(j)
         else
            viou(lp,j) = fjm1*viou(lp,j) - f(j)*r(j)*viu(j) +
     &                   dcl*r(j)*f(j)
         end if
         if (rrc .lt. 3*one) j3rc = j
   20 continue
      dcl = dcl/2
c
c   Start the iteration loop to find cl.
c
      eviae = ev(i)
      devold = zero
      do 60 j = 1, 100
c
         call dsolv2(j,2,id,1,norb,ncore,nops)
         dev = eviae - ev(i)
c
c    The abs(dev-devold) condition was added to eliminate
c    division by zero errors in the calculation of
c    dcl = -dev*dcl / (dev-devold).
c
         if (((abs(dev).lt.small2).or.(abs(dev-devold).lt.small3)) .and.
     &       (j.ne.1)) then
c
            go to 70
c
         else
            if (j .gt. 20 .or. abs(dev) .lt. 0.001D0) then
c
c                 Use newton-raphson iteration to change cl.
c
               dcl = -dev*dcl/(dev-devold)
            else
               if (dev*dcl .lt. zero) dcl = -dcl/3
            end if
         end if
c
c  njtj  ***  modification end  ***
c
c  Find the new potential.
c
   30    continue
         if (down(i)) then
            do 40 k = 2, nr
               viod(lp,k) = viod(lp,k) + dcl*r(k)*f(k)
   40       continue
         else
            do 50 k = 2, nr
               viou(lp,k) = viou(lp,k) + dcl*r(k)*f(k)
   50       continue
         end if
c
c        Update...
c
         cl = cl + dcl
         devold = dev
c
   60 continue
c
c  End the iteration loop for cl.
c
      call ext(820+lp)
c
c   Find the pseudo-wavefunction.
c
   70 continue
      if (down(i)) then
         do 80 j = 2, nr
            v(j) = (viod(lp,j)+llp/r(j))/r(j) + vid(j)
   80    continue
      else
         do 90 j = 2, nr
            v(j) = (viou(lp,j)+llp/r(j))/r(j) + viu(j)
   90    continue
      end if
c
      call difnrl(0,i,v,arps,br,nops(i),lo(i),so(i),ev(i),iflag)
c
c  Compute delta and gamma.
c
      gamma = abs(ar(j3rc)/arps(j3rc)+ar(j3rc+1)/arps(j3rc+1))/2
      ag = zero
      gg = zero
      ll = 4
      do 100 j = 2, nr
         ag = ag + ll*arps(j)*g(j)*rab(j)
         gg = gg + ll*g(j)*g(j)*rab(j)
         ll = 6 - ll
  100 continue
      ag = ag/3
      gg = gg/3
      delta = sqrt((ag/gg)**2+(1/gamma**2-1)/gg) - ag/gg
c
c     Modify the pseudo-wavefunction and pseudo-potential and
c     add to charge density.
c
      if (down(i)) then
         do 110 j = 2, nr
            arps(j) = gamma*(arps(j)+delta*g(j))
            vod(j) = vod(j) + zo(i)*arps(j)*arps(j)
            if (arps(j) .lt. small .and. r(j) .gt. one) arps(j) = small
            rrp = r(j)/rc(lp)
            gpp = (llp-aa*(2*lp+aa-1)*rrp**aa+(aa*rrp**aa)**2)*g(j)/
     &            r(j)**2
            viod(lp,j) = viod(lp,j) + gamma*delta*
     &                   ((ev(i)-v(j))*g(j)+gpp)*r(j)/arps(j)
  110    continue
      else
         do 120 j = 2, nr
            arps(j) = gamma*(arps(j)+delta*g(j))
            vou(j) = vou(j) + zo(i)*arps(j)*arps(j)
            if (arps(j) .lt. small .and. r(j) .gt. one) arps(j) = small
            rrp = r(j)/rc(lp)
            gpp = (llp-aa*(2*lp+aa-1)*rrp**aa+(aa*rrp**aa)**2)*g(j)/
     &            r(j)**2
            viou(lp,j) = viou(lp,j) + gamma*delta*
     &                   ((ev(i)-v(j))*g(j)+gpp)*r(j)/arps(j)
  120    continue
      end if
c
c  wtrans is called to fourier transform the pseudo
c  wave function and save it to the current plot.dat file.
c
      ist = nint(sign(1.d0,arps(nr-85)))
cGuima modifica
      call potrw(arps,r,nr,lo(i),0,ist,rc(lp))
cGuima fim da modificaccao
      call wtrans(arps,r,nr,lo(i),ist)
c
      write(6,9000) nops(i), il(lp), so(i), ev(i), rc(lp), cl, gamma,
     &  delta
 9000 format(1x,i1,a1,f6.1,5f12.6)
c
      return
c
      end
c
c $Id: excorr.f,v 1.3 1999/02/26 14:26:43 wdpgaara Exp $
c
      subroutine excorr(id,cdd,cdu,cdc,vod,vou,vxc,vc,exc,ec)
c
      implicit none
c
      include 'radial.h'
      include 'param.h'
c
c     Compute the LDA exchange and correlation energy and potential.
c
c     Revised by Alberto Garcia
c
c    The only major modification is that the constants for the
c    ceperly-alder 'ca' method are placed in parameter
c    statements, this was done so non-opt compiliers
c    would minimize the number of calculations.
c
C     .. Parameters ..
c
      double precision tiny_charge
      parameter (tiny_charge=1.d-12)
c
      double precision zero, one, pfive, opf, pnn
      parameter (zero=0.D0,one=1.D0,pfive=.5D0,opf=1.5D0,pnn=.99D0)
      double precision pthree, psevf, c0504
      parameter (pthree=0.3D0,psevf=0.75D0,c0504=0.0504D0)
      double precision c0254, c014, c0406
      parameter (c0254=0.0254D0,c014=0.014D0,c0406=0.0406D0)
      double precision c15p9, c0666, c11p4
      parameter (c15p9=15.9D0,c0666=0.0666D0,c11p4=11.4D0)
      double precision c045, c7p8, c88, c20p592
      parameter (c045=0.045D0,c7p8=7.8D0,c88=0.88D0,c20p592=20.592D0)
      double precision c3p52, c0311, c0014
      parameter (c3p52=3.52D0,c0311=0.0311D0,c0014=0.0014D0)
      double precision c0538, c0096, c096
      parameter (c0538=0.0538D0,c0096=0.0096D0,c096=0.096D0)
      double precision c0622, c004, c0232
      parameter (c0622=0.0622D0,c004=0.004D0,c0232=0.0232D0)
      double precision c1686, c1p3981, c2611
      parameter (c1686=0.1686D0,c1p3981=1.3981D0,c2611=0.2611D0)
      double precision c2846, c1p0529, c3334
      parameter (c2846=0.2846D0,c1p0529=1.0529D0,c3334=0.3334D0)
      double precision con1, con2, con3
      parameter (con1=1.D0/6,con2=0.008D0/3,con3=0.3502D0/3)
      double precision con4, con5, con6
      parameter (con4=0.0504D0/3,con5=0.0028D0/3,con6=0.1925D0/3)
      double precision con7, con8, con9
      parameter (con7=0.0206D0/3,con8=9.7867D0/6,con9=1.0444D0/3)
      double precision con10, con11
      parameter (con10=7.3703D0/6,con11=1.3336D0/3)
C     ..
C     .. Scalar Arguments ..
      double precision vxc, vc, exc, ec
      character id*1
C     ..
C     .. Array Arguments ..
      double precision cdd(nrmax), cdu(nrmax), cdc(nrmax), 
     &                 vod(nrmax), vou(nrmax)
c
C     .. Local Scalars ..
      double precision a0, alb, aln, alp,  be, beta,
     &                 cdsum, ecf, ecp, ect, excf, excp,
     &                 exct, exf, exp_var, ftrd, fz, fzp, pi, rs,
     &                 rslog, sb, sqrs, te, tftm, trd, vcd, vcf,
     &                 vcp, vcu, vxcd, vxcf, vxcp, vxcu, vxf, 
     &                 lda_xpot, x, z
      integer i, ll
C     ..
C     .. External Subroutines ..
      external ext
C     ..
C     .. Intrinsic Functions ..
      intrinsic atan, log, sqrt
C     ..
      logical leqi
      external leqi
C     ..
c
      pi = 4*atan(one)
c
c
      trd = one/3
      ftrd = 4*trd
      tftm = 2**ftrd - 2
      a0 = (4/(9*pi))**trd
c
c      set x-alpha
c
      alp = one
      if (.not. leqi(icorr,'xa')) alp = 2*trd
c
c     Initialize
c
      vxc = zero
      vc = zero
      exc = zero
      ec = zero
c
c      start loop (at the second point...see below)
c
      ll = 4
      do 70 i = 2, nr
c
         cdsum = cdd(i) + cdu(i)
         if (ifcore .ge. 1) cdsum = cdsum + cdc(i)
c
         vxcd = 0.d0
         vxcu = 0.d0
         vcd = 0.d0
         vcu = 0.d0
         exct = 0.d0
         ect = 0.d0
c
cag****!!!!!!         if (cdsum .le. tiny_charge) go to 100
c
         rs = (3*r(i)**2/cdsum)**trd
c
c        Spin variables
c
         z = zero
         fz = zero
         fzp = zero
         if (leqi(id,'s')) then
            z = (cdd(i)-cdu(i))/cdsum
            fz = ((1+z)**ftrd+(1-z)**ftrd-2)/tftm
            fzp = ftrd*((1+z)**trd-(1-z)**trd)/tftm
         end if
c
c      exchange (only use (xa))
c
         lda_xpot = -3*alp/(pi*a0*rs)
         exp_var = 3*lda_xpot/4
c
c        Relativistic correction to exchange
c
         if (leqi(id,'r')) then
            beta = c014/rs
            sb = sqrt(1+beta*beta)
            alb = log(beta+sb)
            lda_xpot = lda_xpot*(-pfive+opf*alb/(beta*sb))
            exp_var = exp_var*(one-opf*((beta*sb-alb)/beta**2)**2)
         end if
c
   60    continue
c
         vxf = 2**trd*lda_xpot
         exf = 2**trd*exp_var
         vcp = zero
         ecp = zero
         vcf = zero
         ecf = zero
c
         if (leqi(icorr,'ca')) then
c          ceperly-alder (ca)
c          The Perdew-Zunger parameterization is used.
c          See Phys. Rev. B 23 5075 (1981).
            if (rs .gt. one) then
               sqrs = sqrt(rs)
               te = one + con10*sqrs + con11*rs
               be = one + c1p0529*sqrs + c3334*rs
               ecp = -c2846/be
               vcp = ecp*te/be
               te = one + con8*sqrs + con9*rs
               be = one + c1p3981*sqrs + c2611*rs
               ecf = -c1686/be
               vcf = ecf*te/be
            else
               rslog = log(rs)
               ecp = (c0622+c004*rs)*rslog - c096 - c0232*rs
               vcp = (c0622+con2*rs)*rslog - con3 - con4*rs
               ecf = (c0311+c0014*rs)*rslog - c0538 - c0096*rs
               vcf = (c0311+con5*rs)*rslog - con6 - con7*rs
            end if
c
         else if (leqi(icorr,'xa')) then
c
c          correlation
c
         else if (leqi(icorr,'wi')) then
c
c          wigner (wi)
            vcp = -(c3p52*rs+c20p592)/(3*(rs+c7p8)**2)
            ecp = -c88/(rs+c7p8)
c
         else if (leqi(icorr,'hl')) then
c          hedin-lundqvist (hl)
            x = rs/21
            aln = log(1+1/x)
            vcp = -c045*aln
            ecp = aln + (x**3*aln-x*x) + x/2 - trd
            if (x .gt. 500*one) ecp = ((con1/x-pthree)/x+psevf)/x
            ecp = -c045*ecp
c
         else if (leqi(icorr,'gl')) then
c          gunnarson-lundqvist-wilkins (gl)
            x = rs/c11p4
            aln = log(1+1/x)
            vcp = -c0666*aln
            ecp = aln + (x**3*aln-x*x) + x/2 - trd
            if (x .gt. 500*one) ecp = ((con1/x-pthree)/x+psevf)/x
            ecp = -c0666*ecp
            x = rs/c15p9
            aln = log(1+1/x)
            vcf = -c0406*aln
            ecf = aln + (x**3*aln-x*x) + x/2 - trd
            if (x .gt. 500*one) ecf = ((con1/x-pthree)/x+psevf)/x
            ecf = -c0406*ecf
c
         else if (leqi(icorr,'bh')) then
c          von barth - hedin (bh)
            x = rs/30
            aln = log(1+1/x)
            vcp = -c0504*aln
            ecp = aln + (x**3*aln-x*x) + x/2 - trd
            if (x .gt. 500*one) ecp = ((con1/x-pthree)/x+psevf)/x
            ecp = -c0504*ecp
            x = rs/75
            aln = log(1+1/x)
            vcf = -c0254*aln
            ecf = aln + (x**3*aln-x*x) + x/2 - trd
            if (x .gt. 500*one) ecf = ((con1/x-pthree)/x+psevf)/x
            ecf = -c0254*ecf
c
         else if (leqi(icorr,'gr')) then
c          von barth - hedin + gradient corrections (gr)
            x = rs/30
            aln = log(1+1/x)
            vcp = -c0504*aln
            ecp = aln + (x**3*aln-x*x) + x/2 - trd
            if (x .gt. 500*one) ecp = ((con1/x-pthree)/x+psevf)/x
            ecp = -c0504*ecp
c
         else
c
            write(6,9050) icorr
 9050       format('error in velect - icorr =',a2,' not implemented')
            call ext(400)
c
         end if
c
         vxcp = lda_xpot + vcp
         vxcf = vxf + vcf
         vxcd = vxcp
         vxcu = vxcp
         excp = exp_var + ecp
         excf = exf + ecf
         vcd = vcp
         vcu = vcp
         exct = excp
         ect = ecp
c
         if (z .ne. zero) then
            vxcd = vxcd + fz*(vxcf-vxcp) + (1-z)*fzp*(excf-excp)
            vxcu = vxcu + fz*(vxcf-vxcp) - (1+z)*fzp*(excf-excp)
            vcd = vcd + fz*(vcf-vcp) + (1-z)*fzp*(ecf-ecp)
            vcu = vcu + fz*(vcf-vcp) - (1+z)*fzp*(ecf-ecp)
            exct = exct + fz*(excf-excp)
            ect = ect + fz*(ecf-ecp)
         end if
c
 100     continue
c
         vod(i) = vod(i) + vxcd
         vou(i) = vou(i) + vxcu
c
c        Add to the integrated value
c
         vxc = vxc + ll*(cdd(i)*vxcd+cdu(i)*vxcu)*rab(i)
         vc = vc + ll*(cdd(i)*vcd+cdu(i)*vcu)*rab(i)
         exc = exc + ll*cdsum*exct*rab(i)
         ec = ec + ll*cdsum*ect*rab(i)
         ll = 6 - ll
c
   70 continue
c
      vxc = vxc/3
      vc = vc/3
      ec = ec/3
      exc = exc/3
c
c     Extrapolate backwards for the first point...
c
      vod(1) = vod(2) - (vod(3)-vod(2))*r(2)/(r(3)-r(2))
      vou(1) = vou(2) - (vou(3)-vou(2))*r(2)/(r(3)-r(2))
c
      return
c
      end
c
c $Id: tm2.f,v 1.5 2002/07/04 18:29:34 wdpgaara Exp $
c
      subroutine tm2(i,wfr,br)
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
      include 'energy.h'
c
      include 'tm2_blk.h'
c
      double precision zero, one
      parameter (zero=0.D0,one=1.D0)
      double precision small, ai
      parameter (small=1.D-12,ai=2*137.0360411D0)
      double precision accuracy
      parameter (accuracy=1.d-10)
c
      integer i
      double precision wfr(nrmax), br(nrmax)
c
      double precision work(5)
c
C     .. Local Scalars ..
      double precision arp, bj1, bj2, bj3, bj4, bj5,
     &                 cdps, ddelta, expd, fdnew, fdold, gamma,
     &                 poly, polyr, r2, rc10, rc9, rp, vj, x1, x2,
     &                 xlamda, rcond
      integer ierr, ist, j, k, ka, ll
      logical spin_down, bracketed
C     ..
C     .. Local Arrays ..
      double precision aj(5,5), bj(5)
      double precision aa(nrmax), aap(nrmax), aapp(nrmax), 
     &                 wspl(3*nrmax)
c
      integer indx(5)
C     ..
      integer isrchfgt
      double precision zbrent, v0pp
      external zbrent, v0pp
c
c AG
      fdold = 0.d0
c
      lp = lo(i) + 1
      ka = lo(i) + 1
      eigv = ev(i)
      spin_down = down(i)
      if (spin_down .and. lo(i) .ne. 0) ka = -lo(i)
c
c     Need this to put ar in a common block...
c
      call scopy(nr,wfr,1,ar,1)
c
c     Reset rc to grid point r(j) such that r(j) <= rc < r(j+1)
c
      jrc = isrchfgt(nr,r,1,rc_input(lp)) - 1
      rc(lp) = r(jrc)
c
c  Find the integrated charge inside rc (1-charge outside).
c
      ll = 2
      if (relativistic) then
         cdrc = -(ar(jrc)*ar(jrc)+br(jrc)*br(jrc))*rab(jrc)
         if (mod(jrc,2) .ne. 0) then
            do 40 k = jrc, 1, -1
               cdrc = cdrc + ll*(ar(k)*ar(k)+br(k)*br(k))*rab(k)
               ll = 6 - ll
   40       continue
         else
            do 50 k = jrc, 4, -1
               cdrc = cdrc + ll*(ar(k)*ar(k)+br(k)*br(k))*rab(k)
               ll = 6 - ll
   50       continue
            cdrc = cdrc - (ar(4)*ar(4)+br(4)*br(4))*rab(4)
            cdrc = cdrc + 9*((ar(1)*ar(1)+br(1)*br(1))*rab(1)+
     &             3*(ar(2)*ar(2)+br(2)*br(2))*rab(2)+
     &             3*(ar(3)*ar(3)+br(3)*br(3))*rab(3)+
     &             (ar(4)*ar(4)+br(4)*br(4))*rab(4))/8
         end if
         cdrc = cdrc/3
      else
         cdrc = -ar(jrc)*ar(jrc)*rab(jrc)
         if (mod(jrc,2) .ne. 0) then
            do 60 k = jrc, 1, -1
               cdrc = cdrc + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
   60       continue
         else
            do 70 k = jrc, 4, -1
               cdrc = cdrc + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
   70       continue
            cdrc = cdrc - ar(4)*ar(4)*rab(4)
            cdrc = cdrc + 9*(ar(1)*ar(1)*rab(1)+3*ar(2)*ar(2)*rab(2)+
     &             3*ar(3)*ar(3)*rab(3)+ar(4)*ar(4)*rab(4))/8
         end if
         cdrc = cdrc/3
      end if
c
c  Find the values for wave(arc), d(wave)/dr(arp), potential(vrc),
c  d(potential)/dr(vrp), and d2(potential)/dr2(vrpp)
c
      rc1 = r(jrc)
      rc2 = rc1*rc1
      rc3 = rc2*rc1
      rc4 = rc2*rc2
      rc5 = rc4*rc1
      rc6 = rc4*rc2
      rc7 = rc4*rc3
      rc8 = rc4*rc4
      rc9 = rc4*rc5
      rc10 = rc4*rc6
c
      arc = ar(jrc)
      arp = br(jrc)
c
      if (relativistic) then
         if (spin_down) then
            arp = ka*ar(jrc)/r(jrc) + (eigv-viod(lp,jrc)/r(jrc)-
     &            vid(jrc)+ai*ai)*br(jrc)/ai
         else
            arp = ka*ar(jrc)/r(jrc) + (eigv-viou(lp,jrc)/r(jrc)-
     &            viu(jrc)+ai*ai)*br(jrc)/ai
         end if
      end if
cag?         arp = arp
      brc = arp/arc
c
c
      do 500 j = 2, nr
         if (spin_down) then
            aa(j) = viod(lp,j)/r(j) + vid(j)
         else
            aa(j) = viou(lp,j)/r(j) + viu(j)
         endif
  500 continue   
c
c      Use splines to compute V' and V'' 
c   
      aa(1) = aa(2) - (aa(3)-aa(2))*r(2)/(r(3)-r(2))
      call splift(r,aa,aap,aapp,nr,wspl,ierr,0,zero,zero,zero,zero)
c       
      vrc = aa(jrc)
      vap = aap(jrc)
      vapp = aapp(jrc)
c
c
c   Set up matrix without the d2(potential(0)/dr2=0 condition
c   to find an initial guess for gamma.
c   Note that the equations in the paper have been modified to
c   account for the use of rydberg units.
c
      delta = zero
c
      bj(1) = log(arc/rc1**lp)
      bj1 = bj(1)
      bj(2) = brc - lp/rc1
      bj2 = bj(2)
      bj(3) = vrc - eigv - 2*lp/rc1*bj2 - bj2**2
      bj3 = bj(3)
      bj(4) = vap + 2*lp/rc2*bj2 - 2*lp/rc1*bj3 - 2*bj2*bj3
      bj4 = bj(4)
      bj(5) = vapp - 4*lp/rc3*bj2 + 4*lp/rc2*bj3 - 2*lp/rc1*bj4 -
     &        2*bj3**2 - 2*bj2*bj4
      bj5 = bj(5)
c
      aj(1,1) = rc2
      aj(1,2) = rc4
      aj(1,3) = rc6
      aj(1,4) = rc8
      aj(1,5) = rc10
      aj(2,1) = 2*rc1
      aj(2,2) = 4*rc3
      aj(2,3) = 6*rc5
      aj(2,4) = 8*rc7
      aj(2,5) = 10*rc9
      aj(3,1) = 2*one
      aj(3,2) = 12*rc2
      aj(3,3) = 30*rc4
      aj(3,4) = 56*rc6
      aj(3,5) = 90*rc8
      aj(4,1) = zero
      aj(4,2) = 24*rc1
      aj(4,3) = 120*rc3
      aj(4,4) = 336*rc5
      aj(4,5) = 720*rc7
      aj(5,1) = zero
      aj(5,2) = 24*one
      aj(5,3) = 360*rc2
      aj(5,4) = 1680*rc4
      aj(5,5) = 5040*rc6
c
c     Use LU decomposition to solve the linear system of
c     equations. We can re-use the decomposition inside the
c     loop. This is not the case if gaussian elimination is
c     used ! (Alberto Garcia, April 1991) See Numerical Recipes, p. 31.
c
      call sgeco(aj,5,5,indx,rcond,work)
      if (rcond .lt. 1.d-7) write(6,*) ' rcond too small:' ,rcond
c
      call sgesl(aj,5,5,indx,bj,0)
c
      gamma = bj(1)
      alpha = bj(2)
      alpha1 = bj(3)
      alpha2 = bj(4)
      alpha3 = bj(5)
c
c  Start iteration loop to find delta, uses false position.
c
      do 150 j = 1, 50
c
c  Generate pseudo wavefunction-note missing factor exp(delta).
c
         do 100 k = 1, jrc
            rp = r(k)
            r2 = rp*rp
            polyr = r2*((((alpha3*r2+alpha2)*r2+alpha1)*r2+alpha)*r2+
     &              gamma)
            ar(k) = rp**lp*exp(polyr)
  100    continue
c
c           Integrate pseudo charge density from r = 0 to rc.
c
         ll = 2
         cdps = -ar(jrc)*ar(jrc)*rab(jrc)
         if (mod(jrc,2) .ne. 0) then
            do 110 k = jrc, 1, -1
               cdps = cdps + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
  110       continue
         else
            do 120 k = jrc, 4, -1
               cdps = cdps + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
  120       continue
            cdps = cdps - ar(4)*ar(4)*rab(4)
            cdps = cdps + 9*(ar(1)*ar(1)*rab(1)+3*ar(2)*ar(2)*rab(2)+
     &             3*ar(3)*ar(3)*rab(3)+ar(4)*ar(4)*rab(4))/8
         end if
         cdps = cdps/3
c
c           Calculate new delta
c
         fdnew = log(cdrc/cdps) - 2*delta
         if (abs(fdnew) .lt. small) go to 160
         if (j .eq. 1) then
            ddelta = -one/2
         else
            ddelta = -fdnew*ddelta/(fdnew-fdold)
         endif
         delta = delta + ddelta
c
         bj(1) = bj1 - delta
         bj(2) = bj2
         bj(3) = bj3
         bj(4) = bj4
         bj(5) = bj5
c
         call sgesl(aj,5,5,indx,bj,0)
c
         gamma = bj(1)
         alpha = bj(2)
         alpha1 = bj(3)
         alpha2 = bj(4)
         alpha3 = bj(5)
c
         fdold = fdnew
c
  150 continue
c
c  End iteration loop for delta.
c
      write(6,9000) lp - 1
      call ext(820+lp)
 9000 format(//'error in pseud2 - nonconvergence in finding',
     &      /' starting delta for angular momentum ',i1)
c
c  Bracket the correct gamma, use gamma and -gamma
c  from above as initial brackets, expands brackets
c  until a root is found..
c
  160 continue
      alpha4 = zero
      x1 = gamma
      x2 = -gamma
c
      call zbrac(v0pp,x1,x2,bracketed)
c
      if ( .not. bracketed) then
         write(6,9010) lp
         call ext(830+lp)
 9010    format(//'Error in zbractk - cannot bracket orbital ',i2)
      end if
c
c  Iteration loop to find correct gamma, uses
c  bisection to find gamma.
c
      gamma = zbrent(v0pp,x1,x2,accuracy)
c
c  Augment charge density and invert schroedinger equation
c  to find new potential.
c
  170 continue
      expd = exp(delta)
         do 180 j = 1, jrc
            r2 = r(j)*r(j)
            poly = r2*(((((alpha4*r2+alpha3)*r2+alpha2)*r2+alpha1)*r2+
     &             alpha)*r2+gamma)
            ar(j) = r(j)**lp*expd*exp(poly)
c
c           For those of us with inferior minds, xlamda = p'(r)/r
c           and the thing goes like this (eq. 23 in TM2, in rydberg):
c
c                                               2
c           V = eigv + 2 p'* l(l+1)/r + p'' + p'   ==>
c
c           V = eigv + (p'/r) * [ l(l+1) + rp' ] + p''
c                        |                  |
c                     "xlamda"         "xlamda*r2"
c
            xlamda = ((((12*alpha4*r2+10*alpha3)*r2+8*alpha2)*r2+
     &               6*alpha1)*r2+4*alpha)*r2 + 2*gamma
c
            vj = eigv + xlamda*(2*lp+xlamda*r2) +
     &           ((((132*alpha4*r2+90*alpha3)*r2+56*alpha2)*r2+
     &           30*alpha1)*r2+12*alpha)*r2 + 2*gamma
c
          if (spin_down) then
            viod(lp,j) = (vj-vid(j))*r(j)
            vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
          else
            viou(lp,j) = (vj-viu(j))*r(j)
            vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
          endif
c
  180    continue
c
         do 190 j = jrc + 1, nr
          if (spin_down) then
            vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
          else
            vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
          endif
  190    continue
c
c  njtj  ***  plotting routines ***
c  potrw is called to save a useful number of points
c  of the pseudowave function to make a plot.  The
c  info is written to the current plot.dat file.
c  wtrans is called to fourier transform the the pseudo
c  wave function and save it to the current plot.dat file.
c
      ist = 1
cGuima modifica
      call potrw(ar,r,nr,lo(i),0,ist,rc(lp))
cGuima fim da modificaccao
      call wtrans(ar,r,nr,lo(i),ist)
c
      write(6,9020) lp, il(lp), so(i), eigv, rc(lp), cdrc, delta
 9020 format(1x,i1,a1,f6.1,5f12.6)
c
      return
c
      end
c
c $Id: ker.f,v 1.2 1997/05/22 17:32:17 wdpgaara Exp $
c
c $Log: ker.f,v $
c Revision 1.2  1997/05/22 17:32:17  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine ker(i,ar,br)
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
      include 'energy.h'
c
      double precision zero, one, pfive, smtol
      parameter (zero=0.D0,one=1.D0,pfive=0.5D0,smtol=1.D-12)
      double precision small, ai
      parameter (small=1.D-12,ai=2*137.0360411D0)
c
      integer i
      double precision ar(nrmax), br(nrmax)
c
      integer ist, iswtch, j, jrc, lp, ka, ll, k
      double precision polyr, eigv, cdrc, rc2, rc3, rc4, arc, arp, 
     &                 brc, vrc, alpha, beta, gamma, delta, cdps,
     &                 fdnew, ddelta, fdold, expd, xlamda, vj
c
      integer isrchfgt
c
      lp = lo(i) + 1
      ka = lo(i) + 1
      eigv = ev(i)
      if (down(i) .and. lo(i) .ne. 0) ka = -lo(i)
c
c     Reset rc to grid point r(j) such that r(j) <= rc < r(j+1)
c
      jrc = isrchfgt(nr,r,1,rc_input(lp)) - 1
      rc(lp) = r(jrc)
c
c  Find the integrated charge inside rc.
c
      ll = 2
      if (relativistic) then
         cdrc = -(ar(jrc)*ar(jrc)+br(jrc)*br(jrc))*rab(jrc)
         if (mod(jrc,2) .ne. 0) then
            do 40 k = jrc, 1, -1
               cdrc = cdrc + ll*(ar(k)*ar(k)+br(k)*br(k))*rab(k)
               ll = 6 - ll
   40       continue
         else
            do 50 k = jrc, 4, -1
               cdrc = cdrc + ll*(ar(k)*ar(k)+br(k)*br(k))*rab(k)
               ll = 6 - ll
   50       continue
            cdrc = cdrc - (ar(4)*ar(4)+br(4)*br(4))*rab(4)
            cdrc = cdrc + 9*((ar(1)*ar(1)+br(1)*br(1))*rab(1)+
     &             3*(ar(2)*ar(2)+br(2)*br(2))*rab(2)+
     &             3*(ar(3)*ar(3)+br(3)*br(3))*rab(3)+
     &             (ar(4)*ar(4)+br(4)*br(4))*rab(4))/8
         end if
         cdrc = cdrc/3
      else
         cdrc = -ar(jrc)*ar(jrc)*rab(jrc)
         if (mod(jrc,2) .ne. 0) then
            do 60 k = jrc, 1, -1
               cdrc = cdrc + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
   60       continue
         else
            do 70 k = jrc, 4, -1
               cdrc = cdrc + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
   70       continue
            cdrc = cdrc - ar(4)*ar(4)*rab(4)
            cdrc = cdrc + 9*(ar(1)*ar(1)*rab(1)+3*ar(2)*ar(2)*rab(2)+
     &             3*ar(3)*ar(3)*rab(3)+ar(4)*ar(4)*rab(4))/8
         end if
         cdrc = cdrc/3
      end if
c
c   The initial values for alpha, beta, gamma and delta.
c
      rc2 = r(jrc)*r(jrc)
      rc3 = r(jrc)*rc2
      rc4 = r(jrc)*rc3
      iswtch = 1
      if (ar(jrc) .lt. zero) iswtch = -1
      arc = iswtch*ar(jrc)
      arp = br(jrc)
c
      if (relativistic) then
         if (down(i)) then
            arp = ka*ar(jrc)/r(jrc) + (eigv-viod(lp,jrc)/r(jrc)-
     &            vid(jrc)+ai*ai)*br(jrc)/ai
         else
            arp = ka*ar(jrc)/r(jrc) + (eigv-viou(lp,jrc)/r(jrc)-
     &            viu(jrc)+ai*ai)*br(jrc)/ai
         end if
      end if
c
      brc = arp/ar(jrc)
      if (down(i)) then
         vrc = viod(lp,jrc)/r(jrc) + vid(jrc)
      else
         vrc = viou(lp,jrc)/r(jrc) + viu(jrc)
      end if
      alpha = (3*log(arc/r(jrc)**lp)-2*(r(jrc)*brc-lp)+
     &        (rc2*vrc+lp*lp-rc2*(eigv+brc*brc))/2)/rc4
      beta = (-8*log(arc/r(jrc)**lp)+5*(r(jrc)*brc-lp)-
     &       (rc2*vrc+lp*lp-rc2*(eigv+brc*brc)))/rc3
      gamma = (6*log(arc/r(jrc)**lp)-3*(r(jrc)*brc-lp)+
     &        (rc2*vrc+lp*lp-rc2*(eigv+brc*brc))/2)/rc2
      delta = zero
c
c  Start the iteration loop to find delta.
c
      do 110 j = 1, 50
c
c  Generate the pseudo-wavefunction (note missing factor exp(delta)).
c
         do 80 k = 1, jrc
            polyr = r(k)*r(k)*((alpha*r(k)+beta)*r(k)+gamma)
            ar(k) = iswtch*r(k)**lp*exp(polyr)
   80    continue
c
c  Integrate  the pseudo charge density from r = 0 to rc.
c
         ll = 2
         cdps = -ar(jrc)*ar(jrc)*rab(jrc)
         if (mod(jrc,2) .ne. 0) then
            do 90 k = jrc, 1, -1
               cdps = cdps + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
   90       continue
         else
            do 100 k = jrc, 4, -1
               cdps = cdps + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
  100       continue
            cdps = cdps - ar(4)*ar(4)*rab(4)
            cdps = cdps + 9*(ar(1)*ar(1)*rab(1)+3*ar(2)*ar(2)*rab(2)+
     &             3*ar(3)*ar(3)*rab(3)+ar(4)*ar(4)*rab(4))/8
         end if
         cdps = cdps/3
c
c  Find the new delta.
c
         fdnew = log(cdrc/cdps) - 2*delta
         if (abs(fdnew) .lt. smtol) go to 120
         if (j .eq. 1) then
            ddelta = pfive
         else
            ddelta = -fdnew*ddelta/(fdnew-fdold)
         end if
         alpha = alpha - 3*ddelta/rc4
         beta = beta + 8*ddelta/rc3
         gamma = gamma - 6*ddelta/rc2
         delta = delta + ddelta
         fdold = fdnew
  110 continue
c
c  End the iteration loop for delta.
c
      call ext(820+lp)
c
c    Augment the charge density and invert schroedinger equation
c  to find new potential.
c
  120 continue
      expd = exp(delta)
      if (down(i)) then
         do 130 j = 1, jrc
            ar(j) = expd*ar(j)
            xlamda = (4*alpha*r(j)+3*beta)*r(j) + 2*gamma
            vj = eigv + xlamda*(2*lp+xlamda*r(j)**2) +
     &           (12*alpha*r(j)+6*beta)*r(j) + 2*gamma
            viod(lp,j) = (vj-vid(j))*r(j)
            vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
  130    continue
         do 140 j = jrc + 1, nr
            vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
  140    continue
      else
         do 150 j = 1, jrc
            ar(j) = expd*ar(j)
            xlamda = (4*alpha*r(j)+3*beta)*r(j) + 2*gamma
            vj = eigv + xlamda*(2*lp+xlamda*r(j)**2) +
     &           (12*alpha*r(j)+6*beta)*r(j) + 2*gamma
            viou(lp,j) = (vj-viu(j))*r(j)
            vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
  150    continue
         do 160 j = jrc + 1, nr
            vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
  160    continue
      end if
      write(6,9000) lo(i)+1, il(lp), so(i), eigv, rc(lp), cdrc, delta
 9000 format(1x,i1,a1,f6.1,5f12.6)
c
c  njtj  ***  plotting routines ***
c  potrw is called to save a usefull number of points
c  of the pseudowave function to make a plot.  The
c  info is written to the current plot.dat file.
c  wtrans is called to fourier transform the the pseudo
c  wave function and save it to the current plot.dat file.
c
      ist = 1
      if (ar(nr-85) .lt. zero) ist = -1
cGuima modifica
      call potrw(ar,r,nr,lo(i),0,ist,rc(lp))
cGuima fim da modificaccao
      call wtrans(ar,r,nr,lo(i),ist)
c
      return
c
      end
c
      subroutine wrapup(pot_id)
c
      implicit none

c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'charge.h'
      include 'elecpot.h'
      include 'energy.h'
      include 'compat.h'
c
c     ecuts now set in compat_params...
c
      double precision zero, tpfive, one
      parameter (zero=0.d0,tpfive=2.5d0, one=1.d0)
c
      integer i, icore, j, jcut, lp, noi, npotd, npotu, ifull
      integer nops(norbmx), position
      character ray(6)*10, title*70, pot_id*40, id*1
c
      double precision zval, zratio, zion, ac, bc, cdcp, tanb, rbold,
     &                 rbnew, pi, ecut, vp2z, fcut, zot, vpsdm,
     &                 vps, rmind, vpsum, rminu, zelu, zeld, zelt,
     &                 viodj, viouj, cc
      double precision rcut(10), v(nrmax)
c
      double precision absval, minabs, maxabs, norm1, norm2
      double precision fourier_area(5), dummy_real
      integer n_channels, lun
cGuima modifica
      logical lexist
      double precision adicao(lmax,nrmax)
cGuima end modifica
      double precision cutoff_function
      external cutoff_function

      logical new_scheme
c
      external logder
c
      pi = 4*atan(one)
c
c     Do not use relativity for what follows
c
      if (polarized) then
         id = 's'
      else
         id = ' '
      endif
c
c  Reset the n quantum numbers to include all valence orbitals.
c  Compute the ratio between the valence charge present and the
c  valence charge of a neutral atom.
c  Transfer pseudo valence charge to charge array
c
      zval = zero
      zratio = zero
c
      do 5 i = 1, norb
         nops(i) = 0
    5 continue
c
      do 10 i = ncp, norb
         nops(i) = lo(i) + 1
         zval = zval + zo(i)
   10 continue
      zion = zval + znuc - zel
      if (zval .ne. zero) zratio = zion/zval
c
      do 20 i = 1, nr
         cdd(i) = vod(i)
         cdu(i) = vou(i)
   20 continue
c
c=====================================================================
c  If a core correction is indicated construct pseudo core charge
c  if cfac < 0 or the valence charge is zero the full core is used
c
      ifcore = job - 1
      if (ifcore .ne. 0) then
         ac = zero
         bc = zero
         cc = zero
         icore = 1
         if (cfac .le. zero .or. zratio .eq. zero) then
            write(6,9000) r(icore), ac, bc, cc
            write(6,'(a)') '(Full core used)'
            call coreq
         else
            if (rcfac .le. zero) then
               do 30 i = nr, 2, -1
                  if (cdc(i) .gt. cfac*zratio*
     &                (cdd(i)+cdu(i))) go to 50
   30          continue
            else
               do 40 i = nr, 2, -1
                  if (r(i) .le. rcfac) go to 50
   40          continue
            end if
   50       continue
            icore = i
C
C--------------------------------------------------------------------
C Choice of methods to compute the core correction:
C
C 1. Traditional 'Froyen-Louie-Cohen' with cdc(r) = Arsin(Br)
C    and value and first-derivative matching. It is the default
c    for LDA calculations. See compat_params.f
C    and input.f for info on how to force its use from the input file.
C
C 2. New by Jose Luis Martins' group, using a Kerker-like exp( ) function
C    and matching also the second derivative. This is the default with
C    the 'mons' compatibility mode for GGA calculations.
C

            new_scheme = .false.
            if (is_gga) then
               write(6,'(a)') 'Note: GGA calculation ==> New CC scheme'
               new_scheme = .true.
            endif

            if ( (use_old_cc) .and. (is_gga)) then

               write(6,'(/,2a,/a)')
     $              'WARNING: Using old-style core corrections',
     $            ' despite this being a GGA calculation.',
     $              'I hope you know what you are doing...'
               new_scheme = .false.
            endif

            if ( (use_new_cc) .and. (.not. is_gga)) then

               write(6,'(/,2a,/a)')
     $              'WARNING: Using new core corrections',
     $            ' despite this being an LDA calculation.',
     $         'Results will not be compatible with older versions.'
               new_scheme = .true.
            endif

            if (.not. new_scheme) then

c           Fit to  cdc(r) = ac*r * sin(bc*r) inside r(icore)
c
c           Find derivative at core radius. Use a five-point formula
c           instead of the old two-point method:
c              cdcp = (cdc(icore+1)-cdc(icore))/(r(icore+1)-r(icore))
c           Since the r abscissae are not equally spaced, we perform
c           the derivative using the chain rule:
c
c           r(i) = a [ exp(b*(i-1)) - 1 ]
c           r(x) = a [ exp(b*x) - 1 ]
c
c           f'(r) = f'(x) / r'(x)
c
c           r'(x) = a*b*exp(b*x) = (r(x)+a)*b
c           r'(i) = a*b*exp[b*(i-1)] = rab(i)
c
c           To compute f'(x), use eq. 25.3.6 
c           (p. 883) in Abramowitz-Stegun, with h=1
c
c           f'(0) = 1/12 f(-2) - 2/3 f(-1) + 2/3 f(1) - 1/12 f(2)
c
cold            cdcp = (cdc(icore+1)-cdc(icore))/(r(icore+1)-r(icore))

            cdcp =   cdc(icore-2) - 8*cdc(icore-1) +
     $             8*cdc(icore+1) -   cdc(icore+2)
            cdcp = cdcp / 12.d0 / rab(icore)
c
c           Now fit ac and bc using the function and derivative
c           information. 
c
c           g(r) = Arsin(Br) ==>  g / (rg'-g) = tanBr/(Br)
c
c           Use an iterative method to find Br and from that ac and bc.
c           Start near Br=2.5 so that we are in an invertible region 
c           of the tanx/x function.
c 
c
            tanb = cdc(icore)/(r(icore)*cdcp-cdc(icore))
            rbold = tpfive
            do 70 i = 1, 50
               rbnew = pi + atan(tanb*rbold)
               if (abs(rbnew-rbold) .lt. .00001D0) then
                  bc = rbnew/r(icore)
                  ac = cdc(icore)/(r(icore)*sin(rbnew))
                  do 60 j = 1, icore
                     cdc(j) = ac*r(j)*sin(bc*r(j))
   60             continue
                  write(6,9000) r(icore), ac, bc, cc
c
                  call coreq
                  go to 80
c
               else
                  rbold = rbnew
               end if
   70       continue
            write(6,9010)
            call ext(830)

            else

c                 Use subroutine provided by JLM to fit
c                 cdc(r) = r^2*exp(ac+bc*r^2+cc*r^4) inside r(icore)
c
                  CALL PCC_EXP(NR,ICORE,AC,BC,CC,R,CDC)
                  write(6,9000) r(icore), ac, bc, cc

c
           endif

           call coreq

         end if
      end if
C---------------------------------------------------------------------
 9000 format(//' Core correction used',/' Pseudo core inside r =',f6.3,
     &      /' ac =',f6.3,' bc =',f6.3,' cc =',f6.3,/)
 9010 format(//' Error in pseudo - nonconvergence in finding ',
     &      /'pseudo-core values')
c
c  End the pseudo core charge.
c======================================================================
c
c  Compute the potential due to pseudo valence charge.
c
c  njtj  ***  NOTE  ***
c  Spin-polarized potentials should be unscreened with
c  spin-polarized valence charge.  This was not
c  done in pseudo and pseudok in earlier versions
c  of this program.
c  njtj  ***  NOTE  ***
c
   80 continue
c
      call Velect(0,1,id,zval)
c
c  Construct the ionic pseudopotential and find the cutoff,
c  ecut should be adjusted to give a reassonable ionic cutoff
c  radius, but should not alter the pseudopotential, ie.,
c  the ionic cutoff radius should not be inside the pseudopotential
c  cutoff radius.
c
c  Note that the cutting off of the pseudopotentials (making
c  them approach -2*Zion/r faster) is not strictly necessary.
c  It might even be argued that it should be left to "client"
c  programs to decide what to do.

cag
c
c     On the issue of plotting:
c
c     For non-relativistic, non-spin-polarized calculations, all
c     the orbitals are considered as "down".
c     For relativistic calculations, the "s" orbitals are considered
c     as "up". The actual things plotted on the files (which only
c     record l, not the up/down character) depend on the order of
c     enumeration of the orbitals. Since these are always "down/up",
c     the potentials plotted are always the "up" ones (except of
c     course for scalar calculations and for "s" states in relativistic
c     calculations, for which the distinction is irrelevant.

      write(6,9020)
 9020 format(/)
      ecut = ecuts
      do 150 i = ncp, norb
         lp = lo(i) + 1
         if (down(i)) then
            do 90 j = 2, nr
               v(j) = viod(lp,j)/r(j) + vid(j)
               viod(lp,j) = viod(lp,j) + (vid(j)-vod(j))*r(j)
               vp2z = viod(lp,j) + 2*zion
               if (abs(vp2z) .gt. ecut) jcut = j
   90       continue
cag
c           Plot screened ionic potential
c
            call potrvs(v,r,nr-120,lo(i))
cag
c           Default cutoff function: f(r)=exp(-5*(r-r_cut)). It damps
c           down the residual of rV+2*Zion.
c           Should be made smoother... Vps ends up with a kink at rcut.
c           Maybe use one of the Vanderbilt generalized gaussians.
cag
            rcut(i-ncore) = r(jcut)
            if (rcut(i-ncore) .lt. rc(lp)) then
               write(6,'(a,2f8.4)') 'Vps rcut point moved out to rc: ',
     $              rcut(i-ncore), rc(lp)
               rcut(i-ncore) = rc(lp)
            endif
            do 100 j = jcut, nr
cag               fcut = exp(-5*(r(j)-r(jcut)))
               fcut = cutoff_function(r(j)-r(jcut))
               viod(lp,j) = -2*zion + fcut*(viod(lp,j)+2*zion)
  100       continue
            do 110 j = 2, nr
               v(j) = viod(lp,j)/r(j)
  110       continue
c
            call potran(lo(i)+1,v,r,nr,zion,dummy_real)
            call potrv(v,r,nr-120,lo(i),zion)
c
         else
            do 120 j = 2, nr
               v(j) = viou(lp,j)/r(j) + viu(j)
               viou(lp,j) = viou(lp,j) + (viu(j)-vou(j))*r(j)
               vp2z = viou(lp,j) + 2*zion
               if (abs(vp2z) .gt. ecut) jcut = j
  120       continue
cag
c           Plot screened ionic potential
c
            call potrvs(v,r,nr-120,lo(i))
cag
            rcut(i-ncore) = r(jcut)
            if (rcut(i-ncore) .lt. rc(lp)) then
               write(6,'(a,2f8.4)') 'Vps rcut point moved out to rc: ',
     $              rcut(i-ncore), rc(lp)
               rcut(i-ncore) = rc(lp)
            endif
            do 130 j = jcut, nr
cag               fcut = exp(-5*(r(j)-r(jcut)))
               fcut = cutoff_function(r(j)-r(jcut))
               viou(lp,j) = -2*zion + fcut*(viou(lp,j)+2*zion)
  130       continue
            do 140 j = 2, nr
               v(j) = viou(lp,j)/r(j)
  140       continue
c
            call potran(lo(i)+1,v,r,nr,zion,fourier_area(lo(i)+1))
            call potrv(v,r,nr-120,lo(i),zion)
c
         end if
c
  150 continue
c
c     Write out the Fourier area for each pseudo channel
c
      call get_unit(lun)
      open(unit=lun,file="FOURIER_AREA",form="formatted",
     $     status="unknown")
      rewind(lun)
c
c     Compute also the minimum, maximum, mean, and root-mean-square.
c
      n_channels = 0
      maxabs = -1.0d0
      minabs = 1.0d10
      norm1 = 0.0d0
      norm2 = 0.0d0
      do j=lo(ncp) + 1, lo(norb) + 1
         absval = fourier_area(j)
         if (absval .gt. maxabs) maxabs = absval
         if (absval .lt. minabs) minabs = absval
         norm1 = norm1 + absval
         norm2 = norm2 + absval*absval
         n_channels =  n_channels + 1
      enddo
      norm1 = norm1 / n_channels
      norm2 = sqrt( norm2 / n_channels)
      write(lun,"(i4)") n_channels
      write(lun,"(5f10.5)") (fourier_area(j),j=lo(ncp)+1,lo(norb)+1)
      write(lun,"(4f10.5)") minabs, maxabs, norm1, norm2
      close(lun)
c
      write(6,9020)
c
c   Convert spin-polarized potentials back to nonspin-polarized
c   by occupation weight(zo).  Assumes core polarization is
c   zero, ie. polarization is only a valence effect.
c
      if (polarized) then
         do 180 i = ncp, norb, 2
            lp = lo(i) + 1
            zot = zo(i) + zo(i+1)
            if (zot .ne. zero) then
               do 160 j = 2, nr
                  viod(lp,j) = (viod(lp,j)*zo(i)+viou(lp,j)*zo(i+1))/zot
                  viou(lp,j) = viod(lp,j)
  160          continue
            else
               do 170 j = 2, nr
                  viod(lp,j) = viod(lp,j)/2 + viou(lp,j)/2
                  viou(lp,j) = viod(lp,j)
  170          continue
            end if
  180    continue
      end if
c
      do 190 i = 2, nr
         vid(i) = vod(i)
         viu(i) = vou(i)
  190 continue
c
c   Test the pseudopotential self consistency.  Spin-polarized
c   is tested as spin-polarized(since up/down potentials are
c   now the same)
c
      call dsolv2(0,1,id,ncp,norb,0,nops)
c
c  Printout the pseudo eigenvalues after cutoff.
c
      write(6,9030) (il(lo(i)+1),rcut(i-ncore),i=ncp,norb)
      write(6,9040) (ev(i),i=ncp,norb)
 9030 format(//' test of eigenvalues',//' rcut =',8(2x,a1,f7.2))
 9040 format(' eval =',8(2x,f8.5))
c
c  Printout the data for potentials.
c
      write(6,9050)
 9050 format(///' l    vps(0)    vpsmin      at r',/)
c
      do 220 i = 1, lmax
         if (indd(i)+indu(i) .eq. 0) go to 220
         if (indd(i) .ne. 0) then
            vpsdm = zero
            do 200 j = 2, nr
               if (r(j) .lt. .00001D0) go to 200
               vps = viod(i,j)/r(j)
               if (vps .lt. vpsdm) then
                  vpsdm = vps
                  rmind = r(j)
               end if
  200       continue
            write(6,9060) il(i), viod(i,2)/r(2), vpsdm, rmind
         end if
         if (indu(i) .ne. 0) then
            vpsum = zero
            do 210 j = 2, nr
               if (r(j) .lt. .00001D0) go to 210
               vps = viou(i,j)/r(j)
               if (vps .lt. vpsum) then
                  vpsum = vps
                  rminu = r(j)
               end if
  210       continue
            write(6,9060) il(i), viou(i,2)/r(2), vpsum, rminu
         end if
 9060    format(1x,a1,3f10.3)
  220 continue
c
c   Print out the energies from etotal. (Valence only...)
c
      call etotal(ncp,norb)
c
c   Compute the logarithmic derivative as a function of energy 
c
      if (logder_radius .gt. 0.d0) call logder(ncp,norb,'PS')
c
c  Find the jobname and date.
c
      ray(1) = 'ATM 3.2.2'
c      call cal_date(ray(2))
c  
      read(pot_id,'(4a10)') (ray(i),i=3,6)
c
c  Encode the title array.
c
      title = ' '
      position = 1
      do 240 i = 1, lmax
         if (indd(i) .eq. 0 .and. indu(i) .eq. 0) go to 240
         zelu = zero
         zeld = zero
         if (indd(i) .ne. 0) then
            noi = no(indd(i))
            zeld = zo(indd(i))
         end if
         if (indu(i) .ne. 0) then
            noi = no(indu(i))
            zelu = zo(indu(i))
         end if
         zelt = zeld + zelu
         if ( .not. polarized) then
            write(title(position:),9070) noi, il(i), zelt, ispp, rc(i)
 9070       format(i1,a1,f5.2,a1,' r=',f5.2,'/')
            position = position + 17
         else
            write(title(position:),9090)
     $                             noi, il(i), zeld, zelu, ispp, rc(i)
 9090       format(i1,a1,f4.2,',',f4.2,a1,f4.2,'/')
            position = position + 17
         end if
  240 continue
c
c  Construct relativistic sum and difference potentials.
c
      if (relativistic) then
         if (indu(1) .eq. 0) go to 260
         indd(1) = indu(1)
         indu(1) = 0
         do 250 j = 2, nr
            viod(1,j) = viou(1,j)
            viou(1,j) = zero
  250    continue
  260    continue
         do 280 i = 2, lmax
            if (indd(i) .eq. 0 .or. indu(i) .eq. 0) go to 280
            do 270 j = 2, nr
               viodj = viod(i,j)
               viouj = viou(i,j)
               viod(i,j) = ((i-1)*viodj+i*viouj)/(2*i-1)
               viou(i,j) = 2*(viouj-viodj)/(2*i-1)
  270       continue
  280    continue
      end if
c
c  Determine the number of  potentials.  Coded them as
c  two digits, where the first digit is the number
c  of down or sum potentials and the second the number of
c  up or difference potentials.
c
      npotd = 0
      npotu = 0
      do 290 i = 1, lmax
         if (indd(i) .ne. 0) npotd = npotd + 1
         if (indu(i) .ne. 0) npotu = npotu + 1
  290 continue
c
c  Write the heading to the current pseudo.dat
c  file (unit=1).
c
      ifull = 0
      if (cfac .le. zero .or. zratio .eq. zero) ifull = 1
      if (ifcore .eq. 1) then
         if (ifull .eq. 0) then
            nicore = 'pcec'
         else
            nicore = 'fcec'
         end if
      else if (ifcore .eq. 2) then
         if (ifull .eq. 0) then
            nicore = 'pche'
         else
            nicore = 'fche'
         end if
      else
         nicore = 'nc  '
      end if
      if (polarized) then
         irel = 'isp'
      else if (relativistic) then
         irel = 'rel'
      else
         irel = 'nrl'
      end if

      open(unit=1,file='VPSIN',status='unknown',form='unformatted')
      rewind 1
      open(unit=2,file='VPSFMT',status='unknown',form='formatted')
      rewind 2
      open(unit=22,file='VTOTALps',status='unknown')
      rewind 22
      write(1) nameat, icorr, irel, nicore, (ray(i),i=1,6), title,
     &         npotd, npotu, nr - 1, a, b, zion
      write(1) (r(i),i=2,nr)
c
      write(2,8005) nameat, icorr, irel, nicore
      write(2,8010) (ray(j),j=1,6), title
      write(2,8015) npotd, npotu, nr-1, a, b, zion
      write(2,8040) 'Radial grid follows'
      write(2,8030) (r(j),j=2,nr)
c
 8000 format(1x,i2)
 8005 format(1x,a2,1x,a2,1x,a3,1x,a4)
 8010 format(1x,6a10,/,1x,a70)
 8015 format(1x,2i3,i5,3g20.12)
 8030 format(4(g20.12))
 8040 format(1x,a)
c
c  Write the potentials to file (unit=1).
c
cGuima inserccao
      inquire(file='adiciona',exist=lexist)
      if(lexist)  open(unit=21,file='adiciona')
cGuima termina inserccao
      do 300 i = 1, lmax
         if (indd(i) .eq. 0) go to 300
cGuima inserccao
         if(lexist)then
            read(21,*)
            read(21,*)
            read(21,*)(adicao(i,j),j=2,nr)
            do j=2,nr
               viod(i,j)=viod(i,j)+adicao(i,j)
            enddo
         endif
cGuima termina inserccao
         write(1) i - 1, (viod(i,j),j=2,nr)
         write(2,8040) 'Down Pseudopotential follows (l on next line)'
         write(2,8000) i-1
         write(2,8030) (viod(i,j),j=2,nr)
cGuima inserccao
         write(22,8040) 'Down potential follows (l on next line)'
         write(22,8000) i-1
         write(22,8030) (viod(i,j)+vid(j)*r(j),j=2,nr)
cGuima termina inserccao
  300 continue
      do 310 i = 1, lmax
         if (indu(i) .eq. 0) go to 310
cGuima inserccao
         if(lexist)then
            read(21,*)
            read(21,*)
            read(21,*)(adicao(i,j),j=2,nr)
            do j=2,nr
               viou(i,j)=viou(i,j)+adicao(i,j)
            enddo
         endif
cGuima termina inserccao
         write(1) i - 1, (viou(i,j),j=2,nr)
         write(2,8040) 'Up Pseudopotential follows (l on next line)'
         write(2,8000) i-1
         write(2,8030) (viou(i,j),j=2,nr)
cGuima inserccao
         write(22,8040) 'Up potential follows (l on next line)'
         write(22,8000) i-1
         write(22,8030) (viou(i,j)+viu(j)*r(j),j=2,nr)
cGuima termina inserccao
 310  continue
c
c  Write the charge densities to units 1 and 3(formatted).
c  Note that this charge density is the "pseudo" one.
c
      write(2,8040) 'Core charge follows'

      if (ifcore .ne. 1) then
         write(1) (zero,i=2,nr)
         write(2,8030) (zero,i=2,nr)
      else
         write(1) (cdc(i),i=2,nr)
         write(2,8030) (cdc(i),i=2,nr)
      end if
      write(1) (zratio*(cdd(i)+cdu(i)),i=2,nr)
      write(2,8040) 'Valence charge follows'
      write(2,8030) (zratio*(cdd(i)+cdu(i)),i=2,nr)

      close(1)
      close(2)
c
      open(unit=3,file='PSCHARGE',form='formatted',status='unknown')
c
c     NOTE: We no longer put "zratio" here!!!
c     (We still do in the ps file for compatibility with PW and
c      SIESTA) 
c     (Only affects plots for ionic configurations)
c     BEAR THIS IN MIND IF YOU ARE USING THE HEURISTIC CORE CORRECTION
c     CRITERION: If you specify a given "pc_weight" in the input file,
c     do not be surprised if the plot does not show rcore
c     in the place you expect it to be.
c
      if (ifcore .ne. 1) then
         do 400 j = 2, nr
            write(3,9900) r(j), cdu(j), cdd(j), zero
 400     continue
      else
         do 410 j = 2, nr
            write(3,9900) r(j), cdu(j), cdd(j), cdc(j)
 410     continue
      endif
c
 9900 format(1x,f15.10,3x,3f15.8)
c
      close(unit=3)
c
      return
c
      end

      double precision function cutoff_function(r)
      implicit none
c
c     Generalized cutoff function
c
      double precision r
c
c     Standard cutoff
c
      cutoff_function = exp(-5.d0*r)

      end










C
c $Id: splift.f,v 1.2 1997/05/22 17:32:32 wdpgaara Exp $
c
c $Log: splift.f,v $
c Revision 1.2  1997/05/22 17:32:32  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:55  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine splift(x,y,yp,ypp,n,w,ierr,isx,a1,b1,an,bn)
C
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES
C  *                   A PRIME CONTRACTOR TO THE
C  *                UNITED STATES DEPARTMENT OF ENERGY
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
C  * OWNED RIGHTS.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
C  * PART IS SAND77-1441.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     WRITTEN BY RONDALL E. JONES
C
C     ABSTRACT
C         SPLIFT FITS AN INTERPOLATING CUBIC SPLINE TO THE N DATA POINT
C         GIVEN IN X AND Y AND RETURNS THE FIRST AND SECOND DERIVATIVES
C         IN YP AND YPP.  THE RESULTING SPLINE (DEFINED BY X, Y, AND
C         YPP) AND ITS FIRST AND SECOND DERIVATIVES MAY THEN BE
C         EVALUATED USING SPLINT.  THE SPLINE MAY BE INTEGRATED USING
C         SPLIQ.  FOR A SMOOTHING SPLINE FIT SEE SUBROUTINE SMOO.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.   X(N), Y(N), YP(N), YPP(N), W(3N)
C
C       --INPUT--
C
C         X    - ARRAY OF ABSCISSAS OF DATA (IN INCREASING ORDER)
C         Y    - ARRAY OF ORDINATES OF DATA
C         N    - THE NUMBER OF DATA POINTS.  THE ARRAYS X, Y, YP, AND
C                YPP MUST BE DIMENSIONED AT LEAST N.  (N .GE. 4)
C         ISX  - MUST BE ZERO ON THE INITIAL CALL TO SPLIFT.
C                IF A SPLINE IS TO BE FITTED TO A SECOND SET OF DATA
C                THAT HAS THE SAME SET OF ABSCISSAS AS A PREVIOUS SET,
C                AND IF THE CONTENTS OF W HAVE NOT BEEN CHANGED SINCE
C                THAT PREVIOUS FIT WAS COMPUTED, THEN ISX MAY BE
C                SET TO ONE FOR FASTER EXECUTION.
C         A1,B1,AN,BN - SPECIFY THE END CONDITIONS FOR THE SPLINE WHICH
C                ARE EXPRESSED AS CONSTRAINTS ON THE SECOND DERIVATIVE
C                OF THE SPLINE AT THE END POINTS (SEE YPP).
C                THE END CONDITION CONSTRAINTS ARE
C                        YPP(1) = A1*YPP(2) + B1
C                AND
C                        YPP(N) = AN*YPP(N-1) + BN
C                WHERE
C                        ABS(A1).LT. 1.0  AND  ABS(AN).LT. 1.0.
C
C                THE SMOOTHEST SPLINE (I.E., LEAST INTEGRAL OF SQUARE
C                OF SECOND DERIVATIVE) IS OBTAINED BY A1=B1=AN=BN=0.
C                IN THIS CASE THERE IS AN INFLECTION AT X(1) AND X(N).
C                IF THE DATA IS TO BE EXTRAPOLATED (SAY, BY USING SPLIN
C                TO EVALUATE THE SPLINE OUTSIDE THE RANGE X(1) TO X(N))
C                THEN TAKING A1=AN=0.5 AND B1=BN=0 MAY YIELD BETTER
C                RESULTS.  IN THIS CASE THERE IS AN INFLECTION
C                AT X(1) - (X(2)-X(1)) AND AT X(N) + (X(N)-X(N-1)).
C                IN THE MORE GENERAL CASE OF A1=AN=A  AND B1=BN=0,
C                THERE IS AN INFLECTION AT X(1) - (X(2)-X(1))*A/(1.0-A)
C                AND AT X(N) + (X(N)-X(N-1))*A/(1.0-A).
C
C                A SPLINE THAT HAS A GIVEN FIRST DERIVATIVE YP1 AT X(1)
C                AND YPN AT Y(N) MAY BE DEFINED BY USING THE
C                FOLLOWING CONDITIONS.
C
C                A1=-0.5
C
C                B1= 3.0*((Y(2)-Y(1))/(X(2)-X(1))-YP1)/(X(2)-X(1))
C
C                AN=-0.5
C
C                BN=-3.0*((Y(N)-Y(N-1))/(X(N)-X(N-1))-YPN)/(X(N)-X(N-1)
C
C       --OUTPUT--
C
C         YP   - ARRAY OF FIRST DERIVATIVES OF SPLINE (AT THE X(I))
C         YPP  - ARRAY OF SECOND DERIVATIVES OF SPLINE (AT THE X(I))
C         IERR - A STATUS CODE
C              --NORMAL CODE
C                 1 MEANS THAT THE REQUESTED SPLINE WAS COMPUTED.
C              --ABNORMAL CODES
C                 2 MEANS THAT N, THE NUMBER OF POINTS, WAS .LT. 4.
C                 3 MEANS THE ABSCISSAS WERE NOT STRICTLY INCREASING.
C
C       --WORK--
C
C         W    - ARRAY OF WORKING STORAGE DIMENSIONED AT LEAST 3N.
C
C     .. Parameters ..
      double precision four
      parameter (four=4.D0)
C     ..
C     .. Scalar Arguments ..
      double precision a1, an, b1, bn
      integer ierr, isx, n
C     ..
C     .. Array Arguments ..
      double precision w(n,3), x(n), y(n), yp(n), ypp(n)
C     ..
C     .. Local Scalars ..
      double precision dnew, dold
      integer i, j, nm1, nm2
C     ..
      if (n .lt. 4) then
         ierr = 2
c
         return
c
      end if
      nm1 = n - 1
      nm2 = n - 2
      if (isx .gt. 0) go to 40
      do 10 i = 2, n
         if (x(i)-x(i-1) .le. 0) then
            ierr = 3
c
            return
c
         end if
   10 continue
C
C     DEFINE THE TRIDIAGONAL MATRIX
C
      w(1,3) = x(2) - x(1)
      do 20 i = 2, nm1
         w(i,2) = w(i-1,3)
         w(i,3) = x(i+1) - x(i)
         w(i,1) = 2*(w(i,2)+w(i,3))
   20 continue
      w(1,1) = four
      w(1,3) = -4*a1
      w(n,1) = four
      w(n,2) = -4*an
C
C     L U DECOMPOSITION
C
      do 30 i = 2, n
         w(i-1,3) = w(i-1,3)/w(i-1,1)
         w(i,1) = w(i,1) - w(i,2)*w(i-1,3)
   30 continue
C
C     DEFINE *CONSTANT* VECTOR
C
   40 continue
      ypp(1) = 4*b1
      dold = (y(2)-y(1))/w(2,2)
      do 50 i = 2, nm2
         dnew = (y(i+1)-y(i))/w(i+1,2)
         ypp(i) = 6*(dnew-dold)
         yp(i) = dold
         dold = dnew
   50 continue
      dnew = (y(n)-y(n-1))/(x(n)-x(n-1))
      ypp(nm1) = 6*(dnew-dold)
      ypp(n) = 4*bn
      yp(nm1) = dold
      yp(n) = dnew
C
C     FORWARD SUBSTITUTION
C
      ypp(1) = ypp(1)/w(1,1)
      do 60 i = 2, n
         ypp(i) = (ypp(i)-w(i,2)*ypp(i-1))/w(i,1)
   60 continue
C
C     BACKWARD SUBSTITUTION
C
      do 70 j = 1, nm1
         i = n - j
         ypp(i) = ypp(i) - w(i,3)*ypp(i+1)
   70 continue
C
C     COMPUTE FIRST DERIVATIVES
C
      yp(1) = (y(2)-y(1))/(x(2)-x(1)) - (x(2)-x(1))*(2*ypp(1)+ypp(2))/6
      do 80 i = 2, nm1
         yp(i) = yp(i) + w(i,2)*(ypp(i-1)+2*ypp(i))/6
   80 continue
      yp(n) = yp(n) + (x(n)-x(nm1))*(ypp(nm1)+2*ypp(n))/6
C
      ierr = 1
c
      return
c
      end
C
c $Id: spliq.f,v 1.2 1997/05/22 17:32:32 wdpgaara Exp $
c
c $Log: spliq.f,v $
c Revision 1.2  1997/05/22 17:32:32  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:55  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine spliq(x,y,yp,ypp,n,xlo,xup,nup,ans,ierr)
C
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES
C  *                   A PRIME CONTRACTOR TO THE
C  *                UNITED STATES DEPARTMENT OF ENERGY
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
C  * OWNED RIGHTS.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
C  * PART IS SAND77-1441.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     THIS ROUTINE WAS WRITTEN BY M. K. GORDON
C
C     ABSTRACT
C
C     SUBROUTINE SPLIQ INTEGRATES A CUBIC SPLINE (GENERATED BY
C     SPLIFT, SMOO, ETC.) ON THE INTERVALS (XLO,XUP(I)), WHERE XUP
C     IS A SEQUENCE OF UPPER LIMITS ON THE INTERVALS OF INTEGRATION.
C     THE ONLY RESTRICTIONS ON XLO AND XUP(*) ARE
C                XLO .LT. XUP(1),
C                XUP(I) .LE. XUP(I+1)   FOR EACH I .
C     ENDPOINTS BEYOND THE SPAN OF ABSCISSAS ARE ALLOWED.
C     THE SPLINE OVER THE INTERVAL (X(I),X(I+1)) IS REGARDED
C     AS A CUBIC POLYNOMIAL EXPANDED ABOUT X(I) AND IS INTEGRATED
C     ANALYTICALLY.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.  X(N), Y(N), YP(N), YPP(N), XUP(NUP), ANS(NUP)
C
C      --INPUT--
C
C        X    - ARRAY OF ABSCISSAS (IN INCREASING ORDER) THAT DEFINE TH
C               SPLINE.  USUALLY X IS THE SAME AS X IN SPLIFT OR SMOO.
C        Y    - ARRAY OF ORDINATES THAT DEFINE THE SPLINE.  USUALLY Y I
C               THE SAME AS Y IN SPLIFT OR AS R IN SMOO.
C        YP   - ARRAY OF FIRST DERIVATIVES OF THE SPLINE AT ABSCISSAS.
C               USUALLY YP IS THE SAME AS YP IN SPLIFT OR R1 IN SMOO.
C        YPP  - ARRAY OF SECOND DERIVATIVES THAT DEFINE THE SPLINE.
C               USUALLY YPP IS THE SAME AS YPP IN SPLIFT OR R2 IN SMOO.
C        N    - THE NUMBER OF DATA POINTS THAT DEFINE THE SPLINE.
C        XLO  - LEFT ENDPOINT OF INTEGRATION INTERVALS.
C        XUP  - RIGHT ENDPOINT OR ARRAY OF RIGHT ENDPOINTS OF
C               INTEGRATION INTERVALS IN ASCENDING ORDER.
C        NUP  - THE NUMBER OF RIGHT ENDPOINTS.  IF NUP IS GREATER THAN
C               1, THEN XUP AND ANS MUST BE DIMENSIONED AT LEAST NUP.
C
C      --OUTPUT--
C
C        ANS -- ARRAY OF INTEGRAL VALUES, THAT IS,
C               ANS(I) = INTEGRAL FROM XLO TO XUP(I)
C        IERR -- ERROR STATUS
C                = 1 INTEGRATION SUCCESSFUL
C                = 2 IMPROPER INPUT - N.LT.4 OR NUP.LT.1
C                = 3 IMPROPER INPUT - ABSCISSAS NOT IN
C                        STRICTLY ASCENDING ORDER
C                = 4 IMPROPER INPUT - RIGHT ENDPOINTS XUP NOT
C                        IN ASCENDING ORDER
C                = 5 IMPROPER INPUT - XLO.GT.XUP(1)
C                = 6 INTEGRATION SUCCESSFUL BUT AT LEAST ONE ENDPOINT
C                        NOT WITHIN SPAN OF ABSCISSAS
C              ** NOTE.  ERRCHK PROCESSES DIAGNOSTICS FOR CODES 2,3,4,5
C
C   CHECK FOR IMPROPER INPUT
C
C     .. Scalar Arguments ..
      double precision xlo
      integer ierr, n, nup
C     ..
C     .. Array Arguments ..
      double precision ans(nup), x(n), xup(nup), y(n), yp(n), ypp(n)
C     ..
C     .. Local Scalars ..
      double precision hdiff, hi, hi2, hi3, hlo, hlo2, hsum, hup, hup2,
     &                 hup3, hup4, psum0, psum1, psum2, psum3, sum,
     &                 sum0, sum1, sum2, sum3
      integer i, j, m, nm1, nm2
C     ..
      ierr = 2
      if (n .lt. 4 .or. nup .lt. 1) then
c
         return
c
      end if
      nm1 = n - 1
      nm2 = n - 2
      ierr = 3
      do 10 i = 1, nm1
         if (x(i) .ge. x(i+1)) then
c
            return
c
         end if
   10 continue
      if (nup .ne. 1) then
         ierr = 4
         do 20 i = 2, nup
            if (xup(i-1) .gt. xup(i)) then
c
               return
c
            end if
   20    continue
      end if
      ierr = 5
      if (xlo .gt. xup(1)) then
c
         return
c
      end if
      ierr = 1
      if (xlo .lt. x(1) .or. xup(nup) .gt. x(n)) ierr = 6
C
C   LOCATE XLO IN INTERVAL (X(I),X(I+1))
C
      do 30 i = 1, nm2
         if (xlo .lt. x(i+1)) go to 40
   30 continue
      i = nm1
   40 continue
      hlo = xlo - x(i)
      hlo2 = hlo*hlo
      hi = x(i+1) - x(i)
      hi2 = hi*hi
      do 50 j = 1, nup
         if (xup(j) .gt. x(i+1) .and. xlo .lt. x(nm1)) go to 60
C
C   COMPUTE SPECIAL CASES OF XUP IN INTERVAL WITH XLO
C
         hup = xup(j) - x(i)
         hsum = hup + hlo
         hdiff = hup - hlo
         hup2 = hup*hup
         sum = (ypp(i+1)-ypp(i))*hsum*hdiff*(hup2+hlo2)/(24*hi)
         sum = sum + ypp(i)*hdiff*(hup2+hlo*hup+hlo2)/6
         sum = sum + yp(i)*hdiff*hsum/2
         sum = sum + y(i)*hdiff
         ans(j) = sum
   50 continue
c
      return
C
C   COMPUTE INTEGRAL BETWEEN XLO AND X(I+1) AS FOUR TERMS IN TAYLOR
C   POLYNOMIAL AND ADVANCE I TO I+1
C
   60 continue
      hdiff = hi - hlo
      hsum = hi + hlo
      sum0 = y(i)*hdiff
      sum1 = yp(i)*hdiff*hsum
      sum2 = ypp(i)*hdiff*(hi2+hi*hlo+hlo2)
      sum3 = (ypp(i+1)-ypp(i))*hdiff*hsum*(hi2+hlo2)/hi
      i = i + 1
C
C   LOCATE EACH XUP(M) IN INTERVAL (X(I),X(I+1))
C
      do 90 m = j, nup
   70    continue
         if (xup(m) .lt. x(i+1) .or. i .eq. nm1) go to 80
C
C   AUGMENT INTEGRAL BETWEEN ABSCISSAS TO INCLUDE INTERVAL
C   (X(I),X(I+1)) AND ADVANCE I TO I+1
C
         hi = x(i+1) - x(i)
         hi2 = hi*hi
         hi3 = hi2*hi
         sum0 = sum0 + y(i)*hi
         sum1 = sum1 + yp(i)*hi2
         sum2 = sum2 + ypp(i)*hi3
         sum3 = sum3 + (ypp(i+1)-ypp(i))*hi3
         i = i + 1
c
         go to 70
C
C   INTEGRAL BETWEEN X(I) AND XUP(M) IS ZERO
C
   80    continue
         if (xup(m) .ne. x(i)) then
C
C   COMPUTE INTEGRAL BETWEEN X(I) AND XUP(M) AND EVALUATE
C   TAYLOR POLYNOMIAL IN REVERSE ORDER
C
            hup = xup(m) - x(i)
            hup2 = hup*hup
            hup3 = hup2*hup
            hup4 = hup3*hup
            hi = x(i+1) - x(i)
            psum0 = y(i)*hup
            psum1 = yp(i)*hup2
            psum2 = ypp(i)*hup3
            psum3 = (ypp(i+1)-ypp(i))*hup4/hi
            sum = (sum3+psum3)/24 + (sum2+psum2)/6
            sum = sum + (sum1+psum1)/2
            sum = sum + (sum0+psum0)
         else
            sum = ((sum3/24+sum2/6)+sum1/2) + sum0
         end if
         ans(m) = sum
   90 continue
c
      return
c
      end
c
      subroutine velect(iter,iconv,id,zelec)
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'charge.h'
      include 'elecpot.h'
      include 'energy.h'
      include 'compat.h'
c
c    velect generates the electronic output potential from
c    the electron charge density.  The ionic part is
c    added in dsolv1/dsolv2.
c
c    NOTE:  Velect can be called with zelec=ZEL (all electron)
c                             or with zelec=ZVAL (valence only)
c
C     .. Parameters ..
c
      double precision tiny_charge
      parameter (tiny_charge=1.d-12)
c
      double precision zero, one, pnn
      parameter (zero=0.D0,one=1.D0,pnn=.99D0)
c
C     ..
C     .. Scalar Arguments ..
      double precision zelec
      integer iconv, iter
      character id*1
C     ..
C     .. Arrays in Common ..
c
      double precision s1(nrmax), s2(nrmax), w(3*nrmax), y(nrmax),
     &                 yp(nrmax), ypp(nrmax)
C     ..
C     .. Local Scalars ..
      double precision a1, an, b1, bn, ehart, xlo, xnorm,
     &                 ec, exc, vc, vxc, pi, dx, dc, ex, xccor
      integer i, ierr, isx, ll, nrm, relflag, nspin
C     ..
C     .. Local Arrays ..
      double precision dens(nrmax,2), vxcarr(nrmax,2)

C     .. External Subroutines ..
      external ext, splift, spliq, excorr
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs
C     ..
      logical leqi
      external leqi
c
C     .. Common blocks ..
      common  s1, s2, w, y, yp, ypp
C     ..
c
      pi = 4.d0 * atan(1.d0)
      ehart = zero
c
c      fit cd/r by splines
c
      y(1) = zero
      do 10 i = 2, nr
         y(i) = (cdd(i)+cdu(i))/r(i)
   10 continue
      if (ifcore .eq. 2) then
         do 20 i = 2, nr
            y(i) = y(i) + cdc(i)/r(i)
   20    continue
      end if
      isx = 0
      a1 = zero
      an = zero
      b1 = zero
      bn = zero
      nrm = nr
      call splift(r,y,yp,ypp,nrm,w,ierr,isx,a1,b1,an,bn)
      if (ierr .ne. 1) then
         write(6,9000) ierr
         call ext(420+ierr)
      end if
 9000 format(1x,'****** Error in splift ierr =',i2)
c
c      compute the integrals of cd/r and cd from
c      r(1)=0 to r(i)
c
      xlo = zero
      call spliq(r,y,yp,ypp,nrm,xlo,r,nrm,s2,ierr)
      if (ierr .ne. 1) then
         write(6,9010) ierr
         call ext(440+ierr)
      end if
 9010 format(1x,'****** Error in spliq ierr =',i2)
      do 30 i = 1, nr
         ypp(i) = r(i)*ypp(i) + 2*yp(i)
         yp(i) = r(i)*yp(i) + y(i)
         y(i) = r(i)*y(i)
   30 continue
      call spliq(r,y,yp,ypp,nrm,xlo,r,nrm,s1,ierr)
      if (ierr .ne. 1) then
         write(6,9020) ierr
         call ext(460+ierr)
      end if
 9020 format(1x,'****** Error in spliq ierr =',i2)
c
c      check normalization
c
      xnorm = zero
      if (zelec .ne. zero) xnorm = zelec/s1(nr)
      if (iter .gt. 3 .and. abs(zelec-s1(nr)) .gt. 0.01D0) then
         if (zelec .lt. s1(nr)+1.0D0) then
            write(6,9030) iter, xnorm
 9030       format(/' warning *** charge density rescaled in',' velect',
     &            /' iteration number',i4,3x,'scaling factor =',f6.3,/)
         else
            xnorm = pnn*xnorm
            write(6,9040) iter, xnorm
 9040       format(/' warning *** charge density partially rescaled in',
     &            ' velect',/' iteration number',i4,3x,
     &            'scaling factor =',f6.3,/)
         end if
      end if
c
c      compute new hartree potential
c      renormalize the charge density
c
      do 40 i = 2, nr
         vod(i) = 2*xnorm*(s1(i)/r(i)+s2(nr)-s2(i))
         vou(i) = vod(i)
         cdd(i) = xnorm*cdd(i)
         cdu(i) = xnorm*cdu(i)
   40 continue
c
c      compute hartree contribution to total energy
c
      if (iconv .eq. 1) then
         ehart = zero
         ll = 4
         do 50 i = 2, nr
            ehart = ehart + ll*(cdd(i)+cdu(i))*vod(i)*rab(i)
            ll = 6 - ll
   50    continue
c
c        Divide by two to account for overcounting of the interactions
c        and by three to finish the Simpson rule calculation...
c
         ehart = 0.5d0*(ehart/3)
      end if
c
c     Add the exchange and correlation potential and calculate
c     the total energy contributions.
c
      if ( leqi(icorr,'gl') .or. leqi(icorr,'hl') .or.
     $     leqi(icorr,'wi') .or. leqi(icorr,'bh') ) then
         use_excorr = .true.
      endif

      if (use_excorr) then

         call excorr(id,cdd,cdu,cdc,vod,vou,vxc,vc,exc,ec)
         etot(4) = ehart
         etot(5) = vxc
         etot(6) = (3*vc-4*ec)
         etot(7) = exc
c
         return

      else
c
c        New XC scheme based on code from Jose Soler and Carlos Balbas
c
c        Compute dens(i,nspin) = density up, density down
c
         do i=2,nr
            if (ispp .eq. 's') then
               dens(i,1) = cdu(i)/(4.d0*pi*r(i)**2)
               dens(i,2) = cdd(i)/(4.d0*pi*r(i)**2)
            else
               dens(i,1) = 0.5d0*(cdu(i) + cdd(i))/(4.d0*pi*r(i)**2)
               dens(i,2) = dens(i,1)
            endif
            if (ifcore .ge. 1) then
               dens(i,1) = dens(i,1) + 0.5d0 * cdc(i)/(4.d0*pi*r(i)**2)
               dens(i,2) = dens(i,2) + 0.5d0 * cdc(i)/(4.d0*pi*r(i)**2)
            endif
         enddo

c
c        Extrapolate the density at r=0  
c
         dens(1,1) = dens(2,1) - (dens(3,1)-dens(2,1))*r(2)/(r(3)-r(2))
         dens(1,2) = dens(2,2) - (dens(3,2)-dens(2,2))*r(2)/(r(3)-r(2))
         if (dens(1,1) .lt. 0.d0) dens(1,1) = 0.d0
         if (dens(1,2) .lt. 0.d0) dens(1,2) = 0.d0
c
c        Define 'relflag' and 'nspin' for the interface ATOMXC
c
         if (ispp .eq. 'r') relflag = 1
         if (ispp .ne. 'r') relflag = 0
         nspin = 2
c
         r(1) = 0.0d0

         is_gga = (leqi(icorr,'pb')        ! PBE
     $              .or. leqi(icorr,'bl')  ! BLYP
     $              .or. leqi(icorr,'rp')  ! RPBE
     $              .or. leqi(icorr,'rv')) ! revPBE

         if (icorr .eq. 'ca') then
            call atomxc('LDA','ca',relflag,nr,nrmax,r,nspin,dens,
     .           ex,ec,dx,dc,vxcarr)      
         elseif(icorr .eq. 'pw') then
            call atomxc('LDA','pw92',relflag,nr,nrmax,r,nspin,dens,
     .           ex,ec,dx,dc,vxcarr)
         elseif(icorr .eq. 'pb') then
            call atomxc('GGA','pbe',relflag,nr,nrmax,r,nspin,dens,   
     .           ex,ec,dx,dc,vxcarr)
         elseif(icorr .eq. 'rp') then
            call atomxc('GGA','rpbe',relflag,nr,nrmax,r,nspin,dens,   
     .           ex,ec,dx,dc,vxcarr)
         elseif(icorr .eq. 'rv') then
            call atomxc('GGA','revpbe',relflag,nr,nrmax,r,nspin,dens,   
     .           ex,ec,dx,dc,vxcarr)
         elseif(icorr .eq. 'bl') then
            call atomxc('GGA','lyp',relflag,nr,nrmax,r,nspin,dens,
     .           ex,ec,dx,dc,vxcarr)
         else
            stop 'XC'
         endif

c
c        Add vxc to total potential and energies
c   
         do i=2,nr
            vou(i) = vou(i) + vxcarr(i,1)
            vod(i) = vod(i) + vxcarr(i,2)
         enddo

clcb
         xccor=0.0d0
         ll = 4
         do i=2,nr
            xccor = xccor + ll * rab(i)*
     .           (vxcarr(i,1)*cdu(i) + vxcarr(i,2)*cdd(i))
            ll = 6 - ll
         enddo

         etot(4) = ehart
         etot(5) = xccor / 3
         etot(6) = xccor - 4*(ex + ec)
         etot(7) = ex + ec 

clcb  
c
c        Obtain total potential at r = 0
c
         vod(1) = vod(2) - (vod(3)-vod(2))*r(2)/(r(3)-r(2))
         vou(1) = vou(2) - (vou(3)-vou(2))*r(2)/(r(3)-r(2))

c     *** lcb-jms modification end ********************
         return

      endif
c
c
      end









C
c $Id: vionic.f,v 1.4 2002/07/08 18:08:26 wdpgaara Exp $
c
c $Log: vionic.f,v $
c Revision 1.4  2002/07/08 18:08:26  wdpgaara
c Re-implemented the "valence charge modification" feature.
c New routine: change_valence.
c Superseded vionic kludge. Watch out for compiler complaints.
c
c Revision 1.3  2002/07/05 18:22:39  wdpgaara
c Fix format of grid parameters
c
c Revision 1.2  1997/05/22 17:32:36  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:55  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine Vionic
c
c  Vionic sets up the ionic potential.
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'charge.h'
c
c  njtj ***  major modifications  ***
c    If a potential does not exist, it is approximated
c    by an existing potential.
c    A nonspin or spin-polarized pseudo test, uses the
c    down(nonspin generation), weighted average(spin-
c    polarized), or averaged(relativistic) potentials.
c    A relativistic pseudo test, must use relativistic
c    generated potentials.  The Schroedinger equation is
c    used to integrate a relativistic pseudo test,
c    not the Dirac equation.
c  njtj  ***  major modifications  ***
c
C     .. Parameters ..
      double precision zero
      parameter (zero=0.D0)
C     ..
C     .. Local Scalars ..
      double precision vdiff, vsum, zion
      integer i, j, loi, npotd, npotu, nrm
      character icorrt*2, namet*2
C     ..
C     .. Local Arrays ..
      integer npd(5), npu(5)
      character ray(6)*10, title(7)*10
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs
C     ..
c
c  2*znuc part
c
      ifcore = 0
      if (job .lt. 4) then
         do 20 i = 1, lmax
            do 10 j = 1, nr
               viod(i,j) = -2*znuc
               viou(i,j) = -2*znuc
   10       continue
   20    continue
      else
c
c  read pseudopotentials from tape1
c
         open(unit=1,file='VPSIN',status='old',form='unformatted')
         rewind 1
         read(1) namet, icorrt, irel, nicore, (ray(i),i=1,6),
     &     (title(i),i=1,7), npotd, npotu, nrm, a, b, zion
c
         if (nicore .eq. 'fcec' .or. nicore .eq. 'pcec') ifcore = 1
         if (nicore .eq. 'fche' .or. nicore .eq. 'pche') ifcore = 2
         nr = nrm + 1
         read(1) (r(i),i=2,nr)
         r(1) = zero
c
c   down potentials (or average relativistic potentials)
c
c njtj  ***  major start  ***
c   if a potential does not exist, it is replaced by the
c   next existing lower angular momentum potential or
c   the next existing higher if no lower exists.
c
         do 30 i = 1, lmax
            npd(i) = 0
   30    continue
         do 40 i = 1, npotd
            read(1) loi, (viod(loi+1,j),j=2,nr)
            viod(loi+1,1) = zero
            npd(loi+1) = 1
   40    continue
         if (npd(1) .eq. 0) then
            do 60 i = 2, lmax
               if (npd(i) .gt. 0) then
                  do 50 j = 1, nr
                     viod(1,j) = viod(i,j)
   50             continue
c
                  go to 70
c
               end if
   60       continue
         end if
   70    continue
         do 90 i = 2, lmax
            if (npd(i) .eq. 0) then
               do 80 j = 1, nr
                  viod(i,j) = viod(i-1,j)
   80          continue
            end if
   90    continue
c
c   up potentials (or spin orbit potentials)
c
         if (npotu .le. 0) go to 170
         do 100 i = 1, lmax
            npu(i) = 0
  100    continue
         do 110 i = 1, npotu
            read(1) loi, (viou(loi+1,j),j=2,nr)
            viou(loi+1,1) = zero
            npu(loi+1) = 1
  110    continue
         if (npu(1) .eq. 0) then
            do 130 i = 2, lmax
               if (npu(i) .gt. 0) then
                  do 120 j = 1, nr
                     viou(1,j) = viou(i,j)
  120             continue
c
                  go to 140
c
               end if
  130       continue
         end if
  140    continue
         do 160 i = 2, lmax
            if (npu(i) .eq. 0) then
               do 150 j = 1, nr
                  viou(i,j) = viou(i-1,j)
  150          continue
            end if
  160    continue
c
c  njtj  ***  major end  ***
c
c
c  core and valence charges
c
  170    continue
         read(1) (cdc(i),i=2,nr)
         cdc(1) = zero
c
c  replace valence charge on tape(valence charge modify)
c  This is a horrible kludge, and it probably has never
c  been used after the Froyen days...
c  It might also violate the Fortran Standard.
c
c  Note the instantaneous job number change needed
c  Leave it as is, but use routine change_valence
c  for a cleaner operation.
c
         if (job .eq. 6) then
            write(1) (cdd(i)+cdu(i),i=2,nr)
            close(1)
c
            return
c
         end if
         read(1) (cdd(i),i=2,nr)

         close(1)

         cdd(1) = zero
c
c  njtj  ***   major start  ***
c   distribute charge as up and down charge
c   generate radial intergration grid
c   set up potentials equal to down potentials for
c   spin-polarized pseudo test of nonspin and relativistic
c   generated potentails.  Construct spin-orbit potentials
c   from relativistic sum and difference potentials and
c   change ispp='r' to ispp=' '.
c
         do 180 i = 1, nr
            rab(i) = (r(i)+a)*b
            cdd(i) = cdd(i)/2
            cdu(i) = cdd(i)
  180    continue
         if (ispp .eq. 's' .and. irel .ne. 'isp') then
            do 200 i = 1, lmax
               do 190 j = 1, nr
                  viou(i,j) = viod(i,j)
  190          continue
  200       continue
         end if
         if (ispp .eq. 'r') then
            ispp = ' '
            if (irel .ne. 'rel') then
               write(6,9000)
 9000          format(//'Pseudopotential is not relativistic!!!!',
     &               /' setting up potentials equal to down!!!',//)
               do 220 i = 1, lmax
                  do 210 j = 1, nr
                     viou(i,j) = viod(i,j)
  210             continue
  220          continue
            else
               do 230 j = 1, nr
                  viou(1,j) = viod(1,j)
  230          continue
               do 250 i = 2, lmax
                  do 240 j = 1, nr
                     vsum = viod(i,j)
                     vdiff = viou(i,j)
                     viod(i,j) = vsum - i*vdiff/2
                     viou(i,j) = vsum + (i-1)*vdiff/2
  240             continue
  250          continue
            end if
         end if
c
c   njtj  ***  major end   ***
c
c
c   printout
c
         write(6,9010) namet, icorrt, irel, nicore,
     &     (ray(i),i=1,6), (title(i),i=1,7)
 9010    format(//1x,a2,2x,a2,2x,a3,2x,a4,
     &         '  pseudopotential read from tape',/1x,2a10,5x,4a10,/1x,
     &         7a10,//)
         if (nameat .ne. namet) write(6,9020) nameat, namet
 9020    format(' input element ',a2,' not equal to element on tape ',
     &         a2,//)
         if (icorr .ne. icorrt) write(6,9030) icorr, icorrt
 9030    format(' input correlation ',a2,
     &         ' not equal to correlation from tape ',a2,//)
         write(6,9040) r(2), nr, r(nr)
 9040    format(' radial grid parameters',//' r(1) = .0 , r(2) =',d8.2,
     &         ' , ... , r(',i4,') =',f6.2,//)
      end if
c
c   add potential from shell charge
c
      if (abs(zsh) .gt. 0.D-5) then
         do 270 i = 1, lmax
            do 260 j = 1, nr
               if (r(j) .ge. rsh) then
                  viod(i,j) = viod(i,j) - 2*zsh
                  viou(i,j) = viou(i,j) - 2*zsh
               else
                  viod(i,j) = viod(i,j) - 2*zsh*r(j)/rsh
                  viou(i,j) = viou(i,j) - 2*zsh*r(j)/rsh
               end if
  260       continue
  270    continue
      end if
c
      return
c
      end
C
c $Id: wtrans.f,v 1.3 1997/05/22 18:05:43 wdpgaara Exp $
c
      subroutine wtrans(vd,r,nr,i,ist)
c
      implicit none
c
      include 'plot.h'
c
c **********************************************************
c *  The wave function is fitted
c *  with a second degree polynomial which is muliplied
c *  with the appropriate functions and then integrated
c *  by parts in taking the fourier transform. 
c **********************************************************
c
c  The potential times r is fitted to the polynominal
c  a + bx + cx^2 at every other point.
c
C     .. Parameters ..
      double precision zero, one
      parameter (zero=0.D0,one=1.D0)
C     ..
C     .. Scalar Arguments ..
      integer i, ist, nr
C     ..
C     .. Array Arguments ..
      double precision r(nr), vd(nr)
C     ..
C     .. Local Scalars ..
      double precision d1, d2, d3, q, q2, r0, rm, rp, v0, vm, vp
      integer j, k
      character*7 filename
C     ..
C     .. Local Arrays ..
      double precision vql(100)
C     ..
C     .. Intrinsic Functions ..
      intrinsic cos, sin
C     ..
      rm = zero
      vm = zero
      do 10 k = 2, nr-1, 2
         r0 = r(k)
         v0 = vd(k)
         rp = r(k+1)
         vp = vd(k+1)
         d1 = 1/((rp-rm)*(r0-rm))
         d2 = 1/((rp-r0)*(rm-r0))
         d3 = 1/((r0-rp)*(rm-rp))
         a(k) = vm*d1 + v0*d2 + vp*d3
         b(k) = -vm*(r0+rp)*d1 - v0*(rm+rp)*d2 - vp*(rm+r0)*d3
         c(k) = vm*r0*rp*d1 + v0*rm*rp*d2 + vp*rm*r0*d3
         rm = rp
         vm = vp
   10 continue
c
c  Find the fourier transform-vql.
c
      do 30 j = 1, 54
         q = one/4*j
         q2 = q*q
         vql(j) = zero
         rm = zero
         do 20 k = 2, nr - 1, 2
            rp = r(k+1)
            vql(j) = vql(j) + (2*a(k)*rp+b(k))/q*sin(q*rp) -
     &               ((a(k)*rp+b(k))*rp+c(k)-2*a(k)/q2)*cos(q*rp) -
     &               (2*a(k)*rm+b(k))/q*sin(q*rm) +
     &               ((a(k)*rm+b(k))*rm+c(k)-2*a(k)/q2)*cos(q*rm)
            rm = rp
   20    continue
         vql(j) = vql(j)/q2
   30 continue
c
c  Print out the transform vql(q) to the current plot.dat
c  file (unit=3) for latter plotting.
c
      write(filename,9900) i
 9900 format('PSWFNQ',i1)
      open(unit=3,file=filename,form='formatted',status='unknown')
c
      do 40 j = 1, 48
         write(3,9000) one/4*j, ist*vql(j)
   40 continue
 9000 format(1x,f7.4,3x,f10.6)
c
      close(unit=3)
cag
c      write(3,9010) i
c 9010 format(1x,'marker fw',i1)
c
      return
c
      end
c
      double precision function v0pp(gamma)
c
      implicit none
c
      include 'radial.h'
      include 'tm2_blk.h'
c
      double precision gamma
c
C     .. Parameters ..
      double precision zero, pfive, one, errmin
      parameter (zero=0.D0,pfive=0.5D0,one=1.D0,errmin=1.D-12)
C     ..
C     .. Local Scalars ..
      double precision bj1, bj2, bj2a, bj3, bj3a, bj4, bj5, cdps,
     &                 ddelta, fdnew, fdold, polyr, r2, rc10, rc11,
     &                 rc12, rc9, rp, rcond
      integer j, k, ll
C     ..
C     .. Local Arrays ..
      double precision aj(5,5), bj(5), work(5)
      integer indx(5)
C     ..
C     .. External Subroutines ..
      external ext
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, exp, log
C     ..

      fdold = 0.d0

      rc9 = rc8*rc1
      rc10 = rc8*rc2
      rc11 = rc8*rc3
      rc12 = rc8*rc4
c
      delta = zero
c
      bj(1) = log(arc/rc1**lp) - gamma*rc2
      bj1 = bj(1)
      bj(2) = brc - lp/rc1 - 2*gamma*rc1
      bj2a = bj(2) + 2*gamma*rc1
      bj2 = bj(2)
      bj(3) = vrc - eigv - 2*lp/rc1*bj2a - bj2a**2 - 2*gamma
      bj3 = bj(3)
      bj3a = bj(3) + 2*gamma
      bj(4) = vap + 2*lp/rc2*bj2a - 2*lp/rc1*bj3a - 2*bj2a*bj3a
      bj4 = bj(4)
      bj(5) = vapp - 4*lp/rc3*bj2a + 4*lp/rc2*bj3a - 2*lp/rc1*bj4 -
     &        2*bj3a**2 - 2*bj2a*bj4
      bj5 = bj(5)
c
      aj(1,1) = rc4
      aj(1,2) = rc6
      aj(1,3) = rc8
      aj(1,4) = rc10
      aj(1,5) = rc12
      aj(2,1) = 4*rc3
      aj(2,2) = 6*rc5
      aj(2,3) = 8*rc7
      aj(2,4) = 10*rc9
      aj(2,5) = 12*rc11
      aj(3,1) = 12*rc2
      aj(3,2) = 30*rc4
      aj(3,3) = 56*rc6
      aj(3,4) = 90*rc8
      aj(3,5) = 132*rc10
      aj(4,1) = 24*rc1
      aj(4,2) = 120*rc3
      aj(4,3) = 336*rc5
      aj(4,4) = 720*rc7
      aj(4,5) = 1320*rc9
      aj(5,1) = 24*one
      aj(5,2) = 360*rc2
      aj(5,3) = 1680*rc4
      aj(5,4) = 5040*rc6
      aj(5,5) = 11880*rc8
c
c     Use LU decomposition (AG, April 1991) See Numerical Recipes.
c
      call sgeco(aj,5,5,indx,rcond,work)
      if (rcond .lt. 1.d-7) write(6,*) ' rcond too small:' ,rcond
c
      call sgesl(aj,5,5,indx,bj,0)
c
      alpha = bj(1)
      alpha1 = bj(2)
      alpha2 = bj(3)
      alpha3 = bj(4)
      alpha4 = bj(5)
c
c   start iteration loop to find delta (with gamma fixed)
c
      do 80 j = 1, 100
c
c   generate pseudo wavefunction-note missing factor exp(delta)
c
         do 30 k = 1, jrc
            rp = r(k)
            r2 = rp*rp
            polyr = r2*(((((alpha4*r2+alpha3)*r2+alpha2)*r2+alpha1)*r2+
     &              alpha)*r2+gamma)
            ar(k) = rp**lp*exp(polyr)
   30    continue
c
c   integrate pseudo charge density from r = 0 to rc
c
         ll = 2
         cdps = -ar(jrc)*ar(jrc)*rab(jrc)
         if (jrc .ne. 2*(jrc/2)) then
            do 40 k = jrc, 1, -1
               cdps = cdps + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
   40       continue
         else
            do 50 k = jrc, 4, -1
               cdps = cdps + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
   50       continue
            cdps = cdps - ar(4)*ar(4)*rab(4)
            cdps = cdps + 9*(ar(1)*ar(1)*rab(1)+3*ar(2)*ar(2)*rab(2)+
     &             3*ar(3)*ar(3)*rab(3)+ar(4)*ar(4)*rab(4))/8
         end if
         cdps = cdps/3
c
c        Calculate new delta (with gamma fixed), uses false position
c
         fdnew = log(cdrc/cdps) - 2*delta
         if (abs(fdnew) .lt. errmin) then
            v0pp = 8*((2*one*(lp-one)+5*one)*alpha+gamma**2)
c
            return
c
         end if
c
         if (j .eq. 1) then
            ddelta = -pfive
         else
            ddelta = -fdnew*ddelta/(fdnew-fdold)
         endif
         delta = delta + ddelta
c
         bj(1) = bj1 - delta
         bj(2) = bj2
         bj(3) = bj3
         bj(4) = bj4
         bj(5) = bj5
c
         call sgesl(aj,5,5,indx,bj,0)
c
         alpha = bj(1)
         alpha1 = bj(2)
         alpha2 = bj(3)
         alpha3 = bj(4)
         alpha4 = bj(5)
c
         fdold = fdnew
c
   80 continue
      v0pp = 1.d-60
      write(6,9000)
 9000 format(//'error in gamfind (aka v0pp) - delta not found')
      call ext(860+lp)
c
      end

c
c $Id: chg_mism.f,v 1.2 1997/05/22 17:32:04 wdpgaara Exp $
c
c $Log: chg_mism.f,v $
c Revision 1.2  1997/05/22 17:32:04  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      double precision function chg_mism(xdelta)
c
c     chg_mism gives the discrepancy (one half the log of their
c     ratio) between the pseudo and all-electron charge densities
c     inside rc when the first polynomial coefficient is xdelta.
c
      implicit none
c
      include 'radial.h'
      include 'nonlinear.h'
      include 'coeffs.h'
      include 'linear.h'
c
      double precision xdelta
c
      double precision ar(nrmax), bj(5)
c
      double precision polyr, rp, r2, cdps
      integer k, ll
c
c     Find the {alpha} set for this xdelta and the current gamma.
c
      call genrhs(gamma,xdelta,bj)
      call sgesl(alin,5,5,indx,bj,0)
c
         alpha = bj(1)
         alpha1 = bj(2)
         alpha2 = bj(3)
         alpha3 = bj(4)
         alpha4 = bj(5)
c
c   Generate the pseudo wavefunction
c   (note that exp(xdelta) is put at the end )
c
         do 30 k = 1, jrc
            rp = r(k)
            r2 = rp*rp
            polyr = r2*(((((alpha4*r2+alpha3)*r2+alpha2)*r2+alpha1)*r2+
     &              alpha)*r2+gamma)
            ar(k) = rp**lp*exp(polyr)
   30    continue
c
c   Integrate pseudo charge density from r = 0 to rc
c
         ll = 2
         cdps = -ar(jrc)*ar(jrc)*rab(jrc)
         if (mod(jrc,2) .ne. 0) then
            do 40 k = jrc, 1, -1
               cdps = cdps + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
   40       continue
         else
            do 50 k = jrc, 4, -1
               cdps = cdps + ll*ar(k)*ar(k)*rab(k)
               ll = 6 - ll
   50       continue
            cdps = cdps - ar(4)*ar(4)*rab(4)
            cdps = cdps + 9*(ar(1)*ar(1)*rab(1)+3*ar(2)*ar(2)*rab(2)+
     &             3*ar(3)*ar(3)*rab(3)+ar(4)*ar(4)*rab(4))/8
         end if
         cdps = cdps/3
c
         chg_mism = log(cdrc/cdps) - 2*xdelta
c
         return
c  
         end
c
      subroutine dsolv2(iter,iconv,id,nfirst,nlast,n_of_core_orbs,nn)
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'charge.h'
      include 'elecpot.h'
      include 'energy.h'
c
c  dsolv2 finds the (non) relativistic wave function using
c  difnrl to integrate the Schroedinger equation or
c  difrel to integrate the Dirac equation.
c  The energy level from the previous iteration are used
c  as initial guesses, and it must therefore be reasonably
c  accurate.
c
C     .. Parameters ..
      double precision zero, smev
      parameter (zero=0.D0,smev=1.D-4)
C     ..
C     .. Scalar Arguments ..
      integer iconv, iter, nfirst, nlast, n_of_core_orbs
      character id*1
C     ..
c
c     Watch out for nn: it is not always equal to 'no'. 
c     In particular, some nn(i) could be zero.
c     The same is true for id...
c
      integer nn(*)
C     ..
C     .. Local Scalars ..
      integer i, iflag, j, llp, lp
C     ..
C     .. External Subroutines ..
      external difnrl, difrel, orban
      logical leqi
      external leqi
C     ..
      double precision ar(nrmax), br(nrmax), v(nrmax),
     &                 orb_charge(nrmax)
C     ..
c
c  Initialize arrays for charge density.
c
      do 10 i = 1, nr
         cdd(i) = zero
         cdu(i) = zero
   10 continue
      if (ifcore .ne. 1) then
         do 30 i = 1, nr
            cdc(i) = zero
   30    continue
      end if
c
c  Start the loop over orbitals.
c  Note that spin zero is treated as down.
c
      do 130 i = nfirst, nlast
         if (nn(i) .le. 0) go to 130
         if (zo(i) .eq. 0.0D0 .and. iconv .eq. 0) go to 130
         if (ev(i) .ge. 0.0D0) ev(i) = -smev
c
c  Set up the potential, set the wave function array to zero-ar.
c
         lp = lo(i) + 1
         llp = lo(i)*lp
         do 40 j = 1, nr
            ar(j) = zero
   40    continue
         if (down(i)) then
            do 50 j = 2, nr
               v(j) = viod(lp,j)/r(j) + vid(j)
   50       continue
         else
            do 60 j = 2, nr
               v(j) = viou(lp,j)/r(j) + viu(j)
   60       continue
         end if
         if (.not. leqi(id,'r')) then
            do 70 j = 2, nr
               v(j) = v(j) + llp/r(j)**2
   70       continue
         end if
c
c  Call the integration routine.
c
         if (.not. leqi(id,'r')) then
            call difnrl(iter,i,v,ar,br,nn(i),lo(i),so(i),ev(i),iflag)
         else
            call difrel(iter,i,v,ar,br,nn(i),lo(i),so(i),ev(i))
         end if
c
c  Add to the charge density.
c
         if (leqi(id,'r')) then
            do 300 j = 1, nr
               orb_charge(j) = zo(i)*(br(j)*br(j)+ar(j)*ar(j))
 300        continue
         else
            do 320 j = 1, nr
               orb_charge(j) = zo(i) * ar(j)*ar(j)
 320        continue
         endif
c
            if (down(i)) then
               do 80 j = 1, nr
                  cdd(j) = cdd(j) + orb_charge(j)
   80          continue
            else
               do 90 j = 1, nr
                  cdu(j) = cdu(j) + orb_charge(j)
   90          continue
            end if

cag            if (ifcore .ne. 1 .and. i .le. n_of_core_orbs) then
            if ( i .le. n_of_core_orbs) then
              do 95 j = 1, nr
                 cdc(j) = cdc(j) + orb_charge(j)
   95         continue
            end if

c
c  Compute various quantitities if last iteration.
c
         if (iconv .eq. 1) call orban(i,id,ar,br,nn(i),lo(i),zo(i),
     &                                so(i),ev(i),ek(i),ep(i))
  130 continue
c
c  End loop over orbitals.
c
      return
c
      end
C
      subroutine etotal(nfirst,nlast)
c
      implicit none
c
      include 'orbital.h'
      include 'param.h'
      include 'energy.h'
c
c  etotal computes the total energy from the
c  electron charge density.
c
c  Only the orbitals between nfirst and nlast are considered.
c
c      etot(i)    i=1,10 contains various contributions to the total
c                 energy.
c                 (1)   sum of eigenvalues ev
c                 (2)   sum of orbital kinetic energies ek
c                 (3)   el-ion interaction from sum of orbital
c                       potential energies ep
c                 (4)   electrostatic el-el interaction  (from velect)
c                 (5)   vxc (exchange-correlation) correction to sum
c                       of eigenvalues                   (from velect)
c                 (6)   3 * vc - 4 * ec
c                       correction term for virial theorem
c                       when correlation is included     (from velect)
c                 (7)   exchange and correlation energy  (from velect)
c                 (8)   kinetic energy from eigenvalues  (1,3,4,5)
c                 (9)   potential energy
c                 (10)  total energy
c
c
c      sum up eigenvalues ev, kinetic energies ek, and
c      el-ion interaction ep
c
C     .. Parameters ..
      double precision zero
      parameter (zero=0.D0)
C     ..
      integer nfirst, nlast
c
C     .. Local Scalars ..
      double precision vsum
      integer i
      character*2 id
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs
C     ..
      etot(1) = zero
      etot(2) = zero
      etot(3) = zero
      do 10 i = nfirst, nlast
         etot(1) = etot(1) + zo(i)*ev(i)
         etot(2) = etot(2) + zo(i)*ek(i)
         etot(3) = etot(3) + zo(i)*ep(i)
   10 continue
c
c   kinetic energy
c
      etot(8) = etot(1) - etot(3) - 2*etot(4) - etot(5)
c
c   potential energy
c
      etot(9) = etot(3) + etot(4) + etot(7)
c
c      total energy
c
      etot(10) = etot(1) - etot(4) - etot(5) + etot(7)
c
c   printout
c
      write(6,9000) nameat
 9000 format(//1x,a2,' output data for orbitals',/1x,28('-'),
     &      //' nl    s      occ',9x,'eigenvalue',4x,'kinetic energy',
     &      6x,'pot energy',/)
c
      id = " "
      do 20 i = nfirst, nlast
         if (i .ge. ncp) id="&v"
         write(6,9010) no(i), il(lo(i)+1), so(i), zo(i), ev(i),
     &     ek(i), ep(i), id
 9010    format(1x,i1,a1,f6.1,f10.4,3f17.8,2x,a2)
   20 continue
c
      write(6,'(a)') '---------------------------- &v'
      write(6,9020) (etot(i),i=1,10)
 9020 format(//' total energies',/1x,14('-'),
     &      //' sum of eigenvalues        =',f18.8,
     &      /' kinetic energy from ek    =',f18.8,
     &      /' el-ion interaction energy =',f18.8,
     &      /' el-el  interaction energy =',f18.8,
     &      /' vxc    correction         =',f18.8,
     &      /' virial correction         =',f18.8,
     &      /' exchange + corr energy    =',f18.8,
     &      /' kinetic energy from ev    =',f18.8,
     &      /' potential energy          =',f18.8,/1x,45('-'),
     &      /' total energy              =',f18.8)
c
      if (job .ge. 4 .or. abs(zsh) .gt. 0.00001D0) return
c
c   virial theorem
c
      vsum = 2*etot(8) + etot(9) + etot(6)
      write(6,9030) 2*etot(8), etot(9), etot(6), vsum
 9030 format(//' virial theorem(nonrelativistic)',/1x,14('-'),
     &      //' kinetic energy  *  2      =',f18.8,
     &      /' potential energy          =',f18.8,
     &      /' virial correction         =',f18.8,/1x,45('-'),
     &      /' virial sum                =',f18.8)
c
      return
c
      end
C
c $Id: ext.f,v 1.2 1997/05/22 17:32:13 wdpgaara Exp $
c
c $Log: ext.f,v $
c Revision 1.2  1997/05/22 17:32:13  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
C
C
      subroutine ext(i)
c
c  Stops program in case of errors or completion.
c
c  i is a stop parameter
c   000-099 main (0 is normal exit)
c   100-199 input
c   200-299 charge
c   300-399 vionic
c   400-499 velect
c   500-599 dsolv1
c   600-699 dsolv2 (including difnrl and difrel)
c   700-799 etotal
c   800-899 pseudo, pseudk, pseudt and pseudv
c
C     .. Scalar Arguments ..
      integer i
C     ..
      if (i .ne. 0) write(6,FMT=9000) i
 9000 format('stop parameter =',i3)
      close(unit=1)
      close(unit=3)
      close(unit=5)
      close(unit=6)
c
      stop
c
      end
c
cGuima modificado
      subroutine Input(maxit)
cGuima termina modificacao
c

      implicit none
cGuima modificado
      integer maxit
ccGuima termina modificaccao
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'compat.h'
      include 'input.h'
c
c  subroutine to read input parameters
c
C     .. Parameters ..
      double precision one, zero, pfive
      parameter (one=1.D0,zero=0.D0,pfive=0.5D0)
c
C     .. Local Scalars ..
      double precision rmax, sc, si, zcore, zd, zu,  zval
      integer i, j, li, ni
      character type*2, flavor*3, compat_str*20
C     ..
C     .. Local Arrays ..
      integer lc(15), nc(15), nomin(0:4)
C     ..
C     .. External Functions ..
      double precision nucl_z
      logical leqi
      external nucl_z, leqi
C     ..
C     .. External Subroutines ..
c      external ext, cal_date
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, exp
C     ..
c
c  data for orbitals:
c
c              1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p,4f,5d,6s,6p
c
      data nc /1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 4, 5, 6, 6/
      data lc /0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1/
c
c     This statement now goes in a block data subprogram
cccc      data il /'s', 'p', 'd', 'f', 'g'/
C     ..
c
      do 10 i = 0, 4
         nomin(i) = 10
   10 continue
      do 20 i = 1, norbmx
         no(i) = 0
         lo(i) = 0
         so(i) = zero
         zo(i) = zero
   20 continue
c
c
c  read the type of calculation and title card
c   job =
c   ae = 0 all electron calculation
c   pg = 1 pseudopotential generation w/o core correction
c   pe = 2 pseudopotential generation w/  core correction exchange
c   ph = 3 pseudopotential generation w/  core correction hartree
c   pt = 4 pseudopotential test
c   pm = 5 pseudopotential test + valence charge modify
c
      job = -1
      read(5,9000,err=999,end=999) type, title
 9000 format(3x,a2,5a10)
c
c  if type = ' ' , no more data, program ends
c
      if (type .eq. 'ae') then
         job = 0
      else if (type .eq. 'pg') then
         job = 1
      else if (type .eq. 'pe') then
         job = 2
      else if (type .eq. 'ph') then
         job = 3
      else if (type .eq. 'pt') then
         job = 4
      else if (type .eq. 'pm') then
         job = 5
      else
         job = -1
c
         return
c
      end if
c
      ifcore = job - 1
c
c  njtj  ***  major modification  start  ***
c  There are seven ways to generate the pseudopotential :
c    flavor = van Vanderbilt
c    flavor = tam Troullier and Martins
c    flavor = ker (yes) Kerker
c    flavor = hsc (no)  Hamann Schluter and Chiang
c    flavor = min (oth) datafile made for minimization
c    flavor = bhs Bachelet, Hamann and Schluter
c    flavor = tm2 Improved Troullier and Martins
c
c
c     Flavor should not be necessary for the tests...
c
      if ((job .gt. 0) .and. (job .lt. 4)) then
         read(5,9010) flavor, logder_radius
 9010    format(8x,a3,f9.3)
         if (leqi(flavor,'tm2')) then
            scheme = 6
         else if (leqi(flavor,'bhs')) then
            scheme = 5
         else if (leqi(flavor,'oth') .or. leqi(flavor,'min')) then
            scheme = 4
         else if (leqi(flavor,'tbk') .or. leqi(flavor,'tam')) then
            scheme = 2
         else if (leqi(flavor,'yes') .or. leqi(flavor,'ker')) then
            scheme = 1
         else if (leqi(flavor,'no ') .or. leqi(flavor,' no')
     &            .or. leqi(flavor,'hsc') ) then
            scheme = 0
         else
            write(6,9020) flavor
            call ext(150)
         end if
      end if
 9020 format(//'error in input - flavor =',a3,' unknown')
c  njtj  ***  major modification end  ***
c
c   read element name, correlation type, polarization flag...
c   ispp = ' ' - nonspin calculation
c   ispp = s  - spin polarized calculation
c   ispp = r  - relativistic calculation
c
c   ... and a compatibility string (obsolete -- to be removed)
c
      read(5,9030) nameat, icorr, ispp, compat_str
 9030 format(3x,a2,3x,a2,a1,1x,a20)
c
      call compat_params(compat_str)
c
      if (.not.( leqi(ispp,'s') .or. leqi(ispp,'r') ) ) ispp = ' '
c
c     Not all the correlation schemes can be used for a
c     spin-polarized calculation.
c
      if ( leqi(ispp,'s') .and.
     &     ( leqi(icorr,'xa') .or. leqi(icorr,'wi') .or.
     &       leqi(icorr,'hl') .or. leqi(icorr,'gr'))
     &   )   ispp = ' '
c
      normal = leqi(ispp,' ')
      relativistic = leqi(ispp,'r')
      polarized = leqi(ispp,'s')
c
c  njtj   ***  major modification start  ***
c   Floating point comparison error modification.
c   Read the atomic number (nuclear charge),
c   shell charge and radius (added to the nuclear potential),
c   and radial grid parameters.
c
      read(5,9040) znuc, zsh, rsh, rmax, aa, bb
 9040 format(6f10.3)
      if (abs(znuc) .le. 0.00001D0) znuc = nucl_z(nameat)
      if (job .lt. 4) then
c
c   set up grid
c
         if (abs(rmax) .lt. 0.00001D0) rmax = rmax_def
         if (abs(aa) .lt. 0.00001D0) aa = aa_def
         if (abs(bb) .lt. 0.00001D0) bb = bb_def
         a = exp(-aa)/znuc
         b = 1/bb
         do 30 i = 1, nrmax
            r(i) = a*(exp(b*(i-1))-1)
            rab(i) = (r(i)+a)*b
            if (r(i) .gt. rmax) go to 40
   30    continue
c
         write(6,9050)
 9050    format(/' error in input - arraylimits',
     &      ' for radial array exceeded',/)
c
         call ext(100)
c
   40    continue
c
         nr = i - 1
c
      end if
c  njtj  ***  major modification end  ***
c
c   read the number of core and valence orbitals
c
c
      read(5,9060) ncore, nval
 9060 format(2i5,2f10.3)
      if (ncore .gt. 15) then
         write(6,9070)
         call ext(101)
      end if
 9070 format(//'error in input - max number of core orbitals','is 15')
c
c   compute occupation numbers and orbital energies for the core
c
      zcore = zero
      if (ncore .eq. 0) go to 70
      sc = zero
      if (ispp .ne. ' ') sc = -pfive
      norb = 0
      do 60 i = 1, ncore
         do 50 j = 1, 2
            if (ispp .eq. ' ' .and. j .eq. 2) go to 50
            norb = norb + 1
            no(norb) = nc(i)
            lo(norb) = lc(i)
            so(norb) = sc
            zo(norb) = 2*lo(norb) + 1
            if (ispp .eq. ' ') zo(norb) = 2*zo(norb)
            if (ispp .eq. 'r') zo(norb) = 2*(lo(norb)+sc) + 1
            zcore = zcore + zo(norb)
            if (abs(zo(norb)) .lt. 0.1D0) norb = norb - 1
            if (ispp .ne. ' ') sc = -sc
   50    continue
   60 continue
      ncore = norb
c
c   for the valence orbitals
c
   70 continue
      if (job .ge. 4) ncore = 0
      norb = ncore
      zval = zero
      if (nval .eq. 0) go to 130
      do 90 i = 1, nval
         read(5,9060) ni, li, zd, zu
         si = zero
         if (ispp .ne. ' ') si = pfive
         do 80 j = 1, 2
            if (ispp .eq. ' ' .and. j .eq. 2) go to 80
            norb = norb + 1
            if (ispp .ne. ' ') si = -si
            no(norb) = ni
            lo(norb) = li
            so(norb) = si
            zo(norb) = zd + zu
            if (ispp .eq. 's') then
               if (si .lt. 0.1D0) then
                  zo(norb) = zd
               else
                  zo(norb) = zu
               end if
            else if (ispp .eq. 'r') then
               zo(norb) = zo(norb)*(2*(li+si)+1)/(4*li+2)
            end if
            zval = zval + zo(norb)
            if (ispp .eq. 'r' .and. li+si .lt. zero) norb = norb - 1
            if (norb .eq. 0) go to 80
            if (nomin(lo(norb)) .gt. no(norb))
     &                               nomin(lo(norb)) = no(norb)
   80    continue
   90 continue
cGuima modificado
      read(5,*)maxit
cGuima termina modificaccao
c
c   Compute 'down' flag and
c   abort if two orbitals are equal
c
      nval = norb - ncore
      do 110 i = 1, norb
         down(i) = (so(i) .lt. 0.1D0)
         do 100 j = 1, norb
            if (i .le. j) go to 100
            if (no(i) .ne. no(j)) go to 100
            if (lo(i) .ne. lo(j)) go to 100
            if (abs(so(i)-so(j)) .gt. 0.001D0) go to 100
            write(6,9080) i
            call ext(110+i)
  100    continue
  110 continue
 9080 format(//'error in input - orbital ',i2,'is already occupied')
c
c   reduce n quantum number if pseudoatom
c
      if (job .ge. 4) then
         do 120 i = 1, nval
            no(i) = no(i) - nomin(lo(i)) + lo(i) + 1
  120    continue
      end if
  130 continue
      zion = znuc - zcore - zval
      zel = zval
      if (job .lt. 4) then
         zel = zel + zcore
      else
         znuc = znuc - zcore
      end if
c
c   find jobname and date and printout.
c
      ray(1) = 'ATM 3.2.2'
      call cal_date(ray(2))
cag
      ncp = ncore + 1
cag
 999  continue
      return
c
      end
c
      block data orb_init
      implicit none
      include 'orbital.h'
      data il /'s', 'p', 'd', 'f', 'g'/
      end

c
      subroutine Header
c
c     Writes out the header information
c     read by subroutine Input.
c
c     Most of the information is in
c     the standard program common blocks.
c     A few things in the input.h blocks.
c
c     Alberto Garcia, July 8, 2002
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'compat.h'
      include 'input.h'
c
      character*3 name
      integer i
      double precision xji, zero
      
      parameter(zero = 0.d0)
c
c   find jobname and date and printout.
c
      ray(1) = 'ATM 3.2.2'
c      call cal_date(ray(2))
c
c   printout
c
      write(6,9090) ray(1), ray(2), title, '&v&d'
 9090 format(1x,a10,a10,1x,5a10,1x,a4/60('-'),/)
      if (job .eq. 0) then
         write(6,9100) nameat
      else if (job .lt. 4) then
         write(6,9110) nameat
      else if (job .eq. 4) then
         write(6,9120) nameat
      else if (job .eq. 5) then
         write(6,9130) nameat
      end if
 9100 format(1x,a2,' all electron calculation ',/1x,27('-'),/)
 9110 format(1x,a2,' pseudopotential generation',/1x,29('-'),/)
 9120 format(1x,a2,' pseudopotential test',/1x,23('-'),/)
 9130 format(1x,a2,' pseudo test + charge mod ',/1x,27('-'),/)
      if (ispp .eq. 'r') then
         write(6,9140)
 9140    format(' r e l a t i v i s t i c ! !',/)
         name = '   '
      else if (ispp .eq. ' ') then
         name = 'non'
      else
         name = '   '
      end if
      write(6,9150) icorr, name
 9150 format(' correlation = ',a2,3x,a3,'spin-polarized',/)
      write(6,9160) znuc, ncore, nval, zel, zion
 9160 format(' nuclear charge             =',f10.6,
     &      /' number of core orbitals    =',i3,
     &      /' number of valence orbitals =',i3,
     &      /' electronic charge          =',f10.6,
     &      /' ionic charge               =',f10.6,//)
      if (zsh .gt. 0.00001D0) write(6,9170) zsh, rsh
 9170 format(' shell charge =',f6.2,' at radius =',f6.2,//)
      write(6,9180)
 9180 format(' input data for orbitals',
     &      //'  i    n    l    s     j     occ',/)
      xji = zero
      do 140 i = 1, norb
         if (ispp .eq. 'r') xji = lo(i) + so(i)
         write(6,9190) i, no(i), lo(i), so(i), xji, zo(i)
 9190    format(1x,i2,2i5,2f6.1,f10.4)
  140 continue
      if (job .lt. 4) write(6,9200) r(2), nr, r(nr), aa, bb
 9200 format(//' radial grid parameters',//' r(1) = .0 , r(2) =',e9.3,
     &      ' , ... , r(',i4,') =',f8.3,/' a =',f7.3,'  b =',f8.3,/)
c
      return
c
      end



C
      subroutine orban(iorb,id,ar,br,n,l,occup,spin,eigv,ekin,epot)
c
      implicit none
c
      include 'radial.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
c      
c  orban is used to analyze and printout data
c  about the orbital.
c
C     .. Parameters ..
      double precision ai, zero
      parameter (ai=2*137.0360411D0,zero=0.D0)
C     ..
C     .. Scalar Arguments ..
      integer iorb, l, n
      double precision spin, occup, eigv, ekin, epot 
      character id*1
C     ..
C     .. Array Arguments ..
      double precision ar(*), br(*)
C     ..
C     .. Local Scalars ..
      double precision ar2, arp, arpm, br2, deni, sa2
      integer i, i90, i99, ka, ll, llp, lp, nextr, nzero
      integer kj, ist
C     ..
C     .. Local Arrays ..
      double precision aextr(10), bextr(10), rextr(10), rzero(10)
C     ..
      logical leqi
      external leqi
C     ..
c
      ka = l + 1
      lp = ka
      if (spin .lt. 0.1D0 .and. l .ne. 0) ka = -l
c
c      compute zeroes and extrema
c
      nzero = 0
      nextr = 0
      rzero(1) = zero
      arp = br(2)
      if (leqi(id,'r')) then
         if (spin .lt. 0.1D0) then
            arp = ka*ar(2)/r(2) + (eigv-viod(lp,2)/r(2)-vid(2)+
     &            ai*ai)*br(2)/ai
         else
            arp = ka*ar(2)/r(2) + (eigv-viou(lp,2)/r(2)-viu(2)+
     &            ai*ai)*br(2)/ai
         end if
      end if
      do 20 i = 3, nr
         if (nextr .ge. n-l) go to 30
c
         if (ar(i)*ar(i-1) .le. zero) then
c
c            zero
c
             nzero = nzero + 1
             rzero(nzero) = (ar(i)*r(i-1)-ar(i-1)*r(i))/
     &                                            (ar(i)-ar(i-1))
         endif
c
         arpm = arp
         arp = br(i)
c
         if (leqi(id,'r')) then
            if (spin .lt. 0.1D0) then
               arp = ka*ar(i)/r(i) + (eigv-viod(lp,i)/r(i)-vid(i)+
     &               ai*ai)*br(i)/ai
            else
               arp = ka*ar(i)/r(i) + (eigv-viou(lp,i)/r(i)-viu(i)+
     &               ai*ai)*br(i)/ai
            end if
         end if
c
         if (arp*arpm .le. zero) then
c
c          extremum
c
           nextr = nextr + 1
           rextr(nextr) = (arp*r(i-1)-arpm*r(i))/(arp-arpm)
           aextr(nextr) = (ar(i)+ar(i-1))/2 -
     &                    (arp**2+arpm**2)*(r(i)-r(i-1))/(4*(arp-arpm))
           bextr(nextr) = br(i)
c
         endif
c
   20 continue
c
c   Find orbital kinetic and potential energy
c   the potential part includes only the interaction with
c   the nuclear part
c
   30 continue
      ekin = br(1)*br(1)*rab(1)
      epot = zero
      sa2 = zero
      lp = l + 1
      llp = l*lp
      ll = 2
      if (2*(nr/2) .eq. nr) ll = 4
      i90 = nr
      i99 = nr
      do 40 i = nr, 2, -1
         ar2 = ar(i)*ar(i)
         br2 = br(i)*br(i)
         deni = ar2
         if (leqi(id,'r')) deni = deni + br2
         ekin = ekin + ll*(br2+ar2*llp/r(i)**2)*rab(i)
         if (spin .lt. 0.1D0) then
            epot = epot + ll*deni*viod(lp,i)*rab(i)/r(i)
         else
            epot = epot + ll*deni*viou(lp,i)*rab(i)/r(i)
         end if
         ll = 6 - ll
         if (sa2 .le. 0.1D0) then
           sa2 = sa2 + deni*rab(i)
           if (sa2 .le. 0.01D0) i99 = i
           i90 = i
         endif
   40 continue
c
      ekin = ekin/3
      epot = epot/3
      if (leqi(id,'r')) ekin = zero
c
c     Printout
c
      write(6,9000) n, l, spin
 9000 format(/' n =',i2,'  l =',i2,'  s =',f4.1)
c
      write(6,9010) 'a extr    ', (aextr(i),i=1,nextr)
      if (leqi(id,'r')) write(6,9010) 'b extr    ', (bextr(i),i=1,nextr)
      write(6,9010) 'r extr    ', (rextr(i),i=1,nextr)
      write(6,9010) 'r zero    ', (rzero(i),i=1,nzero)
      write(6,9010) 'r 90/99 % ', r(i90), r(i99)
c
      if (eigv .eq. zero) then
         if (occup .ne. zero) then
            write(6,9020) occup
         else
            write(6,9030)
         end if
      end if
c
 9010 format(8x,a10,2x,8f8.3)
 9020 format(8x,'WARNING: This orbital is not bound',' and contains ',
     &      f6.4,' electrons!!')
 9030 format(8x,'WARNING:  This orbital is not bound!')
c
c
      if (job .ne. 0 .and. job .ne. 4) return
c
c    Plot valence wavefunctions if AE or PT job
c
      if (iorb .gt. ncore) then
c
c     Plot and make the wavefunction 'upright'
c
         ist = nint(sign(1.d0,ar(nr-85)))
c
         if (job .eq. 4) then
            kj = -1
         else 
            kj = 1
         endif
cGuima modifica
         call potrw(ar,r,nr,l,kj,ist,0.d0)
cGuima fim da modificaccao
      endif
c
      return
c
      end


C
      subroutine potran(i,vd,r,nr,zion,fourier_area)
c
      implicit none
c
      include 'plot.h'
c
c ***********************************************************
c *                                                         *
c *  The potential is fitted with a                         *
c *  second degree polynomial, which is multiplied with the *
c *  appropriate functions and then integrated by parts     *
c *  to find the fourier transform.                         *
c *                                                         *
c ***********************************************************
c
c  The potential times r is fitted to the polynominal
c  a + bx + cx^2 at every other point.
c
C     .. Parameters ..
      double precision zero, one
      parameter (zero=0.D0,one=1.D0)
C     ..
C     .. Scalar Arguments ..
      double precision zion, fourier_area
      integer i, nr
C     ..
C     .. Array Arguments ..
      double precision r(nr), vd(nr)
C     ..
C     .. Local Scalars ..
      double precision d1, d2, d3, q, q2, r0, rm, rp, v0, vline, vm, vp
      integer j, k
      character*9 filename
C     ..
C     .. Local Arrays ..
      double precision vql(100)
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, cos, sin
C     ..
      rm = zero
      vm = 2*zion
      do 10 k = 2, nr-1, 2
         r0 = r(k)
         v0 = r0*vd(k) + 2*zion
         rp = r(k+1)
         vp = rp*vd(k+1) + 2*zion
         d1 = 1/((rp-rm)*(r0-rm))
         d2 = 1/((rp-r0)*(rm-r0))
         d3 = 1/((r0-rp)*(rm-rp))
         a(k) = vm*d1 + v0*d2 + vp*d3
         b(k) = -vm*(r0+rp)*d1 - v0*(rm+rp)*d2 - vp*(rm+r0)*d3
         c(k) = vm*r0*rp*d1 + v0*rm*rp*d2 + vp*rm*r0*d3
         rm = rp
         vm = vp
   10 continue
c
c  Find the fourier transform q^2/4pi/zion*vql. Everything is
c  rescaled  by zion.  Integration is from q=1 to q=24
c
      do 30 j = 1, 94
         q = one/4*j
         q2 = q*q
         vql(j) = zero
         rm = zero
         do 20 k = 2, nr - 1, 2
            rp = r(k+1)
            vql(j) = vql(j) + (2*a(k)*rp+b(k))/q*sin(q*rp) -
     &               ((a(k)*rp+b(k))*rp+c(k)-2*a(k)/q2)*cos(q*rp) -
     &               (2*a(k)*rm+b(k))/q*sin(q*rm) +
     &               ((a(k)*rm+b(k))*rm+c(k)-2*a(k)/q2)*cos(q*rm)
            rm = rp
   20    continue
         vql(j) = vql(j)/2/zion - one
   30 continue
c
c  Print out the transforms( really q^2/(4pi*zion)*v(q) ) to
c  the current plot.dat file (unit=3) for latter plotting.
c
      write(filename,9900) i-1
 9900 format('PSPOTQ',i1)
      open(unit=3,file=filename,form='formatted',status='unknown')
c
      do 40 j = 1, 94
         write(3,9000) one/4*j, vql(j), 0.d0
   40 continue
 9000 format(1x,f7.4,3x,2f10.6)
c
      close(3)
c
cag      write(3,9010) i
c 9010 format(1x,'marker fn',i1)
c
c     Compute the absolute area 
c
      vline = 7*one + 32*abs(vql(1)) + 12*abs(vql(2)) + 32*abs(vql(3)) +
     &        7*abs(vql(4))
      do 50 j = 4, 88, 4
         vline = vline + 7*abs(vql(j)) + 32*abs(vql(j+1)) +
     &           12*abs(vql(j+2)) + 32*abs(vql(j+3)) + 7*abs(vql(j+4))
   50 continue
      fourier_area = vline/90
      write(6,9020) i - 1, fourier_area
 9020 format(1x,'The Fourier(q^2/(4pi*zion)*V(q)) absolute',
     &      ' area for l=',i1,' is ',f10.6)
c
      return
c
      end
c
      subroutine potrv(vd,r,nr,k,zion)
c
c  Step size of 0.01 is adjustable as seen fit to give
c  a reasonable plot.
c
C     .. Scalar Arguments ..
      integer k, nr
      double precision zion
C     ..
C     .. Array Arguments ..
      double precision r(nr), vd(nr)
C     ..
C     .. Local Scalars ..
      double precision step
      integer j
      character filename*7
C     ..
cag
      write(filename,9900) k
 9900 format('PSPOTR',i1)
      open(unit=3,file=filename,form='formatted',status='unknown')
cag
c     Write out r, V(r) and the Coulomb potential (for reference)
c
      step = 0.0D0
      do 10 j = 5, nr
         if (r(j) .ge. step) then
            write(3,9000) r(j), vd(j), -2*zion/r(j)
            step = step + 0.01D0
         end if
   10 continue
 9000 format(1x,f7.4,3x,f10.5,g20.10)
c
      close(unit=3)
cag
      return
c
      end
c
      subroutine potrvs(vd,r,nr,k)
c
c  Generates file for plotting of Screened ionic pseudopotentials
c  Step size of 0.01 is adjustable as seen fit to give
c  a reasonable plot.
c
C     .. Scalar Arguments ..
      integer k, nr
C     ..
C     .. Array Arguments ..
      double precision r(nr), vd(nr)
C     ..
C     .. Local Scalars ..
      double precision step
      integer j
      character filename*10
C     ..
      write(filename,9900) k
 9900 format('SCRPSPOTR',i1)
      open(unit=3,file=filename,form='formatted',status='unknown')

c     Write out r and the screened Vps(r)
c
      step = 0.0D0
      do 10 j = 5, nr
         if (r(j) .ge. step) then
            write(3,9000) r(j), vd(j)
            step = step + 0.01D0
         end if
   10 continue
 9000 format(1x,f7.4,3x,f10.5)
c
      close(unit=3)
      return
c
      end

cGuima inserccao
      subroutine totalw(i)
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
      include 'energy.h'
      character filename*7
      write(filename,9900) i
 9900 format('VTOTAL',i1)
      open(unit=3,file=filename,form='formatted',status='unknown')
 8030 format(4(g20.12))
         write(3,'(a)')' Raios'
         write(3,8030) (r(j),j=2,nr)
         write(3,'(a)') 'Down potential follows (l on next line)'
         write(3,'(i2)') i
         write(3,8030) (viod(i,j)+vid(j)*r(j),j=2,nr)
         write(3,'(a)') 'Up potential follows (l on next line)'
         write(3,'(i2)') i
         write(3,8030) (viou(i,j)+viu(j)*r(j),j=2,nr)
      close(unit=3)
      return
      end
cGuima termina inserccao


c
      subroutine potrw(vd,r,nr,k,kj,ist,rc)
c
c  Step size of 0.01 is adjustable as seen fit to give
c  a reasonalble plot.
c
C     .. Parameters ..
      double precision zero, pzf
cGuima modifica
      parameter (zero=0.D0,pzf=0.0D0)
cGuima fim da modificaccao
C     ..
C     .. Scalar Arguments ..
      integer ist, k, kj, nr
      double precision rc
C     ..
C     .. Array Arguments ..
      double precision r(nr), vd(nr)
C     ..
C     .. Local Scalars ..
      double precision step
      integer j
      character filename*7
cGuima insere
      real*8 med1,medr2,medr4,prev,rprev,carga(nr),selfp(nr),autoen
cGuima termina insere
      logical lexist
C     ..
c
      if (kj .eq. 0) then
         write(filename,9900) k
 9900 format('PSWFNR',i1)
      else if (kj .eq. -1) then
         write(filename,9930) k
 9930 format('PTWFNR',i1)
      else
         write(filename,9910) k
 9910 format('AEWFNR',i1)
         call totalw(k)
      endif
c
      open(unit=3,file=filename,form='formatted',status='unknown')
c
c     Write out r, the wavefunction, and rc (kludge to pass it to
c     the plotting program)
c
      step = zero
      do 10 j = 2, nr
         if (r(j) .ge. step) then
            write(3,9000) r(j), vd(j)*ist, rc
            step = step + pzf
         end if
   10 continue
cGuima modifica
 9000 format(1x,f17.12,3x,e21.12,2x,f8.4)
c  Finds the potential due the number density of the orbital
      carga(1)=0.d0
      selfp(1)=0.d0
      r(1)=0.d0
      prev=0.d0
      do j=2,nr
         carga(j)=carga(j-1)+(r(j)-r(j-1))*0.5d0*(vd(j)**2+prev*r(j-1))
         selfp(j)=selfp(j-1)+(r(j)-r(j-1))*0.5d0*(vd(j)**2/r(j)+prev)
         prev=vd(j)**2/r(j)
      enddo
      selfp(1)=selfp(nr)
      do j=2,nr
         selfp(j)=selfp(nr)-selfp(j)+carga(j)/r(j)
      enddo
cGuima fim da modificaccao
      close(unit=3)
cGuima    Calculo das medias <1>, <r**2>, <r**4>,auto-energia
      med1=0.d0
      medr2=0.d0
      medr4=0.d0
      prev=0.d0
      rprev=0.d0
      do j=2,nr
         med1=med1+(r(j)-rprev)*0.5d0*(vd(j)**2+prev**2)
         medr2=medr2+(r(j)-rprev)*0.5d0*(vd(j)**2*r(j)**2+
     #      prev**2*rprev**2)
         medr4=medr4+(r(j)-rprev)*0.5d0*(vd(j)**2*r(j)**4+
     #      prev**2*rprev**4)
         prev=vd(j)
         rprev=r(j)
      enddo
      prev=0.d0
      rprev=0.d0
      autoen=0.d0
      do j=2,nr
         autoen=autoen+(r(j)-rprev)*0.5d0*(vd(j)**2*selfp(j)+prev)
         prev=vd(j)**2*selfp(j)
         rprev=r(j)
      enddo
      write(6,'(3f15.8,a)')medr2/med1,medr4/med1,autoen/med1**2,
     #      ' <r**2>,<r**4>,self-energy'
      write(27,'(3f15.8,a)')medr2/med1,medr4/med1,autoen/med1**2,
     #      ' <r**2>,<r**4>,self-energy'
      do j=2,nr;write(27,*)vd(j);enddo
cGuima   Termina o c�lculo das medias.
c
cag
      return
c
      end


C
      subroutine prdiff(nconf,econf,jobold)
      implicit none
c
c   Prints out the energy differences between
c   different atomic configurations.
c
c   njtj  ***  modifications  ***
c     econf is able to handle larger numbers
c     of configurations.
c   njtj  ***  modifications  ***
c
c
C     .. Scalar Arguments ..
      integer nconf
      integer jobold  
C     ..
C     .. Array Arguments ..
      double precision econf(100)
C     ..
      integer nconfmax
      parameter (nconfmax = 20)

C     .. Local Scalars ..
      integer i, j, iu, nconf_ae, n_excitations
      double precision absval, maxabs, norm1, norm2
      logical ae_found
C     ..
C     .. Local Arrays ..
      double precision econf_ae(nconfmax,nconfmax)

C     ..
      write(6,9000) (i,i=1,nconf)
      do 10 i = 1, nconf
         write(6,9010) i, (econf(i)-econf(j),j=1,i)
   10 continue
 9000 format(/' &d total energy differences in series',
     $      //,' &d',2x,9i9)
 9010 format(' &d',1x,i2,1x,9f9.4)
      write(6,'(/,a,/)') '*----- End of series ----* spdfg &d&v'
c
c     Write external files for easier interfacing
c
      if (jobold .eq. 0) then      ! All-electron series just done

         call get_unit(iu)
         open(iu,file='AE_ECONF',form='formatted',status='unknown')
         rewind(iu)
         write(iu,"(i4)") nconf
         write(iu,*) ((econf(i)-econf(j),j=1,i-1),i=2,nconf)
         close(iu)

      else if (jobold .eq. 4) then   ! Ps test series just done

         call get_unit(iu)
         open(iu,file='PT_ECONF',form='formatted',status='unknown')
         rewind(iu)
         write(iu,"(i4)") nconf
         write(iu,*) ((econf(i)-econf(j),j=1,i-1),i=2,nconf)
         close(iu)
c
c        Check whether there is an AE_ECONF file and compute the differences
c
         inquire(file="AE_ECONF",exist=ae_found)
         if (ae_found) then
            call get_unit(iu)
            open(iu,file='AE_ECONF',form='formatted',status='old')
            rewind(iu)
            read(iu,*) nconf_ae
            if (nconf .ne. nconf_ae) STOP "nconf_ae_pt"
            if (nconf .ge. nconfmax) STOP "nconf_ae_overflow"
            read(iu,*) ((econf_ae(i,j),j=1,i-1),i=2,nconf)
            close(iu)
c
c           Write a new file
c
            call get_unit(iu)
            open(iu,file='ECONF_DIFFS',
     $             form='formatted',status='unknown')
            rewind(iu)
c
c           Compute AE-PT differences in excitation energies,
c           and print them, together with the maximum, mean_abs,
c           and root-mean-square values.
c
            n_excitations = 0
            maxabs = -1.0d0
            norm1 = 0.0d0
            norm2 = 0.0d0
            do i=2,nconf
               do j=1,i-1
                  econf_ae(i,j) = econf_ae(i,j) - (econf(i)-econf(j))
                  absval = abs(econf_ae(i,j))
                  if (absval .gt. maxabs) maxabs = absval
                  norm1 = norm1 + absval
                  norm2 = norm2 + absval*absval
                  n_excitations =  n_excitations + 1
               enddo
            enddo
            norm1 = norm1 / n_excitations
            norm2 = sqrt( norm2 / n_excitations)
            write(iu,"(i4)") n_excitations
            write(iu,"(40f10.5)") ((econf_ae(i,j),j=1,i-1),i=2,nconf)
            write(iu,"(3f10.5)") maxabs, norm1, norm2
            close(iu)
            
         endif ! found AE_ECONF
      endif
c
      end
c
c $Id: string.f,v 1.3 1999/02/26 14:26:47 wdpgaara Exp $
c
c $Log: string.f,v $
c Revision 1.3  1999/02/26 14:26:47  wdpgaara
c Cosmetic changes.
c
c Revision 1.2  1997/05/22  17:32:33  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:55  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      SUBROUTINE CHRLEN(STRING,NCHAR,LCHAR)
C
C***********************************************************************
C
C  CHRLEN accepts a STRING of NCHAR characters and returns LCHAR,
C  the length of the string up to the last nonblank, nonnull.
C
      CHARACTER STRING*(*)
      Integer nchar, lchar
c
      Integer i, ncopy 
C
      NCOPY=NCHAR
      IF(NCOPY.LE.0)NCOPY=LEN(STRING)
C
      DO 10 I=1,NCOPY
        LCHAR=NCOPY+1-I
        IF(STRING(LCHAR:LCHAR).NE.' '.AND.
     *     STRING(LCHAR:LCHAR).NE.CHAR(0))RETURN
 10     CONTINUE
      LCHAR=0
      RETURN
      END
c
      subroutine chrcap(string,nchar)
C
C***********************************************************************
C
C  CHRCAP accepts a STRING of NCHAR characters and replaces
C  any lowercase letters by uppercase ones.
C
C
C     .. Scalar Arguments ..
      integer nchar
      character string*(*)
C     ..
C     .. Local Scalars ..
      integer i, itemp, ncopy
C     ..
C     .. Intrinsic Functions ..
      intrinsic char, ichar, len, lge, lle
C     ..
      ncopy = nchar
      if (ncopy .le. 0) ncopy = len(string)
      do 10 i = 1, ncopy
C
         if (lge(string(i:i),'a') .and. lle(string(i:i),'z')) then
            itemp = ichar(string(i:i)) + ichar('A') - ichar('a')
            string(i:i) = char(itemp)
         end if
   10 continue
c
      return
c
      end
c
      logical function leqi(strng1,strng2)
C
C***********************************************************************
C
C  Case-insensitive lexical equal-to comparison
C
C
C     .. Scalar Arguments ..
      character strng1*(*), strng2*(*)
C     ..
C     .. Local Scalars ..
      integer i, len1, len2, lenc
      character s1*1, s2*1
C     ..
C     .. External Subroutines ..
      external chrcap
C     ..
C     .. Intrinsic Functions ..
      intrinsic len, min
C     ..
      len1 = len(strng1)
      len2 = len(strng2)
      lenc = min(len1,len2)
C
      leqi = .FALSE.
      do 10 i = 1, lenc
         s1 = strng1(i:i)
         s2 = strng2(i:i)
         call chrcap(s1,1)
         call chrcap(s2,1)
         if (s1 .ne. s2) return
   10 continue
C
      if (len1 .gt. lenc .and. strng1(lenc+1:len1) .ne. ' ') return
      if (len2 .gt. lenc .and. strng2(lenc+1:len2) .ne. ' ') return
      leqi = .TRUE.
c
      return
c
      end
c
      subroutine loc_des(message)
c
c     Processes message to locate the delimiters % and $ used
c     in the warnp routine.
c
c     Alberto Garcia, Feb 1, 1991
c
C     .. Scalar Arguments ..
      character message*(*)
C     ..
C     .. Scalars in Common ..
      integer form_length
      character form_spec*200
C     ..
C     .. Local Scalars ..
      integer dol_pos, pct_pos
      character work*200
C     ..
C     .. Intrinsic Functions ..
      intrinsic index
C     ..
C     .. Common blocks ..
      common /fordes/form_spec, form_length
C     ..
      pct_pos = index(message,'%')
c
      form_spec = form_spec(1:form_length)//message(1:(pct_pos-1))//
     &            ''','
      form_length = form_length + (pct_pos-1) + 2
      work = message(pct_pos+1:)
c
      dol_pos = index(work,'$')
      form_spec = form_spec(1:form_length)//work(1:dol_pos-1)//','''
      form_length = form_length + (dol_pos-1) + 2
c
c        Return the rest of message
c
      message = ' '
      message = work(dol_pos+1:)
c
      return
c
      end
c
      SUBROUTINE ZBRAC(FUNC,X1,X2,SUCCES)
C
C
C     .. Parameters ..
      INTEGER NTRY
      PARAMETER (NTRY=50)
      DOUBLE PRECISION FACTOR
      PARAMETER (FACTOR=1.6D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X1,X2
      LOGICAL SUCCES
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION F1,F2
      INTEGER J
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      IF (X1.EQ.X2) THEN
          WRITE (*,FMT=*) 'You have to guess an initial range'
          STOP 'range'

      END IF

      F1 = FUNC(X1)
      F2 = FUNC(X2)
      SUCCES = .TRUE.
      DO 11 J = 1,NTRY
c
c     AG: avoid overflow
c
          IF (sign(1.d0,F1)*sign(1.d0,F2) .LT. 0.D0) RETURN
          IF (ABS(F1).LT.ABS(F2)) THEN
              X1 = X1 + FACTOR* (X1-X2)
              F1 = FUNC(X1)

          ELSE
              X2 = X2 + FACTOR* (X2-X1)
              F2 = FUNC(X2)
          END IF

   11 CONTINUE
      SUCCES = .FALSE.
      RETURN

      END
      SUBROUTINE ZBRAK(FX,X1,X2,N,XB1,XB2,NB)
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION X1,X2
      INTEGER N,NB
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION XB1(1),XB2(1)
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION FX
      EXTERNAL FX
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DX,FC,FP,X
      INTEGER I,NBB
C     ..
      NBB = NB
      NB = 0
      X = X1
      DX = (X2-X1)/N
      FP = FX(X)
      DO 11 I = 1,N
          X = X + DX
          FC = FX(X)
c
c     AG: avoid overflow
c
          IF (sign(1.d0,FC)*sign(1.d0,FP) .LT. 0.D0) THEN
              NB = NB + 1
              XB1(NB) = X - DX
              XB2(NB) = X
          END IF

          FP = FC
          IF (NBB.EQ.NB) RETURN
   11 CONTINUE
      RETURN

      END
      DOUBLE PRECISION FUNCTION ZBRENT(FUNC,X1,X2,TOL)
C
C
C     .. Parameters ..
      INTEGER ITMAX
      PARAMETER (ITMAX=100)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.D-8)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION TOL,X1,X2
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM
      INTEGER ITER
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN,SIGN
C     ..
      A = X1
      B = X2
      FA = FUNC(A)
      FB = FUNC(B)
c
c     AG: avoid overflow
c
      IF (sign(1.d0,FB)*sign(1.d0,FA) .GT. 0.D0) THEN
          WRITE (*,FMT=*) 'Root must be bracketed for ZBRENT.'
          STOP 'range'
      END IF

      FC = FB
      DO 11 ITER = 1,ITMAX
c
c     AG: avoid overflow
c
          IF (sign(1.d0,FB)*sign(1.d0,FC) .GT. 0.D0) THEN
CAG          IF (FB*FC.GT.0.D0) THEN
              C = A
              FC = FA
              D = B - A
              E = D
          END IF

          IF (ABS(FC).LT.ABS(FB)) THEN
              A = B
              B = C
              C = A
              FA = FB
              FB = FC
              FC = FA
          END IF

          TOL1 = 2.D0*EPS*ABS(B) + 0.5D0*TOL
          XM = .5D0* (C-B)
          IF (ABS(XM).LE.TOL1 .OR. FB.EQ.0.D0) THEN
              ZBRENT = B
              RETURN

          END IF

          IF (ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
              S = FB/FA
              IF (A.EQ.C) THEN
                  P = 2.D0*XM*S
                  Q = 1.D0 - S

              ELSE
                  Q = FA/FC
                  R = FB/FC
                  P = S* (2.D0*XM*Q* (Q-R)- (B-A)* (R-1.D0))
                  Q = (Q-1.D0)* (R-1.D0)* (S-1.D0)
              END IF

              IF (P.GT.0.D0) Q = -Q
              P = ABS(P)
              IF (2.D0*P.LT.MIN(3.D0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
                  E = D
                  D = P/Q

              ELSE
                  D = XM
                  E = D
              END IF

          ELSE
              D = XM
              E = D
          END IF

          A = B
          FA = FB
          IF (ABS(D).GT.TOL1) THEN
              B = B + D

          ELSE
              B = B + SIGN(TOL1,XM)
          END IF

          FB = FUNC(B)
   11 CONTINUE
      WRITE (*,FMT=*) 'ZBRENT exceeding maximum iterations.'
      ZBRENT = B
      STOP 'ITER'

      END

      DOUBLE PRECISION FUNCTION RTBIS(FUNC,X1,X2,XACC)
C
C
C     .. Parameters ..
      INTEGER JMAX
      PARAMETER (JMAX=40)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X1,X2,XACC
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DX,F,FMID,XMID
      INTEGER J
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      FMID = FUNC(X2)
      F = FUNC(X1)
c
c     AG: avoid overflow
c
      IF (sign(1.d0,F)*sign(1.d0,FMID) .GE. 0.D0) THEN
CAG      IF (F*FMID.GE.0.D0) THEN
          WRITE (*,FMT=*) 'Root must be bracketed for bisection.'
          STOP 'RTBIS'

      END IF

      IF (F.LT.0.D0) THEN
          RTBIS = X1
          DX = X2 - X1

      ELSE
          RTBIS = X2
          DX = X1 - X2
      END IF

      DO 11 J = 1,JMAX
          DX = DX*.5D0
          XMID = RTBIS + DX
          FMID = FUNC(XMID)
          IF (FMID.LE.0.D0) RTBIS = XMID
          IF (ABS(DX).LT.XACC .OR. FMID.EQ.0.D0) RETURN
   11 CONTINUE
      WRITE (*,FMT=*) 'too many bisections'
      STOP 'ITER'
      END
      SUBROUTINE BRAC(FUNC,X1,X2,SUCCES)
C
c     This is ZBRAC frm Numerical Recipes.
c     Changed name to brac to avoid recursion...
c
C
C     .. Parameters ..
      INTEGER NTRY
      PARAMETER (NTRY=50)
      DOUBLE PRECISION FACTOR
      PARAMETER (FACTOR=1.6D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X1,X2
      LOGICAL SUCCES
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION F1,F2
      INTEGER J
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      IF (X1.EQ.X2) THEN
          WRITE (*,FMT=*) 'You have to guess an initial range'
          STOP 'range'

      END IF

      F1 = FUNC(X1)
      F2 = FUNC(X2)
      SUCCES = .TRUE.
      DO 11 J = 1,NTRY
c
c     AG: avoid overflow
c
          IF (sign(1.d0,F1)*sign(1.d0,F2) .LT. 0.D0) RETURN
CAG          IF (F1*F2.LT.0.D0) RETURN
          IF (ABS(F1).LT.ABS(F2)) THEN
              X1 = X1 + FACTOR* (X1-X2)
              F1 = FUNC(X1)

          ELSE
              X2 = X2 + FACTOR* (X2-X1)
              F2 = FUNC(X2)
          END IF

   11 CONTINUE
      SUCCES = .FALSE.
      RETURN

      END
C
      DOUBLE PRECISION FUNCTION BRENT(FUNC,X1,X2,TOL)
C
c     This is ZBRENT from Numerical Recipes.
c     Changed name to avoid recursion...
c
C
C     .. Parameters ..
      INTEGER ITMAX
      PARAMETER (ITMAX=100)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.D-8)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION TOL,X1,X2
C     ..
C     .. Function Arguments ..
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM
      INTEGER ITER
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN,SIGN
C     ..
      A = X1
      B = X2
      FA = FUNC(A)
      FB = FUNC(B)
c
c     AG: avoid overflow
c
      IF (sign(1.d0,FB)*sign(1.d0,FA) .GT. 0.D0) THEN
CAG      IF (FB*FA.GT.0.D0) THEN
          WRITE (*,FMT=*) 'Root must be bracketed for BRENT.'
          STOP 'range'

      END IF

      FC = FB
      DO 11 ITER = 1,ITMAX
c
c     AG: avoid overflow
c
          IF (sign(1.d0,FB)*sign(1.d0,FC) .GT. 0.D0) THEN
CAG          IF (FB*FC.GT.0.D0) THEN
              C = A
              FC = FA
              D = B - A
              E = D
          END IF

          IF (ABS(FC).LT.ABS(FB)) THEN
              A = B
              B = C
              C = A
              FA = FB
              FB = FC
              FC = FA
          END IF

          TOL1 = 2.D0*EPS*ABS(B) + 0.5D0*TOL
          XM = .5D0* (C-B)
          IF (ABS(XM).LE.TOL1 .OR. FB.EQ.0.D0) THEN
              BRENT = B
              RETURN

          END IF

          IF (ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
              S = FB/FA
              IF (A.EQ.C) THEN
                  P = 2.D0*XM*S
                  Q = 1.D0 - S

              ELSE
                  Q = FA/FC
                  R = FB/FC
                  P = S* (2.D0*XM*Q* (Q-R)- (B-A)* (R-1.D0))
                  Q = (Q-1.D0)* (R-1.D0)* (S-1.D0)
              END IF

              IF (P.GT.0.D0) Q = -Q
              P = ABS(P)
              IF (2.D0*P.LT.MIN(3.D0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
                  E = D
                  D = P/Q

              ELSE
                  D = XM
                  E = D
              END IF

          ELSE
              D = XM
              E = D
          END IF

          A = B
          FA = FB
          IF (ABS(D).GT.TOL1) THEN
              B = B + D

          ELSE
              B = B + SIGN(TOL1,XM)
          END IF

          FB = FUNC(B)
   11 CONTINUE
      WRITE (*,FMT=*) 'BRENT exceeding maximum iterations.'
      BRENT = B
      RETURN

      END
c
c $Id: genrhs.f,v 1.2 1997/05/22 17:32:13 wdpgaara Exp $
c
c $Log: genrhs.f,v $
c Revision 1.2  1997/05/22 17:32:13  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine genrhs(gamma,delta,y)
c
c     Generates the right-hand sides for set of five linear equations
c     that determines alpha(i) given gamma and delta.
c
      implicit none
c
      include 'nonlinear.h'
c
      double precision gamma, delta
      double precision y(5)
c
      double precision a2, a3
c
      y(1) = log(arc/rc1**lp) - gamma * rc2 - delta
      y(2) = brc - lp/rc1 - 2 * gamma * rc1
      a2 = y(2) + 2*gamma*rc1
      y(3) = vrc - eigv - 2*lp/rc1 * a2 - a2**2 - 2*gamma
      a3 = y(3) + 2*gamma
      y(4) = vap + 2*lp/rc2 * a2 - 2*lp/rc1 * a3 - 2*a2*a3
      y(5) = vapp - 4*lp/rc3 * a2 + 4*lp/rc2 * a3 - 2*lp/rc1 * y(4) -
     &       2 * a3**2 - 2 * a2 * y(4)
c
      return
c
      end
c
      subroutine logder(nfirst,nlast,flag)
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
      include 'energy.h'
c
      include 'ode_blk.h'
      include 'ode_path.h'
c
c     Computes the logarithmic derivative (LD) as a function of energy.
c     The all-electron and pseudopotential results can be compared
c     to assess the transferability of the pseudopotential.
c
c     Note:  Non-relativistic wave functions are used.
c
c     Alberto Garcia, April 27, 1991
c
c-----
c     October 1995: Removed output to unit 6.
c
c     Number of energy values at which the LD is computed and
c     half energy range around the eigenvalue (~2 rydberg)
c
      integer npt_energ
      double precision half_range
      parameter (npt_energ=51,half_range=2.d0)
      double precision eps
      parameter (eps=1.d-7)
c
c     Orbitals to be processed: nfirst <= i <= nlast
c
      integer nfirst, nlast
c
c     Flag: 'AE' for all-electron
c           'PS' for pseudo
c
      character flag*2
c
      double precision y(2)
c
      integer i, j, jr0, k, lp, llp, nok, nbad
      double precision emin, emax, eigv, step, r0, ld, h1, hmin
      character filename*9
c
      external ode, rkqc
c
c    Ode integration operational parameters
c
      h1 = 0.1d0
      hmin = 0.d0
      kmax = 200
      dxsav = 0.1d0
c
ccccc      write(6,'(/,1x,a,1x,a2,/)') 'Logarithmic derivatives', flag
c
      do 100 i = nfirst, nlast
c
        eigv = ev(i)
        emin = eigv - half_range
        emax = eigv + half_range
        step = 2 * half_range / (npt_energ - 1)
c
         lp = lo(i) + 1
         llp = lo(i)*lp
c
         write(filename,9900) flag, 'LOGD', lo(i)
 9900    format(a2,a4,i1)
         open(unit=3,file=filename,form='formatted',status='unknown')
c
c        Radius at which the LD will be calculated. Fix it at
c        a grid point. (Numerical Recipes, p.90)
c
         r0 = logder_radius
         call locate(r,nr,r0,jr0)
         r0 = r(jr0)
c
ccccccc         write(6,9000) i, no(i), lo(i), so(i), lo(i)+so(i), r0
 9000    format(/,1x,i2,2i5,2f6.1,2x,'r0:',f5.3,/)
c
         if (down(i)) then
            do 50 j = 2, nr
               v(j) = viod(lp,j)/r(j) + vid(j) + llp/r(j)**2
   50       continue
         else
            do 60 j = 2, nr
               v(j) = viou(lp,j)/r(j) + viu(j) + llp/r(j)**2
   60       continue
         end if
c
        do 30 k = 0, npt_energ - 1
c
          energ = emin + step * k
c
c         Call the integration routine. Adaptive stepsize.
c
c         y(1): wavefunction   ( = r R )
c         y(2): first derivative
c
c         Initial values (at r(2), near r=0 and so R ~ r^l )
c
          y(1) = r(2) ** lp
          y(2) = lp * r(2)**(lp-1)
c
          call odeint(y,2,r(2),r0,eps,h1,hmin,nok,nbad,ode,rkqc)
c
c         Write solution to file 20
c
c$$$          do 200 j = 1, kmax
c$$$             write(20,8000)  xp(j),(yp(m,j),m=1,2)
c$$$  200     continue
c$$$ 8000     format(1x,f10.6,4x,2(2x,f14.6))
c$$$c 
c$$$          write(20,'(/,1x,a,/)') '---------------------'
c
          ld = y(2) / y(1)
ccccccccc          write(6,9100) energ, y(2), y(1), ld, nok, nbad
          write(3,8000) energ, ld, eigv
 8000     format(1x,f12.6,3x,f12.6,3x,f12.6)
 9100     format(1x,f12.6,3x,2(2x,f10.6),3x,f10.6,2(2x,i4))
c
   30   continue
c
          close(3)
c
  100 continue  
c
      return
c
      end

c
      subroutine ode(x,y,dxdy)
c
      implicit none
c
      include 'radial.h'
      include 'ode_blk.h'
c
c     Returns the derivatives (right hand sides) for the
c     system of differential equations:
c
c     dy1/dx = y2
c
c     dy2/dx = ( V(x) - energy ) y1
c
c     appropriate for the solution of the Schroedinger equation.
c
      double precision x
      double precision y(2), dxdy(2)
      double precision potx, delta_pot
c
c     Interpolate (three points used)
c
      call hunt(r,nr,x,jint)
      call polint(r(jint),v(jint),3,x,potx,delta_pot)
c
      dxdy(1) = y(2)
      dxdy(2) = ( potx - energ ) * y(1)
c
      return
c
      end
c
      subroutine denplot
c
c  Prints the charge density
c     
      include 'radial.h'
      include 'charge.h'
      include 'param.h'
c
      double precision pi
      parameter (pi=3.141592653589d0)
c
      double precision fx, step, delta
      integer j
c
c     Minimum step-size for plotting. This cuts down on file size.
c     Set it to zero to recover old behavior.
c
      parameter (delta = 0.000d0)
ccc      parameter (delta = 0.005d0)
c
c     Specify name according to type of job.
c
      if (job .eq. 4) then
         open(unit=3,file='PTCHARGE',form='formatted',status='unknown')
      else
         open(unit=3,file='AECHARGE',form='formatted',status='unknown')
      endif
c
c     Keep for backwards compatibility
c
      open(unit=66,file='CHARGE',form='formatted',status='unknown')
c
      open(unit=4,file='RHO',form='formatted',status='unknown')
c
c     Write out r, cdu, cdd and cdc
c
      step = 0.0d0
      do 10 j = 2, nr
         if (r(j) .ge. step) then
            write(3,9000) r(j), cdu(j), cdd(j), cdc(j)
            write(66,9000) r(j), cdu(j), cdd(j), cdc(j)
            fx = 1.d0 / ( 4.d0 * pi * r(j)**2)
            write(4,9000) r(j), fx*cdu(j), fx*cdd(j), fx*cdc(j)
            step = step + delta
         endif
   10 continue
 9000 format(1x,f15.10,3x,3f18.5)
c
      close(unit=3)
      close(unit=66)
      close(unit=4)
c
      return
c
      end

c
      subroutine saxpy(n,sa,sx,incx,sy,incy)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOP FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
C
C     .. Scalar Arguments ..
      double precision sa
      integer incx, incy, n
C     ..
C     .. Array Arguments ..
      double precision sx(*), sy(*)
C     ..
C     .. Local Scalars ..
      integer i, ix, iy, m, mp1
C     ..
C     .. Intrinsic Functions ..
      intrinsic mod
C     ..
      if (n .le. 0) return
      if (sa .eq. 0.0D0) return
      if (incx .eq. 1 .and. incy .eq. 1) go to 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      if (incy .lt. 0) iy = (-n+1)*incy + 1
      do 10 i = 1, n
         sy(iy) = sy(iy) + sa*sx(ix)
         ix = ix + incx
         iy = iy + incy
   10 continue
c
      return
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 continue
      m = mod(n,4)
      if (m .eq. 0) go to 40
      do 30 i = 1, m
         sy(i) = sy(i) + sa*sx(i)
   30 continue
      if (n .lt. 4) return
   40 continue
      mp1 = m + 1
      do 50 i = mp1, n, 4
         sy(i) = sy(i) + sa*sx(i)
         sy(i+1) = sy(i+1) + sa*sx(i+1)
         sy(i+2) = sy(i+2) + sa*sx(i+2)
         sy(i+3) = sy(i+3) + sa*sx(i+3)
   50 continue
c
      return
c
      end
c
      double precision function sdot(n,sx,incx,sy,incy)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
C
C     .. Scalar Arguments ..
      integer incx, incy, n
C     ..
C     .. Array Arguments ..
      double precision sx(*), sy(*)
C     ..
C     .. Local Scalars ..
      double precision stemp
      integer i, ix, iy, m, mp1
C     ..
C     .. Intrinsic Functions ..
      intrinsic mod
C     ..
      stemp = 0.0D0
      sdot = 0.0D0
      if (n .le. 0) return
      if (incx .eq. 1 .and. incy .eq. 1) go to 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      if (incy .lt. 0) iy = (-n+1)*incy + 1
      do 10 i = 1, n
         stemp = stemp + sx(ix)*sy(iy)
         ix = ix + incx
         iy = iy + incy
   10 continue
      sdot = stemp
c
      return
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 continue
      m = mod(n,5)
      if (m .eq. 0) go to 40
      do 30 i = 1, m
         stemp = stemp + sx(i)*sy(i)
   30 continue
      if (n .lt. 5) go to 60
   40 continue
      mp1 = m + 1
      do 50 i = mp1, n, 5
         stemp = stemp + sx(i)*sy(i) + sx(i+1)*sy(i+1) +
     &           sx(i+2)*sy(i+2) + sx(i+3)*sy(i+3) + sx(i+4)*sy(i+4)
   50 continue
   60 continue
      sdot = stemp
c
      return
c
      end
c
      subroutine scopy(n,sx,incx,sy,incy)
C
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
C
C     .. Scalar Arguments ..
      integer incx, incy, n
C     ..
C     .. Array Arguments ..
      double precision sx(*), sy(*)
C     ..
C     .. Local Scalars ..
      integer i, ix, iy, m, mp1
C     ..
C     .. Intrinsic Functions ..
      intrinsic mod
C     ..
      if (n .le. 0) return
      if (incx .eq. 1 .and. incy .eq. 1) go to 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      if (incy .lt. 0) iy = (-n+1)*incy + 1
      do 10 i = 1, n
         sy(iy) = sx(ix)
         ix = ix + incx
         iy = iy + incy
   10 continue
c
      return
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 continue
      m = mod(n,7)
      if (m .eq. 0) go to 40
      do 30 i = 1, m
         sy(i) = sx(i)
   30 continue
      if (n .lt. 7) return
   40 continue
      mp1 = m + 1
      do 50 i = mp1, n, 7
         sy(i) = sx(i)
         sy(i+1) = sx(i+1)
         sy(i+2) = sx(i+2)
         sy(i+3) = sx(i+3)
         sy(i+4) = sx(i+4)
         sy(i+5) = sx(i+5)
         sy(i+6) = sx(i+6)
   50 continue
c
      return
c
      end
c
      subroutine sscal(n,sa,sx,incx)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO 1.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
C
C     .. Scalar Arguments ..
      double precision sa
      integer incx, n
C     ..
C     .. Array Arguments ..
      double precision sx(*)
C     ..
C     .. Local Scalars ..
      integer i, m, mp1, nincx
C     ..
C     .. Intrinsic Functions ..
      intrinsic mod
C     ..
      if (n .le. 0) return
      if (incx .eq. 1) go to 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      nincx = n*incx
      do 10 i = 1, nincx, incx
         sx(i) = sa*sx(i)
   10 continue
c
      return
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 continue
      m = mod(n,5)
      if (m .eq. 0) go to 40
      do 30 i = 1, m
         sx(i) = sa*sx(i)
   30 continue
      if (n .lt. 5) return
   40 continue
      mp1 = m + 1
      do 50 i = mp1, n, 5
         sx(i) = sa*sx(i)
         sx(i+1) = sa*sx(i+1)
         sx(i+2) = sa*sx(i+2)
         sx(i+3) = sa*sx(i+3)
         sx(i+4) = sa*sx(i+4)
   50 continue
c
      return
c
      end
      subroutine tridib(n,eps1,d,e,e2,lb,ub,m11,m,w,ind,ierr,rv4,rv5)
c
      integer i,j,k,l,m,n,p,q,r,s,ii,m1,m2,m11,m22,tag,ierr,isturm
      double precision d(n),e(n),e2(n),w(m),rv4(n),rv5(n)
      double precision u,v,lb,t1,t2,ub,xu,x0,x1,eps1,tst1,tst2,epslon
      integer ind(m)
c
c     this subroutine is a translation of the algol procedure bisect,
c     num. math. 9, 386-393(1967) by barth, martin, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 249-256(1971).
c
c     this subroutine finds those eigenvalues of a tridiagonal
c     symmetric matrix between specified boundary indices,
c     using bisection.
c
c     on input
c
c        n is the order of the matrix.
c
c        eps1 is an absolute error tolerance for the computed
c          eigenvalues.  if the input eps1 is non-positive,
c          it is reset for each submatrix to a default value,
c          namely, minus the product of the relative machine
c          precision and the 1-norm of the submatrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is arbitrary.
c
c        m11 specifies the lower boundary index for the desired
c          eigenvalues.
c
c        m specifies the number of eigenvalues desired.  the upper
c          boundary index m22 is then obtained as m22=m11+m-1.
c
c     on output
c
c        eps1 is unaltered unless it has been reset to its
c          (last) default value.
c
c        d and e are unaltered.
c
c        elements of e2, corresponding to elements of e regarded
c          as negligible, have been replaced by zero causing the
c          matrix to split into a direct sum of submatrices.
c          e2(1) is also set to zero.
c
c        lb and ub define an interval containing exactly the desired
c          eigenvalues.
c
c        w contains, in its first m positions, the eigenvalues
c          between indices m11 and m22 in ascending order.
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc..
c
c        ierr is set to
c          zero       for normal return,
c          3*n+1      if multiple eigenvalues at index m11 make
c                     unique selection impossible,
c          3*n+2      if multiple eigenvalues at index m22 make
c                     unique selection impossible.
c
c        rv4 and rv5 are temporary storage arrays.
c
c     note that subroutine tql1, imtql1, or tqlrat is generally faster
c     than tridib, if more than n/4 eigenvalues are to be found.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      tag = 0
      xu = d(1)
      x0 = d(1)
      u = 0.0d0
c     .......... look for small sub-diagonal entries and determine an
c                interval containing all the eigenvalues ..........
      do 40 i = 1, n
         x1 = u
         u = 0.0d0
         if (i .ne. n) u = dabs(e(i+1))
         xu = dmin1(d(i)-(x1+u),xu)
         x0 = dmax1(d(i)+(x1+u),x0)
         if (i .eq. 1) go to 20
         tst1 = dabs(d(i)) + dabs(d(i-1))
         tst2 = tst1 + dabs(e(i))
         if (tst2 .gt. tst1) go to 40
   20    e2(i) = 0.0d0
   40 continue
c
      x1 = n
      x1 = x1 * epslon(dmax1(dabs(xu),dabs(x0)))
      xu = xu - x1
      t1 = xu
      x0 = x0 + x1
      t2 = x0
c     .......... determine an interval containing exactly
c                the desired eigenvalues ..........
      p = 1
      q = n
      m1 = m11 - 1
      if (m1 .eq. 0) go to 75
      isturm = 1
   50 v = x1
      x1 = xu + (x0 - xu) * 0.5d0
      if (x1 .eq. v) go to 980
      go to 320
   60 if (s - m1) 65, 73, 70
   65 xu = x1
      go to 50
   70 x0 = x1
      go to 50
   73 xu = x1
      t1 = x1
   75 m22 = m1 + m
      if (m22 .eq. n) go to 90
      x0 = t2
      isturm = 2
      go to 50
   80 if (s - m22) 65, 85, 70
   85 t2 = x1
   90 q = 0
      r = 0
c     .......... establish and process next submatrix, refining
c                interval by the gerschgorin bounds ..........
  100 if (r .eq. m) go to 1001
      tag = tag + 1
      p = q + 1
      xu = d(p)
      x0 = d(p)
      u = 0.0d0
c
      do 120 q = p, n
         x1 = u
         u = 0.0d0
         v = 0.0d0
         if (q .eq. n) go to 110
         u = dabs(e(q+1))
         v = e2(q+1)
  110    xu = dmin1(d(q)-(x1+u),xu)
         x0 = dmax1(d(q)+(x1+u),x0)
         if (v .eq. 0.0d0) go to 140
  120 continue
c
  140 x1 = epslon(dmax1(dabs(xu),dabs(x0)))
      if (eps1 .le. 0.0d0) eps1 = -x1
      if (p .ne. q) go to 180
c     .......... check for isolated root within interval ..........
      if (t1 .gt. d(p) .or. d(p) .ge. t2) go to 940
      m1 = p
      m2 = p
      rv5(p) = d(p)
      go to 900
  180 x1 = x1 * (q - p + 1)
      lb = dmax1(t1,xu-x1)
      ub = dmin1(t2,x0+x1)
      x1 = lb
      isturm = 3
      go to 320
  200 m1 = s + 1
      x1 = ub
      isturm = 4
      go to 320
  220 m2 = s
      if (m1 .gt. m2) go to 940
c     .......... find roots by bisection ..........
      x0 = ub
      isturm = 5
c
      do 240 i = m1, m2
         rv5(i) = ub
         rv4(i) = lb
  240 continue
c     .......... loop for k-th eigenvalue
c                for k=m2 step -1 until m1 do --
c                (-do- not used to legalize -computed go to-) ..........
      k = m2
  250    xu = lb
c     .......... for i=k step -1 until m1 do -- ..........
         do 260 ii = m1, k
            i = m1 + k - ii
            if (xu .ge. rv4(i)) go to 260
            xu = rv4(i)
            go to 280
  260    continue
c
  280    if (x0 .gt. rv5(k)) x0 = rv5(k)
c     .......... next bisection step ..........
  300    x1 = (xu + x0) * 0.5d0
         if ((x0 - xu) .le. dabs(eps1)) go to 420
         tst1 = 2.0d0 * (dabs(xu) + dabs(x0))
         tst2 = tst1 + (x0 - xu)
         if (tst2 .eq. tst1) go to 420
c     .......... in-line procedure for sturm sequence ..........
  320    s = p - 1
         u = 1.0d0
c
         do 340 i = p, q
            if (u .ne. 0.0d0) go to 325
            v = dabs(e(i)) / epslon(1.0d0)
            if (e2(i) .eq. 0.0d0) v = 0.0d0
            go to 330
  325       v = e2(i) / u
  330       u = d(i) - x1 - v
            if (u .lt. 0.0d0) s = s + 1
  340    continue
c
         go to (60,80,200,220,360), isturm
c     .......... refine intervals ..........
  360    if (s .ge. k) go to 400
         xu = x1
         if (s .ge. m1) go to 380
         rv4(m1) = x1
         go to 300
  380    rv4(s+1) = x1
         if (rv5(s) .gt. x1) rv5(s) = x1
         go to 300
  400    x0 = x1
         go to 300
c     .......... k-th eigenvalue found ..........
  420    rv5(k) = x1
      k = k - 1
      if (k .ge. m1) go to 250
c     .......... order eigenvalues tagged with their
c                submatrix associations ..........
  900 s = r
      r = r + m2 - m1 + 1
      j = 1
      k = m1
c
      do 920 l = 1, r
         if (j .gt. s) go to 910
         if (k .gt. m2) go to 940
         if (rv5(k) .ge. w(l)) go to 915
c
         do 905 ii = j, s
            i = l + s - ii
            w(i+1) = w(i)
            ind(i+1) = ind(i)
  905    continue
c
  910    w(l) = rv5(k)
         ind(l) = tag
         k = k + 1
         go to 920
  915    j = j + 1
  920 continue
c
  940 if (q .lt. n) go to 100
      go to 1001
c     .......... set error -- interval cannot be found containing
c                exactly the desired eigenvalues ..........
  980 ierr = 3 * n + isturm
 1001 lb = t1
      ub = t2
      return
      end
      subroutine tinvit(nm,n,d,e,e2,m,w,ind,z,
     x                  ierr,rv1,rv2,rv3,rv4,rv6)
c
      integer i,j,m,n,p,q,r,s,ii,ip,jj,nm,its,tag,ierr,group
      double precision d(n),e(n),e2(n),w(m),z(nm,m),
     x       rv1(n),rv2(n),rv3(n),rv4(n),rv6(n)
      double precision u,v,uk,xu,x0,x1,eps2,eps3,eps4,norm,order,epslon,
     x       pythag
      integer ind(m)
c
c     this subroutine is a translation of the inverse iteration tech-
c     nique in the algol procedure tristurm by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
c
c     this subroutine finds those eigenvectors of a tridiagonal
c     symmetric matrix corresponding to specified eigenvalues,
c     using inverse iteration.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        e2 contains the squares of the corresponding elements of e,
c          with zeros corresponding to negligible elements of e.
c          e(i) is considered negligible if it is not larger than
c          the product of the relative machine precision and the sum
c          of the magnitudes of d(i) and d(i-1).  e2(1) must contain
c          0.0d0 if the eigenvalues are in ascending order, or 2.0d0
c          if the eigenvalues are in descending order.  if  bisect,
c          tridib, or  imtqlv  has been used to find the eigenvalues,
c          their output e2 array is exactly what is expected here.
c
c        m is the number of specified eigenvalues.
c
c        w contains the m eigenvalues in ascending or descending order.
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc.
c
c     on output
c
c        all input arrays are unaltered.
c
c        z contains the associated set of orthonormal eigenvectors.
c          any vector which fails to converge is set to zero.
c
c        ierr is set to
c          zero       for normal return,
c          -r         if the eigenvector corresponding to the r-th
c                     eigenvalue fails to converge in 5 iterations.
c
c        rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (m .eq. 0) go to 1001
      tag = 0
      order = 1.0d0 - e2(1)
      q = 0
c     .......... establish and process next submatrix ..........
  100 p = q + 1
c
      do 120 q = p, n
         if (q .eq. n) go to 140
         if (e2(q+1) .eq. 0.0d0) go to 140
  120 continue
c     .......... find vectors by inverse iteration ..........
  140 tag = tag + 1
      s = 0
c
      do 920 r = 1, m
         if (ind(r) .ne. tag) go to 920
         its = 1
         x1 = w(r)
         if (s .ne. 0) go to 510
c     .......... check for isolated root ..........
         xu = 1.0d0
         if (p .ne. q) go to 490
         rv6(p) = 1.0d0
         go to 870
  490    norm = dabs(d(p))
         ip = p + 1
c
         do 500 i = ip, q
  500    norm = dmax1(norm, dabs(d(i))+dabs(e(i)))
c     .......... eps2 is the criterion for grouping,
c                eps3 replaces zero pivots and equal
c                roots are modified by eps3,
c                eps4 is taken very small to avoid overflow ..........
         eps2 = 1.0d-3 * norm
         eps3 = epslon(norm)
         uk = q - p + 1
         eps4 = uk * eps3
         uk = eps4 / dsqrt(uk)
         s = p
  505    group = 0
         go to 520
c     .......... look for close or coincident roots ..........
  510    if (dabs(x1-x0) .ge. eps2) go to 505
         group = group + 1
         if (order * (x1 - x0) .le. 0.0d0) x1 = x0 + order * eps3
c     .......... elimination with interchanges and
c                initialization of vector ..........
  520    v = 0.0d0
c
         do 580 i = p, q
            rv6(i) = uk
            if (i .eq. p) go to 560
            if (dabs(e(i)) .lt. dabs(u)) go to 540
c     .......... warning -- a divide check may occur here if
c                e2 array has not been specified correctly ..........
            xu = u / e(i)
            rv4(i) = xu
            rv1(i-1) = e(i)
            rv2(i-1) = d(i) - x1
            rv3(i-1) = 0.0d0
            if (i .ne. q) rv3(i-1) = e(i+1)
            u = v - xu * rv2(i-1)
            v = -xu * rv3(i-1)
            go to 580
  540       xu = e(i) / u
            rv4(i) = xu
            rv1(i-1) = u
            rv2(i-1) = v
            rv3(i-1) = 0.0d0
  560       u = d(i) - x1 - xu * v
            if (i .ne. q) v = e(i+1)
  580    continue
c
         if (u .eq. 0.0d0) u = eps3
         rv1(q) = u
         rv2(q) = 0.0d0
         rv3(q) = 0.0d0
c     .......... back substitution
c                for i=q step -1 until p do -- ..........
  600    do 620 ii = p, q
            i = p + q - ii
            rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
            v = u
            u = rv6(i)
  620    continue
c     .......... orthogonalize with respect to previous
c                members of group ..........
         if (group .eq. 0) go to 700
         j = r
c
         do 680 jj = 1, group
  630       j = j - 1
            if (ind(j) .ne. tag) go to 630
            xu = 0.0d0
c
            do 640 i = p, q
  640       xu = xu + rv6(i) * z(i,j)
c
            do 660 i = p, q
  660       rv6(i) = rv6(i) - xu * z(i,j)
c
  680    continue
c
  700    norm = 0.0d0
c
         do 720 i = p, q
  720    norm = norm + dabs(rv6(i))
c
         if (norm .ge. 1.0d0) go to 840
c     .......... forward substitution ..........
         if (its .eq. 5) go to 830
         if (norm .ne. 0.0d0) go to 740
         rv6(s) = eps4
         s = s + 1
         if (s .gt. q) s = p
         go to 780
  740    xu = eps4 / norm
c
         do 760 i = p, q
  760    rv6(i) = rv6(i) * xu
c     .......... elimination operations on next vector
c                iterate ..........
  780    do 820 i = ip, q
            u = rv6(i)
c     .......... if rv1(i-1) .eq. e(i), a row interchange
c                was performed earlier in the
c                triangularization process ..........
            if (rv1(i-1) .ne. e(i)) go to 800
            u = rv6(i-1)
            rv6(i-1) = rv6(i)
  800       rv6(i) = u - rv4(i) * rv6(i-1)
  820    continue
c
         its = its + 1
         go to 600
c     .......... set error -- non-converged eigenvector ..........
  830    ierr = -r
         xu = 0.0d0
         go to 870
c     .......... normalize so that sum of squares is
c                1 and expand to full order ..........
  840    u = 0.0d0
c
         do 860 i = p, q
  860    u = pythag(u,rv6(i))
c
         xu = 1.0d0 / u
c
  870    do 880 i = 1, n
  880    z(i,r) = 0.0d0
c
         do 900 i = p, q
  900    z(i,r) = rv6(i) * xu
c
         x0 = x1
  920 continue
c
      if (q .lt. n) go to 100
 1001 return
      end
C  LINPACKS.FOR  02 December 1989
C  The single precision version of LINPACK.
C
C***********************************************************************
C
      subroutine sgeco(a,lda,n,ipvt,rcond,z)
C
C     SGECO FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW SGECO BY SGESL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW SGECO BY SGEDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW SGECO BY SGEDI.
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       REAL(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK SGEFA
C     BLAS SAXPY,SDOT,SSCAL,SASUM
C     FORTRAN ABS,AMAX1,SIGN
C
C     INTERNAL VARIABLES
C
C
C
C     COMPUTE 1-NORM OF A
C
C     .. Scalar Arguments ..
      double precision rcond
      integer lda, n
C     ..
C     .. Array Arguments ..
      double precision a(lda,*), z(*)
      integer ipvt(*)
C     ..
C     .. Local Scalars ..
      double precision anorm, ek, s, sm, t, wk, wkm, ynorm
      integer info, j, k, kb, kp1, l
C     ..
C     .. External Functions ..
      double precision sasum, sdot
      external sasum, sdot
C     ..
C     .. External Subroutines ..
      external saxpy, sgefa, sscal
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, dmax1, sign
C     ..
      anorm = 0.0D0
      do 10 j = 1, n
         anorm = dmax1(anorm,sasum(n,a(1,j),1))
   10 continue
C
C     FACTOR
C
      call sgefa(a,lda,n,ipvt,info)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      ek = 1.0D0
      do 20 j = 1, n
         z(j) = 0.0D0
   20 continue
      do 100 k = 1, n
         if (z(k) .ne. 0.0D0) ek = sign(ek,-z(k))
         if (abs(ek-z(k)) .le. abs(a(k,k))) go to 30
         s = abs(a(k,k))/abs(ek-z(k))
         call sscal(n,s,z,1)
         ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = abs(wk)
         sm = abs(wkm)
         if (a(k,k) .eq. 0.0D0) go to 40
         wk = wk/a(k,k)
         wkm = wkm/a(k,k)
c
         go to 50
c
   40    continue
         wk = 1.0D0
         wkm = 1.0D0
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
         do 60 j = kp1, n
            sm = sm + abs(z(j)+wkm*a(k,j))
            z(j) = z(j) + wk*a(k,j)
            s = s + abs(z(j))
   60    continue
         if (s .ge. sm) go to 80
         t = wkm - wk
         wk = wkm
         do 70 j = kp1, n
            z(j) = z(j) + t*a(k,j)
   70    continue
   80    continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0D0/sasum(n,z,1)
      call sscal(n,s,z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + sdot(n-k,a(k+1,k),1,z(k+1),1)
         if (abs(z(k)) .le. 1.0D0) go to 110
         s = 1.0D0/abs(z(k))
         call sscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0D0/sasum(n,z,1)
      call sscal(n,s,z,1)
C
      ynorm = 1.0D0
C
C     SOLVE L*V = Y
C
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call saxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (abs(z(k)) .le. 1.0D0) go to 130
         s = 1.0D0/abs(z(k))
         call sscal(n,s,z,1)
         ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0D0/sasum(n,z,1)
      call sscal(n,s,z,1)
      ynorm = s*ynorm
C
C     SOLVE  U*Z = V
C
      do 160 kb = 1, n
         k = n + 1 - kb
         if (abs(z(k)) .le. abs(a(k,k))) go to 150
         s = abs(a(k,k))/abs(z(k))
         call sscal(n,s,z,1)
         ynorm = s*ynorm
  150    continue
         if (a(k,k) .ne. 0.0D0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.0D0) z(k) = 1.0D0
         t = -z(k)
         call saxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
C     MAKE ZNORM = 1.0
      s = 1.0D0/sasum(n,z,1)
      call sscal(n,s,z,1)
      ynorm = s*ynorm
C
      if (anorm .ne. 0.0D0) rcond = ynorm/anorm
      if (anorm .eq. 0.0D0) rcond = 0.0D0
c
      return
c
      end
      double precision function epslon (x)
      double precision x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing floating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to 
c            the accuracy used in floating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying 
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger floating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      end
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
c
      double precision function sasum(n,sx,incx)
C
C     TAKES THE SUM OF THE ABSOLUTE VALUES.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
C
C     .. Scalar Arguments ..
      integer incx, n
C     ..
C     .. Array Arguments ..
      double precision sx(*)
C     ..
C     .. Local Scalars ..
      double precision stemp
      integer i, m, mp1, nincx
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, mod
C     ..
      sasum = 0.0D0
      stemp = 0.0D0
      if (n .le. 0) return
      if (incx .eq. 1) go to 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      nincx = n*incx
      do 10 i = 1, nincx, incx
         stemp = stemp + abs(sx(i))
   10 continue
      sasum = stemp
c
      return
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 continue
      m = mod(n,6)
      if (m .eq. 0) go to 40
      do 30 i = 1, m
         stemp = stemp + abs(sx(i))
   30 continue
      if (n .lt. 6) go to 60
   40 continue
      mp1 = m + 1
      do 50 i = mp1, n, 6
         stemp = stemp + abs(sx(i)) + abs(sx(i+1)) + abs(sx(i+2)) +
     &           abs(sx(i+3)) + abs(sx(i+4)) + abs(sx(i+5))
   50 continue
   60 continue
      sasum = stemp
c
      return
c
      end
c
      subroutine sgefa(a,lda,n,ipvt,info)
C
C     SGEFA FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.
C
C     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SAXPY,SSCAL,ISAMAX
C
C     INTERNAL VARIABLES
C
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
C     .. Scalar Arguments ..
      integer info, lda, n
C     ..
C     .. Array Arguments ..
      double precision a(lda,*)
      integer ipvt(*)
C     ..
C     .. Local Scalars ..
      double precision t
      integer j, k, kp1, l, nm1
C     ..
C     .. External Functions ..
      integer isamax
      external isamax
C     ..
C     .. External Subroutines ..
      external saxpy, sscal
C     ..
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
C
C        FIND L = PIVOT INDEX
C
         l = isamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         if (a(l,k) .eq. 0.0D0) go to 40
C
C           INTERCHANGE IF NECESSARY
C
         if (l .eq. k) go to 10
         t = a(l,k)
         a(l,k) = a(k,k)
         a(k,k) = t
   10    continue
C
C           COMPUTE MULTIPLIERS
C
         t = -1.0D0/a(k,k)
         call sscal(n-k,t,a(k+1,k),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
         do 30 j = kp1, n
            t = a(l,j)
            if (l .eq. k) go to 20
            a(l,j) = a(k,j)
            a(k,j) = t
   20       continue
            call saxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30    continue
c
         go to 50
c
   40    continue
         info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0D0) info = n
c
      return
c
      end
c
      integer function isamax(n,sx,incx)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
C
C     .. Scalar Arguments ..
      integer incx, n
C     ..
C     .. Array Arguments ..
      double precision sx(*)
C     ..
C     .. Local Scalars ..
      double precision smax
      integer i, ix
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs
C     ..
      isamax = 0
      if (n .lt. 1) return
      isamax = 1
      if (n .eq. 1) return
      if (incx .eq. 1) go to 30
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      ix = 1
      smax = abs(sx(1))
      ix = ix + incx
      do 20 i = 2, n
         if (abs(sx(ix)) .le. smax) go to 10
         isamax = i
         smax = abs(sx(ix))
   10    continue
         ix = ix + incx
   20 continue
c
      return
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   30 continue
      smax = abs(sx(1))
      do 40 i = 2, n
         if (abs(sx(i)) .le. smax) go to 40
         isamax = i
         smax = abs(sx(i))
   40 continue
c
      return
c
      end
c
      subroutine sgesl(a,lda,n,ipvt,b,job)
C
C     SGESL SOLVES THE REAL SYSTEM
C     A * X = B  OR  TRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY SGECO OR SGEFA.
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE OUTPUT FROM SGECO OR SGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SGECO OR SGEFA.
C
C        B       REAL(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
C        OR SGEFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SAXPY,SDOT
C
C     INTERNAL VARIABLES
C
C
C     .. Scalar Arguments ..
      integer job, lda, n
C     ..
C     .. Array Arguments ..
      double precision a(lda,*), b(*)
      integer ipvt(*)
C     ..
C     .. Local Scalars ..
      double precision t
      integer k, kb, l, nm1
C     ..
C     .. External Functions ..
      double precision sdot
      external sdot
C     ..
C     .. External Subroutines ..
      external saxpy
C     ..
      nm1 = n - 1
      if (job .ne. 0) go to 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
      if (nm1 .lt. 1) go to 30
      do 20 k = 1, nm1
         l = ipvt(k)
         t = b(l)
         if (l .eq. k) go to 10
         b(l) = b(k)
         b(k) = t
   10    continue
         call saxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20 continue
   30 continue
C
C        NOW SOLVE  U*X = Y
C
      do 40 kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/a(k,k)
         t = -b(k)
         call saxpy(k-1,t,a(1,k),1,b(1),1)
   40 continue
c
      go to 100
c
   50 continue
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
      do 60 k = 1, n
         t = sdot(k-1,a(1,k),1,b(1),1)
         b(k) = (b(k)-t)/a(k,k)
   60 continue
C
C        NOW SOLVE TRANS(L)*X = Y
C
      if (nm1 .lt. 1) go to 90
      do 80 kb = 1, nm1
         k = n - kb
         b(k) = b(k) + sdot(n-k,a(k+1,k),1,b(k+1),1)
         l = ipvt(k)
         if (l .eq. k) go to 70
         t = b(l)
         b(l) = b(k)
         b(k) = t
   70    continue
   80 continue
   90 continue
  100 continue
c
      return
c
      end
c
      integer function isrchfgt(n,array,inc,target)
C
C     Returns the location of the first element in a real ARRAY
C          that is greater than a real TARGET.  Returns N+1 if
C          TARGET is not found and 0 if N < 1.
C     Martin J. McBride.  7/17/85.
C     General Electric CRD, Information System Operation.
C
c
C     .. Scalar Arguments ..
      double precision target
      integer inc, n
C     ..
C     .. Array Arguments ..
      double precision array(*)
C     ..
C     .. Local Scalars ..
      integer i, ix
C     ..
      isrchfgt = 0
      if (n .lt. 1) return
      ix = 1
      if (inc .lt. 0) ix = (-inc)*(n-1) + 1
      do 10 i = 1, n
         if (array(ix) .gt. target) go to 20
         ix = ix + inc
   10 continue
   20 continue
      isrchfgt = i
c
      return
c
      end
      SUBROUTINE HUNT(XX,N,X,JLO)
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
      INTEGER JLO,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION XX(N)
C     ..
C     .. Local Scalars ..
      INTEGER INC,JHI,JM
      LOGICAL ASCND
C     ..
      ASCND = XX(N) .GT. XX(1)
      IF (JLO.LE.0 .OR. JLO.GT.N) THEN
          JLO = 0
          JHI = N + 1
          GO TO 3

      END IF

      INC = 1
      IF (X.GE.XX(JLO) .EQV. ASCND) THEN
    1     JHI = JLO + INC
          IF (JHI.GT.N) THEN
              JHI = N + 1

          ELSE IF (X.GE.XX(JHI) .EQV. ASCND) THEN
              JLO = JHI
              INC = INC + INC
              GO TO 1

          END IF

      ELSE
          JHI = JLO
    2     JLO = JHI - INC
          IF (JLO.LT.1) THEN
              JLO = 0

          ELSE IF (X.LT.XX(JLO) .EQV. ASCND) THEN
              JHI = JLO
              INC = INC + INC
              GO TO 2

          END IF

      END IF

    3 IF (JHI-JLO.EQ.1) RETURN
      JM = (JHI+JLO)/2
      IF (X.GT.XX(JM) .EQV. ASCND) THEN
          JLO = JM

      ELSE
          JHI = JM
      END IF

      GO TO 3

      END
      SUBROUTINE LOCATE(XX,N,X,J)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
      INTEGER J,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION XX(N)
C     ..
C     .. Local Scalars ..
      INTEGER JL,JM,JU
C     ..
      JL = 0
      JU = N + 1
   10 IF (JU-JL.GT.1) THEN
          JM = (JU+JL)/2
          IF ((XX(N).GT.XX(1)) .EQV. (X.GT.XX(JM))) THEN
              JL = JM

          ELSE
              JU = JM
          END IF

          GO TO 10

      END IF

      J = JL
      RETURN

      END
      SUBROUTINE ODEINT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS,
     +                  RKQC)
C
      include 'ode_path.h'
C     ..
C     .. Parameters ..
      INTEGER MAXSTP
      PARAMETER (MAXSTP=10000)
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION TINY
      PARAMETER (TINY=1.D-30)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS,H1,HMIN,X1,X2
      INTEGER NBAD,NOK,NVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION YSTART(NVAR)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL DERIVS,RKQC
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION H,HDID,HNEXT,X,XSAV
      INTEGER I,NSTP
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DYDX(NMAX),Y(NMAX),YSCAL(NMAX)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN
C     ..
      X = X1
      H = SIGN(H1,X2-X1)
      NOK = 0
      NBAD = 0
      KOUNT = 0
      DO 11 I = 1,NVAR
          Y(I) = YSTART(I)
   11 CONTINUE
      XSAV = X - DXSAV*TWO
      DO 16 NSTP = 1,MAXSTP
          CALL DERIVS(X,Y,DYDX)
          DO 12 I = 1,NVAR
              YSCAL(I) = ABS(Y(I)) + ABS(H*DYDX(I)) + TINY
   12     CONTINUE
          IF (KMAX.GT.0) THEN
              IF (ABS(X-XSAV).GT.ABS(DXSAV)) THEN
                  IF (KOUNT.LT.KMAX-1) THEN
                      KOUNT = KOUNT + 1
                      XP(KOUNT) = X
                      DO 13 I = 1,NVAR
                          YP(I,KOUNT) = Y(I)
   13                 CONTINUE
                      XSAV = X
                  END IF

              END IF

          END IF

          IF ((X+H-X2)* (X+H-X1).GT.ZERO) H = X2 - X
          CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
          IF (HDID.EQ.H) THEN
              NOK = NOK + 1

          ELSE
              NBAD = NBAD + 1
          END IF

          IF ((X-X2)* (X2-X1).GE.ZERO) THEN
              DO 14 I = 1,NVAR
                  YSTART(I) = Y(I)
   14         CONTINUE
              IF (KMAX.NE.0) THEN
                  KOUNT = KOUNT + 1
                  XP(KOUNT) = X
                  DO 15 I = 1,NVAR
                      YP(I,KOUNT) = Y(I)
   15             CONTINUE
              END IF

              RETURN

          END IF

          IF (ABS(HNEXT).LT.HMIN) THEN
              WRITE (*,FMT=*) 'ODEINT - Stepsize smaller than minimum.'
              RETURN

          END IF

          H = HNEXT
   16 CONTINUE
      WRITE (*,FMT=*) 'ODEINT - Too many steps.'
      RETURN

      END
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
C
C     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=10)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION DY,X,Y
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION XA(N),YA(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DEN,DIF,DIFT,HO,HP,W
      INTEGER I,M,NS
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION C(NMAX),D(NMAX)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      NS = 1
      DIF = ABS(X-XA(1))
      DO 11 I = 1,N
          DIFT = ABS(X-XA(I))
          IF (DIFT.LT.DIF) THEN
              NS = I
              DIF = DIFT
          END IF

          C(I) = YA(I)
          D(I) = YA(I)
   11 CONTINUE
      Y = YA(NS)
      NS = NS - 1
      DO 13 M = 1,N - 1
          DO 12 I = 1,N - M
              HO = XA(I) - X
              HP = XA(I+M) - X
              W = C(I+1) - D(I)
              DEN = HO - HP
              IF (DEN.EQ.0.D0) THEN
                  WRITE (*,FMT=*) 'POLINT - DEN=0.0'
                  RETURN

              END IF

              DEN = W/DEN
              D(I) = HP*DEN
              C(I) = HO*DEN
   12     CONTINUE
          IF (2*NS.LT.N-M) THEN
              DY = C(NS+1)

          ELSE
              DY = D(NS)
              NS = NS - 1
          END IF

          Y = Y + DY
   13 CONTINUE
      RETURN

      END
      SUBROUTINE RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
C
C     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=10)
      DOUBLE PRECISION FCOR
      PARAMETER (FCOR=.0666666667D0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.D0)
      DOUBLE PRECISION SAFETY
      PARAMETER (SAFETY=0.9D0)
      DOUBLE PRECISION ERRCON
      PARAMETER (ERRCON=6.D-4)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS,HDID,HNEXT,HTRY,X
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DYDX(NMAX),Y(NMAX),YSCAL(NMAX)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL DERIVS
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ERRMAX,H,HH,PGROW,PSHRNK,XSAV
      INTEGER I
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DYSAV(NMAX),YSAV(NMAX),YTEMP(NMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL RK4
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX
C     ..
      PGROW = -0.20D0
      PSHRNK = -0.25D0
      XSAV = X
      DO 11 I = 1,N
          YSAV(I) = Y(I)
          DYSAV(I) = DYDX(I)
   11 CONTINUE
      H = HTRY
    1 HH = 0.5D0*H
      CALL RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X = XSAV + HH
      CALL DERIVS(X,YTEMP,DYDX)
      CALL RK4(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X = XSAV + H
      IF (X.EQ.XSAV) THEN
          WRITE (*,FMT=*) 'Stepsize not significant in RKQC.'
          RETURN

      END IF

      CALL RK4(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX = 0.D0
      DO 12 I = 1,N
          YTEMP(I) = Y(I) - YTEMP(I)
          ERRMAX = MAX(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
   12 CONTINUE
      ERRMAX = ERRMAX/EPS
      IF (ERRMAX.GT.ONE) THEN
          H = SAFETY*H* (ERRMAX**PSHRNK)
          GO TO 1

      ELSE
          HDID = H
          IF (ERRMAX.GT.ERRCON) THEN
              HNEXT = SAFETY*H* (ERRMAX**PGROW)

          ELSE
              HNEXT = 4.D0*H
          END IF

      END IF

      DO 13 I = 1,N
          Y(I) = Y(I) + YTEMP(I)*FCOR
   13 CONTINUE
      RETURN

      END
      SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
C
C
C     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=10)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION H,X
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DYDX(N),Y(N),YOUT(N)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL DERIVS
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION H6,HH,XH
      INTEGER I
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DYM(NMAX),DYT(NMAX),YT(NMAX)
C     ..
      HH = H*0.5D0
      H6 = H/6.D0
      XH = X + HH
      DO 11 I = 1,N
          YT(I) = Y(I) + HH*DYDX(I)
   11 CONTINUE
      CALL DERIVS(XH,YT,DYT)
      DO 12 I = 1,N
          YT(I) = Y(I) + HH*DYT(I)
   12 CONTINUE
      CALL DERIVS(XH,YT,DYM)
      DO 13 I = 1,N
          YT(I) = Y(I) + H*DYM(I)
          DYM(I) = DYT(I) + DYM(I)
   13 CONTINUE
      CALL DERIVS(X+H,YT,DYT)
      DO 14 I = 1,N
          YOUT(I) = Y(I) + H6* (DYDX(I)+DYT(I)+2.D0*DYM(I))
   14 CONTINUE
      RETURN

      END
c
      subroutine compat_params(str)
c
c     Some internal parameters tend to change over time... This routine
c     assigns values to them based on a set of directives, specified
c     at the top of the input file as:
c     
c     %define COMPAT_UCB
c     %define NEW_CC
c  
c

c       
c     COMPAT_UCB :  Revert to the standard circa 1990 UCB values. Note
c                that these correspond to the first released version
c                of Jose Luis Martins code, not to the old Froyen
c                version (that would be 'froyen' --to be implemented)
c     (The default is: to use a denser grid up to larger radii.
c                Use a larger value for the ps's ecuts.
c                Use the Soler-Balbas XC package)
c
c     NEW_CC     : New core-correction scheme
c     OLD_CC     : Old core-correction scheme  (see wrapup.f)
c
c     The default is to use the new CC scheme only for GGA calculations.
c
c
c     For compatibility with an interim scheme using strings in the
c     input file, this routine accepts an argument "str". If not
c     empty, the user is warned that those strings now carry no
c     weight.
c

      character*(*) str

      include 'compat.h'

      logical leqi, defined
      external leqi, defined
c
      if (str .ne. " ") then
         write(6,'(a,a20)')
     $        '** WARNING: Compatibility string obsolete: ',
     $        str
         stop 'COMPAT'
      endif

      if (defined('COMPAT_UCB')) then

         write(6,'(a)') '*** UCB compatibility mode ***'
         aa_def = 6.d0
         bb_def = 40.d0
         rmax_def = 80.d0
         ecuts = 1.2d-4
         use_excorr = .true.
c
      else

         aa_def = 6.d0
         bb_def = 80.d0
         rmax_def = 120.d0 
         ecuts = 1.d-3
         use_excorr = .false.

      endif
c
c     Flag to use the old excorr subroutine.
c
      if (defined('USE_OLD_EXCORR')) use_excorr = .true.
c
c     Avoid cutting off ionic unscreened pseudopotentials
c
      if (defined('NO_PS_CUTOFFS')) ecuts = 0.d0
c
      use_old_cc = defined('OLD_CC')
      use_new_cc = defined('NEW_CC')

      end

      SUBROUTINE PCC_EXP(NR,ICORE,AC,BC,CC,R,CDC)

c mmga
c mmga   M.M.G. Alemany, January 2000  
c mmga
c
c     constructs the partial core correction. Core charge,
c     first and second derivatives are conserved at the point
c     icore. Interpolation is done with Lagrange formula.
c     The core function is exp(ac+bc*r**2+cc*r**4).

      IMPLICIT none


      integer icore, nr
      double precision ac, bc, cc
      double precision R(NR),CDC(NR)

      integer nn
      PARAMETER (NN = 5)
      double precision DR_K(-NN:NN),DDR_K(-NN:NN)
      double precision CDC_SCA(-NN:NN)

      integer i, in1, in2, in, jn, kn, ln
      double precision f1, f2, dr, cdcp, ddr, ddcdc, cdcpp

      IF ( (ICORE - NN).LT.1 .OR. (ICORE + NN).GT.NR ) THEN
        WRITE(6,2017)
        CALL EXT(830)
      ENDIF

      DO 2000   I = ICORE-NN,ICORE+NN
        CDC_SCA(I-ICORE)= LOG(CDC(I)) - 2.D0 * LOG(R(I))
 2000 CONTINUE

      IN1 =-NN
      IN2 = NN

      DO 2011 IN = IN1,IN2
        IF (IN.EQ.0) THEN
          DR_K(IN) = 0.D0
          DO 2012 JN = IN1,IN2
            IF (JN.NE.0) DR_K(IN) = DR_K(IN) + 1.D0/(0 - JN)
 2012     CONTINUE
        ELSE
          F1 = 1.D0
          F2 = 1.D0
          DO 2014 JN = IN1,IN2
            IF (JN.NE.IN .AND. JN.NE.0) F1 = F1 * (0  - JN)
            IF (JN.NE.IN)               F2 = F2 * (IN - JN)
 2014     CONTINUE
          DR_K(IN) = F1 / F2
        ENDIF
 2011 CONTINUE

      DR = 0.D0
      DO 2016 IN = IN1,IN2
        DR = DR + DR_K(IN) * R(ICORE+IN)
 2016 CONTINUE

      CDCP = 0.D0
      DO 2001 IN = IN1,IN2
        CDCP = CDCP + DR_K(IN) * CDC_SCA(IN)
 2001 CONTINUE
      CDCP = CDCP / DR

      DO 2002 IN = IN1,IN2
        DDR_K(IN)= 0.D0
        IF (IN.EQ.0) THEN
          DO 2003 JN = IN1,IN2
            IF (JN.NE.0) THEN
              DO 2004 KN = JN+1,IN2
                IF (KN.NE.0) DDR_K(IN) = DDR_K(IN) + 
     .                                   2.D0/((0 - JN)*(0 - KN))
 2004         CONTINUE
            ENDIF
 2003     CONTINUE
        ELSE
          F2 = 1.D0
          DO 2005 JN = IN1,IN2
            IF (JN.NE.IN) THEN
              F2 = F2 * (IN - JN)
              DO 2006 KN = JN+1,IN2
                IF (KN.NE.IN) THEN
                  F1 = 2.D0
                  DO 2007 LN = IN1,IN2
                    IF (LN.NE.IN .AND. LN.NE.JN .AND. LN.NE.KN)
     .                  F1 = F1 * (0 - LN)
 2007             CONTINUE
                  DDR_K(IN) = DDR_K(IN) + F1                    
                ENDIF
 2006         CONTINUE
            ENDIF
 2005     CONTINUE
          DDR_K(IN) = DDR_K(IN) / F2
        ENDIF
 2002 CONTINUE

      DDR   = 0.D0
      DDCDC = 0.D0
      DO 2010 IN = IN1,IN2
        DDR   = DDR   + DDR_K(IN) * R(ICORE+IN)
        DDCDC = DDCDC + DDR_K(IN) * CDC_SCA(IN)
 2010 CONTINUE

      CDCPP = (DDCDC - DDR * CDCP) / DR**2



      CC = R(ICORE) * CDCPP - CDCP
      CC = CC / (8.D0 * R(ICORE)**3)

      BC = CDCP - 4.D0 * CC * R(ICORE)**3
      BC = BC / (2.D0 * R(ICORE))

      AC = CDC_SCA(0) - BC*R(ICORE)**2 - CC*R(ICORE)**4


      DO 2009 I = 1,ICORE
          CDC(I)= R(I)*R(I) * 
     .            EXP( (CC*R(I)*R(I) + BC) * R(I)*R(I) + AC ) 
 2009 CONTINUE


 2017 FORMAT(//,' error in pcc_exp - ',/,
     . ' derivatives at r(icore) not calculated')


      RETURN

      END
C
      subroutine change_valence
c
c     Generates a new pseudopotential file
c     with a modified valence charge in it.
c     This might be useful ...
c
      include 'radial.h'
      include 'param.h'
      include 'ion.h'
      include 'charge.h'
c
C     .. Local Scalars ..
      double precision zion
      integer i, j, loi, npotd, npotu, nrm
      character icorrt*2, namet*2
C     ..
C     .. Local Arrays ..
      character ray(6)*10, title(7)*10
C     ..


      open(unit=1,file='VPSIN',status='old',form='unformatted')
      open(unit=2,file='VPS_NEW_VALENCE',
     $        status='unknown',form='unformatted')
      rewind 1
      rewind 2
      read(1) namet, icorrt, irel, nicore, (ray(i),i=1,6),
     &     (title(i),i=1,7), npotd, npotu, nrm, a, b, zion
c
      read(1) (r(i),i=2,nr)
      ray(6) = 'NEWVALENCE'
      write(2) namet, icorrt, irel, nicore, (ray(i),i=1,6),
     &     (title(i),i=1,7), npotd, npotu, nrm, a, b, zion
c
      write(2) (r(i),i=2,nr)
c
c   down potentials (or average relativistic potentials)
c
      do 40 i = 1, npotd
         read(1) loi, (viod(loi+1,j),j=2,nr)
         write(2) loi, (viod(loi+1,j),j=2,nr)
 40   continue
c
c   up potentials (or spin orbit potentials)
c
      if (npotu .gt. 0) then
         do 110 i = 1, npotu
            read(1) loi, (viou(loi+1,j),j=2,nr)
            write(2) loi, (viou(loi+1,j),j=2,nr)
  110    continue
      endif
c     
c  core and valence charges
c
      read(1) (cdc(i),i=2,nr)
      write(2) (cdc(i),i=2,nr)

c     Write the actual valence charge we
c     have in the program's structures
c     Without the infamous "zratio" ...
c
      write(2) (cdd(i)+cdu(i),i=2,nr)

      close(1)
      close(2)

      end






c
      subroutine coreq
c
c  Compute and print the fourier transform of the pseudocore charge
c
c  Note: It is assumed that the function is reasonably smooth... otherwise
c        the Gauss-Legendre quadrature will not work!
c
      include 'radial.h'
      include 'charge.h'
c
      double precision pi
      parameter (pi=3.141592653589d0)
      double precision tol
      parameter ( tol = 1.d-6 )
c
c     Number of abscissae
c
      integer ngauss
      parameter (ngauss=100)
      double precision cored(ngauss), t(ngauss), w(ngauss)
c
c     How many points at what spacing? Should go up to q=20, say
c
      double precision delql
      parameter (delql = 0.05d0)
      integer nql
      parameter (nql = 400)

      integer j, nrpnew
      double precision rmin, rmax, zc, q, dcq

      double precision divdif
      external divdif
c
c     Fourier transform the core charge
c
c     find smallest r above which 4pi*r^2*core charge is negligible
c
      do j = nr, 1, -1
         if (abs(cdc(j)) .gt. tol) go to 130
      enddo
 130  continue
      nrpnew = j + 1
      if (nrpnew .gt. nr) nrpnew = nr
c
      rmin = r(1)
      rmax = r(nrpnew)
c
c     Generate abscissas and weights for Gauss-Legendre integration
c
      call gauleg(rmin,rmax,t,w,ngauss)
c
c     Total core charge    (q=0)
c
c     We evaluate the function to be integrated at the calculated
c     abscissas (and multiply by the weight for convenience)
c     The Cern library function divdif ( interpolating to third
c     order ) is used.
c
      norder = 3

      zc = 0.d0
      do  j = 1, ngauss
         cored(j) = w(j)*divdif(cdc,r,nr,t(j),norder)
         zc = zc + cored(j)
      enddo
      write(6,'(a,f8.4)') 'Total pseudocore charge: ', zc

      call get_unit(iu)
      open(iu,file='COREQ',form='formatted',status='unknown')
      rewind(iu)

      write(iu,'(f8.2,4x,f12.4)') 0.d0, zc
c
c      Rest of the Fourier transform
c
      do  k = 1, nql
         q = delql*k
         dcq = 0.d0
         do j = 1, ngauss
            dcq = dcq + cored(j)*sin(q*t(j))/t(j)
         enddo
         dcq = dcq/q
         write(iu,'(f8.2,4x,f12.4)') q, dcq
      enddo
      close(iu)

      end


      subroutine get_unit(lun)

C     Get an available Fortran unit number

      integer lun

      integer i
      logical unit_used

      do i = 10, 99
         lun = i
         inquire(lun,opened=unit_used)
         if (.not. unit_used) return
      enddo
      stop 'NO LUNS'
      end
c
      subroutine gauleg(x1,x2,x,w,n)
C
C $Id: gauleg.f,v 1.1 2000/02/10 17:56:11 wdpgaara Exp $
C
C $Log: gauleg.f,v $
C Revision 1.1  2000/02/10 17:56:11  wdpgaara
C Implement Fourier transform of core charge.
C
C Revision 1.1.1.1  1997/01/07 08:37:19  wdpgaara
C PS fourier transform package
C
c Revision 1.3  1991/12/13  23:32:01  alberto
c More twiddling
c
c Revision 1.2  1991/12/13  23:12:50  alberto
c Cosmetic changes only
c
c Revision 1.1  1991/12/13  22:56:11  alberto
c Initial revision
c
C     Taken from Numerical Recipes (W.H. Press et al., 1986, p.125)
C
c     Given the lower and upper limits of integration x1 and x2, 
c     and given n, this routine returns arrays x and w of length
c     n, containing the abscissas and weights of the Gauss-Legendre
c     quadrature formula.
c
C     .. Parameters ..
      double precision eps
      parameter (eps=3.D-14)
C     ..
C     .. Scalar Arguments ..
      double precision x1, x2
      integer n
C     ..
C     .. Array Arguments ..
      double precision w(n), x(n)
C     ..
C     .. Local Scalars ..
      double precision p1, p2, p3, pp, xl, xm, z, z1
      integer i, j, m
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, cos
C     ..
c
c     The roots are symmetric in the interval, so we only have to
c     find half of them.
c
      m = (n+1)/2
      xm = 0.5D0*(x2+x1)
      xl = 0.5D0*(x2-x1)
c
c     Loop over the desired roots
c
      do 30 i = 1, m
         z = cos(3.141592654D0*(i-.25D0)/(n+.5D0))
c
c        Starting with the above approximation to the Ith root, we 
c        enter the main loop of refinement by Newton's method
c
   10    continue
         p1 = 1.0D0
         p2 = 0.0D0
c
c        Look up the recurrence relation to get the Legendre 
c        polynomial evaluated at z
c
         do 20 j = 1, n
            p3 = p2
            p2 = p1
            p1 = ((2.0D0*j-1.0D0)*z*p2-(j-1.0D0)*p3)/j
   20    continue
c
c        P1 is now the desired Legendre polynomial. We next compute
c        pp, its derivative, by a standard relation involving also p2,
c        the polynomial of one lower order.
c
         pp = n*(z*p1-p2)/(z*z-1.0D0)
         z1 = z
c
c        Newton's method
c
         z = z1 - p1/pp
         if (abs(z-z1) .gt. eps) go to 10
c
c        Scale the root to the desired interval and put in its
c        symmetric counterpart
c
         x(i) = xm - xl*z
         x(n+1-i) = xm + xl*z
c
c       Compute the weight and its symmetric counterpart
c
         w(i) = 2.0D0*xl/((1.0D0-z*z)*pp*pp)
         w(n+1-i) = w(i)
c
   30 continue
c
      return
c
      end
C $Id: divdif.f,v 1.2 2002/09/26 14:32:16 wdpgaara Exp $
c $Log: divdif.f,v $
c Revision 1.2  2002/09/26 14:32:16  wdpgaara
c Replace print by write
c
c Revision 1.1  2000/02/10 17:56:11  wdpgaara
c Implement Fourier transform of core charge.
c
c Revision 1.1.1.1  1997/01/07 08:37:19  wdpgaara
c PS fourier transform package
c
c Revision 1.3  1991/12/13  23:45:02  alberto
c More cosmetic changes
c
c Revision 1.2  1991/12/13  23:41:49  alberto
c Cosmetic changes
c
      double precision function divdif(f,x,n,z,m)
c
      implicit double precision (a-h,o-z)
c
C     CERN LIBRARY PROGRAM NO E-105.
C     REVISED VERSION JULY 1973.
C     PURPOSE = TO INTERPOLATE IN TABLE OF GIVEN FUNCTION VALUES WHICH
C               ARE STORED AFTER INCREASING OR DECREASING VALUES OF THE
C               ARGUMENTS.NEWTONS GENERAL INTERPOLATION FORMULA IS USED.
C     PARAMETERS ( IN LIST ).
C     F       = THE ARRAY OF THE GIVEN FUNCTION VALUES.F(K)=F(X(,)).
C     X       = THE ARRAY OF GIVEN ARGUMENTS.
C     N       = DIMENSION OF THE ARRAYS F AND X,I.D.THE NUMBER OF POINTS
C               TABLE.
C     Z       = ARGUMENT FOR WHICH THE INTERPOLATION IS WANTED.
C     M       = ORDER OF INTERPOLATION.
C
C     PARAMETERS ( IN COMMON BLOCK / DIVCOF / ).
C     MM      = THE NUMBER OF ELEMENTS STORED IN THE FOLLOWING ARRAYS
C               (MM=M+1).
C     ARG     = AN ARRAY USED FOR STORING THE ARGUMENTS USED IN THE IN-
C               TERPOLATION.
C     VAL     = AN ARRAY USED FOR STORING THE FUNCTION VALUES USED IN
C               THE INTERPOLATION.
C     COF     = AN ARRAY USED FOR STORING THE COEFFICIENTS IN NEWTONS
C               INTERPOLATION FORMULA.
C
      integer n, m
      double precision z, zero
      double precision f(n), x(n)
      integer mm, nmax
c
      common /divcof/ arg(11), val(11), cof(11)
      common /divint/ mm
      data zero, mmax/0.0D0, 10/
c
C
C     INTERNAL PARAMETER.
C     MMAX    = THE MAXIMUM ORDER OF INTERPOLATION PERMITTED.THE DIMEN-
C               SIONS OF THE ARRAYS ARG , VAL AND COF IN THE COMMON
C               BLOCK / DIVCOF / SHOULD BE MMAX+1.
C
   10 continue
      if ((z-x(1))*(x(n)-z) .ge. zero) go to 30
C     Z-VALUE OUTSIDE RANGE,PRINT ERROR MESSAGE.
   20 continue
      write(*,9000) z
      divdif = zero
c
      return
C
   30 continue
      if ((m.le.(n-1)) .and. (m.le.mmax)) go to 40
      mm = m
      if (m .gt. (n-1)) m = n - 1
      if (m .gt. mmax) m = mmax
C     REQUIRED ORDER OF INTERPOLATION TOO HIGH.PRINT ERROR MESSAGE AND
C     REDUCE ORDER.
      write(*,9010) mm, m
C
C     START ACTUAL CALCULATION.
C     COMPUTE POINTER,IPOINT,FOR THE LEFT BOUNDARY OF THE INTERVAL IN
C     WHICH WE HAVE Z I.D. Z IN THE INTERVAL X(IPOINT),X(IPOINT+1).
   40 continue
      cof1 = z - x(1)
      do 50 i = 2, n
         ipoint = i - 1
         cof2 = z - x(i)
         if (cof1*cof2 .le. zero) go to 60
         cof1 = cof2
   50 continue
C     CONSTRUCT TABLE TO BE USED IN THE INTERPOLATION.
   60 continue
      il = ipoint
      iu = il + 1
      jl = 1
      ju = 1
      mm = m + 1
      do 80 i = 1, mm
         i1 = 1
         i2 = 1
         if ((jl.eq.0) .or. (ju.eq.0)) go to 70
         cof1 = dabs(z-x(il))
         cof2 = dabs(x(iu)-z)
         if (cof1 .gt. cof2) i1 = 0
         if (i1 .eq. 1) i2 = 0
   70    continue
         if ((jl.eq.0) .or. (i1.eq.0)) ii = iu
         if ((ju.eq.0) .or. (i2.eq.0)) ii = il
         arg(i) = x(ii)
         cof(i) = f(ii)
         val(i) = f(ii)
         if ((jl.eq.1) .and. (i1.eq.1)) il = il - 1
         if ((ju.eq.1) .and. (i2.eq.1)) iu = iu + 1
         if (il .lt. 1) jl = 0
         if (iu .gt. n) ju = 0
   80 continue
C
      do 100 i = 1, m
         do 90 j = i, m
            index = m + 1 + i - j
            jndex = index - i
            cof(index) = (cof(index)-cof(index-1))/
     1                   (arg(index)-arg(jndex))
   90    continue
  100 continue
C
      sum = cof(m+1)
      do 110 i = 1, m
         index = m + 1 - i
         sum = (z-arg(index))*sum + cof(index)
  110 continue
C
      divdif = sum
c
      return
C
 9000 format(//5x,'*** ERROR MESSAGE FUNCTION DIVDIF *** , ARGUMENT Z ='
     1      ,e21.14,' OUTSIDE RANGE.',//)
 9010 format(//5x,
     1'*** ERROR MESSAGE FUNCTION DIVDIF *** , ORDER OF INTERPOLATION M
     2=',i3,' TOO HIGH,M IS REDUCED TO ',i2,//)
C     OPPO
      end
c
c     Copyright (c) 1993, 1995, 1997 
c     Alberto Garcia, wdpgaara@lg.ehu.es
c
c     Redistribution and use, with or without modification, are 
c     permitted provided that the above copyright notice is retained.
c
c Symbols ------
c
c     This package filters the input stream, detecting and acting upon
c     directives of the form
c
c     %define NAME
c     %delete NAME
c     %NAME = value
c     %show
c
c     and ignoring comment lines (with a '#' in column one).
c
c     The only user-callable routines are:
c
c        check_directives ( call check_directives(unit_no: integer) )
c        getline    ( call getline(unit_no: integer, line: string) )
c        set_value  ( new_value = set_value(NAME: string, default: string))
c        defined    ( if (defined('NAME')) ...   )
c        insert     [experimental...]
c
c
c     Check_directives reads comment and directive lines, if present,
c     leaving the file ready for further processing.
c
c     Getline will provide the caller with the next input line that is
c     not a comment or a directive.
c
c     Set_value will return the value associated with the name NAME, or
c     a default value (which MUST BE PASSED AS A STRING) if NAME is not
c     in the symbol table.
c
c     %define and %delete are implemented as special cases of
c     the more general form
c
c     %NAME = value
c
c     by simply associating a value of '1' to the symbol to be %define(d)
c     and a value of '0' to the symbol to be %delete(d), and the routine
c     Defined simply calls Set_value.
c
c     
      subroutine check_directives(unit_no)
      
      implicit none

c     Reads any lines containing comments or directives, and readies
c     the file connected to unit_no for reading the next block of
c     information

      integer unit_no
      character*10 dummy_string

      call getline(unit_no,dummy_string)
      if (dummy_string .ne. '#') backspace(unit_no)

      return
      end
c
      subroutine getline(unit_no,string)

      implicit none

c     Provides the caller with a line that does not contain comments
c     or directives.
c
c     If it encounters the end of file, returns '#'
c
      integer unit_no
      character*(*) string

      character*132 line
      external directive

c 'line' is used as an internal file (it should be big enough to
c        accomodate requests of different sizes).

c Special lines:
c
c ---> Comment lines:
c      The character '#' must appear in column 1
c ---> A line containing '%' in column 1 is treated as a 
c      directive.
c
c      Blank lines are skipped
c

  10  continue

      read(unit_no,'(a132)',end=999) line

c  ...skip comment lines ( those include blank lines )

      if (line .eq. ' ' .or. line(1:1) .eq. '#') then
        write(6,'(a78)') line
c
        go to 10
c
      endif
c
c  ...process directives
c
      if (line(1:1) .eq. '%') then
c
	 write(6,'(a78)') line
         call directive(line(2:))
c
         go to 10
c
      endif
c
      string = line
      RETURN
 999  continue
      string = '#'
      return
c
      end   
  
c
      subroutine directive(str)

      implicit none

      character*(*) str
c
      character*30 name, val_str
      integer eqloc
      
      logical success
      logical leqi, put_pair
      external leqi, put_pair
c
      success = .true.

      eqloc = index(str,'=')
      
      if (eqloc .ne. 0) then
         name = str(1:eqloc-1)
         val_str = str(eqloc+1:)
         success = put_pair(name,val_str)
      else 
         name = str(8:)
         if (leqi(str(1:6),'define')) then
            success = put_pair(name,'1.')
         else if (leqi(str(1:6),'delete')) then
            success = put_pair(name,'0.')
         else if (leqi(str(1:4),'show')) then
            call print_pairs
         else
            write(6,'(/,a,1x,a,/)') ' Unrecognized directive:', str
         endif

      endif
      
      if (.not. success) write(6,'(/,a,2x,a,/)') 
     &     'Set full... Could not process ',name
      
      return
      end
c
      logical function defined(name)
c
      implicit none
c
      character*(*) name
      double precision set_value
      external set_value
c
      defined = (nint(set_value(name,'0')) .eq. 1)
c
      return
c
      end
c
      subroutine insert(module,name)
c
      implicit none
c
      character*(*) name, module
      logical success
      logical put_pair
c
      success = put_pair(name,'1.')
      write(6,*) ' <---- defining ', name, ' in ',module, success
c
      return
c
      end
c
      logical function put_pair(name,num_str)
	
      implicit none

      character*(*) name, num_str

      include 'set2.h'
      
      double precision value
      integer i

      logical leqi
      external leqi, get_real

      call get_real(num_str,value)

      do 10 i = 1, nels

         if (leqi(name,el_name(i))) then
            put_pair = .true.
            el_val(i) = value
cag            write(6,*) name, ' reset. New value: ',num_str
            return
         endif
 10      continue

         if (nels .ne. nmax) then
            nels = nels + 1
            el_name(nels) = name
            el_val(i) = value
            put_pair = .true.
cag            write(6,*) name, ' set to: ', num_str
         endif

         return

         end

         double precision function set_value(name,default_str)

         implicit none
         
         character*(*) name, default_str
         double precision default_value

         include 'set2.h'

         logical leqi
         external leqi
	 external get_real

         integer i
c
         call get_real(default_str,default_value)
         set_value = default_value

         do 10 i=1, nels
            if (leqi(el_name(i),name)) set_value = el_val(i)
 10      continue

         return

         end
         
c     
      subroutine get_real(string,value)
         
c     Alberto Garcia, April 2, 1995
c     
c     It turns out that compiler behavior is not uniform when it comes
c     to interpreting real input format statements. This routine takes
c     the string "string" and extracts the value it represents as a real
c     number. String can contain the letters e,d,E,D signifying an
c     exponent, and embedded blanks everywhere. 
c     The bn descriptor is used in case some compiler's default behavior is
c     to pad the mantissa or exponent fields with significant zeroes.
c     
c     I have tested the routine against all kinds of *reasonable* strings.
c     Please let me know if you find any bugs.
c     
         implicit none
         
         character*(*) string
         double precision value
c     
c     80 should be enough...
c     
         character*80 mantissa, exponent
c     
         integer exp_loc, exp_val
c     
c     Paranoid section.
c     Make sure that the routine still works if the case of the source
c     is changed! Sorry, ASCII only...
c  
      character*1 low_e, cap_e, low_d, cap_d
      
      cap_d = char(68)
      cap_e = char(69)
      low_d = char(100)
      low_e = char(101)

cdebug      write(6,*) cap_d, cap_e, low_d, low_e
c
c     We are counting on finding only one exponent... (reasonable?)
c

      exp_loc = index(string,cap_e) + index(string,cap_d) +
     +          index(string,low_e) + index(string,low_d)

      if (exp_loc .ne. 0) then
	mantissa = string(1:exp_loc-1)
	exponent = string(exp_loc+1:)
      else
	mantissa = string
      endif
c
c    BN means "Treat blanks as no-significant"
c
	read(mantissa,'(bn,f80.0)') value
	
	if (exp_loc .ne. 0) then
		read(exponent,'(bn,i80)') exp_val
		value = value * 10.d0 ** exp_val
	endif

	return

	End
c
      subroutine print_pairs
c
      implicit none
c
      include 'set2.h'
c
      integer i
c
      write(6,9000) nels, nmax
 9000 format(/,1x,'Symbol set contains ',i2,' pairs,',
     +         ' out of a maximum of ',i3,':')
c 
      write(6,9010) (el_name(i),el_val(i),i=1,nels)
 9010 format(1x,a30,2x,g25.15)
      write(6,'(/)')
c
      return
c
      end    








      subroutine prversion
c
c     Simple routine to print the version string. Could be extended to
c     provide more information, if needed.
c
c     Alberto Garcia, Feb. 23, 1998
c
      implicit none
      include 'version.h'

      write(6,'(a)') version

      end
      SUBROUTINE ATOMXC( FUNCTL, AUTHOR, IREL,
     .                   NR, MAXR, RMESH, NSPIN, DENS,
     .                   EX, EC, DX, DC, VXC )

C *******************************************************************
C Finds total exchange-correlation energy and potential for a
C spherical electron density distribution.
C This version implements the Local (spin) Density Approximation and
C the Generalized-Gradient-Aproximation with the 'explicit mesh 
C functional' method of White & Bird, PRB 50, 4954 (1994).
C Gradients are 'defined' by numerical derivatives, using 2*NN+1 mesh
C   points, where NN is a parameter defined below
C Coded by L.C.Balbas and J.M.Soler. December 1996. Version 0.5.
C ************************* INPUT ***********************************
C CHARACTER*(*) FUNCTL : Functional to be used:
C              'LDA' or 'LSD' => Local (spin) Density Approximation
C                       'GGA' => Generalized Gradient Corrections
C                                Uppercase is optional
C CHARACTER*(*) AUTHOR : Parametrization desired:
C     'CA' or 'PZ' => LSD Perdew & Zunger, PRB 23, 5075 (1981)
C           'PW91' => GGA Perdew & Wang, JCP, 100, 1290 (1994) 
C           'PW92' => LSD Perdew & Wang, PRB, 45, 13244 (1992). This is
C                     the local density limit of the next:
C            'PBE' => GGA Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)
C            'LYP' => GGA Becke-Lee-Yang-Parr (see subroutine blypxc)
C                     Uppercase is optional
C INTEGER IREL         : Relativistic exchange? (0=>no, 1=>yes)
C INTEGER NR           : Number of radial mesh points
C INTEGER MAXR         : Physical first dimension of RMESH, DENS and VXC
C REAL*8  RMESH(MAXR)  : Radial mesh points
C INTEGER NSPIN        : NSPIN=1 => unpolarized; NSPIN=2 => polarized
C REAL*8  DENS(MAXR,NSPIN) : Total (NSPIN=1) or spin (NSPIN=2) electron
C                            density at mesh points
C ************************* OUTPUT **********************************
C REAL*8  EX              : Total exchange energy
C REAL*8  EC              : Total correlation energy
C REAL*8  DX              : IntegralOf( rho * (eps_x - v_x) )
C REAL*8  DC              : IntegralOf( rho * (eps_c - v_c) )
C REAL*8  VXC(MAXR,NSPIN) : (Spin) exch-corr potential
C ************************ UNITS ************************************
C Distances in atomic units (Bohr).
C Densities in atomic units (electrons/Bohr**3)
C Energy unit depending of parameter EUNIT below
C ********* ROUTINES CALLED *****************************************
C GGAXC, LDAXC
C *******************************************************************

C Next line is nonstandard but may be suppressed
      IMPLICIT NONE

C Argument types and dimensions
      CHARACTER*(*)     FUNCTL, AUTHOR
      INTEGER           IREL, MAXR, NR, NSPIN
      DOUBLE PRECISION  DENS(MAXR,NSPIN), RMESH(MAXR), VXC(MAXR,NSPIN)
      DOUBLE PRECISION  DC, DX, EC, EX

C Internal parameters
C NN    : order of the numerical derivatives: the number of radial 
C          points used is 2*NN+1
C MR    : must be equal or larger than NR
C MSPIN : must be equal or larger than NSPIN (4 for non-collinear spin)
      INTEGER NN, MR, MSPIN
      PARAMETER ( NN    =    5 )
      PARAMETER ( MR    = 2048 )
      PARAMETER ( MSPIN =    4 )

C Fix energy unit:  EUNIT=1.0 => Hartrees,
C                   EUNIT=0.5 => Rydbergs,
C                   EUNIT=0.03674903 => eV
      DOUBLE PRECISION EUNIT
      PARAMETER ( EUNIT = 0.5D0 )

C DVMIN is added to differential of volume to avoid division by zero
      DOUBLE PRECISION DVMIN
      PARAMETER ( DVMIN = 1.D-12 )

C Local variables and arrays
      LOGICAL
     .  GGA
      INTEGER
     .  IN, IN1, IN2, IR, IS, JN
      DOUBLE PRECISION
     .  AUX(MR), D(MSPIN), DECDD(MSPIN), DECDGD(3,MSPIN),
     .  DEXDD(MSPIN), DEXDGD(3,MSPIN),
     .  DGDM(-NN:NN), DGIDFJ(-NN:NN), DRDM, DVOL, 
     .  DVCDN(MSPIN,MSPIN), DVXDN(MSPIN,MSPIN),
     .  EPSC, EPSX, F1, F2, GD(3,MSPIN), PI
      EXTERNAL
     .  GGAXC, LDAXC

C Set GGA switch
      IF ( FUNCTL.EQ.'LDA' .OR. FUNCTL.EQ.'lda' .OR.
     .     FUNCTL.EQ.'LSD' .OR. FUNCTL.EQ.'lsd' ) THEN
        GGA = .FALSE.
      ELSEIF ( FUNCTL.EQ.'GGA' .OR. FUNCTL.EQ.'gga') THEN
        GGA = .TRUE.
        IF (MR.LT.NR) STOP 'ATOMXC: Parameter MR too small'
      ELSE
        WRITE(6,*) 'ATOMXC: Unknown functional ', FUNCTL
        STOP
      ENDIF

C Initialize output
      EX = 0
      EC = 0
      DX = 0
      DC = 0
      DO 20 IS = 1,NSPIN
        DO 10 IR = 1,NR
          VXC(IR,IS) = 0
   10   CONTINUE
   20 CONTINUE

C Get number pi
      PI = 4 * ATAN(1.D0)

C Loop on mesh points
      DO 140 IR = 1,NR

C       Find interval of neighbour points to calculate derivatives
        IN1 = MAX(  1, IR-NN ) - IR
        IN2 = MIN( NR, IR+NN ) - IR

C       Find weights of numerical derivation from Lagrange
C       interpolation formula
        DO 50 IN = IN1,IN2
          IF (IN .EQ. 0) THEN
            DGDM(IN) = 0
            DO 30 JN = IN1,IN2
              IF (JN.NE.0) DGDM(IN) = DGDM(IN) + 1.D0 / (0 - JN)
   30       CONTINUE
          ELSE
            F1 = 1
            F2 = 1
            DO 40 JN = IN1,IN2
              IF (JN.NE.IN .AND. JN.NE.0) F1 = F1 * (0  - JN)
              IF (JN.NE.IN)               F2 = F2 * (IN - JN)
   40       CONTINUE
            DGDM(IN) = F1 / F2
          ENDIF
   50   CONTINUE

C       Find dr/dmesh
        DRDM = 0
        DO 60 IN = IN1,IN2
          DRDM = DRDM + RMESH(IR+IN) * DGDM(IN)
   60   CONTINUE

C       Find differential of volume. Use trapezoidal integration rule
        DVOL = 4 * PI * RMESH(IR)**2 * DRDM
C       DVMIN is a small number added to avoid a division by zero
        DVOL = DVOL + DVMIN
        IF (IR.EQ.1 .OR. IR.EQ.NR) DVOL = DVOL / 2
        IF (GGA) AUX(IR) = DVOL

C       Find the weights for the derivative d(gradF(i))/d(F(j)), of
C       the gradient at point i with respect to the value at point j
        IF (GGA) THEN
          DO 80 IN = IN1,IN2
            DGIDFJ(IN) = DGDM(IN) / DRDM
   80     CONTINUE
        ENDIF

C       Find density and gradient of density at this point
        DO 90 IS = 1,NSPIN
          D(IS) = DENS(IR,IS)
   90   CONTINUE
        IF (GGA) THEN
          DO 110 IS = 1,NSPIN
            GD(1,IS) = 0
            GD(2,IS) = 0
            GD(3,IS) = 0
            DO 100 IN = IN1,IN2
              GD(3,IS) = GD(3,IS) + DGIDFJ(IN) * DENS(IR+IN,IS)
  100       CONTINUE
  110     CONTINUE
        ENDIF

C       Find exchange and correlation energy densities and their 
C       derivatives with respect to density and density gradient
        IF (GGA) THEN
          CALL GGAXC( AUTHOR, IREL, NSPIN, D, GD,
     .                EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
        ELSE
          CALL LDAXC( AUTHOR, IREL, NSPIN, D, EPSX, EPSC, DEXDD, DECDD,
     .                DVXDN, DVCDN )
        ENDIF

C       Add contributions to exchange-correlation energy and its
C       derivatives with respect to density at all points
        DO 130 IS = 1,NSPIN
          EX = EX + DVOL * D(IS) * EPSX
          EC = EC + DVOL * D(IS) * EPSC
          DX = DX + DVOL * D(IS) * (EPSX - DEXDD(IS))
          DC = DC + DVOL * D(IS) * (EPSC - DECDD(IS))
          IF (GGA) THEN
            VXC(IR,IS) = VXC(IR,IS) + DVOL * ( DEXDD(IS) + DECDD(IS) )
            DO 120 IN = IN1,IN2
              DX= DX - DVOL * DENS(IR+IN,IS) * DEXDGD(3,IS) * DGIDFJ(IN)
              DC= DC - DVOL * DENS(IR+IN,IS) * DECDGD(3,IS) * DGIDFJ(IN)
              VXC(IR+IN,IS) = VXC(IR+IN,IS) + DVOL *
     .               (DEXDGD(3,IS) + DECDGD(3,IS)) * DGIDFJ(IN)
  120       CONTINUE
          ELSE
            VXC(IR,IS) = DEXDD(IS) + DECDD(IS)
          ENDIF
  130   CONTINUE

  140 CONTINUE

C Divide by volume element to obtain the potential (per electron)
      IF (GGA) THEN
        DO 160 IS = 1,NSPIN
          DO 150 IR = 1,NR
            DVOL = AUX(IR)
            VXC(IR,IS) = VXC(IR,IS) / DVOL
  150     CONTINUE
  160   CONTINUE
      ENDIF

C Divide by energy unit
      EX = EX / EUNIT
      EC = EC / EUNIT
      DX = DX / EUNIT
      DC = DC / EUNIT
      DO 180 IS = 1,NSPIN
        DO 170 IR = 1,NR
          VXC(IR,IS) = VXC(IR,IS) / EUNIT
  170   CONTINUE
  180 CONTINUE

      END


      SUBROUTINE EXCHNG( IREL, NSP, DS, EX, VX )

C *****************************************************************
C  Finds local exchange energy density and potential
C  Adapted by J.M.Soler from routine velect of Froyen's 
C    pseudopotential generation program. Madrid, Jan'97. Version 0.5.
C **** Input ******************************************************
C INTEGER IREL    : relativistic-exchange switch (0=no, 1=yes)
C INTEGER NSP     : spin-polarizations (1=>unpolarized, 2=>polarized)
C REAL*8  DS(NSP) : total (nsp=1) or spin (nsp=2) electron density
C **** Output *****************************************************
C REAL*8  EX      : exchange energy density
C REAL*8  VX(NSP) : (spin-dependent) exchange potential
C **** Units ******************************************************
C Densities in electrons/Bohr**3
C Energies in Hartrees
C *****************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.D0,ONE=1.D0,PFIVE=.5D0,OPF=1.5D0,C014=0.014D0)
      DIMENSION DS(NSP), VX(NSP)

       PI=4*ATAN(ONE)
       TRD = ONE/3
       FTRD = 4*TRD
       TFTM = 2**FTRD-2
       A0 = (4/(9*PI))**TRD

C      X-alpha parameter:       
       ALP = 2 * TRD

       IF (NSP .EQ. 2) THEN
         D1 = MAX(DS(1),0.D0)
         D2 = MAX(DS(2),0.D0)
         D = D1 + D2
         IF (D .LE. ZERO) THEN
           EX = ZERO
           VX(1) = ZERO
           VX(2) = ZERO
           RETURN
         ENDIF
         Z = (D1 - D2) / D
         FZ = ((1+Z)**FTRD+(1-Z)**FTRD-2)/TFTM
         FZP = FTRD*((1+Z)**TRD-(1-Z)**TRD)/TFTM 
       ELSE
         D = DS(1)
         IF (D .LE. ZERO) THEN
           EX = ZERO
           VX(1) = ZERO
           RETURN
         ENDIF
         Z = ZERO
         FZ = ZERO
         FZP = ZERO
       ENDIF
       RS = (3 / (4*PI*D) )**TRD
       VXP = -(3*ALP/(2*PI*A0*RS))
       EXP_VAR = 3*VXP/4
       IF (IREL .EQ. 1) THEN
         BETA = C014/RS
         SB = SQRT(1+BETA*BETA)
         ALB = LOG(BETA+SB)
         VXP = VXP * (-PFIVE + OPF * ALB / (BETA*SB))
         EXP_VAR = EXP_VAR * (ONE-OPF*((BETA*SB-ALB)/BETA**2)**2) 
       ENDIF
       VXF = 2**TRD*VXP
       EXF = 2**TRD*EXP_VAR
       IF (NSP .EQ. 2) THEN
         VX(1) = VXP + FZ*(VXF-VXP) + (1-Z)*FZP*(EXF-EXP_VAR)
         VX(2) = VXP + FZ*(VXF-VXP) - (1+Z)*FZP*(EXF-EXP_VAR)
         EX    = EXP_VAR + FZ*(EXF-EXP_VAR)
       ELSE
         VX(1) = VXP
         EX    = EXP_VAR
       ENDIF
      END




      SUBROUTINE GGAXC( AUTHOR, IREL, NSPIN, D, GD,
     .                  EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )

C Finds the exchange and correlation energies at a point, and their
C derivatives with respect to density and density gradient, in the
C Generalized Gradient Correction approximation.
C Lengths in Bohr, energies in Hartrees
C Written by L.C.Balbas and J.M.Soler, Dec'96. Version 0.5.
C Modified by V.M.Garcia-Suarez to include non-collinear spin. June 2002

      IMPLICIT          NONE
      CHARACTER*(*)     AUTHOR
      INTEGER           IREL, NSPIN, NS, IS, IX
      DOUBLE PRECISION  THETA, PHI, D(NSPIN), DECDD(NSPIN),
     .                  DECDGD(3,NSPIN), DEXDD(NSPIN), DEXDGD(3,NSPIN),
     .                  EPSC, EPSX, GD(3,NSPIN),
     .                  DD(2), DTOT, DPOL, GDTOT(3), GDPOL(3),
     .                  GDD(3,2), TINY, DECDN(2), DEXDN(2),
     .                  VPOL, DECDGN(3,2), DEXDGN(3,2),
     .                  VGPOLX, VGPOLC, C2, S2, ST, CP, SP

      PARAMETER ( TINY = 1.D-12 )

      IF (NSPIN .EQ. 4) THEN
C       Find eigenvalues of density matrix (up and down densities
C       along the spin direction)
C       Note: D(1)=D11, D(2)=D22, D(3)=Real(D12), D(4)=Im(D12)
        NS = 2
        DTOT = D(1) + D(2)
        DPOL = SQRT( (D(1)-D(2))**2 + 4.D0*(D(3)**2+D(4)**2) )
        DD(1) = 0.5D0 * ( DTOT + DPOL ) 
        DD(2) = 0.5D0 * ( DTOT - DPOL )
        THETA = ACOS((D(1)-D(2))/(DPOL+TINY))
        C2 = COS(THETA/2)
        S2 = SIN(THETA/2)
        ST = SIN(THETA)
        PHI = ATAN(-D(4)/(D(3)+TINY))
        CP = COS(PHI)
        SP = SIN(PHI)
C       Find diagonal elements of the gradient
        DO 10 IX = 1,3
          GDD(IX,1) = GD(IX,1)*C2**2 + GD(IX,2)*S2**2 +
     .                2.d0*C2*S2*(GD(IX,3)*CP - GD(IX,4)*SP)
          GDD(IX,2) = GD(IX,1)*S2**2 + GD(IX,2)*C2**2 -
     .                2.d0*C2*S2*(GD(IX,3)*CP - GD(IX,4)*SP)
   10   CONTINUE 
      ELSE
        NS = NSPIN
        DO 20 IS = 1,NSPIN
cag       Avoid negative densities
          DD(IS) = max(D(IS),0.0d0)
          DO 30 IX = 1,3
            GDD(IX,IS) = GD(IX,IS)
   30     CONTINUE
   20   CONTINUE
      ENDIF

      IF (AUTHOR.EQ.'PBE' .OR. AUTHOR.EQ.'pbe') THEN
        CALL PBEXC( IREL, NS, DD, GDD,
     .              EPSX, EPSC, DEXDN, DECDN, DEXDGN, DECDGN )
cag
      ELSE IF (AUTHOR.EQ.'RPBE' .OR. AUTHOR.EQ.'rpbe') THEN
        CALL RPBEXC( IREL, NS, DD, GDD,
     .              EPSX, EPSC, DEXDN, DECDN, DEXDGN, DECDGN )
      ELSE IF (AUTHOR.EQ.'REVPBE' .OR. AUTHOR.EQ.'revpbe') THEN
        CALL REVPBEXC( IREL, NS, DD, GDD,
     .              EPSX, EPSC, DEXDN, DECDN, DEXDGN, DECDGN )
      ELSE IF (AUTHOR.EQ.'LYP'.OR.AUTHOR.EQ.'lyp') THEN
        CALL BLYPXC(NSPIN,D,GD,EPSX,EPSC,dEXdn,dECdn,dEXdgn,dECdgn)
cag
      ELSEIF (AUTHOR.EQ.'PW91' .OR. AUTHOR.EQ.'pw91') THEN
        CALL PW91XC( IREL, NS, DD, GDD,
     .               EPSX, EPSC, DEXDN, DECDN, DEXDGN, DECDGN )
      ELSE
        WRITE(6,*) 'GGAXC: Unknown author ', AUTHOR
        STOP
      ENDIF

      IF (NSPIN .EQ. 4) THEN
C       Find dE/dD(ispin) = dE/dDup * dDup/dD(ispin) +
C                           dE/dDdown * dDown/dD(ispin)
        VPOL  = (DEXDN(1)-DEXDN(2)) * (D(1)-D(2)) / (DPOL+TINY)
        DEXDD(1) = 0.5D0 * ( DEXDN(1) + DEXDN(2) + VPOL )
        DEXDD(2) = 0.5D0 * ( DEXDN(1) + DEXDN(2) - VPOL )
        DEXDD(3) = (DEXDN(1)-DEXDN(2)) * D(3) / (DPOL+TINY)
        DEXDD(4) = (DEXDN(1)-DEXDN(2)) * D(4) / (DPOL+TINY)
        VPOL  = (DECDN(1)-DECDN(2)) * (D(1)-D(2)) / (DPOL+TINY)
        DECDD(1) = 0.5D0 * ( DECDN(1) + DECDN(2) + VPOL )
        DECDD(2) = 0.5D0 * ( DECDN(1) + DECDN(2) - VPOL )
        DECDD(3) = (DECDN(1)-DECDN(2)) * D(3) / (DPOL+TINY)
        DECDD(4) = (DECDN(1)-DECDN(2)) * D(4) / (DPOL+TINY)
C       Gradient terms
        DO 40 IX = 1,3
          DEXDGD(IX,1) = DEXDGN(IX,1)*C2**2 + DEXDGN(IX,2)*S2**2
          DEXDGD(IX,2) = DEXDGN(IX,1)*S2**2 + DEXDGN(IX,2)*C2**2
          DEXDGD(IX,3) = 0.5D0*(DEXDGN(IX,1) - DEXDGN(IX,2))*ST*CP
          DEXDGD(IX,4) = 0.5D0*(DEXDGN(IX,2) - DEXDGN(IX,1))*ST*SP
          DECDGD(IX,1) = DECDGN(IX,1)*C2**2 + DECDGN(IX,2)*S2**2
          DECDGD(IX,2) = DECDGN(IX,1)*S2**2 + DECDGN(IX,2)*C2**2
          DECDGD(IX,3) = 0.5D0*(DECDGN(IX,1) - DECDGN(IX,2))*ST*CP
          DECDGD(IX,4) = 0.5D0*(DECDGN(IX,2) - DECDGN(IX,1))*ST*SP
   40   CONTINUE
      ELSE
        DO 60 IS = 1,NSPIN
          DEXDD(IS) = DEXDN(IS)
          DECDD(IS) = DECDN(IS)
          DO 50 IX = 1,3
            DEXDGD(IX,IS) = DEXDGN(IX,IS)
            DECDGD(IX,IS) = DECDGN(IX,IS)
   50     CONTINUE
   60   CONTINUE
      ENDIF

      END


      SUBROUTINE LDAXC( AUTHOR, IREL, NSPIN, D, EPSX, EPSC, VX, VC,
     .                  DVXDN, DVCDN )

C ******************************************************************
C Finds the exchange and correlation energies and potentials, in the
C Local (spin) Density Approximation.
C Written by L.C.Balbas and J.M.Soler, Dec'96.
C Non-collinear spin added by J.M.Soler, May'98
C *********** INPUT ************************************************
C CHARACTER*(*) AUTHOR : Parametrization desired:
C     'CA' or 'PZ' => LSD Perdew & Zunger, PRB 23, 5075 (1981)
C           'PW92' => LSD Perdew & Wang, PRB, 45, 13244 (1992)
C                     Uppercase is optional
C INTEGER IREL     : Relativistic exchange? (0=>no, 1=>yes)
C INTEGER NSPIN    : NSPIN=1 => unpolarized; NSPIN=2 => polarized;
C                    NSPIN=4 => non-collinear polarization
C REAL*8  D(NSPIN) : Local (spin) density. For non-collinear
C                    polarization, the density matrix is given by:
C                    D(1)=D11, D(2)=D22, D(3)=Real(D12), D(4)=Im(D12)
C *********** OUTPUT ***********************************************
C REAL*8 EPSX, EPSC : Exchange and correlation energy densities
C REAL*8 VX(NSPIN), VC(NSPIN) : Exchange and correlation potentials,
C                               defined as dExc/dD(ispin)
C REAL*8 DVXDN(NSPIN,NSPIN)  :  Derivative of exchange potential with
C                               respect the charge density, defined 
C                               as DVx(spin1)/Dn(spin2)
C REAL*8 DVCDN(NSPIN,NSPIN)  :  Derivative of correlation potential
C                               respect the charge density, defined 
C                               as DVc(spin1)/Dn(spin2)
C *********** UNITS ************************************************
C Lengths in Bohr, energies in Hartrees
C ******************************************************************

      IMPLICIT          NONE
      CHARACTER*(*)     AUTHOR
      INTEGER           IREL, NSPIN
      DOUBLE PRECISION  D(NSPIN), EPSC, EPSX, VX(NSPIN), VC(NSPIN),
     .                  DVXDN(NSPIN,NSPIN), DVCDN(NSPIN,NSPIN)

      INTEGER           IS, NS, ISPIN1, ISPIN2
      DOUBLE PRECISION  DD(2), DPOL, DTOT, TINY, VCD(2), VPOL, VXD(2)

      PARAMETER ( TINY = 1.D-12 )

      IF (NSPIN .EQ. 4) THEN
C       Find eigenvalues of density matrix (up and down densities
C       along the spin direction)
C       Note: D(1)=D11, D(2)=D22, D(3)=Real(D12), D(4)=Im(D12)
        NS = 2
        DTOT = D(1) + D(2)
        DPOL = SQRT( (D(1)-D(2))**2 + 4.D0*(D(3)**2+D(4)**2) )
        DD(1) = 0.5D0 * ( DTOT + DPOL )
        DD(2) = 0.5D0 * ( DTOT - DPOL )
      ELSE
        NS = NSPIN
        DO 10 IS = 1,NSPIN
cag       Avoid negative densities
          DD(IS) = max(D(IS),0.0d0)
   10   CONTINUE
      ENDIF


      DO ISPIN2 = 1, NSPIN
        DO ISPIN1 = 1, NSPIN
          DVXDN(ISPIN1,ISPIN2) = 0.D0
          DVCDN(ISPIN1,ISPIN2) = 0.D0
        ENDDO
      ENDDO

      IF ( AUTHOR.EQ.'CA' .OR. AUTHOR.EQ.'ca' .OR.
     .     AUTHOR.EQ.'PZ' .OR. AUTHOR.EQ.'pz') THEN
        CALL PZXC( IREL, NS, DD, EPSX, EPSC, VXD, VCD, DVXDN, DVCDN )
      ELSEIF ( AUTHOR.EQ.'PW92' .OR. AUTHOR.EQ.'pw92' ) THEN
        CALL PW92XC( IREL, NS, DD, EPSX, EPSC, VXD, VCD )
      ELSE
        WRITE(6,*) 'LDAXC: Unknown author ', AUTHOR
        STOP
      ENDIF

      IF (NSPIN .EQ. 4) THEN
C       Find dE/dD(ispin) = dE/dDup * dDup/dD(ispin) +
C                           dE/dDdown * dDown/dD(ispin)
        VPOL  = (VXD(1)-VXD(2)) * (D(1)-D(2)) / (DPOL+TINY)
        VX(1) = 0.5D0 * ( VXD(1) + VXD(2) + VPOL )
        VX(2) = 0.5D0 * ( VXD(1) + VXD(2) - VPOL )
        VX(3) = (VXD(1)-VXD(2)) * D(3) / (DPOL+TINY)
        VX(4) = (VXD(1)-VXD(2)) * D(4) / (DPOL+TINY)
        VPOL  = (VCD(1)-VCD(2)) * (D(1)-D(2)) / (DPOL+TINY)
        VC(1) = 0.5D0 * ( VCD(1) + VCD(2) + VPOL )
        VC(2) = 0.5D0 * ( VCD(1) + VCD(2) - VPOL )
        VC(3) = (VCD(1)-VCD(2)) * D(3) / (DPOL+TINY)
        VC(4) = (VCD(1)-VCD(2)) * D(4) / (DPOL+TINY)
      ELSE
        DO 20 IS = 1,NSPIN
          VX(IS) = VXD(IS)
          VC(IS) = VCD(IS)
   20   CONTINUE
      ENDIF
      END



      SUBROUTINE PBEXC( IREL, NSPIN, DENS, GDENS,
     .                  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Perdew-Burke-Ernzerhof Generalized-Gradient-Approximation.
C Ref: J.P.Perdew, K.Burke & M.Ernzerhof, PRL 77, 3865 (1996)
C Written by L.C.Balbas and J.M.Soler. December 1996. Version 0.5.
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER NSPIN          : Number of spin polarizations (1 or 2)
C REAL*8  DENS(NSPIN)    : Total electron density (if NSPIN=1) or
C                           spin electron density (if NSPIN=2)
C REAL*8  GDENS(3,NSPIN) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(NSPIN)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( DENS(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(NSPIN)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( DENS(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,NSPIN): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,NSPIN): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C EXCHNG, PW92C
C ********************************************************************

      IMPLICIT          NONE
      INTEGER           IREL, NSPIN
      DOUBLE PRECISION  DENS(NSPIN), DECDD(NSPIN), DECDGD(3,NSPIN),
     .                  DEXDD(NSPIN), DEXDGD(3,NSPIN), GDENS(3,NSPIN)

C Internal variables
      INTEGER
     .  IS, IX

      DOUBLE PRECISION
     .  A, BETA, D(2), DADD, DECUDD, DENMIN, 
     .  DF1DD, DF2DD, DF3DD, DF4DD, DF1DGD, DF3DGD, DF4DGD,
     .  DFCDD(2), DFCDGD(3,2), DFDD, DFDGD, DFXDD(2), DFXDGD(3,2),
     .  DHDD, DHDGD, DKFDD, DKSDD, DPDD, DPDZ, DRSDD, 
     .  DS(2), DSDD, DSDGD, DT, DTDD, DTDGD, DZDD(2), 
     .  EC, ECUNIF, EX, EXUNIF,
     .  F, F1, F2, F3, F4, FC, FX, FOUTHD,
     .  GAMMA, GD(3,2), GDM(2), GDMIN, GDMS, GDMT, GDS, GDT(3),
     .  H, HALF, KAPPA, KF, KFS, KS, MU, PHI, PI, RS, S,
     .  T, THD, THRHLF, TWO, TWOTHD, VCUNIF(2), VXUNIF(2), ZETA

C Lower bounds of density and its gradient to avoid divisions by zero
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( GDMIN  = 1.D-12 )

C Fix some numerical parameters
      PARAMETER ( FOUTHD=4.D0/3.D0, HALF=0.5D0,
     .            THD=1.D0/3.D0, THRHLF=1.5D0,
     .            TWO=2.D0, TWOTHD=2.D0/3.D0 )

C Fix some more numerical constants
      PI = 4 * ATAN(1.D0)
      BETA = 0.066725D0
      GAMMA = (1 - LOG(TWO)) / PI**2
      MU = BETA * PI**2 / 3
      KAPPA = 0.804D0

C Translate density and its gradient to new variables
      IF (NSPIN .EQ. 1) THEN
        D(1) = HALF * DENS(1)
        D(2) = D(1)
        DT = MAX( DENMIN, DENS(1) )
        DO 10 IX = 1,3
          GD(IX,1) = HALF * GDENS(IX,1)
          GD(IX,2) = GD(IX,1)
          GDT(IX) = GDENS(IX,1)
   10   CONTINUE
      ELSE
        D(1) = DENS(1)
        D(2) = DENS(2)
        DT = MAX( DENMIN, DENS(1)+DENS(2) )
        DO 20 IX = 1,3
          GD(IX,1) = GDENS(IX,1)
          GD(IX,2) = GDENS(IX,2)
          GDT(IX) = GDENS(IX,1) + GDENS(IX,2)
   20   CONTINUE
      ENDIF
      GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
      GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
      GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
      GDMT = MAX( GDMIN, GDMT )

C Find local correlation energy and potential
      CALL PW92C( 2, D, ECUNIF, VCUNIF )

C Find total correlation energy
      RS = ( 3 / (4*PI*DT) )**THD
      KF = (3 * PI**2 * DT)**THD
      KS = SQRT( 4 * KF / PI )
      ZETA = ( D(1) - D(2) ) / DT
      ZETA = MAX( -1.D0+DENMIN, ZETA )
      ZETA = MIN(  1.D0-DENMIN, ZETA )
      PHI = HALF * ( (1+ZETA)**TWOTHD + (1-ZETA)**TWOTHD )
      T = GDMT / (2 * PHI * KS * DT)
      F1 = ECUNIF / GAMMA / PHI**3
      F2 = EXP(-F1)
      A = BETA / GAMMA / (F2-1)
      F3 = T**2 + A * T**4
      F4 = BETA/GAMMA * F3 / (1 + A*F3)
      H = GAMMA * PHI**3 * LOG( 1 + F4 )
      FC = ECUNIF + H

C Find correlation energy derivatives
      DRSDD = - (THD * RS / DT)
      DKFDD =   THD * KF / DT
      DKSDD = HALF * KS * DKFDD / KF
      DZDD(1) =   1 / DT - ZETA / DT
      DZDD(2) = - (1 / DT) - ZETA / DT
      DPDZ = HALF * TWOTHD * ( 1/(1+ZETA)**THD - 1/(1-ZETA)**THD )
      DO 40 IS = 1,2
        DECUDD = ( VCUNIF(IS) - ECUNIF ) / DT
        DPDD = DPDZ * DZDD(IS)
        DTDD = (- T) * ( DPDD/PHI + DKSDD/KS + 1/DT )
        DF1DD = F1 * ( DECUDD/ECUNIF - 3*DPDD/PHI )
        DF2DD = (- F2) * DF1DD
        DADD = (- A) * DF2DD / (F2-1)
        DF3DD = (2*T + 4*A*T**3) * DTDD + DADD * T**4
        DF4DD = F4 * ( DF3DD/F3 - (DADD*F3+A*DF3DD)/(1+A*F3) )
        DHDD = 3 * H * DPDD / PHI
        DHDD = DHDD + GAMMA * PHI**3 * DF4DD / (1+F4)
        DFCDD(IS) = VCUNIF(IS) + H + DT * DHDD

        DO 30 IX = 1,3
          DTDGD = (T / GDMT) * GDT(IX) / GDMT
          DF3DGD = DTDGD * ( 2 * T + 4 * A * T**3 )
          DF4DGD = F4 * DF3DGD * ( 1/F3 - A/(1+A*F3) ) 
          DHDGD = GAMMA * PHI**3 * DF4DGD / (1+F4)
          DFCDGD(IX,IS) = DT * DHDGD
   30   CONTINUE
   40 CONTINUE

C Find exchange energy and potential
      FX = 0
      DO 60 IS = 1,2
        DS(IS)   = MAX( DENMIN, 2 * D(IS) )
        GDMS = MAX( GDMIN, 2 * GDM(IS) )
        KFS = (3 * PI**2 * DS(IS))**THD
        S = GDMS / (2 * KFS * DS(IS))
        F1 = 1 + MU * S**2 / KAPPA
        F = 1 + KAPPA - KAPPA / F1
c
c       Note nspin=1 in call to exchng...
c
        CALL EXCHNG( IREL, 1, DS(IS), EXUNIF, VXUNIF(IS) )
        FX = FX + DS(IS) * EXUNIF * F

        DKFDD = THD * KFS / DS(IS)
        DSDD = S * ( -(DKFDD/KFS) - 1/DS(IS) )
        DF1DD = 2 * (F1-1) * DSDD / S
        DFDD = KAPPA * DF1DD / F1**2
        DFXDD(IS) = VXUNIF(IS) * F + DS(IS) * EXUNIF * DFDD

        DO 50 IX = 1,3
          GDS = 2 * GD(IX,IS)
          DSDGD = (S / GDMS) * GDS / GDMS
          DF1DGD = 2 * MU * S * DSDGD / KAPPA
          DFDGD = KAPPA * DF1DGD / F1**2
          DFXDGD(IX,IS) = DS(IS) * EXUNIF * DFDGD
   50   CONTINUE
   60 CONTINUE
      FX = HALF * FX / DT

C Set output arguments
      EX = FX
      EC = FC
      DO 90 IS = 1,NSPIN
        DEXDD(IS) = DFXDD(IS)
        DECDD(IS) = DFCDD(IS)
        DO 80 IX = 1,3
          DEXDGD(IX,IS) = DFXDGD(IX,IS)
          DECDGD(IX,IS) = DFCDGD(IX,IS)
   80   CONTINUE
   90 CONTINUE

      END




      SUBROUTINE PW91XC( IREL, NSPIN, DENS, GDENS,
     .                  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Perdew-Wang91 Generalized-Gradient-Approximation.
C Ref: JCP 100, 1290 (1994)
C Written by J.L. Martins  August 2000
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER NSPIN          : Number of spin polarizations (1 or 2)
C REAL*8  DENS(NSPIN)    : Total electron density (if NSPIN=1) or
C                           spin electron density (if NSPIN=2)
C REAL*8  GDENS(3,NSPIN) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(NSPIN)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( DENS(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(NSPIN)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( DENS(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,NSPIN): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,NSPIN): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C EXCHNG, PW92C
C ********************************************************************

      IMPLICIT          NONE
      INTEGER           IREL, NSPIN
      DOUBLE PRECISION  DENS(NSPIN), DECDD(NSPIN), DECDGD(3,NSPIN),
     .                  DEXDD(NSPIN), DEXDGD(3,NSPIN), GDENS(3,NSPIN)

C Internal variables
      INTEGER
     .  IS, IX
      DOUBLE PRECISION
     .  A, BETA, D(2), DADD, DECUDD, DENMIN, 
     .  DF1DD, DF2DD, DF3DD, DF4DD, DF1DGD, DF3DGD, DF4DGD,
     .  DFCDD(2), DFCDGD(3,2), DFDD, DFDGD, DFXDD(2), DFXDGD(3,2),
     .  DHDD, DHDGD, DKFDD, DKSDD, DPDD, DPDZ, DRSDD, 
     .  DS(2), DSDD, DSDGD, DT, DTDD, DTDGD, DZDD(2), 
     .  EC, ECUNIF, EX, EXUNIF,
     .  F, F1, F2, F3, F4, FC, FX, FOUTHD,
     .  GAMMA, GD(3,2), GDM(2), GDMIN, GDMS, GDMT, GDS, GDT(3),
     .  H, HALF, KF, KFS, KS, PHI, PI, RS, S,
     .  T, THD, THRHLF, TWO, TWOTHD, VCUNIF(2), VXUNIF(2), ZETA
     
      DOUBLE PRECISION F5, F6, F7, F8, ASINHS
      DOUBLE PRECISION DF5DD,DF6DD,DF7DD,DF8DD
      DOUBLE PRECISION DF1DS, DF2DS, DF3DS, DFDS, DF7DGD

C Lower bounds of density and its gradient to avoid divisions by zero
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( GDMIN  = 1.D-12 )

C Fix some numerical parameters
      PARAMETER ( FOUTHD=4.D0/3.D0, HALF=0.5D0,
     .            THD=1.D0/3.D0, THRHLF=1.5D0,
     .            TWO=2.D0, TWOTHD=2.D0/3.D0 )

C Fix some more numerical constants
      PI = 4 * ATAN(1.D0)
      BETA = 15.75592D0 * 0.004235D0
      GAMMA = BETA**2 / (2.0D0 * 0.09D0)

C Translate density and its gradient to new variables
      IF (NSPIN .EQ. 1) THEN
        D(1) = HALF * DENS(1)
        D(2) = D(1)
        DT = MAX( DENMIN, DENS(1) )
        DO 10 IX = 1,3
          GD(IX,1) = HALF * GDENS(IX,1)
          GD(IX,2) = GD(IX,1)
          GDT(IX) = GDENS(IX,1)
   10   CONTINUE
      ELSE
        D(1) = DENS(1)
        D(2) = DENS(2)
        DT = MAX( DENMIN, DENS(1)+DENS(2) )
        DO 20 IX = 1,3
          GD(IX,1) = GDENS(IX,1)
          GD(IX,2) = GDENS(IX,2)
          GDT(IX) = GDENS(IX,1) + GDENS(IX,2)
   20   CONTINUE
      ENDIF
      GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
      GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
      GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
      GDMT = MAX( GDMIN, GDMT )

C Find local correlation energy and potential
      CALL PW92C( 2, D, ECUNIF, VCUNIF )

C Find total correlation energy
      RS = ( 3 / (4*PI*DT) )**THD
      KF = (3 * PI**2 * DT)**THD
      KS = SQRT( 4 * KF / PI )
      S = GDMT / (2 * KF * DT)
      T = GDMT / (2 * KS * DT)
      ZETA = ( D(1) - D(2) ) / DT
      ZETA = MAX( -1.D0+DENMIN, ZETA )
      ZETA = MIN(  1.D0-DENMIN, ZETA )
      PHI = HALF * ( (1+ZETA)**TWOTHD + (1-ZETA)**TWOTHD )
      F1 = ECUNIF / GAMMA / PHI**3
      F2 = EXP(-F1)
      A = BETA / GAMMA / (F2-1)
      F3 = T**2 + A * T**4
      F4 = BETA/GAMMA * F3 / (1 + A*F3)
      F5 = 0.002568D0 + 0.023266D0*RS + 7.389D-6*RS**2
      F6 = 1.0D0 + 8.723D0*RS + 0.472D0*RS**2 + 0.07389D0*RS**3
      F7 = EXP(-100.0D0 * S**2 * PHI**4)
      F8 =  15.75592D0*(0.001667212D0 + F5/F6 -0.004235D0 + 
     .          3.0D0*0.001667212D0/7.0D0)
      H = GAMMA * PHI**3 * LOG( 1 + F4 ) + F8 * T**2 * F7
      FC = ECUNIF + H

C Find correlation energy derivatives
      DRSDD = - THD * RS / DT
      DKFDD =   THD * KF / DT
      DKSDD = HALF * KS * DKFDD / KF
      DZDD(1) =   1 / DT - ZETA / DT
      DZDD(2) = - 1 / DT - ZETA / DT
      DPDZ = HALF * TWOTHD * ( 1/(1+ZETA)**THD - 1/(1-ZETA)**THD )
      DO 40 IS = 1,2
        DECUDD = ( VCUNIF(IS) - ECUNIF ) / DT
        DPDD = DPDZ * DZDD(IS)
        DTDD = - T * ( DPDD/PHI + DKSDD/KS + 1/DT )
        DSDD = - S * ( DPDD/PHI + DKFDD/KF + 1/DT )
        DF1DD = F1 * ( DECUDD/ECUNIF - 3*DPDD/PHI )
        DF2DD = - F2 * DF1DD
        DADD = - A * DF2DD / (F2-1)
        DF3DD = (2*T + 4*A*T**3) * DTDD + DADD * T**4
        DF4DD = F4 * ( DF3DD/F3 - (DADD*F3+A*DF3DD)/(1+A*F3) )
        DF5DD = (0.023266D0 + 2.0D0*7.389D-6*RS)*DRSDD
        DF6DD = (8.723D0 + 2.0D0*0.472D0*RS
     .            + 3.0D0*0.07389D0*RS**2)*DRSDD
        DF7DD = -200.0D0 * S * PHI**4 * DSDD * F7
     .         -100.0D0 * S**2 * 4.0D0* PHI**3 * DPDD * F7
        DF8DD = 15.75592D0 * DF5DD/F6 - 15.75592D0*F5*DF6DD / F6**2
        DHDD = 3 * H * DPDD / PHI
        DHDD = DHDD + GAMMA * PHI**3 * DF4DD / (1+F4)
        DHDD = DHDD + DF8DD * T**2 * F7
        DHDD = DHDD + F8 * 2*T*DTDD *F7
        DHDD = DHDD + F8 * T**2 * DF7DD
        
        DFCDD(IS) = VCUNIF(IS) + H + DT * DHDD
        DO 30 IX = 1,3
          DTDGD = (T / GDMT) * GDT(IX) / GDMT
          DSDGD = (S / GDMT) * GDT(IX) / GDMT
          DF3DGD = DTDGD * ( 2 * T + 4 * A * T**3 )
          DF4DGD = F4 * DF3DGD * ( 1/F3 - A/(1+A*F3) ) 
          DF7DGD = -200.0D0 * S * PHI**4 * DSDGD * F7
          DHDGD = GAMMA * PHI**3 * DF4DGD / (1+F4)
          DHDGD = DHDGD + F8 * 2*T*DTDGD *F7 + F8 * T**2 *DF7DGD
          DFCDGD(IX,IS) = DT * DHDGD
   30   CONTINUE
   40 CONTINUE

C Find exchange energy and potential
      FX = 0
      DO 60 IS = 1,2
        DS(1) = MAX( DENMIN, 2 * D(IS) )
        GDMS = MAX( GDMIN, 2 * GDM(IS) )
        KFS = (3 * PI**2 * DS(1))**THD
        S = GDMS / (2 * KFS * DS(1))
        F4 = SQRT(1.0D0 + (7.7956D0*S)**2)
        ASINHS = LOG(7.7956D0*S + F4)
        F1 = 1.0D0 + 0.19645D0 * S * ASINHS
        F2 = 0.2743D0 - 0.15084D0*EXP(-100.0D0*S*S)
        F3 = 1.0D0 / (F1 + 0.004D0 * S*S*S*S)
        F = (F1 + F2 * S*S ) * F3
     .       
        CALL EXCHNG( IREL, 1, DS, EXUNIF, VXUNIF )
        FX = FX + DS(1) * EXUNIF * F

        DKFDD = THD * KFS / DS(1)
        DSDD = S * ( -DKFDD/KFS - 1/DS(1) )
        DF1DS = 0.19645D0 * ASINHS +
     .    0.19645D0 * S * 7.7956D0 / F4
        DF2DS = 0.15084D0*200.0D0*S*EXP(-100.0D0*S*S)
        DF3DS = - F3*F3 * (DF1DS + 4.0D0*0.004D0 * S*S*S)
        DFDS =  DF1DS * F3 + DF2DS * S*S * F3 + 2.0D0 * S * F2 * F3
     .            + (F1 + F2 * S*S ) * DF3DS   
        DFXDD(IS) = VXUNIF(1) * F + DS(1) * EXUNIF * DFDS * DSDD

        DO 50 IX = 1,3
          GDS = 2 * GD(IX,IS)
          DSDGD = (S / GDMS) * GDS / GDMS
          DFDGD = DFDS * DSDGD
          DFXDGD(IX,IS) = DS(1) * EXUNIF * DFDGD
   50   CONTINUE
   60 CONTINUE
      FX = HALF * FX / DT

C Set output arguments
      EX = FX
      EC = FC
      DO 90 IS = 1,NSPIN
        DEXDD(IS) = DFXDD(IS)
        DECDD(IS) = DFCDD(IS)
        DO 80 IX = 1,3
          DEXDGD(IX,IS) = DFXDGD(IX,IS)
          DECDGD(IX,IS) = DFCDGD(IX,IS)
   80   CONTINUE
   90 CONTINUE

      END



      SUBROUTINE PW92C( NSPIN, DENS, EC, VC )

C ********************************************************************
C Implements the Perdew-Wang'92 local correlation (beyond RPA).
C Ref: J.P.Perdew & Y.Wang, PRB, 45, 13244 (1992)
C Written by L.C.Balbas and J.M.Soler. Dec'96.  Version 0.5.
C ********* INPUT ****************************************************
C INTEGER NSPIN       : Number of spin polarizations (1 or 2)
C REAL*8  DENS(NSPIN) : Local (spin) density
C ********* OUTPUT ***************************************************
C REAL*8  EC        : Correlation energy density
C REAL*8  VC(NSPIN) : Correlation (spin) potential
C ********* UNITS ****************************************************
C Densities in electrons per Bohr**3
C Energies in Hartrees
C ********* ROUTINES CALLED ******************************************
C None
C ********************************************************************

C Next line is nonstandard but may be supressed
      IMPLICIT          NONE

C Argument types and dimensions
      INTEGER           NSPIN
      DOUBLE PRECISION  DENS(NSPIN), EC, VC(NSPIN)

C Internal variable declarations
      INTEGER           IG
      DOUBLE PRECISION  A(0:2), ALPHA1(0:2), B, BETA(0:2,4), C,
     .                  DBDRS, DECDD(2), DECDRS, DECDZ, DENMIN, DFDZ,
     .                  DGDRS(0:2), DCDRS, DRSDD, DTOT, DZDD(2),
     .                  F, FPP0, FOUTHD, G(0:2), HALF, ONE,
     .                  P(0:2), PI, RS, THD, THRHLF, ZETA

C Add tiny numbers to avoid numerical errors
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( ONE    = 1.D0 + 1.D-12 )

C Fix some numerical constants
      PARAMETER ( FOUTHD=4.D0/3.D0, HALF=0.5D0,
     .            THD=1.D0/3.D0, THRHLF=1.5D0 )

C Parameters from Table I of Perdew & Wang, PRB, 45, 13244 (92)
      DATA P      / 1.00d0,     1.00d0,     1.00d0     /
      DATA A      / 0.031091d0, 0.015545d0, 0.016887d0 /
      DATA ALPHA1 / 0.21370d0,  0.20548d0,  0.11125d0  /
      DATA BETA   / 7.5957d0,  14.1189d0,  10.357d0,
     .              3.5876d0,   6.1977d0,   3.6231d0,
     .              1.6382d0,   3.3662d0,   0.88026d0,
     .              0.49294d0,  0.62517d0,  0.49671d0 /

C Find rs and zeta
      PI = 4 * ATAN(1.D0)
      IF (NSPIN .EQ. 1) THEN
        DTOT = MAX( DENMIN, DENS(1) )
        ZETA = 0
        RS = ( 3 / (4*PI*DTOT) )**THD
C       Find derivatives dRs/dDens and dZeta/dDens
        DRSDD = (- RS) / DTOT / 3
        DZDD(1) = 0
      ELSE
        DTOT = MAX( DENMIN, DENS(1)+DENS(2) )
        ZETA = ( DENS(1) - DENS(2) ) / DTOT
        RS = ( 3 / (4*PI*DTOT) )**THD
        DRSDD = (- RS) / DTOT / 3
        DZDD(1) =   (ONE - ZETA) / DTOT
        DZDD(2) = - (ONE + ZETA) / DTOT
      ENDIF

C Find eps_c(rs,0)=G(0), eps_c(rs,1)=G(1) and -alpha_c(rs)=G(2)
C using eq.(10) of cited reference (Perdew & Wang, PRB, 45, 13244 (92))
      DO 20 IG = 0,2
        B = BETA(IG,1) * RS**HALF   +
     .      BETA(IG,2) * RS         +
     .      BETA(IG,3) * RS**THRHLF +
     .      BETA(IG,4) * RS**(P(IG)+1)
        DBDRS = BETA(IG,1) * HALF      / RS**HALF +
     .          BETA(IG,2)                         +
     .          BETA(IG,3) * THRHLF    * RS**HALF +
     .          BETA(IG,4) * (P(IG)+1) * RS**P(IG)
        C = 1 + 1 / (2 * A(IG) * B)
        DCDRS = - ( (C-1) * DBDRS / B )
        G(IG) = (- 2) * A(IG) * ( 1 + ALPHA1(IG)*RS ) * LOG(C)
        DGDRS(IG) = (- 2) *A(IG) * ( ALPHA1(IG) * LOG(C) +
     .                            (1+ALPHA1(IG)*RS) * DCDRS / C )
   20 CONTINUE

C Find f''(0) and f(zeta) from eq.(9)
      C = 1 / (2**FOUTHD - 2)
      FPP0 = 8 * C / 9
      F = ( (ONE+ZETA)**FOUTHD + (ONE-ZETA)**FOUTHD - 2 ) * C
      DFDZ = FOUTHD * ( (ONE+ZETA)**THD - (ONE-ZETA)**THD ) * C

C Find eps_c(rs,zeta) from eq.(8)
      EC = G(0) - G(2) * F / FPP0 * (ONE-ZETA**4) +
     .    (G(1)-G(0)) * F * ZETA**4
      DECDRS = DGDRS(0) - DGDRS(2) * F / FPP0 * (ONE-ZETA**4) +
     .        (DGDRS(1)-DGDRS(0)) * F * ZETA**4
      DECDZ = (- G(2)) / FPP0 * ( DFDZ*(ONE-ZETA**4) - F*4*ZETA**3 ) +
     .        (G(1)-G(0)) * ( DFDZ*ZETA**4 + F*4*ZETA**3 )
      
C Find correlation potential
      IF (NSPIN .EQ. 1) THEN
        DECDD(1) = DECDRS * DRSDD
        VC(1) = EC + DTOT * DECDD(1)
      ELSE
        DECDD(1) = DECDRS * DRSDD + DECDZ * DZDD(1)
        DECDD(2) = DECDRS * DRSDD + DECDZ * DZDD(2)
        VC(1) = EC + DTOT * DECDD(1)
        VC(2) = EC + DTOT * DECDD(2)
      ENDIF

      END



      SUBROUTINE PW92XC( IREL, NSPIN, DENS, EPSX, EPSC, VX, VC )

C ********************************************************************
C Implements the Perdew-Wang'92 LDA/LSD exchange correlation
C Ref: J.P.Perdew & Y.Wang, PRB, 45, 13244 (1992)
C Written by L.C.Balbas and J.M.Soler. Dec'96. Version 0.5.
C ********* INPUT ****************************************************
C INTEGER IREL        : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER NSPIN       : Number of spin polarizations (1 or 2)
C REAL*8  DENS(NSPIN) : Local (spin) density
C ********* OUTPUT ***************************************************
C REAL*8  EPSX       : Exchange energy density
C REAL*8  EPSC       : Correlation energy density
C REAL*8  VX(NSPIN)  : Exchange (spin) potential
C REAL*8  VC(NSPIN)  : Correlation (spin) potential
C ********* UNITS ****************************************************
C Densities in electrons per Bohr**3
C Energies in Hartrees
C ********* ROUTINES CALLED ******************************************
C EXCHNG, PW92C
C ********************************************************************

      IMPLICIT          NONE
      INTEGER           IREL, NSPIN
      DOUBLE PRECISION  DENS(NSPIN), EPSX, EPSC, VC(NSPIN), VX(NSPIN)

      CALL EXCHNG( IREL, NSPIN, DENS, EPSX, VX )
      CALL PW92C( NSPIN, DENS, EPSC, VC )
      END



      SUBROUTINE PZXC( IREL, NSP, DS, EX, EC, VX, VC, DVXDN, DVCDN )

C *****************************************************************
C  Perdew-Zunger parameterization of Ceperley-Alder exchange and 
C  correlation. Ref: Perdew & Zunger, Phys. Rev. B 23 5075 (1981).
C  Adapted by J.M.Soler from routine velect of Froyen's 
C    pseudopotential generation program. Madrid, Jan'97.
C **** Input *****************************************************
C INTEGER IREL    : relativistic-exchange switch (0=no, 1=yes)
C INTEGER NSP     : spin-polarizations (1=>unpolarized, 2=>polarized)
C REAL*8  DS(NSP) : total (nsp=1) or spin (nsp=2) electron density
C **** Output *****************************************************
C REAL*8  EX            : exchange energy density
C REAL*8  EC            : correlation energy density
C REAL*8  VX(NSP)       : (spin-dependent) exchange potential
C REAL*8  VC(NSP)       : (spin-dependent) correlation potential
C REAL*8  DVXDN(NSP,NSP): Derivative of the exchange potential
C                         respect the charge density, 
C                         Dvx(spin1)/Dn(spin2)
C REAL*8  DVCDN(NSP,NSP): Derivative of the correlation potential
C                         respect the charge density, 
C                         Dvc(spin1)/Dn(spin2)
C **** Units *******************************************************
C Densities in electrons/Bohr**3
C Energies in Hartrees
C *****************************************************************

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION DS(NSP), VX(NSP), VC(NSP), 
     .           DVXDN(NSP,NSP), DVCDN(NSP,NSP)

       PARAMETER (ZERO=0.D0,ONE=1.D0,PFIVE=.5D0,OPF=1.5D0,PNN=.99D0)
       PARAMETER (PTHREE=0.3D0,PSEVF=0.75D0,C0504=0.0504D0) 
       PARAMETER (C0254=0.0254D0,C014=0.014D0,C0406=0.0406D0)
       PARAMETER (C15P9=15.9D0,C0666=0.0666D0,C11P4=11.4D0)
       PARAMETER (C045=0.045D0,C7P8=7.8D0,C88=0.88D0,C20P59=20.592D0)
       PARAMETER (C3P52=3.52D0,C0311=0.0311D0,C0014=0.0014D0)
       PARAMETER (C0538=0.0538D0,C0096=0.0096D0,C096=0.096D0)
       PARAMETER (C0622=0.0622D0,C004=0.004D0,C0232=0.0232D0)
       PARAMETER (C1686=0.1686D0,C1P398=1.3981D0,C2611=0.2611D0)
       PARAMETER (C2846=0.2846D0,C1P053=1.0529D0,C3334=0.3334D0)
Cray       PARAMETER (ZERO=0.0,ONE=1.0,PFIVE=0.5,OPF=1.5,PNN=0.99)
Cray       PARAMETER (PTHREE=0.3,PSEVF=0.75,C0504=0.0504) 
Cray       PARAMETER (C0254=0.0254,C014=0.014,C0406=0.0406)
Cray       PARAMETER (C15P9=15.9,C0666=0.0666,C11P4=11.4)
Cray       PARAMETER (C045=0.045,C7P8=7.8,C88=0.88,C20P59=20.592)
Cray       PARAMETER (C3P52=3.52,C0311=0.0311,C0014=0.0014)
Cray       PARAMETER (C0538=0.0538,C0096=0.0096,C096=0.096)
Cray       PARAMETER (C0622=0.0622,C004=0.004,C0232=0.0232)
Cray       PARAMETER (C1686=0.1686,C1P398=1.3981,C2611=0.2611)
Cray       PARAMETER (C2846=0.2846,C1P053=1.0529,C3334=0.3334)

C    Ceperly-Alder 'ca' constants. Internal energies in Rydbergs.
       PARAMETER (CON1=1.D0/6, CON2=0.008D0/3, CON3=0.3502D0/3) 
       PARAMETER (CON4=0.0504D0/3, CON5=0.0028D0/3, CON6=0.1925D0/3)
       PARAMETER (CON7=0.0206D0/3, CON8=9.7867D0/6, CON9=1.0444D0/3)
       PARAMETER (CON10=7.3703D0/6, CON11=1.3336D0/3)
Cray       PARAMETER (CON1=1.0/6, CON2=0.008/3, CON3=0.3502/3) 
Cray       PARAMETER (CON4=0.0504/3, CON5=0.0028/3, CON6=0.1925/3)
Cray       PARAMETER (CON7=0.0206/3, CON8=9.7867/6, CON9=1.0444/3)
Cray       PARAMETER (CON10=7.3703/6, CON11=1.3336/3) 

C      X-alpha parameter:
       PARAMETER ( ALP = 2.D0 / 3.D0 )

C      Other variables converted into parameters by J.M.Soler
       PARAMETER ( TINY = 1.D-6 )
       PARAMETER ( PI   = 3.14159265358979312D0 )
       PARAMETER ( TWO  = 2.0D0 ) 
       PARAMETER ( HALF = 0.5D0 ) 
       PARAMETER ( TRD  = 1.D0 / 3.D0 ) 
       PARAMETER ( FTRD = 4.D0 / 3.D0 )
       PARAMETER ( TFTM = 0.51984209978974638D0 )
       PARAMETER ( A0   = 0.52106176119784808D0 )
       PARAMETER ( CRS  = 0.620350490899400087D0 )
       PARAMETER ( CXP  = (- 3.D0) * ALP / (PI*A0) )
       PARAMETER ( CXF  = 1.25992104989487319D0 )

C      Find density and polarization
       IF (NSP .EQ. 2) THEN
         D1 = MAX(DS(1),ZERO)
         D2 = MAX(DS(2),ZERO)
         D = D1 + D2
         IF (D .LE. ZERO) THEN
           EX = ZERO
           EC = ZERO
           VX(1) = ZERO
           VX(2) = ZERO
           VC(1) = ZERO
           VC(2) = ZERO
           RETURN
         ENDIF
c
c        Robustness enhancement by Jose Soler (August 2002)
c
         Z = (D1 - D2) / D
         IF (Z .LE. -ONE) THEN
           FZ = (TWO**FTRD-TWO)/TFTM
           FZP = -FTRD*TWO**TRD/TFTM
           DFZPDN = FTRD*TRD*TWO**(-ALP)/TFTM
         ELSEIF (Z .GE. ONE) THEN
           FZ = (TWO**FTRD-TWO)/TFTM
           FZP = FTRD*TWO**TRD/TFTM
           DFZPDN = FTRD*TRD*TWO**(-ALP)/TFTM
         ELSE
           FZ = ((ONE+Z)**FTRD+(ONE-Z)**FTRD-TWO)/TFTM
           FZP = FTRD*((ONE+Z)**TRD-(ONE-Z)**TRD)/TFTM 
           DFZPDN = FTRD*TRD*((ONE+Z)**(-ALP) + (ONE-Z)**(-ALP))/TFTM
         ENDIF
       ELSE
         D = DS(1)
         IF (D .LE. ZERO) THEN
           EX = ZERO
           EC = ZERO
           VX(1) = ZERO
           VC(1) = ZERO
           RETURN
         ENDIF
         Z = ZERO
         FZ = ZERO
         FZP = ZERO
       ENDIF
       RS = CRS / D**TRD

C      Exchange
       VXP = CXP / RS
       EXP_VAR = 0.75D0 * VXP
       IF (IREL .EQ. 1) THEN
         BETA = C014/RS
         IF (BETA .LT. TINY) THEN
           SB = ONE + HALF*BETA**2
           ALB = BETA
         ELSE
           SB = SQRT(1+BETA*BETA)
           ALB = LOG(BETA+SB)
         ENDIF
         VXP = VXP * (-PFIVE + OPF * ALB / (BETA*SB))
         EXP_VAR = EXP_VAR *(ONE-OPF*((BETA*SB-ALB)/BETA**2)**2) 
       ENDIF
       VXF = CXF * VXP
       EXF = CXF * EXP_VAR
       DVXPDN = TRD * VXP / D
       DVXFDN = TRD * VXF / D

C      Correlation 
       IF (RS .GT. ONE) THEN  
         SQRS=SQRT(RS)
         TE = ONE+CON10*SQRS+CON11*RS
         BE = ONE+C1P053*SQRS+C3334*RS
         ECP = -(C2846/BE)
         VCP = ECP*TE/BE
         DTEDN = ((CON10 * SQRS *HALF) + CON11 * RS)*(-TRD/D)
         BE2 = BE * BE
         DBEDN = ((C1P053 * SQRS *HALF) + C3334 * RS)*(-TRD/D)
         DVCPDN = -(C2846/BE2)*(DTEDN - 2.0D0 * TE * DBEDN/BE)
         DECPDN = (C2846/BE2)*DBEDN
         TE = ONE+CON8*SQRS+CON9*RS
         BE = ONE+C1P398*SQRS+C2611*RS
         ECF = -(C1686/BE)
         VCF = ECF*TE/BE
         DTEDN = ((CON8 * SQRS * HALF) + CON9 * RS)*(-TRD/D)
         BE2 = BE * BE
         DBEDN = ((C1P398 * SQRS * HALF) + C2611 * RS)*(-TRD/D)
         DVCFDN = -(C1686/BE2)*(DTEDN - 2.0D0 * TE * DBEDN/BE)
         DECFDN = (C1686/BE2)*DBEDN
       ELSE
         RSLOG=LOG(RS)
         ECP=(C0622+C004*RS)*RSLOG-C096-C0232*RS
         VCP=(C0622+CON2*RS)*RSLOG-CON3-CON4*RS
         DVCPDN = (CON2*RS*RSLOG + (CON2-CON4)*RS + C0622)*(-TRD/D)
         DECPDN = (C004*RS*RSLOG + (C004-C0232)*RS + C0622)*(-TRD/D)
         ECF=(C0311+C0014*RS)*RSLOG-C0538-C0096*RS
         VCF=(C0311+CON5*RS)*RSLOG-CON6-CON7*RS
         DVCFDN = (CON5*RS*RSLOG + (CON5-CON7)*RS + C0311)*(-TRD/D)
         DECFDN = (C0014*RS*RSLOG + (C0014-C0096)*RS + C0311)*(-TRD/D)
       ENDIF

       ISP1 = 1
       ISP2 = 2

C      Find up and down potentials
       IF (NSP .EQ. 2) THEN
         EX    = EXP_VAR + FZ*(EXF-EXP_VAR)
         EC    = ECP + FZ*(ECF-ECP)
         VX(1) = VXP + FZ*(VXF-VXP) + (ONE-Z)*FZP*(EXF-EXP_VAR)
         VX(2) = VXP + FZ*(VXF-VXP) - (ONE+Z)*FZP*(EXF-EXP_VAR)
         VC(1) = VCP + FZ*(VCF-VCP) + (ONE-Z)*FZP*(ECF-ECP)
         VC(2) = VCP + FZ*(VCF-VCP) - (ONE+Z)*FZP*(ECF-ECP)

C        Derivatives of exchange potential respect the density

         DVXDN(ISP1,ISP1) =
     .             DVXPDN
     .              +  FZP*(VXF-VXP-EXF+EXP_VAR)*( 2.D0*D2/(D*D) )
     .              +  FZ*(DVXFDN-DVXPDN)+(1-Z)*FZP*(VXF-VXP)/(4.D0*D)
     .              +  (1-Z)*DFZPDN*(EXF-EXP_VAR)*( 2.D0*D2/(D*D) )
         DVXDN(ISP1,ISP2) =
     .                 DVXPDN
     .              +  FZP*(VXF-VXP-EXF+EXP_VAR)*(-2.D0*D1/(D*D) )
     .              +  FZ*(DVXFDN-DVXPDN)+(1-Z)*FZP*(VXF-VXP)/(4.D0*D)
     .              +  (1-Z)*DFZPDN*(EXF-EXP_VAR)*( -2.D0*D1/(D*D) )
         DVXDN(ISP2,ISP1) =
     .                 DVXPDN
     .              +  FZP*(VXF-VXP-EXF+EXP_VAR)*( 2.D0*D2/(D*D) )
     .              +  FZ*(DVXFDN-DVXPDN)-(1+Z)*FZP*(VXF-VXP)/(4.D0*D)
     .              -  (1+Z)*DFZPDN*(EXF-EXP_VAR)*( 2.D0*D2/(D*D) )
         DVXDN(ISP2,ISP2) =
     .                 DVXPDN
     .              +  FZP*(VXF-VXP-EXF+EXP_VAR)*(-2.D0*D1/(D*D) )
     .              +  FZ*(DVXFDN-DVXPDN)-(1+Z)*FZP*(VXF-VXP)/(4.D0*D)
     .              -  (1+Z)*DFZPDN*(EXF-EXP_VAR)*( -2.D0*D1/(D*D) )

C        Derivatives of correlation potential respect the density

         DVCDN(ISP1,ISP1) =
     .                DVCPDN
     .              + FZP*(VCF-VCP-ECF+ECP)*( 2.D0*D2/(D*D) )
     .              + FZ*(DVCFDN-DVCPDN)+ (1-Z)*FZP*(DECFDN-DECPDN)
     .              + (1-Z)*DFZPDN*(ECF-ECP)*( 2.D0*D2/(D*D) )
         DVCDN(ISP1,ISP2) =
     .                DVCPDN
     .              + FZP*(VCF-VCP-ECF+ECP)*(-2.D0*D1/(D*D) )
     .              + FZ*(DVCFDN-DVCPDN)+ (1-Z)*FZP*(DECFDN-DECPDN)
     .              + (1-Z)*DFZPDN*(ECF-ECP)*( -2.D0*D1/(D*D) )
         DVCDN(ISP2,ISP1) =
     .                DVCPDN
     .              + FZP*(VCF-VCP-ECF+ECP)*( 2.D0*D2/(D*D) )
     .              + FZ*(DVCFDN-DVCPDN)- (1+Z)*FZP*(DECFDN-DECPDN)
     .              - (1+Z)*DFZPDN*(ECF-ECP)*( 2.D0*D2/(D*D) )
         DVCDN(ISP2,ISP2) =
     .                DVCPDN
     .              + FZP*(VCF-VCP-ECF+ECP)*(-2.D0*D1/(D*D) )
     .              + FZ*(DVCFDN-DVCPDN)- (1+Z)*FZP*(DECFDN-DECPDN)
     .              - (1+Z)*DFZPDN*(ECF-ECP)*( -2.D0*D1/(D*D) )

       ELSE
         EX    = EXP_VAR
         EC    = ECP
         VX(1) = VXP
         VC(1) = VCP
         DVXDN(1,1) = DVXPDN
         DVCDN(1,1) = DVCPDN
       ENDIF

C      Change from Rydbergs to Hartrees
       EX = HALF * EX
       EC = HALF * EC
       DO 10 ISP = 1,NSP
         VX(ISP) = HALF * VX(ISP)
         VC(ISP) = HALF * VC(ISP)
         DO 5 ISP2 = 1,NSP
           DVXDN(ISP,ISP2) = HALF * DVXDN(ISP,ISP2)
           DVCDN(ISP,ISP2) = HALF * DVCDN(ISP,ISP2)
    5    CONTINUE
   10  CONTINUE
      END

       subroutine blypxc(nspin,dens,gdens,EX,EC,
     .                   dEXdd,dECdd,dEXdgd,dECdgd) 
c ***************************************************************
c Implements Becke gradient exchange functional (A.D. 
c Becke, Phys. Rev. A 38, 3098 (1988)) and Lee, Yang, Parr
c correlation functional (C. Lee, W. Yang, R.G. Parr, Phys. Rev. B
c 37, 785 (1988)), as modificated by Miehlich,Savin,Stoll and Preuss,
c Chem. Phys. Lett. 157,200 (1989). See also Johnson, Gill and Pople,
c J. Chem. Phys. 98, 5612 (1993). Some errors were detected in this
c last paper, so not all of the expressions correspond exactly to those
c implemented here.
c Written by Maider Machado. July 1998.
c **************** INPUT ******************************************** 
c integer nspin          : Number of spin polarizations (1 or 2)
c real*8  dens(nspin)    : Total electron density (if nspin=1) or
c                           spin electron density (if nspin=2)
c real*8  gdens(3,nspin) : Total or spin density gradient
c ******** OUTPUT *****************************************************
c real*8  ex             : Exchange energy density
c real*8  ec             : Correlation energy density
c real*8  dexdd(nspin)   : Partial derivative
c                           d(DensTot*Ex)/dDens(ispin),
c                           where DensTot = Sum_ispin( DENS(ispin) )
c                          For a constant density, this is the
c                          exchange potential
c real*8  decdd(nspin)   : Partial derivative
c                           d(DensTot*Ec)/dDens(ispin),
c                           where DensTot = Sum_ispin( DENS(ispin) )
c                          For a constant density, this is the
c                          correlation potential
c real*8  dexdgd(3,nspin): Partial derivative
c                           d(DensTot*Ex)/d(GradDens(i,ispin))
c real*8  decdgd(3,nspin): Partial derivative
c                           d(DensTot*Ec)/d(GradDens(i,ispin))
c ********* UNITS ****************************************************
c Lengths in Bohr
c Densities in electrons per Bohr**3
c Energies in Hartrees
c Gradient vectors in cartesian coordinates
c ********************************************************************
 
      implicit none
      integer nspin
      double precision  dens(nspin), gdens(3,nspin), EX, EC,
     .                  dEXdd(nspin), dECdd(nspin), dEXdgd(3,nspin),
     .                  dECdgd(3,nspin)

c Internal variables
      integer is,ix,ois
      double precision pi, beta, thd, tthd, thrhlf, half, fothd,
     .                 d(2),gd(3,2),dmin, ash,gdm(2),denmin,dt, 
     .                 g(2),x(2),a,b,c,dd,onzthd,gdmin, 	     
     .                 ga, gb, gc,becke,dbecgd(3,2),
     .                 dgdx(2), dgdxa, dgdxb, dgdxc,dgdxd,dbecdd(2),
     .                 den,omega, domega, delta, ddelta,cf,
     .                 gam11, gam12, gam22, LYPa, LYPb1,
     .                 LYPb2,dLYP11,dLYP12,dLYP22,LYP,
     .                 dd1g11,dd1g12,dd1g22,dd2g12,dd2g11,dd2g22,
     .                 dLYPdd(2),dg11dd(3,2),dg22dd(3,2),
     .                 dg12dd(3,2),dLYPgd(3,2)
  
c Lower bounds of density and its gradient to avoid divisions by zero
      parameter ( denmin=1.d-8 )
      parameter (gdmin=1.d-8)
      parameter (dmin=1.d-5)

c Fix some numerical parameters 
      parameter ( thd = 1.d0/3.d0, tthd=2.d0/3.d0 )
      parameter ( thrhlf=1.5d0, half=0.5d0,
     .            fothd=4.d0/3.d0, onzthd=11.d0/3.d0)

c Empirical parameter for Becke exchange functional (a.u.)
      parameter(beta= 0.0042d0) 

c Constants for LYP functional (a.u.) 
      parameter(a=0.04918d0, b=0.132d0, c=0.2533d0, dd=0.349d0)

       pi= 4*atan(1.d0)
       

c Translate density and its gradient to new variables
      if (nspin .eq. 1) then
        d(1) = half * dens(1)
        d(1) = max(denmin,d(1))
        d(2) = d(1)
        dt = max( denmin, dens(1) )
        do ix = 1,3
          gd(ix,1) = half * gdens(ix,1)    
          gd(ix,2) = gd(ix,1)
        enddo 
      else
        d(1) = dens(1)
        d(2) = dens(2)
        do is=1,2
         d(is) = max (denmin,d(is))
        enddo
        dt = max( denmin, dens(1)+dens(2) )  
        do ix = 1,3
          gd(ix,1) = gdens(ix,1)
          gd(ix,2) = gdens(ix,2)
        enddo
      endif

      gdm(1) = sqrt( gd(1,1)**2 + gd(2,1)**2 + gd(3,1)**2 )
      gdm(2) = sqrt( gd(1,2)**2 + gd(2,2)**2 + gd(3,2)**2 )
 
      do is=1,2
      gdm(is)= max(gdm(is),gdmin)
      enddo

c Find Becke exchange energy
       ga = -thrhlf*(3.d0/4.d0/pi)**thd
      do is=1,2
       if(d(is).lt.dmin) then
        g(is)=ga
       else
        x(is) = gdm(is)/d(is)**fothd
        gb = beta*x(is)**2
        ash=log(x(is)+sqrt(x(is)**2+1)) 
        gc = 1+6*beta*x(is)*ash        
        g(is) = ga-gb/gc
       endif
      enddo

c   Density of energy 
      becke=(g(1)*d(1)**fothd+g(2)*d(2)**fothd)/dt

      
c Exchange energy derivatives
       do is=1,2
        if(d(is).lt.dmin)then
         dbecdd(is)=0.
         do ix=1,3
          dbecgd(ix,is)=0.
         enddo
        else
        dgdxa=6*beta**2*x(is)**2
        ash=log(x(is)+sqrt(x(is)**2+1))
        dgdxb=x(is)/sqrt(x(is)**2+1)-ash
        dgdxc=-2*beta*x(is)
        dgdxd=(1+6*beta*x(is)*ash)**2
        dgdx(is)=(dgdxa*dgdxb+dgdxc)/dgdxd
        dbecdd(is)=fothd*d(is)**thd*(g(is)-x(is)*dgdx(is))
        do ix=1,3
         dbecgd(ix,is)=d(is)**(-fothd)*dgdx(is)*gd(ix,is)/x(is)
        enddo 
        endif
       enddo

c  Lee-Yang-Parr correlation energy
      den=1+dd*dt**(-thd)
      omega=dt**(-onzthd)*exp(-c*dt**(-thd))/den
      delta=c*dt**(-thd)+dd*dt**(-thd)/den
      cf=3.*(3*pi**2)**tthd/10.
      gam11=gdm(1)**2
      gam12=gd(1,1)*gd(1,2)+gd(2,1)*gd(2,2)+gd(3,1)*gd(3,2)
      gam22=gdm(2)**2
      LYPa=-4*a*d(1)*d(2)/(den*dt)
      LYPb1=2**onzthd*cf*a*b*omega*d(1)*d(2)
      LYPb2=d(1)**(8./3.)+d(2)**(8./3.)
      dLYP11=-a*b*omega*(d(1)*d(2)/9.*(1.-3.*delta-(delta-11.)
     .*d(1)/dt)-d(2)**2)
      dLYP12=-a*b*omega*(d(1)*d(2)/9.*(47.-7.*delta)
     .-fothd*dt**2)
      dLYP22=-a*b*omega*(d(1)*d(2)/9.*(1.-3.*delta-(delta-11.)*
     .d(2)/dt)-d(1)**2)

c    Density of energy
      LYP=(LYPa-LYPb1*LYPb2+dLYP11*gam11+dLYP12*gam12
     .+dLYP22*gam22)/dt

c   Correlation energy derivatives
       domega=-thd*dt**(-fothd)*omega*(11.*dt**thd-c-dd/den)
       ddelta=thd*(dd**2*dt**(-5./3.)/den**2-delta/dt)

c   Second derivatives with respect to the density
       dd1g11=domega/omega*dLYP11-a*b*omega*(d(2)/9.*
     . (1.-3.*delta-2*(delta-11.)*d(1)/dt)-d(1)*d(2)/9.*
     . ((3.+d(1)/dt)*ddelta-(delta-11.)*d(1)/dt**2))

       dd1g12=domega/omega*dLYP12-a*b*omega*(d(2)/9.*
     . (47.-7.*delta)-7./9.*d(1)*d(2)*ddelta-8./3.*dt)

      dd1g22=domega/omega*dLYP22-a*b*omega*(1./9.*d(2)
     . *(1.-3.*delta-(delta-11.)*d(2)/dt)-d(1)*d(2)/9.*
     . ((3.+d(2)/dt)*ddelta-(delta-11.)*d(2)/dt**2)-2*d(1))

       
      dd2g22=domega/omega*dLYP22-a*b*omega*(d(1)/9.*
     . (1.-3.*delta-2*(delta-11.)*d(2)/dt)-d(1)*d(2)/9.*
     . ((3+d(2)/dt)*ddelta-(delta-11.)*d(2)/dt**2))
      
 
      dd2g12=domega/omega*dLYP12-a*b*omega*(d(1)/9.*
     . (47.-7.*delta)-7./9.*d(1)*d(2)*ddelta-8./3.*dt)
      
      dd2g11=domega/omega*dLYP11-a*b*omega*(1./9.*d(1)
     . *(1.-3.*delta-(delta-11.)*d(1)/dt)-d(1)*d(2)/9.*
     . ((3.+d(1)/dt)*ddelta-(delta-11.)*d(1)/dt**2)-2*d(2))


        dLYPdd(1)=-4*a/den*d(1)*d(2)/dt*
     . (thd*dd*dt**(-fothd)/den
     . +1./d(1)-1./dt)-2**onzthd*cf*a*b*(domega*d(1)*d(2)*
     . (d(1)**(8./3.)+d(2)**(8./3.))+omega*d(2)*(onzthd*
     . d(1)**(8./3.)+d(2)**(8./3.)))+dd1g11*gam11+
     . dd1g12*gam12+dd1g22*gam22


       dLYPdd(2)=-4*a/den*d(1)*d(2)/dt*(thd*dd*dt**(-fothd)/den
     . +1./d(2)-1./dt)-2**onzthd*cf*a*b*(domega*d(1)*d(2)*
     . (d(1)**(8./3.)+d(2)**(8./3.))+omega*d(1)*(onzthd*
     . d(2)**(8./3.)+d(1)**(8./3.)))+dd2g22*gam22+
     . dd2g12*gam12+dd2g11*gam11


c   Second derivatives with respect to the density gradient
       do is=1,2
        if (is.eq.1)then
          ois=2
        else
          ois=1
        endif
        do ix=1,3
         dg11dd(ix,is)=2*gd(ix,is)
         dg22dd(ix,is)=2*gd(ix,is)
         dg12dd(ix,is)=gd(ix,ois)
         dLYPgd(ix,is)=dLYP11*dg11dd(ix,is)+dLYP12*dg12dd(ix,is)+
     .   dLYP22*dg22dd(ix,is)
        enddo
       enddo

c    Set output arguments
       EX=becke
       EC=LYP
       do is=1,nspin
        dEXdd(is)=dbecdd(is)
        dECdd(is)=dLYPdd(is)
        do ix=1,3
         dEXdgd(ix,is)=dbecgd(ix,is)
         dECdgd(ix,is)=dLYPgd(ix,is)
        enddo
       enddo
       end 
C
C     The following routine was taken from 
C     siesta@uam.es--2004/siesta-devel--reference--1.4--patch-46
C     (i.e., before the changeover to aux_xc by J. Gale)
C
      SUBROUTINE RPBEXC( IREL, nspin, Dens, GDens,
     .                   EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Hammer's RPBE Generalized-Gradient-Approximation (GGA).
C A revision of PBE (Perdew-Burke-Ernzerhof)
C Ref: Hammer, Hansen & Norskov, PRB 59, 7413 (1999) and
C J.P.Perdew, K.Burke & M.Ernzerhof, PRL 77, 3865 (1996)
C
C Written by M.V. Fernandez-Serra. March 2004. On the PBE routine of
C L.C.Balbas and J.M.Soler. December 1996. Version 0.5.
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin          : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(nspin)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,nspin): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C EXCHNG, PW92C
C ********************************************************************

      implicit          none
      INTEGER           IREL, nspin
      double precision  Dens(nspin), DECDD(nspin), DECDGD(3,nspin),
     .                  DEXDD(nspin), DEXDGD(3,nspin), GDens(3,nspin)

C Internal variables
      INTEGER
     .  IS, IX

      double precision
     .  A, BETA, D(2), DADD, DECUDD, DENMIN,
     .  DF1DD, DF2DD, DF3DD, DF4DD, DF1DGD, DF3DGD, DF4DGD,
     .  DFCDD(2), DFCDGD(3,2), DFDD, DFDGD, DFXDD(2), DFXDGD(3,2),
     .  DHDD, DHDGD, DKFDD, DKSDD, DPDD, DPDZ, DRSDD,
     .  DS(2), DSDD, DSDGD, DT, DTDD, DTDGD, DZDD(2),
     .  EC, ECUNIF, EX, EXUNIF,
     .  F, F1, F2, F3, F4, FC, FX, FOUTHD,
     .  GAMMA, GD(3,2), GDM(2), GDMIN, GDMS, GDMT, GDS, GDT(3),
     .  H, HALF, KAPPA, KF, KFS, KS, MU, PHI, PI, RS, S,
     .  T, THD, THRHLF, TWO, TWOTHD, VCUNIF(2), VXUNIF(2), ZETA

C Lower bounds of density and its gradient to avoid divisions by zero
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( GDMIN  = 1.D-12 )

C Fix some numerical parameters
      PARAMETER ( FOUTHD=4.D0/3.D0, HALF=0.5D0,
     .            THD=1.D0/3.D0, THRHLF=1.5D0,
     .            TWO=2.D0, TWOTHD=2.D0/3.D0 )

C Fix some more numerical constants
      PI = 4 * ATAN(1.D0)
      BETA = 0.066725D0
      GAMMA = (1 - LOG(TWO)) / PI**2
      MU = BETA * PI**2 / 3
      KAPPA = 0.804D0

C Translate density and its gradient to new variables
      IF (nspin .EQ. 1) THEN
        D(1) = HALF * Dens(1)
        D(2) = D(1)
        DT = MAX( DENMIN, Dens(1) )
        DO 10 IX = 1,3
          GD(IX,1) = HALF * GDens(IX,1)
          GD(IX,2) = GD(IX,1)
          GDT(IX) = GDens(IX,1)
   10   CONTINUE
      ELSE
        D(1) = Dens(1)
        D(2) = Dens(2)
        DT = MAX( DENMIN, Dens(1)+Dens(2) )
        DO 20 IX = 1,3
          GD(IX,1) = GDens(IX,1)
          GD(IX,2) = GDens(IX,2)
          GDT(IX) = GDens(IX,1) + GDens(IX,2)
   20   CONTINUE
      ENDIF
      GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
      GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
      GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
      GDMT = MAX( GDMIN, GDMT )

C Find local correlation energy and potential
      CALL PW92C( 2, D, ECUNIF, VCUNIF )

C Find total correlation energy
      RS = ( 3 / (4*PI*DT) )**THD
      KF = (3 * PI**2 * DT)**THD
      KS = SQRT( 4 * KF / PI )
      ZETA = ( D(1) - D(2) ) / DT
      ZETA = MAX( -1.D0+DENMIN, ZETA )
      ZETA = MIN(  1.D0-DENMIN, ZETA )
      PHI = HALF * ( (1+ZETA)**TWOTHD + (1-ZETA)**TWOTHD )
      T = GDMT / (2 * PHI * KS * DT)
      F1 = ECUNIF / GAMMA / PHI**3
      F2 = EXP(-F1)
      A = BETA / GAMMA / (F2-1)
      F3 = T**2 + A * T**4
      F4 = BETA/GAMMA * F3 / (1 + A*F3)
      H = GAMMA * PHI**3 * LOG( 1 + F4 )
      FC = ECUNIF + H

C Find correlation energy derivatives
      DRSDD = - (THD * RS / DT)
      DKFDD =   THD * KF / DT
      DKSDD = HALF * KS * DKFDD / KF
      DZDD(1) =   1 / DT - ZETA / DT
      DZDD(2) = - (1 / DT) - ZETA / DT
      DPDZ = HALF * TWOTHD * ( 1/(1+ZETA)**THD - 1/(1-ZETA)**THD )
      DO 40 IS = 1,2
        DECUDD = ( VCUNIF(IS) - ECUNIF ) / DT
        DPDD = DPDZ * DZDD(IS)
        DTDD = (- T) * ( DPDD/PHI + DKSDD/KS + 1/DT )
        DF1DD = F1 * ( DECUDD/ECUNIF - 3*DPDD/PHI )
        DF2DD = (- F2) * DF1DD
        DADD = (- A) * DF2DD / (F2-1)
        DF3DD = (2*T + 4*A*T**3) * DTDD + DADD * T**4
        DF4DD = F4 * ( DF3DD/F3 - (DADD*F3+A*DF3DD)/(1+A*F3) )
        DHDD = 3 * H * DPDD / PHI
        DHDD = DHDD + GAMMA * PHI**3 * DF4DD / (1+F4)
        DFCDD(IS) = VCUNIF(IS) + H + DT * DHDD

        DO 30 IX = 1,3
          DTDGD = (T / GDMT) * GDT(IX) / GDMT
          DF3DGD = DTDGD * ( 2 * T + 4 * A * T**3 )
          DF4DGD = F4 * DF3DGD * ( 1/F3 - A/(1+A*F3) )
          DHDGD = GAMMA * PHI**3 * DF4DGD / (1+F4)
          DFCDGD(IX,IS) = DT * DHDGD
   30   CONTINUE
   40 CONTINUE

C Find exchange energy and potential
      FX = 0
      DO 60 IS = 1,2
        DS(IS)   = MAX( DENMIN, 2 * D(IS) )
        GDMS = MAX( GDMIN, 2 * GDM(IS) )
        KFS = (3 * PI**2 * DS(IS))**THD
        S = GDMS / (2 * KFS * DS(IS))
cea Hammer's RPBE (Hammer, Hansen & Norskov PRB 59 7413 (99)
cea     F1 = DEXP( - MU * S**2 / KAPPA)
cea     F = 1 + KAPPA * (1 - F1)
cea Following is standard PBE
cea     F1 = 1 + MU * S**2 / KAPPA
cea     F = 1 + KAPPA - KAPPA / F1
cea (If revPBE Zhang & Yang, PRL 80,890(1998),change PBE's KAPPA to 1.245)
        F1 = DEXP( - MU * S**2 / KAPPA)
        F = 1 + KAPPA * (1 - F1)

c       Note nspin=1 in call to exchng...

        CALL EXCHNG( IREL, 1, DS(IS), EXUNIF, VXUNIF(IS) )
        FX = FX + DS(IS) * EXUNIF * F

cMVFS   The derivatives of F  also need to be changed for Hammer's RPBE.
cMVFS   DF1DD = 2 * F1 * DSDD  * ( - MU * S / KAPPA)
cMVFS   DF1DGD= 2 * F1 * DSDGD * ( - MU * S / KAPPA)
cMVFS   DFDD  = -1 * KAPPA * DF1DD
cMVFS   DFDGD = -1 * KAPPA * DFDGD

        DKFDD = THD * KFS / DS(IS)
        DSDD = S * ( -(DKFDD/KFS) - 1/DS(IS) )
c       DF1DD = 2 * (F1-1) * DSDD / S
c       DFDD = KAPPA * DF1DD / F1**2
        DF1DD = 2* F1 * DSDD * ( - MU * S / KAPPA)
        DFDD = -1 * KAPPA * DF1DD
        DFXDD(IS) = VXUNIF(IS) * F + DS(IS) * EXUNIF * DFDD

        DO 50 IX = 1,3
          GDS = 2 * GD(IX,IS)
          DSDGD = (S / GDMS) * GDS / GDMS
c         DF1DGD = 2 * MU * S * DSDGD / KAPPA
c         DFDGD = KAPPA * DF1DGD / F1**2
          DF1DGD =2*F1 * DSDGD * ( - MU * S / KAPPA)
          DFDGD = -1 * KAPPA * DF1DGD
          DFXDGD(IX,IS) = DS(IS) * EXUNIF * DFDGD
   50   CONTINUE
   60 CONTINUE
      FX = HALF * FX / DT

C Set output arguments
      EX = FX
      EC = FC
      DO 90 IS = 1,nspin
        DEXDD(IS) = DFXDD(IS)
        DECDD(IS) = DFCDD(IS)
        DO 80 IX = 1,3
          DEXDGD(IX,IS) = DFXDGD(IX,IS)
          DECDGD(IX,IS) = DFCDGD(IX,IS)
   80   CONTINUE
   90 CONTINUE

      END
C
C     The following routine is merely PBEXC above with
C     a change in KAPPA.
C
      SUBROUTINE REVPBEXC( IREL, NSPIN, DENS, GDENS,
     .                  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Perdew-Burke-Ernzerhof Generalized-Gradient-Approximation.
C Ref: J.P.Perdew, K.Burke & M.Ernzerhof, PRL 77, 3865 (1996)
C MODIFIED by Zhang and Yang, PRL 80,890(1998)
C (By just changing the value of KAPPA in PBE...)
C Written by L.C.Balbas and J.M.Soler. December 1996. Version 0.5.
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER NSPIN          : Number of spin polarizations (1 or 2)
C REAL*8  DENS(NSPIN)    : Total electron density (if NSPIN=1) or
C                           spin electron density (if NSPIN=2)
C REAL*8  GDENS(3,NSPIN) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(NSPIN)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( DENS(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(NSPIN)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( DENS(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,NSPIN): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,NSPIN): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C EXCHNG, PW92C
C ********************************************************************

      IMPLICIT          NONE
      INTEGER           IREL, NSPIN
      DOUBLE PRECISION  DENS(NSPIN), DECDD(NSPIN), DECDGD(3,NSPIN),
     .                  DEXDD(NSPIN), DEXDGD(3,NSPIN), GDENS(3,NSPIN)

C Internal variables
      INTEGER
     .  IS, IX

      DOUBLE PRECISION
     .  A, BETA, D(2), DADD, DECUDD, DENMIN,
     .  DF1DD, DF2DD, DF3DD, DF4DD, DF1DGD, DF3DGD, DF4DGD,
     .  DFCDD(2), DFCDGD(3,2), DFDD, DFDGD, DFXDD(2), DFXDGD(3,2),
     .  DHDD, DHDGD, DKFDD, DKSDD, DPDD, DPDZ, DRSDD,
     .  DS(2), DSDD, DSDGD, DT, DTDD, DTDGD, DZDD(2),
     .  EC, ECUNIF, EX, EXUNIF,
     .  F, F1, F2, F3, F4, FC, FX, FOUTHD,
     .  GAMMA, GD(3,2), GDM(2), GDMIN, GDMS, GDMT, GDS, GDT(3),
     .  H, HALF, KAPPA, KF, KFS, KS, MU, PHI, PI, RS, S,
     .  T, THD, THRHLF, TWO, TWOTHD, VCUNIF(2), VXUNIF(2), ZETA

C Lower bounds of density and its gradient to avoid divisions by zero
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( GDMIN  = 1.D-12 )

C Fix some numerical parameters
      PARAMETER ( FOUTHD=4.D0/3.D0, HALF=0.5D0,
     .            THD=1.D0/3.D0, THRHLF=1.5D0,
     .            TWO=2.D0, TWOTHD=2.D0/3.D0 )

C Fix some more numerical constants
      PI = 4 * ATAN(1.D0)
      BETA = 0.066725D0
      GAMMA = (1 - LOG(TWO)) / PI**2
      MU = BETA * PI**2 / 3
      KAPPA = 1.245d0

C Translate density and its gradient to new variables
      IF (NSPIN .EQ. 1) THEN
        D(1) = HALF * DENS(1)
        D(2) = D(1)
        DT = MAX( DENMIN, DENS(1) )
        DO 10 IX = 1,3
          GD(IX,1) = HALF * GDENS(IX,1)
          GD(IX,2) = GD(IX,1)
          GDT(IX) = GDENS(IX,1)
   10   CONTINUE
      ELSE
        D(1) = DENS(1)
        D(2) = DENS(2)
        DT = MAX( DENMIN, DENS(1)+DENS(2) )
        DO 20 IX = 1,3
          GD(IX,1) = GDENS(IX,1)
          GD(IX,2) = GDENS(IX,2)
          GDT(IX) = GDENS(IX,1) + GDENS(IX,2)
   20   CONTINUE
      ENDIF
      GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
      GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
      GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
      GDMT = MAX( GDMIN, GDMT )

C Find local correlation energy and potential
      CALL PW92C( 2, D, ECUNIF, VCUNIF )

C Find total correlation energy
      RS = ( 3 / (4*PI*DT) )**THD
      KF = (3 * PI**2 * DT)**THD
      KS = SQRT( 4 * KF / PI )
      ZETA = ( D(1) - D(2) ) / DT
      ZETA = MAX( -1.D0+DENMIN, ZETA )
      ZETA = MIN(  1.D0-DENMIN, ZETA )
      PHI = HALF * ( (1+ZETA)**TWOTHD + (1-ZETA)**TWOTHD )
      T = GDMT / (2 * PHI * KS * DT)
      F1 = ECUNIF / GAMMA / PHI**3
      F2 = EXP(-F1)
      A = BETA / GAMMA / (F2-1)
      F3 = T**2 + A * T**4
      F4 = BETA/GAMMA * F3 / (1 + A*F3)
      H = GAMMA * PHI**3 * LOG( 1 + F4 )
      FC = ECUNIF + H

C Find correlation energy derivatives
      DRSDD = - (THD * RS / DT)
      DKFDD =   THD * KF / DT
      DKSDD = HALF * KS * DKFDD / KF
      DZDD(1) =   1 / DT - ZETA / DT
      DZDD(2) = - (1 / DT) - ZETA / DT
      DPDZ = HALF * TWOTHD * ( 1/(1+ZETA)**THD - 1/(1-ZETA)**THD )
      DO 40 IS = 1,2
        DECUDD = ( VCUNIF(IS) - ECUNIF ) / DT
        DPDD = DPDZ * DZDD(IS)
        DTDD = (- T) * ( DPDD/PHI + DKSDD/KS + 1/DT )
        DF1DD = F1 * ( DECUDD/ECUNIF - 3*DPDD/PHI )
        DF2DD = (- F2) * DF1DD
        DADD = (- A) * DF2DD / (F2-1)
        DF3DD = (2*T + 4*A*T**3) * DTDD + DADD * T**4
        DF4DD = F4 * ( DF3DD/F3 - (DADD*F3+A*DF3DD)/(1+A*F3) )
        DHDD = 3 * H * DPDD / PHI
        DHDD = DHDD + GAMMA * PHI**3 * DF4DD / (1+F4)
        DFCDD(IS) = VCUNIF(IS) + H + DT * DHDD

        DO 30 IX = 1,3
          DTDGD = (T / GDMT) * GDT(IX) / GDMT
          DF3DGD = DTDGD * ( 2 * T + 4 * A * T**3 )
          DF4DGD = F4 * DF3DGD * ( 1/F3 - A/(1+A*F3) )
          DHDGD = GAMMA * PHI**3 * DF4DGD / (1+F4)
          DFCDGD(IX,IS) = DT * DHDGD
   30   CONTINUE
   40 CONTINUE

C Find exchange energy and potential
      FX = 0
      DO 60 IS = 1,2
        DS(IS)   = MAX( DENMIN, 2 * D(IS) )
        GDMS = MAX( GDMIN, 2 * GDM(IS) )
        KFS = (3 * PI**2 * DS(IS))**THD
        S = GDMS / (2 * KFS * DS(IS))
        F1 = 1 + MU * S**2 / KAPPA
        F = 1 + KAPPA - KAPPA / F1
c
c       Note nspin=1 in call to exchng...
c
        CALL EXCHNG( IREL, 1, DS(IS), EXUNIF, VXUNIF(IS) )
        FX = FX + DS(IS) * EXUNIF * F

        DKFDD = THD * KFS / DS(IS)
        DSDD = S * ( -(DKFDD/KFS) - 1/DS(IS) )
        DF1DD = 2 * (F1-1) * DSDD / S
        DFDD = KAPPA * DF1DD / F1**2
        DFXDD(IS) = VXUNIF(IS) * F + DS(IS) * EXUNIF * DFDD

        DO 50 IX = 1,3
          GDS = 2 * GD(IX,IS)
          DSDGD = (S / GDMS) * GDS / GDMS
          DF1DGD = 2 * MU * S * DSDGD / KAPPA
          DFDGD = KAPPA * DF1DGD / F1**2
          DFXDGD(IX,IS) = DS(IS) * EXUNIF * DFDGD
   50   CONTINUE
   60 CONTINUE
      FX = HALF * FX / DT

C Set output arguments
      EX = FX
      EC = FC
      DO 90 IS = 1,NSPIN
        DEXDD(IS) = DFXDD(IS)
        DECDD(IS) = DFCDD(IS)
        DO 80 IX = 1,3
          DEXDGD(IX,IS) = DFXDGD(IX,IS)
          DECDGD(IX,IS) = DFCDGD(IX,IS)
   80   CONTINUE
   90 CONTINUE

      END
c
          subroutine cal_date(bdate)
c
c    gets the date (day-MONTH-year)
c
          character*10 bdate
          character*3 month(12)
          character*1 dash
          
          integer values(8)
c
      data month/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP',
     &           'OCT','NOV','DEC'/
      data dash/'-'/
c
c       call date_and_time(values=values)

c          write(bdate,102) values(3),dash,month(values(2)),dash,
c     &                     mod(values(1),100)
 102  format(i2,a1,a3,a1,i2.2,' ')

          end
c
      double precision function second() 
c
      real tt
c      call cpu_time(tt)
c      second = dble(tt)
      end
c

