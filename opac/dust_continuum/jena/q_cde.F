      subroutine q_cde(e1,e2,lam,rad,qex,qsc,qabs,g)
      implicit doubleprecision (a-h,o-z)
      real*8 lsol,msol
C            --  Fundamental                constants  --
      parameter (h      =    6.62619650D-27, ! Planck constant
     +           c      = 2.99792456211D+10, ! Velocity of light
     +           ak     =       1.38062D-16, ! Boltzmann constant
     +           rgas   =        8.3143D+7 , ! Gas constant
C            --  Derivatives of fundamental constants --
     +           hdk    = 4.79942639502D-11, ! h/ak
     +           thdcc  = 1.47452759903D-47, ! 2*h/c^2
     +           hcdk   =     1.43883182736D0, ! h*c/ak
C            --  Radiation                  constants  --
     +           stebo  =     5.6696196D-05, ! Stefan-Boltzmann const
C            --  Atomic                     constants  --
     +           amu    =      1.660531D-24, ! atomic mass unit (C12 scale)
C                       =              1/NA, ! NA = Avogadro's number
C            --  General Astronomical       constants  --
     +           year   =         31557600., ! Year of 325.25 days [sec]  
     +           parsec =       3.08561D+18, ! Parsec [cm]
     +           rsol   =       6.95987D+10, ! Solar radius [cm]
     +           lsol   =        3.8268D+33, ! Solar luminosity [ergs/sec]
     +           msol   =        1.9892D+33, ! Solar mass [gr]
C            --  Number                     constants  --
     +           zero   =           0.0D+00, ! 0
     +           half   =           0.5D+00, ! 1/2
     +           one    =           1.0D+00, ! 1
     +           two    =           2.0D+00, ! 2
     +           expmax =         100.0D+00, ! maximum exponent
C            --  Mathematical               constants  --
     +           pi     = 3.1415926535898D0  , ! pi
     +           ln10   = 2.30258509299405D0 ) ! ln(10)
ccccccc      include './include/consts.dust'
      real*8    e1,e2,lam ,rad
      real*8    qex,qsc,qabs,g
      real*8    fact
      complex*8 refrel,epsil
      real*8    x    

c     size parameter
      X = 2.0d0*PI*rad/lam

c     The formula for qabs is valid in the rayleigh limit. The formula is
c     given in Bohren en Huffman p.355-366. The use here is restricted to
c     dust in a vacuum. 

      REFREL = cmplx(e1,e2)
      EPSIL  = REFREL*REFREL
      FACT   = 4.0d0*X/3.0d0
      FACT   = FACT*PI*rad*rad    !calculate the absorbtion coefficient
      qabs   = 2.0d0*IMAG(EPSIL/(EPSIL-1.0d0)*LOG(EPSIL))

      qabs   = FACT*qabs
      qsc    = 0.0d0
      qex    = qabs+qsc
      g      = 1.0d0

      return
      end
c
c



