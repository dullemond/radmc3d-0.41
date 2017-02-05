;==========================================================================
;
;              PARAMETERS FOR THE 2-D AGN TORUS MODELS
;
;==========================================================================

show = 1

;--------------------------------------------------------------------------
;                         The central star
;--------------------------------------------------------------------------
rstar  = 100.d0 * RS    ; Radius of the AGN
mstar  = 1d8 * MS       ; Mass of the AGN
lagn   = 1d11 * LS      ; Luminosity of the AGN
nphot  = 1000000L       ; Nr of photons

;--------------------------------------------------------------------------
;                            Scattering
;--------------------------------------------------------------------------
scat   = 0              ; Include scattering (0=no scattering)?

;--------------------------------------------------------------------------
;                         The spatial grid 
;--------------------------------------------------------------------------
nr     = 178            ; number of points in radius
nt     = 60             ; Number of points in theta for main theta grid
ntex   = 10             ; Extra points in theta 
np     = 360            ; Nr of phi points
rin    = 0.3*pc         ; Inner radius of radial grid (--> tin is computed)
tin    = 0.d0           ; Inner temperature (--> rin is computed)
rout   = 30.0*pc        ; Outer radius of radial grid
hrgrid = 1.2            ; Main grid span in Theta (measured negatively from equator)
hrgmax = 0.9*!pi/2.d0   ; Total grid span in Theta (measured negatively from equator)
nomirror = 1            ; Switch off the equatorial mirroring

;--------------------------------------------------------------------------
;              Grid refinement and smoothing of solution
;--------------------------------------------------------------------------
nlevr  = 3		; higher resolution at level #
nspanr = 2		; number of steps to resolve
nstepr = 2		; number of sub-steps to create
drsm   = 3d-1           ; smoothing radius
rdxsm  = 1d0            ; smoothing level

;--------------------------------------------------------------------------
;                    The Sigma(R) setup parameters
;--------------------------------------------------------------------------
r0     = 10 * pc        ; Pivot of power law definitions:
sig0   = 0.0            ; Sigma at r0 (either sig0 or mdisk)
mdisk  = 4.D6 * MS      ; Mass of the disk (either sig0 or mdisk)
plsig1 =  1.0d0		; Powerlaw for Sigma(R) for R<R0
plsig2 = -12.0d0	; Powerlaw for Sigma(R) for R>R0

;--------------------------------------------------------------------------
;                    Parameters for the blobbiness
;--------------------------------------------------------------------------
nblob  = 200            ; Number of blobs
sblob  = 0.05           ; Dimensionless size of blobs
iseed  = 730911         ; Seed for random number generator
istyle = 10

;--------------------------------------------------------------------------
;           Settings for H(R) (initial value, minimal value)
;--------------------------------------------------------------------------
hrdisk = 0.50           ; Vertical pressure scale height in units of radius
plh    = 0.00           ; Powerlaw for H(R) in setup
cdens  = 1              ; Constant density in Theta direction?

