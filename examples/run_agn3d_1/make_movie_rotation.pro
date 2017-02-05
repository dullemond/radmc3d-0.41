@readradmc
@natconst
;
; Make movie trajectory
;
nframes = 360
phi0    = 0.d0
phi1    = 360.d0
phi     = phi0 + (phi1-phi0)*dindgen(nframes)/(nframes-1.d0)
incl    = 75.
;
; Size of image
;
sizecm  = 14*pc
;
; Write this trajectory to the movie.inp file
;
openw,1,'movie.inp'
printf,1,1
printf,1,nframes
for iframe=0,nframes-1 do begin
   printf,1,0.d0,0.d0,0.d0,sizecm,sizecm,0.d0,incl,phi[iframe]
endfor
close,1
;
; Now call RADMC-3D
;
;spawn,'radmc3d movie npix 200 lambda 10 nofluxcons'
;
; Read the images
;
lgrange = [-13,-8]
xsize=500
ysize=500
window,0,xsize=xsize,ysize=ysize
device,decomposed=0
loadct,3
mpegID = MPEG_OPEN([xsize,ysize],bitrate=10485720.0)

for iframe=1,nframes do begin
   nr=strcompress(string(iframe),/remove_all)
   if iframe lt 10 then nr='0'+nr
   if iframe lt 100 then nr='0'+nr
   if iframe lt 1000 then nr='0'+nr
   a=readimage(ext=nr)
   ;;
   ;; Plot on screen
   ;;
   plotimage,a,/log,lgrange=lgrange,/pc
   ;;
   ;; Read from screen and send to MPEG
   ;;
   image=tvrd(/order,true=1) 
   MPEG_PUT, mpegID,FRAME=iframe, IMAGE=image
   ;; wait,1
endfor
MPEG_SAVE, mpegID, FILENAME='animation_rotation.mpg'
MPEG_CLOSE, mpegID
;
end
