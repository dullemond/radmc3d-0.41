@readradmc
@natconst
;
; Make movie trajectory
;
;nframes = 360
nframes = 300
;nframes = 1
x0      = 0.5*pc
y0      = 10*pc
x       = x0 + dblarr(nframes)
y       = -y0 + 2*y0*dindgen(nframes)/(nframes-1.d0)
;y       = 2*y0*dindgen(nframes)/(nframes-1.d0)
z       = dblarr(nframes)
incl    = 75.*2*!dpi/360.
phi     = 30.*2*!dpi/360.
y1      = y
z1      = z
y       = cos(!dpi/2-incl)*y1 + sin(!dpi/2-incl)*z1
z       = -sin(!dpi/2-incl)*y1 + cos(!dpi/2-incl)*z1 
x1      = x
y1      = y
x       = cos(phi)*x1 + sin(phi)*y1
y       = -sin(phi)*x1 + cos(phi)*y1
;
; Size of image
;
sizerad = 0.8d0
;
; Write this trajectory to the movie.inp file
;
openw,1,'movie.inp'
printf,1,-1
printf,1,nframes
for iframe=0,nframes-1 do begin
   printf,1,0.d0,0.d0,0.d0,sizerad,sizerad,0.d0,x[iframe],z[iframe],z[iframe]
endfor
close,1
;
; Now call RADMC-3D
;
spawn,'radmc3d movie npix 200 lambda 10 nofluxcons'
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
MPEG_SAVE, mpegID, FILENAME='animation_journey.mpg'
MPEG_CLOSE, mpegID
;
end
