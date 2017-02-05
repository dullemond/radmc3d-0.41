@readradmc.pro

set_plot,'x'
device,decomposed=0
loadct,5

common colors,r_orig,g_orig,b_orig,r_curr,g_curr,b_curr

spawn,'nice radmc3d child',unit=iounit

for i = 90, 90, 10 do begin

print, 'generating figure ', i

incl = 90.
phi  = 1.0 * i
npix = 800
ilambda = 145
;zoomau = [-200000.0d0,200000.0d0,-200000.0d0,200000.0d0]
;zoomau = [-100000.0d0,100000.0d0,-100000.0d0,100000.0d0]
zoomau = [-70000.0d0,70000.0d0,-70000.0d0,70000.0d0]

makeimage,iounit=iounit,incl=incl,phi=phi,npix=npix,posang=posang,$
             ifreq=ilambda,nostar=nostar,zoomau=zoomau,$
             nofluxcons=nofluxcons,plottau=plottau
image=readimage(iounit=iounit)
;plotimage,image,/log,lgrange=[-20,-15],/pc

;lgrange = [-20.0,-15.0]
lgrange=[alog10(min(image.image)),alog10(max(image.image))]

imgidx = floor( ( (alog10(image.image>1d-30)-lgrange[0])/(lgrange[1]-lgrange[0]) < 1. > 0. ) * 255 )

zf = 1

image_rgb=dblarr(3,image.nx,image.ny)
image_rgb[0,*,*]=r_orig[imgidx]
image_rgb[1,*,*]=g_orig[imgidx]
image_rgb[2,*,*]=b_orig[imgidx]
q  = rebin(image_rgb,3,zf*image.nx,zf*image.ny)
 
fnum  = string(i,'(i4.4)')
fname = 'movie' + fnum + '.jpg'

write_jpeg,fname,q,true=1

print, 'generated figure ', i

endfor

printf,iounit,'quit'

end
