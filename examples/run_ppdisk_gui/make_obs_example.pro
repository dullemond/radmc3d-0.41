@../../idl/readradmc.pro

spec = readspectrum()

powav  = [0.55, 0.78, 0.9, 1.25, 1.65, 2.25, 3.8, 4.5, 5.8, 8.0, 12.0, 25.0, 60.0, 100., 200., 850., 1300.]
poflux = 10.^interpol(alog10(spec.spectrum), alog10(spec.lambda), alog10(powav))
noise = randomn(seed, n_elements(poflux)) * (0.1 * poflux)
ii    = where(powav gt 8.0)
noise(ii) = noise(ii) * 3.
npoflux = poflux + noise


sowav1  = 1.25 * (2.5/1.25)^(dindgen(200)/199)
soflux1 = 10.^interpol(alog10(spec.spectrum), alog10(spec.lambda), alog10(sowav1))
noise = randomn(seed, n_elements(soflux1)) * (0.2 * soflux1)
snoflux1 = soflux1 + noise

sowav2  = 5.5 * (35/5.5)^(dindgen(380)/279)
soflux2 = 10.^interpol(alog10(spec.spectrum), alog10(spec.lambda), alog10(sowav2))
noise = randomn(seed, n_elements(soflux2)) * (0.1 * soflux2)
snoflux2 = soflux2 + noise


window,0
plot_oo, spec.lambda, spec.spectrum*2.9979d14/spec.lambda, xr=[0.1, 1000]
oplot, powav, npoflux*2.9979d14/powav, psym=4, col=rgb(255,0,0), thick=2
oplot, sowav1, snoflux1*2.9979d14/sowav1, col=rgb(0,255,0)
oplot, sowav2, snoflux2*2.9979d14/sowav2, col=rgb(0,255,0)

openw, 1, 'obs_data_test.dat'
for i=0, 

close, 1

end
