
; File: moster.pro
;  Cre: 2012
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.2 $)

; Plots the Moster relation in our code.

set_plot, 'ps'
device, filename='moster.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0
;window, xsize=1000, ysize=1000
;Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4
!P.charsize = 2.0
!P.thick = 3.0
!P.charthick = 3 

restore, 'halo_template.sav'
halo_data = read_ascii('set36/halos.out', template=stars_template)
stars_data = read_ascii('set36/halos_stars.out', template=stars_template)
halos = halo_data.field001[1:*,99]
stars = stars_data.field001[1:*,99]
mstar = stars/halos 
mss = smooth(mstar, 30, /edge_truncate) 
;plot, halos, mstar, /xlog, /ylog, xrange=[5.0e-1, 1.0e6],yrange=[1.0e-2,1.0e-1] 
plot, halos, mss, /xlog, /ylog, xrange=[5.0e-1, 1.0e6], yrange=[1.0e-3,1.0e-1], xtitle='!NM!Dhalo!N/10!E10!N M!D' + sun + '!N', ytitle='!NM!D*!N/M!Dhalo!N'
xyouts, 1.0e-1, 0.06, 'z = 0', charsize=2.0

gamma = 0.6
beta = 1.5
n = 0.035
m1 = 10.0^11.5
moster = fltarr(270)
for i=0, 269 do begin moster[i]=2.0*n/((halos[i]*1.0e10/m1)^(-beta)+(halos[i]*1.0e10/m1)^gamma)
oplot, halos, moster, linestyle=2, color=2

;---------------------------

halo_data = read_ascii('set22/halos.out', template=stars_template)
stars_data = read_ascii('set22/halos_stars.out', template=stars_template)
halos = halo_data.field001[1:*,99]
stars = stars_data.field001[1:*,99]
mstar = stars/halos 
mss = smooth(mstar, 30, /edge_truncate) 
; oplot, halos, mss, color=2

gamma = 0.9
beta = 0.9
n = 0.02
m1 = 10.0^12.4
moster = fltarr(270)
for i=0, 269 do begin moster[i]=2.0*n/((halos[i]*1.0e10/m1)^(-beta)+(halos[i]*1.0e10/m1)^gamma)
;oplot, halos, moster, linestyle=2, color=2

;---------------------------

halo_data = read_ascii('set33/halos.out', template=stars_template)
stars_data = read_ascii('set33/halos_stars.out', template=stars_template)
halos = halo_data.field001[1:*,99]
stars = stars_data.field001[1:*,99]
mstar = stars/halos 
mss = smooth(mstar, 30, /edge_truncate) 
;oplot, halos, mss, color=4

gamma = 0.9
beta = 0.8
n = 0.015
m1 = 10.0^12.4
moster = fltarr(270)
for i=0, 269 do begin moster[i]=2.0*n/((halos[i]*1.0e10/m1)^(-beta)+(halos[i]*1.0e10/m1)^gamma)
;oplot, halos, moster, linestyle=2, color=4

device, /close_file
set_plot, 'X'

