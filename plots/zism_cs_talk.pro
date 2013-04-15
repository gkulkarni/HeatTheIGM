
; File: zism_cs.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; This file plots cross-section weighted Z_ISM, in contrast to
; zism.pro.

set_plot, 'ps'
device, filename='zism_cs_talk.ps', xsize=7.0, ysize=7.0, $
        /inches, color=1, /HELVETICA, yoffset=1.0
!P.font = 0 

TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
!P.charsize = 1.5 
!P.thick = 3 
!P.charthick = 1 

restore, 'halo_template.sav'
zism_data = read_ascii('set49/avghalos_zism.out', template=stars_template)
n_data = read_ascii('set34/nofmh.out', template=stars_template) 
m_data = read_ascii('set49/avghalos.out', template=stars_template) 

z = zism_data.field001[0,*]
avgz = dblarr(100)
avgz_cs = dblarr(100) 

for i = 0, 99 do begin 

   zism = double(zism_data.field001[1:*, i])
   n = double(n_data.field001[1:*, i])
   m = double(m_data.field001[1:*, i])
   num = n_elements(n) 

   dr = total(n)
   nr = total(n*zism)
   avgz[i] = nr/dr 

   cs = dblarr(num) 
   for j = 0, num-1 do begin 
      cs[j] = double(sigma(m[j],z[i]))
   endfor

   dr = total(n*cs) 
   nr = total(n*cs*zism) 
   avgz_cs[i] = nr/dr

endfor 

plot, z, avgz_cs, xrange=[1,20], yrange=[1.0e-4,1.0e0], $
      /ylog, /xlog, ytickformat='exponent', xtitle='redshift', $
      ytitle='Z / Zsun', xthick=3, ythick=3 

openr, lun, '../data/prochaska03.dat', /get_lun
obsdata_zism = fltarr(2, 109)
readf, lun, obsdata_zism
z_obs = obsdata_zism[0,*]
zism_obs = obsdata_zism[1,*]
oplot, z_obs, zism_obs, psym=7, color=2, thick=4

device, /close_file
set_plot, 'X'

END

