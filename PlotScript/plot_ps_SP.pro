PRO plot_ps_sp

!p.charsize=2
!p.thick=3
loadct,39
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO2_SP/"

namez = ['0.5','0.9','1.5']
namezmean=['0.51', '0.93', '1.57']
nx=['_120','_225','_320']
nx_cell= [120.,225., 320.]
nz_cell= [125.,125.,175.]
cell=8.

nline=[2,2,0,0];,0,0]
nz = n_elements(namez)
nk = n_elements(namek)
ngal_mean=[15.66,4.65,0.84,  15.53,4.66,0.84,  16.31,4.82,0.85,  14.20,4.59,0.73,   12.23,4.29,0.61]
SN = 1./ngal_mean*8.*8.*8.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0

lcol  = [0,95,210,35,135,  120]

xmin = 0.02
xmax = 0.04
window,0
plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,0.2],/yl,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum x k2 [(Mpc)^-5]',/nodata,yra=[4,80],xma=[9,1],yma=[3,1]
for iz=0,2 do begin

   nsim=0
   for i=0,8 do begin ; probleme avec le #9 ...
      
      prod = '_norm'+strtrim(i,2)
      
      for igd=0,4 do begin
         t  = dir + 'PS'+prod+nx[iz]+'_z'+namez[iz]+'_G'+strtrim(igd,2)+'_wngal.txt' 
         p = read_ascii(t, template =  TEMP_POW_SPEC_TXT2) 
         
         if (igd eq 0) then begin
            xs = p.(0)
            p_osc =p.(1)
         endif else p_osc = p_osc + p.(1)
      endfor
      SN = mean((p.(3))[where(p.(0) gt 0.2 and p.(0) lt 0.4)])
      
      p_osc = p_osc / 5.      
      p2fit = p_osc - SN 
      
      if (i eq 0) then pmean = p2fit else  pmean = pmean + p2fit
      oplot,xs,p2fit*xs*xs,col=50+iz*50
         nsim ++
   endfor
   print,'n PS for pmean ',double(nsim)
   pmean = pmean / double(nsim)
   
   if (iz eq 0 ) then begin
      pm = dblarr(n_elements(pmean),3)
      pm[*,iz] = pmean
   endif else pm[*,iz] = pmean
   t  = dir + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   xt=ps.(0)
   pt=ps.(1)
   oplot,xt,pt*xt*xt
   ptinter = interpol(pt,xt,xs)
   if (iz eq 0 ) then begin
      pth = dblarr(n_elements(pmean),3)
      pth[*,iz] = ptinter
   endif else pth[*,iz] = ptinter

endfor
mytext = ['Normalisation 1, z='+namez,'spectre theorique']
legend,mytext,line=0,col=[50+indgen(3)*50,0],box=1,/fill,/right,/bottom,charsize=1.5

write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/ps_norm.jpg' ,tvrd(true=3),true=3

window,2
plot,xs,pth[*,0]/pm[*,0],/xs,/ys,xra=[0.01,0.1],yra=[1,2],/nodata
for iz=0,2 do begin

   oplot,xs,pth[*,iz]/pm[*,iz],psym=-2,symsi=3,col=50+iz*50
   ok = where(xs ge xmin and xs le xmax)
   print,iz,sqrt(mean(pth[ok,iz]/pm[ok,iz])),sqrt(median(pth[ok,iz]/pm[ok,iz])),mean(pth[ok,iz]/pm[ok,iz]),median(pth[ok,iz]/pm[ok,iz])
endfor
oplot,[xmin,xmin],[0,2],li=2
oplot,[xmax,xmax],[0,2],li=2

;stop
write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/ps_norm_ratio.jpg' ,tvrd(true=3),true=3
end
