PRO plot_ps_osc_sp

!p.charsize=2
!p.thick=3
loadct,12
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


!p.multi=0
window,0

plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,0.2],/yl,xtit='wavenumber [Mpc^-1]',ytit='PSxk2 [(Mpc)^-5]',/nodata,yra=[10,85],xma=[9,1],yma=[3,1]

for iz=0,2 do begin
   t  = dir + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   xt=ps.(0)
   pt=ps.(1)
   oplot,xt,pt*xt*xt,col=123

   for i=0,9 do begin
      
      prod = '_nCR'+strtrim(i,2)
      
      for igd=0,4 do begin
         t  = dir + 'PS'+prod+nx[iz]+'_z'+namez[iz]+'_G'+strtrim(igd,2)+'_wngal.txt' 
                                ;  print,t
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
      oplot,xs,p2fit*xs*xs,col=lcol[1],psym=1
      
      delta_k = ((p.(0))[15]-(p.(0))[10]) / 5. ; pour moyenner eventuelle erreur d'arrondi
      Vsurvey = nx_cell[iz]*nx_cell[iz]*nz_cell[iz]*cell*cell*cell * 5.
      CteSigmaPk = 2*!pi/sqrt(delta_k * Vsurvey)
      sig_osc = CteSigmaPk /p.(0)* (p2fit + SN)
      
   endfor
   
   pmean = pmean / double(i)
   if (iz eq 0 ) then begin
      pm = dblarr(n_elements(pmean),3)
      pm[*,iz] = pmean
   endif else pm[*,iz] = pmean

   oplot,xs,pmean*xs*xs,th=1
   sigmean = CteSigmaPk /p.(0)* (pmean + SN)
   if (iz eq 0 ) then begin
      sm = dblarr(n_elements(pmean),3)
      sm[*,iz] = sigmean
   endif else sm[*,iz] = sigmean
  errplot,xs,(pmean-sigmean)*xs*xs,(pmean+sigmean)*xs*xs,th=1
   
endfor
mytext = ['1 simu, 5 grids','<10 simu>','spectre theorique']
legend,mytext,line=[-1,0,0],psym=[1,0,0],col=[lcol[1],0,123],box=1,/fill,/right,/bottom,charsize=1.5

write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/ps_osc.jpg' ,tvrd(true=3),true=3
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
window,2,xs=600,ys=900

!p.multi=[0,1,3]
for iz=0,2 do begin
   plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,0.2],xtit='wavenumber [Mpc^-1]',ytit='PS / < PS > ',/nodata,yra=[0.88+0.03*iz,1.12-0.03*iz],xma=[9,2],yma=[3,1],charsi=3

   n1=0
   n2=0
   n3=0
   ok = where(xs ge 0.02 and xs le 0.2,nok)
   n=10*nok
   for i=0,9 do begin
      
      prod = '_nCR'+strtrim(i,2)
      
      for igd=0,4 do begin
         t  = dir + 'PS'+prod+nx[iz]+'_z'+namez[iz]+'_G'+strtrim(igd,2)+'_wngal.txt' 
       ;  print,t
         p = read_ascii(t, template =  TEMP_POW_SPEC_TXT2) 
         
         if (igd eq 0) then p_osc = p.(1) else p_osc = p_osc + p.(1)
      endfor
      SN = mean((p.(3))[where(p.(0) gt 0.2 and p.(0) lt 0.4)])
      
      p_osc = p_osc / 5.      
      p2fit = p_osc - SN 
      
      oplot,xs,p2fit/pm[*,iz],col=lcol[1],psym=1,symsi=2
      sig_osc = CteSigmaPk /p.(0)* (p2fit + SN)
 
      ok1 = where( abs(p2fit[ok]-pm[ok,iz]) ge sm[ok,iz]*1.,nok1) 
      ok2 = where( abs(p2fit[ok]-pm[ok,iz]) ge sm[ok,iz]*2.,nok2) 
      ok3 = where( abs(p2fit[ok]-pm[ok,iz]) ge sm[ok,iz]*3.,nok3) 
      n1 = n1 + nok1
      n2 = n2 + nok2
      n3 = n3 + nok3
                                ;  print,iz,i,nok1,nok2,nok3
   endfor
   errplot,xs,(pm[*,iz]-sm[*,iz])/pm[*,iz],(pm[*,iz]+sm[*,iz])/pm[*,iz],th=2
   errplot,xs,(pm[*,iz]-sm[*,iz]-sm[*,iz])/pm[*,iz],(pm[*,iz]+sm[*,iz]+sm[*,iz])/pm[*,iz],th=2,li=2
   errplot,xs,(pm[*,iz]-sm[*,iz]-sm[*,iz]-sm[*,iz])/pm[*,iz],(pm[*,iz]+sm[*,iz]+sm[*,iz]+sm[*,iz])/pm[*,iz],th=2,li=1
   print,'stat ',iz,100.-n1/double(n)*100.,100.-n2/double(n)*100.,100.-n3/double(n)*100.

   legend,'z='+namez[iz],box=1,/fill,/right,/top,charsize=1.5

endfor

write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/ps_sigma.jpg' ,tvrd(true=3),true=3
end
