PRO plot_sp


!p.charsize=2
!p.thick=3
loadct,39
restore,'temp_ascii_new.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO2_SP/"

namez = ['0.5','0.9','1.5']
namezmean=['0.51', '0.93', '1.57']
nx=['_120','_225','_320']
nx_cell= [120.,225., 320.]
nz_cell= [125.,125.,175.]
cell=8.
!p.multi=0

lcol  = [0,95,210,35,135,  120]
;plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,0.2],xtit='wavenumber [Mpc^-1]',ytit='<PS> / P_theo',/nodata,yra=[0.5,1.2],xma=[9,1],yma=[3,1]
plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,0.2],xtit='wavenumber [Mpc^-1]',ytit='<PS> / P_theo',/nodata,yra=[0.,2],xma=[9,1],yma=[3,1]
for iz=0,2 do begin
   nsim=0
   for i=0,9 do begin ; probleme avec le #9 ...
      
      prod = '_nCR'+strtrim(i,2)
     ; prod = '_norm'+strtrim(i,2)

      for igd=0,4 do begin
         t  = dir + 'PS'+prod+nx[iz]+'_z'+namez[iz]+'_G'+strtrim(igd,2)+'_wngal.txt' 
         p = read_ascii(t, template =  TEMP_POW_SPEC_5col) 
         
         if (igd eq 0 and i eq 0) then begin
            xs = p.(0)
            p_mean =dblarr(n_elements(xs))
            p_sn   =p_mean
         endif
         
         p_mean = p_mean + p.(1)
         p_sn  = p_sn  + p.(3)
         
         nsim ++
      endfor
   endfor
   print,'N PS = ',nsim
   SN = mean(p_sn[where(xs gt 0.2 and xs lt 0.4)])

   if (iz eq 0) then pm = dblarr(n_elements(xs),3)
   if (iz eq 0) then SNm = dblarr(3)
   p_mean = p_mean - SN
   pm[*,iz] = p_mean / nsim
   SNm[iz] = SN / nsim
   t  = dir + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   xt=ps.(0)
   pt=ps.(1)
   ptinter = interpol(pt,xt,xs)
   oplot,xs,pm[*,iz] / ptinter,col = 50 + iz*50
endfor
legend,'z = '+namez,line=0,col=50+indgen(3)*50,box=1,/fill,/right,/bottom,charsize=1.5
read,xx
;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/ps_sur_pth.jpg' ,tvrd(true=3),true=3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.02,0.2],xtit='wavenumber [Mpc^-1]',ytit='dispersion(10 Pk) / sig_Pk',/nodata,yra=[0,2],xma=[9,1],yma=[3,1]
for iz=0,2 do begin
   for i=0,9 do begin ; probleme avec le #9 ...
      
      prod = '_nCR'+strtrim(i,2)
      
      for igd=0,4 do begin
         t  = dir + 'PS'+prod+nx[iz]+'_z'+namez[iz]+'_G'+strtrim(igd,2)+'_wngal.txt' 
         p = read_ascii(t, template =  TEMP_POW_SPEC_5col) 

         if (i eq 0 and igd eq 0) then p_var = dblarr(n_elements(xs),10)
         p_var[*,i] = p_var[*,i] + p.(1)
      endfor

      p_var[*,i] = (p_var[*,i] / 5. - SNm[iz]) - pm[*,iz]
      pkvar = dblarr(n_elements(xs))
   endfor
   for j = 0, n_elements(xs)-1 do pkvar[j] = stddev(p_var[j,*])

   delta_k = ( (p.(0))[15]-(p.(0))[10] ) / 5. ; pour moyenner eventuelle erreur d'arrondi
   Vsurvey = nx_cell[iz]*nx_cell[iz]*nz_cell[iz]*cell*cell*cell * 5.
   CteSigmaPk = 2*!pi/sqrt(delta_k * Vsurvey)
   sig_osc = CteSigmaPk /p.(0)* (pm[*,iz] + SNm[iz])

   oplot,xs,pkvar / sig_osc,col = 50 + iz*50
   print,xs[10],pm[10,iz], pkvar[10], sig_osc[10],SNm[iz]
   ;read,xx
endfor
legend,'z = '+namez,line=0,col=50+indgen(3)*50,box=1,/fill,/right,/bottom,charsize=1.5
stop
;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/ps_sur_pth_sigma.jpg' ,tvrd(true=3),true=3

END



