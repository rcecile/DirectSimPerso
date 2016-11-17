PRO save_ps_mean_err

namez = ['0.9','1.3','1.8','1.8']
nsuff = ['', '', ' thin', ' thick']
nx= ['_640', '_900', '_1024', '_500']
nameerr = ['', '_err0.03', '_errPpodds']


lcol  = [50,90,250,150]
lpsym=[0,-5,-4]
llin = [0,3,2]

nz = n_elements(namez)
nerr = n_elements(nameerr)
nsim=10

restore,'temp_ascii.sav'   
dir="/sps/lsst/data/rcecile/TJP_BAO_PS/"
dirn="/sps/lsst/data/rcecile/TJP_noBAO_PS/"
myformat = '(F9.7," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",G14.6," ",I6," ",I6," ",G14.6)'
loadct,39
!p.thick=3
!p.charsize=2

ftheo  = dir + 'simu_ps'+nx[0]+'_z'+namez[0]+'_ntpk.txt'
ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
xt=ptheo.(0)
nxt = n_elements(xt)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; CAS SPECTRO : ON MIXE PAR CELLE SIZE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ir=0
suff = nameerr[ir]

plot,[0.01,0.1],[1,1],/xs,/ys,xra=[0.005,0.15],/nodata,yra=[1000,.5e5],/yl,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum'

;;; CELL 8 Mpc

fobs  = dir + 'PS_G2'+nx[0]+'_z'+namez[0]+ nameerr[0]+'_wngal.txt' 
pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_TXT)     
xs = pobs.(0)
nxs = n_elements(xs)
ok = intarr(nxs)
for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))

rapnobao=dblarr(nxs)
for iz = 0,2 do begin
   rap=dblarr(nxs)

   ftheo  = dir + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
   ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
   for j=0,nsim-1 do begin
      fnobao  = dirn + 'PS_G2_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
      p = read_ascii(fnobao, template =  TEMP_POW_SPEC_TXT)     
      rap = rap + (p.(1)-p.(6))/(ptheo.(2))[ok]
   endfor
   rapnobao = rapnobao + rap/nsim
endfor
rapnobao = rapnobao / 3.

for iz = 0,2 do begin
   ; spectre th??orique no osc, no error
   
   ftheo  = dir + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
   ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     

   fobs  = dir + 'PS_G2'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
   pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_TXT)     
                                   ; < spectre simul?? no osc, with error >
   
   p2fit = (pobs.(1)-pobs.(6)) / rapnobao
   
   oplot,xs,p2fit,col=lcol[iz]
   
   fname = dir + 'PS_G2_2fitErrk'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
   print,fname
   openw,lun1, fname, /get_lun
   
   newp1 = p2fit + pobs.(6)
   
   newp7 = (newp1-pobs.(6))/(pobs.(1)-pobs.(6))*pobs.(7) ; erreur "undampee", comme le signal
   print,'  ' ,suff,(pobs.(0))[40],(pobs.(7))[40],newp7[40],(newp1[40]-(pobs.(6))[40])/((pobs.(1)[40]-pobs.(6))[40])
   
   for i=0,n_elements(p2fit)-1 do  printf, lun1, (pobs.(0))[i], $
                                           newp1[i], $ ;spectre
                                           (pobs.(2))[i], (pobs.(3))[i], (pobs.(4))[i], (pobs.(5))[i], (pobs.(6))[i], $
                                           newp7[i], $ ; erreur
                                           (pobs.(8))[i], (pobs.(9))[i], (pobs.(10))[i],format=myformat
   close, lun1
   free_lun, lun1
   oplot,xt,ptheo.(2),li=2
endfor
legend,'z = '+ namez+nsuff,col=lcol,box=1,line=0,/fill,/left,/bottom,charsize=1.5

;;; CELL 16 Mpc
iz=3

fobs  = dir + 'PS_G2'+nx[iz]+'_z'+namez[iz]+ nameerr[0]+'_wngal.txt' 
pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_TXT)     
xs = pobs.(0)
nxs = n_elements(xs)
ok = intarr(nxs)
for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))
rapnobao=dblarr(nxs)

                                ; spectre th??orique no osc, no error
ftheo  = dir + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
for j=0,nsim-1 do begin
   fnobao  = dirn + 'PS_G2_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
   p = read_ascii(fnobao, template =  TEMP_POW_SPEC_TXT)     
   rapnobao = rapnobao + (p.(1)-p.(6))/(ptheo.(2))[ok]
endfor
rapnobao = rapnobao / double(nsim)
                                ; spectre observ??
fobs  = dir + 'PS_G2'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_TXT)     
                                ; < spectre simul?? no osc, with error >

p2fit = (pobs.(1)-pobs.(6)) / rapnobao

oplot,xs,p2fit,col=lcol[iz]

fname = dir + 'PS_G2_2fitErrk'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
print,fname
openw,lun1, fname, /get_lun

newp1 = p2fit + pobs.(6)

newp7 = (newp1-pobs.(6))/(pobs.(1)-pobs.(6))*pobs.(7) ; erreur "undampee", comme le signal
print,'  ' ,suff,(pobs.(0))[40],(pobs.(7))[40],newp7[40],(newp1[40]-(pobs.(6))[40])/((pobs.(1)[40]-pobs.(6))[40])

for i=0,n_elements(p2fit)-1 do  printf, lun1, (pobs.(0))[i], $
                                        newp1[i], $ ;spectre
                                        (pobs.(2))[i], (pobs.(3))[i], (pobs.(4))[i], (pobs.(5))[i], (pobs.(6))[i], $
                                        newp7[i], $ ; erreur
                                        (pobs.(8))[i], (pobs.(9))[i], (pobs.(10))[i],format=myformat
close, lun1
free_lun, lun1
oplot,xt,ptheo.(2),li=2
legend,'z = '+ namez+nsuff,col=lcol,box=1,line=0,/fill,/left,/bottom,charsize=1.5
legend,lerr,li=0,box=1,/fill,/center,/bottom,charsize=1.5
read,xx

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; CAS GAUSS : ON MIXE TOUT cell8 et on exporte pour le dernier
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ir=1
suff = nameerr[ir]


;;; CELL 8 Mpc

plot,[0.01,0.1],[1,1],/xs,/ys,xra=[0.005,0.15],/nodata,yra=[1000,.5e5],/yl,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum'
fobs  = dir + 'PS_G2'+nx[0]+'_z'+namez[0]+ nameerr[0]+'_wngal.txt' 
pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_TXT)     
xs = pobs.(0)
nxs = n_elements(xs)
ok = intarr(nxs)
for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))

rapnobao=dblarr(nxs)
for iz = 0,2 do begin ; on ne melange pas avec thick car pas meme echantillonage en x : galere et comme c est pas terrible ...
   rap=dblarr(nxs)
   
   ftheo  = dir + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
   ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
   for j=0,nsim-1 do begin
      fnobao  = dirn + 'PS_G2_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
      p = read_ascii(fnobao, template =  TEMP_POW_SPEC_TXT)  
      rap = rap + (p.(1)-p.(6))/(ptheo.(2))[ok]
   endfor
   rapnobao = rapnobao + rap/nsim
endfor
rapnobao = rapnobao / 3.

for iz = 0,2 do begin
   
   ; spectre th??orique no osc, no error
   ftheo  = dir + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
   ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
   
                                ; spectre observ??
   fobs  = dir + 'PS_G2'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
   pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_TXT)     

   p2fit = (pobs.(1)-pobs.(6)) / rapnobao

   oplot,xs,p2fit,col=lcol[iz]
   
   fname = dir + 'PS_G2_2fitErrk'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
   print,fname
   openw,lun1, fname, /get_lun
   
   newp1 = p2fit + pobs.(6)
   
   newp7 = (newp1-pobs.(6))/(pobs.(1)-pobs.(6))*pobs.(7) ; erreur "undampee", comme le signal
   print,'  ' ,suff,(pobs.(0))[40],(pobs.(7))[40],newp7[40],(newp1[40]-(pobs.(6))[40])/((pobs.(1)[40]-pobs.(6))[40])

   for i=0,n_elements(p2fit)-1 do  printf, lun1, (pobs.(0))[i], $
                                           newp1[i], $ ;spectre
                                           (pobs.(2))[i], (pobs.(3))[i], (pobs.(4))[i], (pobs.(5))[i], (pobs.(6))[i], $
                                           newp7[i], $ ; erreur
                                           (pobs.(8))[i], (pobs.(9))[i], (pobs.(10))[i],format=myformat
   close, lun1
   free_lun, lun1
oplot,xt,ptheo.(2),li=2
endfor
;;; CELL 16 Mpc
iz=3

xs_8 = xs

fobs  = dir + 'PS_G2'+nx[iz]+'_z'+namez[iz]+ nameerr[0]+'_wngal.txt' 
pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_TXT)     
xs = pobs.(0)
nxs = n_elements(xs)
ok = intarr(nxs)
for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))
rapnobao16 = interpol(rapnobao,xs_8,xs)

                                ; spectre th??orique no osc, no error
ftheo  = dir + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     

                                ; spectre observ??
fobs  = dir + 'PS_G2'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_TXT)     

p2fit = (pobs.(1)-pobs.(6)) / rapnobao16

oplot,xs,p2fit,col=lcol[iz]

fname = dir + 'PS_G2_2fitErrk'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
print,fname
openw,lun1, fname, /get_lun

newp1 = p2fit + pobs.(6)

newp7 = (newp1-pobs.(6))/(pobs.(1)-pobs.(6))*pobs.(7) ; erreur "undampee", comme le signal
print,'  ' ,suff,(pobs.(0))[40],(pobs.(7))[40],newp7[40],(newp1[40]-(pobs.(6))[40])/((pobs.(1)[40]-pobs.(6))[40])

for i=0,n_elements(p2fit)-1 do  printf, lun1, (pobs.(0))[i], $
                                        newp1[i], $ ;spectre
                                        (pobs.(2))[i], (pobs.(3))[i], (pobs.(4))[i], (pobs.(5))[i], (pobs.(6))[i], $
                                        newp7[i], $ ; erreur
                                        (pobs.(8))[i], (pobs.(9))[i], (pobs.(10))[i],format=myformat
close, lun1
free_lun, lun1
oplot,xt,ptheo.(2),li=2
legend,'z = '+ namez+nsuff,col=lcol,box=1,line=0,/fill,/left,/bottom,charsize=1.5
read,xx
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; CAS PODDS : ON MIXE PAR GRILLE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ir = 2
suff = nameerr[ir]

plot,[0.01,0.1],[1,1],/xs,/ys,xra=[0.005,0.15],/nodata,yra=[1000,.5e5],/yl,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum'

for iz = 0,nz-1 do begin
   
   ; spectre th??orique no osc, no error
   ftheo  = dir + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
   ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
                                ; spectre observ??
   fobs  = dir + 'PS_G2'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
   pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_TXT)     
   xs = pobs.(0)
   ok = intarr(nxs)
   for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))
   
                                ; < spectre simul?? no osc, with error >
   rapnobao=dblarr(nxs)
   for j=0,nsim-1 do begin
      fnobao  = dirn + 'PS_G2_'+strtrim(j,2)+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
      p = read_ascii(fnobao, template =  TEMP_POW_SPEC_TXT)     
      rapnobao = rapnobao + ( p.(1)-p.(6) )/(ptheo.(2)[ok])
   endfor
   rapnobao = rapnobao / double(nsim)
   p2fit = (pobs.(1)-pobs.(6)) / (rapnobao) 
   
   oplot,xs,p2fit,col=lcol[iz]
   
   fname = dir + 'PS_G2_2fitErrk'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
   print,fname
   openw,lun1, fname, /get_lun
   
   newp1 = p2fit + pobs.(6)
   
   newp7 = (newp1-pobs.(6))/(pobs.(1)-pobs.(6))*pobs.(7) ; erreur "undampee", comme le signal
   print,'  ' ,suff,(pobs.(0))[40],(pobs.(7))[40],newp7[40],(newp1[40]-(pobs.(6))[40])/((pobs.(1)[40]-pobs.(6))[40])
   
   for i=0,n_elements(p2fit)-1 do  printf, lun1, (pobs.(0))[i], $
                                              newp1[i], $ ;spectre
                                              (pobs.(2))[i], (pobs.(3))[i], (pobs.(4))[i], (pobs.(5))[i], (pobs.(6))[i], $
                                              newp7[i], $ ; erreur
                                              (pobs.(8))[i], (pobs.(9))[i], (pobs.(10))[i],format=myformat
   close, lun1
   free_lun, lun1
  ;    print,'  ',fname
   oplot,xt,ptheo.(2),li=2
endfor
legend,'z = '+ namez+nsuff,col=lcol,box=1,line=0,/fill,/left,/bottom,charsize=1.5

END
