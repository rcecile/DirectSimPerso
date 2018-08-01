PRO plot_sn,saveplot,iroot
;;;;;;;;;; ATTENTION pas tout a jour pour grille z=1.3
lcol  = [95,210,35,135,  120,0]
loadct,12
lroot = ["_superSN5","_testSN_120","_testSN_80"]
ltit=["paper grids","grids 120x120x120","grids 80x80x80"]
if (iroot eq 0) then nx=['_120','_225','_300'] else nx=['','','']

root = lroot[iroot]
print,'root = ',root

!p.charsize=1.5

restore,'temp_ascii_new.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO2_PS/"
dirg='/sps/lsst/data/rcecile/Planck_BAO2_grids/'
namez = ['0.5','0.9','1.3']
namezmean = ['0.51','0.93','1.36']
lpsym=[0,4,0]
llin = [2,0,0]

kmin=[0.1,0.07,0.02]
nz = n_elements(namez)
nk = n_elements(namek)
myxt='wavenumber [Mpc^-1]'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

!p.multi=0
ngal_mean=[15.66,4.65,0.84,  15.53,4.66,0.84,  16.31,4.82,0.85,  14.20,4.59,0.73,   12.23,4.29,0.61]
SN = 1./ngal_mean*8.*8.*8.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; COMPUTE SN
if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_shotnoise'+root+'.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif else window,iroot

nameerr = ['', '_err0.03', '_errP','_errPBDT9','_errPBDT8']
namesf = ['SelFunc_speczmean_nofz.txt','SelFunc_gaussmean_nofz.txt','SelFunc_photzmean_nofz.txt','SelFunc_pbdt9mean_nofz.txt','SelFunc_pbdt8mean_nofz.txt']
SNtoUseC = dblarr(n_elements(nameerr),n_elements(namez),80)


SN0 = dblarr(n_elements(namez))
for isuff=0,n_elements(nameerr)-1 do begin
;for isuff=0,2,2 do begin
   for iz=0,n_elements(namez)-1 do begin      
      suff = nameerr[isuff]
      if (isuff eq 0 and iz eq 0) then plot,[0.001,0.4],[1e4,1e4],/xs,/ys,xra=[0.01,.4],xtit=myxt,ytit='PS of random grids [Mpc^-3]',/nodata,$
                                /yl,xma=[9,1],/xl,tit=ltit[iroot],yra=[25,3.e4];,yma=[4,.5]
      for ig=0,4  do begin
         t  = dir + 'PS'+root+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(ig,2)+'_wngal.txt' 
         p = read_ascii(t, template =  TEMP_POW_SPEC_4col) 
         
         oplot,p.(0),p.(2),th=4,li=llin[iz],psym=lpsym[iz],col=lcol[isuff]
         oplot,[0.01,0.5],[SN[isuff*3 +iz],SN[isuff*3 +iz]],th=1 ;,col=lcol[isuff],li=2
         
         mySN = mean((p.(2))[where(p.(0) gt 0.2 and p.(0) lt 0.4)])
         print,iz,isuff,mySN
      endfor
   endfor
endfor

lerr=['spectroZ','Gaussian 0.03','photoZ','photoZ BDT 90%','photoZ BDT 80%']
legend,[lerr],li=0,col=[lcol],th=3,box=1,/fill,/right,/top,charsize=1.25
legend,'z='+namez,li=[2,0,0],psym=[-3,4,-3],box=1,col=lcol[4],/fill,/center,/top,charsize=1.5,th=3

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif
write_jpeg,'/sps/lsst/users/rcecile/BAO_InOut/sn'+root+'.jpg' ,tvrd(true=3),true=3

stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

nx= ['_120','_225', '_320']
nx_cell= [120.,225., 320.]
nz_cell= [125.,125.,175.]
cell=8.

use_cte = 1
window,2
plot,[0.01,0.4],[1e4,1e4],/xs,/ys,xra=[0.01,.3],xtit=myxt,ytit='ratio PS / <PS spectro>',/nodata,$
                                yra=[0.9,1.4],xma=[9,1],yma=[4,.5],/xl
lsuff=['','_err0.03_00','_err0.03_02','_err0.03_03','_err0.03_04']
lcol  = [95,210,20,120,160,60]

for iz=0,2 do begin
   suff = lsuff[0]
   kmax = 0.01-0.01*iz
   print,iz,kmax
   for ig=0,4  do begin
      t  = dir + 'PS'+root+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(ig,2)+'_wngal.txt' 
      p = read_ascii(t, template =  TEMP_POW_SPEC_4col) 
      if (ig eq 0) then psepc = p.(2) else  psepc += p.(2)
      ok = where(p.(0) gt kmax)
   endfor
   psepc /= double(ig)
   for isuff=0,n_elements(lsuff)-1 do begin
      suff = lsuff[isuff]
      for ig=0,4  do begin
         t  = dir + 'PS'+root+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(ig,2)+'_wngal.txt' 
         p = read_ascii(t, template =  TEMP_POW_SPEC_4col) 
         if (ig eq 0) then perr = p.(2) else perr += p.(2)
      endfor
      perr /= ig
      oplot,p.(0)[ok],perr[ok]/psepc[ok],th=4,li=llin[iz],psym=lpsym[iz],col=lcol[isuff]
   endfor
   
endfor
lerr=['spectroZ','Gaussian 0.00','Gaussian 0.02','Gaussian 0.03','Gaussian 0.04']
legend,[lerr],li=0,col=[lcol],th=3,box=1,/fill,/right,/top,charsize=1.25
legend,'z='+namez,li=[2,0,0],psym=[-3,4,-3],box=1,col=lcol[3],/fill,/center,/top,charsize=1.5,th=3
;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/sn_gauss_ratio.jpg' ,tvrd(true=3),true=3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

window,3
plot,[0.01,0.4],[1e4,1e4],/xs,/ys,xra=[0.02,.3],xtit=myxt,ytit='PS',/nodata,$
                                yra=[40,1500],/yl,xma=[9,1],yma=[4,.5],/xl
lsuff=['','_err0.03_00','_err0.03_02','_err0.03_03','_err0.03_04']
lcol  = [95,210,20,120,160,60]

for iz=0,2 do begin
   kmax = 0.04-0.01*iz
   print,iz,kmax
   
   suff = lsuff[0]
   for isuff=0,n_elements(lsuff)-1 do begin
      suff = lsuff[isuff]
      for ig=0,4  do begin
         t  = dir + 'PS'+root+nx[iz]+'_z'+namez[iz]+suff+'_G'+strtrim(ig,2)+'_wngal.txt' 
         p = read_ascii(t, template =  TEMP_POW_SPEC_4col) 
         ok = where(p.(0) gt kmax)
         oplot,p.(0)[ok],p.(2)[ok],th=4,li=llin[iz],psym=lpsym[iz],col=lcol[isuff]
      endfor
   endfor
   
endfor
legend,[lerr],li=0,col=[lcol],th=3,box=1,/fill,/right,/center,charsize=1.25
legend,'z='+namez,li=[2,0,0],psym=[-3,4,-3],box=1,col=lcol[3],/fill,/left,/bottom,charsize=1.5,th=3
;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/sn_gauss_ps.jpg' ,tvrd(true=3),true=3
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

window,3
plot,[0.01,0.4],[1e4,1e4],/xs,/ys,xra=[0.02,.2],xtit=myxt,ytit='PS',/nodata,$
                                yra=[120,700],/yl,xma=[9,1],yma=[4,.5],/xl
lsuff=['_err0.03_00','_err0.03_02','_err0.03_03','_err0.03_04']
lcol  = [50,100,130,200]

iz=1   
kmax = 0.04-0.01*iz
print,iz,kmax
   
for isuff=0,n_elements(lsuff)-1 do begin
   suff = lsuff[isuff]
   for isim =0,9 do begin
      for ig=0,4  do begin
         t  = dir + 'PS'+root+nx[iz]+'_z'+namez[iz]+suff+'_'+strtrim(isim,2)+'_G'+strtrim(ig,2)+'_wngal.txt' 
         print,t
         p = read_ascii(t, template =  TEMP_POW_SPEC_4col) 
         if (isuff + isim + ig eq 0) then  myp = dblarr(n_elements(p.(0)),n_elements(lsuff))

         if (isim + ig eq 0) then myp[*,isuff] = p.(2) else myp[*,isuff] +=  p.(2)
         oplot,p.(0),p.(2),th=2,col=lcol[isuff]
      endfor
   endfor
endfor
myp /= 50.
for isuff=0,n_elements(lsuff)-1 do oplot,p.(0),myp[*,isuff],th=3

lerr=['Gaussian 0.00','Gaussian 0.02','Gaussian 0.03','Gaussian 0.04']
legend,[lerr],li=0,col=[lcol],th=3,box=1,/fill,/right,/center,charsize=1.25
;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/sn_gauss_ps.jpg' ,tvrd(true=3),true=3

window,4
sn=mean(myp[where(p.(0) gt 0.15),*])

plot,p.(0),1./p.(0)/alog(myp[*,3]/sn),/xs,/ys,xra=[0.02,.2],xtit=myxt,ytit='PS',/nodata,xma=[9,1],yma=[4,.5];,yra=[0,1]
for isuff=0,n_elements(lsuff)-1 do oplot,p.(0),1./p.(0)/alog(myp[*,isuff]/sn),col=lcol[isuff],th=4
legend,[lerr],li=0,col=[lcol],th=3,box=1,/fill,/right,/center,charsize=1.25

stop
end

