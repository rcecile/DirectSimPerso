PRO plot_bao_sigma,saveplot
!p.charsize=1.5

h= 0.679
h3 = h*h*h
s=150.6
restore,'temp_ascii_new.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO2_PS/"

z=[0.5,0.9, 1.5, 1.3]
namez = ['0.5','0.9','1.3']
nx= ['_120','_225', '_300']
nz = n_elements(namez)

nsuff=['','_err0.03','_errP','_errPBDT9','_errPBDT8']
n = n_elements(nsuff)
lerr=['spectroZ','Gaussian 0.03','photoZ','photoZ BDT 90%','photoZ BDT 80%']
myformat ='(F7.1)'
kmin=[4,3,2]
kmax=[20,20,20,20,20,    20,15,15,15,15,      20,12,12,12,12]

kmax=[15,15,15,15,15,    20,15,15,15,15,      20,12,12,12,12]

lcol  = [95,210,35,135,  120,0]
loadct,12

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (saveplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/users/rcecile/Fig/plot_sigma.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif else window,0


plot,[0.5,19.5],[0,1.75],/nodata,/xs,/ys,ytit='Error on the BAO scale [%]',xticknam=replicate(" ",20),xma=[7.5,.5],yma=[.75,.5];,yticklen=1
for iz=0,2 do begin
   print,'iz = ', namez[iz]
   for ierr=0,n-1 do begin
      sk ='_k0.0' +strtrim(kmin[iz],2)+'_0.'+strtrim(kmax[iz*5+ierr],2)
      ffit = dir + 'fit_lfZuccaAllFalse'+ sk +nx[iz]+'_z'+namez[iz]+nsuff[ierr]+'_result.txt'
      pfitw = read_ascii(ffit, template =  TEMP_POW_SPEC_FIT_RES)  
      if (ierr eq 0) then sigzs = (pfitw.(2))[0]
      if (ierr ge 2) then print,'    ierr',nsuff[ierr], '      krange = ',sk,'         Chi2 = ',(pfitw.(4))[0],$
             ' sig = ',(pfitw.(2))[0],'  --------  sigZp/sigZs = ',(pfitw.(2))[0]/sigzs,' diff [sig] = ',((pfitw.(0))[0]-s)/(pfitw.(2))[0]
      idx = iz*(n+2)+ierr+1

      oplot,[idx,idx],[0,(pfitw.(1))[0]]/s*100.,th=10,col=lcol[ierr]
   endfor
endfor

;xyouts, 1,139.6,'grids @ z=0.5',charsize=2
;xyouts, 8,139.6,'grids @ z=0.9',charsize=2
;xyouts,15,139.6,'grids @ z=1.5',charsize=2
legend,'z = '+namez[0],charsize=2,/top,/left,/box
legend,'z = '+namez[1],charsize=2,/top,/center,/box
legend,'z = '+namez[2],charsize=2,/top,/right,/box

legend,lerr,li=0,col=lcol,box=1,/fill,position=[7,1.3],charsize=1.5,th=3

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

end
