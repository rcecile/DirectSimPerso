PRO prepa_hist_zs_zb_interquart

dir='/sps/lsst/data/rcecile/TJP_BAO/'
nslice =70

pinterquart = dblarr(nslice)
pbias = dblarr(nslice)
poutlier = dblarr(nslice)
dinterquart = dblarr(nslice)
dbias = dblarr(nslice)
doutlier = dblarr(nslice)
zstat = dblarr(nslice)

; better to do this part in qlogin mode
for i = 0,nslice-1 do begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   file='cat_G2_Slice'+strtrim(i,2)+'_errPBDT.fits'
   print,file
   m=mrdfits(dir+file,1,h,col=['ZS','ZP'])
   ok = where(m.(1) lt 9)
   t = (m.(1)-m.(0))/(1.+m.(0))
   file='cat_G2_Slice'+strtrim(i,2)+'_errPBDT_cata.fits'
   print,file
   m=mrdfits(dir+file,1,h,col=['ZS','ZP'])
   ok = where(m.(1) lt 9)
   t = [t,((m.(1)-m.(0))/(1.+m.(0)))[ok]]

   RSTAT,t, Med, Hinge1, Hinge2
   dbias[i] = Med
   dinterquart[i] = Hinge2-Hinge1
   out = where(abs(t) gt 3.*dinterquart[i],nout)
   doutlier[i] = double(nout)/double(n_elements(t))*100.

   print,'-------------------------------------------------------------------------------------------------'
   print,'stat ',i,zstat[i],'       bias ',dbias[i],'    sigma ',dinterquart[i],'   outlier ',doutlier[i]
   print,'-------------------------------------------------------------------------------------------------'
   save,zstat,pinterquart,pbias,poutlier,dinterquart,dbias,doutlier,file=dir+'histo_zs_zb_requirement_med.sav'
endfor

END
