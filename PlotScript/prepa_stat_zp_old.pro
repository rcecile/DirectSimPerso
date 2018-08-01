PRO prepa_stat_zp,isuff

dir='/sps/lsst/data/rcecile/Planck_BAO2/'

nslice = 100
suff = ['_errP','_errPBDT9','_errPBDT8']
pinterquart = dblarr(nslice)
pbias = dblarr(nslice)
poutlier = dblarr(nslice)
poutlier08 = dblarr(nslice)
poutlier10 = dblarr(nslice)
poutlier12 = dblarr(nslice)
z_max = 3.

z = findgen(nslice+1)*z_max/nslice

; better to do this part in qlogin mode
for i = 3,nslice-2 do begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   file='cat_AllzOrdZP_Slice'+strtrim(i,2)+suff[isuff]+'.fits'
   print,file
   m=mrdfits(dir+file,1,h,col=['ZS','ZP'])
   t = (m.(1)-m.(0))/(1.+m.(0))
   
   RSTAT,t, Med, Hinge1, Hinge2
   pbias[i] = Med
   pinterquart[i] = Hinge2-Hinge1
   out = where(abs(t) gt 0.15,nout)
   poutlier[i] = double(nout)/double(n_elements(t))*100.

   out = where(abs(t) gt 0.08,nout)
   poutlier08[i] = double(nout)/double(n_elements(t))*100.
   out = where(abs(t) gt 0.10,nout)
   poutlier10[i] = double(nout)/double(n_elements(t))*100.
   out = where(abs(t) gt 0.12,nout)
   poutlier12[i] = double(nout)/double(n_elements(t))*100.

   print,'As a function of z_s --------------------------------------------------'
   print,'-------------------------------------------------------------------------------------------------'
   print,'stat ',i,'       bias ',pbias[i],'    sigma ',pinterquart[i],'   outlier ',poutlier[i],' et ',poutlier08[i],poutlier10[i],poutlier12[i]
   print,'-------------------------------------------------------------------------------------------------'
   save,z,pbias,pinterquart,poutlier,poutlier08,poutlier10,poutlier12,file=dir+'stat_zp_1pzs_all'+suff[isuff]+'.sav'
endfor

END
