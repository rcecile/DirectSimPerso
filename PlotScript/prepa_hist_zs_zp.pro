PRO prepa_plot_concat2

dir='/sps/lsst/data/rcecile/PaperBAO_cell8/'

restore,dir+'histo_zs_zp_concat.sav'
restore,dir+'histo_zs_zp_concat_diag.sav'
i_min = i+1
hdiag = h*0.
h2diag = h*0.

binz=0.05
zmax=3.

for i = i_min,63 do begin
   file='cat_gold_Slice'+strtrim(i,2)+'_errP.fits'
   print,file
   m=mrdfits(dir+file,1,h,col=['zs','zP'])
   h = hist_2d(m.(0),m.(1),bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
   hdiag = h2diag + h

   file='cat_gold_Slice'+strtrim(i,2)+'_errPpodds.fits'
   print,file
   m=mrdfits(dir+file,1,h,col=['zs','zP'])
   h = hist_2d(m.(0),m.(1),bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
   h2diag = h2diag + h

   print,'Save up to slice ',i
   save,h,h2,hdiag,h2diag,x,zpp,zpm,i,file=dir+'histo_zs_zp_concat_diag.sav'
endfor

END

restore,dir+'/histo_zs_zp_concat_diag.sav'
restore,dir+'histo_zs_zp_concat.sav'   
save,h,h2,hdiag,h2diag,x,zpp,zpm,file=dir+'/histo_zs_zp_concat_tot.sav'
