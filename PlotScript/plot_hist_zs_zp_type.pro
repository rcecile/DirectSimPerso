PRO plot_hist_zs_zp_type,doplot,isuff
; to do after rangez_zp.pro

dir='/sps/lsst/data/rcecile/Planck_BAO2/'
loadct,39
ws = 500

suff = ['_errP','_errPBDT9','_errPBDT8','_err0.03']
mytit = [' with BDT 90% cut',' with BDT 80% cut',', no quality cut','with Gaussian error']
mytitt = ['Early','Late','StarBurst']+' type galaxies'
myname= ['Early','Late','StarBurst']
mytitt = myname+' type galaxies'

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

filename='histo_zs_zp_lfZuccaAllFalse_type'+suff[isuff]+'.sav'
print,'restore ',filename
restore,dir+filename

lall = all_hist0+all_hist1+all_hist2

n = n_elements(lall[0,*])
lha0 = all_hist0
lha1 = all_hist1
lha2 = all_hist2
lh0 = all_hist0
lh1 = all_hist1
lh2 = all_hist2
for i=0,n-1 do lha0[i,*] = lh0[i,*]/total(lall[i,*])*100.
for i=0,n-1 do lha1[i,*] = lh1[i,*]/total(lall[i,*])*100.
for i=0,n-1 do lha2[i,*] = lh2[i,*]/total(lall[i,*])*100.
for i=0,n-1 do lh0[i,*] = lh0[i,*]/total(lh0[i,*])*100.
for i=0,n-1 do lh1[i,*] = lh1[i,*]/total(lh1[i,*])*100.
for i=0,n-1 do lh2[i,*] = lh2[i,*]/total(lh2[i,*])*100.

myc = [findgen(10)/10,1+findgen(5)*10]
ic0=70
mycc = [ic0+findgen(10)*15,(ic0+10*15.+findgen(5)*10)]

for it=0,2 do begin
   if it eq 0 then lh = lh0
   if it eq 1 then lh = lh1
   if it eq 2 then lh = lh2
   im=contour(lh,z,z,xtit='Spectroscopic redshift z_s',ytit='Photometric redshift z_p'+mytit[isuff],/fill,POSITION=[0.13,0.2,0.95,0.99],FONT_SIZE=10,$
              RGB_INDICES=mycc, C_VALUE=myc,RGB_TABLE=39, DIM=[ws,ws*1.1],xticklen=1,yticklen=1,xsubticklen=.025,ysubticklen=.025)
   cb = COLORBAR(TARGET=im,POSITION=[0.1,0.065,0.95,0.1],  TITLE='Distribution of z_p(z_s) [% of this type] for  '+mytitt[it],FONT_SIZE=10)

   if it eq 0 then lh = lha0
   if it eq 1 then lh = lha1
   if it eq 2 then lh = lha2
   im=contour(lh,z,z,xtit='Spectroscopic redshift z_s',ytit='Photometric redshift z_p'+mytit[isuff],/fill,POSITION=[0.13,0.2,0.95,0.99],FONT_SIZE=10,$
              RGB_INDICES=mycc, C_VALUE=myc,RGB_TABLE=39, DIM=[ws,ws*1.1],xticklen=1,yticklen=1,xsubticklen=.025,ysubticklen=.025)
   cb = COLORBAR(TARGET=im,POSITION=[0.1,0.065,0.95,0.1],  TITLE='Distribution of z_p(z_s) [% of all] for  '+mytitt[it],FONT_SIZE=10)
   

   if (doplot eq 1) then im.Save, "/sps/lsst/users/rcecile/Fig/stat_zs_zp"+suff[isuff]+".pdf",xmargin=0.2, ymargin=0, /CENTIMETERS, page_size=[6.8, 7.5], height=7.50

endfor

stop
   
END

;w.save, 'publication_vector_output.pdf',xmargin=0, ymargin=0, /CENTIMETERS, page_size=[figWidthCm, figHeightCm], height=figHeightCm
