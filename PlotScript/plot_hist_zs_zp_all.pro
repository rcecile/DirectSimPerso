PRO plot_hist_zs_zp_all,doplot,isuff
; to do after rangez_zp.pro

dir='/sps/lsst/data/rcecile/Planck_BAO2/'
loadct,39

suff = ['_errP','_errPBDT9','_errPBDT8','_err0.03']
mytit = [', no quality cut',' with BDT 90% cut',' with BDT 80% cut','with Gaussian error']
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

filename='histo_zs_zp_lfZuccaAllFalse'+suff[isuff]+'.sav'
print,'restore ',filename
restore,dir+filename

lhp = alog10(double(all_hist))

ws = 660
myc = findgen(12)/2.
mycdiv = 42.

myc = findgen(6)
mycdiv = 50.

im=contour(lhp,z,z,xtit='Spectroscopic redshift z_s',ytit='Photometric redshift z_p'+mytit[isuff],/fill,POSITION=[0.13,0.2,0.95,0.99],FONT_SIZE=15,$
           RGB_INDICES=myc*mycdiv, C_VALUE=myc,RGB_TABLE=39, DIM=[ws,ws*1.1],xticklen=1,yticklen=1,xsubticklen=.025,ysubticklen=.025)
cb = COLORBAR(TARGET=im,POSITION=[0.1,0.065,0.95,0.1],  TITLE='log(Ngal) with redshift bin width = 0.01',FONT_SIZE=15)


if (doplot eq 1) then im.Save, "/sps/lsst/users/rcecile/Fig/stat_zs_zp"+suff[isuff]+".pdf",xmargin=0.2, ymargin=0, /CENTIMETERS, page_size=[6.8, 7.5], height=7.50

stop
   
END

;w.save, 'publication_vector_output.pdf',xmargin=0, ymargin=0, /CENTIMETERS, page_size=[figWidthCm, figHeightCm], height=figHeightCm
