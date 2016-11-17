PRO plot_hist_zs_zp_all,doplot,isuff

dir='/sps/lsst/data/rcecile/TJP_BAO/'
loadct,39

suff = ['_errPpodds','_errPBDT','_errP','_err0.03']
nsuff = ['_errPpodds','_errPBDT','_errP','_errG']
mytit = ['',' with ODDS cut',' with BDT cut','with Gaussian error']
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

restore,dir+'histo_zs_zp'+suff[isuff]+'.sav'
lhp = alog10(double(all_hist))
binz=0.01
zmax=3.
z=findgen(zmax/binz+1)*binz

ws = 600
myc = findgen(12)/2.
mycdiv = 42.

myc = findgen(6)
mycdiv = 50.

im=contour(lhp,z,z,xtit='spectroZ',ytit='photoZ'+mytit[isuff],/fill,POSITION=[0.12,0.2,0.97,0.93],FONT_SIZE=15,$
           RGB_INDICES=myc*mycdiv, C_VALUE=myc,RGB_TABLE=39, DIM=[ws,ws*1.225],xticklen=1,yticklen=1,xsubticklen=.025,ysubticklen=.025)
cb = COLORBAR(TARGET=im,POSITION=[0.1,0.08,0.95,0.12],  TITLE='log(Ngal) with redshift bin width = 0.01',FONT_SIZE=15)

if (doplot eq 1) then im.Save, "/sps/lsst/dev/rcecile/Fig/stat_zs_zp"+nsuff[isuff]+".pdf"
   
 print,total(all_hist)
   restore,dir+'/histo_zs_zp_0.01.sav'
 print,total(hp)
 print,total(hd)

END
