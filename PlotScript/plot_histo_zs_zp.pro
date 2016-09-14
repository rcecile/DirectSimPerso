PRO plot_histo_zs_zp,doplot
; plot_histo_zs_zp,0 --> complete the histogram (cannot be done in qlogin mode)
; plot_histo_zs_zp,1 --> plot
; plot_histo_zs_zp,2 --> save the plot

dir='/sps/lsst/data/rcecile/TJP_BAO/'
   loadct,39

; to do interactvely
if (doplot eq 0) then begin
   restore,dir+'/histo_zs_zp_0.01_diag.sav'
   file='cat_G2_AllSlice_errP_cataLT10.fits'
   print,file
   m=mrdfits(dir+file,1,h,col=['ZS','ZP'])
   h = hist_2d(m.(0),m.(1),bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
   hp = hp + h

   h = histogram( (m.(0)-m.(1))/(1.+m.(0)) ,bin=diff_bin,min=-1*diff_max,max=diff_max)
   hpdiff = hpdiff + h
   
   file='cat_G2_AllSlice_errPpodds_cataLT10.fits'
   print,file
   m=mrdfits(dir+file,1,h,col=['ZS','ZP'])
   h = hist_2d(m.(0),m.(1),bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
   hd = hd + h
   
   h = histogram( (m.(0)-m.(1))/(1.+m.(0)) ,bin=diff_bin,min=-1*diff_max,max=diff_max)
   hddiff = hddiff + h

   print,'Save up to slices CATA',total(hp),total(hd)
   save,hp,hd,hpdiff,hddiff,z,binz,zmax,diff_bin,diff_max,i,file=dir+'histo_zs_zp_0.01.sav'
endif
   s(z-z[i])/(1.+z[i]) gt 0.09),i])
 

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (doplot ge 1) then begin

   restore,dir+'/histo_zs_zp_0.01.sav'
   lhp = alog10(double(hp))
   lhd = alog10(double(hd))
   
   loadct,39
   
   ws = 600
   myc = findgen(12)/2.
   mycdiv = 42.
   
   myc = findgen(6)
   mycdiv = 50.
    
   im=contour(lhp,z,z,xtit='spectroZ',ytit='photoZ',tit='raw photoZ',/fill,POSITION=[0.12,0.2,0.97,0.93],FONT_SIZE=15,$
              RGB_INDICES=myc*mycdiv, C_VALUE=myc,RGB_TABLE=39, DIM=[ws,ws*1.225],xticklen=1,yticklen=1,xsubticklen=.025,ysubticklen=.025)
   cb = COLORBAR(TARGET=im,POSITION=[0.1,0.08,0.95,0.12],  TITLE='log(Ngal) with redshift bin width = 0.01',FONT_SIZE=15)
   
   
   if (doplot eq 2) then im.Save, "/sps/lsst/dev/rcecile/Fig/stat_zs_zp.pdf"
   
   im=contour(lhd,z,z,xtit='spectroZ',ytit='photoZ',tit='photoZ with ODDS cut',/fill,POSITION=[0.12,0.2,0.97,0.93],FONT_SIZE=15,$
              RGB_INDICES=myc*mycdiv, C_VALUE=myc,RGB_TABLE=39, DIM=[ws,ws*1.225],xticklen=1,yticklen=1,xsubticklen=.025,ysubticklen=.025)
   cb = COLORBAR(TARGET=im,POSITION=[0.1,0.08,0.95,0.12],  TITLE='log(Ngal) with redshift bin width = 0.01',FONT_SIZE=15)
   
   if (doplot eq 2) then im.Save, "/sps/lsst/dev/rcecile/Fig/stat_zs_zd.pdf"
 

endif

END

