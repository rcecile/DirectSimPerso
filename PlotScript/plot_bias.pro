PRO plot_bias,doplot,savettot

dir='/sps/lsst/data/rcecile/Planck_BAO2/'
suff = ['_errP','_errPBDT9','_errPBDT8']
loadct,12
lcol  = [0,210,105,  135, 120]
!p.thick=2
!p.charsize=1.5
lcol  = [35,105,210,  135, 120]
lcol  = [35,135, 120]

!x.margin=[8.5,0.5]

!p.multi=[0,1,3]

if (doplot) then begin  
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/users/rcecile/Fig/plot_bias.eps', /PORTRAIT,/COLOR,XSIZE=6.6,YSIZE=6.6,FONT_SIZE=4
endif

ni=[17,35,63,40]

if (savettot eq 1) then begin
   binsize=0.001
   binmax=0.6
   isuff=0
   for i=0,3 do begin
      for isuff=0,2 do begin
         file='cat_lfZuccaAllFalse_zOrdZP_Slice'+strtrim(ni[i],2)+suff[isuff]+'.fits'
         m=mrdfits(dir+file,1,h,col=['ZS','ZP'])
         print,file,mean(m.(1))
         t = (m.(0)-m.(1))/(1.+m.(1))
         h=histogram(t,min=binmax*(-1.),max=binmax,bin=binsize)
         x=findgen(n_elements(h))*binsize-binmax
         plot,x,h,/xs,/ys,xra=[-0.02,0.02]
         if (i eq 0 and isuff eq 0) then ttot=dblarr(n_elements(h),4,3)
         ttot[*,i,isuff] = h
      endfor
   endfor
   save,ttot,x,file=dir+'ttot.sav'
endif

restore,dir+'ttot.sav'
mtext=['photoZ','photoZ, BDT 90% cut','photoZ, BDT 80% cut']

bmin=-0.6
bmax=0.6
bstep=0.0003
myformat = '(F4.2)'
!p.multi=[0,2,2]
myx=3.5
myxmar1=[myx, 0.0, myx, 0.]
myxmar2=[0.0, myx, 0.0, myx]
myymar1=[1.5, 1.5, 4.0, 4.0]
myymar2=[1., 1., -1.5, -1.5]
     
emptytick=replicate(" ",10)

for i=0,3 do begin
   for isuff=0,2 do begin
      h =  ttot[*,i,isuff]
      if (i ge 2) then myxtit='(z_s-z_p)/(1+z_p)' else myxtit=' '
      if (i lt 2) then myxtickname = emptytick else myxtickname=''
      if ( (i mod 2) eq 0) then myytickname = '' else myytickname=emptytick
      if (isuff eq 0) then $
         plot,x,h/max(h),/xs,/ys,xra=[-0.09,0.09],/yl,yra=[0.0501,1],xmar=[myxmar1[i],myxmar2[i]],ymar=[myymar1[i],myymar2[i]],xtit=myxtit,ytickname=myytickname,xtickname=myxtickname
      oplot,x,h/max(h),col=lcol[isuff]
      restore,dir+'statZp_lfZuccaAllFalse'+suff[isuff]+'.sav'
  ;    oplot,[x[where(h eq max(h))],x[where(h eq max(h))]],[0,max(h)*2],th=4
      oplot,[pbias[ni[i]],pbias[ni[i]]],[0.001,1],th=2,col=lcol[isuff]
      print,i,ni[i],'  ',suff[isuff],pbias[ni[i]]
      mys = string(zSLICE[ni[i]],format=myformat)
      legend,'z_p='+mys,charsize=1.25,/fill,/right,/top
   endfor
   oplot,[0,0],[0.001,1],th=2,li=2
   if ( (i mod 2) eq 1) then AXIS, YAXIS=1, ys=1

endfor

;legend,mtext,lin=0,col=lcol,box=1,/fill,/left,/top,charsize=1.5

if (doplot eq 1) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

stop


END
