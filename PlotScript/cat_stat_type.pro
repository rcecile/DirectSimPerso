PRO cat_stat_type,doplot
;cat_stat_type,0 : to compute the histogram (long), must be done first
;cat_stat_type,1 : to do the plot (fast)
loadct,39
!p.charsize=2.
!p.thick=4
lcol  = [0,80]

if (doplot eq 0) then begin

   dir='/sps/lsst/data/rcecile/TJP_BAO/'
   nSlice=70
   
   nonly=dblarr(nSlice)
   n_type=3
   ntype=dblarr(nSlice,n_type)
   z=dblarr(nSlice)

   for i=0,nSlice-1 do begin
      name=dir+"cat_G2_Slice"+strtrim(i,2)+"_ZONLY.fits"
      h = headfits(name,ext=1)
      nn = sxpar(h,'NAXIS2')
      nonly[i] = double(nn)
      PRINT,NAME,NN
   endfor
   
   for i=0,nSlice-1 do begin
      name=dir+"cat_G2_Slice"+strtrim(i,2)+".fits"
      m=mrdfits(name,1,h,col=['TYPE'])
      zmin = sxpar(h,'ZCAT_MIN')
      zmax = sxpar(h,'ZCAT_MAX')
      z[i] = (zmin+zmax)/2.
      for it=0,n_type-1 do begin
         ok = where(fix(m.(0)) eq it+1,nok)
         ntype[i,it] = double(nok)
      endfor
      PRINT,NAME,ntype[i,*]
      m=0
   endfor

   mytot=dblarr(n_type+2)
   mytot[0] = total(nonly)
   for i=0,n_type-1 do mytot[i+2] = total(ntype[*,i])
   mytot[1] = total(mytot[2:2+n_type-1])
   myformat='(E8.2)'
   myout = string(mytot,format=myformat)

   nall=dblarr(nSlice)
   for i=0,nSlice-1 do nall[i] = total(ntype[i,*])

   save,z,nonly,ntype,n_type,nall,myout,file='/sps/lsst/data/rcecile/TJP_BAO/cat_stat_type.save'

    
endif


if (doplot eq 2) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/cat_stat_type.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=8.8,FONT_SIZE=4

endif

if (doplot ge 1) then restore,'/sps/lsst/data/rcecile/TJP_BAO/cat_stat_type.save'

h=dloscom(z)
ang=1600.*8. / h / sqrt(2.)
correct = where(ang lt 120./180*!pi)
nonlyc= nonly
nallc = nall
nonlyc[correct] = nonlyc[correct] / ang[correct]*120./180*!pi
nallc[correct] = nallc[correct] / ang[correct]*120./180*!pi


!p.multi=[0,1,2]
plot,z,alog10(nonly),/xs,/ys,ytit='log10(Ngal)',yra=[5,9.5],xmar=[7,1],ymar=[0,2],xtickname=replicate(" ",10);,xticklen=1
oplot,z,alog10(nall),col=lcol[1]
oplot,z,alog10(nonlyc),li=2
oplot,z,alog10(nallc),col=lcol[1],li=2
what=['No magnitude cut      ',$
      'Golden selection        ']+':  N gal  ='+myout
legend,what,col=[lcol[0],lcol[1]],line=0,box=1,/fill,/left,/bottom,charsize=1.5

plot,z,alog10([1,100]),xtit='Redshift',/xs,/ys,ytit='log10(Ngal)',psym=10,yra=[1.5,7.9],xmar=[7,1],ymar=[4,0],/nodata;,xticklen=1
for i=0,n_type-1 do oplot,z,alog10(ntype[*,i]),col=lcol[1],lin=2-i
what=['Early type    ',$
      'Late type     ',$
      'StarBurst type' ]

legend,what,col=lcol[1],line=[2-indgen(n_type)],box=1,/fill,/center,/bottom,charsize=1.75

if (doplot eq 2) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif



end

