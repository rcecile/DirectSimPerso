PRO cat_hist_type,doplot
;cat_hist_type,0 : to compute the histogram (long), must be done first
;cat_hist_type,1 : to do the plot (fast)
loadct,39
!p.charsize=2.
!p.thick=3
surf = !pi / !pi * 180. / !pi * 180. *60. * 60. ; 1/4 du ciel en arc-min2 + changement de bining
lcol  = [0,80]


suff = 'lfZuccaAllFalse'

dir='/sps/lsst/data/rcecile/Planck_BAO2/'

if (doplot eq 0) then begin

   nSlice=70

   z=findgen(61)/60.*3.
   ztype=dblarr(360,61)
   zngal=dblarr(61)
   z123=dblarr(61,3)

   istart=0
   if (istart gt 0) then restore,'/sps/lsst/data/rcecile/Planck_BAO2/cat_stat_type'+suff+'.save'


   for is=istart,nSlice-1 do begin

      name=dir+"cat_"+suff+"_Slice"+strtrim(is,2)+".fits"
      print,name
      
      m=mrdfits(name,1,h,col=['TYPE','ZS'])
      for it=1,3 do begin
         ok = where(fix(m.(0)) eq it,nok)
         h = histogram( (m.(1))[ok] ,min=0,bin=0.05,max=3)
         z123[*,it-1]  = z123[*,it-1] + 1.d*h
        endfor

      for i=0,60 do zngal[i] = total(z123[i,*])
         
      window,10,xs=600,ys=600
      plot,[0,3],[1e3,2e7],/xs,/ys,psy=10,xra=[0.2,2.45],/yl,/nodata,yra=[0.005,2.5],xticklen=1,yticklen=1
      oplot,z,zngal[*]
      for i=0,2 do oplot,z,z123[*,i]/surf,col=100+i*50
      legend,['All','Early','Late','StarBurst'],line=0,col=[0,100,150,200],box=1,/fill,/right,/top,charsize=1.5
      save,z,ztype,zngal,z123,file='/sps/lsst/data/rcecile/Planck_BAO2/cat_stat_type'+suff+'.save'

   endfor

   stop
    
endif

if (doplot eq 2) then begin

   print,'  '
   for is=0,n_elements(lsuff)-1 do begin

      restore,'/sps/lsst/data/rcecile/Planck_BAO2/cat_stat_type'+lsuff[is]+'.save'
      mytot=dblarr(4)
      for i=0,2 do mytot[i] = total(z123[*,i])/surf ;/0.003*0.05
      mytot[3] = total(zngal)/surf; /0.003*0.05

      print,'  '
      print,lsuff[is],mytot
   endfor
   print,'  '
endif

if (doplot eq 1) then begin
   restore,'/sps/lsst/data/rcecile/Planck_BAO2/cat_stat_type'+suff+'.save'
   window,isuff+0,xs=600,ys=600
   if (max(z123)/surf gt 1) then mymin = 0.005 else mymin= 0.0005
   if (max(z123)/surf gt 1) then mymax = 2.5 else mymax=0.25
   plot,z,zngal[*]/surf,/xs,/ys,xra=[0.2,2.45],/yl,yra=[mymin,mymax],tit=suff,xtit='Redshift',ytit='Ngal',xticklen=1,yticklen=1
   for i=0,2 do oplot,z,z123[*,i]/surf,col=100+i*50
   legend,['All','Early','Late','StarBurst'],line=0,col=[0,100,150,200],box=1,/fill,/right,/top,charsize=1.5
;   write_jpeg,'/sps/lsst/users/rcecile/BAO_InOut/check_ngal'+suff+'.jpg' ,tvrd(true=3),true=3
   window,isuff+2
   
   plot,z,zngal[*],/xs,/ys,xra=[0.2,1.15],/yl,yra=[0.01,1],tit=suff,xtit='Redshift',ytit='Ngal'
   for i=0,2 do oplot,z,z123[*,i]/zngal,col=100+i*50
   legend,['All','Early','Late','StarBurst'],line=0,col=[0,100,150,200],box=1,/fill,/right,/bottom,charsize=1.5
 ;  write_jpeg,'/sps/lsst/users/rcecile/BAO_InOut/check_ngal_frac'+suff+'.jpg' ,tvrd(true=3),true=3
   
   read,xx
   
   restore,'temp_ascii_new.sav'
   all = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_LF_Dahlen_ZONLY.txt' , template =  TEMP_SEL_FUNC)
   gold = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_LF_Dahlen_specz_nofz.txt' , template =  TEMP_SEL_FUNC)
   
   
   nonly = dblarr(n_elements(z))
   nall = nonly
   for i=0,n_elements(z)-2 do begin
      ok = where(all.(0) ge z[i] and all.(0) lt z[i+1])
      nonly[i] = total( (all.(1))[ok] )
      nall[i] = total( (gold.(1)*all.(1))[ok] ) 
   endfor
   
   myformat='(E8.1)'
   ok = where(z ge 0.2 and z le 2.45)
   myout = string([total(nonly[ok]),total(nall[ok])],format=myformat)
   
   !p.multi=[0,1,2]
   
   
;window,0
   plot,z,alog10(smooth(nonly,5)),/xs,/ys,ytit='log10(Ngal)',yra=[5,8.9],xmar=[7,1],ymar=[0,2],xtickname=replicate(" ",10),xra=[0.2,2.45] ;,xticklen=1
   oplot,z,alog10(smooth(nall,5)),col=lcol[1]
   
   what=['No magnitude cut',$
         'Golden sample']+', Ngal ='+myout
   legend,what,col=[lcol[0],lcol[1]],line=0,box=1,/fill,/left,/bottom,charsize=1.75
   
   n_type=3
   plot,z,alog10([1,150]),xtit='Redshift',/xs,/ys,ytit='log10(Ngal)',psym=10,yra=[2,7.9],xmar=[7,1],ymar=[4,0],/nodata,xra=[0.2,2.45] ;,xticklen=1
   for i=0,n_type-1 do oplot,z,alog10(z123[*,i]),col=lcol[1],lin=2-i
   
   what=['Early type    ',$
         'Late type     ',$
         'StarBurst type' ]
   legend,what,col=lcol[1]+150,line=[2-indgen(n_type)],box=1,/fill,/left,/bottom,charsize=1.75
   
;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/check_ngal.jpg' ,tvrd(true=3),true=3
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
   
   !p.multi=0
   
   plot,z, (smooth(nall,5))/surf /0.03*0.05,/xs,/ys,/yl,yra=[0.04,2],xra=[0.2,2.45],xtit='redshift',ytit='Ngal/arc-min2/bin en z 0.05',th=5
   for i=0,n_type-1 do oplot,z,(ztype[*,i])/surf/0.03*0.05,lin=5-i*2
   
   what=['Early type    ',$
         'Late type     ',$
         'StarBurst type' ]
   legend,what,col=0,line=[5-indgen(n_type)*2],box=1,/fill,/right,/top,charsize=1.75
   
;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/check_ngal_zSch.jpg' ,tvrd(true=3),true=3
endif
stop
end
