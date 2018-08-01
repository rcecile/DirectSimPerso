PRO plot_mag_type,isuff,doplot
loadct,39
!p.charsize=2.
!p.thick=3
lcol  = [0,80]

lsuff=['LF_ZuccaAllFixed_allFalse_SED0_1_2345','LF_ZuccaAllFixed_allFalse_SED_superGold_BrightOnly']
suff = lsuff[isuff]

dir='/sps/lsst/data/rcecile/Planck_BAO2/'

if (doplot eq 0) then begin
   
   z=findgen(90)
;          MB   z type
   MTz=dblarr(111,90,3)
   Tz=dblarr(90,3)

   istart=4
   nSlice = 93
   
   for is=istart,nSlice do begin
      
      name=dir+"cat_"+suff+"_zOrd_Slice"+strtrim(is,2)+".fits"
      print,name
      m=mrdfits(name,1,h,col=['ZS','TYPE','MB'])
      z[is] = mean(m.(0))
      for it =1,3 do begin
         ok = where(fix(m.(1)) eq it)
         h = histogram( (m.(2))[ok] ,min=-24,bin=0.1,max=-13)
         MTz[*,is,it-1] = h
         Tz[is,it-1] = total( MTz[*,is,it-1])
      endfor
      
      save,z,MTz,Tz,file='/sps/lsst/data/rcecile/Planck_BAO2/cat_stat_type_mag'+suff+'.save'
      plot,[0,3],[1e3,2e7],/xs,/ys,psy=10,xra=[0.2,2.45],/yl,/nodata,yra=[1e3,max(tz)],xticklen=1,yticklen=1
      for i=0,2 do oplot,z,Tz[*,i],col=100+i*50
      legend,['All','Early','Late','StarBurst'],line=0,col=[0,100,150,200],box=1,/fill,/right,/top,charsize=1.5

   endfor


endif

if (doplot eq 1) then begin
   restore,'/sps/lsst/data/rcecile/Planck_BAO2/cat_stat_type_mag'+suff+'.save'
stop
endif

END
