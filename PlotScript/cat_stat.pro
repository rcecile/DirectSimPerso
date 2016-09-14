PRO cat_stat,cas
;cat_stat,"" : for the plot
;cat_stat,"4" : for check

loadct,39
!p.charsize=2
!p.thick=3

dir='/sps/lsst/data/rcecile/'
if (strlen(cas) eq 0) then dir = dir+"TJP_BAO/" else  dir = dir+"TJP_noBAO/"
if (strlen(cas) eq 0) then cas = "" else  cas = "_"+simu
print,dir

suff=["_ZONLY",""];,"_errp","_errpodds"]
nsuff=n_elements(suff)
lcol  = [0,80, 150, 240, 210]

loadct,39
!p.charsize=2
!p.thick=3
;recopie une tranche sauf les catastrophiques (> 3 sigma de G 0.03), regroupes
nSlice=70
nSlice=65
     
n=dblarr(nSlice,nsuff)
z=dblarr(nSlice+1)
z[nSlice]=10

for is=0,nsuff-1 do begin
   for i=0,nSlice-1 do begin
      name=dir+"cat_G2"+cas+"_Slice"+strtrim(i,2)+suff[is]+".fits"
      h=headfits(name,ext=1) 
      if (is eq 1) then begin
         zmin = sxpar(h,'ZCAT_MIN')
         zmax = sxpar(h,'ZCAT_MAX')
         z[i] = (zmin+zmax)/2.
      endif
      nn = sxpar(h,'NAXIS2')
      n[i,is] = double(nn)
      print,is,i,"",name,nn
   endfor
   if (is le 1 ) then continue
   name=dir+"cat_G2"+cas+"_AllSlice"+suff[is]+"_cataLT10.fits"
   m=mrdfits(name,1)
   for i=0,nSlice-1 do begin
      ok = where(m.(6) ge z[i] and m.(6) lt z[i+1],nok)
      n[i,is] = n[i,is] + double(nok)
      print,'cata ',i,nok
   endfor

endfor

mytot=dblarr(nsuff)
for is=0,nsuff-1 do mytot[is] = total(n[*,is])
myformat='(E8.2)'
myout = string(mytot,format=myformat)

!p.multi=0
plot,z[0:nSlice-1],n[0:nSlice-1,0],/xs,/ys,xtit='Redshift',ytit='Nb of galaxies / slice 80 Mpc thick',psym=10,/yl,/nodata,yra=[1e5,2e9]
for is=0,nsuff-1 do oplot,z,n[*,is],col=lcol[is]
what=['No magnitude cut             ',$
      'Golden selection             ',$
      'Golden selection + photoz    ',$
      'Golden selection + podds cut']+':  N gal  ='+myout
legend,what,col=lcol,line=0,box=1,/fill,/right,/top,charsize=1.5


stop

end

