PRO error_distrib

loadct,39
!p.charsize=2
!p.thick=3
lcol  = [150, 240, 210]

dir='/sps/lsst/data/rcecile/TJP_BAO/'

suff=["_err0.03","_errP","_errPpodds"]
nsuff = n_elements(suff)

hmin=-0.09
hmax=0.09
hbin=0.001
x=findgen((hmax-hmin)/hbin+1)*hbin+hmin
nx=n_elements(x)
hh=dblarr(n_elements(x),nsuff,6)

z=dblarr(6)
for i=10,60,10 do begin
   for is=0,nsuff-1 do begin
      name=dir+"cat_zOrd_Slice"+strtrim(i,2)+suff[is]+".fits"
      m=mrdfits(name,1,hd)
      h=histogram((m.(6)-m.(3))/(1.+m.(3)),min=hmin,max=hmax,bin=hbin)
      h = h/double(h[nx/2+1])
      zmin = sxpar(hd,'ZCAT_MIN')
      zmax = sxpar(hd,'ZCAT_MAX')
      z[i/10-1]=(zmin+zmax)/2.
      plot,x,h,tit=strtrim((zmin+zmax)/2.)
      hh[*,is,i/10-1] = h
   endfor
endfor
print,z

save,z,hh,file='/sps/lsst/data/rcecile/TJP_BAO/error_distrib.sav'

END
