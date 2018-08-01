PRO plot_ps,saveplot

dir="/sps/lsst/data/rcecile/Planck_BAO2_PS/"
root = "_lfZuccaAllFalse"
offset=0

print,'root = ',root
namez = ['0.5','0.9','1.3']
lnamez = ['1.3','0.9','0.5']
nx= ['_120','_225',  '_300']
nx_cell= [120.,225., 300.]
nz_cell= [125.,125.,125.]
cell=8.
namezmean=['0.51', '0.93', '1.36']
lpsym=[0,4,0]
llin = [2,0,0]
lth=[2,6,10]

nameerr = ['', '_err0.03', '_errP','_errPBDT9','_errPBDT8']
lcol  = [95,210,35,135,  120]
lerr=['spectroZ','Gaussian 0.03','photoZ','photoZ BDT 90%','photoZ BDT 80%']

loadct,12

nz = n_elements(namez)
nerr = n_elements(nameerr)
nsim=10

kmax = intarr(n_elements(namez), n_elements(nameerr))
kmax[*,0] = [20,20,15]
kmax[*,1] = [20,20,15]
kmax[*,2] = [20,20,15]
kmax[*,3] = [18,18,12]
kmax[*,4] = [16,16,12]

restore,'/sps/lsst/users/rcecile/BAOProgs/DirectSimPerso/PlotScript/temp_ascii_new.sav'   
myformat = '(F9.7," ",G14.6," ",G14.6," ",G14.6)'
!p.thick=3
!p.charsize=1.5

if (saveplot ) then begin
; current plotting device.
      mydevice = !D.NAME
      SET_PLOT, 'PS'
      DEVICE, FILENAME='/sps/lsst/users/rcecile/Fig/ps'+root+'.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=6.6,FONT_SIZE=5
   endif else window,0


!p.multi=0
;plot,[0.01,0.1],[1,1],/xs,/ys,xra=[0.02,0.19],/nodata,yra=[.2,150],xtit='wavenumber k [Mpc^-1]',ytit='Power spectrum x k^2 [(Mpc)^-5]',xma=[4.5,1],yma=[3.5,.5],ytickn=replicate(" " ,10),/yl,/xl
plot,[0.01,0.1],[1,1],/xs,/ys,xra=[0.02,0.19],/nodata,yra=[30,80000],xtit='wavenumber k [Mpc^-1]',ytit='Power spectrum [(Mpc)^-3]',xma=[4.5,1],yma=[3.5,.5],ytickn=replicate(" " ,10),/yl,/xl

for iz = 0,nz-1 do begin
   t  = dir + 'simu_ps'+nx[iz]+'_z'+namezmean[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   joli=indgen(60)
   joli = joli^1.5
   xt=(ps.(0))[joli]
   pt=(ps.(1))[joli]
;   oplot,xt,pt*xt*xt,li=2
;    oplot,xt,pt,li=llin[iz],psym=lpsym[iz]
    oplot,xt,pt,th=lth[iz]

   for ir = 0, nerr-1 do begin
      suff = nameerr[ir]
      fobs  = dir + 'PS'+root+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
      print,fobs
      
      pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_4col)     
      p2fit = pobs.(1) -  pobs.(2)

;      oplot,pobs.(0),(p2fit)*pobs.(0)*pobs.(0),col=lcol[ir],li=llin[iz],psym=lpsym[iz]
;      oplot,pobs.(0),(p2fit),col=lcol[ir],li=llin[iz],psym=lpsym[iz]
      oplot,pobs.(0),(p2fit),col=lcol[ir],th=lth[iz]

;      errplot,pobs.(0),(p2fit-pobs.(3))*pobs.(0)*pobs.(0),(p2fit+pobs.(3))*pobs.(0)*pobs.(0),col=lcol[ir]
;      errplot,pobs.(0),(p2fit-pobs.(3)),(p2fit+pobs.(3)),col=lcol[ir]
   endfor
endfor
legend,['fiducial model',[lerr]],li=0,col=[0,[lcol]],th=3,box=1,/fill,/left,/bottom,charsize=1.25
;legend,'z='+lnamez,li=[0,0,2],psym=[-3,4,-3],box=1,col=lcol[4],/fill,/right,/top,charsize=1.5,th=3
legend,'z='+lnamez,li=[0,0,0],th=[10,6,2],box=1,col=lcol[4],/fill,/right,/top,charsize=1.5

if (saveplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

stop
!p.multi=[0,3,5]
window,2+offset,xs=1200,ys=900
 xbao = 2*!pi/150

for ir = 0, nerr-1 do begin
   for iz = 0,nz-1 do begin
      suff = nameerr[ir]
      fobs  = dir + 'PS'+root+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
      
      pobs = read_ascii(fobs, template =  TEMP_POW_SPEC_4col)     
      p2fit = pobs.(1) -  pobs.(2)

      if (iz eq 0) then mymin = 0.04 else if (iz eq 1) then mymin = 0.03 else mymin=0.02
      if (iz eq 0 or ir eq 0) then mymax = 0.2 else if (iz eq 1 and ir gt 0) then mymax = 0.15 else mymax = 0.1
      plot,pobs.(0),(p2fit)*pobs.(0)*pobs.(0),col=lcol[ir],/xs,/ys,/xl,xra=[mymin,mymax],xmar=[0,0],yma=[2,1],/yl
      oplot,[xbao,xbao],[.1,50000],th=20,col=230
      oplot,pobs.(0),(p2fit)*pobs.(0)*pobs.(0),col=lcol[ir]

  ;    errplot,pobs.(0),(p2fit-pobs.(3))*pobs.(0)*pobs.(0),(p2fit+pobs.(3))*pobs.(0)*pobs.(0),col=lcol[ir]
   endfor
endfor

;write_jpeg,'/sps/lsst/users/rcecile/BAO_InOut/ps_check.jpg' ,tvrd(true=3),true=3

END
