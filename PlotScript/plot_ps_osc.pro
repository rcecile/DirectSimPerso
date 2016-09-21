PRO plot_ps_osc,saveplot

nsuff=['','_err0.03','_errPpodds']

!p.charsize=2

loadct,39
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/TJP_BAO_PS/"

namez = ['0.9','1.3','1.8','1.8']
nx= ['_640', '_900','_1024', '_500']
gsuff = ['', '', ' thin', ' thick']

nz = n_elements(namez)
nk = n_elements(namek)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0

if (saveplot eq 1) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_ps_wosc.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=5.5,FONT_SIZE=4
endif

lpsym=[0,1,0,0]
lcol= [0,80,150,210]
lpsym=[0,4,0,0]
lth=[2,4,4,4]


!p.multi=0
plot,[0.01,0.1],[1e4,1e4],/xs,/ys,xra=[0.0015,0.19],/yl,xtit='wavenumber [Mpc^-1]',ytit='Power spectrum [(Mpc)^-3]',/nodata,yra=[.5,20e4],xma=[9,1],yma=[3,1]

for iz=0,nz-1 do begin
   
   
   t  = dir + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
   ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
   xt=ps.(0)
   pt=ps.(1)
   oplot,xt,pt,li=iz*2,th=4


   for isuff=0,n_elements(nsuff)-1 do begin

      suff = nsuff[isuff]
      
      nsim=0
      t  = dir + 'PS_G2'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
  ;    if (isuff eq 2 and iz eq 2) then t  = dir + 'PS_G2'+nx[iz]+'_z'+namez[iz]+'_k0.06'+suff+'_wngal.txt' 
      print,t
      check = FILE_TEST(t)
      if (check eq 0 ) then continue
      nsim = nsim +  1
      p = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
      xs = p.(0)
      pref = (p.(1)-p.(6))
      oplot,xs,pref,col=lcol[isuff+1],li=iz*2,th=lth[isuff];,psym=lpsym[isuff+1],symsize=2
     
   endfor


endfor   

mytext=['spectroZ','Gauss 0.03','photoZ with podds cut']
mytext = ['theoretical spectrum with BAO',mytext]
legend,mytext,psym=lpsym,col=lcol,box=1,/fill,/left,/bottom,charsize=1.5,th=4
legend,'z = '+ namez+gsuff,line=indgen(nz)*2,box=1,/fill,/right,/top,charsize=1.5,th=4

if (saveplot eq 1) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
   print,'plot saved'
endif


end
