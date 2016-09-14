PRO plot_osc,iz,doplot

!p.charsize=2
!p.thick=2
!p.symsize=4
h= 0.679
h3 = h*h*h

loadct,39
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/TJP_BAO_PS/"
dirn="/sps/lsst/data/rcecile/TJP_noBAO_PS/"

namez = ['0.7','1.4']
nx= ['_450', '_875']
namez = ['0.9','1.3','1.8']
nx= ['_640', '_900', '_1000']

nz = n_elements(namez)
suffn = '_err0.03'

suff = suffn
lpsym=[0,-5]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=0

ltext ='z = '+ namez
lcol  = [80, 150, 200]
lline= 0

t  = dir + 'PS_G2'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt'
print,t
p = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
xs = p.(0)
obs =  p.(1)-p.(6) 

t  = dir + 'simu_ps'+nx[iz]+'_z'+namez[iz]+'_ntpk.txt'
ps = read_ascii(t, template =  TEMP_POW_SPEC_TH)
xt=ps.(0)
osc = ps.(1)
th = ps.(2)
ok = intarr(n_elements(xs))
for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))

!p.multi=0

pref=dblarr(n_elements(p.(0)))
nsim=0
t  = dir + 'PS_G2_2fitErr'+nx[iz]+'_z'+namez[iz]+suff+'_wngal.txt' 
print,t
p = read_ascii(t, template =  TEMP_POW_SPEC_TXT) 
pref = (p.(1)-p.(6))

legend,['"observations"','10 realizations wo BAOs'],col=[0,lcol[iz]+2*5],li=[0,0],psym=[1,0],box=1,/fill,/right,/top,charsize=1.5



plot,xt/h, (osc/th),/xs,/ys,xra=[0.04,0.15],yra=[0,2],xtit='wavenumber/h [Mpc^-1]',ytit='Power spectrum',tit='z = '+namez[iz]
oplot,xs/h,(obs/pref),psym=-4,symsi=2,col=120,th=3

legend,['theoretical case','corrected "observations"/theoretical no BAO case','"observations"/<10 realizations wo BAOs>'],col=[0,200,120],li=0,psym=[0,5,4],box=1,/fill,/right,/bottom,charsize=1.5

;stop


end
