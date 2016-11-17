PRO cell.pro

g1= randomu(seed,640,640,640,poisson=1)
g10 = rebin(g1,64,64,64)
g8 = rebin(g1,80,80,80)

restore,'/sps/lsst/dev/rcecile/BAOProgs/DirectSimPerso/PipeScripts/temp_ascii.sav'        ; contient dir
t  =  '/sps/lsst/data/rcecile/TJP_BAO_PS/PS_G2_2fitErrk_640_z0.9_wngal.txt'
pcorr = read_ascii(t, template =  TEMP_2FIT)     
xs = double(pcorr.(0))

ps10 = fft(g10,dim=3,/double)
ps10=abs(ps10)
ps10 = ps10*ps10
ps1 = fft(g10,dim=3,/double)
ps1=abs(ps1)
ps1 = ps1*ps1
plot,ps10,/xs,/ys,psym=1,/yl,yra=[1e-15,1e-1]
oplot,ps1,psym=1,col=123

ps1 =  fft_powerspectrum(g1, xs[2]-xs[1],FREQ=xs)
END
