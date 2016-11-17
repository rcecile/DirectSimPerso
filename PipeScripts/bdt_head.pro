PRO bdt_head

href=headfits("/sps/lsst/data/rcecile/TJP_BAO/cat_G2_Slice60.fits",ext=1)

dir='/sps/lsst/PhotozBAO/ricol/SIMU50deg/bao/FITS/'
dir='/sps/lsst/data/rcecile/TJP_BAO/'

list=['H0','OMEGAM0','OMEGAB0','OMEGAR0','OMEGAT0','OMEGADE0','OMEGADK','DE_W0','DE_WA','SIGMA8','N_S','CELL','SEGMSIZE']
listdel=['ZS','ZP','THETA','PHI','TYPE','PODDS','BDT']
for f=0,69 do begin
   name=dir+"file"+strtrim(f,2)+".fits"
   hjs=headfits(name,ext=1)
   print,name,' ok'
   for i=0,n_elements(listdel)-1 do sxdelpar,hjs,listdel[i]
   print,name,' del done'
   for i=0,n_elements(list)-1 do sxaddpar,hjs,list[i], sxpar(href,list[i])
   modfits, name,0,hjs,EXTEN_NO =1
   print,name,' add done'
endfor

end
