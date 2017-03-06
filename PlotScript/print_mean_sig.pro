PRO print_mean_sig

suff='cube_'
dir='/sps/lsst/data/rcecile/Planck_BAO_grids/'
D=['0.5','0.9','1.5']

snx=['160','300','225']
nz=[140,125,150]
nD = n_elements(D)
restore,'temp_ascii.sav'

err=['_errP','_errPBDT9','_errPBDT8']
serr=['_photoZ','_photoZ_BDT90','_photoZ_BDT80']

for ir =0,2 do begin
s = read_ascii(dir+"Sigma_"+serr+".txt", template =  TEMP_SEL_FUNC)


for id = 0,nD-1 do begin
   grid=dir+'grid_'+suff+strtrim(snx[id],2)+'_z'+D[id]+'.fits'
   print,grid
   c=mrdfits(grid,6,h)

endfor

endfor






END
