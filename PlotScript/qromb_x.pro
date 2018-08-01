function qromb_x, fname, a, b, _extra=extra
  common qromb_x_common, funcname
  funcname = fname
  return, qromb('qromb_xfx', a, b, _extra=extra)
end
