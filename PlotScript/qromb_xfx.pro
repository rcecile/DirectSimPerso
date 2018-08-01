
function qromb_xfx, x
  common qromb_x_common, funcname
  return, call_function(funcname, x)*x
end
