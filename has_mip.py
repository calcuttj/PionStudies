def has_mip(e, cut_val):
  
  for dedxs in e.reco_daughter_dEdX:
    these_dedxs = [i for i in dedxs]
    if len( these_dedxs ) < 1: continue 

    avg_dedx = sum( these_dedxs ) / len( these_dedxs )

    if avg_dedx < cut_val: 
      return True

  return False

def has_mip_chi2(e, cut_val):
  
  for chi2 in e.reco_daughter_Chi2_proton:
    if chi2 > 9000.: return True 

    if chi2 > cut_val:
      return True

  return False