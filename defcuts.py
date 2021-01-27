#from ROOT import * 

def defcuts():
  cuts = {
    "PrimaryPion":        "reco_beam_true_byHits_process == \"primary\" && reco_beam_true_byHits_matched && reco_beam_true_byHits_origin == 4 && reco_beam_true_byHits_PDG == 211",
    "PrimaryMuon":        "reco_beam_true_byHits_process == \"primary\" && reco_beam_true_byHits_matched && reco_beam_true_byHits_origin == 4 && reco_beam_true_byHits_PDG == -13",
    "PrimaryProton":      "reco_beam_true_byHits_process == \"primary\" && reco_beam_true_byHits_matched && reco_beam_true_byHits_origin == 4 && reco_beam_true_byHits_PDG == 2212",
    "PrimaryElectron":    "reco_beam_true_byHits_process == \"primary\" && reco_beam_true_byHits_matched && reco_beam_true_byHits_origin == 4 && ( reco_beam_true_byHits_PDG == -11 || reco_beam_true_byHits_PDG == 11 )",

    "Cosmic":             "reco_beam_true_byHits_origin == 2",

    #"PrimaryBeamNotTrig": "reco_beam_true_byHits_process == \"primary\" && !reco_beam_true_byHits_matched && reco_beam_true_byHits_origin == 4",
    "PrimaryBeamNotTrig": "reco_beam_true_byHits_process == \"primaryBackground\"",

    "UpstreamInt":        "(reco_beam_true_byHits_process == \"neutronInelastic\" ||  reco_beam_true_byHits_process == \"protonInelastic\" || reco_beam_true_byHits_process == \"pi-Inelastic\" || reco_beam_true_byHits_process == \"hadElastic\" || reco_beam_true_byHits_process == \"pi+Inelastic\") && reco_beam_true_byHits_origin == 4",
    "NeutronInel":        "reco_beam_true_byHits_process == \"neutronInelastic\"",
    "ProtonInel":         "reco_beam_true_byHits_process == \"protonInelastic\"",

    "Decay":      "reco_beam_true_byHits_process == \"Decay\" && reco_beam_true_byHits_origin == 4",
    "hadElastic":         "reco_beam_true_byHits_process == \"hadElastic\"",
    "PiMinusElastic":     "reco_beam_true_byHits_process == \"pi-Inelastic\"",
    "Other":              "reco_beam_true_byHits_origin == -1"

  }

  cuts["PrimaryPionDecay"]     = cuts["PrimaryPion"] + " && true_beam_endProcess == \"Decay\""

  cuts["PrimaryPionFastScint"] = cuts["PrimaryPion"] + " && true_beam_endProcess == \"FastScintillation\""
  cuts["PrimaryPionInteract"]  = cuts["PrimaryPion"] + " && true_beam_endProcess == \"pi+Inelastic\"" 

  cuts["UpstreamIntToPiPlus"]  = cuts["UpstreamInt"] + " && reco_beam_true_byHits_PDG == 211"
  cuts["UpstreamIntToKaon"]    = cuts["UpstreamInt"] + " && abs(reco_beam_true_byHits_PDG) == 321"
  cuts["UpstreamIntToPiMinus"] = cuts["UpstreamInt"] + " && reco_beam_true_byHits_PDG == -211"
  cuts["UpstreamIntToProton"]  = cuts["UpstreamInt"] + " && reco_beam_true_byHits_PDG == 2212"
  cuts["UpstreamIntToNuc"]     = cuts["UpstreamInt"] + " && reco_beam_true_byHits_PDG  > 2212"
  cuts["UpstreamIntToGamma"]   = cuts["UpstreamInt"] + " && reco_beam_true_byHits_PDG == 22"
  cuts["UpstreamIntToPion"]    = cuts["UpstreamInt"] + " && abs(reco_beam_true_byHits_PDG) == 211"


  return cuts

def testcuts(e):
  if( e.reco_beam_true_byHits_process == "primary" and e.reco_beam_true_byHits_matched and e.reco_beam_true_byHits_origin == 4 and e.reco_beam_true_byHits_PDG == 211 ):
    return "PrimaryPion"
  elif( e.reco_beam_true_byHits_process == "primary" and e.reco_beam_true_byHits_matched and e.reco_beam_true_byHits_origin == 4 and e.reco_beam_true_byHits_PDG == -13 ):
    return "PrimaryMuon"
  elif( e.reco_beam_true_byHits_process == "primary" and e.reco_beam_true_byHits_matched and e.reco_beam_true_byHits_origin ==4 and e.reco_beam_true_byHits_PDG == 2212 ):
    return "PrimaryProton"
  elif( e.reco_beam_true_byHits_process == "primary" and e.reco_beam_true_byHits_matched and e.reco_beam_true_byHits_origin ==4 and abs(e.reco_beam_true_byHits_PDG) == 11 ):
    return "PrimaryElectron"
  elif( e.reco_beam_true_byHits_origin == 2 ):
    return "Cosmic"
  #elif( e.reco_beam_true_byHits_process == "primary" and not e.reco_beam_true_byHits_matched  and e.reco_beam_true_byHits_origin == 4 ):
  elif( e.reco_beam_true_byHits_process == "primaryBackground" ):
    return "PrimaryBeamNotTrig"



  elif( ( e.reco_beam_true_byHits_process == "neutronInelastic" or e.reco_beam_true_byHits_process == "protonInelastic" 
  or    e.reco_beam_true_byHits_process == "pi-Inelastic" or e.reco_beam_true_byHits_process == "hadElastic" 
  or    e.reco_beam_true_byHits_process == "pi+Inelastic" ) and e.reco_beam_true_byHits_origin == 4  and e.reco_beam_true_byHits_PDG == 211 ):
    return "UpstreamIntToPiPlus"
  elif( ( e.reco_beam_true_byHits_process == "neutronInelastic" or e.reco_beam_true_byHits_process == "protonInelastic" 
  or    e.reco_beam_true_byHits_process == "pi-Inelastic" or e.reco_beam_true_byHits_process == "hadElastic" 
  or    e.reco_beam_true_byHits_process == "pi+Inelastic" ) and e.reco_beam_true_byHits_origin == 4  and e.reco_beam_true_byHits_PDG == 321 ):
    return "UpstreamIntToKaon"
  elif( ( e.reco_beam_true_byHits_process == "neutronInelastic" or e.reco_beam_true_byHits_process == "protonInelastic" 
  or    e.reco_beam_true_byHits_process == "pi-Inelastic" or e.reco_beam_true_byHits_process == "hadElastic" 
  or    e.reco_beam_true_byHits_process == "pi+Inelastic" ) and e.reco_beam_true_byHits_origin == 4  and e.reco_beam_true_byHits_PDG == -211 ):
    return "UpstreamIntToPiMinus" 
  elif( ( e.reco_beam_true_byHits_process == "neutronInelastic" or e.reco_beam_true_byHits_process == "protonInelastic" 
  or    e.reco_beam_true_byHits_process == "pi-Inelastic" or e.reco_beam_true_byHits_process == "hadElastic" 
  or    e.reco_beam_true_byHits_process == "pi+Inelastic" ) and e.reco_beam_true_byHits_origin == 4  and e.reco_beam_true_byHits_PDG == 2212 ):
    return "UpstreamIntToProton"
  elif( ( e.reco_beam_true_byHits_process == "neutronInelastic" or e.reco_beam_true_byHits_process == "protonInelastic" 
  or    e.reco_beam_true_byHits_process == "pi-Inelastic" or e.reco_beam_true_byHits_process == "hadElastic" 
  or    e.reco_beam_true_byHits_process == "pi+Inelastic" ) and e.reco_beam_true_byHits_origin == 4  and e.reco_beam_true_byHits_PDG >  2212 ):
    return "UpstreamIntToNuc" 

#  elif( e.reco_beam_true_byHits_process == "pi+Inelastic"  and e.reco_beam_true_byHits_PDG == 22 ):
#    return "PionInelToGamma"
#  elif( ( e.reco_beam_true_byHits_process == "neutronInelastic" or e.reco_beam_true_byHits_process == "protonInelastic" 
#  or    e.reco_beam_true_byHits_process == "pi-Inelastic" or e.reco_beam_true_byHits_process == "hadElastic" ) and e.reco_beam_true_byHits_origin == 4 ):
#    return "UpstreamInt"    
#  elif( e.reco_beam_true_byHits_process == "neutronInelastic" ):
#    return "NeutronInel"
#  elif( e.reco_beam_true_byHits_process == "protonInelastic" ):
#    return"ProtonInel"
  elif( e.reco_beam_true_byHits_process == "Decay" ):
    return "Decay"
  elif( e.reco_beam_true_byHits_origin == -1 ):
    return "Other"

  return "bad"
  
def testcuts_FS(e):

  if not ( e.reco_beam_type == 13 and (e.true_beam_PDG == -13 or e.true_beam_PDG == 211 or e.true_beam_PDG == 2212 or abs(e.true_beam_PDG) == 11) ):
    return "bad"

  if(   e.reco_beam_true_byHits_process == "primary" and e.reco_beam_true_byHits_matched 
  and   e.reco_beam_true_byHits_origin == 4 and e.reco_beam_true_byHits_PDG == 211 
  and   e.true_beam_endProcess == "Decay"):
    return "PrimaryPionDecay"

  elif( e.reco_beam_true_byHits_process == "primary" and e.reco_beam_true_byHits_matched 
  and   e.reco_beam_true_byHits_origin == 4 and e.reco_beam_true_byHits_PDG == 211 
  and   e.true_beam_endProcess == "FastScintillation"):
    return "PrimaryPionFastScint"

  elif( e.reco_beam_true_byHits_process == "primary" and e.reco_beam_true_byHits_matched 
  and   e.reco_beam_true_byHits_origin == 4 and e.reco_beam_true_byHits_PDG == 211 
  and   e.true_beam_endProcess == "pi+Inelastic"):
    return "PrimaryPionInteract"

  elif( e.reco_beam_true_byHits_process == "primary" and e.reco_beam_true_byHits_matched and e.reco_beam_true_byHits_origin == 4 and e.reco_beam_true_byHits_PDG == -13 ):
    return "PrimaryMuon"
  elif( e.reco_beam_true_byHits_process == "primary" and e.reco_beam_true_byHits_matched and e.reco_beam_true_byHits_origin ==4 and e.reco_beam_true_byHits_PDG == 2212 ):
    return "PrimaryProton"
  elif( e.reco_beam_true_byHits_process == "primary" and e.reco_beam_true_byHits_matched and e.reco_beam_true_byHits_origin ==4 and abs(e.reco_beam_true_byHits_PDG) == 11 ):
    return "PrimaryElectron"
  elif( e.reco_beam_true_byHits_origin == 2 ):
    return "Cosmic"
  #elif( e.reco_beam_true_byHits_process == "primary" and not e.reco_beam_true_byHits_matched  and e.reco_beam_true_byHits_origin == 4 ):
  elif( e.reco_beam_true_byHits_process == "primaryBackground" ):
    return "PrimaryBeamNotTrig"

  elif( ( e.reco_beam_true_byHits_process == "neutronInelastic" or e.reco_beam_true_byHits_process == "protonInelastic" 
  or    e.reco_beam_true_byHits_process == "pi-Inelastic" or e.reco_beam_true_byHits_process == "hadElastic" 
  or    e.reco_beam_true_byHits_process == "pi+Inelastic" ) and e.reco_beam_true_byHits_origin == 4  and e.reco_beam_true_byHits_PDG == 211 ):
    return "UpstreamIntToPiPlus"
  elif( ( e.reco_beam_true_byHits_process == "neutronInelastic" or e.reco_beam_true_byHits_process == "protonInelastic" 
  or    e.reco_beam_true_byHits_process == "pi-Inelastic" or e.reco_beam_true_byHits_process == "hadElastic" 
  or    e.reco_beam_true_byHits_process == "pi+Inelastic" ) and e.reco_beam_true_byHits_origin == 4  and e.reco_beam_true_byHits_PDG == 321 ):
    return "UpstreamIntToKaon"
  elif( ( e.reco_beam_true_byHits_process == "neutronInelastic" or e.reco_beam_true_byHits_process == "protonInelastic" 
  or    e.reco_beam_true_byHits_process == "pi-Inelastic" or e.reco_beam_true_byHits_process == "hadElastic" 
  or    e.reco_beam_true_byHits_process == "pi+Inelastic" ) and e.reco_beam_true_byHits_origin == 4  and e.reco_beam_true_byHits_PDG == -211 ):
    return "UpstreamIntToPiMinus" 
  elif( ( e.reco_beam_true_byHits_process == "neutronInelastic" or e.reco_beam_true_byHits_process == "protonInelastic" 
  or    e.reco_beam_true_byHits_process == "pi-Inelastic" or e.reco_beam_true_byHits_process == "hadElastic" 
  or    e.reco_beam_true_byHits_process == "pi+Inelastic" ) and e.reco_beam_true_byHits_origin == 4  and e.reco_beam_true_byHits_PDG == 2212 ):
    return "UpstreamIntToProton"
  elif( ( e.reco_beam_true_byHits_process == "neutronInelastic" or e.reco_beam_true_byHits_process == "protonInelastic" 
  or    e.reco_beam_true_byHits_process == "pi-Inelastic" or e.reco_beam_true_byHits_process == "hadElastic" 
  or    e.reco_beam_true_byHits_process == "pi+Inelastic" ) and e.reco_beam_true_byHits_origin == 4  and e.reco_beam_true_byHits_PDG >  2212 ):
    return "UpstreamIntToNuc"  

#  elif( e.reco_beam_true_byHits_process == "pi+Inelastic"  and e.reco_beam_true_byHits_PDG == 22 ):
#    return "PionInelToGamma"
#  elif( e.reco_beam_true_byHits_process == "neutronInelastic" or e.reco_beam_true_byHits_process == "protonInelastic" 
#  or    e.reco_beam_true_byHits_process == "pi-Inelastic" or e.reco_beam_true_byHits_process == "hadElastic" ):
#    return "UpstreamInt"    
  #elif( e.reco_beam_true_byHits_process == "neutronInelastic" ):
  #  return "NeutronInel"
  #elif( e.reco_beam_true_byHits_process == "protonInelastic" ):
  #  return"ProtonInel"
  elif( e.reco_beam_true_byHits_process == "Decay" ):
    return "Decay"
  elif( e.reco_beam_true_byHits_origin == -1 ):
    return "Other"

  return "bad"

def ang_cut_str(coslow=.93):
  cut = " && ( true_beam_startDirX*reco_beam_trackDirX + true_beam_startDirY*reco_beam_trackDirY + true_beam_startDirZ*reco_beam_trackDirZ > " + str(coslow) + " ) "
  return cut

def pos_cut_str(  xlow=-3., xhigh=7., ylow=-8., yhigh=7., zlow=27.5, zhigh=32.5):
  cut =  " && ( reco_beam_startX - (true_beam_startX + -1.*true_beam_startZ*(true_beam_startDirX/true_beam_startDirZ) ) > " + str(xlow) +" ) "
  cut += " && ( reco_beam_startX - (true_beam_startX + -1.*true_beam_startZ*(true_beam_startDirX/true_beam_startDirZ) ) < " + str(xhigh) +" ) "
  cut += " && ( reco_beam_startY - (true_beam_startY + -1.*true_beam_startZ*(true_beam_startDirY/true_beam_startDirZ) ) > " + str(ylow) +" ) "
  cut += " && ( reco_beam_startY - (true_beam_startY + -1.*true_beam_startZ*(true_beam_startDirY/true_beam_startDirZ) ) < " + str(yhigh) +" ) "
  cut += " && (reco_beam_startZ > " + str(zlow) +" ) "
  cut += " && (reco_beam_startZ < " + str(zhigh) +" ) "
  
  return cut

def data_ang_cut_str(coslow=.93):
  cut = " && ( data_BI_dirX*reco_beam_trackDirX + data_BI_dirY*reco_beam_trackDirY + data_BI_dirZ*reco_beam_trackDirZ > " + str(coslow) + " ) "
  return cut

def data_pos_cut_str( xlow=0., xhigh=10., ylow=-5., yhigh=10., zlow=30., zhigh=35.):
  cut =  " && ( (reco_beam_startX - data_BI_X) > " + str(xlow) +" ) "
  cut += " && ( (reco_beam_startX - data_BI_X) < " + str(xhigh) +" ) "
  cut += " && ( (reco_beam_startY - data_BI_Y) > " + str(ylow) +" ) "
  cut += " && ( (reco_beam_startY - data_BI_Y) < " + str(yhigh) +" ) "
  cut += " && (reco_beam_startZ > " + str(zlow) +" ) "
  cut += " && (reco_beam_startZ < " + str(zhigh) +" ) "
  cut += " && @data_BI_nMomenta.size() == 1 && @data_BI_nTracks.size()"
  
  return cut




def ang_pos_test_cut(e, coslow=.93, xlow=-3., xhigh=7., ylow=-8., yhigh=7., zlow=27.5, zhigh=32.5):
  if (e.true_beam_startDirX*e.reco_beam_trackDirX + e.true_beam_startDirY*e.reco_beam_trackDirY + e.true_beam_startDirZ*e.reco_beam_trackDirZ < coslow): return 0

  if ( e.reco_beam_startX - (e.true_beam_startX + -1.*e.true_beam_startZ*(e.true_beam_startDirX/e.true_beam_startDirZ) ) < xlow ): return 0
                            
  if ( e.reco_beam_startX - (e.true_beam_startX + -1.*e.true_beam_startZ*(e.true_beam_startDirX/e.true_beam_startDirZ) ) > xhigh ): return 0
                            
  if ( e.reco_beam_startY - (e.true_beam_startY + -1.*e.true_beam_startZ*(e.true_beam_startDirY/e.true_beam_startDirZ) ) < ylow ): return 0
                            
  if ( e.reco_beam_startY - (e.true_beam_startY + -1.*e.true_beam_startZ*(e.true_beam_startDirY/e.true_beam_startDirZ) ) > yhigh ): return 0

  if( e.reco_beam_startZ < zlow or e.reco_beam_startZ > zhigh ): return 0

  return 1

def data_ang_pos_test_cut(e, coslow=.93, xlow=0., xhigh=10., ylow=-5., yhigh=10., zlow=30., zhigh=35.):
  
  if (e.data_BI_nMomenta != 1 or e.data_BI_nTracks != 1 ): 
    #print e.event, e.data_BI_nMomenta, e.data_BI_nTracks
    return 0

  if ( (e.data_BI_dirX*e.reco_beam_trackDirX + e.data_BI_dirY*e.reco_beam_trackDirY + e.data_BI_dirZ*e.reco_beam_trackDirZ) < coslow): 
    #print e.event, "cos", (e.data_BI_dirX*e.reco_beam_trackDirX + e.data_BI_dirY*e.reco_beam_trackDirY + e.data_BI_dirZ*e.reco_beam_trackDirZ) 
    return 0

  if ( (e.reco_beam_startX - e.data_BI_X ) < xlow ): 
    #print e.event, "X", (e.reco_beam_startX - e.data_BI_X ), e.reco_beam_startX, e.data_BI_X 
    return 0

  if ( (e.reco_beam_startX - e.data_BI_X ) > xhigh ): 
    #print e.event, "X", (e.reco_beam_startX - e.data_BI_X ), e.reco_beam_startX, e.data_BI_X 
    return 0

  if ( (e.reco_beam_startY - e.data_BI_Y ) < ylow ): 
    #print e.event, "Y", (e.reco_beam_startY - e.data_BI_Y ), e.reco_beam_startY, e.data_BI_Y 
    return 0

  if ( (e.reco_beam_startY - e.data_BI_Y ) > yhigh ): 
    #print e.event, "Y", (e.reco_beam_startY - e.data_BI_Y ), e.reco_beam_startY, e.data_BI_Y 
    return 0

  if( e.reco_beam_startZ < zlow or e.reco_beam_startZ > zhigh ): 
    #print e.event, "Z",  e.reco_beam_startZ
    return 0

  return 1

