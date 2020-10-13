from ROOT import *
from set_style import *
from add_stack import *
from scale_stack import *

from argparse import ArgumentParser as ap

parser = ap()

parser.add_argument( "-d", type=str, help='Data file' )
parser.add_argument( "-m", type=str, help='MC file' ) 
parser.add_argument( "-t", type=str, help='Title of stack/hist' )
parser.add_argument( "--Xmin", type=float, help='legend positions', default = -1. )
parser.add_argument( "--Xmax", type=float, help='legend positions', default = -1. )
parser.add_argument( "--Ymin", type=float, help='legend positions', default = -1. )
parser.add_argument( "--Ymax", type=float, help='legend positions', default = -1. )
parser.add_argument( "-s", type=str, help='Name of output plot file', default='stack_try.pdf' )
parser.add_argument( "-l", type=int, help='Cut level', default = 0 )
parser.add_argument( "--listNames", type=int, help='List contents of MC file and exit', default = 0 )

args = parser.parse_args()

fMC   = TFile(args.m) 
if args.listNames:
  for i in fMC.GetListOfKeys(): print i.GetName()
  exit()


gROOT.SetBatch(1)
gStyle.SetOptStat(0)

fData = TFile(args.d)

c1 = TCanvas("c1", "c1", 500, 400)
c1.SetTicks()

if args.l == 0:
  data_hist = fData.Get(args.t + "hist")
  mc_hist   = fMC.Get(args.t + "stack")

elif args.l == 1:
  data_hist = fData.Get("first_cut_dir/" + args.t + "_ang_pos_cut")
  mc_hist   = fMC.Get("first_cut_dir/" + args.t + "stack_ang_pos")

elif args.l == 2:
  data_hist = fData.Get("second_cut_dir/" + args.t + "_ang_pos_endZ_cut")
  mc_hist   = fMC.Get("second_cut_dir/" + args.t + "stack_ang_pos_len")

elif args.l == 3:
  data_hist = fData.Get("third_cut_dir/" + args.t + "_ang_pos_endZ_chi2_cut")
  mc_hist   = fMC.Get("third_cut_dir/" + args.t + "stack_ang_pos_chi2")

else: exit()

data_hist.Sumw2()

#data_hist.Scale( add_stack(mc_hist) / data_hist.Integral() )
mc_hist = scale_stack(mc_hist, data_hist)

mc_hist.SetMaximum( max([mc_hist.GetMaximum(), data_hist.GetMaximum()]) )

mc_hist.Draw("HIST")
set_style( mc_hist, data_hist.GetXaxis().GetTitle(), "")
data_hist.Draw("same Pe")

leg = fMC.Get("leg")
leg.AddEntry(data_hist, "Data", "lp")

if( args.Xmin > -1. and args.Xmax > -1. and args.Ymin > -1. and args.Ymax > -1. ):
  leg.SetX1( args.Xmin )
  leg.SetX2( args.Xmax )
  leg.SetY1( args.Ymin )
  leg.SetY2( args.Ymax )
leg.Draw("same")
  




'''elif args.l == 1:
  data_hist = fData.Get("first_cut_dir/" + args.t + "_ang_pos_cut")
  mc_hist   = fMC.Get("first_cut_dir/" + args.t + "stack_ang_pos")

  data_hist.Sumw2()

  #data_hist.Scale( add_stack(mc_hist) / data_hist.Integral() )
  mc_hist = scale_stack(mc_hist, data_hist)

  mc_hist.Draw("HIST")
  set_style( mc_hist, data_hist.GetXaxis().GetTitle(), "")
  data_hist.Draw("same Pe")
  
  leg = fMC.Get("leg")
  leg.AddEntry(data_hist, "Data", "lp")

  if( args.Xmin > -1. and args.Xmax > -1. and args.Ymin > -1. and args.Ymax > -1. ):
    leg.SetX1( args.Xmin )
    leg.SetX2( args.Xmax )
    leg.SetY1( args.Ymin )
    leg.SetY2( args.Ymax )
  leg.Draw("same")
  
elif args.l == 2:
  data_hist = fData.Get("second_cut_dir/" + args.t + "_ang_pos_endZ_cut")
  mc_hist   = fMC.Get("second_cut_dir/" + args.t + "stack_ang_pos_len")

  data_hist.Sumw2()


  #data_hist.Scale( add_stack(mc_hist) / data_hist.Integral() )
  mc_hist = scale_stack(mc_hist, data_hist)

  mc_hist.Draw("HIST")
  set_style( mc_hist, data_hist.GetXaxis().GetTitle(), "")
  data_hist.Draw("same Pe")
  
  leg = fMC.Get("leg")
  leg.AddEntry(data_hist, "Data", "lp")

  if( args.Xmin > -1. and args.Xmax > -1. and args.Ymin > -1. and args.Ymax > -1. ):
    leg.SetX1( args.Xmin )
    leg.SetX2( args.Xmax )
    leg.SetY1( args.Ymin )
    leg.SetY2( args.Ymax )
  leg.Draw("same")
  
elif args.l == 3:
  data_hist = fData.Get("third_cut_dir/" + args.t + "_ang_pos_endZ_chi2_cut")
  mc_hist   = fMC.Get("third_cut_dir/" + args.t + "stack_ang_pos_chi2")

  data_hist.Sumw2()


  #data_hist.Scale( add_stack(mc_hist) / data_hist.Integral() )
  mc_hist = scale_stack(mc_hist, data_hist)

  mc_hist.Draw("HIST")
  set_style( mc_hist, data_hist.GetXaxis().GetTitle(), "")
  data_hist.Draw("same Pe")
  
  leg = fMC.Get("leg")
  leg.AddEntry(data_hist, "Data", "lp")

  if( args.Xmin > -1. and args.Xmax > -1. and args.Ymin > -1. and args.Ymax > -1. ):
    leg.SetX1( args.Xmin )
    leg.SetX2( args.Xmax )
    leg.SetY1( args.Ymin )
    leg.SetY2( args.Ymax )
  leg.Draw("same")
 

else: exit() '''

c1.SaveAs(args.s)


