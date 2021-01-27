import ROOT as RT
import sys
from set_style import * 
from add_stack import *
from scale_stack import *
from argparse import ArgumentParser as ap

parser = ap()

parser.add_argument( "-d", type=str, help='Data file' )
parser.add_argument( "-m", type=str, help='MC file' ) 
parser.add_argument( "-s", type=str, help='Name of stack' )
parser.add_argument( "--min", type=float, help='Minimum of user range', default = -1. )
parser.add_argument( "--max", type=float, help='Maximum of user range', default = -1. )
parser.add_argument( "--Xmin", type=float, help='legend positions', default = -1. )
parser.add_argument( "--Xmax", type=float, help='legend positions', default = -1. )
parser.add_argument( "--Ymin", type=float, help='legend positions', default = -1. )
parser.add_argument( "--Ymax", type=float, help='legend positions', default = -1. )
parser.add_argument( "-p", type=str, help="Title", default = "")
parser.add_argument( "-l", type=int, help="Set Log?", default=0 )

parser.add_argument( "-t", type=str, help='Name of output plot file', default='stack_try.pdf' )
parser.add_argument( "-r", type=int, help='Rebin?', default = 0 )
parser.add_argument( "--listNames", type=int, help='List stacks?', default=0)

args = parser.parse_args()


do_range = True
if( args.min == -1. and args.max == -1. ): do_range = False

fMC = RT.TFile(args.m)
if args.listNames:
  for i in fMC.GetListOfKeys(): print(i.GetName())
  exit()

stack_title = args.s 

RT.gROOT.SetBatch(1)
RT.gStyle.SetOptStat(0)

fData = RT.TFile(args.d)
#data_stack = fData.Get(stack_title)

mc_stack = fMC.Get(stack_title)


#mc_hist = mc_stack.GetHists().At(0).Clone()
#for i in range(1,mc_stack.GetNhists()):
#  mc_hist.Add(mc_stack.GetHists().At(i))

data_hist = fData.Get(stack_title).GetHists().At(0).Clone("")
#data_hist = fData.Get(stack_title)


#print data_hist.Integral()
#print add_stack( mc_stack )


#data_hist.Scale( mc_hist.Integral() / data_hist.Integral() )
data_hist.Sumw2()
#data_hist.Scale( add_stack( mc_stack ) / data_hist.Integral() )
mc_stack = scale_stack(mc_stack, data_hist)

if( args.r ):
  data_hist.Rebin(args.r)
  data_hist.Scale( 1. / args.r)
  
  new_hists = []
  for i in range(0, mc_stack.GetNhists()):
    new_hists.append( mc_stack.GetHists().At(i).Clone() )
    new_hists[i].Rebin(args.r) 
    new_hists[i].Scale( 1. / args.r )

  mc_stack = RT.THStack()
  for h in new_hists:
    mc_stack.Add( h )


  
c1 = RT.TCanvas("c1", "c1", 500, 400)
c1.SetTicks()
if args.l: c1.SetLogy()
mc_stack.Draw()
set_style( mc_stack, args.p, "" )
if do_range:
  mc_stack.GetXaxis().SetRangeUser(args.min, args.max)

#print data_hist.GetMaximum(), mc_stack.GetMaximum()
if data_hist.GetMaximum() > mc_stack.GetMaximum():
  mc_stack.SetMaximum(1.1*data_hist.GetMaximum())

mc_stack.Draw("HIST")
data_hist.SetMarkerStyle(20)
data_hist.SetMarkerSize(.5)
data_hist.Draw("same PE")
leg = fMC.Get("leg")
leg.AddEntry(data_hist, "Data", "lp")

if( args.Xmin > -1. and args.Xmax > -1. and args.Ymin > -1. and args.Ymax > -1. ):
  leg.SetX1( args.Xmin )
  leg.SetX2( args.Xmax )
  leg.SetY1( args.Ymin )
  leg.SetY2( args.Ymax )
leg.Draw("same")
c1.SaveAs(args.t)
