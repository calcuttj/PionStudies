import ROOT as RT
def scale_stack(s, h):
  new_stack = RT.THStack()
  hists = [i.Clone() for i in s.GetHists()]

  stack_total = 0
  for hist in hists:
    stack_total += hist.Integral()

  for hist in hists:
    hist.Scale(h.Integral() / stack_total)
    new_stack.Add(hist)
  
  return new_stack
