#!/usr/bin/python26

from getdist import plots, MCSamples, loadMCSamples
import getdist

ch_path=[]
while 1:
  ch_path.append(raw_input('Path to chains or q to quit: ').strip())
  if ch_path[-1]=='q':
    ch_path=ch_path[:-1]
    break

ch=[]
for p in ch_path:
  ch.append(loadMCSamples(p))

pars=ch[0].getParamNames()
if len(ch)>1:
  print 'Select parameters based on first chain entry; will look for matching parameter names in subsequent entries'
print 'Parameters are: '
print pars
plotpars=[]
while 1:
  p=raw_input('Parameter to plot or q to quit: ').strip()
  if p=='q': break
  plotpars.append(pars.parWithName(p).label)

plotpar_names=['plot'+str(i) for i in range(len(plotpars))]

do_markers=raw_input('Enter true values (y/n): ').strip()
markers={}
if 'y' in do_markers.lower():
  print 'Enter true value or non-numerical value to skip'
  for p in plotpars:
    m=raw_input('True value for '+p+': ')
    try:
      markers[p]=float(m)
    except:
      pass

for i, samps in enumerate(ch):
  parsi=samps.getParamNames()
  for p_label in plotpars:
    found=False
    for j, pj in enumerate(parsi.names):
      if pj.label==p_label:
        found=True
        p_name=plotpar_names[plotpars.index(p_label)]
        if p_label[0]=='M' and 'M_{\odot}' in p_label:
          # Add a derived parameter for mass to get around labelling issues
          M=eval('samps.getParams().'+pj.name)
          M1=samps.ranges.getLower(pj.name)
          M2=samps.ranges.getUpper(pj.name)
          p_label_new=p_label.replace('M_{\odot}', '10^{14} M_{\odot}')
          samps.addDerived(M*1e-14, name=p_name, label=p_label_new)
          samps.updateBaseStatistics()
          # Adding the derived parameter doesn't retain the range info
          if M1 is None:
            M1='N'
          else:
            M1*=1e-14
          if M2 is None:
            M2='N'
          else:
            M2*=1e-14
          samps.setRanges({p_name: [M1, M2]})
          if p_label in markers.keys() and i==0:
            markers[p_label]*=1e-14
        elif p_label[0]=='S' and 'Jy' in p_label:
          # Convert source fluxes to mJy
          S=eval('samps.getParams().'+pj.name)
          S1=samps.ranges.getLower(pj.name)
          S2=samps.ranges.getUpper(pj.name)
          samps.addDerived(S*1e3, name=p_name, label=p_label.replace('Jy', 'mJy'))
          samps.updateBaseStatistics()
          # Adding the derived parameter doesn't retain the range info
          if S1 is None:
            S1='N'
          else:
            S1*=1e3
          if S2 is None:
            S2='N'
          else:
            S2*=1e3
          samps.setRanges({p_name: [S1, S2]})
          if p_label in markers.keys() and i==0:
            markers[p_label]*=1e3
        else:
          pj.name=p_name
    if not found:
      print 'Warning: '+p_label+' not found in '+ch_path[i]
  samps.setParamNames(parsi)

leg=[]
for c in ch:
  leg.append(c.getName())

g=plots.getSubplotPlotter(width_inch=8)
g.settings.axes_fontsize=8
g.settings.alpha_filled_add=0.4
g.triangle_plot(ch, plotpar_names,
    filled_compare=True, 
    legend_labels=leg, 
    legend_loc='upper right')
for ipar, p in enumerate(plotpars):
  if p in markers.keys():
    ax=g.subplots[ipar,ipar]
    ax.axvline(x=markers[p], color='k')
    for jpar, p2 in enumerate(plotpars[ipar+1:]):
      if p2 in markers.keys():
        jplot=ipar+jpar+1
        ax=g.subplots[jplot,ipar]
        ax.plot(markers[p], markers[p2], '*k')
      
g.export(ch_path[0]+'_tri.png')

