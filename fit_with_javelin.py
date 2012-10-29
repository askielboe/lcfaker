# Generate lightCurves
from FakeLight import LightContAndLine

lc = LightContAndLine(2000)

t = lc.lightCurveCont.time
f = lc.lightCurveCont.flux

# from numpy import asarray
# 
# t = asarray(t)
# f = asarray(f)
# 
# nBins = 100
# binSize = (max(t)-min(t))/nBins
# newBins = [i*binSize for i in range(nBins)]
# newBins = asarray(newBins)
# 
# from rebin import rebin
# rebin(t,f,newBins, interp_kind='piecewise_constant')
# 
# from histogram import histogram
# 
# t = 't',asarray(t)
# f = asarray(f)
# h = histogram( 'h', [t], f)

# lc.lightCurveLine.lag_luminosity(lc.lightCurveCont, 100., 1.0, 0.5)
# lc.trim()
# lc.lightCurveLine.smooth(1.5)

lc.reprocess(200, 1.0, 1.5, 100., 1.0, 0.519)

# Should be equivalent to:
lc.lightCurveLine.lag_luminosity(lc.lightCurveCont, 100., 1.0, 0.519)
lc.trim()
lc.rbin(200)

lc.lightCurveLine.smooth(1.5)


###

lc.observeIntervals([250,400,650,800],[45,55,60,65])

lc.plot('o')

lc.saveToTxt()

# Run JAVELIN
from javelin.zylc import get_data
from javelin.lcmodel import Cont_Model, Rmap_Model

c = get_data(["lc_cont.txt"])
cmod = Cont_Model(c)
cmod.do_mcmc()

cy = get_data(["lc_cont.txt","lc_line.txt"], names=["Continuum", "Line"])
cymod = Rmap_Model(cy)
cymod.do_mcmc(conthpd=cmod.hpd)

cymod.show_hist()

cymod.get_hpd()
cymodhpd = cymod.hpd

par_best = cymodhpd[1,:]
print(par_best)

javdata_best =  cymod.do_pred(par_best)
javdata_best.plot(set_pred=True, obs=cy)

