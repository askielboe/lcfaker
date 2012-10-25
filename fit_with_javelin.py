# # Generate lightCurves
# from FakeLight import LightContAndLine
# 
# lc = LightContAndLine(1000, 1.0, 5.0, 100.0, 0.5)
# lc.bin(100)
# lc.observe(5, 50)
# lc.plot('o')
# lc.saveToTxt()

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

