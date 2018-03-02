import pyneb as pn

def getTemDen_helper(inQ, outQ, atom, lev_i1, lev_j1, lev_i2, lev_j2,
         wave1, wave2, maxError, method, log, start_x, end_x, to_eval, nCut, maxIter,NLevels=None):
    ''' multiprocessing helper for getTemDen'''
    thisAtom = pn.Atom(atom=atom,NLevels=NLevels)
    for inputs in iter(inQ.get, 'STOP'):
        jobid = inputs[0]
        int_ratio = inputs[1]
        tem = inputs[2]
        den = inputs[3]
        res = thisAtom._getTemDen_1(int_ratio, tem=tem, den=den, lev_i1=lev_i1, lev_j1=lev_j1, lev_i2=lev_i2, lev_j2=lev_j2,
                  wave1=wave1, wave2=wave2, maxError=maxError, method=method, log=log, start_x=start_x, 
                  end_x=end_x, to_eval=to_eval, nCut=nCut, maxIter=maxIter)
        outQ.put((jobid, res))
    return
