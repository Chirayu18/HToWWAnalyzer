import dask
import dask_awkward as dak
import awkward as ak
def applyTriggerPaths(events):
    trigmask = (events.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ ) | (events.HLT.Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ )
    trig_events = events[trigmask]
    trigLeps = events.TrigObj[(events.TrigObj.id == 11) | (events.TrigObj.id == 13)]
    trigLeps = trigLeps[dak.argsort(trigLeps.pt, ascending = False, axis=1)]
    lead_trig_lep = trigLeps[:,0]
    sublead_trig_lep = trigLeps[:,1]
    pt_req_mask = (lead_trig_lep.pt>25) & (sublead_trig_lep.pt>15)
    return trig_events[pt_req_mask] , lead_trig_lep, sublead_trig_lep
     
