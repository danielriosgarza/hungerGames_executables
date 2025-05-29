# -*- coding: utf-8 -*-
"""
Created on Thu May 29 11:34:11 2025

@author: drgarza
"""

def getTransitionFunction(statement):
    
    def tf(metObj):
        if statement == '""':
            return 1.0
        return eval(statement)
    return tf
    


def update_Bh_Params(parameter, new_value, reactor):
    
    #1
    if parameter == 'xa_mumax':
        reactor.microbiome.subpopD['xa'].mumax = new_value
        return reactor
        

    #2
    elif parameter == 'xb_mumax':
        reactor.microbiome.subpopD['xb'].mumax = new_value
        return reactor
    
    #3
    elif parameter == 'xa_pHopt':
        reactor.microbiome.subpopD['xa'].pHopt = new_value
        return reactor
    
    #4
    elif parameter == 'xb_pHopt':
        reactor.microbiome.subpopD['xb'].pHopt = new_value
        return reactor
    
    #5
    elif parameter == 'xa_pHalpha':
        reactor.microbiome.subpopD['xa'].pHalpha = new_value
        return reactor
    
    #6
    elif parameter == 'xb_pHalpha':
        reactor.microbiome.subpopD['xb'].pHalpha = new_value
        return reactor
    
    #7
    elif parameter == 'xa_k_s1':
        reactor.microbiome.subpopD['xa'].feedingTerms[0].monodKs[reactor.microbiome.subpopD['xa'].feedingTerms[0].metIDs.index('trehalose')] = new_value
        return reactor
    
    #8
    elif parameter == 'xa_k_s2':
        reactor.microbiome.subpopD['xa'].feedingTerms[1].monodKs[reactor.microbiome.subpopD['xa'].feedingTerms[1].metIDs.index('pyruvate')] = new_value
        return reactor
    
    #9
    elif parameter == 'xb_k_s3':
        reactor.microbiome.subpopD['xb'].feedingTerms[0].monodKs[reactor.microbiome.subpopD['xb'].feedingTerms[0].metIDs.index('glucose')] = new_value
        return reactor
    
    #10
    elif parameter == 'xb_k_s4':
        reactor.microbiome.subpopD['xb'].feedingTerms[0].monodKs[reactor.microbiome.subpopD['xb'].feedingTerms[0].metIDs.index('glutamate')] = new_value
        return reactor
    
    #11
    elif parameter == 'xb_k_s2':
        reactor.microbiome.subpopD['xb'].feedingTerms[1].monodKs[reactor.microbiome.subpopD['xb'].feedingTerms[1].metIDs.index('pyruvate')] = new_value
        return reactor
    
    #12
    elif parameter == 'xa_g_s1':
        reactor.microbiome.subpopD['xa'].feedingTerms[0].yields[reactor.microbiome.subpopD['xa'].feedingTerms[0].metIDs.index('trehalose')] = new_value
        return reactor
    
    #13
    elif parameter == 'xa_g_s2':
        reactor.microbiome.subpopD['xa'].feedingTerms[1].yields[reactor.microbiome.subpopD['xa'].feedingTerms[1].metIDs.index('pyruvate')] = new_value
        return reactor
    
    #14
    elif parameter == 'xa_g_s6_s1':
        reactor.microbiome.subpopD['xa'].feedingTerms[0].yields[reactor.microbiome.subpopD['xa'].feedingTerms[0].metIDs.index('acetate')] = new_value
        return reactor
    
    #15
    elif parameter == 'xa_g_s5_s1':
        reactor.microbiome.subpopD['xa'].feedingTerms[0].yields[reactor.microbiome.subpopD['xa'].feedingTerms[0].metIDs.index('lactate')] = new_value
        return reactor
    
    #16
    elif parameter == 'xa_g_s6_s2':
        reactor.microbiome.subpopD['xa'].feedingTerms[1].yields[reactor.microbiome.subpopD['xa'].feedingTerms[1].metIDs.index('acetate')] = new_value
        return reactor
    
    #17
    elif parameter == 'xb_g_s3':
        reactor.microbiome.subpopD['xb'].feedingTerms[0].yields[reactor.microbiome.subpopD['xb'].feedingTerms[0].metIDs.index('glucose')] = new_value
        return reactor
    
    #18
    elif parameter == 'xb_g_s4':
        reactor.microbiome.subpopD['xb'].feedingTerms[0].yields[reactor.microbiome.subpopD['xb'].feedingTerms[0].metIDs.index('glutamate')] = new_value
        return reactor
    
    #19
    elif parameter == 'xb_g_s2':
        reactor.microbiome.subpopD['xb'].feedingTerms[1].yields[reactor.microbiome.subpopD['xb'].feedingTerms[1].metIDs.index('pyruvate')] = new_value
        return reactor
    
    #20
    elif parameter == 'xb_g_s6_s3_s4':
        reactor.microbiome.subpopD['xb'].feedingTerms[0].yields[reactor.microbiome.subpopD['xb'].feedingTerms[0].metIDs.index('acetate')] = new_value
        return reactor
    
    #21
    elif parameter == 'xb_g_s6_s2':
        reactor.microbiome.subpopD['xb'].feedingTerms[1].yields[reactor.microbiome.subpopD['xb'].feedingTerms[1].metIDs.index('acetate')] = new_value
        return reactor
    
    
    #22
    elif parameter == 'z1_r':
        reactor.microbiome.bacteria['bh'].connections['xa'][0] = ('xb', getTransitionFunction(f"(4.572019894399926e-08**1/(4.572019894399926e-08 + metObj.metD['trehalose'].concentration**1.0))"), new_value)
        return reactor
    
    #23
    elif parameter == 'z2_r':
        reactor.microbiome.bacteria['bh'].connections['xb'][0] = ('xa', getTransitionFunction(f"(metObj.metD['trehalose'].concentration**50.0)/(metObj.metD['trehalose'].concentration**50.0 + 0.21**50.0)"), new_value)
        return reactor
    
    #24
    elif parameter == 'z3_r':
        reactor.microbiome.bacteria['bh'].connections['xa'][1] = ('xc', getTransitionFunction('""'), new_value)
        return reactor
    
    #25
    elif parameter == 'z4_r':
        reactor.microbiome.bacteria['bh'].connections['xb'][1] = ('xc', getTransitionFunction(f"(0.0001293346**2.9055456987/(0.0001293346**2.9055456987 + ((metObj.metD['glucose'].concentration + metObj.metD['glutamate'].concentration - ((metObj.metD['glucose'].concentration - metObj.metD['glutamate'].concentration)**2)**0.5)/2)**2.905545698700455))"), new_value)
        return reactor
    
    #26
    elif parameter == 'z5_r':
        reactor.microbiome.bacteria['bh'].connections['xc'][0] = ('xd', getTransitionFunction('""'), new_value)
        return reactor
    
    #27
    elif parameter == 'z1_l_s1':
        current = reactor.microbiome.bacteria['bh'].connections['xa'][0][:]
        reactor.microbiome.bacteria['bh'].connections['xa'][0] = (current[0], getTransitionFunction(f"({new_value}**1/({new_value} + metObj.metD['trehalose'].concentration**1.0))"), current[2])
        return reactor
    
    #28
    elif parameter == 'z2_l_s1':
        current = reactor.microbiome.bacteria['bh'].connections['xb'][0][:]
        reactor.microbiome.bacteria['bh'].connections['xb'][0] = (current[0], getTransitionFunction(f"(metObj.metD['trehalose'].concentration**50.0)/(metObj.metD['trehalose'].concentration**50.0 + {new_value}**50.0)"), current[2])
        return reactor
    
    #29
    elif parameter == 'z4_l_s3_s4':
        current = reactor.microbiome.bacteria['bh'].connections['xb'][1][:]
        reactor.microbiome.bacteria['bh'].connections['xb'][1] = (current[0], getTransitionFunction(f"({new_value}**2.9055456987/({new_value}**2.9055456987 + ((metObj.metD['glucose'].concentration + metObj.metD['glutamate'].concentration - ((metObj.metD['glucose'].concentration - metObj.metD['glutamate'].concentration)**2)**0.5)/2)**2.905545698700455))"), current[2])
        return reactor
    
    #30
    elif parameter == 'z1_h_s1':
        current = reactor.microbiome.bacteria['bh'].connections['xa'][0][:]
        reactor.microbiome.bacteria['bh'].connections['xa'][0] = (current[0], getTransitionFunction(f"(4.572019894399926e-08**{new_value}/(4.572019894399926e-08 + metObj.metD['trehalose'].concentration**{new_value}))"), current[2])
        return reactor
    
    #31
    elif parameter == 'z2_h_s1':
        current = reactor.microbiome.bacteria['bh'].connections['xb'][0][:]
        reactor.microbiome.bacteria['bh'].connections['xb'][0] = (current[0], getTransitionFunction(f"(metObj.metD['trehalose'].concentration**{new_value})/(metObj.metD['trehalose'].concentration**{new_value} + 0.21**{new_value})"), current[2])
        return reactor
    
    #32
    elif parameter == 'z4_h_s3_s4':
        current = reactor.microbiome.bacteria['bh'].connections['xb'][0][:]
        reactor.microbiome.bacteria['bh'].connections['xb'][0] = (current[0], getTransitionFunction(f"(metObj.metD['trehalose'].concentration**{new_value})/(metObj.metD['trehalose'].concentration**{new_value} + 0.21**{new_value})"), current[2])
        return reactor
   
    