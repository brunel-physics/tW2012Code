from ROOT import *
from random import random

def makeLeptonCuts(cutFlowMode,runType,tree,returnLeptons = False,fillTreeFlag=False,fillTree=0,weight=1.0,btagUD=0,plotList=[],dataset="",plotConf=0,LESSyst = 0):
    btagUp = True if btagUD == 1 else False
    btagDown = True if btagUD == 2 else False
    selectedEles =[]
    selectedMuons = []
    numEles = 0
    #Electron ID!"
#    print "Electron"
    for i in range(tree.numElePF2PAT):
        #Do electron ID here - TODO make cuts a config file?
#        print tree.elePF2PATPT[i]
#        print getLeptonPt(tree,i,True,LESSyst)
        if getLeptonPt(tree,i,True,LESSyst) < 20: continue
        if fabs(tree.elePF2PATEta[i]) > 2.5 : continue
        if fabs(tree.elePF2PATBeamSpotCorrectedTrackD0[i]) > 0.04 : continue
#        if fabs(tree.elePF2PATTrackD0[i]) >= 0.04: continue
        if tree.elePF2PATMissingInnerLayers[i] > 1: continue
        if not tree.elePF2PATPhotonConversionVeto[i] : continue
#        print tree.elePF2PATPT[i] , tree.elePF2PATMVA[i] , tree.elePF2PATComRelIsodBeta[i]
        if tree.elePF2PATMVA[i] < 0.5 : continue
        if tree.elePF2PATComRelIsoRho[i]/getLeptonPt(tree,i,True,LESSyst) > 0.15: continue
        numEles += 1
        selectedEles.append(i)
    #Loose electron ID
    numLooseEles = 0
    for i in range(tree.numElePF2PAT):
        if getLeptonPt(tree,i,True,LESSyst) < 10: continue
        if fabs(tree.elePF2PATEta[i]) > 2.5: continue
        if tree.elePF2PATComRelIsoRho[i] / getLeptonPt(tree,i,True,LESSyst)  > 0.15 : continue
        if tree.elePF2PATMVA[i] < 0.5: continue
        numLooseEles +=1
    #Muon id!
    numMus = 0
#    print "muons. That is all."
    for i in range(tree.numMuonPF2PAT):
#        print tree.muonPF2PATPt[i]
#        print getLeptonPt(tree,i,False,LESSyst,"Pt")
        if getLeptonPt(tree,i,False,LESSyst,"Pt") < 20: continue
        if fabs(tree.muonPF2PATEta[i]) > 2.5: continue
        if not tree.muonPF2PATGlobalID[i] and not tree.muonPF2PATTrackID[i]: continue
        if tree.muonPF2PATComRelIsodBeta[i] > 0.2: continue
        numMus += 1
        selectedMuons.append(i)
    #loose muon id
    numLooseMuos = 0
    for i in range(tree.numMuonPF2PAT):
        if not tree.muonPF2PATGlobalID[i] and not tree.muonPF2PATTrackID[i]: continue
        if getLeptonPt(tree,i,False,LESSyst,"Pt") < 10: continue
        if fabs(tree.muonPF2PATEta[i]) > 2.5: continue
        if tree.muonPF2PATComRelIsodBeta[i] > 0.2: continue
        numLooseMuos +=1

    if not (numEles + numMus ==2):
        return false
#    print numEles, numMus
    lepton1 = 0
    lepton2 = 0
    if cutFlowMode=="ee" and numEles == 2 and numLooseMuos == 0 and numLooseEles==2:
        lepton1 = TLorentzVector(getLeptonPt(tree,selectedEles[0],True,LESSyst,"PX"),getLeptonPt(tree,selectedEles[0],True,LESSyst,"PY"),getLeptonPt(tree,selectedEles[0],True,LESSyst,"PZ"),getLeptonPt(tree,selectedEles[0],True,LESSyst,"E"))
        lepton2 = TLorentzVector(getLeptonPt(tree,selectedEles[1],True,LESSyst,"PX"),getLeptonPt(tree,selectedEles[1],True,LESSyst,"PY"),getLeptonPt(tree,selectedEles[1],True,LESSyst,"PZ"),getLeptonPt(tree,selectedEles[1],True,LESSyst,"E"))
        if tree.elePF2PATCharge[selectedEles[0]] * tree.elePF2PATCharge[selectedEles[1]] > 0.: return false
    elif cutFlowMode=="emu" and numEles ==1 and numMus == 1 and numLooseEles == 1 and numLooseMuos == 1:
        lepton1 = TLorentzVector(getLeptonPt(tree,selectedEles[0],True,LESSyst,"PX"),getLeptonPt(tree,selectedEles[0],True,LESSyst,"PY"),getLeptonPt(tree,selectedEles[0],True,LESSyst,"PZ"),getLeptonPt(tree,selectedEles[0],True,LESSyst,"E"))
        lepton2 = TLorentzVector(getLeptonPt(tree,selectedMuons[0],False,LESSyst,"PX"),getLeptonPt(tree,selectedMuons[0],False,LESSyst,"PY"),getLeptonPt(tree,selectedMuons[0],False,LESSyst,"PZ"),getLeptonPt(tree,selectedMuons[0],False,LESSyst,"E"))
        if tree.elePF2PATCharge[selectedEles[0]] * tree.muonPF2PATCharge[selectedMuons[0]] > 0.: return false
    elif cutFlowMode=="mumu" and numMus == 2 and numLooseMuos == 2 and numLooseEles==0:
        lepton1 = TLorentzVector(getLeptonPt(tree,selectedMuons[0],False,LESSyst,"PX"),getLeptonPt(tree,selectedMuons[0],False,LESSyst,"PY"),getLeptonPt(tree,selectedMuons[0],False,LESSyst,"PZ"),getLeptonPt(tree,selectedMuons[0],False,LESSyst,"E"))
        lepton2 = TLorentzVector(getLeptonPt(tree,selectedMuons[1],False,LESSyst,"PX"),getLeptonPt(tree,selectedMuons[1],False,LESSyst,"PY"),getLeptonPt(tree,selectedMuons[1],False,LESSyst,"PZ"),getLeptonPt(tree,selectedMuons[1],False,LESSyst,"E"))
        if tree.muonPF2PATCharge[selectedMuons[0]] * tree.muonPF2PATCharge[selectedMuons[1]] > 0.: return false
    else:
        return False
    if  fillTree: 
        fillTree.Fill(1,weight)
#        print plotList
        for plot in plotList.keys():
            if int(plotConf.get(plot,'cutStage')) < 1:
                plotList[plot]['lepSel'][dataset].Fill(float(eval(plotConf.get(plot,'fillExp'))),weight)

    if (lepton1+lepton2).M() > 20 and (cutFlowMode == "emu" or ((lepton1+lepton2).M() <81 or (lepton1+lepton2).M() > 101)) and not returnLeptons:
        if fillTree: fillTree.Fill(2,weight)
        for plot in plotList.keys():
            if int(plotConf.get(plot,'cutStage')) < 2:
#                print plotConf.get(plot,'fillExp')
                plotList[plot]['lepMass'][dataset].Fill(float(eval(plotConf.get(plot,'fillExp'))),weight)

        return True
    elif returnLeptons:
        return (lepton1,lepton2)
    return False

def makeJetCuts(tree,jerUD,jesUD,jecUnc):
    numJets = 0
    for i in range(tree.numJetPF2PAT):
        (jetPx,jetPy) = calcJERJetCorr(tree,i,jerUD)
        jesUncertainty = jecUnc.getUncertainty(sqrt(jetPx*jetPx + jetPy * jetPy),tree.jetPF2PATEta[i],jesUD)
        if jesUD == 1:
            jetPx *= (1+jesUncertainty)
            jetPy *= (1+jesUncertainty)
        if jesUD == 2:
            jetPx *= (1-jesUncertainty)
            jetPy *= (1-jesUncertainty)
        if sqrt(jetPx*jetPx + jetPy * jetPy) < 30:continue
        if fabs(tree.jetPF2PATEta[i]) > 2.5:continue
#        if tree.jetPF2PATdRClosestLepton[i] < 0.3: continue
        if tree.jetPF2PATNConstituents[i] < 2: continue
        if tree.jetPF2PATNeutralHadronEnergyFractionCorr[i] >= 0.99 or tree.jetPF2PATNeutralEmEnergyFractionCorr[i] >= 0.99: continue
        if fabs(tree.jetPF2PATEta[i]) < 2.4 and (tree.jetPF2PATChargedEmEnergyFraction[i] >= 0.99 or tree.jetPF2PATChargedHadronEnergyFraction[i] <= 0. or tree.jetPF2PATChargedMultiplicity[i] < 0.) : continue
        
        numJets += 1
#        print tree.jetPF2PATBDiscriminator[i]
    if numJets > 0:
        return True
    else:
        return False

def jetRegion(tree,runType,jerUD=0,jesUD=0,jesUnc=0,bTagUD=0,loose=False): #Returns a number that is the signal/control region the event lies in.
    "Returns the signal/control region of tree. This only makes sense if it has passed the other cuts already."
    numJets = 0
    numTags = 0
    tightJets = []
    bTags = []
    for i in range(tree.numJetPF2PAT):
        (jetPx,jetPy) = (tree.jetPF2PATPx[i],tree.jetPF2PATPy[i])
        if runType == 'mc':
            (jetPx,jetPy) = calcJERJetCorr(tree,i,jerUD)
            jesUncertainty = jesUnc.getUncertainty(sqrt(jetPx*jetPx + jetPy * jetPy),tree.jetPF2PATEta[i],jesUD)
            if jesUD == 1:
                jetPx *= (1+jesUncertainty)
                jetPy *= (1+jesUncertainty)
            if jesUD == 2:
                jetPx *= (1-jesUncertainty)
                jetPy *= (1-jesUncertainty)
        jetPt = sqrt(jetPx*jetPx + jetPy * jetPy)
        if sqrt(jetPx*jetPx + jetPy * jetPy) < 30:continue
        if fabs(tree.jetPF2PATEta[i]) > 2.4:continue
#        if tree.jetPF2PATdRClosestLepton[i] < 0.3: continue
        if tree.jetPF2PATNConstituents[i] < 2: continue
        if tree.jetPF2PATNeutralHadronEnergyFractionCorr[i] >= 0.99 or tree.jetPF2PATNeutralEmEnergyFractionCorr[i] >= 0.99: continue
        if fabs(tree.jetPF2PATEta[i]) < 2.4 and (tree.jetPF2PATChargedEmEnergyFraction[i] >= 0.99 or tree.jetPF2PATChargedHadronEnergyFraction[i] <= 0. or tree.jetPF2PATChargedMultiplicity[i] < 0.) : continue

        numJets += 1
        tightJets.append(i)
        if tree.jetPF2PATBDiscriminator[i] > 0.679:
            SF = random()
            btagSF = 0.726981*((1.+(0.253238*jetPt))/(1.+(0.188389*jetPt)))
            if bTagUD == 1:
                ptmin = [20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600]
                ptmax = [30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800]
                for j in range(len(ptmin)):
                    if jetPt > ptmin[j] and jetPt < ptmax[j]:
                        SFb_error = [  0.0554504,   0.0209663,   0.0207019,   0.0230073,   0.0208719,   0.0200453,   0.0264232,   0.0240102,   0.0229375,   0.0184615,   0.0216242,   0.0248119,   0.0465748,   0.0474666,   0.0718173,   0.0717567 ]
                        btagSF += SFb_error[j]
            if bTagUD == 2:
                ptmin = [20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600]
                ptmax = [30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800]
                for j in range(len(ptmin)):
                    if jetPt > ptmin[j] and jetPt < ptmax[j]:
                        SFb_error = [  0.0554504,   0.0209663,   0.0207019,   0.0230073,   0.0208719,   0.0200453,   0.0264232,   0.0240102,   0.0229375,   0.0184615,   0.0216242,   0.0248119,   0.0465748,   0.0474666,   0.0718173,   0.0717567 ]
                        btagSF -= SFb_error[j]
            if SF < btagSF or runType == 'data':
                numTags+=1
                bTags.append(i)
            
    #Here I guess do some loose jet vetos? Dependant on whether it has one jet or whatever.

    #The only ones of these that are interesting at the moment are 1, 4 and 5.
#    print tightJets,bTags
    return (tightJets,bTags)
    if numJets == 0 and numTags == 0:
        return 0
    if loose:
        return 2
    if numJets == 1 and numTags == 0:
        return 1
    if numJets == 1 and numTags == 1:
        if looseJets(tree) == 1: #If there are extra loose jets ignore the event.
            return 2
        else:
            return 0
    if numJets == 2 and numTags == 0:
        return 3
    if numJets == 2 and numTags == 1:
        return 4 #Should I put in some sort of loose veto in these two as well?
    if numJets == 2 and numTags == 2:
        return 5
    return 6

def numBJets(tree,runType):
    "Returns the number of events that are classified as loose in the event"
    numLooseJets = 0
    for i in range(tree.numJetPF2PAT):
        if tree.jetPF2PATPt[i] < 30:continue
        if fabs(tree.jetPF2PATEta[i]) > 2.5:continue
#        if tree.jetPF2PATdRClosestLepton[i] < 0.3: continue
        if tree.jetPF2PATNConstituents[i] < 2: continue
        if tree.jetPF2PATNeutralHadronEnergyFractionCorr[i] >= 0.99 or tree.jetPF2PATNeutralEmEnergyFractionCorr[i] >= 0.99: continue
        if fabs(tree.jetPF2PATEta[i]) < 2.4 and (tree.jetPF2PATChargedEmEnergyFraction[i] >= 0.99 or tree.jetPF2PATChargedHadronEnergyFraction[i] <= 0. or tree.jetPF2PATChargedMultiplicity[i] < 0.) : continue
        if tree.jetPF2PATBDiscriminator[i] > 0.679:
            SF = random()
            btagSF = 0.726981*((1.+(0.253238*tree.jetPF2PATPt[i]))/(1.+(0.188389*tree.jetPF2PATPt[i])))
            if SF < btagSF or runType == 'data':
                numLooseJets += 1
                
    return numLooseJets
    
def looseJets(tree,runType,jerUD=0,jesUD=0,jecUnc=0):
    "Returns the number of events that are classified as loose in the event"
    numLooseJets = 0
    for i in range(tree.numJetPF2PAT):
        (jetPx,jetPy) = (tree.jetPF2PATPx[i],tree.jetPF2PATPy[i])
        if runType == 'mc':
            (jetPx,jetPy) = calcJERJetCorr(tree,i,jerUD)
            jesUncertainty = jecUnc.getUncertainty(sqrt(jetPx*jetPx + jetPy * jetPy),tree.jetPF2PATEta[i],jesUD)
            if jesUD == 1:
                jetPx *= (1+jesUncertainty)
                jetPy *= (1+jesUncertainty)
            if jesUD == 2:
                jetPx *= (1-jesUncertainty)
                jetPy *= (1-jesUncertainty)
        jetPt = sqrt(jetPx*jetPx + jetPy * jetPy)
        if sqrt(jetPx*jetPx + jetPy * jetPy) < 20:continue
        if fabs(tree.jetPF2PATEta[i]) > 2.4:continue
#        if tree.jetPF2PATdRClosestLepton[i] < 0.3: continue
        if tree.jetPF2PATNConstituents[i] < 2: continue
        if tree.jetPF2PATNeutralHadronEnergyFractionCorr[i] >= 0.99 or tree.jetPF2PATNeutralEmEnergyFractionCorr[i] >= 0.99: continue
        if fabs(tree.jetPF2PATEta[i]) < 2.4 and (tree.jetPF2PATChargedEmEnergyFraction[i] >= 0.99 or tree.jetPF2PATChargedHadronEnergyFraction[i] <= 0. or tree.jetPF2PATChargedMultiplicity[i] < 0.) : continue
        if tree.jetPF2PATBDiscriminator[i] > 0.679:
            SF = random()
            btagSF = 0.726981*((1.+(0.253238*jetPt))/(1.+(0.188389*jetPt)))
            if SF < btagSF or runType == 'data':
                numLooseJets+=1
        
#        numLooseJets += 1

    return numLooseJets
        
#can experiment with different cuts here. This is where cuts should go - stuff like leptons and jets get their own id sections too.
def makeCuts(cutFlowMode,runType,tree,systematics,fillTreeFlag=False,fillTree=0,weight=1.0,jecUnc=0,plotList=[],dataset="",plotConf=0,LESFlag = 0): 
    "This calls the other cut functions and decides1 if the event should be processed"
    btagDirection = 0
    if systematics[4]:
        btagDirection = 1
    if systematics[5]:
        btagDirection = 2
    if not makeLeptonCuts(cutFlowMode,runType,tree,False,fillTreeFlag,fillTree,weight,plotList=plotList,dataset=dataset,plotConf=plotConf,LESSyst=LESFlag): return False
    #    print "Selected leptons"
    eventMetx = tree.metPF2PATPx
    eventMety = tree.metPF2PATPy
#    if sqrt(eventMety * eventMety + eventMetx * eventMetx)  > 45. :print "\n",eventMetx,eventMety,sqrt(eventMety * eventMety + eventMetx * eventMetx)
    if runType == 'mc':
        if systematics[8] or systematics[9]:
            uncmetx = tree.metPF2PATPx
            uncmety = tree.metPF2PATPy
            for i in range(tree.numElePF2PAT):
                uncmetx += tree.elePF2PATPX[i]
                uncmety += tree.elePF2PATPX[i]
            for i in range(tree.numMuonPF2PAT):
                uncmetx += tree.muonPF2PATPX[i]
                uncmety += tree.muonPF2PATPY[i]
            for i in range(tree.numJetPF2PAT):
                uncmetx += tree.jetPF2PATPx[i]
                uncmety += tree.jetPF2PATPy[i]
            if systematics[8]:
                eventMetx += uncmetx*0.1
                eventMety += uncmety*0.1
            if systematics[9]:
                eventMetx -= uncmetx*0.1
                eventMety -= uncmety*0.1

        if systematics[0] or systematics[1]:
            metUD = 1 if systematics[0] else 2
            (eventMetx,eventMety) = calcMETforJER(eventMetx,eventMety,tree,metUD)

        if systematics[2] or systematics[3]:
            jesUD = 1 if systematics[2] else 2
            (eventMetx,eventMety) = jecUnc.getMetAfterJESUnc(eventMetx,eventMety,tree,jesUD)
#    if sqrt(eventMety * eventMety + eventMetx * eventMetx)  > 45. :print eventMetx,eventMety,sqrt(eventMety * eventMety + eventMetx * eventMetx)
    eventMET = sqrt(eventMety * eventMety + eventMetx * eventMetx)
    if (not cutFlowMode=="emu") and eventMET < 50: return False
    if fillTreeFlag:
        fillTree.Fill(3,weight)
        (lepton1,lepton2) = makeLeptonCuts(cutFlowMode,runType,tree,True,False,LESSyst=LESFlag)
        for plot in plotList.keys():
            if int(plotConf.get(plot,'cutStage')) < 3:
                plotList[plot]['MET'][dataset].Fill(float(eval(plotConf.get(plot,'fillExp'))),weight)
        #    print "MET"
    jerUD = 0
    jesUD = 0

    if runType == 'mc':
        if systematics[0]:
            jerUD = 1
        if systematics[1]:
            jerUD = 2

        if systematics[2]:
            jesUD = 1
        if systematics[3]:
            jesUD = 2
            
    if not makeJetCuts(tree,jerUD,jesUD,jecUnc):return False
#    if cutFlowMode == "emu" and tree.elePF2PATPT[0] + tree.muonPF2PATPt[0] + tree.metPF2PATEt + tree.jetPF2PATPt[0] < 60.:return false
#    print "jet cuts!"
    #else: return True
    return True

def passesTrigger(tree,cutMode,dataset):
    "Checks that the data tree has at least one of the trigger requirements"
    if cutMode == 'emu':
        if tree.HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4 == 0: return false
        if tree.HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5 == 0: return false
        if tree.HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6 == 0: return false
        if tree.HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7 == 0: return false
        if tree.HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8 == 0: return false
        if "Prompt" in dataset and tree.HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9 == 0: return false
        return true
    if cutMode == 'ee':
        if tree.HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15 == 1: return true
        if tree.HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16 == 1: return true
        if tree.HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17 == 1: return true
        if tree.HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18 == 1: return true
#        if "Prompt" in dataset and tree.HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v19 == 1: return true
#        if tree.HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v19 == 1: return true
        
    if cutMode == 'mumu':
        if tree.HLT_Mu17_Mu8_v14 == 0: return false
        if tree.HLT_Mu17_Mu8_v15 == 0: return false
        if tree.HLT_Mu17_Mu8_v16 == 0: return false
        if tree.HLT_Mu17_Mu8_v17 == 0: return false
        if tree.HLT_Mu17_Mu8_v18 == 0: return false
        if tree.HLT_Mu17_Mu8_v19 == 0: return false
        if "Prompt" in dataset and tree.HLT_Mu17_Mu8_v20 == 0: return false
        if "Prompt" in dataset and tree.HLT_Mu17_Mu8_v21 == 0: return false
        if "Prompt" in dataset and not (tree.HLT_Mu17_Mu8_v20 == 1 or tree.HLT_Mu17_Mu8_v21 == 1 or tree.HLT_Mu17_Mu8_v19 == 1 or tree.HLT_Mu17_Mu8_v18 == 1 or tree.HLT_Mu17_Mu8_v17 == 1 or tree.HLT_Mu17_Mu8_v16 == 1 or tree.HLT_Mu17_Mu8_v15 == 1 or tree.HLT_Mu17_Mu8_v14 == 1): print "Doesn't have any of the obvious triggers - might be missing triggers"
#        if tree.HLT_Mu17_Mu8_v20 == 0: return false
#        if tree.HLT_Mu17_Mu8_v21 == 0: return false
        return true
    return false
        
def makeLooseCuts(cutFlowMode,runType,tree): #can experiment with different cuts here. This is where cuts should go - stuff like leptons and jets get their own id sections too.
    "This calls the other cut functions and decides if the event should be processed"
    if not makeLeptonCuts(cutFlowMode,tree): return False
    if jetRegion(tree,runType,True) == 2:return True
    else: return False

def calcMETforJER(initMetx,initMety,tree,metUD):
    "Recalculate the MET of the event using JER up or down"
    for i in range(tree.numJetPF2PAT):
        if tree.genJetPF2PATPT[i] > 15.0:
            initMetx += tree.jetPF2PATPx[i]
            initMety += tree.jetPF2PATPy[i]
            (jetPx,jetPy) = calcJERJetCorr(tree,i,metUD)
            initMetx -= jetPx
            initMety -= jetPy
    return (initMetx,initMety)

def calcJERJetCorr(tree,index,metUD):
    "Recalculates jet pt for JER systematic"
    if metUD == 0 or tree.genJetPF2PATPX[index] < -998.:
        return (tree.jetPF2PATPx[index],tree.jetPF2PATPy[index])
    oldCorrFact = 0.
    newCorrFact = 0.
    if tree.jetPF2PATEta[index] <= 1.1:
        oldCorrFact = 0.066
        newCorrFact = -0.006 if metUD == 1 else 0.136
    elif tree.jetPF2PATEta[index] <= 1.7:
        oldCorrFact = 0.191
        newCorrFact = 0.129 if metUD == 1 else 0.251
    elif tree.jetPF2PATEta[index] <= 2.3:
        oldCorrFact = 0.096
        newCorrFact = 0.011 if metUD == 1 else 0.176
    elif tree.jetPF2PATEta[index] <= 5.0:
        oldCorrFact = 0.166
        newCorrFact = -0.033 if metUD == 1 else 0.356
    else:
        oldCorrFact = 0.166
        newCorrFact = -0.033 if metUD == 1 else 0.356
    newPx =(1+newCorrFact)* (tree.jetPF2PATPx[index] + tree.genJetPF2PATPX[index] * oldCorrFact)/(1+oldCorrFact) - newCorrFact*tree.genJetPF2PATPX[index]
    newPy =(1+newCorrFact)* (tree.jetPF2PATPy[index] + tree.genJetPF2PATPY[index] * oldCorrFact)/(1+oldCorrFact) - newCorrFact*tree.genJetPF2PATPY[index]
    return (newPx, newPy)

def ptsysForPlotting(tree,cutMode,runType):
    jets = numberOfJets(tree,True)
    (lepton1,lepton2) = makeLeptonCuts(cutMode,runType,tree,True)
    sumJetPx = 0.
    sumJetPy = 0.
    for i in range(len(jets)):
        sumJetPx += tree.jetPF2PATPx[i]
        sumJetPy += tree.jetPF2PATPy[i]
    ptx = lepton1.Px() + lepton2.Px() + tree.metPF2PATPx + sumJetPx
    pty = lepton1.Py() + lepton2.Py() + tree.metPF2PATPy + sumJetPy
    return sqrt(ptx * ptx + pty + pty)

def htForPlotting(tree,cutMode,runType):
    "Gives the Ht of the event, defined as the sum of pt of the leptons, MET and all jets that pass jetID"
    jets = numberOfJets(tree,True)
    (lepton1,lepton2) = makeLeptonCuts(cutMode,runType,tree,True)
    sumJetPt = 0.
    for i in range(len(jets)):
        sumJetPt += tree.jetPF2PATPt[jets[i]]
    return lepton1.Pt() + lepton2.Pt() + tree.metPF2PATEt + sumJetPt


def leadingJetPt(tree):
    "Returns the pt of the first jet that passes jet id criteria."
    for i in range(tree.numJetPF2PAT):
        if tree.jetPF2PATPt[i] < 30:continue
        if fabs(tree.jetPF2PATEta[i]) > 2.5:continue
#        if tree.jetPF2PATdRClosestLepton[i] < 0.3: continue
        if tree.jetPF2PATNConstituents[i] < 2: continue
        if tree.jetPF2PATNeutralHadronEnergyFractionCorr[i] >= 0.99 or tree.jetPF2PATNeutralEmEnergyFractionCorr[i] >= 0.99: continue
        if fabs(tree.jetPF2PATEta[i]) < 2.4 and (tree.jetPF2PATChargedEmEnergyFraction[i] >= 0.99 or tree.jetPF2PATChargedHadronEnergyFraction[i] <= 0. or tree.jetPF2PATChargedMultiplicity[i] < 0.) : continue
        
        return tree.jetPF2PATPt[i]
    return -1

def numberOfJets(tree,returnJets=False):

    numJets = 0
    jets = []
    for i in range(tree.numJetPF2PAT):
        if tree.jetPF2PATPt[i] < 30:continue
        if fabs(tree.jetPF2PATEta[i]) > 2.5:continue
#        if tree.jetPF2PATdRClosestLepton[i] < 0.3: continue
        if tree.jetPF2PATNConstituents[i] < 2: continue
        if tree.jetPF2PATNeutralHadronEnergyFractionCorr[i] >= 0.99 or tree.jetPF2PATNeutralEmEnergyFractionCorr[i] >= 0.99: continue
        if fabs(tree.jetPF2PATEta[i]) < 2.4 and (tree.jetPF2PATChargedEmEnergyFraction[i] >= 0.99 or tree.jetPF2PATChargedHadronEnergyFraction[i] <= 0. or tree.jetPF2PATChargedMultiplicity[i] < 0.) : continue
       
        numJets += 1
        jets.append(i)
        
    if not returnJets:
        return numJets
    else:
        return jets

def getLeptonPt(tree,index,ele,LES,direction="PT"):
    """
    The important thing here is the direction - 0 is for pt and 1,2,3 are for px,y and z respectively.
    LES is true if there needs to be an adjustment. False means it just returns the pt in the tree.
    ele is true for an electron, false for a muon
    I've now generalised this method so that it just takes an argument direction and scales that.
    This means it can now be used for energy as well as anything else you arbitrarily want to scale.
    """
    directionMap = ["PT","PX","PY","PZ"]
    scale = 1.0
    if LES:
        if ele:
            if fabs(tree.elePF2PATEta[index]) < 1.4442:
                scale = 1.005 if LES == 1 else 0.995
            else:
                scale = 1.01 if LES ==1 else 0.99
        else:
            scale = 1.002 if LES == 1 else 0.998
    if ele:
        return scale * eval("tree.elePF2PAT"+direction+"["+str(index)+"]")
    else:
        return scale * eval("tree.muonPF2PAT"+direction+"["+str(index)+"]")
