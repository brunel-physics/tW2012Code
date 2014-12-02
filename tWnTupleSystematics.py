from ROOT import *
import sys
from tWnTupleCutsSyst import *
from optparse import OptionParser
import ConfigParser
from setTDRStyle import setTDRStyle
from jetCorrectionUncertainty import JetCorrectionUncertainty
from pileupReweighting import PileupReweighting
from histogramPlotter import HistogramPlotter
#from leptonSF import leptonSF

def getZWeight(cutMode,tree):
    if cutMode == "mumu":
        if tree.metPF2PATPt < 10.: return 0.8841
        if tree.metPF2PATPt < 20.: return 0.9386
        if tree.metPF2PATPt < 30.: return 1.0131
        if tree.metPF2PATPt < 40.: return 1.1012
        if tree.metPF2PATPt < 50.: return 1.1850
        if tree.metPF2PATPt < 60.: return 1.2500
        else: return 1.3071
    if cutMode == "ee":
        if tree.metPF2PATPt < 10.: return 0.9215
        if tree.metPF2PATPt < 20.: return 0.9608
        if tree.metPF2PATPt < 30.: return 1.0246
        if tree.metPF2PATPt < 40.: return 1.0964
        if tree.metPF2PATPt < 50.: return 1.1633
        if tree.metPF2PATPt < 60.: return 1.2529
        else: return 1.2194
    if cutMode == "emu":
        if tree.metPF2PATPt < 10.: return 0.9028
        if tree.metPF2PATPt < 20.: return 0.9497
        if tree.metPF2PATPt < 30.: return 1.0189
        if tree.metPF2PATPt < 40.: return 1.0988
        if tree.metPF2PATPt < 50.: return 1.17415
        if tree.metPF2PATPt < 60.: return 1.25145
        else: return 1.26325

def pileupSystWeight(pileupWeight,puA,puUD):
    ibin = 0
    for i in range(puA.GetNbinsX()):
        if pileupWeight == puA.GetBinContent(i):
            ibin = i
    p0_minus = [ 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. ]
    p1_minus = [-0.677786, -0.619614, -0.49465, -0.357963, -0.238359, -0.110002, 0.0348629, 0.191263, 0.347648, 0.516615, 0.679646, 0.836673, 0.97764, 1.135, 1.29922, 1.42467, 1.55901, 1.61762, 1.67275, 1.96008]
    p2_minus = [ 0.526164, 0.251816, 0.11049, 0.026917, -0.0464692, -0.087022, -0.0931581, -0.0714295, -0.0331772, 0.0347473, 0.108658, 0.193048, 0.272314, 0.376357, 0.4964, 0.58854, 0.684959, 0.731063, 0.760044, 1.02386]
    p1_expoM = [ 1.63363e-03, 6.79290e-04, 3.69900e-04, 2.24349e-04, 9.87156e-06]
    p2_expoM = [ 2.64692, 3.26585, 3.53229, 4.18035, 5.64027]
    p0_plus = [ 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. ]
    p1_plus = [ -0.739059, -0.594445, -0.477276, -0.359707, -0.233573, -0.103458, 0.0373401, 0.176571, 0.337617, 0.499074, 0.675126, 0.840522, 1.00917, 1.15847, 1.23816, 1.44271, 1.52982, 1.46385, 1.5802, 0.988689]
    p2_plus = [ 0.208068, 0.130033, 0.0850356, 0.0448344, 0.000749832, -0.0331347, -0.0653281, -0.0746009, -0.0800667, -0.0527636, -0.00402649, 0.103338, 0.261261, 0.491084, 0.857966, 1.19495, 1.75071, 2.65559, 3.35433, 5.48835]
    p1_expoP = [ 1.42463e-01, 4.18966e-02, 1.12697e-01, 1.66197e-01, 1.50768e-01]
    p2_expoP = [ 1.98758, 2.27217, 2.26799, 2.38455, 2.52428]

    Shift = 0
    if puUD == 1:
        Shift = 0.6
    else:
        Shift = -0.6

    Pweight = 1.

    if ibin >= 0 and ibin < 20:
        if Shift < 0.:
            Pweight = p0_minus[ibin] + p1_minus[ibin]*Shift + p2_minus[ibin]*Shift*Shift
        else:
            Pweight = p0_plus[ibin] + p1_plus[ibin]*Shift + p2_plus[ibin]*Shift*Shift
    elif ibin < 25.:
        if Shift < 0.:
            Pweight = p1_expoM[ibin-20]*exp(p2_expoM[ibin-20]*Shift)
        else:
            Pweight = p1_expoP[ibin-20]*exp(p2_expoP[ibin-20]*Shift)
    else:
        Pweight = 0.
#    print ibin, pileupWeight,Pweight
    return Pweight
    
    
def main():
    parser = OptionParser("usage: %prog [options] arg")
    parser.add_option("-c","--config",dest="configFile",help="Config file containing list of files for each dataset")
    parser.add_option("-m","--mode",dest="cutMode",help="Cut mode to be used on the full dataset")
    parser.add_option("-o","--outFolder",dest="outFolder",default="",help="Name of the folder to save the cut flows (root and png) to")
    parser.add_option("-e","--postfix",dest="postfix",default="ANotHelpfulPostfix",help="A sensible name that will allow you to identify the cut flows")
    parser.add_option("-s","--sysOut",dest="sysOut",help="Folder the one bin histograms will be put into")
    parser.add_option("-p","--postOBH",dest="obhPostfix",default="",help = "Post fix on the one bin histograms. This means the systematic and up or down.")
    parser.add_option("-d","--data",action="store_true",dest="dataFlag",default=False,help="Only run on data")
    parser.add_option("-a","--saveSkim",dest="skimSaveFile",default="",help="If you want to save the skims, set this option to the out-folder on data0.")
    parser.add_option("-q","--mc",action="store_true",dest="mcFlag",default=False,help="Only run on MC")
    parser.add_option("-r","--reweightPU",action="store_true",dest="reweightPileup",default=False,help="Use new distribution to reweight pileup in MC")
    parser.add_option("-x","--plotConfig",dest="plotConfig",help="A configuration containing plot information")
    parser.add_option("-l","--lumiSpecified",action="store_true",dest="lumiSpecified",default=False,help="Set this flag if running over a specific run range and have lumi information for the data in the config.")
    parser.add_option("-k","--recordEventNumbers",action="store_true",default=False,dest="eventDumpFlag",help="Set this flag to record the event run,lumi and number for each region of passing events.")
    parser.add_option("-t","--ttbarReweight",action="store_true",default=False,dest="ttbarFlag",help="Set this flag for doing the top pt reweighting on ttbar. Will also produce systematic samples.")
    parser.add_option("--JERPlus",action="store_true",dest="jerUFlag",default=False,help="JER up systematic")
    parser.add_option("--JERMinus",action="store_true",dest="jerDFlag",default=False,help="JER down systematic")
    parser.add_option("--JESPlus",action="store_true",dest="jesUFlag",default=False,help="JES up systematic")
    parser.add_option("--JESMinus",action="store_true",dest="jesDFlag",default=False,help="JES down systematic")
    parser.add_option("--btagsUp",action="store_true",dest="btagUFlag",default=False,help="b-tag SF up")
    parser.add_option("--btagsDown",action="store_true",dest="btagDFlag",default=False,help="b-tag SF down")
    parser.add_option("--pileupUp",action="store_true",dest="puUFlag",default=False,help="PU up")
    parser.add_option("--pileupDown",action="store_true",dest="puDFlag",default=False,help="PU down")
    parser.add_option("--metUp",action="store_true",dest="metUFlag",default=False,help="unclustered MET Up")
    parser.add_option("--metDown",action="store_true",dest="metDFlag",default=False,help="Unclustered MET Down")
    parser.add_option("--lepSFUp",action="store_true",dest="lepSFUFlag",default=False,help="Lepton scale Factor Up")
    parser.add_option("--lepSFDown",action="store_true",dest="lepSFDFlag",default=False,help="Lepton Scale Factor Down")
    parser.add_option("--zPlusJetsUp",action="store_true",dest="zPlusJetsUp",default=False,help="Scale up lepton scale factor")
    parser.add_option("--zPlusJetsDown",action="store_true",dest="zPlusJetsDown",default=False,help="Scale down z+jets scale factor")
    parser.add_option("--LESUp",action="store_true",dest="LESUp",default=False,help="Scale up lepton energies")
    parser.add_option("--LESDown",action="store_true",dest="LESDown",default=False,help="Scale Down lepton energies")

    (options, args) = parser.parse_args()
    systs = [options.jerUFlag,options.jerDFlag,options.jesUFlag,options.jesDFlag,options.btagUFlag,options.btagDFlag,options.puUFlag,options.puDFlag,options.metUFlag,options.metDFlag,options.lepSFUFlag,options.lepSFDFlag,options.zPlusJetsUp,options.zPlusJetsDown,options.LESUp,options.LESDown]
    print systs
    
    datasets = []
    config = ConfigParser.RawConfigParser()

    plotConf = ConfigParser.RawConfigParser()
    plots = []
#    print options.plotConfig
    if options.plotConfig:
        plotConf.read(options.plotConfig)
        plots = plotConf.sections()

    #config file will contain a list of sections that each correspond to a dataset
    if options.configFile:
        config.read(options.configFile)
        datasets = config.sections()


    sampleNamesForPlots = ['tW','zPlusJets','ttbar','other', 'data']
    cutFlowStages = ["lepSel","lepMass","MET","1jet","btag","Ht","1j1t","2j1t","2j2t"]
    plotList = {}

    if options.eventDumpFlag:
        eventDumpFiles = {}
        for dataset in datasets:
            eventDumpFiles[dataset] = {}
            for region in ['1j1t','2j1t','2j2t']:
                eventDumpFiles[dataset][region] = open("eventDumpLists/"+dataset+"_"+region + "_eventDump.txt","w")

    for plot in plots:
        plotList[plot] = {}
        for stage in cutFlowStages:
            plotList[plot][stage] = {}
            for sample in sampleNamesForPlots:
                plotList[plot][stage][sample] = TH1F(sample + "_" + plot + "_" + stage, sample + "_" + plot + "_" + stage, int(plotConf.get(plot,"bins")),float(plotConf.get(plot,"minX")),float(plotConf.get(plot,"maxX")) )

    #Initialise the one bin histograms here.
    regionNames = ['1j1t','2j1t','2j2t']    
    sampleNames = ['twdr','tt','other','DATA']
    systematicName = ""
    lesFlag = 0
    upDown = ""
    if systs[0]:
        systematicName = "JER"
        upDown = "plus"
    elif systs[1]:
        systematicName = "JER"
        upDown = "minus"
    elif systs[2]:
        systematicName = "JES"
        upDown = "plus"
    elif systs[3]:
        systematicName = "JES"
        upDown = "minus"
    elif systs[4]:
        systematicName = "BtagSF"
        upDown = "plus"
    elif systs[5]:
        systematicName = "BtagSF"
        upDown = "minus"
    elif systs[6]:
        systematicName = "PU"
        upDown = "plus"
    elif systs[7]:
        systematicName = "PU"
        upDown = "minus"
    elif systs[8]:
        systematicName = "UnclusteredMET"
        upDown = "plus"
    elif systs[9]:
        systematicName = "UnclusteredMET"
        upDown = "minus"
    elif systs[10]:
        systematicName = "lepSF"
        upDown = "plus"
    elif systs[11]:
        systematicName =  "lepSF"
        upDown = "minus"
    elif systs[12]:
        systematicName = "zJetsSF"
        upDown = "plus"
    elif systs[13]:
        systematicName = "zJetsSF"
        upDown = "minus"
    elif systs[14]:
        systematicName = "LES"
        upDown = "plus"
        lesFlag = 1
    elif systs[15]:
        systematicName = "LES"
        lesFlag = 2
        upDown = "minus"
    oneBinHistos = {}
    histoPostfix = ""
    if upDown == "plus" or upDown == "minus":
        histoPostfix = "__" + systematicName + "__" + upDown
    if not options.obhPostfix == "":
        histoPostfix = options.obhPostfix
    for sample in sampleNames:
        oneBinHistos[sample] = {}
        for region in regionNames:
            oneBinHistos[sample][region] = TH1F(options.cutMode + region + "__" + sample + histoPostfix,options.cutMode + region + "__" + sample + histoPostfix,1,0.5,1.5)
            #oneBinHistos[sample][region].Sumw2()

    oneBinHistosFilled = {}
    for sample in sampleNames:
        oneBinHistosFilled[sample] = False
            
    cutFlows = {}
    for dataset in ['data','tW','ttbar','zPlusJets','other']:
        cutFlows[dataset]=TH1F(dataset+"_cutFlow_"+options.cutMode,dataset+"_cutFlow_"+options.cutMode,6,1,7)

    regionPlots = {}
    for dataset in ['data','tW','ttbar','zPlusJets','other']:
        regionPlots[dataset]=TH1F(dataset+"_regionPlot_"+options.cutMode,dataset+"_regionPlot_"+options.cutMode,3,1,4)

    #Add in some extra one bin histograms only if doing top pt reweighting
    ttbarSystHist = {}
    if options.ttbarFlag:
        for region in regionNames:
            ttbarSystHist[region] = {}
            ttbarSystHist[region]["plus"] = TH1F(options.cutMode + region + "__tt__topPt__plus",options.cutMode + region + "__tt__topPt__plus",1,0.5,1.5)
            ttbarSystHist[region]["minus"] = TH1F(options.cutMode + region + "__tt__topPt__minus",options.cutMode + region + "__tt__topPt__minus",1,0.5,1.5)

    totalLumi = 0.
    lumiA = 0.
    lumiB = 0.
    lumiC = 0.

    lumiA = 808.472 + 82.136

    if not options.lumiSpecified:
        if options.cutMode == "mumu":
            lumiB = 4426.
        else:
            lumiB = 4429.
    
        if options.cutMode == "emu":
            lumiC = 495.003 + 6401.    
        else:
            lumiC = 486.186 + 6396.
        if options.cutMode == "ee":
            lumiC = 486.168 + 6401.
    #If running over a specific nrun range set the lumi here so MC is properly normalised.oith three, Rockstar has the flexibility to move between them to pace the both story and the gameplay. Seeing Michael beat up a bunch of gang members with a baseball bat might not gel too well with the reformed-gangster-trying-to-be-good that we see in his cutscenes, but Trevor after a few too many whiskies? I can definitely imagine that.
        
    totalLumi = lumiA + lumiB + lumiC

#    for dataset in datasets:
#        if config.get(dataset,'runType') == 'mc':continue
#        totalLumi += float(config.get(dataset,"luminosity")) 
#        if 'A13Jul' in dataset or 'A06Aug' in dataset:
#            lumiA += float(config.get(dataset,"luminosity")) 
#        if 'B13Jul' in dataset:
#            lumiB += float(config.get(dataset,"luminosity"))
#        if 'C24Aug' in dataset or 'CPrompt' in dataset:
#            lumiC += float(config.get(dataset,"luminosity")) 

    #Set up reweighting of pileup if using the data-driven distributions
    pileupA, pileupB, pileupC = 0 ,0,0#These two are used by both reweighting and systematic options, so are left outside the if statements.
    mcName = "pileupHistos/pileup_MC_Summer12.root"
    getNPV = PileupReweighting("pileupHistos/run2012A_13Jul.root","pileupHistos/pileup_MC_Summer12.root")
    if options.reweightPileup:
        pileupA = PileupReweighting("pileupHistos/systematics/run2012A_13Jul.root",mcName)
        pileupB = PileupReweighting("pileupHistos/systematics/run2012B_13Jul.root",mcName)
        pileupC = PileupReweighting("pileupHistos/systematics/run2012C_v2.root",mcName)


            #dummy variables for pileup systematic.

    pileupAUp, pileupBUp, pileupCUp,pileupADown, pileupBDown, pileupCDown = 0,0,0,0,0,0
    

    if systs[6] or systs[7]: #set up pile-up histograms. Will need these to get npv later on.
        #Initialise pileup reweighting stuff.
        pileupAUp = PileupReweighting("pileupHistos/systematics/run2012A_13Jul_Up.root",mcName)
        pileupBUp = PileupReweighting("pileupHistos/systematics/run2012B_13Jul_Up.root",mcName)
        pileupCUp = PileupReweighting("pileupHistos/systematics/run2012C_v2_Up.root",mcName)
        pileupADown = PileupReweighting("pileupHistos/systematics/run2012A_13Jul_Down.root",mcName)
        pileupBDown = PileupReweighting("pileupHistos/systematics/run2012B_13Jul_Down.root",mcName)
        pileupCDown = PileupReweighting("pileupHistos/systematics/run2012C_v2_Down.root",mcName)
        
    jecUnc = 0
    jecUnc = JetCorrectionUncertainty("JESscaleFactors.txt")


#This is here in case I ned the eta dependent scale factor. Hopefully I won't, as the other way is far easier.
#    lepSFs = 0
#    dataFilelepsf = ""
#    if options.cutMode == "ee":
#        dataFilelepsf = "trigger_SF_ee.root"
##    elif options.cutMode == "emu":
 #       dataFilelepsf = "trigger_SF_emu.root"
 #   elif options.cutMode == "mumu":
 #       dataFilelepsf = "trigger_SF_mumu.root"
 #   
 #   lepSFs = leptonSF(dataFilelepsf)

    leptonScaleFactor = 1.0
    if options.cutMode == "ee":
        leptonScaleFactor *= (0.9623 * 0.9623 * 0.975)
        if systs[10]:
            leptonScaleFactor += 0.021
        if systs[11]:
            leptonScaleFactor -= 0.021
    elif options.cutMode == "emu":
        leptonScaleFactor *= (0.9623 * 0.999 * 0.953)
        if systs[10]:
            leptonScaleFactor += .017
        if systs[11]:
            leptonScaleFactor -= .017
    elif options.cutMode == "mumu":
        leptonScaleFactor *= (0.999 * 0.999 * 0.965)
        if systs[10]:
            leptonScaleFactor += 0.022
        if systs[11]:
            leptonScaleFactor -= 0.022

    foundCount = {}
    foundCountWeighted = {}
    for dataset in datasets:
        if options.dataFlag and config.get(dataset,"runType") == "mc" : continue
        if options.mcFlag and config.get(dataset,"runType") == "data" : continue

        #Here I will build a TChain for each dataset I guess, coz it looks like that might make this work better.
        tree = TChain("tree")

        files = open(config.get(dataset,"fileName"),"r")
        
        for file in files:
            tree.Add(file[:-1])

        if not options.skimSaveFile == "":
            outFile = TFile("/data0/tW2012/" + options.skimSaveFile + options.cutMode + "_" + dataset +  "_fullSkim.root","RECREATE")
        outTree = tree.CloneTree(0)

        weight = 1.0
        numberOfAcceptedEvents = 0
        
        ttbarWeight = 1.0
        #If processing top pt reweighted events, find the overall weight here.
        if options.ttbarFlag:
            ttbarFiles = open(config.get(dataset,"fileName"),"r")
            ttbarTopWeight = 0.
            totalEvents = 0
            for file in ttbarFiles:
                tempFile = TFile(file[:-1],"READ")
                ttbarTopWeight += (tempFile.Get("topPtWeightSum")).GetBinContent(1)
                totalEvents += (tempFile.Get("eventcount")).GetBinContent(1)
            ttbarWeight = ttbarTopWeight / float(totalEvents)
            print "Average ttbar weight: {0}".format(ttbarWeight)
        

        print dataset
        if config.get(dataset,"runType") == 'mc':
            if 'tW' in dataset or 'tbarW' in dataset:
                regionPlotName = 'tW'
                sampleName  = 'twdr'
            elif 'zPlusJets' in dataset:
                regionPlotName = 'zPlusJets'
                sampleName = 'other'
            elif 'ttbar' in dataset:
                regionPlotName = 'ttbar'
                sampleName = 'tt'
            else :
                regionPlotName = 'other'
                sampleName = 'other'
            weight = totalLumi * float(config.get(dataset,"crossSection")) / (float(config.get(dataset,"totalEvents")) * ttbarWeight)
            print "Dataset {0} contains {1} events, cross sections of {2} and therefore a weight of {3}".format(dataset,config.get(dataset,"totalEvents"),config.get(dataset,"crossSection"),weight)
            oneBinHistosFilled[sampleName] = True
        if config.get(dataset,"runType") == 'data':
            weight = 1.0
            regionPlotName = 'data'
            sampleName = 'DATA'
            oneBinHistosFilled["DATA"] = True

        runType = config.get(dataset,"runType")
        cutFlowMode = options.cutMode

        foundCount[dataset] = 0
        foundCountWeighted[dataset]=0

        for i in range(tree.GetEntries()):
            sys.stdout.write("\rProcessing event: {0} [{1:.1%}]    Found {2} events, weighted at {3:.2f}".format(str(i) +"/" + str(tree.GetEntries()),(i * 1.) / tree.GetEntries(),foundCount[dataset],foundCountWeighted[dataset]))
            sys.stdout.flush()
            tree.GetEntry(i)
            eventWeight = weight
            if (config.get(dataset,"runType")) == 'mc':
                pileupweight = 1.0
                #This is decide whether to use the nominal pileup distributions (the ones I've been using) or the data-driven ones I got from Danny.
                if not options.reweightPileup:
                    pileupweight = (tree.PileUpWeightRunA * lumiA + tree.PileUpWeightRunB * lumiB + tree.PileUpWeightRunC * lumiC) / totalLumi
                else:
                    npv = getNPV.getNpv(tree.PileUpWeightRunA)
                    pileupweight = (pileupA.weight(npv) * lumiA + pileupB.weight(npv) * lumiB + pileupC.weight(npv) * lumiC) / totalLumi
                puUD = 0
                if systs[6] or systs[7]:
                    if systs[6]:
                        pileupweight = (pileupAUp.weight(npv) * lumiA + pileupBUp.weight(npv) * lumiB + pileupCUp.weight(npv) * lumiC) / totalLumi
                    else:
                        pileupweight = (pileupADown.weight(npv) * lumiA + pileupBDown.weight(npv) * lumiB + pileupCDown.weight(npv) * lumiC) / totalLumi
#                    pileupweight *= pileupSystWeight(tree.PileUpWeightRunA,puA,puUD)
                eventWeight *= pileupweight
                eventWeight *= leptonScaleFactor
                if options.ttbarFlag:
                    eventWeight *= tree.topPtReweight
                if 'zPlusJets' in dataset:
                    zWeight = getZWeight(options.cutMode,tree) - 1
#                    print "\n",getZWeight(options.cutMode,tree),zWeight, eventWeight,
                    if systs[12]:
                        zWeight *= 2.
                    if systs[13]:
                        zWeight = 0.
#                    print zWeight,
                    eventWeight *= 1 + zWeight
#                    print eventWeight
            if config.get(dataset,"runType") == 'data' and not passesTrigger(tree,options.cutMode,dataset): continue
#            print plotList
#            print plots
            if makeCuts(options.cutMode,config.get(dataset,"runType"),tree,systs,True,cutFlows[regionPlotName],eventWeight,jecUnc,plotList,regionPlotName,plotConf,LESFlag=lesFlag):
                
                outTree.Fill()
                jerUD = 0
                if systs[0]:
                    jerUD = 1
                if systs[1]:
                    jerUD = 2
                jesUD = 0
                if systs[2]:
                    jesUD = 1
                if systs[3]:
                    jesUD = 2
                bTagUD = 0
                if systs[4]:
                    bTagUD = 1
                if systs[5]:
                    bTagUD = 2
                (jets,tags)=jetRegion(tree,config.get(dataset,"runType"),jerUD,jesUD,jecUnc,bTagUD = bTagUD)
                (lepton1,lepton2) = makeLeptonCuts(options.cutMode,config.get(dataset,"runType"),tree,True,LESSyst=lesFlag)
                pxsys = 0.
                pysys = 0.
                htSys = 0.
                for i in range(len(jets)):
                    jetPx = tree.jetPF2PATPx[i]
                    jetPy = tree.jetPF2PATPy[i]
                    if config.get(dataset,"runType") == 'mc':
                        (jetPx,jetPy) = calcJERJetCorr(tree,jets[i],jerUD)
                        #jetPy = calcJERJetCorr(tree,jets[i],jerUD)[1]
                        jesUncertainty = jecUnc.getUncertainty(sqrt(jetPx*jetPx + jetPy * jetPy),tree.jetPF2PATEta[i],jesUD)
                        if jesUD == 1:
                            jetPx *= (1+jesUncertainty)
                            jetPy *= (1+jesUncertainty)
                        elif jesUD == 2:
                            jetPx *= (1-jesUncertainty)
                            jetPy *= (1-jesUncertainty)
                    pxsys += jetPx
                    pysys += jetPy
                    htSys += sqrt(jetPx * jetPx + jetPy * jetPy)
                pxsys += lepton1.Px() + lepton2.Px()
                pysys += lepton1.Py() + lepton2.Py()
                htSys += lepton1.Pt() + lepton2.Pt()
                eventMetx = tree.metPF2PATPx
                eventMety = tree.metPF2PATPy
                if config.get(dataset,"runType") == 'mc':
                    if systs[8] or systs[9]:
                        uncmetx = tree.metPF2PATPx
                        uncmety = tree.metPF2PATPy
                        for i in range(tree.numElePF2PAT):
                            uncmetx += getLeptonPt(tree,i,True,LES=lesFlag,direction="PX")
                            uncmety += getLeptonPt(tree,i,True,LES=lesFlag,direction="PY")
                        for i in range(tree.numMuonPF2PAT):
                            uncmetx += getLeptonPt(tree,i,False,LES=lesFlag,direction="PX")
                            uncmety += getLeptonPt(tree,i,False,LES=lesFlag,direction="PY")
                        for i in range(tree.numJetPF2PAT):
                            uncmetx += tree.jetPF2PATPx[i]
                            uncmety += tree.jetPF2PATPy[i]
                    if systs[8]:
                        eventMetx += uncmetx*0.1
                        eventMety += uncmety*0.1
                    if systs[9]:
                        eventMetx -= uncmetx*0.1
                        eventMety -= uncmety*0.1

                    if systs[0] or systs[1]:
                        if systs[0]:
                            metUD = 1
                        else:
                            metUD = 2
                        (eventMetx,eventMety) = calcMETforJER(eventMetx,eventMety,tree,metUD)

                    if systs[2] or systs[3]:
                        if systs[2]:
                            jesUD = 1
                        else:
                            jesUD = 2
                        (eventMetx,eventMety) = jecUnc.getMetAfterJESUnc(eventMetx,eventMety,tree,jesUD)
    
                pxsys += eventMetx
                pysys += eventMety
                htSys += sqrt(eventMetx * eventMetx + eventMety * eventMety)

                #Fill in the one bin histograms.
                if len(jets) == 1 and looseJets(tree,config.get(dataset,"runType"),jerUD,jesUD,jecUnc) == 1 and len(tags) and (htSys > 160. or not options.cutMode == "emu") :
                    if options.eventDumpFlag:
                        eventDumpFiles[dataset]["1j1t"].write(str(tree.eventRun)+","+str(tree.eventLumiblock)[:-2]+","+str(tree.eventNum)+"\n")
                    oneBinHistos[sampleName]['1j1t'].Fill(1.,eventWeight)
                    foundCount[dataset]+=1
                    foundCountWeighted[dataset]+=eventWeight
                    regionPlots[regionPlotName].Fill(1,eventWeight)
                    if options.ttbarFlag:
                        ttbarSystHist['1j1t']["minus"].Fill(1.,eventWeight/tree.topPtReweight)
                        tempWeight = tree.topPtReweight - 1
                        ttbarSystHist['1j1t']["plus"].Fill(1.,eventWeight*((1 + (2. * tempWeight))/(1+tempWeight)))
#                        print "   ",eventWeight, tree.topPtReweight, eventWeight*((1 + (2. * tempWeight))/(1+tempWeight)), eventWeight/tree.topPtReweight
                    for plot in plots:
                        plotList[plot]['1j1t'][regionPlotName].Fill(float(eval(plotConf.get(plot,'fillExp'))),eventWeight)

                if len(jets) == 2 and len(tags) == 1 and (htSys > 160. or not options.cutMode == "emu"):
                    if options.eventDumpFlag:
                        eventDumpFiles[dataset]["2j1t"].write(str(tree.eventRun)+","+str(tree.eventLumiblock)[:-2]+","+str(tree.eventNum)+"\n")
                    oneBinHistos[sampleName]['2j1t'].Fill(1.,eventWeight)
                    regionPlots[regionPlotName].Fill(2,eventWeight)
                    if options.ttbarFlag:
                        ttbarSystHist['2j1t']["minus"].Fill(1.,eventWeight/tree.topPtReweight)
                        tempWeight = tree.topPtReweight - 1
                        ttbarSystHist['2j1t']["plus"].Fill(1.,eventWeight*((1 + (2. * tempWeight))/(1+tempWeight)))
                    for plot in plots:
                        plotList[plot]['2j1t'][regionPlotName].Fill(float(eval(plotConf.get(plot,'fillExp'))),eventWeight)
                if len(jets) == 2 and len(tags) == 2 and (htSys > 160. or not options.cutMode == "emu"):
                    if options.eventDumpFlag:
                        eventDumpFiles[dataset]["2j2t"].write(str(tree.eventRun)+","+str(tree.eventLumiblock)[:-2]+","+str(tree.eventNum)+"\n")
                    oneBinHistos[sampleName]['2j2t'].Fill(1.,eventWeight)
                    regionPlots[regionPlotName].Fill(3,eventWeight)
                    if options.ttbarFlag:
                        ttbarSystHist['2j2t']["minus"].Fill(1.,eventWeight/tree.topPtReweight)
                        tempWeight = tree.topPtReweight - 1
                        ttbarSystHist['2j2t']["plus"].Fill(1.,eventWeight*((1 + (2. * tempWeight))/(1+tempWeight)))
                    for plot in plots:
                        plotList[plot]['2j2t'][regionPlotName].Fill(float(eval(plotConf.get(plot,'fillExp'))),eventWeight)
                if not len(jets) == 1: continue
                if not looseJets(tree,config.get(dataset,"runType"),jerUD,jesUD,jecUnc) == 1: continue
                cutFlows[regionPlotName].Fill(4,eventWeight)
                for plot in plots:
                    plotList[plot]['1jet'][regionPlotName].Fill(float(eval(plotConf.get(plot,'fillExp'))),eventWeight)
                if not len(tags) == 1:continue
                cutFlows[regionPlotName].Fill(5,eventWeight)
                for plot in plots:
                    plotList[plot]['btag'][regionPlotName].Fill(float(eval(plotConf.get(plot,'fillExp'))),eventWeight)
                #                if not options.cutMode == "emu" or lepton1.Pt() + lepton2.Pt() + tree.metPF2PATEt + tree.jetPF2PATPt[jets[0]] > 160.:
                if not options.cutMode == "emu" or htSys > 160.:
                    numberOfAcceptedEvents += 1
                    cutFlows[regionPlotName].Fill(6,eventWeight)
                    for plot in plots:
                        plotList[plot]['Ht'][regionPlotName].Fill(float(eval(plotConf.get(plot,'fillExp'))),eventWeight)

#        print numberOfAcceptedEvents
#        for i in range(7):
#            print cutFlows[regionPlotName].GetBinContent(i)
        if not options.skimSaveFile == "":
            outFile.Write()
            #            inputFile.Close()
            #        outFile.Write()
            outFile.Close()

    for sample in sampleNames:
        if not oneBinHistosFilled[sample]: continue
        for region in regionNames:
            print oneBinHistos[sample][region],
            print oneBinHistos[sample][region].GetBinContent(1),
            print oneBinHistos[sample][region].GetBinError(1),
            oneBinHistos[sample][region].Sumw2()
            print oneBinHistos[sample][region].GetBinError(1)
            oneBinHistos[sample][region].SaveAs(options.sysOut + options.cutMode + region + "__" + sample + histoPostfix+".root")

    #Save thye ttbar histograms, if necessary
    if options.ttbarFlag:
        for region in regionNames:
            for key in ttbarSystHist[region].keys():
                ttbarSystHist[region][key].SaveAs(options.sysOut + options.cutMode + region + "__tt__topPt__"+key+".root")
        

    if systs[2] or systs[3]:
        jecUnc.savePlot(options.outFolder + options.cutMode)
    #These could probably end up in configuration file at some point. Would
    #make my life easier in the long run.

    distributionOutFile = TFile(options.outFolder + "rawInfo" + options.cutMode + options.postfix + ".root","RECREATE")
    distributionOutFile.cd()
    for dataset in sampleNamesForPlots:
        cutFlows[dataset].Write()
    distributionOutFile.Close()
    
    distNames = ["data","other","zPlusJets","ttbar","tW"]
    colourScheme = {"other":kGreen-3,"zPlusJets":kAzure-2,"ttbar":kRed+1}
    legendOrder = ["data","tW","ttbar","zPlusJets","other"]
    legMap = {"data":("Data","p"),"tW":("tW","f"),"ttbar":('t#bar{t}','f'),"zPlusJets":('Z/#gamma*+jets',"f"),"other":("Other","f")}
    #legMap = {"other":("Other","f"),"zPlusJets":('Z/#gamma*+jets',"f"),"ttbar":('t#bar{t}','f'),"tW":("tW","f"),"data":("Data","p")}
    labelText1 = "CMS Preliminary, #sqrt{s} = 8TeV"
    modeStringForPlot = ""
    if options.cutMode == "ee":
        modeStringForPlot = "ee"
    if options.cutMode == "emu":
        modeStringForPlot = "e#mu"
    if options.cutMode == "mumu":
        modeStringForPlot = "#mu#mu"

    labelText2 = "12.2 fb^{-1}, " + modeStringForPlot + " channel"

    #create the plotting object
    cutFlowPlotter = HistogramPlotter(distNames,colourScheme,legMap,labelText1,labelText2, options.outFolder, options.postfix, options.cutMode, legendOrderList = legendOrder)
    xBinNames = ["","Lepton Selection","m_{ll}","E_{T}^{miss}","1 jet","b-tag","H_{T}"]
    cutFlowPlotter.makePlot(cutFlows,"CutFlow", "Cut Stage", xBinLabels = xBinNames)
    xBinForRegions = ["","1j1t","2j1t","2j2t"]
    cutFlowPlotter.makePlot(regionPlots,"RegionPlot","Region Plot", xBinLabels = xBinForRegions)
    
    for plot in plots:
        for stage in cutFlowStages:
            cutFlowPlotter.makePlot(plotList[plot][stage], plot, plotConf.get(plot,"xAxisLabel"),subLabel = stage)
    
    for key in foundCount.keys():
        print key, foundCount[key], foundCountWeighted[key]

def someOtherFunction():    
    gROOT.SetBatch()
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)
    gStyle.SetErrorX(0)
    setTDRStyle()

    labelcms  = TPaveText(0.12,0.88,0.5,0.94,"NDCBR")
    labelcms.SetTextAlign(12)
    labelcms.SetTextSize(0.045)
    labelcms.SetFillColor(kWhite)
    labelcms.AddText("CMS Preliminary, #sqrt{s} = 8 TeV")
    labelcms.SetBorderSize(0)

    labelcms2  = TPaveText(0.12,0.85,0.5,0.88,"NDCBR")
    labelcms2.SetTextAlign(12)
    labelcms2.SetTextSize(0.045)
    labelcms2.SetFillColor(kWhite)

    if options.cutMode == "ee":
        labelcms2.AddText("12.2 fb^{-1}, ee channel")
    if options.cutMode == "emu":
        labelcms2.AddText("12.2 fb^{-1}, e#mu channel")
    if options.cutMode == "mumu":
        labelcms2.AddText("12.2 fb^{-1}, #mu#mu channel")

#    labelcms2.AddText("12.2 fb^{-1}, ee,e#mu,#mu#mu channels")

    labelcms2.SetBorderSize(0)

    gStyle.SetPalette(1)
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetCanvasColor(kWhite)
    gStyle.SetCanvasDefH(600)
    gStyle.SetCanvasDefW(600)
    gStyle.SetLabelFont(18,"")
    
    gStyle.SetTitleXOffset(1.2)
    gStyle.SetTitleYOffset(1.2)
    
    mcDists = THStack("mcDists", "MC Distributions")
    cutFlows['other'].SetFillColor(kGreen-3)
    cutFlows['other'].SetLineColor(kBlack)
    mcDists.Add(cutFlows['other'])
    cutFlows['zPlusJets'].SetFillColor(kAzure-2)
    cutFlows['zPlusJets'].SetLineColor(kBlack)
    mcDists.Add(cutFlows['zPlusJets'])
    cutFlows['ttbar'].SetFillColor(kRed+1)
    cutFlows['ttbar'].SetLineColor(kBlack)
    mcDists.Add(cutFlows['ttbar'])
    cutFlows['tW'].SetLineColor(kBlack)
    mcDists.Add(cutFlows['tW'])
    canvy = TCanvas("Control Region Plots", "Control Region Plots")
    canvy.cd()
    leg = TLegend(0.7,0.7,0.94,0.94)
    leg.SetFillStyle(1001)
    leg.SetFillColor(kWhite)
    leg.SetBorderSize(1)
    leg.AddEntry(cutFlows['data'], 'Data', 'p')
    leg.AddEntry(cutFlows['tW'], 'tW','f')
    leg.AddEntry(cutFlows['ttbar'], 't#bar{t}','f')
    leg.AddEntry(cutFlows['zPlusJets'], 'Z/#gamma*+jets','f')
    leg.AddEntry(cutFlows['other'], 'Other','f')
    cutFlows['data'].SetMarkerStyle(20)
    cutFlows['data'].SetMarkerSize(1.2)
    cutFlows['data'].SetLineWidth(1)
    cutFlows['data'].SetMarkerColor(kBlack)
    cutFlows['data'].SetLineColor(kBlack)
    max = TMath.Max(mcDists.GetMaximum(), cutFlows['data'].GetMaximum())
    mcDists.Draw('')
    mcDists.SetMaximum(max * 1.2)
    mcDists.SetMinimum(0)
    mcDists.GetXaxis().SetBinLabel(1,"Lepton Selection")
    mcDists.GetXaxis().SetBinLabel(2,"m_{ll}")
    mcDists.GetXaxis().SetBinLabel(3,"E_{T}^{miss}")
    mcDists.GetXaxis().SetBinLabel(4,"1 jet")
    mcDists.GetXaxis().SetBinLabel(5,"b-tag")
    mcDists.GetXaxis().SetBinLabel(6,"H_{T}")

    mcDists.GetYaxis().SetTitle("events / 12.2 fb^{-1}")
    mcDists.GetYaxis().CenterTitle()
    cutFlows['data'].Draw("e, sames")
    leg.Draw()
    labelcms.Draw()
    labelcms2.Draw()

    canvy.SaveAs(options.outFolder  + "cutFlow_" + options.cutMode + options.postfix + ".png")
    canvy.SaveAs(options.outFolder  + "cutFlow_" + options.cutMode + options.postfix + ".pdf")
    canvy.SaveAs(options.outFolder  + "cutFlow_" + options.cutMode + options.postfix + ".root")
    mcDists.SetMinimum(1)
    mcDists.SetMaximum(10*max)
    canvy.SetLogy()

    canvy.SaveAs(options.outFolder  + "cutFlow_" + options.cutMode + options.postfix + "_log.png")
    canvy.SaveAs(options.outFolder  + "cutFlow_" + options.cutMode + options.postfix + "_log.pdf")
    canvy.SaveAs(options.outFolder  + "cutFlow_" + options.cutMode + options.postfix + "_log.root")

    
    canvy.SetLogy(0)

    hextra =TH1F( cutFlows['other'].Clone())
    hextra.Add(cutFlows['zPlusJets'])
    hextra.Add(cutFlows['ttbar'])
    hextra.Add(cutFlows['tW'])
    hextra.Sumw2()
    setex2 = TExec("setex2","gStyle.SetErrorX(0.5)")
    setex2.Draw()
    hextra.Sumw2()
    GE = TGraphAsymmErrors(hextra)
    
    GE.SetFillColor(28)
    GE.SetFillStyle(3018)
    GE.SetMarkerSize(0)
    GE.SetLineWidth(0)
    GE.SetLineColor(kWhite)
    leg.AddEntry(GE,"uncertainty","f")
    mcDists.SetMaximum(max*1.2)
    mcDists.SetMinimum(0)
    GE.Draw("sames, e2")
    setex1 = TExec("setex1","gStyle.SetErrorX(0)")
    setex1.Draw()
    cutFlows['data'].Draw("e, sames")

    canvy.SaveAs(options.outFolder  + "error_cutFlow_" + options.cutMode + options.postfix + ".png")
    canvy.SaveAs(options.outFolder  + "error_cutFlow_" + options.cutMode + options.postfix + ".pdf")
    canvy.SaveAs(options.outFolder  + "error_cutFlow_" + options.cutMode + options.postfix + ".root")

    mcDists.SetMaximum(max * 10)
    mcDists.SetMinimum(1)
    canvy.SetLogy()

    canvy.SaveAs(options.outFolder  + "error_cutFlow_" + options.cutMode + options.postfix + "_log.png")
    canvy.SaveAs(options.outFolder  + "error_cutFlow_" + options.cutMode + options.postfix + "_log.pdf")
    canvy.SaveAs(options.outFolder  + "error_cutFlow_" + options.cutMode + options.postfix + "_log.root")
    

#            outTree.Close()


if __name__ == "__main__":
    main()
