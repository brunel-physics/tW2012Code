from ROOT import *

class PileupReweighting:
    def __init__(self, dataFileName, mcFileName):
        self.dataFile = TFile(dataFileName,"READ")
        dataPU = self.dataFile.Get("pileup")
        
        self.mcFile = TFile(mcFileName,"READ")
        mcPU = self.mcFile.Get("pileup")

        self.pu = dataPU.Clone()
        self.pu.Scale(1.0/self.pu.Integral())
        mcPU.Scale(1.0/mcPU.Integral())
        self.pu.Divide(mcPU)
        self.number = 14

    def weight(self,npv):
        bin = self.pu.GetXaxis().FindBin(npv)
        return self.pu.GetBinContent(bin)

    def getNpv(self,weight):
        bin = 0
        for i in range(self.pu.GetNbinsX()):
            if weight == self.pu.GetBinContent(i):
                bin = i
        return bin

