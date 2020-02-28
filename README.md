# MuMuTauTauTreelizerThis tool is used to treelize MiniAOD (or reMiniAOD) samples. The output files contain vectors of different objects (eg. muons, electrons, taus) and flat branches of object counters (eg. number of vertices, event weights etc.).# Introduction for setting up the environment:$ export SCRAM_ARCH=slc6_amd64_gcc700$ cmsrel CMSSW_10_2_16$ cd CMSSW_10_2_16/src$ git cms-init$ git cms-merge-topic cms-egamma:EgammaPostRecoTools # Recipe for implanting latest Egamma ID dependence$ git cms-merge-topic -u cms-tau-pog:CMSSW_10_2_X_tau-pog_DeepTau2017v2p1_nanoAOD # Recipe for implanting lastest deepTauID$ git clone https://github.com/Fengwangdong/MuMuTauTauTreeMaker.git$ scram b -j4NOTE: When executing the scripts below, one needs to customize several items accordingly:[1] isMC: 1 = MC; 0 = data.[2] tauCluster: 0 = slimmedTaus; 1 = slimmedTausBoosted; 2 = slimmedTausMuonCleaned; 3 = slimmedTausElectronCleaned; 4 = slimmedTausNewID (deepTauID).[3] inputFiles: path + input file (file by file, wildcards do not work).# Run the diMuon + ditau treelizer:$ cd MuMuTauTauTreeMaker/MuTauTreelizer/test$ cmsRun runDiMuDiTau_cfg.py isMC=1(0) tauCluster=0(1/2/3/4) inputFiles=/PATH/file1.root inputFiles=/PATH/file2.root# Run the ZMuMu inclusive treelizer:$ cd MuMuTauTauTreeMaker/MuTauTreelizer/test$ cmsRun runZMuMuInclusive_cfg.py isMC=1(0) tauCluster=0(1/2/3/4) inputFiles=/PATH/file1.root inputFiles=/PATH/file2.root# Run the ZMuTau (tau_mu + tau_h) treelizer:$ cd MuMuTauTauTreeMaker/MuTauTreelizer/test$ cmsRun runZTauMuTauHad_cfg.py isMC=1(0) tauCluster=0(1/2/3/4) inputFiles=/PATH/file1.root inputFiles=/PATH/file2.root