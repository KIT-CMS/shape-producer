# -*- coding: utf-8 -*-

from cutstring import *
from estimation_methods import EstimationMethod, SStoOSEstimationMethod, ABCDEstimationMethod
from estimation_methods_2016 import DataEstimation as DataEstimation2016
from estimation_methods_2016 import WEstimationWithQCD as WEstimationWithQCD2016
from estimation_methods_2016 import QCDEstimationWithW as QCDEstimationWithW2016
from systematics import *
from histogram import *
from systematic_variations import *
from era import log_query


class DataEstimation(DataEstimation2016):
    #    def get_cuts(self):
    #        return Cuts(Cut("run <= 300676", "rereco_equivalent"))

    pass


class HTTEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="HTT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer17MiniAOD")

    def get_weights(self):
        return Weights(
            Weight("((gen_match_2 == 5)*0.95 + (gen_match_2 != 5))",
                   "hadronic_tau_sf"),
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
            Weight("puweight", "puweight"),
            #Weight("eventWeight", "eventWeight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "(^GluGluHToTauTau.*125.*|^VBFHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

    def get_cuts(self):
        if self.channel.name == "mt":
            return Cuts(
                Cut("(trg_singlemuon==1 && pt_1>29 && pt_2>30)",
                    "trg_singlemuon"))
        else:
            return Cuts()


class ggHEstimation(HTTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="ggH",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer17MiniAOD")

    def get_files(self):
        query = {
            "process": "^GluGluHToTauTau.*125.*",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

    def get_cuts(self):
        if self.channel.name == "mt":
            return Cuts(
                Cut("(trg_singlemuon==1 && pt_1>29 && pt_2>30)",
                    "trg_singlemuon"))
        else:
            return Cuts()


class qqHEstimation(HTTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="qqH",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer17MiniAOD")

    def get_files(self):
        query = {
            "process": "^VBFHToTauTau.*125.*",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

    def get_cuts(self):
        if self.channel.name == "mt":
            return Cuts(
                Cut("(trg_singlemuon==1 && pt_1>29 && pt_2>30)",
                    "trg_singlemuon"))
        else:
            return Cuts()


class VHEstimation(HTTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(HTTEstimation, self).__init__(
            name="VH",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer17MiniAOD")

    def get_files(self):
        query = {
            "process": "(^W(minus|plus)HToTauTau.*125.*|^ZHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }

        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class TTEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TT",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAOD")

    def get_weights(self):
        return Weights(

            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
            Weight("0.5", "tt_stitching_weight"),

            # Weights for corrections
            Weight("topPtReweightWeight", "topPtReweightWeight"),
            Weight("((gen_match_2 == 5)*0.95 + (gen_match_2 != 5))",
                   "hadronic_tau_sf"),
            #Weight("eventWeight", "eventWeight"),
            Weight("puweight", "puweight"),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^TT",
            "data": False,
            "campaign": self._mc_campaign,
            #"version": "v1" # to be used if only one inclusive sample is desired
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class EWKEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(EWKEstimation, self).__init__(
            name="EWK",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAOD")

    def get_files(self):
        query = {
            "process": "^EWKZ",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

    def get_weights(self):
        return Weights(
            Weight("((gen_match_2 == 5)*0.95 + (gen_match_2 != 5))",
                   "hadronic_tau_sf"),
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
            Weight("puweight", "puweight"),
            #Weight("eventWeight", "eventWeight"),
            self.era.lumi_weight)


class TTLEstimationMT(TTEstimation):
    # L refering to a prompt t-quark to lepton decay as opposed to t->tau->lepton (important for embedded events)
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTL",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAOD")

    def get_cuts(self):
        return Cuts(
            Cut("gen_match_2==5", "genmatch"),
            Cut("!((gen_match_1==4)&&(gen_match_2==5))",
                "ttbar->tau tau veto for embedded events"))


class TTLEstimationET(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTL",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAOD")

    def get_cuts(self):
        return Cuts(
            Cut("gen_match_2==5", "genmatch"),
            Cut("!((gen_match_1==3)&&(gen_match_2==5))",
                "ttbar->tau tau veto for embedded events"))


class TTLEstimationTT(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTL",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAOD")

    def get_cuts(self):
        return Cuts(
            Cut("((gen_match_1 == 5) && (gen_match_2 == 5))*!((gen_match_1 == 5) && (gen_match_2 == 5))",
                "Empty Process")
        )  # All ttbar->real tau events are vetoed for embedded events

class TTLEstimationEM(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTL",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAOD")

    def get_cuts(self):
        return Cuts(
            Cut("!((gen_match_1 == 3) && (gen_match_2 == 4))",
                "ttbar->tau tau veto for embedded events")
        )

class QCDEstimationET(SStoOSEstimationMethod):
    def __init__(self,
                 era,
                 directory,
                 channel,
                 bg_processes,
                 data_process,
                 friend_directory=None,
                 extrapolation_factor=1.0):
        super(QCDEstimationET, self).__init__(
            name="QCD",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            bg_processes=bg_processes,
            data_process=data_process,
            extrapolation_factor=extrapolation_factor)


class QCDEstimation_SStoOS_MTETEM(SStoOSEstimationMethod):
    def __init__(self,
                 era,
                 directory,
                 channel,
                 bg_processes,
                 data_process,
                 extrapolation_factor=1.0):
        super(QCDEstimation_SStoOS_MTETEM, self).__init__(
            name="QCD",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            bg_processes=bg_processes,
            data_process=data_process,
            extrapolation_factor=extrapolation_factor)


class QCDEstimation_ABCD_TT_ISO2(ABCDEstimationMethod):
    def __init__(self, era, directory, channel, bg_processes, data_process):
        super(QCDEstimation_ABCD_TT_ISO2, self).__init__(
            name="QCD",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            bg_processes=bg_processes,
            data_process=data_process,
            AC_cut_names=
            [  # cuts applied in AC, which should be removed in the BD control regions
                "tau_2_iso"
            ],
            BD_cuts=[  # cuts to be applied instead of cuts removed above
                Cut("byTightIsolationMVArun2v1DBoldDMwLT_2<0.5", "tau_2_iso"),
                #Cut("byMediumIsolationMVArun2v1DBoldDMwLT_2<0.5", "tau_2_iso"),
                Cut("byLooseIsolationMVArun2v1DBoldDMwLT_2>0.5",
                    "tau_2_iso_loose")
            ],
            AB_cut_names=
            [  # cuts applied in AB, which should be removed in the CD control regions
                "os"
            ],
            CD_cuts=[  # cuts to be applied instead of cuts removed above
                Cut("q_1*q_2>0", "ss")
            ])


class QCDEstimation_ABCD_TT_ISO2_TRANSPOSED(ABCDEstimationMethod):
    def __init__(self, era, directory, channel, bg_processes, data_process):
        super(QCDEstimation_ABCD_TT_ISO2_TRANSPOSED, self).__init__(
            name="QCD",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            bg_processes=bg_processes,
            data_process=data_process,
            AB_cut_names=
            [  # cuts applied in AB, which should be removed in the CD control regions
                "tau_2_iso"
            ],
            CD_cuts=[  # cuts to be applied instead of cuts removed above
                Cut("byTightIsolationMVArun2v1DBoldDMwLT_2<0.5", "tau_2_iso"),
                #Cut("byMediumIsolationMVArun2v1DBoldDMwLT_2<0.5", "tau_2_iso"),
                Cut("byLooseIsolationMVArun2v1DBoldDMwLT_2>0.5",
                    "tau_2_iso_loose")
            ],
            AC_cut_names=
            [  # cuts applied in AC, which should be removed in the BD control regions
                "os"
            ],
            BD_cuts=[  # cuts to be applied instead of cuts removed above
                Cut("q_1*q_2>0", "ss")
            ])


class QCDEstimation_ABCD_TT_ISO1(ABCDEstimationMethod):
    def __init__(self, era, directory, channel, bg_processes, data_process):
        super(QCDEstimation_ABCD_TT_ISO1, self).__init__(
            name="QCD",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            bg_processes=bg_processes,
            data_process=data_process,
            AC_cut_names=
            [  # cuts applied in AC, which should be removed in the BD control regions
                "tau_1_iso"
            ],
            BD_cuts=[  # cuts to be applied instead of cuts removed above
                Cut("byTightIsolationMVArun2v1DBoldDMwLT_1<0.5", "tau_1_iso"),
                #Cut("byMediumIsolationMVArun2v1DBoldDMwLT_1<0.5", "tau_1_iso"),
                Cut("byLooseIsolationMVArun2v1DBoldDMwLT_1>0.5",
                    "tau_1_iso_loose")
            ],
            AB_cut_names=
            [  # cuts applied in AB, which should be removed in the CD control regions
                "os"
            ],
            CD_cuts=[  # cuts to be applied instead of cuts removed above
                Cut("q_1*q_2>0", "ss")
            ])


class VVEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(VVEstimation, self).__init__(
            name="VV",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAOD")

    def get_files(self):
        query = {
            "process": "(WW|ZZ|WZ)",  # Query for Di-Boson samples
            "data": False,
            "generator": "^pythia8",
            "campaign": self._mc_campaign
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)

        query = {
            "process": "STtW",  # Query for Single-Top samples
            "data": False,
            "generator": "powheg\-pythia8",
            "campaign": self._mc_campaign
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)

        log_query(self.name, query, files)
        return self.artus_file_names(files)

    def get_weights(self):
        return Weights(
            Weight("((gen_match_2 == 5)*0.95 + (gen_match_2 != 5))",
                   "hadronic_tau_sf"),
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
            Weight("puweight", "puweight"),
            #Weight("eventWeight", "eventWeight"),
            self.era.lumi_weight)

    def get_cuts(self):
        if self.channel.name == "mt":
            return Cuts(
                Cut("(trg_singlemuon==1 && pt_1>29 && pt_2>30)",
                    "trg_singlemuon"))
        else:
            return Cuts()


class DYJetsToLLEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(DYJetsToLLEstimation, self).__init__(
            name="DYJetsToLL",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAOD")

    def get_weights(self):
        weightstring = Weights(
            Weight("((((genbosonmass >= 50.0 && (npartons == 0 || npartons >= 4))*2.455936181) + ((genbosonmass >= 50.0 && npartons == 1)*0.5608870881) + ((genbosonmass >= 50.0 && npartons == 2)*0.5745263806) + ((genbosonmass >= 50.0 && npartons == 3)*0.617450628))+ ((genbosonmass < 50.0)*numberGeneratedEventsWeight*crossSectionPerEventWeight))",
                "stitching_weight"),
            # Weights for corrections
            Weight("gen_match_2 == 5","hadronic_tau_sf"),
            Weight("idisoweight_1*triggerweight_1", "muon_sf"),
            self.era.lumi_weight)
        if self.channel.name == "em":
            weightstring = Weights(
                Weight( "((((genbosonmass >= 50.0 && (npartons == 0 || npartons >= 4))*2.455936181) + ((genbosonmass >= 50.0 && npartons == 1)*0.5608870881) + ((genbosonmass >= 50.0 && npartons == 2)*0.5745263806) + ((genbosonmass >= 50.0 && npartons == 3)*0.617450628))+ ((genbosonmass < 50.0)*numberGeneratedEventsWeight*crossSectionPerEventWeight))",
                    "stitching_weight"),
                # Weights for corrections
                Weight("gen_match_2 == 5","hadronic_tau_sf"),
                Weight("idisoweight_1*idisoweight_2", "muon_sf"),
                self.era.lumi_weight)
        return weightstring
            

    def get_files(self):
        query = {
            "process": "(DYJetsToLL_M10to50|DY.?JetsToLL_M50)",
            #"process": "(DYJetsToLL_M10to50|DYJetsToLL_M50)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            #"version": "v1" # to be used if only one inclusive sample is desired
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ZTTEstimation(DYJetsToLLEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(DYJetsToLLEstimation, self).__init__(
            name="ZTT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAOD")

    def get_cuts(self):
        if self.channel.name in {"mt", "et"}:
                ztt_genmatch_cut = Cut("gen_match_2==5", "ztt_genmatch")
        elif self.channel.name == "tt":
            ztt_genmatch_cut = Cut("(gen_match_1==5) && (gen_match_2==5)",
                                    "ztt_genmatch")
        elif self.channel.name == "em":
            ztt_genmatch_cut = Cut("(gen_match_1>2) && (gen_match_2>3)",
                                    "ztt_genmatch")
        return Cuts(ztt_genmatch_cut)


class ZLEstimationMT(ZTTEstimation):
    def get_cuts(self):
        return Cuts(Cut("gen_match_2<5", "zl_genmatch_mt"))


class ZTTEstimationTT(ZTTEstimation):
    def get_cuts(self):
        return Cuts(
            Cut("(gen_match_1==5 && gen_match_2==5)", "ztt_genmatch_tt"))

class ZTTEstimationEM(ZTTEstimation):
    def get_cuts(self):
        return Cuts(
            Cut("(gen_match_1==3 && gen_match_2==4)", "ztt_genmatch_em"))


class ZJEstimationMT(ZTTEstimation):
    def get_cuts(self):
        return Cuts(Cut("gen_match_2==6", "zj_genmatch_mt"))


class ZLEstimationMTSM(ZLEstimationMT):
    def get_weights(self):
        ztt_weights = super(ZLEstimationMT, self).get_weights()
        return ztt_weights + Weights(
            Weight(
                "(((decayMode_2 == 0)*0.75) + ((decayMode_2 == 1 || decayMode_2 == 2)*1.0) + ((decayMode_2 == 10)*1.0))","decay_mode_reweight"))


class ZLLEstimation(DYJetsToLLEstimation):
    def __init__(self, era, directory, channel):
        super(DYJetsToLLEstimation, self).__init__(
            name="ZLL",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAOD")

    def get_cuts(self):
        zll_genmatch_cut = Cut("1 == 1", "zll_genmatch")
        if self.channel.name in ["mt", "et"]:
            zll_genmatch_cut = Cut("gen_match_2!=5", "zll_genmatch")
        elif self.channel.name == "tt":
            zll_genmatch_cut = Cut("(gen_match_1!=5) || (gen_match_2!=5)",
                                   "zll_genmatch")
        elif self.channel.name == "em":
            zll_genmatch_cut = Cut("(gen_match_1<3) || (gen_match_2<4)",
                                   "zll_genmatch")
        return Cuts(zll_genmatch_cut)


class ZTTEmbeddedEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(ZTTEmbeddedEstimation, self).__init__(
            name="ZTT",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign=None)

    def get_weights(self):
        if self.channel.name in {"mt", "et"}:
            return Weights(
                Weight("generatorWeight", "simulation_sf"),
                Weight("muonEffTrgWeight", "scale_factor"),
                Weight("idweight_1*triggerweight_1*isoweight_1", "leptopn_sf"))
        elif self.channel.name == "tt":
            return Weights(
                Weight("generatorWeight", "simulation_sf"),
                Weight("muonEffTrgWeight", "scale_factor"),
                Weight("(idweight_1*isoweight_1*idweight_2*isoweight_2*(triggerweight_1*(triggerweight_1<=1.8)+(triggerweight_1>1.8))*(triggerweight_2*(triggerweight_2<=1.8)+(triggerweight_2>1.8))", "leptopn_sf"))
        elif self.channel.name == "em":
            return Weights(
                Weight("generatorWeight", "simulation_sf"),
                Weight("muonEffTrgWeight", "scale_factor"),
                # no trigger sf yet
                Weight("idweight_1*isoweight_1*idweight_2*isoweight_2", "leptopn_sf"))

    def get_files(self):
        query = {"process": "Embedding2017(B|C|D|E|F)", "embedded": True}
        #query = {"process" : "Embedding2017(B|C|D|E|F)", "embedded" : True}
        if self.channel.name == "mt":
            query["campaign"] = "MuTauFinalState"
            query["scenario"] = ".*v2"
        elif self.channel.name == "et":
            query["campaign"] = "ElTauFinalState"
            query["scenario"] = ".*v2"
        elif self.channel.name == "tt":
            query["campaign"] = "TauTauFinalState"
            query["scenario"] = ".*v2"
        elif self.channel.name == "em":
            query["campaign"] = "ElMuFinalState"
            query["scenario"] = ".*v2"
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

    def get_cuts(self):
        if self.channel.name in ["mt", "et"]:
            ztt_genmatch_cut = Cut("gen_match_2==5","ztt_genmatch")
        elif self.channel.name == "tt":
            ztt_genmatch_cut = Cut("(gen_match_1==5) && (gen_match_2==5)",
                                   "ztt_genmatch")
        elif self.channel.name == "em":
            ztt_genmatch_cut = Cut("(gen_match_1>2) && (gen_match_2>3)",
                                   "ztt_genmatch")
        return Cuts(ztt_genmatch_cut)


class ZTTEmbeddingEstimation_ScaledToMC(EstimationMethod):
    def __init__(self, era, directory, channel, embedding_process,
                 ttbar_tautau_mc_process, z_tautau_mc_process):
        super(ZTTEmbeddingEstimation_ScaledToMC, self).__init__(
            name="ZTT",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign=None)
        self._embedding_process = copy.deepcopy(embedding_process)
        self._ttbar_tautau_mc_process = copy.deepcopy(ttbar_tautau_mc_process)
        self._z_tautau_mc_process = copy.deepcopy(z_tautau_mc_process)

    def create_root_objects(self, systematic):
        yield_category = copy.deepcopy(systematic.category)
        yield_category._variable = None

        shape_category = copy.deepcopy(systematic.category)
        shape_category._name += "_unscaled"

        root_objects = []
        systematic._embedding_systematics = []
        shape_systematic = Systematic(
            category=shape_category,
            process=self._embedding_process,
            analysis=systematic.analysis,
            era=self.era,
            variation=systematic.variation,
            mass=125)
        systematic._embedding_systematics.append(shape_systematic)
        shape_systematic.create_root_objects()
        root_objects += shape_systematic.root_objects

        for process in [
                self._embedding_process, self._ttbar_tautau_mc_process,
                self._z_tautau_mc_process
        ]:
            s = Systematic(
                category=yield_category,
                process=process,
                analysis=systematic.analysis,
                era=self.era,
                variation=systematic.variation,
                mass=125)
            systematic._embedding_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects
        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_embedding_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _embedding_systematics needed for embedding scaled to MC estimation.",
                systematic.name)
            raise Exception

        for s in systematic._embedding_systematics:
            s.do_estimation()

        shapes = [s.shape for s in systematic._embedding_systematics]

        # embedding shape
        embedding_shape = shapes[0]

        # scale factor = MC(TTT + ZTT) yield / embedding yield
        sf = (shapes[2].result + shapes[3].result) / shapes[1].result
        print "Scale factor", sf

        # scaling shape
        embedding_shape.result.Scale(sf)

        embedding_shape.name = systematic.name
        return embedding_shape


class WEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(WEstimation, self).__init__(
            name="W",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIISummer17MiniAOD")

    def get_weights(self):
        return Weights(

            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            #Weight("numberGeneratedEventsWeight","numberGeneratedEventsWeight"), # to be used only for one inclusive sample
            #Weight("crossSectionPerEventWeight","crossSectionPerEventWeight"), # to be used only for one inclusive sample
            Weight(
                "(((npartons == 0 || npartons >= 5)*2.36006270e-3) + ((npartons == 1)*2.34817764e-4) + ((npartons == 2)*1.31144867e-4) + ((npartons == 3)*1.39177532e-4) + ((npartons == 4)*6.46064804e-5))",
                "wj_stitching_weight"),

            # Weights for corrections
            # Weight("1.5*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5))",
            #        "hadronic_tau_sf"),
            Weight("((gen_match_2 == 5)*0.95 + (gen_match_2 != 5))",
                   "hadronic_tau_sf"),
            Weight("puweight", "puweight"),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "W.?JetsToLNu",
            #"process": "WJetsToLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

    def get_cuts(self):
        if self.channel.name == "mt":
            return Cuts(
                Cut("(trg_singlemuon==1 && pt_1>30 && pt_2>30)",
                    "trg_singlemuon"))
        else:
            return Cuts()


class TTTEstimation(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTT",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAOD")

    def get_cuts(self):
        return Cuts(
            Cut("(gen_match_1 > 2 && gen_match_1 < 6) && (gen_match_1 > 2 && gen_match_2 < 6)",
                "gen_match_genuine_taus"))


class TTTEstimationMT(TTTEstimation):
    def get_cuts(self):
        return Cuts(Cut("gen_match_2==5 && gen_match_1 == 4", "ttt_genmatch_mt"))


class TTJEstimationMT(TTTEstimation):
    def get_cuts(self):
        return Cuts(Cut("gen_match_2!=3", "ttj_genmatch_mt"))

class TTTEstimationEM(TTTEstimationMT):
    def get_cuts(self):
        return Cuts(Cut("gen_match_1 == 3 && gen_match_2 == 4", "ttt_genmatch_em"))
    pass

class TTJEstimationEM(TTTEstimation):
    def get_cuts(self):
        return Cuts(Cut("gen_match_2!=3", "ttj_genmatch_mt"))

class TTJEstimation(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTJ",
            folder="nominal",
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAOD")

    def get_cuts(self):
        return Cuts(
            Cut("!(gen_match_1 > 2 && gen_match_1 < 6) && (gen_match_1 > 2 && gen_match_2 < 6)",
                "gen_match_genuine_taus"))


class QCDEstimationMT(QCDEstimationET):
    pass

class QCDEstimationEM(QCDEstimationET):
    pass



class ZLEstimationETSM(ZLEstimationMT):
    def get_weights(self):
        ztt_weights = super(ZLEstimationMT, self).get_weights()
        return ztt_weights + Weights(
            # Weight(
            #     "1.5*(((decayMode_2 == 0)*0.98) + ((decayMode_2 == 1 || decayMode_2 == 2)*1.2) + ((decayMode_2 == 10)*1.0))",
            #     "decay_mode_reweight"))
            Weight(
                "(((decayMode_2 == 0)*0.98) + ((decayMode_2 == 1 || decayMode_2 == 2)*1.2) + ((decayMode_2 == 10)*1.0))",
                "decay_mode_reweight"))


# et is equivalent to mt
class ZJEstimationET(ZJEstimationMT):
    pass


class TTTEstimationET(TTTEstimationMT):
    pass


class TTJEstimationET(TTJEstimationMT):
    pass


class ZJEstimationTT(ZJEstimationMT):
    def get_cuts(self):
        return Cuts(
            Cut("(gen_match_2 == 6 || gen_match_1 == 6)", "zj_genmatch_tt"))


class ZLEstimationTT(ZLEstimationMT):
    def get_cuts(self):
        return Cuts(
            Cut("(gen_match_1<6 && gen_match_2<6 &&! (gen_match_1==5 && gen_match_2==5))",
                "zl_genmatch_tt"))

class ZJEstimationEM(ZJEstimationMT):
    pass


class ZLEstimationEM(ZTTEstimation):
    def get_cuts(self):
        return Cuts(Cut("gen_match_2<4 || gen_match_1<3", "zl_genmatch_em"))



class TTTEstimationTT(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAOD")

    def get_cuts(self):
        return Cuts(
            Cut("((gen_match_1 == 5) && (gen_match_2 == 5))",
                "select ttbar->tau tau events"))


class TTJEstimationTT(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAOD")

    def get_cuts(self):
        return Cuts(
            Cut("gen_match_2==5", "genmatch"),
            Cut("!((gen_match_1==5)&&(gen_match_2==5))",
                "ttbar->tau tau veto for embedded events"))


class QCDEstimationTT(ABCDEstimationMethod):
    def __init__(self,
                 era,
                 directory,
                 channel,
                 bg_processes,
                 data_process,
                 friend_directory=None):
        super(QCDEstimationTT, self).__init__(
            name="QCD",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            bg_processes=bg_processes,
            data_process=data_process,
            AC_cut_names=
            [  # cuts to be removed to include region for shape derivation
                "tau_2_iso"
            ],
            BD_cuts=
            [  # cuts to be applied to restrict to region for shape derivation
                Cut("byTightIsolationMVArun2v1DBoldDMwLT_2<0.5", "tau_2_iso"),
                Cut("byLooseIsolationMVArun2v1DBoldDMwLT_2>0.5",
                    "tau_2_iso_loose")
            ],
            AB_cut_names=
            [  # cuts to be removed to include region for the determination of the extrapolation derivation
                "os"
            ],
            CD_cuts=
            [  # cuts to be applied to restrict to region for the determination of the extrapolation derivation
                Cut("q_1*q_2>0", "ss")
            ])


class WEstimationWithQCD(WEstimationWithQCD2016):
    pass


class QCDEstimationWithW(QCDEstimationWithW2016):
    pass


class TTTTEstimationMT(TTEstimation):
    # true tt->tautau
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTTT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAOD")

    def get_cuts(self):
        return Cuts(
            Cut("((gen_match_1 == 4) && (gen_match_2 == 5))",
                "select ttbar->tau tau events"))


class TTTTEstimationET(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTTT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAOD")

    def get_cuts(self):
        return Cuts(
            Cut("((gen_match_1 == 3) && (gen_match_2 == 5))",
                "select ttbar->tau tau events"))


class TTTEstimationTT(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(TTEstimation, self).__init__(
            name="TTT",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAOD")

    def get_cuts(self):
        return Cuts(
            Cut("((gen_match_1 == 5) && (gen_match_2 == 5))",
                "select ttbar->tau tau events"))
