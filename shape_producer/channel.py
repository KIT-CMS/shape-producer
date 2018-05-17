# -*- coding: utf-8 -*-

from cutstring import Cuts, Cut
import logging
logger = logging.getLogger(__name__)
"""
"""


class Channel(object):
    @property
    def cuts(self):
        return self._cuts

    @property
    def name(self):
        return self._name


class EMSM(Channel):
    def __init__(self):
        self._name = "em"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("iso_1<0.15", "ele_iso"), Cut("iso_2<0.2", "muon_iso"),
            Cut("nbtag==0", "bveto"), Cut("diLepMetMt<60.0", "diLepMetMt"),
            Cut("pZetaMissVis>-35.0", "pzeta"), Cut("q_1*q_2<0", "os"),
            Cut("(trg_electronmuon==1 && pt_1>26 && pt_2>25)",
                "trg_electronmuon"))


class EESM(Channel):
    def __init__(self):
        self._name = "ee"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("iso_1<0.1 && iso_2<0.1", "ele_iso"), Cut("q_1*q_2<0", "os"),
            Cut("(trg_singleelectron==1 && pt_1>26 && pt_2>26)",
                "trg_singleelectron"))


class MMSM(Channel):
    def __init__(self):
        self._name = "mm"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("iso_1<0.15 && iso_2<0.15", "muon_iso"), Cut(
                "q_1*q_2<0", "os"),
            Cut("(trg_singlemuon==1 && pt_1>25 && pt_2>25)", "trg_singlemuon"))


class MT(Channel):
    def __init__(self):
        self._name = "mt"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonTight3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_2>0.5", "tau_iso"),
            #Cut("byTightCombinedIsolationDeltaBetaCorr3Hits_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "muon_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("trg_singlemuon==1", "trg_singlemuon"))


class MTSM(Channel):
    def __init__(self):
        self._name = "mt"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonTight3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "muon_iso"), Cut("q_1*q_2<0", "os"),
            Cut("mt_1<50", "m_t"),
            Cut("((trg_singlemuon==1 && pt_1>20 && pt_2>30))",
                "trg_singlemuoncross"))


class MTSM_17(Channel):
    def __init__(self):
        self._name = "mt"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonTight3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "muon_iso"), Cut("q_1*q_2<0", "os"),
            Cut("mt_1<50", "m_t"),
            Cut("pt_1>30 && pt_2>30", "pt_cut"),
            Cut("trg_singlemuon==1","trg_singlemuoncross")
        )  # in samples trg_singlemuon is "HLT_IsoMu27_v:29.0"


class ET(Channel):
    def __init__(self):
        self._name = "et"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronTightMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_2>0.5", "tau_iso"),
            #Cut("byTightCombinedIsolationDeltaBetaCorr3Hits_2>0.5", "tau_iso"),
            Cut("iso_1<0.1", "ele_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("trg_singleelectron==1", "trg_singleelectron"))


class ETSM(Channel):
    def __init__(self):
        self._name = "et"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronTightMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_2>0.5", "tau_iso"),
            Cut("iso_1<0.1", "ele_iso"), Cut("q_1*q_2<0", "os"),
            Cut("mt_1<50", "m_t"),
            Cut("(trg_singleelectron==1 && pt_1>26 && pt_2>30)",
                "trg_singleelectron"))


class ETSM_17(Channel):
    def __init__(self):
        self._name = "et"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronTightMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_2>0.5", "tau_iso"),
            Cut("iso_1<0.1", "ele_iso"), Cut("q_1*q_2<0", "os"),
            Cut("mt_1<50", "m_t"),
            Cut("(trg_singleelectron==1 && pt_1>37 && pt_2>30)",
                "trg_singleelectron"))


class TT(Channel):
    def __init__(self):
        self._name = "tt"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_1>0.5", "tau_1_iso"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_2>0.5", "tau_2_iso"),
            Cut("q_1*q_2<0", "os"), 
            Cut("trg_doubletau==1", "trg_doubletau"))


class TTSM(Channel):
    def __init__(self):
        self._name = "tt"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_1>0.5", "tau_1_iso"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_2>0.5", "tau_2_iso"),
            Cut("q_1*q_2<0", "os"), Cut("pt_tt>50", "pt_h"),
            Cut("(pt_1>50 && pt_2>40)", "trg_doubletau"))


class TTSM_17(Channel):
    def __init__(self):
        self._name = "tt"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_1>0.5", "tau_1_iso"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_2>0.5", "tau_2_iso"),
            Cut("q_1*q_2<0", "os"), 
            Cut("pt_tt>50", "pt_h"),
            Cut("(pt_1>50 && pt_2>40)", "trg_doubletau"))


class EM(Channel):
    def __init__(self):
        self._name = "em"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("iso_1<0.15", "ele_iso"),
            Cut("iso_2<0.2", "muon_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("trg_muonelectron_lowptmu==1", "trg_muonelectron"),
            Cut("pt_1>20 && pt_2 > 20", "pt_cut"))


class PU(Channel):
    def __init__(self):
        self._name = "pu"
        self._cuts = Cuts()


# collection of channels an analysis can be ran on
class Channels(object):
    def __init__(self, name):
        self._name = name
        self._channels = []

    def add(self, channel):
        self._channels.append(channel)

    @property
    def name(self):
        return self._name
