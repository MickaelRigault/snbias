#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This modules manages the I/O to access the data"""

import os
import sys
import numpy as np
from astropy import constants,table
from astrobject.utils.tools import load_pkl,dump_pkl

#########################################
#                                       #
# Created by: M. Rigault                #
# Purpose: Build basic functionnality   #
#          to study SN bias             #
#                                       #
#########################################

__all__ = ["get_local_hostdata"]

DATAPATH = os.path.dirname(__file__)+"/data/"


################################
#                              #
#  General IO                  #
#                              #
################################
def format_snname(name):
    """
    This methods enable to get a commun formating of the snnames:
    - SN+YEAR+upper cap for single letter or lower cap for multiple letter
    - if SNfactory: SNF+YYYYMMDD-NNN
    - if PTF: PTFYY+lower cap for multiple letter
    """
    if (name.lower()).startswith("sn") is False:
        return name
    
    if len(name)==7 :
        nameDB = name.upper()
    else:
        nameDB = name.replace('sn','SN')
    return nameDB

################################
#                              #
#  Cepheid Data                #
#                              #
################################
CEPHEIDSNE_RIESS16=DATAPATH+"/Riess16/cepheid_sne_distances.dat"
SNEIA_RIESS16=DATAPATH+"/Riess16/scolnic_snedata.dat"
_FIXED_SNEIA_RIESS16=DATAPATH+"/Riess16/scolnic_snedata_rest.dat"
def get_riess16_cepheidsn():
    """ Astropy Table of the Table 5 of Riess et al. 2016. Cepheid SN data information """
    return table.Table.read(CEPHEIDSNE_RIESS16, format="csv", fast_reader=False, comment="#")

def get_riess16_sndata():
    """ SN Ia data used by Riess et al. 2016, from D. Scolnic (http://kicp.uchicago.edu/~dscolnic/Supercal/supercal_vH0.fitres; 7th of April 2016)
    Base on "supercal.fitres" header:
    abs(x1)<3, abs(c)<0.3, fitprob>0.01, pkmjderr<2.0, x1err<1.0, and some reasonable 3-4 sigma cut on Hubble Residual (mures) outliers

    Do not use stone and colfax
    """
    maintable = table.Table.read(SNEIA_RIESS16, format="ascii", guess=False)
    endtable = table.Table.read(_FIXED_SNEIA_RIESS16, format="ascii", guess=False)
    
    return table.vstack([maintable,endtable])

################################
#                              #
#  SN Data                     #
#                              #
################################
SNF_DATA_META      = DATAPATH+"/SNf_SNeIa/META.pkl"
SNF_DATA_HUBBLIZER = DATAPATH+"/SNF-0203-ALLEG2a_SNeIa_hubble.pkl"
SNF_DATA_PHRENO    = DATAPATH+"/phrenology_2015_06_15_ALLAIRE_AtMax.pkl"

def get_snf_data(kind):
    """
    Access the SNfactory data. kind is phreno/hubblizer/idr.
    A formatted dictonnary is returned
    """
    # -------------------
    # - Hubblizer Data
    # -------------------
    if kind=="hubblizer":
        dico =  load_pkl(SNF_DATA_HUBBLIZER)
        data = {}
        for name,d in dico.items():
            if "hubblizer.mBfit" not in d.keys():
                continue
            data[name]={"zcmb":d['host.zcmb'],
                        "mu":d['hubblizer.mBfit'],
                        "mu.err":d['hubblizer.mBfit.err'],
                        "HR":d['hubblizer.dmfit_corr'],
                        "HR_err":d['hubblizer.dmfit_corr.err'],
                        "oHR":d['hubblizer.dmfit_orig'],
                        "oHR_err":d['hubblizer.dmfit_orig.err'],
                        "x1":d["salt2.X1"],
                        "color":d["salt2.Color"]}
        return data
    
    # -------------------
    # - Base IDR
    # -------------------
    elif kind in ["meta","idr"]:
        data = load_pkl(SNF_DATA_META)
        # -----------
        # - Formating
        data_={}
        for name,d in data.items():
            data_[name] = {}
            for k in d.keys():
                if len(k.split("."))>1:
                    data_[name]["".join(k.split(".")[1:]).replace('err',"_err").lower()] = d[k]
                else:
                    data_[name][k.lower()] = d[k]

        del data
        return data_
    # -------------------
    # - Phrenology
    # -------------------
    elif kind in ["phreno","phrenology"]:
        data=load_pkl(SNF_DATA_PHRENO)
        data_ = {}
        for name,d in data.items():
            data_[name]={}
            for k in d.keys():
                if k is None: continue # Who did this dico ??!!
                data_[name][k.split("phrenology.")[-1].replace(".err","_err")] = d[k]
        del data
        return data_
    else:
        raise ValueError("Unknown snf data 'kind' %s"%kind)

################################
#                              #
#  Host Data                   #
#                              #
################################
# ======================= #
# =  Local SNf Data     = #
# ======================= #
LOCAL_SNF_DATA=DATAPATH+"/SNf_localhost/localhost_idr.pkl"

def get_local_hostdata():
    """ This is the data created by the SNfactory local galaxy stellar free fits"""
    return load_pkl(LOCAL_SNF_DATA)

# ======================= #
# =  Global Data        = #
# ======================= #

def get_list_globalhost_sources():
    return [_source.split("gethost_")[-1].split("_data")[0]
                      for _source in dir(sys.modules[__name__])
                      if _source.startswith("gethost")]

def get_global_hostdata(sourcename,**kwargs):
    """
    """
    sourcename = sourcename.lower()
    if sourcename.lower() not in get_list_globalhost_sources():
        raise ValueError("unknown global host source '%s' "%sourcename+"\n"+"known sources: "+", ".join(GLOBAL_HOSTSOURCES))
    return eval("gethost_%s_data(**kwargs)"%sourcename)

# -------------------- #
# -    INPUT FILES   - #
# -------------------- #
_fileBetoule14      = DATAPATH+"jla_lcparams.txt"
_fileAsiogoSrc      = DATAPATH+"AsiagoSource_asu.fit"
_fileAsiogoDico     = DATAPATH+"Asiago_dicoBasic.pkl"
_fileRiessH09morpho = DATAPATH+"Riess_morphology_classification.dat"
_fileSilverman12Src = DATAPATH+"Silverman_2012_SNclassification_HostType.fits"
_fileSilverman12    = DATAPATH+"Silverman_2012_SNclassification_HostType.pkl"
_fileChildress13    = DATAPATH+"MJC_compile_SNdata.pkl"
_fileChildressSNf   = DATAPATH+"Childress_2013_snf_final_host_metadata.pkl"
_fileNeill09Scr1    = DATAPATH+"Neill09_table1.dat"
_fileNeill09Scr2    = DATAPATH+"Neill09_table2.dat"
_fileConley11       = DATAPATH+"Conley_2011_data.txt"
_fileKelly09        = DATAPATH+"Kelly09_data.dat"
# -- own stuff
_fileMorphoCompil   = DATAPATH+"morphology_R14_new_SNe.dat"

# ========================== #
# =    Get Data            = #
# ========================== #
# ----------------- #
# - Asiago        - #
# ----------------- #
def gethost_asiago_data():
    """
    This is the Asiago Catalogue containing the marphology of the Host.
    Most of the time, this is classification is ok, but still be careful
    """
    return load_pkl(_fileAsiogoDico)

# ----------------- #
# - Riess Data    - #
# ----------------- #    
def gethost_riess_hicken09_morphology_data():
    """
    Data from Adam Riess having the morphology classification
    of the Hicken 09 host data.
    """
    data = [l for l in open(_fileRiessH09morpho).read().splitlines() if l[0] != "#"]
    dico = {}
    for l in data:
        A = l.split()
        if A[0].startswith("0") or A[0].startswith("1"):
            name = "SN20%s"%A[0]
        else:
            name = "SN19%s"%A[0]
        if A[1] == "N/A":
            IDm = None
        else:
            IDm = np.int(A[1])
        dico_ = {
            "object":name,
            "IDmorpho":IDm,
            "hostname":A[2],
            "morphology":_morphoHubbleID2morpho_(IDm)
            }
        dico[dico_['object']] = dico_
        
    return dico

# -------------------------- #
# - Betoule 2014 JLA Data  - #
# -------------------------- #
def get_betoule14_data():
    return gethost_betoule14_data()

def gethost_betoule14_data():
    """
    = The JLA data = 
    """
    set_names = ["SNLS","SDSS","LowZ","HST"]
    datas = [l for l in open(_fileBetoule14).read().splitlines() if l[0] != "#"]
    dicoB = {}
    for data in datas:
        name,zcmb,zhel,dz,mb,dmb,x1,dx1,color,dcolor,_3rdvar,d3rdvar,cov_m_s,cov_m_c,cov_s_c,set_ = data.split()
        dico_ = {}
        dico_["object"]     = format_snname(name)
        dico_['JLAname']    = name
        dico_['zcmb']       = np.float(zcmb)
        dico_['zhelio']     = np.float(zhel)
        dico_['zhelio_err'] = np.float(dz)
        dico_['zcmb_err']   = np.float(dz)

        dico_['mu_b']     = np.float(mb)
        dico_['mu_b_err'] = np.float(dmb)

        dico_['x1']        = np.float(x1)
        dico_['x1_err']    = np.float(dx1)
        dico_['color']     = np.float(color)
        dico_['color_err'] = np.float(dcolor)

        dico_['mass']      = np.float(_3rdvar)
        dico_['mass_err']  = np.float(d3rdvar)

        dico_['cov_m_s']   = np.float(cov_m_s)
        dico_['cov_m_c']   = np.float(cov_m_c)
        dico_['cov_s_c']   = np.float(cov_s_c)
        dico_['setID']     = np.int(set_)
        dico_['sample']    = set_names[dico_['setID']-1]
        
        dicoB[dico_["object"]] = dico_
        
    return dicoB

# ------------------------ #
# - Silverman 2012 Data  - #
# ------------------------ #
def gethost_silverman12_data():
    """
    """
    return load_pkl(_fileSilverman12)


# ----------------------- #
# Childress Data        - #
# ----------------------- #
def gethost_childress13_data(formated_dict=True):
    """
    Data from Childress et al. 2013.
    
    ID: 1 = SNLS
        2 = SDSS
        3 = SNfactory
        4 = Low-Z
    """
    sampleID = ["SNLS","SDSS","SNf","LowZ"]
    dico = load_pkl(_fileChildress13)
    if formated_dict is False:
        return dico
    # - Lets Format it -- #
    dicoF = {}
    for name,d in dico.items():
        dicoF_ = {}
        dicoF_["mass"]      = d["Mass"]
        dicoF_["mass_err"]  = d["eMass"]
        dicoF_["object"]    = d["name"]
        dicoF_["x1"]        = d["x1"]
        dicoF_["x1_err"]    = d["ex1"]
        dicoF_["color"]     = d["c"]
        dicoF_["color_err"] = d["ec"]
        dicoF_["HR"]        = d["HR"]
        dicoF_["HR_err"]    = d["eHR"]
        dicoF_["oHR"]       = d["oHR"]
        dicoF_["oHR_err"]   = d["eoHR"]
        dicoF_["zcmb"]      = d["z"]
        dicoF_["zcmb_err"]  = d["ez"]
        dicoF_["setID"]     = int(d["ID"])
        dicoF_["sample"]    = sampleID[dicoF_["setID"]-1]
        dicoF_["object"]    = name
        dicoF[name] = dicoF_
        
    return dicoF

def gethost_childress_snf_data():
    """
    This is the SNfactory global host detailed dico
    """
    def _key_to_val_(key_,val_,dico_):
        """
        """
        if len(val_)==3:
            dico_[key_]        = val_[1]
            dico_[key_+"_err"] = np.nanmean(np.abs([val_[1]-val_[0],val_[1]-val_[2]]))
        else:
            dico_[key_]        = np.NaN
            dico_[key_+"_err"] = np.NaN
                            
    dico = load_pkl(_fileChildressSNf)
    for d in dico.values():
        ssfr  = np.asarray(d['ssfr']).copy()
        mass  = np.asarray(d['mass']).copy()
        metal = np.asarray(d['metal']).copy()
        _key_to_val_("ssfr",ssfr,d)
        _key_to_val_("mass",mass,d)
        _key_to_val_("metal",metal,d)
            
    return dico

# ----------------- #
# - Neill 2009    - #
# ----------------- #
def gethost_neill09_data():
    """
    
    from Neill et al 2009 (http://cdsads.u-strasbg.fr/abs/2009ApJ...707.1449N)
    Info: stretch errors for SNe with s < 0.7 have been multiplied by 3.0
    *Status* : L - low stretch: s < 0.80,
               R - red: C > 0.7,
               H - in Hubble flow: z > 0.0133,
               C - eligible for cosmology fitting
               
    """
    dico1 = _Neilldico_from_table1_()
    dico2 = _Neilldico_from_table2_()
    Keys = ["Host_Type","Host_mass","Host_ssfr","Host_age"]
    for d in dico1.values():
        Host_name = d["hostname"]
        if Host_name not in dico2.keys():
            for Key in Keys:
                d[Key] = np.NaN
                if Key !="Host_Type":
                    d[Key+"_err"] = np.NaN
        else:
            dH = dico2[Host_name]
            if dH["Type"] == "PEC":
                dH["Type"] = "11"

            d["morphoVaucouleurID"]    = np.round(np.float(dH["Type"]))
            d["morphology"]    = _morphoDeVaucouleurID2morpho_(d["morphoVaucouleurID"])
            d["mass"]          = dH["Mass_"]
            d["mass_err"]      = np.nanmean([dH["Mass_"]-dH["Mass_m"],dH["Mass_p"]-dH["Mass_"]])
            d["ssfr"]          = dH["sSFR_"]
            d["ssfr_err"]      = np.nanmean([dH["sSFR_"]-dH["sSFR_m"],dH["sSFR_p"]-dH["sSFR_"]])
            d["hostage"]       = dH["Age_"]
            d["hostage_err"]   = np.nanmean([dH["Age_"]-dH["Age_m"],dH["Age_p"]-dH["Age_"]])
            d["ssfr"]       = dH["sSFR_"]
            d["ssfr_err"]   = np.nanmean([dH["sSFR_"]-dH["sSFR_m"],dH["sSFR_p"]-dH["sSFR_"]])
            
        d["zcmb"] = 10**d["cz"] / constants.c.value
        
    return dico1

# ----------------- #
# - Kellly 2009   - #
# ----------------- #
def gethost_kelly09_data():
    """
    """
    dico = {}
    data = [l for l in open(_fileKelly09) if l[0] !="#"]
    Keys = ["object","zcmb","ibandRadius","Mass_m","Mass","Mass_p"]
    for d in data:
        dico_ = {}
        values = d.split()
        for k,Key in enumerate(Keys):
            if k==0:
                dico_[Key] = values[k]
            else:
                dico_[Key] = np.float(values[k])
        dico[values[0]] = dico_

    dicoK = {}
    for d in dico.values():
        dK_ = {}
        dK_["mass"] = d["Mass"]
        dK_["mass_err"] = np.nanmean([d["Mass"]-d["Mass_m"],d["Mass_p"]-d["Mass"]])
        for K in Keys[:3]:
            dK_[K] =d[K]
        dicoK[d['object']] = dK_

    return dicoK

# ----------------- #
# - Conley 2011   - #
# ----------------- #
def gethost_conley11_data():
    """
    """
    datas = [l.split() for l in open(_fileConley11) if l[0]!="#"]
    dico = {}
    for d in datas:
        dico_ = {}
        name = d[0]
        dico_["zcmb"] , dico_["zcmb_err"],dico_['mB'],dico_['mB_err'],\
          dico_['stretch'],dico_['stretch_err'], \
          dico_['color'],dico_['color_err'],dico_['mass'] = np.asarray(d[1:-2],dtype="float")
        dico_['filters'],obj_ = d[-2:]
          
        nameDB = format_snname(name)
        dico_["mass_err"] = np.NaN
        dico_['object']   = nameDB
        dico[nameDB] = dico_

    return dico

# ----------------- #
# - Own Stuff     - #
# ----------------- #
def gethost_our_morphology_data():
    """
    This is a compilation for some SNe~Ia (nearby once)
      It is mainly based on Asiago, but Greg Aldering made some
    """
    dico = {}
    datas = [l for l in open(_fileMorphoCompil).read().splitlines() if l[0] != "#"]
    for l in datas:
        A = l.split()
        dico_ = {
            "object":A[0],
            "zcmb":np.float(A[1]),
            "inPS1": np.bool(np.int(A[2])),
            "inUn3": np.bool(np.int(A[3])),
            "inK14": np.bool(np.int(A[4])),
            "inR14": np.bool(np.int(A[5])),
            "Fchart":np.bool(np.int(A[6])),
            "Ra":np.float(A[7].split(",")[0]),"Dec":np.float(A[8].split(",")[0]),
            "morphology":A[9]
            }
        dico[dico_["object"]] = dico_
    return dico

# ======================= #
# =  Tools Data         = #
# ======================= #
def _morphoHubbleID2morpho_(ID):
    """
    """
    if ID is None:
        return None

    IDmorpho = np.asarray([1,2,3,4,5,6,7,8,9,10,11])
    Morpho   = np.asarray(["E","E/S0","S0","S0a","Sa","Sab","Sb","Sbc","Sc","Scd","Sd/Irr"])
    i = np.argwhere(IDmorpho == ID)
    if len(i)==0:
        return None
    
    return Morpho[i[0]][0]

def _morphoDeVaucouleurID2morpho_(ID):
    """
    """
    if ID is None:
        return None

    IDmorpho = np.linspace(-6,11,18)
    Morpho   = np.asarray(["E","E","E/S0","E/S0","S0","S0a","S0a", # 0
                           "Sa","Sab","Sb","Sbc","Sc","Sc","Sc","Scd","Sd","Sd/Irr","Irr"])
    i = np.argwhere(IDmorpho == ID)
    if len(i)==0:
        return None
    
    return Morpho[i[0]][0]

# ---------------------------
# - Internal for creation 
def _dump_Asiago(dicoAsiago):
    """
    """
    dump_pkl(dicoAsiago,_fileAsiogoDico)
    
def _create_AsiagoDico(dump_it=False):
    """
    """
    fitsfile = pyfits.open(_fileAsiogoSrc)
    Columns  = fitsfile[1].data

    dicoAsiago = {}
    for i in range(len( Columns['SN'])):
        dico_ = {}
        dico_["object"]     = "SN"+Columns['SN'][i]
        dico_["sntype"]     = Columns['Type'][i]
        dico_["hostname"]   = Columns['Galaxy'][i]
        dico_["morphology"] = Columns['MType'][i]
        dicoAsiago[dico_["object"]] = dico_

    if dump_it:
        _dump_Asiago(dicoAsiago)
    else:
        return dicoAsiago
    
def _dump_Silverman2012_hostmorpho_snclass(dicoSilverman):
    """
    """
    dump_pkl(dicoSilverman,_fileSilverman12)
    
def _create_Silverman2012_hostmorpho_snclass(dump_it=False):
    """
    """
    head,data = pyfits.open(_fileSilverman12Src)
    
    Ncol = data.header['TFIELDS']
    Param = [data.header['TTYPE%d'%(i+1)] for i in range(Ncol)]
    dico = {}
    for d in data.data:
        d_ = {}
        for key,v in zip(Param,d):
            try:
                d_[key] = np.float(v)
            except:
                d_[key] = v.replace(" ","")
        d_['object'] = d_["SimbadName"]
        d_['Ra']     = d_["_RA"]
        d_['Dec']    = d_["_DE"]
        d_['morphology'] = d_['MType']
        d_['hostname'] = d_['Gal']
        d_['sntype'] = d_['SNtype']
        [d_.pop(k) for k in
         ["_RA","_DE","MType","Gal","SNtype"]]
        
        dico[d_['object']] = d_

    if dump_it:
        dump_Silverman2012_hostmorpho_snclass(dico)
    else:
        return dico

def _Neilldico_from_table1_():
    """
    """
    dico = {}
    data = [l for l in open(_fileNeill09Scr1) if l[0] !="#"]
    Keys = ["object","hostname","cz","stretch","m_b_max","Statut"]
    for d in data:
        dico_ = {}
        SNname_Host_cz_s_mb_statut = d.split()
        SNname = "SN%s"%SNname_Host_cz_s_mb_statut[0]
        for i,Key in enumerate(Keys[1:]):
            if Key == "cz":
                dico_[Key] = np.float(SNname_Host_cz_s_mb_statut[i+1])
            elif Key in ["stretch","m_b_max"]:
                val,valerr = SNname_Host_cz_s_mb_statut[i+1].split(';')
                dico_[Key] = np.float(val)
                dico_[Key+"_err"] = np.float(valerr)
            else:
                dico_[Key] = SNname_Host_cz_s_mb_statut[i+1]
        dico_['object'] = SNname
        dico[SNname] = dico_
    return dico

def _Neilldico_from_table2_():
    """
    """
    dico = {}
    data = [l for l in open(_fileNeill09Scr2) if l[0] !="#"]
    Keys = ["hostname","Type","Age_m","Age_","Age_p","Mass_m","Mass_","Mass_p",
            "sSFR_m","sSFR_","sSFR_p","ebmv"]
    for d in data:
        dico_ = {}
        data_split = d.split()
        for i,Key in enumerate(Keys):
            if i<2: # host name and Host tyepe (Vaucouleur)
                dico_[Key] = data_split[i]
            else:
                dico_[Key] = np.float(data_split[i])
        
        dico[data_split[0]] = dico_

    return dico


