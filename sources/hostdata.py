#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" This module gather the host data sources """

import numpy as np
from ..io import get_local_hostdata,get_global_hostdata
from .basesources import BaseSources,SourceCollection

__all__ = ["get_localhost_source","get_globalhost_source"]

def get_localhost_source(**kwargs):
    """ """
    return LocalHostSource(**kwargs)

def get_globalhost_source(sourcelist=["childress13","betoule14","neill09",
                                      "childress_snf","our_morphology","asiago",
                                      "silverman12","conley11","kelly09",
                                      "riess_hicken09_morphology"],
                                      **kwargs):
    """ """
    return GlobalHostSource(**kwargs)


def get_global_source(sourcename):
    """This method enable to get individual global data souce"""
    source_ = BaseSources(empty=True)
    datasource = get_global_hostdata(sourcename)
    source_.create(data=datasource)
    return source_
    

# ======================================= #
#                                         #
#     Local SNfactory Data                #
#                                         #
# ======================================= #
class LocalHostSource( BaseSources ):
    """This instance gather the information about the local SNfactory data"""
    
    def create(self):
        """
        """
        self._properties["data"] = get_local_hostdata() 
        self._update_data_()

    # =========================== #
    # = Internal                = #
    # =========================== #
    
    # =========================== #
    # = Properties              = #
    # =========================== #
    # ------------------ #
    # - Masking        - #
    # ------------------ #
    @property
    def mask_Iae(self,cut_pIae=0.5,mask=None):
        """This mask tells which sn of the snnames is a Iae (pIae>cut_pIae)"""
        return (np.asarray(self.get_key("probaIae",mask=mask))>cut_pIae)

    @property
    def mask_Iaa(self,cut_pIae=0.5,mask=None):
        """This mask tells which sn of the snnames is a Iae (pIae>cut_pIae)"""
        # - Remark not based on mask_Iae to avoid NaN issues
        return (np.asarray(self.get_key("probaIae",mask=mask))<=cut_pIae)
    

# ======================================= #
#                                         #
#     Global  Data                        #
#                                         #
# ======================================= #
class GlobalHostSource( SourceCollection ):
    """
    """
    def create(self,sources=["childress13","betoule14","neill09",
                             "childress_snf","our_morphology","asiago",
                             "silverman12","conley11","kelly09",
                             "riess_hicken09_morphology"],
                        build=True,snnames=None):
        """
        """
        for source in sources:
            self.add_source(get_global_source(source),source)

        if build:
            self.build(snnames=snnames)

    def build(self, keys=["mass","metal","ssfr","morphology","zcmb"],
              snnames=None,sourcemask=None,exclude_asiagonames=True):
        """
        This methods enable to build the self.data dictionnary.
        the sourcemask is apply to self.sourcerank (remove mask data)
        exclude_asiagonames to avoid looping over the large asiago list...
        """
        if not self.has_sources():
            raise AttributeError("No source loaded.")
        # ---------------------
        # - Loop building Data
        usedsources = self.sourcerank if sourcemask is None else self.sourcerank[~sourcemask]
        
        used_names = snnames if snnames is not None else\
          np.unique(np.concatenate([self.sources[sourcename].snnames
                                for sourcename in usedsources
                                if not (sourcename == "asiago" and exclude_asiagonames)]))
        
        super(GlobalHostSource,self).build(keys,snnames=used_names,sourcemask=sourcemask)
        
# ======================================= #
#                                         #
#     SN Data                             #
#                                         #
# ======================================= #
