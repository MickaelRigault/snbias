#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" This module gather the sn data sources """

import warnings

import numpy as np
from ..io import get_snf_data
from .basesources import BaseSources,SourceCollection

from astrobject.utils.decorators import _autogen_docstring_inheritance

__all__ = ["get_snf_source"]

def get_snf_source(**kwargs):
    return SNfSource(**kwargs)

# ======================================= #
#                                         #
#     Local SNfactory Data                #
#                                         #
# ======================================= #
class SNfSource( SourceCollection ):
    """This instance gather the information about the local SNfactory data"""
    
    
    def create(self,sources=["idr","hubblizer","phreno"],
                        build=False,snnames=None):
        """
        """
        for source in sources:
            try:
                self.add_source(BaseSources(data=get_snf_data(source)),source)
            except:
                warnings.warn("Unknown SNF source: '%s'. Skipped"%source)
                
        self._update_()
        usednames = snnames if snnames is not None else self.all_avialable_names
        
        if build:
            self.build(keys=["x1","color","zcmb","HR","oHR","EWCaIIHK",
                             "EWSiII4000","EWSiII5972","EWSiII6355",
                             "Rsjb","vSiII_4128","vSiII_6355"],
                       snnames=self.all_avialable_names)
            
    # =========================== #
    # = Properties              = #
    # =========================== #
    @property
    def maskgood(self):
        """This return False if the SN is not 'good' as defined in the idr"""
        return self.masktraining * self.maskvalidation
    
    @property
    def masktraining(self):
        """The returns False if the SN is not in the 'training' idr subset"""
        if "idr" not in self.sourcerank:
            raise AttributeError("the ìdr` source is not loaded.")
        return np.asarray([self.get_source_sn_key("idr",snname,"subset") != "training"
                           for snname in self.snnames],dtype="bool")

    @property
    def maskvalidation(self):
        """The returns False if the SN is not in the 'training' idr subset"""
        if "idr" not in self.sourcerank:
            raise AttributeError("the ìdr` source is not loaded.")
        return np.asarray([self.get_source_sn_key("idr",snname,"subset") != "validation"
                           for snname in self.snnames],dtype="bool")
    
        
