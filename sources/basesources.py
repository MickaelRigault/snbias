#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import warnings

from astropy.table import Table
from astrobject.astrobject.baseobject import BaseObject
from astrobject.astrobject.collection import BaseCollection
from astrobject.utils import statbox
from astrobject.utils.tools import kwargs_update
from astrobject.utils.decorators import _autogen_docstring_inheritance




def merge_sourcecollections(sources,force_idmerging=True,snnames=None,verbose=True):
    """
    This enables to merge two SourceCollection into a single one.
    set snnames to only build the new sourcecollection on those supernovae
    
    """
    NewCollection = SourceCollection(empty=True)
    buildin_keys=[]
    buildin_snnames=[]
    for i,s in enumerate(sources):
        # -----------------
        # - Test Input
        if "__nature__" not in dir(s) or s.__nature__ != "Source":
            raise TypeError("the source #%d is not a `snbias Source`"%i)
        # -- This is a SourceCollection
        if "add_source" in dir(s):
            buildin_keys+=s.keys.tolist()
            buildin_snnames+=s.snnames.tolist()
            for srank in s.sourcerank:
                if verbose: print "adding %s"%srank 
                NewCollection.add_source(s.sources[srank],srank,force_it=force_idmerging)
        # -- This is a BaseSource
        else:
            warnings.warn("the source #%d is not a `snbias SourceCollection so has associated name`"%i)
            
    NewCollection._update_()
    if verbose: print "-> build the new SourceCollection"
    keys_ = np.unique([k for k in buildin_keys if "_source" not in k and "_err" not in k])
    snnames_ = snnames if snnames is not None else np.unique(buildin_snnames)
    if verbose: print "-> %d object to load and %d keys: "%(len(snnames_),len(keys_))+"\n"+", ".join(keys_.tolist())
    NewCollection.build(keys=keys_,snnames=snnames_)
    return NewCollection

# ======================================= #
#                                         #
#     Mother of Sources                   #
#                                         #
# ======================================= #
class BaseSources( BaseObject ):
    """
    This modules has the fundamental methods that all sources will share
    """
    __nature__ = "Source"
    _error_index = "_err"
    _properties_keys = ["data"]
    _derived_properties_keys = ["keys","snnames"]

    # =========================== #
    # = Construction            = #
    # =========================== #
    def __init__(self,empty=False,create=True,**kwargs):
        """Load the instance. Set empty to True to get an empty instance.
        Set create to False to avoid loading the data."""
        self.__build__()
        if empty:
            return
        if "create" not in dir(self):
            warnings.warn("No `create` method in the instance. Define one to be able to load data.")
            return
        if create:
            self.create(**kwargs)

    def create(self,data=None):
        """
        """
        if data is None:
            warnings.warn("No data to create the instance")
            return
        
        self._properties["data"] = data
        self._update_data_()
        
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def create_key(self,newkey,keyname,newkey_err=None,source=None,
                   force_it=False):
        """
        ** IMPORTANT the newkey must be ordered like self.snnames **
        """
        # -----------------------
        # -- Test Input   
        if len(newkey) != len(self.snnames):
            raise ValueError("the given new key must have the size of self.snnames")
        if newkey_err is not None and len(newkey_err) != len(newkey):
            raise ValueError("newkey and newkey_err must have the same size or set newkey to None")
        elif "__iter__" not in dir(newkey_err):
            if newkey_err is None:
                newkey_err = np.NaN
            newkey_err = [newkey_err]*len(self.snnames)
        if newkey in self.keys and not force_it:
            raise ValueError("the given keyname already exist. Set force_it to True to overwrite")
        
        if "__iter__" not in dir(source):
            source = [source]*len(self.snnames)
        elif len(source) != len(newkey):
            raise ValueError("the given source name must either be a single string or have the same size as newkey")
        
        # -----------------------
        # -- Test Input   
        for i,name in enumerate(self.snnames):
            self.data[name][keyname] = newkey[i]
            self.data[name][keyname+self._error_index] = newkey_err[i]
            self.data[name][keyname+"_source"] = source[i]
            
        self._update_data_()
            
    # ----------------- #
    # - Getter        - #
    # ----------------- #
    def get_sn_index(self,snname,mask=None):
        """return the position of the given supernova in the snnames[~mask] array"""
        if "__iter__" in dir(snname):
            return [self.get_sn_index(_name) for _name in snname]

        if mask is None:
            return np.where(self.snnames == snname)
        return np.where(self.snnames[~mask] == snname)
    
    def get_sn_key(self,snname,key,keysoftexit=False):
        """return the value associated to the given key for given supernova"""
        self._test_key_(key,softexit=keysoftexit)
        # - Several SNe
        if "__iter__" in dir(snname):
            return np.asarray([self.data[_name][key] if self._test_sn_(_name,softexit=True) and\
                     key in self.data[_name].keys() else np.NaN
                     for _name in snname ])
        # - One SN
        else: # this is a single SN
            if not self._test_sn_(snname,softexit=True) or key not in self.data[snname].keys():
                return np.NaN
            return self.data[snname][key]
        
    def get_key(self,key,mask=None):
        """this methods returns the value associated to each sn in snnames[~mask]"""
        if mask is None:
            return self.get_sn_key(self.snnames,key)
        return self.get_sn_key(self.snnames[~mask],key)
    
    # ----------------- #
    # - Plotter       - #
    # ----------------- #
    def show_key_vs_key(self,xkey,ykey,ckey=None,
                        mask=None,savefile=None,
                        ax=None,show=True,logx=False,logy=False,
                        printstat=False,
                        **kwargs):
        """
        """
        from astrobject.utils.mpladdon import figout,errorscatter
        import matplotlib.pyplot as mpl
        self._plot = {}
        
        # --------------
        # - Da Settings
        if ax is None:
            fig = mpl.figure(figsize=[8,8])
            ax  = fig.add_axes([0.1,0.1,0.8,0.8])
            ax.set_xlabel(r"$\mathrm{%s}$"%xkey.replace("_","\ "),fontsize="x-large")
            ax.set_ylabel(r"$\mathrm{%s}$"%ykey.replace("_","\ "),fontsize="x-large")
        elif "hist" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")
        else:
            fig = ax.figure
            
        ecolor = kwargs.pop("ecolor","0.7")
        prop = kwargs_update({"marker":"H","lw":1.5,
                              "zorder":5,"s":100,
                              "c":mpl.cm.Blues(0.6,0.6),
                              "edgecolors":mpl.cm.Blues(1.,1.)},
                              **kwargs)
        # ------------
        # - Da Data
        x = self.get_key(xkey,mask=mask)
        dx = self.get_key(xkey+self._error_index,mask=mask) \
          if xkey+self._error_index in self.keys else None
        y = self.get_key(ykey,mask=mask)
        dy = self.get_key(ykey+self._error_index,mask=mask) \
          if ykey+self._error_index in self.keys else None
        flagok = (x==x) * (y==y)
        if ckey is not None:
            c =  self.get_key(ckey,mask=mask)
            flagok = flagok * (c==c)
            prop["c"] = c[flagok]

        if printstat:
            intro =" Statistic for %s vs. %s"%(xkey,ykey)
            if ckey is not None:
                intro+=" having %s-values"%(ckey)
            print "".center(80,"*")
            print ("  %s  "%intro).center(80,"*")
            print "".center(80,"*")
            print "%d objects"%len(x[flagok])
            foldpc = 10
            st = statbox.get_kfolded_significance("spearman_rank_coef",[x[flagok],y[flagok]],foldpc=10,nsample=1000)
            st50,st5,st95 = np.percentile(st,[50,5,95])
            print "spearman-rank significance level:  %.2f (+%.2f -%.2f) (50%% -5,+95 based on kfolding-%d%%) )"%(st50,st5,st95,foldpc)
            
        # ------------
        # - Da Plot
        pl = ax.scatter(x[flagok],y[flagok],**prop)
        if dx is not None or dy is not None:
            err = ax.errorscatter(x[flagok],y[flagok],dx[flagok],dy[flagok],
                                  ecolor=ecolor)
        else:
            err = None
        # -----------------
        # - log or not log
        if logx:
            ax.set_xscale("log")
        if logy:
            ax.set_yscale("log")
        self._plot['ax'] = ax
        self._plot['fig'] = fig
        self._plot['plot'] = {"scatter":pl,
                              "error":err}
        self._plot['prop'] = prop
        fig.figout(savefile=savefile,show=show)
        return self._plot
    
    def show_key_hist(self,key,mask=None,log=False,savefile=None,
                      ax=None,show=True,**kwargs):
        """
        This methods enables to display an histogram of the given key-values

        Parameters:
        -----------

        key: [string]               the key you want to display (using self.get_key)

        - options -
        
        savefile: [None/string]     give a name (without extention) where the figure will
                                    be saved

        mask: [bool-array]          mask the snnames sources (self.snnames[~mask] used)
        
        log: [bool]                 display the log of the key-values

        ax: [None/ matplotlib.Axes] give the axis where the plot will be made

        show: [bool]                if not saved (savefile=None) the plot will
                                    be shown except if this is False

        -- **kwargs goes to matplotlib histogram's method --

        Return
        ------
        dictionnary having the ploting information (also recorded in self._plot
        
        """
        from astrobject.utils.mpladdon import figout
        import matplotlib.pyplot as mpl
        self._plot = {}
        
        # --------------
        # - Da Settings
        if ax is None:
            fig = mpl.figure(figsize=[8,8])
            ax  = fig.add_axes([0.1,0.1,0.8,0.8])
            ax.set_xlabel(r"$\mathrm{%s}$"%key,fontsize="x-large")
        elif "hist" not in dir(ax):
            raise TypeError("The given 'ax' most likely is not a matplotlib axes. "+\
                             "No imshow available")
        else:
            fig = ax.figure

        prop = kwargs_update({"histtype":"step","lw":2,"fill":True,
                              "fc":mpl.cm.Blues(0.6,0.4),"ec":mpl.cm.Blues(0.8,0.8)},
                              **kwargs)
        # ------------
        # - Da Plot
        x = np.log10(self.get_key(key,mask=mask)) if log else self.get_key(key,mask=mask)
        pl = ax.hist(x[x==x],**prop)
        self._plot['ax'] = ax
        self._plot['fig'] = fig
        self._plot['plot'] = pl
        self._plot['prop'] = prop
        fig.figout(savefile=savefile,show=show)
        return self._plot
    # =========================== #
    # = Internal                = #
    # =========================== #
    def _update_data_(self):
        """ Run this if you change the data"""
        # - These is saved for speed reasons
        self._derived_properties["snnames"] = np.sort(self.data.keys())
        self._derived_properties["keys"] = np.unique(np.concatenate([d.keys() for d in self.data.values()]))
        
    # ---------------- #
    # - Tests        - #
    # ---------------- #
    def _test_key_(self,key, softexit=False):
        """raise an exception (or return False if softexit is True) if the key is not known"""
        if key not in self.keys:
            if softexit:
                return False
            raise ValueError("'%s' is not a known key, these are: "%key+", ".join(self.keys.tolist()))
        return True
    
    def _test_sn_(self,snname, softexit=False):
        """raise an exception (or return False if softexit is True)if the sn is not known"""
        if snname not in self.snnames:
            if softexit:
                return False
            raise ValueError("'%s' is not a known supernovae"%snname)
        return True
        
    # =========================== #
    # = Properties              = #
    # =========================== #
    # ------------------- #
    # - Data In         - #
    # ------------------- #
    @property
    def data(self):
        """The data source of the instance. Dictionnary"""
        if self._properties["data"] is None:
            self._properties["data"] = {}
        return self._properties["data"]

    def has_data(self):
        return len(self.data.keys())>0

    # -------------------- #
    # - Associated Info  - #
    # -------------------- #
    @property
    def snnames(self):
        """list of all the supernova"""
        return self._derived_properties["snnames"]
    
    @property
    def keys(self):
        """list of the keys shared by all SNe Ia entries"""
        if not self.has_data():
            return None
        return self._derived_properties["keys"]
    
# ======================================= #
#                                         #
#     Mother of Sources                   #
#                                         #
# ======================================= #
class SourceCollection( BaseSources, BaseCollection ):
    """
    """
    def __build__(self):
        """
        """
        for new_prop in ["handler","sourcerank"]:
            self._properties_keys.append(new_prop)

        super(SourceCollection,self).__build__()
        
    
    def add_source(self,newsource,id_,rank=-1,
                   force_it=False,keys=None):
        """
        Add to the current SourceCollection a new source.
        This will not change the self.snnames but potentially new objects
        will be visible in self.all_avialable_names. The same goes for self.keys
        even though you can direction add a 'keys' in the current instance by giving
        it (string or array of string) in the input argument 'keys'
        """
        if self.has_data() and id_ in self.list_sources and not force_it:
            raise ValueError("the id %s is already used. change it or set force_it to True to overwrite it.")

        self._handler[id_] = newsource
        # --------------
        # source ranking
        if rank == -1:
            self.sourcerank.append(id_)
        else:
            if rank<0:
                rank+=1
            self.sourcerank.insert(rank,id_)
            
        if keys is not None:
            if "__iter__" not in dir(keys): keys=[keys]
            [self.build_key(key_) for key_ in keys
            if "_err" not in key_ and "_source" not in key_]
            
    def add_sourcecollection(self,sourcecollection,force_it=True):
        """
        """
        if "__nature__" not in dir(sourcecollection) or sourcecollection.__nature__ != "Source":
            raise TypeError("the given source is  a `snbias SourceCollection`")
        if "add_source" not in dir(sourcecollection):
            raise TypeError("the given source is a Source not a SourceCollection. use `self.add_source`")
        
        [self.add_source(sourcecollection.sources[srank],srank,force_it=force_it)
         for srank in sourcecollection.sourcerank] 
        [self.build_key(key_) for key_ in sourcecollection.keys
         if "_err" not in key_ and "_source" not in key_]
        self._update_()
        
    @_autogen_docstring_inheritance(BaseCollection.remove,"BaseCollection.remove")
    def remove(self,id_):
        #
        # Add the rank value removal
        #
        self.sourcerank.pop(self.get_source_rank(id_))
        super(SourceCollection,self).remove(id_)
        
        
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def create(self):
        raise NotImplementedError("SourceCollection basic create function is not avialable.")
    
    # ----------------- #
    # - BUILDING DATA - #     
    # ----------------- #
    def build(self, keys, snnames, sourcemask=None):
        """
        This methods enable to build the self.data dictionnary.
        the sourcemask is apply to self.sourcerank (remove mask data)
        exclude_asiagonames to avoid looping over the large asiago list...
        """
        if not self.has_sources():
            raise AttributeError("No source loaded.")
        # ---------------------
        # - Initializing
        for name in snnames:
            self.data[name] = {}
        self._update_data_()
        # ---------------------
        # - Building the keys
        [self.build_key(key_,sourcemask=sourcemask) for key_ in keys]

    def build_key(self,key,sourcemask=None):
        """This methods enables to add the given key to the data dictionnary"""
        
        if not self.has_data():
            raise AttributeError("no current data dictionnary.")

        values_,sources = self.fetch_sn_key_sources(self.snnames,key,sourcemask=sourcemask)
        values = values_.T
        
        for i,name_ in enumerate(self.snnames):
            index = np.where(values[i]==values[i])[0]
            if len(index)==0:
                self.data[name_][key] = np.NaN
                self.data[name_][key+"_source"] = None
                self.data[name_][key+"_err"] = np.NaN
            else:
                source = sources[index[0]]
                self.data[name_][key] = values[i][index[0]]
                self.data[name_][key+"_source"] = source
                self.data[name_][key+"_err"] = self.sources[source].get_sn_key(name_,
                                                key+self.sources[source]._error_index,keysoftexit=True)

        self._update_data_()
    # ----------------- #
    # - Getter        - #     
    # ----------------- #
    def fetch_sn_key_sources(self,snname,key,sourcemask=None):
        """
        """
        sources = self.get_key_sources(key,mask=sourcemask)
            
        data = np.asarray([self.sources[source].get_sn_key(snname,key, keysoftexit=True) for source in sources])
        return data,sources
    
    def get_source_sn_key(self,sourcename,snname,key):
        """
        """
        self._test_source_(sourcename)
        return self.sources[sourcename].get_sn_key(snname,key,keysoftexit=True)

    # ------------------
    # - Source IO
    def get_source(self,sourcename):
        """returns a of the requesteed source"""
        self._test_source_(sourcename)
        return self._handler[sourcename].copy()

    def get_key_sources(self,key,mask=None,snname=None):
        """This module returns the ranklist of the sources containing the given key.
        Add the snname cut to check if the given sn is there known."""
        sources = self.sourcerank if mask is None else self.sourcerank[~mask]
        
        return [source_ for source_ in sources
                if key in self.sources[source_].keys and
                (snname is None or snname in self.sources[source_].snnames)]
    
    # ----------------- #
    # - Ranking       - #
    # ----------------- #
    def get_source_rank(self,sourcename):
        """give the source id and this returns its rank"""
        self._test_source_(sourcename)
        return np.argwhere(np.asarray(self.sourcerank) == sourcename)[0][0]
            
        
    # =========================== #
    # = Properties              = #
    # =========================== #
    @property
    def sourcerank(self):
        if self._properties["sourcerank"] is None:
            self._properties["sourcerank"] = []
        return self._properties["sourcerank"]

    @property
    def sources(self):
        return self._handler

    @property
    def all_avialable_keys(self):
        return np.sort(np.unique(np.concatenate([s_.keys for s_ in self.sources.values()])))

    @property
    def all_avialable_names(self):
        return np.sort(np.unique(np.concatenate([s_.snnames for s_ in self.sources.values()])))
