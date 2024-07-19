#!/usr/bin/env python3

# S. I. Cooper, Apr. 2024
# modified from https://github.com/bellan/VVXAnalysis/blob/c40b91f1147a784f3748d5bf2dd3ae1d28540931/Combine/scripts/tableYields.py
#
# Usage example:
# python3 scripts/tablesFromDatacard.py 2016preVFP/eejj/eejj_11jul2024_bdt_LQToDEle/scaled/tmp_card_file_LQToDEle.txt 2016postVFP/eejj/eejj_11jul2024_bdt_LQToDEle/scaled/tmp_card_file_LQToDEle.txt 2017/eejj/eejj_11jul2024_bdt_LQToDEle/scaled/tmp_card_file_LQToDEle.txt 2018/eejj/eejj_11jul2024_bdt_LQToDEle/scaled/tmp_card_file_LQToDEle.txt > limits11July_allRunII_LQToDEle/tablesFromDatacard.txt

from argparse import ArgumentParser
from optparse import OptionParser
from fnmatch import fnmatch
from string import Template
from math import sqrt, isnan
from ctypes import c_double
import os
import sys
import pandas as pd

import ROOT

from HiggsAnalysis.CombinedLimit import DatacardParser

if not 'LQANA' in os.environ:
    raise RuntimeError('Need to define $LQANA environment variable to point to your git clone of the rootNtupleAnalyzerV2 repo.')
sys.path.append(os.getenv("LQANA").rstrip("/") + "/scripts/")

from combineCommon import SeparateDatacards


class MyTemplate(Template):
    idpattern = r'(?-i:[a-zA-Z][a-zA-Z]*)'  # ignore underscores and numerals in the identifier name


class EventYield:
    '''
    A class that handles correctly the statistical error on the event yield (sum in quadrature)
    '''
    def __init__(self, val=0, err=0, hasErr=True):
        self.val = val
        self.err = err
        self.hasErr = hasErr

    def __add__(self, other):
        new = EventYield(self.val, self.err)
        new += other
        return new

    def __iadd__(self, other):
        self.val += other.val
        self.err = sqrt(self.err**2 + other.err**2)
        return self

    def __repr__(self):
        if self.hasErr:
            return '(' + str(self.val) + '+-' + str(self.err) + ')'
        else:
            return '(' + str(self.val) + ')'

    def to_string(self, fmt='%.4g'):
        if self.hasErr:
            return (fmt+' pm '+fmt) %(self.val, self.err)
        else:
            return (fmt) %(self.val)

    def __iter__(self):
        return (i for i in (self.val, self.err))


def uniq(orig):
    '''
    Create a list of unique elements respecting the order in the original list
    '''
    out = []
    seen = set()
    for e in orig:
        if(e in seen): continue
        out.append(e)
        seen.update(e)
    return out


def get_mass_from_signalName(signalName):
    return signalName.split("_")[-1].split("M")[-1] # LQ_M300


def sum_yields(*card_yields, **kwargs):
    '''
    Sums EventYields stored in dictionaries, for each process and bin
    '''
    out = {proc_name: {} for y in card_yields for proc_name in y.keys()}

    for proc_name in out.keys():
        out[proc_name] = {}
        bin_names = uniq(bin_name for card_yield in card_yields for bin_name in card_yield[proc_name].keys())
        # print("INFO: card_yield proc_name = '{}'".format(proc_name))
        # print("INFO: got bin names='{}'".format(bin_names))

        for bin_name in bin_names:
            # print("INFO: bin name='{}'".format(bin_name))
            # for card_yield in card_yields:
            #     print("INFO: this card has bin names = '{}'".format(card_yield[proc_name].keys()))
            bin_tot_y = sum((card_yield[proc_name][bin_name] for card_yield in card_yields), EventYield())
            out[proc_name][bin_name] = bin_tot_y

    return out


def find_in_shapeMap(string, shapeMap):
    '''
    Tries to find a matching pattern in a shapeMap and returns the object associated to it
    (usually another shapeMap (which is really just a dictionary) or a list
    '''
    for pattern, data in shapeMap.items():
        if fnmatch(string, pattern):
            return data
    raise LookupError('No pattern matching "%s" in shapeMap. Keys are: %s' %(string, shapeMap.keys()))



def create_dataframe(data, signalName, unblind=False):
    # Remove data_obs
    series_data = data.pop('data_obs', None)

    df = pd.DataFrame(data)
    # df.fillna(value=EventYield(0,0), inplace=True)
    # On lxplus, pandas version is 1.2.2, and fillna is bugged (does not accept object values)
    df = df.applymap(lambda x: x if not(isinstance(x, float) and isnan(x)) else EventYield(0,0))

    # Transpose: rows = samples, columns = years
    df = df.transpose()

    # Total yield of MC for each year
    df2 = df.drop([signalName])
    df.loc['Total BG'] = df2.sum()

    if unblind:
        if(series_data is not None):
            logging.debug('data series: %s', series_data)
            df.loc['Data'] = series_data
        else:
            logging.warning('data series is None')

    # Compute total yield
    df['Total'] = df.sum(axis=1)
    df = df[['Total']]

    df = df.transpose()  # transpose back for printing
    df.rename({signalName: 'LQ signal'}, axis=1, inplace=True)
    # sigYield = EventYield(int(signalName.split('M')[-1]), 0, False)
    # df.insert(0, "M_LQ", str(int(signalName.split('M')[-1])))
    df.insert(0, "M_LQ", int(signalName.split('M')[-1]))
    return df


def print_yield(df, float_format='%.4g', **kwargs):
    '''
    Prints to stdout a LaTeX table of event yields using Pandas DataFrame's to_latex()
    '''
    # Create list of bins that contains unique elements but mantains the order in which they were in data
    # this could not be archieved with a list comprehension (no uniqueness)
    # nor with sorting a set (the starting order would not be respeted)
    # bin_names = []
    # for _, proc_data in data.items():
    #     for bin_name in proc_data:
    #         if bin_name not in bin_names:
    #             bin_names.append(bin_name)

    # Compute total yield
    # total = {bin_name: [0, 0] for bin_name in bin_names}
    # for proc_name, proc_data in data.items():
    #     for bin_name, proc_bin_data in proc_data.items():
    #         total[bin_name] += proc_bin_data

    # print(df)
    dfPrint = df.applymap(lambda x: x.to_string(fmt=float_format) if isinstance(x, EventYield) else str(x))
    dfPrint = dfPrint.applymap(lambda x: x.replace('pm', u"\u00B1"))
    renaming = {}
    renaming['ZJet_amcatnlo_ptBinned_IncStitch'] = 'Z+jets'
    renaming['TTTo2L2Nu'] = 'ttbar'
    renaming['QCDFakes_DATA'] = 'QCD (data)'
    renaming['DIBOSON_nlo'] = 'Diboson'
    renaming['SingleTop'] = 'ST'
    dfPrint.rename(axis='columns', mapper=renaming, inplace=True)
    out_string = dfPrint.to_markdown(floatfmt=float_format.replace('%', ''), tablefmt="fancy_grid", index=False)
    print(out_string)

    # Convert to string using the format supplied by command line args
                   # TTTo2L2Nu     │ QCDFakes_DATA   │ DIBOSON_nlo    │ SingleTop
    formatters = [lambda x: x.to_string(fmt=float_format) if isinstance(x, EventYield) else str(x)]*len(df.columns)
    zjetName = r'$Z$+jets' # r'$\PZ$+jets'
    ttbarName = r'$t\bar{t}$' # r'$\PQt\PAQt$'
    out_string = df.to_latex(formatters=formatters, index=False)\
                   .replace('pm', r'$\pm$')\
                   .replace('ZJet\_amcatnlo\_ptBinned\_IncStitch', zjetName)\
                   .replace('TTTo2L2Nu', ttbarName)\
                   .replace('QCDFakes\_DATA', r'QCD (data)')\
                   .replace('DIBOSON\_nlo', r'Diboson')\
                   .replace('SingleTop', r'ST')
    print(out_string)


def get_yield(card, unblind=False, **kwargs):
    '''
    Retrieve the histograms in the files pointed by the shapeMap in the card
    and use their integral and error to construct a dictionary that maps
    {process: {bin: [yield +- error]}}
    '''
    out = {}
    tf_handles = {}
    total = {}
    for bin_name, bin_data in card.exp.items():
        shapeMap_bin = find_in_shapeMap(bin_name, card.shapeMap)
        # print('DEBUG: shapeMap for bin {} = {}'.format(bin_name, shapeMap_bin))
        mapping = {'CHANNEL': bin_name}
        mapping.update({'MASS': get_mass_from_signalName(card.signals[0])})

        proc_names = list(bin_data.keys())
        if(unblind):
            proc_names.append('data_obs')

        for proc_name in proc_names:
            shapeMap_proc = find_in_shapeMap(proc_name, shapeMap_bin)
            # print('DEBUG: shapeMap for process {} = {}'.format(proc_name, shapeMap_proc))
            mapping.update({'PROCESS': proc_name})
            # filepath, rootpath, _ = shapeMap_proc
            filepath, rootpath = shapeMap_proc
            filepath = MyTemplate(filepath).substitute(mapping)
            rootpath = MyTemplate(rootpath).substitute(mapping)
            # print('DEBUG post sub  filepath={}  rootpath={}'.format(filepath, rootpath))

            if not filepath in tf_handles.keys():
                tf = ROOT.TFile(filepath)
                if(not tf or not tf.IsOpen()):
                    raise FileNotFoundError(filepath)
                else:
                    tf_handles[filepath] = tf
                    # print('DEBUG - opened', filepath)
            h = tf_handles[filepath].Get(rootpath)
            error = c_double(0.)
            if(h):
                integral = h.IntegralAndError(0,-1, error)
            else:
                print('WARN Could not get "{}" from {}'.format(rootpath, filepath))
                integral = 0

            out.setdefault(proc_name, {}).setdefault(bin_name, EventYield()) 
            out[proc_name][bin_name] += EventYield(integral, error.value)

            # logging.debug('\t%s %24s = %5.3f +. %5.3f' %(bin_name, proc_name, integral, error.value))

    for filepath, handle in tf_handles.items():
        # print('DEBUG Closing', filepath)
        handle.Close()

    return out



def parse_args():
    parser = ArgumentParser(description='Reads one or more datacards, finds the histograms referenced by their shapeMaps and formats the expected and observed events in a LaTeX table')
    parser.add_argument('datacards', nargs='+')
    parser.add_argument('-u', '--unblind', action='store_true', help='print the "Observation" row')
    parser.add_argument(      '--float-format', default='%.4g', help='Format string used for floats (default: %(default)s)')
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')
    args, unknown = parser.parse_known_args()

    parser_datacard = OptionParser()
    DatacardParser.addDatacardParserOptions(parser_datacard)
    args_datacard = parser_datacard.parse_args(unknown)
    return args, args_datacard


def main():
    args, args_datacard = parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)

    cardsByMass = {}
    for idx, datacard in enumerate(args.datacards):
        print('INFO: reading', datacard)
        # each card here is a combined card
        massList, tmpCardsByMass, _ = SeparateDatacards(datacard, idx, "/tmp/")
        for mass in massList:
            if int(mass) not in cardsByMass.keys():
                cardsByMass[mass] = []
            dcard = tmpCardsByMass[mass]
            with open(dcard) as f:
                cardsByMass[mass].append(DatacardParser.parseCard(f, options=args_datacard[0]))

    print('INFO: number of cards (for LQ_M{}) ='.format(list(cardsByMass)[0]), len(list(cardsByMass.values())[0]))

    # cards[0].print_structure()
    yieldsByMass = {}
    for mass, cards in cardsByMass.items():
        if not mass in yieldsByMass.keys():
            yieldsByMass[mass] = []
        for card in cards:
            yieldsByMass[mass].append(get_yield(card, **vars(args)))

    completeDataFrame = None
    for mass, evtYields in yieldsByMass.items():
        signalName = "LQ_M{}".format(mass)
        tot_yield = create_dataframe(sum_yields(*evtYields), signalName)
        if completeDataFrame is None:
            completeDataFrame = tot_yield
        else:
            # completeDataFrame = completeDataFrame.append(tot_yield)
            completeDataFrame = pd.concat([completeDataFrame, tot_yield])
        # print('INFO: mass={}, tot_yield={}'.format(mass, tot_yield))
    # print("INFO: args=", args)
    print_yield(completeDataFrame, **vars(args))

    return 0

if __name__ == '__main__':
    main()
