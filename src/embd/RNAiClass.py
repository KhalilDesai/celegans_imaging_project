"""
Port of development/RNAiClass.py (Created on Jul 8, 2016, @author: Becky): the parameter-loading
and PAD/distance parts.

Changes (all output-neutral; docs/DECISIONS_AND_ISSUES.md):
* the four hardcoded normalisation dicts of setOrigin come from ``_rnai_constants`` (generated
  from the original source and tested bit-identical; R1);
* ``pNamesUse`` is built with ``dtype=object`` (ragged arrays raise in numpy >= 1.24; B12);
* methods that need embryo data or images raise NotImplementedError("image analysis not yet
  ported"): refresh_params_data, refresh_embs, getEmbryos, splitByMovement, calcAverageParams,
  asignUse, load_from_pickle, getPrev, getEmbsUsed, getRepEmb, plot*, showMaxProj*,
  makeTalignMovies, fix_head_params_RNAi, ControlClass.defineFolders/getEmbryos,
  get_control_origin, write_controlAvg_mysql;
* load_params_data falls back to the pickle/embryo path only if ``version`` is NULL, exactly as the
  original. The data layer always supplies a version (A4), so that path is never reached.
"""
import numpy as np

from . import _rnai_constants as _K
from . import db_utils_embryos, db_utils_genes
from ._compat import NOT_PORTED, not_ported
from .embdFunc import (get_param_names, get_rnai_revision_version, getNDAvgVector, getNDDistance,
                       split_params)
from .fileLookup import FOLDER_RNAI_DATA
from .params_use_assign import get_params_use
from .varLookup import FOLDER_IN, PARAM_NAMES, printLog

PARAM_NAMES_USE_GLS, PARAM_NAMES_USE_MS = get_params_use()


class RNAiClass(object):
    """
    Class that takes all embryos from a given RNAi in both GS and MS (rules and attributes that apply to RNAi condition)
    """

    def __init__(self, RNAi, verb=True, check_version=False):
        '''
        Constructor
        '''
        self.RNAi = RNAi  # EMBD number... i.e. 132
        self.verb = verb
        conn, curs = db_utils_embryos.initialize_db()
        self.id = db_utils_genes.get_id_from_embd(self.RNAi, curs)
        self.check_version = check_version
        if check_version:
            self.__version__ = get_rnai_revision_version()
        else:
            self.__version__ = None
        self.GLSfolder = FOLDER_IN + 'cropped/EMBD{0:04}/GLS/'.format(self.RNAi)
        self.MSfolder = FOLDER_IN + 'cropped/EMBD{0:04}/MS/'.format(self.RNAi)
        self.fileName = FOLDER_RNAI_DATA + 'EMBD{0:04}'.format(self.RNAi)
        self.paramsGLS = {}  # average, good parameters GLS
        self.paramsGLSstd = {}  # std, good parameters GLS
        self.paramsMS = {}  # average, good parameters MS
        self.paramsMSstd = {}  # std, good parameters MS
        self.paramsEmbsMS = {}
        self.paramsEmbsGLS = {}
        for pN in PARAM_NAMES:  # adds parameter and value, for each pN, to params dictionary
            self.paramsGLS[pN] = np.nan
            self.paramsGLSstd[pN] = np.nan
            self.paramsMS[pN] = np.nan
            self.paramsMSstd[pN] = np.nan
            self.paramsEmbsMS[pN] = np.nan
            self.paramsEmbsGLS[pN] = np.nan
        self.paramsUseGLS = []
        self.paramsUseMS = []
        self.GLSUseEmbs = []
        self.MSUseEmbs = []
        self.GLSMoveEmbs, self.GLSNoMoveEmbs = [], []
        self.MSMoveEmbs, self.MSNoMoveEmbs = [], []
        self.prev_GLS, self.prev_MS = 0, 0
        self.posND = np.ones(2 * len(PARAM_NAMES)) * np.nan
        self.stdV = np.zeros(2 * len(PARAM_NAMES))
        self.pNamesUse = np.array([PARAM_NAMES_USE_GLS, PARAM_NAMES_USE_MS], dtype=object)
        self.origin = None
        self.paramNorms = None
        self.weights = None
        self.NDFlag = False
        self.label = 'EMBD{0:04}'.format(self.RNAi)
        self.refreshed = False
        self.load_params_data()
        self.setOrigin()

    def save_params_data(self):  # makes PICKLE
        conn, curs = db_utils_embryos.initialize_db()
        columns = {}
        for key in self.paramsGLS:
            columns[key + '_GLS'] = self.paramsGLS[key]
        for key in self.paramsMS:
            columns[key + '_MS'] = self.paramsMS[key]
        columns['prev_GLS'] = self.prev_GLS
        columns['prev_MS'] = self.prev_MS
        columns['version'] = self.__version__
        db_utils_embryos.update_row(self.id, columns, curs, 'genes')
        conn.commit()
        conn.close()

        # with open(self.fileName,
        #           'wb') as output:  # this means open file (filename) as the file to save to (output) and do XYZ with it and then close it
        #     pickle.dump(self.__version__, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.paramsGLS, output,
        #                 pickle.HIGHEST_PROTOCOL)  # saves ave GLS params in first arg to fileName (output) leave 3rd arg alone.
        #     pickle.dump(self.paramsGLSstd, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.paramsMS, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.paramsMSstd, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.paramsUseMS, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.paramsUseGLS, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.paramsEmbsMS, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.paramsEmbsGLS, output, pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.prev, output, pickle.HIGHEST_PROTOCOL)

    def load_params_data(self):
        # try:
        if self.verb:
            printLog('loading from MySQL...')
        conn, curs = db_utils_embryos.initialize_db()
        p_names = get_param_names() + ['prev_GLS', 'prev_MS', 'version']
        params = db_utils_embryos.get_row(self.id, curs, 'genes', p_names)
        self.paramsGLS, self.paramsMS, other = split_params(params)
        self.paramsUseGLS = [key for key in self.paramsGLS if not np.isnan(self.paramsGLS[key])]
        self.paramsUseMS = [key for key in self.paramsMS if not np.isnan(self.paramsMS[key])]
        self.prev_GLS = other['prev_GLS']
        self.prev_MS = other['prev_MS']
        ver = other['version']
        if ver is None:
            self.load_from_pickle()
            self.save_params_data()
        if self.check_version and self.__version__ != ver:
            printLog('old version, updating RNA')
            self.refresh_params_data()
        else:
            self.__version__ = ver
            # except Exception as e:
            #     printLog('Failed to load from mysql: {}'.format(e))
            #     self.refresh_params_data()

    def load_from_pickle(self):
        raise NotImplementedError(NOT_PORTED + ' (RNAiClass.load_from_pickle: genes row has no version; '
                                  'the original would reprocess embryos from images)')

    refresh_params_data = not_ported('RNAiClass.refresh_params_data')
    refresh_embs = not_ported('RNAiClass.refresh_embs')
    asignUse = not_ported('RNAiClass.asignUse')
    getEmbryos = not_ported('RNAiClass.getEmbryos')
    splitByMovement = not_ported('RNAiClass.splitByMovement')
    calcAverageParams = not_ported('RNAiClass.calcAverageParams')

    def setOrigin(self, origin=None, norms=None, weights=None):
        # PORT: the four ~300-entry dict literals of the original (development/RNAiClass.py, setOrigin,
        # lines 256-778) live in _rnai_constants.py, generated by tools/extract_rnai_constants.py and tested
        # bit-identical. The original also contains a commented-out, unused block of per-parameter KS-test
        # weights (`KStests = {...}` and its normalisation, lines 780-917); it is not reproduced here. Origin = control mean
        # (controlOrigin); norms = std across RNAi conditions (rnaiStd); controlStd only feeds self.stdV.
        self.controlOrigin = dict(_K.controlOrigin)
        self.controlStd = dict(_K.controlStd)
        self.rnaiMean = dict(_K.rnaiMean)
        self.rnaiStd = dict(_K.rnaiStd)

        if weights is None:
            self.weights = {}
            for key in self.controlOrigin:
                self.weights[key] = 1.
        else:
            self.weights = weights

        if origin is None:
            self.origin = self.controlOrigin
        else:
            self.origin = origin
        if norms is None:
            self.paramNorms = self.rnaiStd
        else:
            self.paramNorms = norms

    def setNDPosition(self):  # paramsWT and paramsNorm (i.e. std for WT or RNAi) are dictionaries
        if self.origin is None:
            self.setOrigin()

        for i in range(len(PARAM_NAMES)):  # iterates through the index of parameter names to analyze
            pN = PARAM_NAMES[i]  # pN is the string value associated with index i
            if pN in self.paramsUseGLS and pN in PARAM_NAMES_USE_GLS:  # checks to see if parameter name is in paramsUseGLS(good usable parameter RNAiclass dictionary)
                if self.paramNorms[pN + '_GLS'] == 0:
                    self.posND[i] = self.weights[pN + '_GLS'] * (self.paramsGLS[pN] - self.origin[pN + '_GLS'])
                else:
                    self.posND[i] = self.weights[pN + '_GLS'] * (self.paramsGLS[pN] - self.origin[pN + '_GLS']) / \
                                    self.paramNorms[
                                        pN + '_GLS']  # calculates normalized value for each parameter and stores it in array posND
                    self.stdV[i] = self.controlStd[pN + '_GLS'] / self.paramNorms[pN + '_GLS']
            if pN in self.paramsUseMS and pN in PARAM_NAMES_USE_MS:  # checks to see if parameter name is in paramsUseGLS(good usable parameter RNAiclass dictionary)
                if self.paramNorms[pN + '_MS'] == 0:
                    self.posND[i + len(PARAM_NAMES)] = self.weights[pN + '_MS'] * (
                        self.paramsMS[pN] - self.origin[pN + '_MS'])
                else:
                    self.posND[i + len(PARAM_NAMES)] = self.weights[pN + '_MS'] * (
                        self.paramsMS[pN] - self.origin[pN + '_MS']) / self.paramNorms[
                                                           pN + '_MS']  # calculates normalized value for each parameter and stores it in array posND
                    self.stdV[i + len(PARAM_NAMES)] = self.controlStd[pN + '_MS'] / self.paramNorms[pN + '_MS']
        self.NDFlag = True

    #         self.posND[self.stdV>0.5]=0

    def getDims(self):  # returns the total number of considered dimensions
        return len([pN for pN in self.paramsUseGLS if pN in PARAM_NAMES_USE_GLS] +
                   [pN for pN in self.paramsUseMS if pN in PARAM_NAMES_USE_MS])

    def getOverlapDims(self, RNAi):  #
        n = 0  # overlap
        pNGLS = [pN for pN in self.paramsUseGLS if pN in PARAM_NAMES_USE_GLS]
        pNMS = [pN for pN in self.paramsUseMS if pN in PARAM_NAMES_USE_MS]
        pNGLSRNAi = [pN for pN in RNAi.paramsUseGLS if pN in PARAM_NAMES_USE_GLS]
        pNMSRNAi = [pN for pN in RNAi.paramsUseMS if pN in PARAM_NAMES_USE_MS]
        for pN in pNGLS:  # counts the number of overlapping parameters between RNAi_1 and RNAi_2 in GLS
            if pN in pNGLSRNAi:
                n += 1  # adds 1 to the value of n
        for pN in pNMS:  # counts the number of overlapping parameters between RNAi_1 and RNAi_2 in MS
            if pN in pNMSRNAi:
                n += 1  # adds 1 to the value of n
        return n

    def getDist2Zero(self):
        if not self.NDFlag:
            self.setNDPosition()
        return getNDDistance(self.posND, 'control', self.pNamesUse)

    def getDistance(self, RNAi):
        if not self.NDFlag:
            self.setNDPosition()
        if not RNAi.NDFlag:
            RNAi.setNDPosition()
        return getNDDistance(self.posND, RNAi.posND, self.pNamesUse)

    def getMiddle(self, RNAi):
        if self.NDFlag:
            return getNDAvgVector([self.posND, RNAi.posND])
        else:
            printLog('WARNING: origin is not set yet')
            return np.nan

    def getPAD(self, RNAi):
        """
        calculates PAD between two RNAi conditions
        :param RNAi:
        :return:
        """
        dist = self.getDistance(RNAi)
        if np.isinf(dist):
            return 0.
        pad = 1. - dist / (self.getDist2Zero() + RNAi.getDist2Zero())
        if pad > 0:
            return pad
        else:
            return 0

            #     def getPAD(self, RNAi): #calculates PAD between two RNAi conditions
            #         if not self.NDFlag: self.setNDPosition()
            #         if not RNAi.NDFlag: RNAi.setNDPosition()
            #         avgV = embdFunc.getNDAvgVector((np.abs(self.posND),np.abs(RNAi.posND)))
            #         dist = embdFunc.getNDDistance(self.posND/avgV/2, RNAi.posND/avgV/2, self.pNamesUse)
            #         return 1.-dist/np.sqrt(len(PARAM_NAMES_USE_GLS)+len(PARAM_NAMES_USE_MS))

    def get_cosine_dist(self, RNAi):
        """
        takes in two RNAi objects and calculates the cosine distance between the two conditions (vectors in
        n-dimensional space)- cosine of 0 is 1 (maximally similar), a cosine of 90deg is 0.
        :param RNAi: RNAi object
        :return: cosine distance
        """
        v1 = self.posND
        v2 = RNAi.posND

        N = len(self.pNamesUse[0]) + len(self.pNamesUse[1])  # total number of parameters considered GLS + MS
        vOverlap = np.invert(np.isnan(
            v1 * v2))  # overlapping vector, checks whether a parameter (coordinate) is not nan for each considered group (RNAi)
        if sum(vOverlap) == 0:
            return np.inf  # if there are no overlapping coordinates, returns infinity
        dotp_overlap = np.dot(v1[vOverlap],v2[vOverlap])  # finds the dot product for overlapping coordinates within the
        # vector (this is equiv. to x1*x2 + y1*y2 + z1*z2...)
        i1 = np.sum(np.invert(np.isnan(v1))) - np.sum(
            vOverlap)  # sums all coordinates that are not equal to nan for v1 and subtracts overlap
        i2 = np.sum(np.invert(np.isnan(v2))) - np.sum(
            vOverlap)  # sums all coordinates that are not equal to nan for v2 and subtracts overlap
        weight = 1/ (1. * N / (N - i1 - i2))  # uses same weighting as used for distance calculation, but inverts it
        # since higher values have opposite meaning in euc distance vs. cosine dis
        print("i1 = {0}, i2 = {1}, overlap = {2}, weight = {3}".format(i1, i2, np.sum(vOverlap), weight))
        # cosine_dist = dotp_overlap / (self.getDist2Zero() * RNAi.getDist2Zero())
        cosine_dist = dotp_overlap / (np.linalg.norm(v1[vOverlap]) * np.linalg.norm(v2[vOverlap]))

        adjusted_cos_dist = cosine_dist * weight
        if adjusted_cos_dist > 0:
            return adjusted_cos_dist
        else:
            return 0

    def get_norm_dist(self, r2):
        """
        calculates PAD and normalized distance between two RNAi conditions
        :param r2: RNAiClass object
        :return: PAD and distance to r2 corrected by number of dimensions
        """
        return np.exp(-self.getDistance(r2) / len(PARAM_NAMES_USE_GLS + PARAM_NAMES_USE_MS))

    def getSimParams(self, RNAi):
        if not self.NDFlag: self.setNDPosition()
        if not RNAi.NDFlag: RNAi.setNDPosition()
        diff = np.abs(self.posND - RNAi.posND)
        inds = np.argsort(diff)
        names = np.array([pN + '_GLS' for pN in PARAM_NAMES] + [pN + '_MS' for pN in PARAM_NAMES])
        names = names[inds]
        posOrdered = self.posND[inds]
        bestPars, bestParVals = [], []
        i = 0
        while len(bestPars) < 10 and i < names.size:
            pN = names[i]
            if pN[-3:] == 'GLS' and pN[:-4] in PARAM_NAMES_USE_GLS:
                bestPars.append(pN)
                bestParVals.append(posOrdered[i])
            elif pN[-3:] == '_MS' and pN[:-3] in PARAM_NAMES_USE_MS:
                bestPars.append(pN)
                bestParVals.append(posOrdered[i])
            i += 1
        return bestPars, bestParVals

    def getDiffParams(self, RNAi):
        if not self.NDFlag: self.setNDPosition()
        if not RNAi.NDFlag: RNAi.setNDPosition()
        diff = np.abs(self.posND - RNAi.posND)
        inds = np.argsort(diff)[::-1]
        names = np.array([pN + '_GLS' for pN in PARAM_NAMES] + [pN + '_MS' for pN in PARAM_NAMES])
        names = names[inds]
        posOrdered = self.posND[inds]
        worstPars, worstParVals = [], []
        i = 0
        while len(worstPars) < 10 and i < names.size:
            pN = names[i]
            if pN[-3:] == 'GLS' and pN[:-4] in PARAM_NAMES_USE_GLS:
                worstPars.append(pN)
                worstParVals.append(posOrdered[i])
            elif pN[-3:] == '_MS' and pN[:-3] in PARAM_NAMES_USE_MS:
                worstPars.append(pN)
                worstParVals.append(posOrdered[i])
            i += 1
        return worstPars, worstParVals

    getPrev = not_ported('RNAiClass.getPrev')
    getEmbsUsed = not_ported('RNAiClass.getEmbsUsed')
    getRepEmb = not_ported('RNAiClass.getRepEmb')
    plotSigmFitAll = not_ported('RNAiClass.plotSigmFitAll')
    plotSigmFitAvg = not_ported('RNAiClass.plotSigmFitAvg')
    plotCentMassAll = not_ported('RNAiClass.plotCentMassAll')
    plotCentMassAvg = not_ported('RNAiClass.plotCentMassAvg')
    plotMOIAll = not_ported('RNAiClass.plotMOIAll')
    plotMOIAvg = not_ported('RNAiClass.plotMOIAvg')
    showMaxProjMS = not_ported('RNAiClass.showMaxProjMS')
    showMaxProjGLS = not_ported('RNAiClass.showMaxProjGLS')
    makeTalignMovies = not_ported('RNAiClass.makeTalignMovies')
    fix_head_params_RNAi = not_ported('RNAiClass.fix_head_params_RNAi')

    def get_human_gene_from_ortholist_embd(self, ce_id):
        raise NotImplementedError('needs ortholist_results_embd.csv, which is not in the public data')


class ControlClass(RNAiClass):
    def __init__(self, verb=True):
        RNAiClass.__init__(self, 0, verb)

    asignUse = not_ported('ControlClass.asignUse')
    defineFolders = not_ported('ControlClass.defineFolders')
    getEmbryos = not_ported('ControlClass.getEmbryos')


get_control_origin = not_ported('RNAiClass.get_control_origin')
write_controlAvg_mysql = not_ported('RNAiClass.write_controlAvg_mysql')
