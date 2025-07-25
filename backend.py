import matplotlib.pyplot as plt
import numpy as np
import pyopenms
from scipy import signal
from scipy.interpolate import CubicSpline

__fragment_types__ = {"a": "N", "b": "N", "c": "N", "x": "C", "y": "C", "z": "C"}
__fragment_mass_calculation__ = {"x": pyopenms.Residue.ResidueType.XIon,
                                 "y": pyopenms.Residue.ResidueType.YIon,
                                 "z": pyopenms.Residue.ResidueType.ZIon,
                                 "a": pyopenms.Residue.ResidueType.AIon,
                                 "b": pyopenms.Residue.ResidueType.BIon,
                                 "c": pyopenms.Residue.ResidueType.CIon}


def sorted_range(seq):  #simple function to sort fragment indices during DFS
    return sorted(range(len(seq)), key=seq.__getitem__)


def merge_dicts(dict1, dict2):
    #implementation of the dict merge, needed for merging of dicts of fragments in series
    res_dict = dict1.copy()
    for k in dict2:
        if k in res_dict:
            res_dict[k] += dict2[k]
        else:
            res_dict[k] = dict2[k]
    return res_dict


def moving_average(vec, N=3):
    # calculator of moving average values for calculation of the no-fragment penalty (gamma values according to paper)
    ma_vec = np.concatenate(
        (vec[::-1], vec, vec[::-1]))  # this is needed to correctly calculate moving averages on borders
    ma_vec.astype(float)
    res = np.convolve(ma_vec, np.ones(N) / N, mode='valid')
    lvec = len(vec)
    return res[lvec:2 * lvec - 1]


def remove_isotopes(iso, cutoff=0.2):  #removal of the low-abundant isotopes from the isotope structure of ion
    mx = 0
    for i in range(len(iso)):
        if iso[i][1] > mx:
            mx = iso[i][1]
    rem = []
    for i in range(len(iso)):
        if iso[i][1] > mx * cutoff:
            rem.append([iso[i][0], iso[i][1]])
    return np.array(rem)


def merge_intervals(set_intervals_1, set_intervals_2):  #implementation of interval merge function
    all_indexes = set()
    for k in set_intervals_1:
        if k[0] not in all_indexes:
            all_indexes.add(k[0])
        if k[1] not in all_indexes:
            all_indexes.add(k[1])
    for k in set_intervals_2:
        if k[0] not in all_indexes:
            all_indexes.add(k[0])
        if k[1] not in all_indexes:
            all_indexes.add(k[1])
    all_indexes = list(all_indexes)
    all_indexes.sort()
    intervals = []
    for i in range(len(all_indexes) - 1):
        x = (all_indexes[i] + all_indexes[i + 1]) // 2
        inside = False
        for k in set_intervals_1:
            if k[0] <= x < k[1]:
                inside = True
        for k in set_intervals_2:
            if k[0] <= x < k[1]:
                inside = True
        if inside:
            intervals.append([all_indexes[i], all_indexes[i + 1]])
    for i in reversed(range(1, len(intervals))):
        if intervals[i - 1][1] == intervals[i][0]:
            intervals[i - 1][1] = intervals[i][1]
            intervals.pop(i)
    return intervals


def interval_check(x, interval_set):
    res = False
    for interval in interval_set:
        if (x >= interval[0]) and (x < interval[1]):
            res = True
    return res


class Ppt:
    #main class that describes the chain subunit, with overloaded comparison functions
    def __init__(self, name, seq):
        self.name = name    #subunit alias, like Ub for Ubiquitin subunit
        self.seq = seq      #subunit seq, like MQIF... for Ub
        self.length = len(self.seq) #length of the subunit polypeptide

    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name


class S:
    # main class for Series fragments, like c(Ub)+2Ub, with all implemented machinery for calculation of fragments
    def __init__(self, f_type, fragmenting_ppt, intact_ppts, f_name):
        self.seqs_of_intact_ppts = None
        self.seq_of_fragmenting_unit = None
        self.f_type = f_type        #Principal type of fragment: N-terminal or C-terminal, important for correct calculation of the fragment sets
        self.f_name = f_name        #Types of fragments in Biemann's notation a,b,c,x,y,z; comparing to f_type, only needed for naming fragments
        self.fragmenting_ppt = fragmenting_ppt #Ppt corresponding to Central subunit undergoing fragmentation
        self.intact_ppts = intact_ppts  #Set of Ppt objects corresponding to intact chains in series

    def __str__(self):  #method for generation of str alias of the fragment series that shows in IT captions
        if len(self.intact_ppts) == 0:
            return self.f_name + "(" + str(self.fragmenting_ppt) + ")"
        else:
            return str(self.intact_ppts) + "+" + self.f_name + "(" + str(self.fragmenting_ppt) + ")"

    def __repr__(self): #method for generation of str alias of the fragment series that shows in IT captions
        return str(self.intact_ppts) + "+" + self.f_name + "(" + str(self.fragmenting_ppt) + ")"

    def __hash__(self): #hash calculator that is needed for correct operation of dicts and sets of fragment series
        if len(self.intact_ppts) == 0:
            for_calc = (self.f_name, self.fragmenting_ppt, 0)
        else:
            intact_ppts_keys = []
            for k in self.intact_ppts:
                intact_ppts_keys.append(str(k) + "_" + str(self.intact_ppts[k]))
            intact_ppts_keys.sort()
            for_calc = (self.f_name, self.fragmenting_ppt, (tuple(intact_ppts_keys)))
        return hash(for_calc)

    def __eq__(self, other): #comparator of the fragment series identity
        return self.__hash__() == other.__hash__()

    def generate_sequences(self): #generates AA sequences for intact and central subunits for use in fragment_seq_generator function
        self.seqs_of_intact_ppts = []
        self.seq_of_fragmenting_unit = self.fragmenting_ppt.seq
        for ppt in self.intact_ppts:
            for k in range(self.intact_ppts[ppt]):
                self.seqs_of_intact_ppts.append(ppt.seq)

    def fragment_seq_generator(self, k):    #generate seq input for pyopenms calculation of fragment isotope masses
        fr_unit_length = len(self.seq_of_fragmenting_unit)
        if self.f_type == "N":
            return "".join(self.seqs_of_intact_ppts) + self.seq_of_fragmenting_unit[:k]
        elif self.f_type == "C":
            return "".join(self.seqs_of_intact_ppts) + self.seq_of_fragmenting_unit[fr_unit_length - k:fr_unit_length]

    def fragment_mass_calculator(self, k, z, max_num_of_isotopes=20): #calculator of isotope masses using pyopenms
        result = []
        seq = pyopenms.AASequence.fromString(self.fragment_seq_generator(k))
        frag = seq.getFormula(__fragment_mass_calculation__[self.f_name], z)
        isotopes = frag.getIsotopeDistribution(pyopenms.CoarseIsotopePatternGenerator(max_num_of_isotopes))
        for iso in isotopes.getContainer():
            result.append([iso.getMZ() / z, iso.getIntensity()])
        return result


class I:    #Class related to the specific isoform and all representations of it

    V: list[Ppt]

    def __init__(self, name, init_formula, ppts):
        self.name = name    #isoform name
        self.ppts = ppts    #set of Ppt instances associated with subunits from which the particular isoform is constructed
        self.init_formula = init_formula    #encoding formula from the config file, "str" format, with aliases matching the ppts record
        self.V = [""]       #vertices in the graph-based representation
        self.E = []         #edges in the graph-based representation
        self.parsing(init_formula[0:len(init_formula)], 0)
        self.S_set = {}     #set of fragment series represented by S instances associated with the particular isoform

    def generate_dict_with_all_ppts_in_I(self): # put in self.all_units dict all Ppts presented in V as keys and number of their occurrence as values
        all_units = {}
        for k in self.V:
            if k not in all_units:
                all_units[k] = 1
            else:
                all_units[k] += 1
        self.all_units = all_units

    def compute_reachability(self):   #DFS implementation: for every edge from E calculate set of vertices V that can be reached from e
        havebeen_matrix = {}
        for e in self.E:
            havebeen = [0] * len(self.V)
            stack = [e[1]]
            while len(stack) > 0:
                cur_nV = stack.pop(-1)
                havebeen[cur_nV] = 1
                for cur_e in self.E:
                    if cur_e[0] == cur_nV:
                        stack.append(cur_e[1])
            hb_dict = {}
            for k in range(len(havebeen)):
                if havebeen[k] == 1:
                    if self.V[k] not in hb_dict:
                        hb_dict[self.V[k]] = 1
                    else:
                        hb_dict[self.V[k]] += 1
            havebeen_matrix[e] = hb_dict
        self.reachability = havebeen_matrix

    def generate_set_of_complimentary_intact_ppts(self, d, v):
        #function needed for calculation of the C-terminal series based on the complimentary N-terminal series
        result = self.all_units.copy()
        result[v] -= 1
        for k in d:
            result[k] -= d[k]
        return result

    def generate_series_of_fragment(self, f_name):
        #main function for; generate series of fragments for isoform of the "f_name" type using the procedure described in paper
        self.generate_dict_with_all_ppts_in_I()
        self.compute_reachability()
        f_type = __fragment_types__[f_name]
        S_set = []
        for ind_v in range(len(self.V)):
            conj_e = []
            conj_ne = []
            for e in self.E:
                if e[0] == ind_v:
                    conj_e.append(e)
                    conj_ne.append(e[2])
            indexes_sorted = sorted_range(conj_ne)
            conj_e_s = []
            conj_ne_s = []
            length_v = self.V[ind_v].length
            for i in indexes_sorted:
                conj_e_s.append(conj_e[i])
                conj_ne_s.append(conj_ne[i])
            if f_type == "N":
                prev_i = 0
                k = 0
                new_S_F = self.V[ind_v]
                new_S_ftype = f_type
                new_S_fname = f_name
                d = {}
                while k < len(conj_ne_s):
                    new_S_indexes = (prev_i, conj_ne_s[k])
                    prev_i = conj_ne_s[k]
                    new_S_D = d.copy()
                    s = S("N", new_S_F, new_S_D, f_name)
                    S_set.append((s, (new_S_indexes)))
                    d = merge_dicts(d, self.reachability[conj_e_s[k]])
                    k += 1
                new_S_indexes = (prev_i, length_v)
                new_S_D = d.copy()
                s = S("N", new_S_F, new_S_D, f_name)
                S_set.append((s, (new_S_indexes)))
            if f_type == "C":
                prev_i = 0
                k = 0
                new_S_F = self.V[ind_v]
                new_S_ftype = f_type
                new_S_fname = f_name
                d = {}
                while k < len(conj_ne_s):
                    new_S_indexes = (length_v - conj_ne_s[k], length_v - prev_i)
                    prev_i = conj_ne_s[k]
                    new_S_D = d.copy()
                    complimentary_S_D = self.generate_set_of_complimentary_intact_ppts(new_S_D, self.V[ind_v])
                    d = merge_dicts(d, self.reachability[conj_e_s[k]])
                    s = S("C", new_S_F, complimentary_S_D, f_name)
                    S_set.append((s, (new_S_indexes)))
                    k += 1
                new_S_indexes = (0, length_v - prev_i)
                new_S_D = d.copy()
                complimentary_S_D = self.generate_set_of_complimentary_intact_ppts(new_S_D, self.V[ind_v])
                s = S("C", new_S_F, complimentary_S_D, f_name)
                S_set.append((s, (new_S_indexes)))
        final_S_set = self.S_set.copy()
        for q in S_set:
            real_value = False
            for k in final_S_set:
                if q[0] == k:
                    real_value = k
            if real_value == False:
                final_S_set[q[0]] = [q[1]]
            else:
                new_intervals = merge_intervals(final_S_set[real_value], [q[1]])
                final_S_set[real_value] = new_intervals
        self.S_set = final_S_set


    def parsing(self, cur_formula, parent): #dynamic programming recursion function for reading the encoding bracketed sequence
        """

        :type cur_formula: str
        """
        commas = []
        q = 0
        for i in range(len(cur_formula)):
            if cur_formula[i] == "(":
                q += 1
            elif cur_formula[i] == ")":
                q -= 1
            if (cur_formula[i] == ",") and (q == 0):
                commas.append(i)
        if len(commas) == 0:
            v_name = cur_formula
            self.V[parent] = self.ppts[v_name]
        else:
            v_name = cur_formula[:commas[0]]
            self.V[parent] = self.ppts[v_name]
            for k in range(0, len(commas)):
                if k == len(commas) - 1:
                    operating_formula = cur_formula[commas[k]:]
                else:
                    operating_formula = cur_formula[commas[k]:commas[k + 1]]
                j = 0
                while operating_formula[j] != "(":
                    j += 1
                linkage = operating_formula[1:j]
                child_sequence = operating_formula[j + 1:len(operating_formula) - 1]
                self.V.append("")
                self.E.append((parent, len(self.V) - 1, int(linkage)))
                self.parsing(child_sequence, len(self.V) - 1)


class InternalMassSpectrum:
    # class associated with internal representation of the mass spectrum
    def __init__(self, high_mz, bin_size):
        self.peaks = None
        self.mz = np.arange(0, int(high_mz / bin_size), dtype=float) * bin_size
        self.intensities = np.zeros(len(self.mz))
        self.parentZ = -1

    def find_peaks(self, height):   #peak picking height is percentage cutoff of the maximum intensity in spectra
        peaks, _ = signal.find_peaks(self.intensities, height=height)
        self.peaks = peaks

    def load_calibrated_MS_spectra_from_txt(self, filename, parentZ): #import of x-y MS spectrum from txt data
        try:
            data = np.genfromtxt(filename, delimiter="\t")
            x = data[:, 0]
        except IndexError:
            data = np.genfromtxt(filename, delimiter=",")
            x = data[:, 0]
        y = data[:, 1]
        cs = CubicSpline(x, y)
        mz = self.mz
        ys = cs(mz)
        ys = np.array(ys)
        for i in range(len(mz)):
            if mz[i] <= x[0]:
                ys[i] = 0
        ys = ys / np.max(ys)
        self.parentZ = parentZ
        self.intensities = ys

    def peak_detection(self, isotopes_list, tol=1e-6, enter_mz_values=False):
        #function for detection of the isotopes from the supplemented list
        #isotopes_list - output from pyopenms coarse isotope calculator
        #tol - maximum deviation (in Th values) for the isotopes in theoretical and experimental data for matching
        #enter_mz_values - obsolete, for debugging purposes
        masses = self.mz[self.peaks]
        intensities = self.intensities[self.peaks]
        i = 0
        j = 0
        max_int = 0.0
        isotope_identified = list(isotopes_list * 0)
        if isotopes_list[0] > masses[-1]:
            return isotope_identified
        i = np.argmin(np.abs(masses - isotopes_list[0])) - 10
        while i < len(masses):
            if masses[i] > isotopes_list[j] + tol:
                j += 1
            if j >= len(isotopes_list):
                break
            if np.abs(masses[i] - isotopes_list[j]) < tol:
                isotope_identified[j] = 1.0
                max_int = max(max_int, intensities[i])
            i += 1
        if enter_mz_values:
            return (isotopes_list, isotope_identified, max_int)
        else:
            return isotope_identified


class IsoformTable:
    #class associated with the isoform table, main component for the analysis
    #contains the information about series of the fragments, and intervals of every series for every isoform
    #contains methods for assessment of which fragments for every series
    #once MS data is processed, every series is supplemented with the list of identified fragments based on the CP valeus

    def __init__(self, Isoforms):
        """

        :type Isoforms: list(I)
        """
        self.series = []    #list of series S presented in IT
        self.table_of_fragments = []  # table with theoretical fragment intervals for every isoform  #

        #in three following dicts keys are series
        self.CP_values = {} # dict with CP values for fragments. values are lists of lists (2D matrix) with CP value of each fragment of each charge state
        self.Q = {} # dict with indicators of the presence of specific fragment in the spectrum. values are associated with the specific fragment indicator (1 if present, 0 if not)#        #
        self.gamma = {} #information on the gamma values for the fragments. values are lists with elements associated with fragment in series

        self.isoforms_names = [] #names of isoforms presented in IT
        self.ticks_array = None # x-ticks for bayesian plot
        self.report_strings = [] # content for table in bayesian plot window


        for isoform in Isoforms:    #unpacking of fragment series data for every isoform and combining them into the one set
            self.isoforms_names.append(isoform.name)
            for serie in isoform.S_set:
                in_series = -1
                for i in range(len(self.series)):
                    if self.series[i] == serie:
                        in_series = i
                if in_series != -1:
                    self.table_of_fragments[in_series][isoform.name] = isoform.S_set[serie]
                else:
                    self.series.append(serie)
                    self.table_of_fragments.append({isoform.name: isoform.S_set[serie]})

        for k in range(len(self.series)):
            self.series[k].generate_sequences()

        self.dict_of_series = {} #re-packing information of the table of fragments in the dict format

        for ind in range(len(self.series)):
            self.dict_of_series[self.series[ind]] = self.table_of_fragments[ind]

        for isoform in self.isoforms_names: #filling gaps in table with empty intervals
            for serie in self.series:
                if isoform not in self.dict_of_series[serie]:
                    self.dict_of_series[serie][isoform] = [(0, 1)]

    def continuity(self, isotopes):
        #supplemntal function for calculation of the longest isotopic series
        interim_res = 0.0
        res = 0.0
        for k in isotopes:
            if k == 1:
                interim_res += 1.0
            if k == 0:
                res = max(res, interim_res)
                interim_res = 0.0
        res = max(res, interim_res)
        return res / len(isotopes)


    def annotate_fragment(self, ms: InternalMassSpectrum, mass_accuracy, maxZ):
        #function that takes mass spectrum (in the internal format) and perform annotation of peaks and calcualtion of CP values
        serie: S
        for s in range(len(self.series)):
            serie = self.series[s]
            number_of_fragments = len(serie.seq_of_fragmenting_unit)
            table_of_identified_fragments = np.zeros((ms.parentZ, number_of_fragments))
            for z in range(1, maxZ):
                for k in range(1, number_of_fragments):
                    b = serie.fragment_mass_calculator(k, z)
                    b = remove_isotopes(b)
                    identified_isotopes = ms.peak_detection(b[:, 0], mass_accuracy)
                    table_of_identified_fragments[z][k] = self.continuity(identified_isotopes)
            self.CP_values[self.series[s]] = table_of_identified_fragments

    def evaluation_function(self, serie, CP2=1.0, CP1=0.25):
        #functon that uses CP values with CP1 and CP2 cutoffs to calculate indicators for presence of fragments and write them down in Q
        result = []
        for k in range(len(self.CP_values[serie][0])):
            vec_x = self.CP_values[serie][:, k]
            res = 0
            for x in vec_x:
                if x >= CP1:
                    res += x
            if res >= CP2:
                result.append(1)
            else:
                result.append(0)
        self.Q[serie] = result

    def gamma_compute(self, serie, N=5): #computation of gamma function
        q = self.Q[serie]
        gamma = moving_average(q, N)
        self.gamma[serie] = gamma
        for isoform in self.isoforms_names:
            for serie in self.series:
                if isoform not in self.dict_of_series[serie]:
                    self.dict_of_series[serie][isoform] = [(0, 1)]

    def compute_Pmatrix(self):
        #creation of the dict P for storage of conditional probabilities matrixes for every series
        #P has series as keys and values as matrixes
        #every matrix P has rows associated with isoforms, and columns with the fragment from the corresponsing series
        #these matrixes are filled with 1 and 0 values; 1 if the fragment can be generated by the isoform and 0 for opposite
        self.P = {}
        for serie in self.series:
            self.P[serie] = {}
            fragments = np.arange(len(serie.seq_of_fragmenting_unit))
            for isoform in self.dict_of_series[serie]:
                res = []
                for fragment in fragments:
                    if interval_check(fragment, self.dict_of_series[serie][isoform]):
                        res.append(1.0)
                        #res.append(0.5 + alpha / 2.0)
                    else:
                        res.append(0.0)
                        #res.append(0.5 - alpha / 2.0)
                self.P[serie][isoform] = np.array(res)

    def compute_Pt_matrix(self, alpha, consider_gamma=False, beta=0.05):
        #computation of the conditional probability matrixes dict Pt based on the P template using alpha, beta, and gamma values
        self.Pt = {}
        for serie in self.series:
            current_P_dict_fragment = self.P[serie]
            Q = self.Q[serie]
            if consider_gamma:
                gamma = self.gamma[serie]
            P = []
            for isoform_name in self.isoforms_names:
                P.append(current_P_dict_fragment[isoform_name])
            Pt = P.copy()
            P = np.array(P)
            Pt = np.array(Pt)
            for i in range(len(Q) - 1):
                n = len(P[:, i])
                if Q[i] == 0:
                    if not consider_gamma:
                        Pt[:, i] = P[:, i] * 0.0 + 1.0
                    else:
                        Pt[:, i] = 0.5 - gamma[i] * beta * (P[:, i] - 1 / 2)
                else:
                    Pt[:, i] = 0.5 + alpha * (P[:, i] - 1 / 2)
                Pt[:, i] = Pt[:, i] / np.sum(Pt[:, i])
            self.Pt[serie] = Pt

    def compute_probabilies(self):  #bayesian calculator
        r = []
        m = len(self.isoforms_names)
        r.append(np.zeros(m) + 1.0 / m)
        self.ticks_array = []
        self.report_strings = []
        for serie in self.series:
            Pt = self.Pt[serie]
            for i in range(len(Pt[0]) - 1):
                rnew = r[-1] * Pt[:, i]
                rnew = rnew / np.sum(rnew)
                r.append(rnew)
                rstring = [str(len(r)), str(serie), str(i)]
                for k in range(len(rnew)):
                    rstring.append("{:.3f}".format(rnew[k]))
                self.ticks_array.append(str(serie) + ",k=" + str(i))
                self.report_strings.append("\t".join(rstring))
        self.r_matrix = np.array(r)
        self.final_probs = [self.isoforms_names, rnew]

def main():
    return 0

if __name__ == "__main__":
    main()

