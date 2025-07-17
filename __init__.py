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


def sorted_range(seq):
    return sorted(range(len(seq)), key=seq.__getitem__)


def dicts_uniter(dict1, dict2):
    """
    :type dict1: dict
    :type dict2: dict
    """
    res_dict = dict1.copy()
    for k in dict2:
        if k in res_dict:
            res_dict[k] += dict2[k]
        else:
            res_dict[k] = dict2[k]
    return res_dict


def moving_average(vec, N=3):
    ma_vec = np.concatenate((vec[::-1],vec,vec[::-1]))
    ma_vec.astype(float)
    res = np.convolve(ma_vec, np.ones(N)/N, mode='valid')
    lvec = len(vec)
    return res[lvec:2*lvec-1]



def remove_isotopes(iso, cutoff=0.2):
    mx = 0
    for i in range(len(iso)):
        if iso[i][1] > mx:
            mx = iso[i][1]
    rem = []
    for i in range(len(iso)):
        if iso[i][1] > mx * cutoff:
            rem.append([iso[i][0], iso[i][1]])
    return np.array(rem)


def interval_overlapper(set_intervals_1, set_intervals_2):
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
            if x >= k[0] and x < k[1]:
                inside = True
        if inside:
            intervals.append([all_indexes[i], all_indexes[i + 1]])
    # print(intervals)
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
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        self.length = len(self.seq)

    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name


class S:
    def __init__(self, f_type, fragmenting_ppt, intact_ppts, f_name):
        """

        :type intact_ppts: dict
        """
        self.f_type = f_type
        self.f_name = f_name
        self.fragmenting_ppt = fragmenting_ppt
        self.intact_ppts = intact_ppts

    def __str__(self):
        if len(self.intact_ppts) == 0:
            return self.f_name + "(" + str(self.fragmenting_ppt) + ")"
        else:
            return str(self.intact_ppts) + "+" + self.f_name + "(" + str(self.fragmenting_ppt) + ")"

    def __repr__(self):
        return str(self.intact_ppts) + "+" + self.f_name + "(" + str(self.fragmenting_ppt) + ")"

    def __hash__(self):
        if len(self.intact_ppts) == 0:
            for_calc = (self.f_name, self.fragmenting_ppt, 0)
        else:
            intact_ppts_keys = []
            for k in self.intact_ppts:
                intact_ppts_keys.append(str(k) + "_" + str(self.intact_ppts[k]))
            intact_ppts_keys.sort()
            for_calc = (self.f_name, self.fragmenting_ppt, (tuple(intact_ppts_keys)))
        # print("hash of",for_calc, "is", hash(for_calc))
        return hash(for_calc)

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def generate_sequences(self):
        self.seqs_of_intact_ppts = []
        self.seq_of_fragmeting_unit = self.fragmenting_ppt.seq
        for ppt in self.intact_ppts:
            for k in range(self.intact_ppts[ppt]):
                self.seqs_of_intact_ppts.append(ppt.seq)

    def fragment_seq_generator(self, k):
        fr_unit_length = len(self.seq_of_fragmeting_unit)
        if self.f_type == "N":
            return "".join(self.seqs_of_intact_ppts) + self.seq_of_fragmeting_unit[:k]
        elif self.f_type == "C":
            return "".join(self.seqs_of_intact_ppts) + self.seq_of_fragmeting_unit[fr_unit_length - k:fr_unit_length]

    def fragment_mass_calculator(self, k, z, max_num_of_isotopes=20):
        result = []
        seq = pyopenms.AASequence.fromString(self.fragment_seq_generator(k))
        frag = seq.getFormula(__fragment_mass_calculation__[self.f_name], z)
        isotopes = frag.getIsotopeDistribution(pyopenms.CoarseIsotopePatternGenerator(max_num_of_isotopes))
        for iso in isotopes.getContainer():
            result.append([iso.getMZ() / z, iso.getIntensity()])
        return result


class I:
    V: list[Ppt]

    def __init__(self, name, init_formula, ppts):
        self.name = name
        self.ppts = ppts
        self.init_formula = init_formula
        self.V = [""]
        self.E = []
        self.parsing(init_formula[0:len(init_formula)], 0)
        self.S_set = {}

    def generate_dict_with_all_ppts_in_I(self):
        all_units = {}
        for k in self.V:
            if k not in all_units:
                all_units[k] = 1
            else:
                all_units[k] += 1
        self.all_units = all_units

    def compute_reachebility(self):
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
        self.reachebility = havebeen_matrix

    def generate_set_of_complimentary_intact_ppts(self, d, v):
        result = self.all_units.copy()
        result[v] -= 1
        for k in d:
            result[k] -= d[k]
        return result

    def generate_series_of_fragment(self, fragname):
        self.generate_dict_with_all_ppts_in_I()
        self.compute_reachebility()
        fragtype = __fragment_types__[fragname]
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
            if fragtype == "N":
                prev_i = 0
                k = 0
                new_S_F = self.V[ind_v]
                new_S_ftype = fragtype
                new_S_fname = fragname
                d = {}
                while k < len(conj_ne_s):
                    new_S_indexes = (prev_i, conj_ne_s[k])
                    prev_i = conj_ne_s[k]
                    new_S_D = d.copy()
                    s = S("N", new_S_F, new_S_D, fragname)
                    S_set.append((s, (new_S_indexes)))
                    d = dicts_uniter(d, self.reachebility[conj_e_s[k]])
                    k += 1
                new_S_indexes = (prev_i, length_v)
                new_S_D = d.copy()
                s = S("N", new_S_F, new_S_D, fragname)
                S_set.append((s, (new_S_indexes)))
            if fragtype == "C":
                prev_i = 0
                k = 0
                new_S_F = self.V[ind_v]
                new_S_ftype = fragtype
                new_S_fname = fragname
                d = {}
                while k < len(conj_ne_s):
                    new_S_indexes = (length_v - conj_ne_s[k], length_v - prev_i)
                    prev_i = conj_ne_s[k]
                    new_S_D = d.copy()
                    complimentary_S_D = self.generate_set_of_complimentary_intact_ppts(new_S_D, self.V[ind_v])
                    d = dicts_uniter(d, self.reachebility[conj_e_s[k]])
                    s = S("C", new_S_F, complimentary_S_D, fragname)
                    S_set.append((s, (new_S_indexes)))
                    k += 1
                new_S_indexes = (0, length_v - prev_i)
                new_S_D = d.copy()
                complimentary_S_D = self.generate_set_of_complimentary_intact_ppts(new_S_D, self.V[ind_v])
                s = S("C", new_S_F, complimentary_S_D, fragname)
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
                new_intervals = interval_overlapper(final_S_set[real_value], [q[1]])
                final_S_set[real_value] = new_intervals
        self.S_set = final_S_set

    def parsing(self, cur_formula, parent):
        """

        :type cur_formula: str
        """
        # print(cur_formula, parent)
        commas = []
        q = 0
        for i in range(len(cur_formula)):
            if cur_formula[i] == "(":
                q += 1
            elif cur_formula[i] == ")":
                q -= 1
            if (cur_formula[i] == ",") and (q == 0):
                commas.append(i)
        # print("COMMAS",commas)
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
    def __init__(self, high_mz, bin_size):
        self.peaks = None
        self.mz = np.arange(0, int(high_mz / bin_size), dtype=float) * bin_size
        self.intensities = np.zeros(len(self.mz))
        self.parentZ = -1

    def find_peaks(self, height):
        peaks, _ = signal.find_peaks(self.intensities, height=height)
        self.peaks = peaks

    def load_calibrated_MS_spectra_from_txt(self, filename, parentZ):
        try:
            data = np.genfromtxt(filename, delimiter="\t")
            x = data[:, 0]
        except IndexError:
            data = np.genfromtxt(filename, delimiter=",")
            x = data[:, 0]
        y = data[:, 1]
        #plt.plot(x,y)
        #plt.show()
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
        masses = self.mz[self.peaks]
        intensities = self.intensities[self.peaks]
        i = 0
        j = 0
        max_int = 0.0
        isotope_identified = list(isotopes_list * 0)
        #print(masses)
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
                    #print(max_int, intensities[i])
                    max_int = max(max_int, intensities[i])
            i += 1
        #print(isotopes_list, isotope_identified)
        if enter_mz_values:
            return (isotopes_list, isotope_identified, max_int)
        else:
            return isotope_identified

    def mass_calibrator(self, series, mass_accuracy=0.01, max_Z=5, cutoff=0.6, intensity_cutoff=0.4, iterations=2):
        theoretical_mz = []
        experimental_mz = []
        ms = self
        for s in range(len(series)):
            serie = series[s]
            number_of_fragments = len(serie.seq_of_fragmeting_unit)
            table_of_identified_fragments = np.zeros((ms.parentZ, number_of_fragments))
            for z in range(1, max_Z):
                for k in range(1, number_of_fragments):
                    b = serie.fragment_mass_calculator(k, z)
                    b = remove_isotopes(b, cutoff)
                    identified_isotopes = ms.peak_detection(b[:, 0], mass_accuracy, enter_mz_values=True)
                    for i in range(len(identified_isotopes)):
                        if identified_isotopes[i] != 0:
                            if identified_isotopes[i][1] > intensity_cutoff:
                                theoretical_mz.append(b[:,0][i])
                                experimental_mz.append(identified_isotopes[i][0])
                                #print(serie, k, z, b[:,0][i], identified_isotopes[i][0])
        theormz = theoretical_mz
        devmz = np.array(theoretical_mz)-np.array(experimental_mz)
        fig, ax = plt.subplots(2,1)
        for k in range(iterations):
            param = np.polyfit(theormz, devmz, 1)
            polynom = np.poly1d(param)
            deviations_from_approx = devmz-polynom(devmz)
            ax[0].plot(theormz,devmz,"o")
            ax[0].plot(ms.mz,polynom(ms.mz),"-")
            #print(np.std(deviations_from_approx))
            new_theor = []
            new_devmz = []
            for i in range(len(theormz)):
                if np.abs(deviations_from_approx[i]) < 2*np.std(deviations_from_approx):
                    new_theor.append(theormz[i])
                    new_devmz.append(devmz[i])
            theormz = np.array(new_theor)
            devmz = np.array(new_devmz)
        param = np.polyfit(theormz, devmz, 1)
        polynom = np.poly1d(param)
        new_mz = ms.mz+polynom(ms.mz)
        ax[1].plot(ms.mz, ms.intensities)
        ax[0].set_xlim(min(ms.mz),max(ms.mz))
        cs = CubicSpline(new_mz, ms.intensities)
        ys = cs(self.mz)
        ys = np.array(ys)
        ys = ys / np.max(ys)
        self.intensities = ys
        ax[1].plot(self.mz, self.intensities)
        plt.show()

class IsoformTable:
    def __init__(self, Isoforms):
        """

        :type Isoforms: list(I)
        """
        self.ticks_array = None
        self.series = []
        self.table_of_fragments = []
        self.Q = {}
        self.ksi = {}
        self.identified_fragments = {}
        self.isoforms_names = []
        self.report_strings = []
        for isoform in Isoforms:
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
        self.dict_of_series = {}
        for ind in range(len(self.series)):
            self.dict_of_series[self.series[ind]] = self.table_of_fragments[ind]

        for isoform in self.isoforms_names:
            for serie in self.series:
                if isoform not in self.dict_of_series[serie]:
                    self.dict_of_series[serie][isoform] = [(0,1)]

    def continuity(self, isotopes):
        interim_res = 0.0
        res = 0.0
        for k in isotopes:
            if k == 1:
                interim_res += 1.0
            if k == 0:
                res = max(res, interim_res)
                interim_res = 0.0
        res = max(res, interim_res)
        #print("cont", isotopes, res, res / len(isotopes))
        return res / len(isotopes)


    #def correct_mass(self, ms:InternalMassSpectrum, mass_accuracy):

    def annotate_fragment(self, ms: InternalMassSpectrum, mass_accuracy, maxZ):
        serie: S
        for s in range(len(self.series)):
            serie = self.series[s]
            number_of_fragments = len(serie.seq_of_fragmeting_unit)
            table_of_identified_fragments = np.zeros((ms.parentZ, number_of_fragments))
            for z in range(1, maxZ):
                for k in range(1, number_of_fragments):
                    b = serie.fragment_mass_calculator(k, z)
                    b = remove_isotopes(b)
                    identified_isotopes = ms.peak_detection(b[:, 0], mass_accuracy)
                    #print(k, z, identified_isotopes)
                    # c = remove_intensities(b)
                    # print(str(serie),z, k, self.continuity(identified_isotopes))
                    table_of_identified_fragments[z][k] = self.continuity(identified_isotopes)
            self.identified_fragments[self.series[s]] = table_of_identified_fragments




        '''   
        plt.plot(theoretical_mz, np.array(theoretical_mz)-np.array(experimental_mz), "ro")
        x_ref = np.arange(min(ms.mz), max(ms.mz))
        
        y_ref = polynom(x_ref)
        plt.plot(x_ref, y_ref,"b")
        plt.xlim(min(ms.mz), max(ms.mz))
        plt.ylim(-0.003, 0.003)

        plt.show()
        '''

    def evaluation_function_old(self, serie, cutoff=1.0, sgo=0.25):
        result = []
        for k in range(len(self.identified_fragments[serie][0])):
            vec_x = self.identified_fragments[serie][:, k]
            interim_res = 0
            res = 0
            for x in vec_x:
                if x >= sgo:
                    interim_res += x
                else:
                    res = max(interim_res, res)
                    interim_res = 0
            if res >= cutoff:
                result.append(1)
            else:
                result.append(0)
        self.Q[serie] = result

    def evaluation_function(self, serie, cutoff=1.0, sgo=0.25):
        result = []
        for k in range(len(self.identified_fragments[serie][0])):
            vec_x = self.identified_fragments[serie][:, k]
            #print(vec_x)
            res = 0
            for x in vec_x:
                if x >= sgo:
                    res += x
            if res >= cutoff:
                result.append(1)
            else:
                result.append(0)
        self.Q[serie] = result

    def ksi_compute(self, serie, N=5):
        q = self.Q[serie]
        ksi = moving_average(q, N)
        self.ksi[serie] = ksi

    #    def interval_deconvolution(self):

        for isoform in self.isoforms_names:
            for serie in self.series:
                if isoform not in self.dict_of_series[serie]:
                    self.dict_of_series[serie][isoform] = [(0,1)]

    def compute_Pmatrix(self):
        self.P = {}
        all_isoforms = self.isoforms_names
        #print(all_isoforms)
        #print(self.table_of_fragments)
        for serie in self.series:
            self.P[serie] = {}
            fragments = np.arange(len(serie.seq_of_fragmeting_unit))
            print(self.dict_of_series)
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



    def compute_Pt_matrix(self, alpha,consider_ksi=False, beta=0.05):
        self.Pt = {}
        for serie in self.series:
            print(serie, self.P[serie])
            current_P_dict_fragment = self.P[serie]
            Q = self.Q[serie]
            if consider_ksi:
                ksi = self.ksi[serie]
            P = []
            for isoform_name in self.isoforms_names:
                P.append(current_P_dict_fragment[isoform_name])
            Pt = P.copy()
            P = np.array(P)
            Pt = np.array(Pt)
            for i in range(len(Q)-1):
                n = len(P[:, i])
                if Q[i] == 0:
                    if not consider_ksi:
                        Pt[:, i] = P[:, i] * 0.0 + 1.0
                    else:
                        Pt[:, i] = 0.5 - ksi[i]*beta*(P[:, i] - 1/2)
                else:
                    Pt[:, i] = 0.5+alpha*(P[:, i] - 1/2)
                Pt[:, i] = Pt[:, i] / np.sum(Pt[:, i])
            self.Pt[serie] = Pt
            #print(alpha)
            #print(str(serie))
            #print("P", P)
            #print("Pt", Pt)
            #print("ksi", ksi)
            #print(len("ksi"), ksi)

    def compute_probabilies(self):
        r = []
        m = len(self.isoforms_names)
        r.append(np.zeros(m) + 1.0 / m)
        self.ticks_array = []
        for serie in self.series:
            Pt = self.Pt[serie]
            for i in range(len(Pt[0]) - 1):
                rnew = r[-1] * Pt[:, i]
                rnew = rnew / np.sum(rnew)
                r.append(rnew)
                rstring = [str(len(r)),str(serie), str(i)]
                for k in range(len(rnew)):
                    rstring.append("{:.3f}".format(rnew[k]))
                self.ticks_array.append(str(serie)+",k="+str(i))
                self.report_strings.append("\t".join(rstring))
        self.r_matrix = np.array(r)
        self.final_probs = [self.isoforms_names, rnew]

    def fragment_plotter(self):
        total_number_of_series = len(self.series)
        rows_number = int(round(np.sqrt(total_number_of_series)))
        axs = []
        fig = plt.figure()
        nrs = rows_number
        ncs = total_number_of_series // rows_number + 1
        for k in range(total_number_of_series):
            ax = fig.add_subplot(ncs, nrs, k + 1)
            ax.imshow(self.identified_fragments[self.series[k]], aspect="equal")
            ax.set_title(str(self.series[k]))
            if self.series[k].f_type == "C":
                ax.invert_xaxis()
            axs.append(ax)
        # plt.tight_layout()
        plt.show()

    def bayesian_show(self):
        for k in range(len(self.r_matrix[0])):
            plt.plot(self.r_matrix[:,k]+k*0.001, "o--",label=self.isoforms_names[k],markersize=2)
        plt.ylim(0,1)
        plt.legend()
        plt.show()





def main_old():
    path_to_file = "k6k6k48_prox_1273_etd5_2mtorr.txt"
    operating_mass_spectra = InternalMassSpectrum(3000, 0.0001)
    operating_mass_spectra.load_calibrated_MS_spectra_from_txt(path_to_file, 20)
    operating_mass_spectra.find_peaks(height=0.002)
    Ubq = Ppt("Ub", "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
    ppts = {"Ub": Ubq}
    test_I = I("test_isoform1", "Ub,48(Ub, 48(Ub))", ppts)
    test_I.generate_series_of_fragment("c")
    test_I.generate_series_of_fragment("z")
    IT = IsoformTable([test_I])
    IT.annotate_fragment(operating_mass_spectra, mass_accuracy=0.003)
    for s in IT.series:
        IT.evaluation_function(s, cutoff=1.0, sgo=0.25)
    IT.compute_Pmatrix(alpha=0.05)
    IT.compute_Pt_matrix()
    IT.compute_probabilies()
    # IT.fragment_plotter()
    return 0


def main():
    path_to_file = "ub3_k6_k48br.txt"
    operating_mass_spectra = InternalMassSpectrum(2000, 0.001)
    operating_mass_spectra.load_calibrated_MS_spectra_from_txt(path_to_file, parentZ=25)
    operating_mass_spectra.find_peaks(height=0.002)
    Ubq = Ppt("Ub", "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
    H2 = Ppt("Ub", "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
    ppts = {"Ub": Ubq}
    #I_K6 = I("K6", "Ub,6(Ub, 6(Ub))", ppts)
    I_array = []
    I_array.append(I("k48k63", "Ub,48(Ub), 63(Ub)", ppts))
    I_array.append(I("k11k48", "Ub,11(Ub), 48(Ub)", ppts))
    I_array.append(I("k6k48", "Ub,6(Ub), 48(Ub)", ppts))
#    I_array.append(I("k48linear", "Ub,48(Ub,48(Ub))", ppts))
    I_array.append(I("k6linear", "Ub,6(Ub,6(Ub))", ppts))
#    I_array.append(I("k63linear", "Ub,63(Ub,63(Ub))", ppts))
    for i in range(len(I_array)):
        I_array[i].generate_series_of_fragment("c")
        I_array[i].generate_series_of_fragment("z")
    IT = IsoformTable(I_array)
    #for i in range(len(IT.series)):
    #    print(i, IT.series[i])
    #print(IT.table_of_fragments)
    IT.series[1].seqs_of_intact_ppts[0] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGRQLEDGRTLSDYNIQRESTLHLVLRLRGG"
    IT.series[2].seqs_of_intact_ppts[0] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGRQLEDGRTLSDYNIQRESTLHLVLRLRGG"
    IT.series[2].seqs_of_intact_ppts[1] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGRQLEDGRTLSDYNIQRESTLHLVLRLRGG"
    IT.series[3].seqs_of_intact_ppts[0] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGRQLEDGRTLSDYNIQRESTLHLVLRLRGG"
    IT.series[3].seqs_of_intact_ppts[1] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGRQLEDGRTLSDYNIQRESTLHLVLRLRGG"
    IT.series[4].seqs_of_intact_ppts[0] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGRQLEDGRTLSDYNIQRESTLHLVLRLRGG"
    IT.annotate_fragment(operating_mass_spectra, mass_accuracy=0.003)
#    IT.series[]
    for s in IT.series:
        IT.evaluation_function(s, cutoff=1.0, sgo=0.25)
    IT.fragment_plotter()
    IT.compute_Pmatrix(alpha=0.05)
    IT.compute_Pt_matrix()
    IT.compute_probabilies()
    IT.bayesian_show()
    #print(IT.Q)

#    operating_mass_spectra.mass_calibrator(IT.series, mass_accuracy=0.003, max_Z=10, cutoff=0.1, intensity_cutoff=0.2, iterations=5)
    # IT.fragment_plotter()
    return 0

def main_H2():
    path_to_file = "H2B_Ub.txt"
    operating_mass_spectra = InternalMassSpectrum(2000, 0.001)
    operating_mass_spectra.load_calibrated_MS_spectra_from_txt(path_to_file, parentZ=20)
    operating_mass_spectra.find_peaks(height=0.002)
    H2B_seq = "AKSAPAPKKGSKKAVTKTQKKDGKKRRKTRKESYAIYVYKVLKQVHPDTGISSKAMSIMNSFVNDVFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSAK"
    H2B = Ppt("H2B", H2B_seq)
    Ub = Ppt("Ub", "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")

    ppts = {"H2B": H2B, "Ub": Ub}
    I_list = []
    for i in range(len(H2B_seq)):
        if H2B_seq[i] == "K":
            I_list.append(I("K"+str(i+1),"H2B,"+str(i+1)+"(Ub)", ppts))
    for i in range(len(I_list)):
        I_list[i].generate_series_of_fragment("c")
        I_list[i].generate_series_of_fragment("z")
    IT = IsoformTable(I_list)
    for i in range(len(IT.series)):
        print(i, IT.series[i])
    #IT.series[1].seqs_of_intact_ppts[0] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGRQLEDGRTLSDYNIQRESTLHLVLRLRGG"
    #IT.series[2].seqs_of_intact_ppts[0] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGRQLEDGRTLSDYNIQRESTLHLVLRLRGG"
    #IT.series[2].seqs_of_intact_ppts[1] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGRQLEDGRTLSDYNIQRESTLHLVLRLRGG"
    #IT.series[3].seqs_of_intact_ppts[0] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGRQLEDGRTLSDYNIQRESTLHLVLRLRGG"
    #IT.series[3].seqs_of_intact_ppts[1] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGRQLEDGRTLSDYNIQRESTLHLVLRLRGG"
    #IT.series[4].seqs_of_intact_ppts[0] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGRQLEDGRTLSDYNIQRESTLHLVLRLRGG"
    IT.annotate_fragment(operating_mass_spectra, mass_accuracy=0.003)
#    IT.series[]
    for s in IT.series:
        IT.evaluation_function(s, cutoff=1.5, sgo=0.25)
    IT.fragment_plotter()
    IT.compute_Pmatrix(alpha=0.05)
    IT.compute_Pt_matrix()
    IT.compute_probabilies()
    IT.bayesian_show()
    print(IT.r_matrix[-1,:])
    #print(IT.Q)

#    operating_mass_spectra.mass_calibrator(IT.series, mass_accuracy=0.003, max_Z=10, cutoff=0.1, intensity_cutoff=0.2, iterations=5)
    # IT.fragment_plotter()
    return 0

def main_new():
    path_to_file = "yuh1k6.txt"
    operating_mass_spectra = InternalMassSpectrum(2000, 0.001)
    operating_mass_spectra.load_calibrated_MS_spectra_from_txt(path_to_file, parentZ=25)
    operating_mass_spectra.find_peaks(height=0.002)
    UbK6R = Ppt("UbK6R", "MQIFVRTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
    Ub = Ppt("Ub", "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
    Yuh = Ppt("Yuh1", "DQFVLNVIKENVQTFSTGQSEAPEATAAKKAHYITYVEENGGIFEL")
    ppts = {"Yuh1": Yuh, "Ub": Ub, "UbK6R":UbK6R}
    I_list = []
    I_list.append(I("Prot","Yuh1,9(Ub,6(UbK6R))", ppts))
    for i in range(len(I_list)):
        I_list[i].generate_series_of_fragment("c")
        I_list[i].generate_series_of_fragment("z")
    IT = IsoformTable(I_list)
    for i in range(len(IT.series)):
        print(i, IT.series[i])
    IT.annotate_fragment(operating_mass_spectra, mass_accuracy=0.003)
#    IT.series[]
    for s in IT.series:
        IT.evaluation_function(s, cutoff=1.5, sgo=0.25)
    IT.fragment_plotter()

def main_k48():
    path_to_file = "yuh1k48.txt"
    operating_mass_spectra = InternalMassSpectrum(2000, 0.001)
    operating_mass_spectra.load_calibrated_MS_spectra_from_txt(path_to_file, parentZ=25)
    operating_mass_spectra.find_peaks(height=0.002)
    UbK6R = Ppt("UbK48R", "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGRQLEDGRTLSDYNIQKESTLHLVLRLRGG")
    Ub = Ppt("Ub", "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
    Yuh = Ppt("Yuh1", "DQFVLNVIKENVQTFSTGQSEAPEATAAKKAHYITYVEENGGIFEL")
    ppts = {"Yuh1": Yuh, "Ub": Ub, "UbK48R":UbK6R}
    I_list = []
    I_list.append(I("Prot","Yuh1,9(Ub,48(UbK48R))", ppts))
    for i in range(len(I_list)):
        I_list[i].generate_series_of_fragment("c")
        I_list[i].generate_series_of_fragment("z")
    IT = IsoformTable(I_list)
    for i in range(len(IT.series)):
        print(i, IT.series[i])
    IT.annotate_fragment(operating_mass_spectra, mass_accuracy=0.003)
#    IT.series[]
    for s in IT.series:
        IT.evaluation_function(s, cutoff=1.5, sgo=0.25)
    IT.fragment_plotter()

def test_main():
    IT= IsoformTable()
    x= [1.0, 1.0]
    print(IT.continuity(x))

def test():
    Q = np.array([1.0,0.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0])
    print(moving_average(Q,5))



if __name__ == "__main__":
    test()

# X,48(Ub,6(Ub,48(Ub),63(Ub))),63(Ub,29(Ub))

'''
    def penalty_compute(self, serie, n):
        Q_array = self.Q[serie]
        for k in range(len(Q_array)):
            Q_m = np.concatenate((Q_array[::-1], Q_array, Q_array[::-1]))


            vec_x = self.identified_fragments[serie][:, k]
            print(vec_x)
            res = 0
            for x in vec_x:
                if x >= sgo:
                    res += x
            if res >= cutoff:
                result.append(1)
            else:
                result.append(0)
        self.Q[serie] = result

'''