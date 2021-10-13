import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
#####################################
class MGF:
    def __init__(self):
        
        self.spectra = []
        self.filename = None

        #filter attributes
        self.filter_peaks = None
        self.pos_spectra = []
        self.neg_spectra = []

        #difference attributes (Tree)
        #self.diff_root = None
        #self.diff_count = 0
        
        #difference attributes (List)
        self.diff_list = []
        self.unique_diffs = []
        self.diff_counts = []
        self.diff_spectra = []


    def add_spectra(self, spectra):
        self.spectra += spectra

    def save(self, filename = False):
        """saves a dictionary generated from a.mgf file back into a .mgf file

        Keyword Arguments:
        mgf -- dictionary generated from a .mgf file
        filename -- filename to give the new .mgf file
        """
        if filename == False:
            filename_mgf = self.filename[:-3] + ("_PROCESSED.mgf")
        else:
            filename_mgf = str(filename) + ".mgf"
        #Write to new mgf file
        #Blank file
        with open(filename_mgf, "w") as f:
            f.write("")
        #Append each line to file
        with open(filename_mgf, "ab") as f:
            for entry in self.spectra:
                f.write(bytearray(b"BEGIN IONS\n"))
                f.write(bytearray("TITLE=index=" + entry.title + "\n", "utf-8"))
                f.write(bytearray("PEPMASS=" + entry.pepmass + "\n", "utf-8"))
                f.write(bytearray("CHARGE=" + entry.charge + "+\n", "utf-8"))
                f.write(bytearray("SCANS=" + entry.scans + "\n", "utf-8"))
                f.write(bytearray("RTINSECONDS=" + entry.rtinseconds + "\n", "utf-8"))
                for peak in entry.spectra:
                    f.write(bytearray(str(peak[0]) + " " + str(peak[1]) + "\n", "utf-8"))
                f.write(b"END IONS\n\n")

    def filter(self, mz_list, intensity_tol = 0.2, wobble_tol = 0.01, ms2_peaks = 10, checks = 1):
        """Takes a .mgf file converted into a dictionary by extract_mgf() and extracts all spectra with peaks of a desired m/z value and 
        intensity relative to the max intensity in the spectra
        
        Keyword Arguments:
        mgf_list -- dictionary generated from a .mgf file
        mz_list -- a list of desired m/z values to check spectra for
        intensity_tol -- intensity relative to the max peak in each spectra that a peak must be above to be           detected
        wobble_tol -- value above and below the specified m/z values peaks can be detected
        checks -- the number of peaks that a spectra must have with the specified m/z value to be detected, a value below 1 will require 1 or more peaks to be present and a value greater than the amount of specified peaks will require all peaks to be present
        """
        #Extract spectra with desired peaks
        if type(mz_list) is float or type(mz_list) is int:
            mz_list = [mz_list]
        if checks < 1:
            checks = 1
        else:
            pass
        if checks > len(mz_list):
            checks = len(mz_list)
        else:
            pass
        entries = 0
        neg_entries = 0
        new_mgf = []
        new_mgf_neg = []
        if len(self.spectra) >= ms2_peaks:
            for entry in self.spectra:
                if entry.check_peaks(mz_list = mz_list, intensity_tol = intensity_tol, wobble_tol = wobble_tol) >= checks:
                    new_mgf.append(entry)
                    entries += 1
                else:
                    new_mgf_neg.append(entry)
                    neg_entries += 1
        else:
            pass
        self.pos_spectra = MGF()
        self.pos_spectra.spectra = new_mgf
        self.pos_spectra.filename = self.filename[:-4] + "_pos.mgf"
        self.neg_spectra = MGF()
        self.neg_spectra.spectra = new_mgf_neg
        self.neg_spectra.filename = self.filename[:-4] + "_neg.mgf"
        ###DEBUG TRACKING
        print(len(new_mgf))
        print(len(new_mgf)/len(self.spectra))
        self.filter_peaks = mz_list

    def summary(self):
        string = ""
        string += "Spectra: {spectra}\n".format(spectra = len(self.spectra))
        av_pep = sum(float(x.pepmass) for x in self.spectra)/len(self.spectra) 
        string += "Average Precursor Mass: {av}\n".format(av = av_pep)
        av_charge = sum(float(x.charge) for x in self.spectra)/len(self.spectra)
        string += "Average Charge: {charge}\n\n".format(charge = av_charge)
        if self.filter_peaks is None:
            string += "Not yet filtered\n\n"
        else:
            string += "Filter: {filter}\n".format(filter = ", ".join(map(str, self.filter_peaks)))
            #string += Tolerance/Checks:
            string += "Positive Spectra: {pos_spec}\n".format(pos_spec = len(self.pos_spectra.spectra))
            string += "Negative Spectra: {neg_spec}\n".format(neg_spec = len(self.neg_spectra.spectra))
            string += "Ratio: {pos_percent:.2f}% positive\n".format(pos_percent = len(self.pos_spectra.spectra)/(len(self.neg_spectra.spectra) + len(self.pos_spectra.spectra)) * 100)
            p_m = sum(float(x.pepmass) for x in self.pos_spectra.spectra)/len(self.pos_spectra.spectra) 
            string += "Average Precursor Mass (Positive Spectra): {pos_mass:.2f}\n".format(pos_mass = p_m)
            pos_av_charge = sum(float(x.charge) for x in self.pos_spectra.spectra)/len(self.pos_spectra.spectra)
            string += "Average Charge (Positive (Spectra)): {charge}\n\n".format(charge = pos_av_charge)
            n_m = sum(float(x.pepmass) for x in self.neg_spectra.spectra)/len(self.neg_spectra.spectra) 
            string += "Average Precursor Mass (Negative Spectra): {neg_mass:.2f}\n".format(neg_mass = n_m)
            neg_av_charge = sum(float(x.charge) for x in self.neg_spectra.spectra)/len(self.neg_spectra.spectra)
            string += "Average Charge (Negative Spectra): {charge}\n\n".format(charge = neg_av_charge)
        if len(self.unique_diffs) != 0: 
            string += "Mass Differences: {diffs}\n".format(diffs = len(self.unique_diffs))
        else:
            string += "Differences not calculated"
        return string

    def diff(self, diff_tol = 0.001, intensity_tol = 0.2):
        return self.list_diff(self, diff_tol = diff_tol, intensity_tol = 0.2)

    def list_diff(self, diff_tol = 0.001, intensity_tol = 0.2, mass_tol = 0):
        """ calculates all mass differences within each spectra in the file

        Keyword Arguments:
        mgf_dict -- dictionary generated from a .mgf file
        diff_tol -- the difference in mz for two peaks to be considered different
        """
        #TODO: Frequency cutoff
        #TODO: Max spectra intensity threshold
        #iterate through mgf finding all differences and inserting them into diff_list
        #calculate mass differences
        print("start")
        print(datetime.now())
        for entry in self.spectra:
            valid_peaks = []
            spectra = entry.spectra
            int_thresh = float(max(spectra, key = lambda x:float(x[1]))[1]) * intensity_tol
            for pair in spectra:
                if float(pair[1]) > int_thresh and float(pair[0]) > mass_tol:
                    valid_peaks.append(pair[0])
                else:
                    pass                
            for row in range(len(valid_peaks) - 1):
                mass = valid_peaks[row]
                for column in range(row + 1, len(valid_peaks)):
                    self.diff_list.append(float(valid_peaks[column]) - float(mass))
        diff_list = self.diff_list
        if len(diff_list) == 0:
            print("No valid peaks")
            return
        print(diff_list[0])
        print("diffs calced")
        print(datetime.now())
        #sort diff_list
        diff_list.sort()
        print("sorted")
        print(datetime.now())
        print(len(diff_list))
        #iterate through diff list and make a count of all diffs
        current = diff_list[0]
        count = 0
        max_mass = 0
        for mass in diff_list:
            if mass > max_mass:
                max_mass = mass
            else:
                pass
            if abs(mass - current) <= diff_tol:
                count += 1
            else:
                self.diff_counts.append(count)
                self.unique_diffs.append(current)
                current = mass
                count = 1
        print("counted")
        #filter by fraction of max peak
        print(datetime.now())
    
    def find_diff(self, diff, tol = 0.001):
        found = []
        for current in range(len(self.unique_diffs)):
            if abs(self.unique_diffs[current] - diff) < tol:
                found.append([self.unique_diffs[current], self.diff_counts[current], self.diff_spectra[current]])
            else:
                pass
        return found

    def top_diffs(self, count = 100):
        index_tops = sorted(range(len(self.diff_counts)), key = lambda i: self.diff_counts[i], reverse = True)[0:count]
        return index_tops

    def precursur_shifts(self, target, tol = 0.01):
        if self.filter_peaks is None:
            print("MGF has not been filtered")
            return None
        #Check if filtered###############
        pos_precursors = []
        neg_precursors = []
        H = 1.00783
        for entry in self.pos_spectra.spectra:
            pre_mass = (float(entry.pepmass) * float(entry.charge)) - (float(entry.charge) * H)
            pos_precursors.append(pre_mass)
        for entry in self.neg_spectra.spectra:
            pre_mass = (float(entry.pepmass) * float(entry.charge)) - (float(entry.charge) * H)
            neg_precursors.append(pre_mass)
        #prec_grid = [[None]*len(neg_precursors) for i in range(len(pos_precursors))]
        print(len(pos_precursors))
        print(len(neg_precursors))
        #retrieve spectra
        hit_indexes = []
        counter = 0
        #base =     22304944
        #third =    16033
        #make lists in format [(index, difference), (index, difference), etc]
        pos_enum = [(ind, item) for ind, item in enumerate(pos_precursors)]
        neg_enum = [(ind, item) for ind, item in enumerate(neg_precursors)] 
        #sort lists by difference 
        pos_enum.sort(key = lambda x: x[1])
        neg_enum.sort(key = lambda x: x[1])
        #set starting values
        #positive pepmass list indexer
        i = 0
        #counter for tracking number of comparisons
        counter = 0
        #set initial lower bound index for comparisons
        bottom = 0
        while i < len(pos_enum):
            #j is negative pepmass list indexer
            #set as lower bound for searching
            j = bottom
            while j < len(neg_enum):
                #tracking comparisons
                counter += 1
                #calculate differences
                res = pos_enum[i][1] - neg_enum[j][1]
                #track minimum check value
                if res > target + tol:
                    #if result is greater than target range, shift the boundary to match it
                    bottom = j
                #if under target range break
                elif res < target - tol:
                    #if result is less than target range, stop comparisons and move to next item in positive list
                    break
                #record hits
                elif res > target - tol and res < target + tol:
                    hit_indexes.append((pos_enum[i][0], neg_enum[j][0]))
                j += 1
            i += 1
        print("Comparisons = {counter}".format(counter = counter))
        print(len(hit_indexes))
        pos_spectra_hits = []
        neg_spectra_hits = []
        for pair in hit_indexes:
            pos_spectra_hits.append(self.pos_spectra.spectra[pair[0]])
            neg_spectra_hits.append(self.neg_spectra.spectra[pair[1]])
        return pos_spectra_hits, neg_spectra_hits


class Spectra:
    def __init__(self, filename, title, pepmass, charge, scans, rt, spectra, rtsecs = True, pepmass_int = None):
        self.filename = filename
        self.title = title
        self.pepmass = pepmass
        self.charge = charge
        self.scans = scans
        self.rtinseconds = rt   
        self.spectra = spectra
        #track whether rt is in seconds or minutes
        self.rtsecs = rtsecs 
        #store pepmass intensity if present
        pepmass_int = pepmass_int

    def __str__(self):
        #return whatever i want it to print
        #construct string
        return "info stuff"

    def show(self):
        print("File: " + self.filename)
        print("Title: " + self.title)
        print("Pepmass: " + self.pepmass)
        print("Charge: " + self.charge)
        print("Scans: " + self.scans)
        print("RT in Seconds: " + self.rtinseconds)
        print("Spectra: ")
        for i in self.spectra:
            print(i[0] + " , " + i[1])

    def check_peaks(self, mz_list, intensity_tol, wobble_tol):
        """Checks whether a given spectra contains peaks of specified m/z values
        
        Keyword Arguments:
        mz_list -- a list of desired m/z values to check spectra for
        intensity_tol -- intensity relative to the max peak in each spectra that a peak must be above to be detected
        wobble_tol -- value above and below the specified m/z values peaks can be detected
        """
        if type(mz_list) is float or type(mz_list) is int:
            mz_list = [mz_list]
        #Check if list?
        mz_presence = []
        mz_occur = 0
        max_peak = 0
        for mz in mz_list:
            mz_peak = 0
            for peak in self.spectra:
                #check for max
                if len(peak) != 2:
                    print(peak)
                    break
                if peak[1] == "IONS":
                    return False
                if float(peak[1]) > max_peak:
                    max_peak = float(peak[1])
                else:
                    pass
                #Tolerance for peak wobble?
                if float(peak[0]) - mz < wobble_tol and float(peak[0]) - mz >= 0 and float(peak[1]) > mz_peak:
                    mz_peak = float(peak[1])
                else:
                    pass
            if mz_peak > max_peak*intensity_tol:
                mz_presence.append(True)
                mz_occur += 1
            else:
                mz_presence.append(False)
            return mz_occur

        def max_intensity():
            return max(self.spectra[1].spectra, key = lambda x:float(x[1]))[1]
###############################################################
def read_mgf(mgf_filename):
    #initialise MGF object
    new_mgf = MGF()
    new_mgf.filename = mgf_filename

    title = None
    rt = None
    rtsec = True
    charge = None
    scans = None
    pepmass = None
    pepmass_int = None
    
    mgf_file = open(mgf_filename,  "r")
    mgf_split = str(mgf_file.read()).split("END IONS")
    mgf_file.close()
    mgf_list = []
    spectra_list = []
#    begin = False
#    start = 0 
#    while begin is False:
#        if 
    for i in range(len(mgf_split) - 1):
        mgf_part = mgf_split[i].split("\n")
        if i != 0:
            mgf_part = mgf_part[2:]
        spec_check = False
        line = 1
        while spec_check is False:
            if "=" in mgf_part[line]:
                line += 1
            else:
                spec_check = True
        spectra = mgf_part[line:-1]
        if spectra[0][0] == "BEGIN":
            print(i)
        for j in range(len(spectra)):
            spectra[j] = spectra[j].split(" ") 
        header = mgf_part[1:line]
        for j in header:
            word = ""
            letter = 0
            while j[letter] != "=":
                word = word + j[letter]
                letter += 1
            if word == "TITLE":
                title = j[6:]
            elif word == "RTINSECONDS":
                rt = j[12:]
                rtsecs = True
            elif word == "RTINMINUTES":
                rt = j[12:]
                rtsecs = False
            elif word == "PEPMASS":
                pepmass = j[8:]
                if len(pepmass.split()) != 1:
                    split_pepmass = pepmass.split()
                    pepmass = split_pepmass[0]
                    pepmass_int = split_pepmass[1]
                else:
                    pass
            elif word == "CHARGE":
                charge = j[7:-1]
            elif word == "SCANS":
                scans = j[6:]
        new_spectra = Spectra(filename = mgf_filename, title = title, pepmass = pepmass, charge = charge, scans= scans, rt = rt, spectra = spectra, rtsecs = rtsecs, pepmass_int = pepmass_int)
        spectra_list.append(new_spectra)
    new_mgf.add_spectra(spectra_list)
    return new_mgf

def list_to_mgf(spectra_list, filename = "new"):
    new_mgf = MGF()
    new_mgf.filename = filename
    new_mgf.add_spectra(spectra_list)
    return new_mgf

def read_glycans(glycan_filename, mgf, new_filename = "_glycan_counts.txt"):
    with open(glycan_filename, "r") as f:
        glycan_list = []
        for line in f:
            glycan_list.append(line.split(" "))
    if len(mgf.unique_diffs) == 0:
        mgf.list_diff()
    else:
        pass
    output_string = ""
    for i in glycan_list:
        peaks = mgf.find_diff(float(i[2]))
        count = sum(x[1] for x in peaks)
        output_string += " ".join(i)[:-1] + " {}\n".format(str(count))
    with open(new_filename, "w") as f:
        f.write(output_string)