import matplotlib.pyplot as plt
import scipy as sci
from scipy.stats import skew
from scipy.stats import kurtosis

import numpy as np
import pandas as pd

def strip_iso(spectrum):
    """
    ##REDUNDANT??##

    Strips isotopic peaks from spectrum

    Keyword Arguments
    spectrum -- array of spectral peaks to be stripped of isotopic peaks

    Returns
    Array of spectral peaks stripped of identified isotopic peaks
    """
    #TODO
    #Maybe pick largest intensity instead of leftmost?
    #reconstruct as recursive to check every peak
    #must see at least X peaks (start with 3)
    #doubly charged testing: 4475
    stripped_spec = []
    base = spectrum[-1]
    diff_list = []
    charge_list = []
    for current_peak in range(len(spectrum) - 2, -1, -1):
        diff = round(float(base[0]) - float(spectrum[current_peak][0]), 1)
        diff_list.append(diff)
        if diff < 1:
            
            pass
            #check next peak
        elif diff == 1:
            #set current peak as base
            base = spectrum[current_peak]
            #check next peak
            pass
        else: 
            #record base peak as main peak
            stripped_spec = [base] + stripped_spec
            #set current peak as base peak
            base = spectrum[current_peak]
    diff_list = set(diff_list)
    print(diff_list)
    return stripped_spec

def find_peaks(spectrum, targets):
    """
    Checks a spectrum for peaks at target m/z values, and returns whether those peaks are present and the intensities of those peaks

    Keyword Arguments
    spectrum -- array of spectral peaks to be searched
    targets -- array of m/z values to check spectrum for
    """
    #check for presence of specific peaks
    #record intensities
    #[126.055]
    """ HexNAc-C2H6O3':126.055,'HexNAc-CH6O3':138.055,'HexNAc-C2H4O2':144.065,'Hex-H2O':145.0501,'HexNAc-2H2O':168.066,
    'HexNAc-H2O':186.076,'HexNAc':204.087,'NeuAc-H2O':274.092,'NeuAc':292.103,'HexHexNAc':366.140 """
    # form arrays to record hits and intensities
    hits = [0]*len(targets)
    intensities = [0]*len(targets)
    #round spectral m/z values and input target values to 1 decimal place
    round_targets = [round(float(x), 1) for x in targets]
    round_spec = [round(float(x[0]), 1) for x in spectrum]
    #check spectrum for target peaks and record presence and intensity of peaks
    for peak in range(len(round_targets)):
        if round_targets[peak] in round_spec:
            #record hit
            hits[peak] = 1
            #record intensity
            intensities[peak] = spectrum[round_spec.index(round_targets[peak])][1]
        else:
            pass
    return hits, intensities

def intensity_strip(spectrum, intensity):
    """
    strip spectrum of all peaks with intensity less than target proportion of highest intensity peak
    and return stripped spectrum

    Keyword Arguments
    spectrum -- array of spectral peaks to be stripped low intensity peaks
    intensity -- float value between 0 and 1 specifying proportion of maximum inntensity to strip peaks below 
    """
    return [x for x in spectrum if float(x[1]) >= float(float(max(spectrum, key = lambda x:float(x[1]))[1])*intensity)]

def yo_shift(spectrum, yo):
    """
    search spectra for a peak with mz value 17 below Y0 ion peak and return boolean value of its presence

    Keyword Arguments
    spectrum -- array of spectral peaks to search
    yo -- m/z value of Y0 peak
    """
    #if a peak exists at Y0 - 17 rounded to 2 decimal places, return True
    if (round(float(yo), 2) - 17) in [round(float(x[0]), 2) for x in spectrum]:
        return True
    else:
        return False

def check_ratio(spectrum):
    """
    checks ratio of the intensities of oxonium ion peaks to determine probable N- or O-linkage of glycan in spectrum
    and returns the ratio and peak intensities of oxonium ions

    Keyword Arguments
    spectrum -- array of spectral peaks to search
    """
    #handle tolerances better

    # Calculated Ratio:
    #(138.055 + 168.066)/(126.055 + 144.065)

    # round peaks in spectrum to 1 decimal place
    round_spec = [round(float(x[0]), 1) for x in spectrum]
    # the 4 oxonium ion peaks to calculate ratios between
    a = 138.055
    b = 168.066
    c = 126.055
    d = 144.065
    peaks_list = [a,b,c,d]
    peak_ints = [0,0,0,0]
    # check intensities of oxonium ion peaks
    for i in range(len(peaks_list)):
        if round(peaks_list[i], 1) in round_spec:
            peak_ints[i] = float(spectrum[round_spec.index(round(peaks_list[i], 1))][1])
        else:
            peak_ints[i] = 0
    #if ion c and d are not present, return 0 for ratio and intensity
    if peak_ints[2] + peak_ints[3] == 0:
        return 0, 0
    #else, return the ratio (a + b)/(c + d)
    else:
        return (peak_ints[0] + peak_ints[1])/(peak_ints[2] + peak_ints[3]), peak_ints

def calc_yo(mgf_spectra, glyco_frame):
    """
    calculates Y0 ion mass of a spectrum in an mgf entry based on the peptide mass, charge, and the reported glycan mass
     of the scan, equation is (raw_pepmass * charge) - ((charge - 0) * H) and returns singly charged Y0 peak mass, # may remove to search for all charges

    Keyword Arguments
    mgf_spectra -- Spectra object to calculate Y0 of
    glyco_frame -- dataframe containing masses of known glycans in "Glycan mass" column and scans of corresponding 
                    spectrum in "Scan" column
    """
    # retrieve scan number
    scan = mgf_spectra.scans
    #if no "Glycan mass" column, return None
    if glyco_frame["Glycan mass"][glyco_frame["Scan"] == int(scan)].empty is True:
        return None
    # else return corresponding glycan mass from glyco_frame
    else:
        glyc_mass = float(glyco_frame["Glycan mass"][glyco_frame["Scan"] == int(scan)])
    #calculate precursor peptide mass by subtracting glycan mass from raw pepmass
    raw_pepmass = calc_pepmass(mgf_spectra)
    calc = raw_pepmass - glyc_mass
    return calc

def calc_pepmass(mgf_spectra):
    """
    calculates precursor peptide mass an mgf entry based on the peptide mass, charge, and the reported glycan mass
     of the scan, equation is (raw_pepmass * charge) - ((charge - 0) * H) and returns singly charged precursor peptide mass, # may remove to search for all charges

    Keyword Arguments
    mgf_spectra -- Spectra object to calculate Y0 of
    glyco_frame -- dataframe containing masses of known glycans in "Glycan mass" column and scans of corresponding 
                    spectrum in "Scan" column
    """
    raw_pepmass = float(mgf_spectra.pepmass)
    charge = float(mgf_spectra.charge)
    # mass of Hydrogen used to calculate mass
    H = 1.00784
    calc = (raw_pepmass * charge) - ((charge - 0) * H)
    return calc

def is_iso(spectrum, peak):
    """
    Determines if target peak in a spectrum is an isotopic peak based on repeating m/z shifts of the same mass
    and returns identified base isotopic peak if target peak is found to be an isotope

    Keyword Arguments
    spectrum -- array of spectral peaks to search
    peak -- index of peak in spectrum to check
    """
    peaks = [peak]
    # calculate difference between target peak and previous peak in spectrum
    diff = round(float(spectrum[peak][0]) - float(spectrum[peak - 1][0]), 1)
    # if difference is greater than one, peak is not isotopic, return original peak, else check for isotopes
    if diff > 1:
        return peak, None
    else:
        # check peak to left of original peak
        peak -= 1
        peaks.append(peak)
        # set 'check_peak' to the peak to the left of 'peak' 
        check_peak = peak - 1
        #continue traversing peaks to the left until the difference between 2 peaks > 1
        while check_peak >= 0:
            # calculate difference between current peak and left peak
            next_diff = round(float(spectrum[peak][0]) - float(spectrum[check_peak][0]), 1)
            # if calculated difference = previous difference, add peak to list and continue traversing left
            if next_diff == diff:
                peaks.append(check_peak)
                peak = check_peak
                check_peak = peak - 1
            # else, check if chain of equal differences are 3 or longer
            else:
                # if longer than 2, check difference between peaks to determine whether isotopes and calculate charge of peaks
                if len(peaks) >= 3:
                    #if differences round to 0, register as no isotope and return original peak
                    if diff == 0:
                        return peaks[0], None
                    # else, return leftmost peak and charge as 1//diff
                    else:
                        return peaks[-1], 1//diff
                # if 2 or shorter, can't determine if isotopes, return original peak 
                else:
                    return peaks[0], None
        # if check_peak reaches 0, return non-isotopic
        return peaks[0], None
        
def iso_strip(spectrum):
    """
    Strips isotopic peaks from spectrum

    Keyword Arguments
    spectrum -- array of spectral peaks to be stripped of isotopic peaks

    Returns
    Array of spectral peaks stripped of identified isotopic peaks"""
    stripped_spec = []
    stripped_charges = []
    # set peak index to last peak in spectrum
    peak = len(spectrum) - 1
    # traverse left through peaks in spectrum, excluding leftmost peak
    while peak >= 1:
        #check if current peak is isotopic, and return current peak or leftmost isotopic peak
        # charge is None if not isotopic, or estimated int if isotopic
        peak, charge = is_iso(spectrum, peak)
        #add peak to left of new spectrum
        stripped_spec = [spectrum[peak]] + stripped_spec
        # record charge if found
        stripped_charges = [charge] + stripped_charges
        peak -= 1
    # return spectrum stripped of isotopes and charges of each peak if found
    return stripped_spec, stripped_charges

def yo_from_glycs(spectrum, glycan_data, tol = 0.01):
    """
    Returns list of probable Y0 peaks in spectrum by subtracting known possible glycan masses from raw peptide mass of the spectrum

    Keyword Arguments
    spectrum -- Spectra object to search
    glycan_data -- List of glycan masses
    tol -- m/z error tolerance for mass matching

    Optimisations:
        - Check charges
        - Work out whats happening with misses
    """
    ##TODO## Account for charged diffs
        #Charged prec and glycan? only prec/only glycan?

    #calculate precursor mass
    H = 1.00784
    pepmass = calc_pepmass(spectrum)
    ###subtract glycan massses and check for corresponding peak
    ###sort and stop if calc < 0  
    # retrieve list of known glycan masses (may retrieve from dataframe in future verision)
    glycan_masses = glycan_data
    candidate_peaks = []
    #cull low intensity and isotopic peaks
    spec = intensity_strip(spectrum.spectra, 0.01)
    spec = iso_strip(spec)[0]
    #
    for glyc_mass in glycan_masses:
        #calculate uncharged Y0 value
        calc_mass = pepmass - glyc_mass
        #calculate charged mass values for charges 1-5 
        for charge in range(1,5):
            # add charge * H to mass and divide by charge for charged mass
            charge_mass = (calc_mass + ((charge - 0) * H))/charge 
            #check spectrum for peaks with m/z = charged mass
            for peak in range(len(spec)):
                # if charged mass within tolerance value of peak, add to candidate list
                if float(spec[peak][0]) >= charge_mass - tol and float(spec[peak][0]) <= charge_mass + tol:
                    candidate_peaks.append((peak, charge, charge_mass, glyc_mass))
                else:
                    pass
    return candidate_peaks

def peak_pick(diffs, glycan_list):
    """
    From list of difference chains and list of glycan masses, pick most likely y0 peak

    Keyword Arguments
    diffs -- list of chains of differences between peaks where each difference is equal to the mass of oxonium ions (passed from glycan_traverse_and_validate function)
    glycan_list -- list of glycan masses (Currently unused)
    """
    #Priority order: Chain length, Validation, Leftmost peak, 204 presence
    # [(18, '1157.50317', 8)]
    #TODO:
    #validate against glycan
    ########UNFINISHED################
    #NEEDS TO BE AFTER CANDIDATE PICKING
    #validated_diffs = [peak for peak in diffs if peak[0] in glycan_list]
    #sort validated list by 

    #sort differences by m/z value in ascending order
    sorted_diffs = sorted(diffs, key = lambda x: x[2], reverse=True)
    longest_chain = 0
    #set initial pick to longest chain
    pick = sorted_diffs[0]
    # loop through difference chains
    for i in sorted_diffs:
        #set difference chain
        chain = i[0]
        #check if chain length is shorter than current longest chain
        if len(chain) - 1 < longest_chain:
            #break loop if chain is shorter than current longest chain 
            break
        else:
            #else continue
            pass
        #set position to end of chain
        chain_pos = len(chain) - 1
        #set end of chain to final diff in chain list
        #need this diff to be 203 representing HexNAc end of chain
        chain_end = round(chain[chain_pos][-1], 0)
        #loop back through diff jumps until 203 (HexNAc) jump is found
        while chain_end != 203:
            #move to previous position
            chain_pos -= 1
            #check if at start of list
            if chain_pos == -1:
                chain_pos = 0
                chain_end = 203
            else:
                #set end of chain to new position
                chain_end = round(chain[chain_pos][-1], 0)
        #Check if position of 203 diff is = or > than current longest chain
        #if =, check which peak is leftmost and use that
        if chain_pos == longest_chain:
            #retrieve current longest chain peak m/z
            current_mark = pick[0][longest_chain][0]
            #retrieve new candidate chain peak m/z
            candidate_mark = chain[chain_pos][0]
            #if new candidate m/z < currently selected peak, use new candidate
            if current_mark >= candidate_mark:
                pick = i
                longest_chain = chain_pos
            else:
                pass
        #if longer, replace candidate chain and longest chain length
        elif chain_pos > longest_chain:
            pick = i
            # print(pick)
            longest_chain = chain_pos
        #else pass
        else:
            pass
    # return whole picked chain and candidate Y0 peak index based on picked chain
    return pick, pick[0][longest_chain][1]

def glycan_traverse_and_validate(spec, glycan_data, tol = 0.01):
    """
    Returns probable Y0 candidate ---------

    Keyword Arguments
    spec -- Spectra object to search
    glycan_data -- List of glycan masses
    tol -- m/z error tolerance for mass matching

    Optimisations
        - Prioritise picking: Length O, HexNAc base X , validation list O, leftmost X
        - Search 1-4+ charge for each diff individually rather than 1-4+ charge for all diffs in traversal
    """
    #find list of candidates using yo_from_glycs
    candidate_peaks = yo_from_glycs(spec, glycan_data, tol)
    glycands = [x[0] for x in candidate_peaks]
    spectrum = iso_strip(intensity_strip(spec.spectra, tol))[0]
    #HexNAc, Hex, Fuc, HexNAcHex
    #find peptide peak by identifying mass shifts corresponding to glycan masses
    #list of oxonium ion masses and corresponding list of identities
    glycan_masses = [203.0794, 146.0579, 162.0528, 365.1322]
    glycan_ids = ["HexNAc", "Fucose", "Hex", "HexNAcHex"]
    #round glycan masses to 1 decimal
    round_glycans = [round(float(x),1) for x in glycan_masses]
    #set max possible mass for glycan to bound stop diff calcs when diffs exceed max
    max_mass = max(round_glycans)
    # set variable to track if a diff chain has been found
    any_found = False
    all_diffs = []
    all_cands = []
    #search until no mass shifts corresponding to glycan masses are found with a lower peak
    #start counter iterates backwards through list of peaks
    #currently searching last 10 peaks
    for start in range(0,11):
        #multiply by charge to account for charged peaks
        #calculating 1-5
        for charge in range(1,5):
            some_diffs = []
            #retrieve initial peak mass
            yo = len(spectrum) - 1 - start
            #add initial peak to begin checks from
            #peaks will be added to this list as they are located to be at target masses
            check_list = [yo]
            # counter for number of peaks checked for target mass shifts
            peak_check = 0
            # keep checking peaks until at the end of check list and no more chained peaks are found
            # check peaks right to left
            while peak_check < len(check_list):
                # set mass of first peak to check
                peak = check_list[peak_check]
                # iterate through peaks to the right of current target peak, looking for target mass shifts
                for next_peak in range(check_list[peak_check] - 1, -1, -1):
                    #calculate mass shift between target peak and current next_peak
                    diff = round((float(spectrum[peak][0]) - float(spectrum[next_peak][0])) * charge, 1)
                    # stop search if difference between the two peaks is > the maximum glycan mass in file
                    if diff > max_mass:
                        break
                    # else continue searching
                    else:
                        #check if diff is in glycan mass list
                        if diff in round_glycans:
                            #add pair of peaks, charge, and mass shift to charge_diffs
                            some_diffs.append((peak, next_peak, charge, diff))
                            #record candidate y0 peak change
                            #continue checking to continue chain of mass shifts
                            if next_peak < yo:
                                yo = next_peak
                            #record peak to start next loop from if is not already included
                            if next_peak in check_list:
                                pass
                            else:
                                check_list.append(next_peak)
                            #record glycan found in this spectrum
                            any_found = True
                        else:
                            pass
                peak_check += 1
            # record full chain of mass shifts, charge searched on, and the length of mass shift chain
            all_diffs.append([some_diffs, charge, len(some_diffs)])
            #if picked peak is next peak in spectra, likely no diff was found
            #record as normal
            #possibly implement special handling of case in future
            if yo == len(spectrum) - 1 - start:
                all_cands.append((yo, spectrum[yo][0], len(some_diffs)))
            #else peak was found
            #record peak found
            else:
                all_cands.append((yo, spectrum[yo][0], len(some_diffs)))
    # if no valid mass shifts are found, return empty list
    # possibly update to clearer flag later
    if any_found is False:
        return []
    # if a valid mass shift has been found
    else:
        #pick peak
        picked_peak = peak_pick(all_diffs, glycands)
        return(spectrum[picked_peak[1]], picked_peak[0])
