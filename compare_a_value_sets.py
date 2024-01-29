import numpy as np
import fortranformat as ff
import scipy.stats as ss
from PyAstronomy import pyasl

OSCILLATOR_CONVERSION_CONST = 6.67025177e13 #This number is from Drake's handbook on atomic physics.

def wl_conversion(wavelengths_in_nm):
    #some data sets have wl > 200nm in air, for some reason.
    #angstrom and vacuum malarky.
    num = len(wavelengths_in_nm)

    #briefly convert to angs because of weird astronomer things
    wavelengths_in_nm = 10.0 * wavelengths_in_nm
    counter = 0
    for ii in range(0,num):
        if wavelengths_in_nm[ii] >= 2000.0:
            counter += 1
            wavelengths_in_nm[ii] = pyasl.airtovac2(wavelengths_in_nm[ii], mode="edlen53", precision=1e-9, maxiter=30)
    wavelengths_in_nm = 0.10 * wavelengths_in_nm
    print("converted ",counter," wavelengths to vac")
    return wavelengths_in_nm

def convert_loggf_to_avalue(wavelength_array_nm,loggf_array,upper_weight_array):
    num = len(wavelength_array_nm)
    avalues = np.zeros(num)
    #from drakes book again. might need to be lower weights here...
    avalues = np.power(10.0,loggf_array) / upper_weight_array
    avalues = avalues * OSCILLATOR_CONVERSION_CONST / np.power(wavelength_array_nm,2)

    return avalues

class dataset:
    def __init__(self,eupper_wavenumber,elower_wavenumber,wavelength_nm,avalue=[],loggf = [],j_upper=[],j_lower=[]):
        self.eupper = eupper_wavenumber
        self.elower = elower_wavenumber
        self.avalue = avalue
        self.loggf = loggf
        self.wavelength = wavelength_nm
        self.j_upper = j_upper
        self.j_lower = j_lower

    #converts your wavelengths if they are in air.
    def convert_wavelengths_to_vaccuum(self):   
        self.wavelength = wl_conversion(self.wavelength)
    
    def calculate_a_values(self):
        loggf = self.loggf

        assert (len(loggf) == len(self.j_upper)),"length of jupper and loggf not the same"

        if len(loggf) == 0:
            print("no loggfs found to convert to a values. have you given me the loggfs?")
            return 
        elif len(self.j_upper) == 0:
            print("no upper j's found (lower j's not needed here)")
        else:
            weights = 2.0 * self.j_upper + 1.0
            print("converting loggf to avalue, assuming wl is in nm and in vac.")
            self.avalue = convert_loggf_to_avalue(self.wavelength,self.loggf,weights)

    def calculate_loggf(self):
        assert False,"not yet implemented"
        print("not yet implemented")
        return 

class compared_data:
    def __init__(self,comparison):
        self.base_eu = comparison.eupper_base_found
        self.base_el = comparison.elower_base_found
        self.base_wl = comparison.foundwlbase
        self.comp_eu = comparison.eupper_comp_found
        self.comp_el = comparison.elower_comp_found
        self.comp_wl = comparison.foundwlcomp
        self.base_a  = comparison.avalue_base_found
        self.comp_a  = comparison.avalue_comp_found
        self.leastsquares_a = np.sum(np.power((comparison.avalue_comp_found-comparison.avalue_base_found),2))
        self.spearman = ss.stats.spearmanr(comparison.avalue_base_found,comparison.avalue_comp_found)

class comparison:
    def __init__(self,base_set,compared_dataset):
        number_to_be_compared = len(base_set.wavelength)
        #print(base_set.wavelength)
        #print(number_to_be_compared)
        self.base = base_set
        self.compared = compared_dataset
        self.number_to_be_compared = number_to_be_compared
        self.wavelengths_not_found = []
        self.wavelength_indices_found = []
        self.compared_positions_of_base_wls_in_compared_dataset = []
        self.eupper_base_found = []
        self.elower_base_found = []
        self.eupper_comp_found = []
        self.elower_comp_found = []
        self.avalue_comp_found = []
        self.avalue_base_found = []
        self.compared_a_values_raw = np.zeros(number_to_be_compared)
        self.compared_wl_raw = -1.0*np.ones(number_to_be_compared)
        self.compared_eupper_raw = -1.0*np.ones(number_to_be_compared)
        self.compared_elower_raw = -1.0*np.ones(number_to_be_compared)
        self.foundwlbase = []
        self.foundwlcomp = []
        self.compared_data = []

    def calculate_comparison(self,tolerance_wl,tolerance_wavenumber):

        #tolerance_wl = 0.01
        compared_wavelengths = self.compared.wavelength
        compared_uppers = self.compared.eupper
        compared_lowers = self.compared.elower
        compared_avalues = self.compared.avalue

        base_wavelengths = self.base.wavelength
        base_uppers = self.base.eupper
        wavelength_indices_found = []
        wavelengths_not_found = []

        base_wavelengths_found = []
        comp_wavelengths_found = []
        print("there are",self.number_to_be_compared,'wavelengths to be compared')
        for jj in range(0,self.number_to_be_compared):
            current_wavelength = base_wavelengths[jj]
            current_upper = base_uppers[jj]
            possible_wavelengths_indices = np.argwhere(np.abs(compared_wavelengths-current_wavelength)<tolerance_wl)
            #print(possible_wavelengths_indices)
            match_found = False

            if len(possible_wavelengths_indices) > 0:
                possible_wavelengths_indices = np.concatenate(possible_wavelengths_indices,axis=0)

                possible_uppers = compared_uppers[possible_wavelengths_indices]

                index_order = np.argsort(np.abs(current_upper - possible_uppers))
                possible_wavelengths_indices = possible_wavelengths_indices[index_order]
                possible_uppers = possible_uppers[index_order]

                for ii in possible_wavelengths_indices:
                    test_upper = compared_uppers[ii]
                    
                    if np.abs(test_upper - current_upper) < tolerance_wavenumber:
                        match_found = True 
                        self.compared_positions_of_base_wls_in_compared_dataset.append(ii)
                        wavelength_indices_found.append(jj)
                        comp_wavelengths_found.append(compared_wavelengths[ii])
                        base_wavelengths_found.append(current_wavelength)
                        self.compared_wl_raw[jj] = compared_wavelengths[ii]
                        self.compared_a_values_raw[jj] = compared_avalues[ii]
                        self.compared_eupper_raw[jj] = compared_uppers[ii]
                        self.compared_elower_raw[jj] = compared_lowers[ii]

                        break 
                if match_found == False:
                    print("requested upper ",current_upper," with possible comparable uppers ",possible_uppers )
                    print("wavelength found ",current_wavelength, "but no matching upper found.")
            else:
                print("wavelength ",current_wavelength,' not found in compared data set')
                wavelengths_not_found.append(current_wavelength)

        #print(positions_of_connor_wls_in_compared_data_set)
        wavelength_indices_found = np.array(wavelength_indices_found)
        wavelengths_not_found = np.array(wavelengths_not_found)
        self.wavelength_indices_found = wavelength_indices_found
        self.wavelengths_not_found = wavelengths_not_found
        print("total wavelengths not found: ",len(wavelengths_not_found))
        print("wavelengths found for comparisons:",len(wavelength_indices_found))
        self.eupper_base_found = self.base.eupper[wavelength_indices_found]
        self.elower_base_found = self.base.elower[wavelength_indices_found]
        self.eupper_comp_found = self.compared.eupper[self.compared_positions_of_base_wls_in_compared_dataset]
        self.elower_comp_found = self.compared.elower[self.compared_positions_of_base_wls_in_compared_dataset]
        self.avalue_comp_found = self.compared.avalue[self.compared_positions_of_base_wls_in_compared_dataset]
        self.avalue_base_found = self.base.avalue[wavelength_indices_found]
        self.foundwlbase = base_wavelengths_found
        self.foundwlcomp = comp_wavelengths_found
        self.compared_data = compared_data(self)
        

    def write_out_data_file(self,filename,label1,label2,ignore_null=False,sort_by_a_value_ratio=False):
        file = open(filename,'w')
        format_string = 'F10.2,1X,F10.2,2ES10.2,F10.1,1X,F10.1,F10.1,1X,F10.1'
        CON = label1 
        KUR = label2

        header = "#    "+CON+"wl      "+KUR+"wl      "+CON+"a      "+KUR+"a     "+CON+"el      "+CON+"eu     "+KUR+"el      "+KUR+"eu" +'\n'

        file.write(header)
        base_data_wl    = self.base.wavelength
        compared_wl     = self.compared_wl_raw
        base_avalue     = self.base.avalue
        compared_avalue = self.compared_a_values_raw
        base_el         = self.base.elower
        base_eu         = self.base.eupper
        comp_el         = self.compared_elower_raw
        comp_eu         = self.compared_eupper_raw

        line = ff.FortranRecordWriter(format_string)



        if sort_by_a_value_ratio:
            avalue_ratio = base_avalue/compared_avalue
            sorted_indices = np.argsort(avalue_ratio)
            base_data_wl    = base_data_wl[sorted_indices]  
            compared_wl     = compared_wl[sorted_indices]  
            base_avalue     = base_avalue[sorted_indices]  
            compared_avalue = compared_avalue[sorted_indices]  
            base_el         = base_el[sorted_indices]  
            base_eu         = base_eu[sorted_indices]  
            comp_el         = comp_el[sorted_indices]  
            comp_eu         = comp_eu[sorted_indices]  



        for jj in range(0,len(base_data_wl)):
            array = [base_data_wl[jj],compared_wl[jj],base_avalue[jj],compared_avalue[jj],base_el[jj],base_eu[jj],comp_el[jj],comp_eu[jj]]
            #string_to_be_written = str(wl_connor[jj]) + " " +str(kurucz_a_values_raw[jj]) +'\n'
            if not(ignore_null and (compared_wl[jj] == -1.0)):
                file.write(line.write(array))
                file.write("\n")

        file.close()

def read_kurucz(kurucz_path,wavelength_convert=True,calculate_a_values=True):
    format_string = 'F11.4,F7.3,F6.2,F12.3,F5.1,1X,A10,F12.3,F5.1,1X,A10,F6.2,F6.2,F6.2,A4,I2,I2,I3,F6.3,I3,F6.3,I5,I5,A10,I5,I5'
    print("reading in Kurucz. Note that Kurucz data has wls in air (nm) and transitions strengths in log gf.")
    f = open(kurucz_path,'r')
    f_opened = f.readlines()
    f.close()
    numlines = len(f_opened)
    
    wavelengths_air_nm = np.zeros(numlines)
    loggf_kurucz = np.zeros(numlines)
    upper_j = np.zeros(numlines)
    lower_j = np.zeros(numlines)
    e_upper = np.zeros(numlines)
    e_lower = np.zeros(numlines)
    reader = ff.FortranRecordReader(format_string)
    for ii in range(0,numlines):
        array = reader.read(f_opened[ii])
        wavelengths_air_nm[ii] = array[0]
        loggf_kurucz[ii] = array[1]

        el = np.abs(array[3])
        eu = np.abs(array[6])

        if eu > el:
            lower_j[ii] = array[4]
            upper_j[ii] = array[7]
            e_upper[ii] = eu
            e_lower[ii] = el
        else:
            upper_j[ii] = array[4]
            lower_j[ii] = array[7]
            e_upper[ii] = el
            e_lower[ii] = eu
    dataclass = dataset(e_upper,e_lower,wavelength_nm=wavelengths_air_nm,loggf=loggf_kurucz,j_upper=upper_j,j_lower=lower_j)

    if wavelength_convert:
        dataclass.convert_wavelengths_to_vaccuum()
    
    if calculate_a_values:
        #if not wavelength_convert:
        #    print("You have selected to not convert the wavelengths to vac.")
        #    print("This is required to calcualte a avalues, so I am overriding your choice, belligerent user.")
        #    dataclass.convert_wavelengths_to_vaccuum()
        dataclass.calculate_a_values()

    return dataclass
        #e_lower[ii] = np.min(np.abs(array[3],np.abs(array[6])))
        #e_upper[ii] = np.max(np.abs(array[3],np.abs(array[6])))
