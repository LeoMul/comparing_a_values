import numpy as np
import fortranformat as ff
class dataset:
    def __init__(self,eupper,elower,avalue,wavelength):
        self.eupper = eupper
        self.elower = elower
        self.avalue = avalue
        self.wavelength = wavelength

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
        self.leastsquares_a = np.sum(np.power(comparison.avalue_comp_found-comparison.avalue_base_found,2))


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
                for ii in possible_wavelengths_indices:
                    test_upper = compared_uppers[ii]
                    if np.abs(test_upper -current_upper) < tolerance_wavenumber:
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
                    print("wavelength found ",current_wavelength, "but no matching upper found.")
                    print("requested upper ",current_upper,"with possible comparable uppers ",possible_uppers )
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
        

    def write_out_data_file(self,filename,label1,label2,ignore_null):
        file = open(filename,'w')
        format_string = 'F10.2,1X,F10.2,2ES10.2,F10.1,1X,F10.1,F10.1,1X,F10.1'
        CON = label1 
        KUR = label2

        header = "#    "+CON+"wl      "+KUR+"wl      "+CON+"a      "+KUR+"a     "+CON+"el      "+CON+"eu     "+KUR+"el      "+KUR+"eu" +'\n'

        file.write(header)
        base_data_wl = self.base.wavelength
        line = ff.FortranRecordWriter(format_string)

        for jj in range(0,len(base_data_wl)):
            array = [base_data_wl[jj],self.compared_wl_raw[jj],self.base.avalue[jj],self.compared_a_values_raw[jj],self.base.elower[jj],self.base.eupper[jj],self.compared_elower_raw[jj],self.compared_eupper_raw[jj]]
            #string_to_be_written = str(wl_connor[jj]) + " " +str(kurucz_a_values_raw[jj]) +'\n'
            if not(ignore_null and (self.compared_wl_raw[jj] == -1.0)):
                file.write(line.write(array))
                file.write("\n")

        file.close()
