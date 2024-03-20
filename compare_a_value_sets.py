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
    def __init__(self,eupper_wavenumber,elower_wavenumber,wavelength_nm,avalue=[],loggf = [],j_upper=[],j_lower=[],csf_upper=[],csf_lower=[]):
        self.eupper     = eupper_wavenumber
        self.elower     = elower_wavenumber
        self.avalue     = avalue
        self.loggf      = loggf
        self.wavelength = wavelength_nm
        self.j_upper    = j_upper
        self.j_lower    = j_lower
        self.csf_upper  = csf_upper
        self.csf_lower  = csf_lower

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
    
    def reduce(self,e_upper_limit):
        indices_wanted = np.argwhere(self.eupper < e_upper_limit)
        #print(indices_wanted)

        if len(indices_wanted) > 0:
            initial = len(self.eupper)
            indices_wanted = np.concatenate(indices_wanted)
            self.eupper     = self.eupper[indices_wanted] 
            self.elower     = self.elower[indices_wanted] 
             
            
            self.wavelength = self.wavelength[indices_wanted] 
            if len(self.j_upper) > 0:
                self.j_upper    = self.j_upper[indices_wanted] 
            if len(self.j_lower) > 0:
                self.j_lower    = self.j_lower[indices_wanted]
            if len(self.avalue) > 0:
                self.loggf    = self.loggf[indices_wanted]
            if len(self.loggf) > 0:
                self.avalue    = self.avalue[indices_wanted]
            if len(self.csf_upper) > 0:
                self.csf_upper    = self.csf_upper[indices_wanted] 
            if len(self.j_lower) > 0:
                self.csf_lower    = self.csf_lower[indices_wanted]
            #self.avalue     = self.avalue[indices_wanted]
            #self.loggf      = self.loggf[indices_wanted] 
            
            final = len(self.eupper)
            print("reduced data set from ",initial," to ",final," lines.")

        else:
            print("no upper levels below this threshold, dataset not changing.")
 

        


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

        base_ju = self.base.j_upper
        base_jl = self.base.j_lower
        comp_ju = self.compared.j_upper        
        comp_jl = self.compared.j_lower

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
                    #print(len(base_ju))


                    if np.abs(test_upper - current_upper) < tolerance_wavenumber:
                        
                        #this is so messy. i am sorry for 
                        #my future self who has to debug this
                        #print(len(base_ju))

                        if (len(base_ju) > 0) and (len(base_jl) > 0) and (len(comp_ju) > 0) and (len(comp_jl) > 0) :

                                if (base_ju[jj] == comp_ju[ii]) and (base_jl[jj] == comp_jl[ii]):
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
                        else:
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
    
    def sort(self,sorting = ''):

        base_data_wl    = self.base.wavelength
        compared_wl     = self.compared_wl_raw
        base_avalue     = self.base.avalue
        compared_avalue = self.compared_a_values_raw
        base_el         = self.base.elower
        base_eu         = self.base.eupper
        comp_el         = self.compared_elower_raw
        comp_eu         = self.compared_eupper_raw

        csf_upper = self.base.csf_upper
        csf_lower = self.base.csf_lower
        j_upper = self.base.j_upper
        j_lower = self.base.j_lower 
        avalue_ratio = base_avalue/compared_avalue

        temp_array_for_my_sanity = base_el 
        temp_array_for_my_sanity = np.concatenate((temp_array_for_my_sanity,base_eu)) 

        levels = (np.unique(temp_array_for_my_sanity))
        
        print(levels)

        if sorting == 'avalue_ratio':
            sorted_indices = np.argsort(avalue_ratio)
        elif sorting == 'wavelength':
            sorted_indices = np.argsort(base_data_wl)
        elif sorting == 'base_el':

            sorted_indices = np.argsort(base_el)
            base_data_wl    = base_data_wl[sorted_indices]  
            compared_wl     = compared_wl[sorted_indices]  
            base_avalue     = base_avalue[sorted_indices]  
            compared_avalue = compared_avalue[sorted_indices]  
            base_el         = base_el[sorted_indices]  
            base_eu         = base_eu[sorted_indices]  
            comp_el         = comp_el[sorted_indices]  
            comp_eu         = comp_eu[sorted_indices]   
            csf_lower_copy = csf_lower
            csf_upper_copy = csf_upper
            j_lower_copy = j_lower
            j_upper_copy = j_upper

            csf_lower = []  
            csf_upper = [] 
            j_lower = []
            j_upper = []
            for jj in sorted_indices:
                csf_lower.append(csf_lower_copy[jj])
                csf_upper.append(csf_upper_copy[jj])
                j_lower.append(j_lower_copy[jj])
                j_upper.append(j_upper_copy[jj])
                #sorted_indices = np.array([])

            unique_energies,counts = np.unique(base_el,return_counts=True)
            #print(unique_energies,counts,np.shape(counts))
            offset = 0
            print(sorted_indices)
            for ii in range(0,np.shape(counts)[0]-1):
        
                #print(offset,offset + counts[ii])
                
                new_indices = np.argsort(base_eu[offset:offset + counts[ii]]) + offset

                #print(base_eu[new_indices])

                sorted_indices[offset:offset + counts[ii]] = new_indices
                offset += counts[ii]

                #print(base_el_temp[counts[ii-1]-1:counts[ii]])
            print(sorted_indices)




        else:
            sorted_indices = np.arange(0,len(base_data_wl),1)
        
        base_data_wl    = base_data_wl[sorted_indices]  
        compared_wl     = compared_wl[sorted_indices]  
        base_avalue     = base_avalue[sorted_indices]  
        compared_avalue = compared_avalue[sorted_indices]  
        base_el         = base_el[sorted_indices]  
        base_eu         = base_eu[sorted_indices]  
        comp_el         = comp_el[sorted_indices]  
        comp_eu         = comp_eu[sorted_indices]  

        #print(sorted_indices)

        #j_lower         =  j_lower  [sorted_indices]  
        #j_upper         =  j_upper  [sorted_indices]  

        csf_lower_copy = csf_lower
        csf_upper_copy = csf_upper
        j_lower_copy = j_lower
        j_upper_copy = j_upper

        csf_lower = []  
        csf_upper = [] 
        j_lower = []
        j_upper = []
        for jj in sorted_indices:
            csf_lower.append(csf_lower_copy[jj])
            csf_upper.append(csf_upper_copy[jj])
            j_lower.append(j_lower_copy[jj])
            j_upper.append(j_upper_copy[jj])
        
        self.base.wavelength         = base_data_wl      
        self.compared_wl_raw         = compared_wl        
        self.base.avalue             = base_avalue        
        self.compared_a_values_raw   = compared_avalue    
        self.base.elower             = base_el               
        self.base.eupper             = base_eu              
        self.compared_elower_raw     = comp_el                  
        self.compared_eupper_raw     = comp_eu    
        self.base.j_lower = j_lower
        self.base.j_upper = j_upper
        self.base.csf_lower = csf_lower
        self.base.csf_upper = csf_upper            

        return 0


    def write_out_data_file(self,filename,label1,label2,ignore_null=False,sorting=''):
        file = open(filename,'w')
        format_string = 'F10.2,1X,F10.2,2ES10.2,F11.2,1X,F11.2,F11.2,1X,F10.1,5X,a11,1X,,A10,1X,F3.1,4X,A11,1x,A10,1X,F3.1'
        CON = label1 
        KUR = label2

        header = "#    "+CON+"wl      "+KUR+"wl      "+CON+"a      "+KUR+"a     "+CON+"el      "+CON+"eu     "+KUR+"el      "+KUR+"eu" +'\n'

        file.write(header)
        line = ff.FortranRecordWriter(format_string)

        base_data_wl    = self.base.wavelength
        compared_wl     = self.compared_wl_raw
        base_avalue     = self.base.avalue
        compared_avalue = self.compared_a_values_raw
        base_el         = self.base.elower
        base_eu         = self.base.eupper
        comp_el         = self.compared_elower_raw
        comp_eu         = self.compared_eupper_raw
        csf_upper = self.base.csf_upper
        csf_lower = self.base.csf_lower
        j_upper = self.base.j_upper
        j_lower = self.base.j_lower 
        avalue_ratio = base_avalue/compared_avalue

        temp_array_for_my_sanity = base_el 
        temp_array_for_my_sanity = np.concatenate((temp_array_for_my_sanity,base_eu)) 

        levels = (np.unique(temp_array_for_my_sanity))
        
        print(levels)

        if sorting == 'avalue_ratio':
            sorted_indices = np.argsort(avalue_ratio)
        elif sorting == 'wavelength':
            sorted_indices = np.argsort(base_data_wl)
        elif sorting == 'base_el':

            sorted_indices = np.argsort(base_el)
            base_data_wl    = base_data_wl[sorted_indices]  
            compared_wl     = compared_wl[sorted_indices]  
            base_avalue     = base_avalue[sorted_indices]  
            compared_avalue = compared_avalue[sorted_indices]  
            base_el         = base_el[sorted_indices]  
            base_eu         = base_eu[sorted_indices]  
            comp_el         = comp_el[sorted_indices]  
            comp_eu         = comp_eu[sorted_indices]   
            csf_lower_copy = csf_lower
            csf_upper_copy = csf_upper
            j_lower_copy = j_lower
            j_upper_copy = j_upper

            csf_lower = []  
            csf_upper = [] 
            j_lower = []
            j_upper = []
            for jj in sorted_indices:
                csf_lower.append(csf_lower_copy[jj])
                csf_upper.append(csf_upper_copy[jj])
                j_lower.append(j_lower_copy[jj])
                j_upper.append(j_upper_copy[jj])
                #sorted_indices = np.array([])

            unique_energies,counts = np.unique(base_el,return_counts=True)
            #print(unique_energies,counts,np.shape(counts))
            offset = 0
            print(sorted_indices)
            for ii in range(0,np.shape(counts)[0]-1):
        
                #print(offset,offset + counts[ii])
                
                new_indices = np.argsort(base_eu[offset:offset + counts[ii]]) + offset

                #print(base_eu[new_indices])

                sorted_indices[offset:offset + counts[ii]] = new_indices
                offset += counts[ii]

                #print(base_el_temp[counts[ii-1]-1:counts[ii]])
            print(sorted_indices)




        else:
            sorted_indices = np.arange(0,len(base_data_wl),1)


        base_data_wl    = base_data_wl[sorted_indices]  
        compared_wl     = compared_wl[sorted_indices]  
        base_avalue     = base_avalue[sorted_indices]  
        compared_avalue = compared_avalue[sorted_indices]  
        base_el         = base_el[sorted_indices]  
        base_eu         = base_eu[sorted_indices]  
        comp_el         = comp_el[sorted_indices]  
        comp_eu         = comp_eu[sorted_indices]  

        #print(sorted_indices)

        #j_lower         =  j_lower  [sorted_indices]  
        #j_upper         =  j_upper  [sorted_indices]  

        csf_lower_copy = csf_lower
        csf_upper_copy = csf_upper
        j_lower_copy = j_lower
        j_upper_copy = j_upper

        csf_lower = []  
        csf_upper = [] 
        j_lower = []
        j_upper = []

        if len(csf_lower)>0:

            for jj in sorted_indices:
                csf_lower.append(csf_lower_copy[jj])
                csf_upper.append(csf_upper_copy[jj])
                j_lower.append(j_lower_copy[jj])
                j_upper.append(j_upper_copy[jj])



        for jj in range(0,len(base_data_wl)):
            lower_index = np.argwhere(base_el[jj] == levels)
            upper_index = np.argwhere(base_eu[jj] == levels)

            if len(lower_index) > 0: 
                lower_index = lower_index[0][0] + 1
            else:
                lower_index = 'indfail'

            if len(upper_index) > 0: 
                upper_index = upper_index[0][0] + 1
            else:
                lower_index = 'indfail'

            array = [base_data_wl[jj],compared_wl[jj],base_avalue[jj],compared_avalue[jj],base_el[jj],base_eu[jj],comp_el[jj],comp_eu[jj]]
            #string_to_be_written = str(wl_connor[jj]) + " " +str(kurucz_a_values_raw[jj]) +'\n'
            if len(csf_lower)>0:

                lower_array = ['level'+str(lower_index) ,csf_lower[jj],j_lower[jj]]
                upper_array = ['level'+str(upper_index) ,csf_upper[jj],j_upper[jj]]
                array.extend(lower_array)
    
                array.extend(upper_array)
            
            if not(ignore_null and (compared_wl[jj] == -1.0)):
                file.write(line.write(array))
                file.write("\n")
            

        file.close()

#def write_out_many_datas(heading_string_array,list_of_comparison_classes):
#    #assumes all the comparison classes have  the same base data. messy for now.
#    #further assumes they have all been sorted the same way 
#
#    num_comparisons = len(heading_string_array)
#    assert(num_comparisons == len(list_of_comparison_classes))
#
#    format_string = 'F10.2,'+ str(num_comparisons)+'ES10.2,F11.2,1X,F11.2,5X,a11,1X,,A10,1X,F3.1,4X,A11,1x,A10,1X,F3.1'
#    base_data_wl    = list_of_comparison_classes[0].base.wavelength
#
#    base_data_wl    = list_of_comparison_classes[0].base.wavelength
#    compared_wl     = list_of_comparison_classes[0].compared_wl_raw
#    base_avalue     = list_of_comparison_classes[0].base.avalue
#    compared_avalue = list_of_comparison_classes[0].compared_a_values_raw
#    base_el         = list_of_comparison_classes[0].base.elower
#    base_eu         = list_of_comparison_classes[0].base.eupper
#    comp_el         = list_of_comparison_classes[0].compared_elower_raw
#    comp_eu         = list_of_comparison_classes[0].compared_eupper_raw
#    csf_upper = list_of_comparison_classes[0].base.csf_upper
#    csf_lower = list_of_comparison_classes[0].base.csf_lower
#    j_upper = list_of_comparison_classes[0].base.j_upper
#    j_lower = list_of_comparison_classes[0].base.j_lower 
#    
#    temp_array_for_my_sanity = base_el 
#    temp_array_for_my_sanity = np.concatenate((temp_array_for_my_sanity,base_eu)) 
#    levels = (np.unique(temp_array_for_my_sanity))
#
#    temp_array_for_my_sanity = base_el 
#    temp_array_for_my_sanity = np.concatenate((temp_array_for_my_sanity,base_eu)) 
#
#    levels = (np.unique(temp_array_for_my_sanity))
#
#    for jj in range(0,len(base_data_wl)):   
#            lower_index = np.argwhere(base_el[jj] == levels)
#            upper_index = np.argwhere(base_eu[jj] == levels)
#
#            if len(lower_index) > 0: 
#                lower_index = lower_index[0][0] + 1
#            else:
#                lower_index = 'indfail'
#
#            if len(upper_index) > 0: 
#                upper_index = upper_index[0][0] + 1
#            else:
#                lower_index = 'indfail'
#
#            array = [base_data_wl[jj],compared_wl[jj],base_avalue[jj],compared_avalue[jj],base_el[jj],base_eu[jj],comp_el[jj],comp_eu[jj]]
#            #string_to_be_written = str(wl_connor[jj]) + " " +str(kurucz_a_values_raw[jj]) +'\n'
#            lower_array = ['level'+str(lower_index) ,csf_lower[jj],j_lower[jj]]
#            upper_array = ['level'+str(upper_index) ,csf_upper[jj],j_upper[jj]]
#            array.extend(lower_array)
#
#            array.extend(upper_array)
#            
#            if not(ignore_null and (compared_wl[jj] == -1.0)):
#                file.write(line.write(array))
#                file.write("\n")
#            
#
#        file.close()
#    return 0 






def read_kurucz(kurucz_path,wavelength_convert=True,calculate_a_values=True):
    format_string = 'F11.4,F7.3,F6.2,F12.3,F5.1,1X,A10,F12.3,F5.1,1X,A10,F6.2,F6.2,F6.2,A4,I2,I2,I3,F6.3,I3,F6.3,I5,I5,A10,I5,I5'
    print("reading in Kurucz. Note that Kurucz data has wls in air (nm) and transitions strengths in log gf.")
    f = open(kurucz_path,'r')
    f_opened = f.readlines()
    f.close()
    numlines = len(f_opened)
    if wavelength_convert == False:
        print("WARNING: wavelength conversion is off - assumes your wavelength data is in ritz vac.")
    wavelengths_air_nm = np.zeros(numlines)
    loggf_kurucz = np.zeros(numlines)
    upper_j = np.zeros(numlines)
    lower_j = np.zeros(numlines)
    e_upper = np.zeros(numlines)
    e_lower = np.zeros(numlines)
    csf_label_lower = []
    csf_label_upper = []
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
            csf_label_lower.append(array[5])
            csf_label_upper.append(array[8])


        else:
            upper_j[ii] = array[4]
            lower_j[ii] = array[7]
            e_upper[ii] = el
            e_lower[ii] = eu
            csf_label_lower.append(array[8])
            csf_label_upper.append(array[5])
    dataclass = dataset(e_upper,e_lower,wavelength_nm=wavelengths_air_nm,loggf=loggf_kurucz,j_upper=upper_j,j_lower=lower_j,csf_lower=csf_label_lower,csf_upper=csf_label_upper)

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
