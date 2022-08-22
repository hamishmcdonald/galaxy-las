import math
import csv
import os
import laspy
import numpy
import scipy
import astropy
import traceback
import timeline

#constants and their units
WEIN_CONSTANT = 2897771.9 #nM K
BOLTZMANN_CONSTANT = 1.380649e-23 #J/(K)
PLANCK_CONSTANT = 6.62607015e-34 #J s
LIGHT_SPEED_CONSTANT = 29979245800 #cm/(s)

def main():
    #local file location of GaiaSource files
    # the files can be found for download here:
    #http://cdn.gea.esac.esa.int/Gaia/gedr3/gaia_source/
    file_directory = 'D:\GaiaSourceTest'

    #interate through GaiaSource files and open them provided they are csv files
    for gaia_file in os.listdir(file_directory):
        if gaia_file.endswith('.csv'):
            try:
                with open(file_directory + '/' + gaia_file) as current_csv:
                    
                    #create las header then create las file using header
                    header = laspy.LasHeader(version="1.4", point_format=2)
                    galaxy_data = laspy.LasData(header)

                    #create extra dimensions in las file for storing meta data
                    galaxy_data.add_extra_dims([
                    laspy.ExtraBytesParams(name="solution_id", type="uint64"),
                    laspy.ExtraBytesParams(name="designation", type="uint64"),
                    laspy.ExtraBytesParams(name="source_id", type="uint64"),
                    ])

                    #create arrays for temporary storage of data to be added to the las file
                    x, y, z = [], [], []
                    red, green, blue = [], [], []
                    solution_id, designation, source_id= [], [], []
                   
                   #create csv reader with optional fieldname capabilites
                   # more information on the fieldnames used in the GaiaSource files can be found here:
                   #https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/
                   # sec_dm_main_source_catalogue/ssec_dm_gaia_source.html
                    current_csv_reader = csv.DictReader(current_csv)

                    row_number = 1

                    #interate through each star in the current GaiaSource file
                    for row in current_csv_reader:
                        try:
                            row_number += 1

                            #calculate x, y, z coordinates of the star using parallax, galactic longitude and latitude
                            # more information on the formulas can be found here:
                            #https://en.wikipedia.org/wiki/Galactic_coordinate_system
                            if row['parallax'] == '':
                                raise Exception("no parallax")
                            else:
                                x_value = math.cos(float(row['b'])) * math.cos(float(row['l'])) / float(row['parallax'])
                                y_value = math.cos(float(row['b'])) * math.sin(float(row['l'])) / float(row['parallax'])
                                z_value = math.sin(float(row['b'])) / float(row['parallax'])
                                
                            #calculate peak wavelength of light emitted from the star
                            if not row['nu_eff_used_in_astrometry'] == '':
                                peak_wavelength_value = 1 / float(row['nu_eff_used_in_astrometry']) * 1000
                            elif not row['pseudocolour'] == '':
                                peak_wavelength_value = 1 / float(row['pseudocolour']) * 1000
                            else:
                                raise Exception("no nu_eff_used_in_astronomy or pseudocolour")
                         
                            #calculate temperature of the star with Wien's law formula using displacement constant and peak 
                            # wavelength
                            #https://www.omnicalculator.com/physics/wiens-law
                            temperature = WEIN_CONSTANT / peak_wavelength_value

                            print(temperature)
                            
                            #calculate rgb values of star using temperature
                            #red_value, green_value, blue_value = timeline.rgb_from_T(temperature, False, 255, False)
                            #https://en.wikipedia.org/wiki/CIE_1931_color_space
                            if 675 < temperature < 250000:
                                red_value, green_value, blue_value = calculate_rgb(temperature)
                            else:
                                raise Exception("temperature outside of normal range: " + str(temperature))
                            

                            #store unique source indentifiers and designations
                            solution_id_value = int(row['solution_id'])
                            designation_value = int(row['designation'][11:])
                            source_id_value = int(row['source_id'])

                            #add values to associated temporary array, done separately incase calculations produce an 
                            # exception
                            x.append(x_value)
                            y.append(y_value)
                            z.append(z_value)
                            red.append(red_value)
                            green.append(green_value)
                            blue.append(blue_value)
                            solution_id.append(solution_id_value)
                            designation.append(designation_value)
                            source_id.append(source_id_value)
                            
                            #print to console if no exceptions occured for the star
                            print("data in row " + str(row_number) + " successfully added")

                        #print to console if an exception occured for the star
                        except Exception as row_exception:
                            print("Exception occured in row " + str(row_number) + " in file " + gaia_file + ": ", 
                                row_exception)
                            traceback.print_exc()
                            print("\n")

                #add temporary arrays to the las file
                galaxy_data.x = x
                galaxy_data.y = y
                galaxy_data.z = z
                galaxy_data.red = red
                galaxy_data.green = green
                galaxy_data.blue = blue
                galaxy_data.solution_id = solution_id
                galaxy_data.designation = designation
                galaxy_data.source_id = source_id

                #write las file to local storage
                galaxy_data.write(gaia_file + ".las")

                #print to console if no exceptions occured for the GaiaSource file 
                print(gaia_file + " was successfully converted")

            #print to console if an exception occured for the GaiaSource file
            except Exception as file_exception:
                print("Exception occured in file " + gaia_file + ": ", file_exception)
                traceback.print_exc()

def calculate_rgb(temperature):
    #wavelength in cm
    wavelength_cm = numpy.linspace(3.5e-5, 8e-5, 1e-5)
    #wavelength in nm
    lam = numpy.linspace(350, 800, 100)
    #spectral radiance of a blackbody in cgs units
    x = PLANCK_CONSTANT * LIGHT_SPEED_CONSTANT / wavelength_cm / BOLTZMANN_CONSTANT / (temperature)
    B = (2 * PLANCK_CONSTANT * LIGHT_SPEED_CONSTANT ** 2 / wavelength_cm ** 5 / (numpy.exp(x) - 1))

    #Color matching functions
    lamcie,xbar,ybar,zbar = cie()
    
    #Interpolate to same axis
    B = numpy.interp(lamcie, lam, B)                 

    #Tristimulus values
    X = scipy.integrate.simps(B * xbar, lamcie)
    Y = scipy.integrate.simps(B * ybar, lamcie)
    Z = scipy.integrate.simps(B * zbar, lamcie)
    XYZ = numpy.array([X,Y,Z])

    x = X / sum(XYZ)
    y = Y / sum(XYZ)
    z   = 1 - x - y

    Y   = 1
    X   = (Y/y) * x
    Z   = (Y/y) * z
    XYZ = numpy.array([X,Y,Z])

    # Matrix for Wide RGB D65 conversion
    XYZ2RGB = numpy.array([[ 1.656492, -0.354851, -0.255038],
                           [-0.707196,  1.655397,  0.036152],
                           [ 0.051713, -0.121364,  1.011530]])
    #Map XYZ to RGB
    RGB = numpy.dot(XYZ2RGB,XYZ)

    #adjust gamma value of RGB colour
    for i,color in enumerate(RGB):
        if color <= 0.0031308:
            RGB[i] = 12.92 * color
        else:
            RGB[i] = (1 + 0.055) * color ** (1 / 2.4) - 0.055

    # RGB = RGB / np.array([0.9505, 1., 1.0890])  #Scale so that Y of "white" (D65) is (0.9505, 1.0000, 1.0890)
    maxRGB = max(RGB.flatten())
    #Normalize to 1 if there are values above
    if maxRGB > 1: RGB = RGB / maxRGB
    #Clip negative values          
    RGB = RGB.clip(min=0)                       

    #Normalize to number of colors
    RGB = 255 * RGB                               
    return RGB

def cie():
    """
    Color matching functions. Columns are wavelength in nm, and xbar, ybar,
    and zbar, are the functions for R, G, and B, respectively.
    """
    lxyz = numpy.array([[380., 0.0014, 0.0000, 0.0065],
                        [385., 0.0022, 0.0001, 0.0105],
                        [390., 0.0042, 0.0001, 0.0201],
                        [395., 0.0076, 0.0002, 0.0362],
                        [400., 0.0143, 0.0004, 0.0679],
                        [405., 0.0232, 0.0006, 0.1102],
                        [410., 0.0435, 0.0012, 0.2074],
                        [415., 0.0776, 0.0022, 0.3713],
                        [420., 0.1344, 0.0040, 0.6456],
                        [425., 0.2148, 0.0073, 1.0391],
                        [430., 0.2839, 0.0116, 1.3856],
                        [435., 0.3285, 0.0168, 1.6230],
                        [440., 0.3483, 0.0230, 1.7471],
                        [445., 0.3481, 0.0298, 1.7826],
                        [450., 0.3362, 0.0380, 1.7721],
                        [455., 0.3187, 0.0480, 1.7441],
                        [460., 0.2908, 0.0600, 1.6692],
                        [465., 0.2511, 0.0739, 1.5281],
                        [470., 0.1954, 0.0910, 1.2876],
                        [475., 0.1421, 0.1126, 1.0419],
                        [480., 0.0956, 0.1390, 0.8130],
                        [485., 0.0580, 0.1693, 0.6162],
                        [490., 0.0320, 0.2080, 0.4652],
                        [495., 0.0147, 0.2586, 0.3533],
                        [500., 0.0049, 0.3230, 0.2720],
                        [505., 0.0024, 0.4073, 0.2123],
                        [510., 0.0093, 0.5030, 0.1582],
                        [515., 0.0291, 0.6082, 0.1117],
                        [520., 0.0633, 0.7100, 0.0782],
                        [525., 0.1096, 0.7932, 0.0573],
                        [530., 0.1655, 0.8620, 0.0422],
                        [535., 0.2257, 0.9149, 0.0298],
                        [540., 0.2904, 0.9540, 0.0203],
                        [545., 0.3597, 0.9803, 0.0134],
                        [550., 0.4334, 0.9950, 0.0087],
                        [555., 0.5121, 1.0000, 0.0057],
                        [560., 0.5945, 0.9950, 0.0039],
                        [565., 0.6784, 0.9786, 0.0027],
                        [570., 0.7621, 0.9520, 0.0021],
                        [575., 0.8425, 0.9154, 0.0018],
                        [580., 0.9163, 0.8700, 0.0017],
                        [585., 0.9786, 0.8163, 0.0014],
                        [590., 1.0263, 0.7570, 0.0011],
                        [595., 1.0567, 0.6949, 0.0010],
                        [600., 1.0622, 0.6310, 0.0008],
                        [605., 1.0456, 0.5668, 0.0006],
                        [610., 1.0026, 0.5030, 0.0003],
                        [615., 0.9384, 0.4412, 0.0002],
                        [620., 0.8544, 0.3810, 0.0002],
                        [625., 0.7514, 0.3210, 0.0001],
                        [630., 0.6424, 0.2650, 0.0000],
                        [635., 0.5419, 0.2170, 0.0000],
                        [640., 0.4479, 0.1750, 0.0000],
                        [645., 0.3608, 0.1382, 0.0000],
                        [650., 0.2835, 0.1070, 0.0000],
                        [655., 0.2187, 0.0816, 0.0000],
                        [660., 0.1649, 0.0610, 0.0000],
                        [665., 0.1212, 0.0446, 0.0000],
                        [670., 0.0874, 0.0320, 0.0000],
                        [675., 0.0636, 0.0232, 0.0000],
                        [680., 0.0468, 0.0170, 0.0000],
                        [685., 0.0329, 0.0119, 0.0000],
                        [690., 0.0227, 0.0082, 0.0000],
                        [695., 0.0158, 0.0057, 0.0000],
                        [700., 0.0114, 0.0041, 0.0000],
                        [705., 0.0081, 0.0029, 0.0000],
                        [710., 0.0058, 0.0021, 0.0000],
                        [715., 0.0041, 0.0015, 0.0000],
                        [720., 0.0029, 0.0010, 0.0000],
                        [725., 0.0020, 0.0007, 0.0000],
                        [730., 0.0014, 0.0005, 0.0000],
                        [735., 0.0010, 0.0004, 0.0000],
                        [740., 0.0007, 0.0002, 0.0000],
                        [745., 0.0005, 0.0002, 0.0000],
                        [750., 0.0003, 0.0001, 0.0000],
                        [755., 0.0002, 0.0001, 0.0000],
                        [760., 0.0002, 0.0001, 0.0000],
                        [765., 0.0001, 0.0000, 0.0000],
                        [770., 0.0001, 0.0000, 0.0000],
                        [775., 0.0001, 0.0000, 0.0000],
                        [780., 0.0000, 0.0000, 0.0000]])
    return lxyz.T

if __name__ == '__main__':
    main()