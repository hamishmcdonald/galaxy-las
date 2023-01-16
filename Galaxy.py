import math
import csv
import os
import laspy
import traceback
import timeline

#local file location of GaiaSource files
#the files can be found for download here:
#http://cdn.gea.esac.esa.int/Gaia/gedr3/gaia_source/
file_directory = 'D:\GaiaSource'

#constants and their units
WEIN_CONSTANT = 2897771.9 #nM K
BOLTZMANN_CONSTANT = 1.380649e-23 #J/(K)
PLANCK_CONSTANT = 6.62607015e-34 #J s
LIGHT_SPEED_CONSTANT = 29979245800 #cm/(s)

#matrix of polynomial co-efficients for temperature to rgb conversion
POLYNOMIAL_COEFFICIENTS = [[202.7407364067034, -26.810279254200008, 17.92772958298324, - 9.45648666312952, 1.6236384714807253], 
[218.60561513620232, 52.17192901392309, -23.061733295188585, 54.04865637275284, -51.058826067299464], 
[216.63554612112824, -19.0450479334475, 9.609252176140112, -6.312170402591599, 3.146919073190075], 
[177.94711834639443, 110.39882919514417, -34.69680828555997, 13.575839801777828, -12.414595228016992]]

def main():
    #iterate through GaiaSource files and open them provided they are csv files
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
                   #more information on the fieldnames used in the GaiaSource files can be found here:
                   #https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/
                   #sec_dm_main_source_catalogue/ssec_dm_gaia_source.html
                    current_csv_reader = csv.DictReader(current_csv)

                    row_number = 1
                    largest_temperature = 0

                    #iterate through each star in the current GaiaSource file
                    for row in current_csv_reader:
                        try:
                            row_number += 1

                            #calculate cartesian coordinates and rgb colorisation of star (row)
                            x_value, y_value, z_value = calculateCartesian(row)
                            red_value, green_value, blue_value = calculateRGB(calculateTemperature(row))

                            #store unique source indentifiers and designations of star (row)
                            solution_id_value = int(row['solution_id'])
                            designation_value = int(row['designation'][11:])
                            source_id_value = int(row['source_id'])

                            #add values to associated temporary array, performed separately in case calculations produce an 
                            #exception
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

#calculate x, y, z coordinates of the star using parallax, galactic longitude and latitude
#more information on the formulas can be found here:
#https://en.wikipedia.org/wiki/Galactic_coordinate_system
def calculateCartesian(row):
    if row['parallax'] == '':
        raise Exception("no parallax")
    else:
        x_value = math.cos(float(row['b'])) * math.cos(float(row['l'])) / float(row['parallax'])
        y_value = math.cos(float(row['b'])) * math.sin(float(row['l'])) / float(row['parallax'])
        z_value = math.sin(float(row['b'])) / float(row['parallax'])
    
    return x_value, y_value, z_value

#calculate peak wavelength of light emitted from the star then calculate temperature of the star with Wien's law formula 
#using displacement constant and peak wavelength
#https://www.omnicalculator.com/physics/wiens-law
def calculateTemperature(row):
    if not row['nu_eff_used_in_astrometry'] == '':
        peak_wavelength_value = 1 / float(row['nu_eff_used_in_astrometry']) * 1000
    elif not row['pseudocolour'] == '':
        peak_wavelength_value = 1 / float(row['pseudocolour']) * 1000
    else:
        raise Exception("no nu_eff_used_in_astronomy or pseudocolour")

    temperature = WEIN_CONSTANT / peak_wavelength_value

    return temperature

#calculate rgb values of star using temperature
#https://en.wikipedia.org/wiki/CIE_1931_color_space
# polynomial equations approximating data collected from timeline.py using getRGB.py
def calculateRGB(t):
    if 0 <= t <= 15000:

        #calculate red value
        if t <= 5705:
            red_value = 255
        elif 5705 < t:
            red_value = calculatePolynomial(t, 0)

        #calcluate green value
        if t <= 665:
            green_value = 0
        elif 665 < t <= 5705:
            green_value = calculatePolynomial(t, 1)
        elif 5705 < t <= 6145:
            green_value = 255
        elif 6145 < t:
            green_value = calculatePolynomial(t, 2)

        #calculate blue value
        if t <= 1395:
            blue_value = 0
        elif 1395 < t <= 6145:
            blue_value = calculatePolynomial(t, 3)
        elif 6145 < t:
            blue_value = 255

    else:
        raise Exception("temperature outside of normal range: " + str(t))

    rgb_value = red_value, green_value, blue_value

    #normalise rgb value
    for value in rgb_value:
        if value > 255:
            value = 255
        if value < 0:
            value = 0
    
    return rgb_value

#calculate polynomial
def calculatePolynomial(t, polynomial):
    colour_value = 0
    for i, coefficient in enumerate(POLYNOMIAL_COEFFICIENTS[polynomial]):
        colour_value += coefficient * t ** i

    return colour_value

if __name__ == '__main__':
    main()