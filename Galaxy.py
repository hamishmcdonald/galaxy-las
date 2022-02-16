import math
import csv
import os
import laspy
import traceback
import timeline

def main():
    file_directory = 'C:\\Users\\HamishMcDonald\\OneDrive - Euclideon PTY LTD\\Desktop\\gaia_source'

    for gaia_file in os.listdir(file_directory):
        if gaia_file.endswith('.csv'):
            try:
                with open(file_directory + '/' + gaia_file) as current_csv:
                
                    header = laspy.LasHeader(version="1.4", point_format=2)
                    galaxy_data = laspy.LasData(header)

                    galaxy_data.add_extra_dims([
                    laspy.ExtraBytesParams(name="peak_wavelength_param", type="float64"),
                    laspy.ExtraBytesParams(name="nu_eff_used_in_astronomy", type="float64"),
                    laspy.ExtraBytesParams(name="pseudocolour", type="float64"), 
                    laspy.ExtraBytesParams(name="rgb", type="3int16"),
                    laspy.ExtraBytesParams(name="solution_id", type="uint64"),
                    laspy.ExtraBytesParams(name="designation", type="uint64"),
                    laspy.ExtraBytesParams(name="source_id", type="uint64"),
                    laspy.ExtraBytesParams(name="parallax", type="float64"),
                    laspy.ExtraBytesParams(name="pm", type="float64"),
                    laspy.ExtraBytesParams(name="phot_g_mean_mag", type="float64"),
                    laspy.ExtraBytesParams(name="phot_g_mean_flux", type="float64"),
                    ])

                    x, y, z = [], [], []
                    red, green, blue = [], [], []
                    solution_id, designation, source_id= [], [], []
                    peak_wavelength_array, nu_eff_used_in_astronomy, pseudocolour = [], [], []
                    parallax, pm, rgb = [], [], []
                    phot_g_mean_mag, phot_g_mean_flux  = [], []
                   
                    current_csv_reader = csv.DictReader(current_csv)

                    row_number = 0

                    for row in current_csv_reader:
                        try:
                            row_number += 1

                            x_value = math.cos(float(row['b'])) * math.cos(float(row['l'])) / float(row['parallax'])
                            y_value = math.cos(float(row['b'])) * math.sin(float(row['l'])) / float(row['parallax'])
                            z_value = math.sin(float(row['b'])) / float(row['parallax'])

                            if not row['nu_eff_used_in_astrometry'] == '':
                                peak_wavelength_value = 1 / float(row['nu_eff_used_in_astrometry']) * 1000
                                nu_eff_used_in_astronomy_value = float(row['nu_eff_used_in_astrometry'])
                                pseudocolour_value = None
                            elif not row['pseudocolour'] == '':
                                peak_wavelength_value = 1 / float(row['pseudocolour']) * 1000
                                nu_eff_used_in_astronomy_value = None
                                pseudocolour_value = float(row['pseudocolour'])
                            else:
                                raise Exception("no nu_eff_used_in_astronomy or pseudocolour")
                         
                            temperature = 2897771.9 / peak_wavelength_value
                            
                            red_value, green_value, blue_value = timeline.rgb_from_T(temperature, False, 255, False)

                            rgb_value = [red_value, green_value, blue_value]

                            solution_id_value = int(row['solution_id'])
                            designation_value = int(row['designation'][11:])
                            source_id_value = int(row['source_id'])
                            parallax_value = float(row['parallax'])
                            pm_value = float(row['pm'])
                            phot_g_mean_mag_value = float(row['phot_g_mean_mag'])
                            phot_g_mean_flux_value= float(row['phot_g_mean_flux'])

                            x.append(x_value)
                            y.append(y_value)
                            z.append(z_value)
                            peak_wavelength_array.append(peak_wavelength_value)
                            nu_eff_used_in_astronomy.append(nu_eff_used_in_astronomy_value)
                            pseudocolour.append(pseudocolour_value)
                            red.append(red_value)
                            green.append(green_value)
                            blue.append(blue_value)
                            rgb.append(rgb_value)
                            solution_id.append(solution_id_value)
                            designation.append(designation_value)
                            source_id.append(source_id_value)
                            parallax.append(parallax_value)
                            pm.append(pm_value)
                            phot_g_mean_mag.append(phot_g_mean_mag_value)
                            phot_g_mean_flux.append(phot_g_mean_flux_value)
                            
                            print("data in row " + str(row_number) + " successfully added")

                        except Exception as row_exception:
                            print("Exception occured in row " + str(row_number) + " in file " + gaia_file + ": ", row_exception)
                            traceback.print_exc()

                galaxy_data.x = x
                galaxy_data.y = y
                galaxy_data.z = z
                galaxy_data.peak_wavelength_param = peak_wavelength_array
                galaxy_data.nu_eff_used_in_astronomy = nu_eff_used_in_astronomy
                galaxy_data.pseudocolour = pseudocolour
                galaxy_data.red = red
                galaxy_data.green = green
                galaxy_data.blue = blue
                galaxy_data.rgb = rgb
                galaxy_data.solution_id = solution_id
                galaxy_data.designation = designation
                galaxy_data.source_id = source_id
                galaxy_data.parallax = parallax
                galaxy_data.pm = pm
                galaxy_data.phot_g_mean_mag = phot_g_mean_mag
                galaxy_data.phot_g_mean_flux = phot_g_mean_flux

                galaxy_data.write(gaia_file + ".las")

                print(gaia_file + " was successfully converted")

            except Exception as file_exception:
                print("Exception occured in file " + gaia_file + ": ", file_exception)
                traceback.print_exc()

if __name__ == '__main__':
    main()