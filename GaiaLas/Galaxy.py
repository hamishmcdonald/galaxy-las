import math
import csv
import os
import laspy
import traceback

#local file location of GaiaSource files
#the files can be found for download here:
#http://cdn.gea.esac.esa.int/Gaia/gedr3/gaia_source/
file_directory = 'E:\GaiaSource'

#constants and their units
WEIN_CONSTANT = 2897771.9 #nM K
BOLTZMANN_CONSTANT = 1.380649e-23 #J/(K)
PLANCK_CONSTANT = 6.62607015e-34 #J s
LIGHT_SPEED_CONSTANT = 29979245800 #cm/(s)

#matrix of polynomial co-efficients for temperature to rgb conversion
POLYNOMIAL_COEFFICIENTS = [
    [202.7407364067034, -26.810279254200008, 17.92772958298324, -9.45648666312952, 1.6236384714807253], 
    [218.60561513620232, 52.17192901392309, -23.061733295188585, 54.04865637275284, -51.058826067299464], 
    [216.63554612112824, -19.0450479334475, 9.609252176140112, -6.312170402591599, 3.146919073190075], 
    [177.94711834639443, 110.39882919514417, -34.69680828555997, 13.575839801777828, -12.414595228016992]
    ]

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
                            #red_value, green_value, blue_value = calculateRGB(calculateTemperature(row))
                            red_value, green_value, blue_value = retrieveRGB(calculateTemperature(row))

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
                            #pass

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
                galaxy_data.write(gaia_file + '.las')

                #print to console if no exceptions occured for the GaiaSource file 
                print(gaia_file + " was successfully converted")

            #print to console if an exception occured for the GaiaSource file
            except Exception as file_exception:
                print(gaia_file + ": ", file_exception)
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

        print(t, red_value, green_value, blue_value)
    else:
        raise Exception("temperature outside of normal range: " + str(t))

    rgb_value = [red_value, green_value, blue_value]

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

def retrieveRGB(t):
    #round temperature to the nearest 100
    t = round(t / 100) * 100

    RGB = {
        0 : [255.0, 0.0, 0.0],
        100 : [255.0, 0.0, 0.0],
        200 : [255.0, 0.0, 0.0],
        300 : [255.0, 0.0, 0.0],
        400 : [255.0, 0.0, 0.0],
        500 : [255.0, 0.0, 0.0],
        600 : [255.0, 0.0, 0.0],
        700 : [255.0, 30.707195158993482, 0.0],
        800 : [255.0, 62.84631865698403, 0.0],
        900 : [255.0, 82.38814221954293, 0.0],
        1000 : [255.0, 97.57751272824949, 0.0],
        1100 : [255.0, 110.27819077186376, 0.0],
        1200 : [255.0, 121.27772322322069, 0.0],
        1300 : [255.0, 131.00195842275375, 0.0],
        1400 : [255.0, 139.71586404984836, 0.11050373094119768],
        1500 : [255.0, 147.59994831496417, 18.644786247062925],
        1600 : [255.0, 154.78513693782884, 30.60807369940283],
        1700 : [255.0, 161.37080136250643, 40.77624591188531],
        1800 : [255.0, 167.43497226726078, 50.10731982367809],
        1900 : [255.0, 173.0405484608862, 58.91408920062816],
        2000 : [255.0, 178.23928830746172, 67.33215692226692],
        2100 : [255.0, 183.07449467341854, 75.42990482515403],
        2200 : [255.0, 187.5828900300578, 83.2456894518547],
        2300 : [255.0, 191.79596786796577, 90.80309972758808],
        2400 : [255.0, 195.74099313122167, 98.11800548192707],
        2500 : [255.0, 199.4417600895473, 105.20207706752826],
        2600 : [255.0, 202.91917804363274, 112.06463602929267],
        2700 : [255.0, 206.19173192828404, 118.71366101565178],
        2800 : [255.0, 209.2758500985908, 125.15634478698563],
        2900 : [255.0, 212.18620195725364, 131.39940437782826],
        3000 : [255.0, 214.9359416544312, 137.44925243526643],
        3100 : [255.0, 217.53690970629796, 143.31208953713326],
        3200 : [255.0, 219.99980132690428, 148.99395145762813],
        3300 : [255.0, 222.33430810591696, 154.50073101398092],
        3400 : [255.0, 224.54923810751725, 159.83818595519116],
        3500 : [255.0, 226.65261832658052, 165.01193959663362],
        3600 : [255.0, 228.65178259284136, 170.02747809396328],
        3700 : [255.0, 230.55344737767754, 174.89014657445878],
        3800 : [255.0, 232.3637774733693, 179.60514534222736],
        3900 : [255.0, 234.08844314068605, 184.17752677696132],
        4000 : [255.0, 235.73267002876855, 188.612193194808],
        4100 : [255.0, 237.3012829410251, 192.91389573834647],
        4200 : [255.0, 238.79874433730635, 197.0872342520505],
        4300 : [255.0, 240.22918831510043, 201.13665804392585],
        4400 : [255.0, 241.59645069285713, 205.06646741095014],
        4500 : [255.0, 242.9040957207876, 208.88081580173264],
        4600 : [255.0, 244.1554398640468, 212.5837124959702],
        4700 : [255.0, 245.35357303659254, 216.17902569159276],
        4800 : [255.0, 246.50137760857132, 219.67048590390402],
        4900 : [255.0, 247.60154546366704, 223.0616895946868],
        5000 : [255.0, 248.65659334385973, 226.35610296219053],
        5100 : [255.0, 249.66887668612821, 229.5570658346544],
        5200 : [255.0, 250.64060212776877, 232.66779562033054],
        5300 : [255.0, 251.57383883332682, 235.69139127587437],
        5400 : [255.0, 252.47052877597042, 238.63083726250522],
        5500 : [255.0, 253.33249608890205, 241.48900746566997],
        5600 : [255.0, 254.16145558764882, 244.26866905919266],
        5700 : [255.0, 254.95902055140263, 246.9724862992151],
        5800 : [254.2753552826387, 255.0, 248.89371634285885],
        5900 : [253.54242494781823, 255.0, 250.7213949411821],
        6000 : [252.84034354247356, 255.0, 252.49732151851293],
        6100 : [252.16738962434155, 255.0, 254.22349337372583],
        6200 : [250.6355871967779, 254.1013696493253, 255.0],
        6300 : [248.43374496729214, 252.49086341010553, 255.0],
        6400 : [246.32601853230557, 250.9435715487046, 255.0],
        6500 : [244.3068069034624, 249.45598408533232, 255.0],
        6600 : [242.3709347025114, 248.02484400407218, 255.0],
        6700 : [240.51361283542045, 246.64712476575556, 255.0],
        6800 : [238.73040343768372, 245.3200101952792, 255.0],
        6900 : [237.01718855947016, 244.04087645426077, 255.0],
        7000 : [235.37014213345006, 242.80727584959166, 255.0],
        7100 : [233.78570483088797, 241.61692226209115, 255.0],
        7200 : [232.26056146481875, 240.4676780080743, 255.0],
        7300 : [230.79162064441437, 239.35754197105797, 255.0],
        7400 : [229.3759964232882, 238.2846388617041, 255.0],
        7500 : [228.01099171753992, 237.24720948201144, 255.0],
        7600 : [226.694083297704, 236.24360188516292, 255.0],
        7700 : [225.42290818314427, 235.27226333571403, 255.0],
        7800 : [224.19525128846723, 234.33173298628398, 255.0],
        7900 : [223.00903418968545, 233.42063519684277, 255.0],
        8000 : [221.862304893604, 232.53767343132546, 255.0],
        8100 : [220.75322850755964, 231.68162467380702, 255.0],
        8200 : [219.6800787185275, 230.8513343130223, 255.0],
        8300 : [218.64123000097175, 230.04571144973636, 255.0],
        8400 : [217.635150481867, 229.26372458647973, 255.0],
        8500 : [216.66039539924398, 228.50439766356135, 255.0],
        8600 : [215.71560109756032, 227.76680640913128, 255.0],
        8700 : [214.7994795093067, 227.05007497447355, 255.0],
        8800 : [213.91081307763326, 226.3533728287053, 255.0],
        8900 : [213.04845007952545, 225.67591188971792, 255.0],
        9000 : [212.21130031324526, 225.01694387054022, 255.0],
        9100 : [211.39833111746267, 224.37575782238872, 255.0],
        9200 : [210.60856369278875, 223.75167785752168, 255.0],
        9300 : [209.84106969933717, 223.14406103665712, 255.0],
        9400 : [209.0949681065388, 222.55229540718528, 255.0],
        9500 : [208.36942227374226, 221.97579817971243, 255.0],
        9600 : [207.6636372421942, 221.4140140316449, 255.0],
        9700 : [206.976857220837, 220.86641352756791, 255.0],
        9800 : [206.30836325000283, 220.3324916471156, 255.0],
        9900 : [205.65747102856437, 219.81176641186545, 255.0],
        10000 : [205.02352889141898, 219.3037776035543, 255.0],
        10100 : [204.40591592537618, 218.80808556658914, 255.0],
        10200 : [203.80404021258323, 218.3242700884425, 255.0],
        10300 : [203.21733719158814, 217.85192935207598, 255.0],
        10400 : [202.64526812700373, 217.3906789550351, 255.0],
        10500 : [202.08731867952264, 216.94015099031103, 255.0],
        10600 : [201.5429975687355, 216.49999318447584, 255.0],
        10700 : [201.0118353218505, 216.06986808896923, 255.0],
        10800 : [200.49338310198542, 215.64945232074794, 255.0],
        10900 : [199.98721161023272, 215.23843584882223, 255.0],
        11000 : [199.492910056175, 214.83652132347567, 255.0],
        11100 : [199.0100851919579, 214.44342344522153, 255.0],
        11200 : [198.5383604054254, 214.05886837078003, 255.0],
        11300 : [198.07737486817993, 213.6825931535697, 255.0],
        11400 : [197.626782734753, 213.3143452164002, 255.0],
        11500 : [197.18625238937673, 212.95388185422908, 255.0],
        11600 : [196.75546573711335, 212.60096976500648, 255.0],
        11700 : [196.33411753635022, 212.2553846067787, 255.0],
        11800 : [195.9219147698959, 211.9169105793561, 255.0],
        11900 : [195.51857605212035, 211.58534002897636, 255.0],
        12000 : [195.12383106977464, 211.26047307450602, 255.0],
        12100 : [194.7374200542964, 210.9421172538295, 255.0],
        12200 : [194.35909328357377, 210.63008718916896, 255.0],
        12300 : [193.988610611283, 210.32420427016888, 255.0],
        12400 : [193.62574102205264, 210.02429635366, 255.0],
        12500 : [193.27026221083278, 209.73019747909075, 255.0],
        12600 : [192.92196018495915, 209.44174759868685, 255.0],
        12700 : [192.58062888751397, 209.15879232146168, 255.0],
        12800 : [192.246069840674, 208.8811826702582, 255.0],
        12900 : [191.9180918078358, 208.60877485106138, 255.0],
        13000 : [191.59651047338357, 208.3414300338663, 255.0],
        13100 : [191.28114813904722, 208.07901414443802, 255.0],
        13200 : [190.97183343586528, 207.82139766633978, 255.0],
        13300 : [190.66840105083475, 207.56845545264727, 255.0],
        13400 : [190.37069146739282, 207.32006654680558, 255.0],
        13500 : [190.07855071892533, 207.07611401211636, 255.0],
        13600 : [189.7918301545582, 206.83648476937975, 255.0],
        13700 : [189.51038621652717, 206.6010694422414, 255.0],
        13800 : [189.2340802284738, 206.36976220982461, 255.0],
        13900 : [188.96277819405177, 206.14246066625438, 255.0],
        14000 : [188.6963506052702, 205.9190656867019, 255.0],
        14100 : [188.43467226003446, 205.69948129960198, 255.0],
        14200 : [188.1776220883793, 205.48361456471866, 255.0],
        14300 : [187.9250829869209, 205.27137545674762, 255.0],
        14400 : [187.67694166108183, 205.0626767541714, 255.0],
        14500 : [187.43308847467245, 204.85743393309167, 255.0],
        14600 : [187.19341730643555, 204.6555650657848, 255.0],
        14700 : [186.95782541318457, 204.45699072373785, 255.0],
        14800 : [186.72621329918988, 204.26163388493944, 255.0],
        14900 : [186.49848459148575, 204.06941984520992, 255.0],
        15000 : [186.27454592079093, 203.88027613336965, 255.0]
    }
  
    return RGB.get(t)

if __name__ == '__main__':
    main()