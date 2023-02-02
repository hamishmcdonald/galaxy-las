import numpy
import timeline
import matplotlib
import csv

temperatures = range(0, 15001, 100)
temperature_values, red_values, green_values, blue_values = [], [], [], []
for temperature in temperatures:

    red_value, green_value, blue_value = timeline.rgb_from_T(temperature, False, 255, False)

    round(red_value, 2)
    round(green_value, 2)
    round(blue_value, 2)

    print(str(temperature) + " : [" + str(red_value) + ", " + str(green_value) + ", " + str(blue_value) + "],")

    temperature_values.append(temperature)
    red_values.append(red_value)
    green_values.append(green_value)
    blue_values.append(blue_value)

#red higher than 5710
red_high_polynomial = numpy.polynomial.Polynomial.fit(temperature_values[1141:], red_values[1141:], 4)
#green in between 
green_lower_polynomial = numpy.polynomial.Polynomial.fit(temperature_values[133:1143], green_values[133:1143], 4)
#green higher than
green_higher_polynomial = numpy.polynomial.Polynomial.fit(temperature_values[1229:], green_values[1229:], 4)
#blue in between
blue_middle_polynomial = numpy.polynomial.Polynomial.fit(temperature_values[279:1231], blue_values[279:1231], 4) 
   
print(red_high_polynomial, green_lower_polynomial, green_higher_polynomial, blue_middle_polynomial)

#produce a matrix of values where each row corresponds to an x value in the data and each column corresponds to a 
#coefficient in the final formula (the total number of columns is degree + 1)
def polynomial_best_fit(x, y, degree):
    matrix = []
    for x_value in x:
        row = []
        for j in range(degree, 0, -1):
            row.append(x_value ^ j)
        row.append(int(1))
        matrix.append(row) 
