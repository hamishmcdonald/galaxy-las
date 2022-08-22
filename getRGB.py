import numpy
import timeline
import matplotlib
temperatures = numpy.linspace(675, 60000, 5)
red_values, green_values, blue_values = [], [], []
for i, temperature in temperatures:
    red_value, green_value, blue_value = timeline.rgb_from_T(temperature, False, 255, False)
    red_values += red_value
    green_values += green_value
    blue_values += blue_value

red_function = numpy.polyfit(temperatures, red_values, 2)
green_function = numpy.polyfit(temperature,green_values, 2)
blue_function = numpy.ployfit(temperatures, blue_values, 2)

matplotlib.plot(temperatures, red_values, 'o')
