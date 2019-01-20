import csv
from decimal import *
import math

# Boltzmann constant in eV / K
k_B = 0.00008617385
e = Decimal(2.71828)

# dictionary of parameters
#parameters = {
#    'Q4' : 2,
#    'W4' : 3,
#    'Q6' : 4,
#    'W6' : 5,
#    'Temp' : 6,
#    'Iteration' : 7
#}

inputFileStart = 'li_ko_5_f.all'
inputFile = ''.join([inputFileStart, '.qw46HV'])
outputFile = ''.join([inputFileStart, '.freeen'])

outputFileOpen = open(outputFile, 'w')


walkers = Decimal(1536)
walk_frac = walkers/(walkers+1)

N_temp = 60
start_temp = 5
delta_temp = 5

linesQW46HV = [line.rstrip('\n') for line in open(inputFile)]
noLines = len(linesQW46HV)-3

# creates a list of a list of all variables
jList = []
for line in linesQW46HV:
    varList = line.split(' ')
    jList.append(varList)
del jList[noLines+2]
del jList[noLines+1]
del jList[noLines]


# determines an energy modifier such that at the maximum temperature and iteration count the lowest partition function is ~1
def mod_finder(jList, noLines, N_temp, delta_temp, k_B, walk_frac):

    last_iteration = Decimal(jList[noLines-1][7])

    temp = start_temp+(N_temp*delta_temp)

    beta = Decimal(1 / (k_B * temp))
    w_i = (walk_frac**last_iteration) - ((walk_frac)**(last_iteration+1))
    boltzmann_factor_ideal = int(1/w_i)
    ener_mod = Decimal(math.log(boltzmann_factor_ideal))/(-beta)

    return ener_mod

ener_mod = mod_finder(jList, noLines, N_temp, delta_temp, k_B, walk_frac)


temp = start_temp
tempCount = 0

while tempCount < N_temp:                                                           # iterates through temperatures
    temp = start_temp + (delta_temp * tempCount)

    Z_fcc = 0
    Z_hhcc = 0
    Z_hcp = 0
    Z_bcc = 0
    Z_rest = 0
    Z_total = 0


    done = 0                                                                        # variable to check completion
    beta = Decimal(1 / (k_B * temp))

    # iterates through configurations
    for row in jList:
        print(row)
        if row[6] == '':
            done = 1
        else:
            if float(row[6]) > 1000 or done == 1:
                pass
            # finds the phase a point belongs to
            else:
                Q4 = float(row[2])
                W4 = float(row[3])
                Q6 = float(row[4])
                W6 = float(row[5])
                enthalpy = Decimal(row[0])+ener_mod
                iteration = Decimal(row[7])
                w_i = (walk_frac**iteration) - ((walk_frac)**(iteration+1))     # Nested Sampling Weight
                boltzmann_factor = e**(-beta*enthalpy)
                if Q4 < 0.2 and Q4 > 0.152:                                     # Tests which phase the point belongs too
                    Z_fcc = Z_fcc + (w_i * boltzmann_factor)                    # Adds partition function to total
                else:
                    if Q4 < 0.152 and Q4 > 0.132:
                        Z_hhcc = Z_hhcc + (w_i * boltzmann_factor)
                    else:
                        if Q4 < 20001 and Q4 > 20000:
                            Z_hcp = Z_hcp + (w_i * boltzmann_factor)
                        else:
                            if Q4 < 0.09345 and Q4 > 0.02:
                                Z_bcc = Z_bcc + (w_i * boltzmann_factor)
                            else:
                                Z_rest = Z_rest + (w_i * boltzmann_factor)


    # finds the log
    Z_fcc = int(Z_fcc)
    if Z_fcc > 0:
        ZLog_fcc = math.log(Z_fcc)
    else:
        ZLog_fcc = 0


    Z_hhcc = int(Z_hhcc)
    if Z_hhcc > 0:
        ZLog_hhcc = math.log(Z_hhcc)
    else:
        ZLog_hhcc = 0


    Z_hcp = int(Z_hcp)
    if Z_hcp > 0:
        ZLog_hcp = math.log(Z_hcp)
    else:
        ZLog_hcp = 0


    Z_bcc = int(Z_bcc)
    if Z_bcc > 0:
        ZLog_bcc = math.log(Z_bcc)
    else:
        ZLog_bcc = 0


    Z_rest = int(Z_rest)
    if Z_rest > 0:
        ZLog_rest = math.log(Z_rest)
    else:
        ZLog_rest = 0


    Z_total = Z_fcc + Z_hhcc + Z_hcp + Z_bcc + Z_rest

    Z_total = int(Z_total)
    if Z_total > 0:
        LogSum = math.log(Z_total)
    else:
        LogSum = 0

    # prints the output
    line = ''.join([str(temp), ' ', str(LogSum), ' ', str(ZLog_fcc), ' ', str(ZLog_hhcc), ' ', str(ZLog_hcp), ' ', str(ZLog_bcc), ' ', '-Infinity', ' ', str(ZLog_rest)])
    outputFileOpen.write(line + '\n')

    tempCount = tempCount + 1

outputFileOpen.close()