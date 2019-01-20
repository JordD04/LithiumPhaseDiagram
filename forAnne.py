N_temp = 60
start_temp = 5
delta_temp = 5

def constants(N_temp, start_temp, delta_temp):
    N_temp = str(N_temp)
    start_temp = str(start_temp)
    delta_temp = str(delta_temp)

    constants_string = ''.join([N_temp, ' ', start_temp, ' ', delta_temp])

    return constants_string

constants_string = constants(N_temp, start_temp, delta_temp)

constants_list = constants_string.split(' ')

N_temp = int(constants_list[0])
start_temp = int(constants_list[1])
delta_temp = int(constants_list[2])

