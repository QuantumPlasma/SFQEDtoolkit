import struct
import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams.update({'font.size': 28})

path = os.getcwd()
print('From: ' + path)

bytes_per_double = 8
bytes_per_int = 4

root_folder = './'
file_name = 'transition.dat'


################## single file reading
file = open(root_folder + file_name, 'rb')

number_of_points = struct.unpack('=I', file.read(bytes_per_int))[0]
print('1:',number_of_points)

bytes_to_read = bytes_per_double * 2 * number_of_points
#print(bytes_to_read)
buffer_bytes = file.read(bytes_to_read)
format = '=' + 'd' * number_of_points * 2
unpacked_data = struct.unpack(format,buffer_bytes)
#unpacked_data = struct.unpack_from(format,buffer_bytes,4)
#print(unpacked_data)

x_diff = np.array(unpacked_data[::2])
#print('x: ')
#for coeff in x_diff:
#    print(coeff)

y_diff = np.array(unpacked_data[1::2])
#print('y: ')
#for coeff in y_diff:
#    print(coeff)

print('rel diff value = ', y_diff[2])

x_diff = x_diff[1::]
y_diff = y_diff[1::]

#difference
#print plot
fig_rel, ax_rel = plt.subplots(figsize=(10,10))
ax_rel.plot(x_diff, y_diff,'g-',label='rel')

minimum = np.amin(y_diff)
maximumum = np.amax(y_diff)
quantity = (maximumum - minimum) / 5.
print(minimum, maximumum, quantity)
ax_rel.set_yticks(np.arange(minimum, maximumum + quantity, quantity))

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

#plt.title('Relative difference')
plt.xlabel(r'$\chi_e$')
plt.ylabel('y - axis')
#plt.legend()
plt.savefig(root_folder + 'plotted_rel', bbox_inches='tight')




plt.show()
file.close()
