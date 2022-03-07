import struct
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams.update({'font.size': 22})

bytes_per_double = 8
bytes_per_int = 4
root_folder = './'
file_name = 'transition.dat'
data_ratio = 4./5.

################################## single file reading
file = open(root_folder + file_name, 'rb')
#file = open(root_folder + 'rel_diff_nrgs' + file_suffix, 'rb')

points_x = struct.unpack('=I', file.read(bytes_per_int))[0]
(a, b) = struct.unpack('=dd', file.read(bytes_per_double * 2))

points_y = struct.unpack('=I', file.read(bytes_per_int))[0]
(c, d) = struct.unpack('=dd', file.read(bytes_per_double * 2))

#used to reshape
new_shape = (points_x, points_y)

total_points = points_x * points_y

bytes_to_read = bytes_per_double * total_points
buffer_bytes = file.read(bytes_to_read)
format = '=' + 'd' * total_points
unpacked_data = struct.unpack(format,buffer_bytes)

func_diff = np.reshape(np.array(unpacked_data), new_shape)
#remove 0 evals in both direction (notice we have not transposed yet)
func_diff = func_diff[1::,1::]


###################################### plotting
fig2, ax3 = plt.subplots(1,1,figsize=(20,10))

field_min = func_diff.min()
field_max = func_diff.max()
quantity = (field_max - field_min) / 7.
print('rel diff NaN check: ', np.isnan(field_min) or np.isnan(field_max))

# print(func_diff.shape)
# count = 0
# for ciccio in range(0,199):
#     for pippo in range(0,199):
#         if np.isnan(func_diff[ciccio,pippo]):
#             count = count + 1
#             #print(func_diff[ciccio,pippo])

# print(count, ' / ', 199*199)

a2 = ax3.imshow(func_diff.transpose(), extent=[a, b, c, d], origin='lower', cmap='gist_rainbow')

ax3.set_aspect(data_ratio/ax3.get_data_ratio())
plt.colorbar(a2, ticks=np.arange(field_min, field_max + quantity, quantity))

#divider = make_axes_locatable(ax3)
#cax = divider.append_axes("right", size="5%", pad=0.10)
#plt.colorbar(a2, cax=cax, ticks=np.arange(field_min, field_max + quantity, quantity))

# cb1 = fig.colorbar(a1, ax=ax2)
# cb1.ax.tick_params(axis='both', which='major', pad=15)
ax3.set_title(r'$F$')
ax3.set_xlabel(r'$\chi_\gamma$')
ax3.set_ylabel(r'$random\,\, number$')

#plt.zlabel('z - axis')
# plt.legend()
plt.savefig(root_folder + 'relative.png', bbox_inches='tight')


plt.show()
