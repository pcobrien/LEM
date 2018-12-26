import numpy as np
import matplotlib.pyplot as plt
import sys

def grid_diffusion_periodic(z0, z, D, dt, dx2, dy2):
	# Compute diffusion for grid using input parameters
	
	grid_size = z.shape[0]
	ind = range(0, grid_size)
	ind_up = np.roll(ind, -1)
	ind_down = np.roll(ind, 1)
	
	z[0:,0:] = z0[0:,0:] + D * dt * ( (z0[ind_up, 0:] - 2*z0[0:, 0:] + z0[ind_down, 0:])/dx2 + (z0[0:, ind_up] - 2*z0[0:, 0:] + z0[0:, ind_down])/dy2 )

	# Replace the "old" grid with the post-diffusion grid
	z0 = z.copy()
	
	return z0, z


def grid_diffusion(z0, z, D, dt, dx2, dy2):
	z[1:-1, 1:-1] = z0[1:-1, 1:-1] + D * dt * ( (z0[2:, 1:-1] - 2*z0[1:-1, 1:-1] + z0[:-2, 1:-1])/dx2 + (z0[1:-1, 2:] - 2*z0[1:-1, 1:-1] + z0[1:-1, :-2])/dy2 )
	z[0,:] = z[1,:]
	z[:,0] = z[:,1]
	z[-1,:] = z[-2,:]
	z[:,-1] = z[:,-2]
	z0 = z.copy()
	return z0, z


def crater_richardson(xcenter_ind, ycenter_ind, radius, grid, dgrid_x, dgrid_y):
	# function to generate parabolic crater as outlined in Richardson 2009
	current_grid = grid.copy()

	# Set constants of crater profile
	diameter = 2.0*radius
	depth = 0.2*diameter
	rim_height = 0.2*depth
	
	#print 'Depth: ', depth
	#print 'Rim height: ', rim_height
	
	crater_volume = (np.pi/2.0)*(depth)*(radius**2)		# May be used later to constrain ejecta blanket profile
	
	ejecta_limit = 0.01	# Arbitrary limit for height at which the ejecta blanket is cut off
	ejecta_volume = 0
	r_ejecta = ( (ejecta_limit/rim_height)**(-1.0/3.0))*radius
	
	grid_size = grid.shape[0]
	center = int(grid_size/2.0)

	# Set up grid and calculate the distance from the center of the crater to everywhere on the grid
	Y, X = np.ogrid[:grid_size, :grid_size]
	
	dx = abs(X - xcenter_ind)
	dy = abs(Y - ycenter_ind)
	dx[dx > grid_size/2.0] = grid_size - dx[dx > grid_size/2.0]
	dy[dy > grid_size/2.0] = grid_size - dy[dy > grid_size/2.0]
	dx = dx*dgrid_x
	dy = dy*dgrid_y
	dist_from_center = np.hypot(dx, dy)
	
	# Grid pixels covered by the crater
	crater_mask = dist_from_center <= radius

	# Grid pixels covered by the ejecta blanket
	ejecta_mask = (dist_from_center > radius) & (dist_from_center <= r_ejecta)

	# Inheritance parameter, set as a constant, see Howard 2007, I=0 --> crater rim horizontal, I=1 --> crater rim parallel to pre-existing slope
	I_i = 0.9

	weighted_grid = grid.copy()
	outside_mask = dist_from_center > radius
	weighted_grid[outside_mask] = current_grid[outside_mask]*(1.0/dist_from_center[outside_mask])
	E_r = np.average(weighted_grid[crater_mask]) + np.average(weighted_grid[~crater_mask])
	
	# Crater elevation profile
	delta_H_crater = (( ( ( (dx)**2 + (dy)**2 ) * (rim_height + depth) ) / radius**2 ) - depth)
	
	# Ejecta elevation profile
	with np.errstate(divide='ignore'):		# Divide by zero at r=0 but we don't care about that point
		delta_H_ejecta = rim_height*( ( (np.hypot(dx, dy))/(radius))**(-3.0))
	
	# Lognormal noise elevation profile
	noise_profile = np.random.lognormal(size=(grid_size, grid_size))/100.0
	
	# Inheritance matrices - determines how the crater is integrated into the existing grid
	G_grid = (1.0 - I_i)*np.ones((grid_size, grid_size))
	min_mask = G_grid > delta_H_ejecta/rim_height
	G_grid[min_mask] = delta_H_ejecta[min_mask]/rim_height
	
	crater_inh_profile = (1.0 - I_i*((dist_from_center/radius)**2.0))*(E_r - current_grid)
	ejecta_inh_profile = G_grid*(E_r - current_grid)
	
	delta_E_crater = delta_H_crater + crater_inh_profile + noise_profile
	delta_E_ejecta = delta_H_ejecta + ejecta_inh_profile + noise_profile
	
	# Add calculated elevations to the grid at the corresponding pixels
	grid[crater_mask] += delta_E_crater[crater_mask]
	grid[ejecta_mask] += delta_E_ejecta[ejecta_mask]
	
	return grid

'''
D = 5.5e-5
D2 = 2.0*D

grid_size = 512
center = int(grid_size/2)
grid_width = 10000.0
dgrid_x = grid_width/grid_size
dgrid_y = grid_width/grid_size
dx2 = dgrid_x**2
dy2 = dgrid_y**2
dt = (dx2 * dy2 / (2 * D * (dx2 + dy2)))
dt2 = (dx2 * dy2 / (2 * D2 * (dx2 + dy2)))
model_time = 4e9 # 1 Gyr
nsteps = int(model_time/dt)
print nsteps
nsteps = 10

grid = np.zeros((grid_size, grid_size))
z0 = grid
z = np.empty((grid_size, grid_size))

z0_2 = z0.copy()
z_2 = z.copy()

z0 = crater_richardson(center, center, 2500.0, z0, dgrid_x, dgrid_y)
z0_2 = crater_richardson(center, center, 2500.0, z0, dgrid_x, dgrid_y)

slice = range(0, 512)
slice_plot = (slice-center*np.ones(len(slice)))*dgrid_x

for i in range(nsteps):
	z0, z = grid_diffusion(z0, z, D, dt2, dx2, dy2)
	z0_2, z_2 = grid_diffusion(z0_2, z_2, D2, dt2, dx2, dy2)
	
	plt.figure()
	plt.plot(slice_plot, z0[256, slice])
	plt.plot(slice_plot, z0_2[256, slice])
	plt.xlabel('Distance from crater center')
	plt.ylabel('Elevation')
	plt.title('Cross-section profile after %6.3f Myr' % int(i*dt2/1e6))
	fname = 'crater_evol_5km_%03d.png' % i
	plt.savefig(fname)
sys.exit()

# Plot a slice across the middle of the grid
plt.figure()
#plt.plot(slice_plot, fresh_slice)
plt.plot(slice_plot, z0[256, slice])
plt.plot(slice_plot, z0_2[256, slice])
plt.xlabel('Distance from crater center')
plt.ylabel('Elevation')
plt.title('Cross-section profile')
plt.show()
'''


