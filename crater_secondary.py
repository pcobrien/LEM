import sys
import numpy as np
import matplotlib.pyplot as plt


def crater_secondary(xcenter_ind, ycenter_ind, diameter, grid, dgrid_x, dgrid_y):
	# function to generate parabolic crater as outlined in Richardson 2009, this function for craters whose centers are off-grid but some part of which overlaps the grid
	
	#if (0 <= xcenter_ind <= grid.shape[0]) & (0 <= xcenter_ind <= grid.shape[0]):
	current_grid = grid.copy()

	# Set constants of crater profile
	radius = diameter/2.0
	depth = 0.11*diameter
	rim_height = 0.2*depth

	crater_volume = (np.pi/2.0)*(depth)*(radius**2)		# May be used later to constrain ejecta blanket profile

	ejecta_limit = 0.01	# Arbitrary limit for height at which the ejecta blanket is cut off
	ejecta_volume = 0
	r_ejecta = ( (ejecta_limit/rim_height)**(-1.0/3.0))*radius

	grid_size = grid.shape[0]

	# Set up grid and calculate the distance from the center of the crater to everywhere on the grid
	Y, X = np.ogrid[:grid_size, :grid_size]

	dx = abs(X - xcenter_ind)
	dy = abs(Y - ycenter_ind)
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

	ret_grid = grid.copy()

	return ret_grid

#####
def crater_secondary_small(xcenter_ind, ycenter_ind, diameter, grid, dgrid_x, dgrid_y):
	# function to generate parabolic crater as outlined in Richardson 2009, this function for craters whose centers are off-grid but some part of which overlaps the grid
	
	#if (0 <= xcenter_ind <= grid.shape[0]) & (0 <= xcenter_ind <= grid.shape[0]):
	current_grid = grid.copy()

	# Set constants of crater profile
	radius = diameter/2.0
	depth = 0.11*diameter
	rim_height = 0.2*depth

	crater_volume = (np.pi/2.0)*(depth)*(radius**2)		# May be used later to constrain ejecta blanket profile

	ejecta_limit = 0.001	# Arbitrary limit for height at which the ejecta blanket is cut off
	ejecta_volume = 0
	r_ejecta = ( (ejecta_limit/rim_height)**(-1.0/3.0))*radius
	r_ejecta = 3.0*radius

	grid_size = grid.shape[0]

	# Set up grid and calculate the distance from the center of the crater to everywhere on the grid
	Y, X = np.ogrid[:grid_size, :grid_size]

	dx = abs(X - xcenter_ind)
	dy = abs(Y - ycenter_ind)
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

	delta_E_crater = delta_H_crater + crater_inh_profile# + noise_profile
	delta_E_ejecta = delta_H_ejecta + ejecta_inh_profile# + noise_profile

	# Add calculated elevations to the grid at the corresponding pixels
	grid[crater_mask] += delta_E_crater[crater_mask]
	grid[ejecta_mask] += delta_E_ejecta[ejecta_mask]

	ret_grid = grid.copy()

	return ret_grid

'''
grid_size = 512
z0 = np.zeros((512, 512))
z = z0.copy()

res = 19.53125
z0 = crater_secondary(256, 256, 500.0, z0, res, res)


z02 = np.zeros((512, 512))
z2 = z0.copy()

res = 19.53125
z02 = crater_secondary_small(256, 256, 500.0, z02, res, res)


# Plot a slice across the middle of the grid
slc = range(0, grid_size)
plt.figure()
plt.plot(slc, z0[256, slc])
plt.plot(slc, z02[256, slc])
plt.xlabel('Distance')
plt.ylabel('Elevation')
plt.title('cross-section profile')
plt.show()
'''
