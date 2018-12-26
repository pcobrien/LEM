import numpy as np
import matplotlib.pyplot as plt
import gc

def pi_group_scale(d_i, v_i, rho_i=2500.0, target='moon'):
	a_i = d_i/2.0
	# a_i = impactor radius in meters
	# v_i = impactor velocity in m/s
	# rho_i = impactor density in kg/m^3
	# target = target body

	#a_i = np.array(a_i)
	#v_i = np.array(v_i)

	if target == 'moon':
		K_1 = 0.2
		mu = 0.55
		Y_bar = 7.6*(10**6)		# Pa
		rho_t = 2250.0	# kg/m^3
		g = 1.668		# m/s^2
		D_tr = 1.350*9.81/g
		b = D_tr**(-0.18)
		D_crit = ((1.18/b)**(1.0/0.18))*1000.0
	elif target == 'ceres':
		K_1 = 0.2
		mu = 0.55
		Y_bar = 0.1*(10**6)		# Pa
		rho_t = 2250.0	# kg/m^3
		g = 0.27
		D_tr = 1.350*9.81/g
		b = D_tr**(-0.18)
		D_crit = ((1.18/b)**(1.0/0.18))*1000.0
	
	m_i = (4.0/3.0)*np.pi*rho_i*(a_i**3.0)		# m_i = impactor mass in kg 
	
	pi_1 = (m_i/rho_t)
	pi_2 = (g*a_i)/(v_i**2)
	pi_3 = (rho_t/rho_i)**(-1.0/3.0)
	pi_4 = ((Y_bar)/(rho_t*(v_i**2.0)))**((2.0+mu)/2.0)
	exp = -(3.0*mu)/(2.0+mu)
		
	V = K_1 * pi_1 * ( (pi_2*pi_3 + pi_4)**exp)
	D_t = ((24.0*V)/np.pi)**(1.0/3.0)
	D_f = 1.18*D_t
	'''
	if D_t > D_crit:
		#print 'Complex crater'
		D_f = b*(D_t**(1.18))
	else:
		#print 'Simple crater'
		D_f = 1.18*D_t
	'''
	# Return final crater diameter in m
	return D_f
