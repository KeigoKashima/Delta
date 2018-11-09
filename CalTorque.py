import numpy as np
from sympy import *
import sympy
import math

##Calculate torques for each mortors by given forces


##-------------definition of variables-----------------------------------------

var('E1 F1 G1 E2 F2 G2 E3 F3 G3')
var('x y z l L a')

##-----------------------function-----------------------------------------------

def ti(*arg):
	"""
	 	Put arguments into the equation of ti

			>>>ti(E1,F1,G1)
			{-F1 + sqrt(E1**2 + F1**2 - G1**2)}/(G1- E1)
	"""
	return ((-arg[1] + sqrt(arg[0]**2 + arg[1]**2 - arg[2]**2))/(arg[2] - arg[0]))


##-----------Definiton of Ei,Fi,Gi------------------------------------------
E1 = -2*L*(x + np.sqrt(3)*(y + 2*a))
F1 = -4*L*z
G1 = x**2 + y**2 + z**2 + 4*a*y +4*a**2 + 4*L**2 - 4*l**2

E2 = 4*L*(x - np.sqrt(3)*a)
F2 = -4*L*z
G2 = x**2 + y**2 + z**2 - 2*np.sqrt(3)*a*x - 2*a*y + 4*a**2 + 4*L**2 - 4*l**2

E3 = -2*L*(x - np.sqrt(3)*(y + 2*a))
F3 = -4*L*z
G3 = x**2 + y**2 + z**2 + 2*np.sqrt(3)*a*x - 2*a*y + 4*a**2 + 4*L**2 - 4*l**2

##------------Calculation----------------------------------------------------
t1  = ti(E1,F1,G1)
t2  = ti(E2,F2,G2)
t3  = ti(E3,F3,G3)

## d/dx (2*arctan(ti)) = 2/(1+ti^2) * dti/dx
f1x = diff(t1, x) * (2/(1+t1**2))
f1y = diff(t1, y) * (2/(1+t1**2))
f1z = diff(t1, z) * (2/(1+t1**2))
f2x = diff(t2, x) * (2/(1+t2**2))
f2y = diff(t2, y) * (2/(1+t2**2))
f2z = diff(t2, z) * (2/(1+t2**2))
f3x = diff(t3, x) * (2/(1+t3**2))
f3y = diff(t3, y) * (2/(1+t3**2))
f3z = diff(t3, z) * (2/(1+t3**2))

J = Matrix([
[f1x,f2x,f3x],
[f1y,f2y,f3y],
[f1z,f2z,f3z]
])


##-----------Substitution-----------

s_x = 0
s_y = 0
s_z = 0
s_l	= 0.20/2 			# 2L > l
s_L = 0.20/2 			# 2L + 2l = 0.4
s_a = 0.1					#np.sqrt(3)*s_l - np.sqrt(4*s_L**2 - s_l**2)
F = Matrix([[0,0,-10.0]]).T

print('[2l,2L,a]=',2*s_l,2*s_L,s_a)
print('F=',F[0],F[1],F[2])
print('z','x','tau[0]','tau[1]','tau[2]')
while s_z < 0.4:
	s_y = 0
	while s_y < s_a:		##Until the hand reach the edge

		J1 = J.subs([(x,s_x),(y,s_y),(z,s_z),(l,s_l),(L,s_L),(a,s_a)])
		J2 = J1.inv()		##Inverse
		tau = J2*F

		##----------Output---------
		if (type(tau[0])==Float)and(type(tau[1])==Float)and(type(tau[2])==Float):
			if s_y > s_a - 0.02:
				print(s_z,s_y,tau[0],tau[1],tau[2])
		else:
			J1 = J.subs([(x,s_x),(y,s_y-0.01),(z,s_z),(l,s_l),(L,s_L),(a,s_a)])
			J2 = J1.inv()		##Invers
			tau = J2*F
			if(s_y == 0):
				print(s_z,s_y,tau[0],tau[1],tau[2])
			else:
				print(s_z,s_y-0.01,tau[0],tau[1],tau[2])
			break

		s_y = s_y + 0.01
	s_z = s_z + 0.01
