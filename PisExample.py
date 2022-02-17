#!/usr/bin/python

# this class just generates a quick text output from the Pislib.py library of a scalar potential and magnetic fields.
#
# accidently expanded this to include a fitting set that should be split off??? M. Mcmcrea 2021_02_04
#
# created by M. McCrea 2020_04_01


import matplotlib.pyplot as plt
#plt.switch_backend('agg')#can be useful in Spyder for different plot display if an error message occurs on plt.show()

import time

from Pislib import *
from PislibHelperFunctions import *


#Setting G_lm combinations
l = np.array([0,1,1,3,3,3,4,5])
m = np.array([0,0,0,2,1,3,0,0])
Glm=np.array([1,2,1,1,2,3,4,5])

#Creating arbitrary function combination from Pislib.py
fx ,fy ,fz = CreateArbGlm(l, m, Glm)
print("type(fx) = ", type(fx))
print("fx=" , fx)
print("fy=" , fy)
print("fz=" , fz)

#defining variables for sympy symbolic evaluation
x=sym.Symbol('x')
y=sym.Symbol('y')
z=sym.Symbol('z')

print("fx(1,0,0) = " , fx.evalf(subs={x: 1, y:0 , z:0}))
print("fy(0,1,0) = " , fy.evalf(subs={x: 0, y:1 , z:0}))
print("fz(0,0,1) = " , fz.evalf(subs={x: 0, y:0 , z:1}))

print("fx(2,0,0) = " , fx.evalf(subs={x: 2, y:0 , z:0}))
print("fy(0,2,0) = " , fy.evalf(subs={x: 0, y:2 , z:0}))
print("fz(0,0,2) = " , fz.evalf(subs={x: 0, y:0 , z:2}))

#creating lambda functions for faster execution if further symbolic manipulation is not required after function creation.
lfx ,lfy ,lfz = CreateArbGlmLambdafied(l, m, Glm)

print("It can be seen that lambdified functiosn evaluate to identical values, but different formating due to being a different type of object.")
print("lfx(1,0,0) = " , lfx(1,0,0))
print("lfy(0,1,0) = " , lfy(0,1,0))
print("lfz(0,0,1) = " , lfz(0,0,1))

print("lfx(2,0,0) = " , lfx(2,0,0))
print("lfy(0,2,0) = " , lfy(0,2,0))
print("lfz(0,0,2) = " , lfz(0,0,2))

print("The types are given below as:")
print("type(fx) = ", type(fx))
print("type(lfx) = ", type(lfx))


#this works
print("lfz(0,0,2)+lfy(0,0,2) = " , lfz(0,0,2)+lfy(0,0,2))
#this works
t = fz+fy
print("t = " , t)
#this gives an error message, lambdified equations cannot be added symbolically.
#print("lfz+lfy = " , lfz+lfy)

print("\n\nlambdified functions are useful as they execute faster.  See below for a comparison showing the lambdified functions are roughly 100 times faster than the symbolic functions in evaluation")
time_start = time.time()
for i in range(500):
  fx.evalf(subs={x: i, y:0 , z:0})
time_end = time.time()
print("total execution time of symbolic functions:" , time_end-time_start, " seconds")

time_start = time.time()
for i in range(500):
  lfx(i,0,0)
time_end = time.time()
print("total execution time of lambdified functions:" , time_end-time_start, " seconds")


print("\n\n This next line of code is a plotting example of the Glm functions from HelperFunctions.  Look for f00.pdf as the output in the same folder for how it works.")
PlotArbGlm(l,m,Glm)

print("\n\n The lambdified functions will operate over numpy arrays while the symbolic function will not.")
value = 10
xx = np.array(range(value))
yy = xx
zz = xx
print("xx = \n" , xx)
print("yy = \n" , yy)
print("zz = \n" , zz)

print("lfx(xx,yy,zz) = \n" , lfx(xx,yy,zz))

#if uncommented this will give an error message.
#print("symbolically: fx(xx,yy,zz) = \n" , fx.evalf(subs={x: xx, y:yy , z:zz}))

print("An issue with using functions in python is that if the function definition does not include a variable then only a scalar will be returned.  This becomes an issue when using l=0, m=0 for the x or y components for example.")
#Setting G_lm combinations
l = np.array([0])
m = np.array([0])
Glm=np.array([1])

fx0 ,fy0 ,fz0 = CreateArbGlm(l, m, Glm)
print("type(fx) = ", type(fx))
print("fx0=" , fx0)
print("fy0=" , fy0)
print("fz0=" , fz0)
lfx0 ,lfy0 ,lfz0 = CreateArbGlmLambdafied(l, m, Glm)

print("lfx0(xx,yy,zz) = \n" , lfx0(xx,yy,zz))
print("lfz0(xx,yy,zz) = \n" , lfz0(xx,yy,zz))



print("\n\nAnd for an example of creating a combination up to an arbitrary order l for all m:")
Bx,By,Bz = CreateGlmToOrder(1,[-0.2,0.5,-0.8,-2,-3,0.25,4,6])
print("Bx = " , Bx)
print("By = " , By)
print("Bz = " , Bz)