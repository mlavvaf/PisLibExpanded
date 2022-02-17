#!/usr/bin/python

# calculates Sigma and Pi from

# http://inspirehep.net/record/1703751

# using symbolic algebra

# Wed Feb 12 15:14:12 CST 2020 Jeff and Mark
# Tue May 19 14:00:00 CST 2020 Mark - added cylindrical field  variables

# scalarpotential class supplies a variety of useful (self) functions.

# Sigma is the scalar potential

#Pix = \partial_x Sigma
#Piy = \partial_y Sigma
#Piz = \partial_z Sigma

#$\vec{B}(x,y,z) = Pix\hat{x} + Piy\hat{y} + Piz\hat{z}$

#Of note is that the choice of l for any m=0 the field along the z-axis has the shape of x^l

from sympy import assoc_legendre
from sympy import cos,sin,Abs,factorial
from sympy import sqrt
from sympy import trigsimp
from sympy import Derivative
from sympy import lambdify

import sympy as sym

class scalarpotential:
    def __init__(self,ell=2,m=0,Glm=1):

        self.ell=ell
        self.m=m
        self.Glm = Glm
        
        r=sym.Symbol('r',real=True,positive=True)
        theta=sym.Symbol('theta',real=True,positive=True)
        phi=sym.Symbol('phi',real=True,positive=True)

        def c(ell,m):
            if(m>=0):
                return (factorial(ell-1)*(-2)**Abs(m))/factorial(ell+Abs(m))*cos(m*phi)
            else:
                return (factorial(ell-1)*(-2)**Abs(m))/factorial(ell+Abs(m))*sin(Abs(m)*phi)

        legendre_ell=ell+1
        
        Sigma=c(legendre_ell,m)*r**legendre_ell*assoc_legendre(legendre_ell,Abs(m),cos(theta))
        Sigma=sym.simplify(Sigma)
        
        
        Sigma=Sigma.subs({Abs(sin(theta)):sin(theta)})
        Sigma=Sigma.expand(trig=True)
        self.Sigma_spherical=Sigma
        
        #coordinate conversions for cartesian
        x=sym.Symbol('x')
        y=sym.Symbol('y')
        z=sym.Symbol('z')
        
        Sigma=Sigma.subs({r*cos(theta):z,
                          r*sin(theta)*sin(phi):y,
                          r*sin(theta)*cos(phi):x})
        Sigma=Sigma.subs({r**2*sin(theta)**2:x**2+y**2})
        Sigma=Sigma.subs({r**2:x**2+y**2+z**2})
        Sigma=sym.simplify(Sigma.expand())
        self.Sigma=Sigma
        # G = sym.Symbol('G')
        # self.Sigma=G*self.Sigma
        # self.Sigma=self.Sigma.subs({G:Glm})
        
        #cartesian magnetic field calculation
        Pix=Derivative(Sigma,x)
        Piy=Derivative(Sigma,y)
        Piz=Derivative(Sigma,z)

        Pix=Pix.doit()
        Piy=Piy.doit()
        Piz=Piz.doit()

        Pix=sym.simplify(Pix.expand())
        Piy=sym.simplify(Piy.expand())
        Piz=sym.simplify(Piz.expand())
        
        Pix = Pix
        Piy = Piy
        Piz = Piz
        
        self.Pix=Pix
        self.Piy=Piy
        self.Piz=Piz
        
        #lambdify changes the formula to a python lambda function for faster execution
        self.fPix=lambdify([x,y,z],Pix,"numpy")
        self.fPiy=lambdify([x,y,z],Piy,"numpy")
        self.fPiz=lambdify([x,y,z],Piz,"numpy")

        #cylindrical coordinates conversions
        rho=sym.Symbol('rho')
        phi=sym.Symbol('phi')
        z=sym.Symbol('z')

        cyPiRho = cos(phi)*Pix.subs({x:rho*cos(phi),y:rho*sin(phi)}) +sin(phi)*Piy.subs({x:rho*cos(phi),y:rho*sin(phi)})
        cyPiPhi = (-sin(phi)*Pix.subs({x:rho*cos(phi),y:rho*sin(phi)}) +cos(phi)*Piy.subs({x:rho*cos(phi),y:rho*sin(phi)}))
        cyPiZ = Piz.subs({x:rho*cos(phi),y:rho*sin(phi)})
        
        cyPiRho=sym.simplify(cyPiRho.expand(trig=True))
        cyPiPhi=sym.simplify(cyPiPhi.expand(trig=True))
        cyPiZ=sym.simplify(cyPiZ.expand(trig=True))

        self.cyPiRho=cyPiRho
        self.cyPiPhi=cyPiPhi
        self.cyPiZ=cyPiZ
        
        #lambdify changes the formula to a python lambda function for faster execution
        self.fcyPiRho=lambdify([rho,phi,z],cyPiRho,"numpy")
        self.cyPiPhi=lambdify([rho,phi,z],cyPiPhi,"numpy")
        self.fcyPiZ=lambdify([rho,phi,z],cyPiZ,"numpy")
        
        #uncomment these line to see equations printed on creation
        # print('Sigma(%i,%i)='%(ell,m) , self.Sigma)
                
        # print('Sigma_spherical(%i,%i)='%(ell,m) , self.Sigma_spherical)

        # print("Pix = " , Pix)
        # print("Piy = " , Piy)
        # print("Piz = " , Piz)
        
        # print("cyPiRho = " , cyPiRho)
        # print("cyPiPhi = " , cyPiPhi)
        # print("cyPiZ = " , cyPiZ)