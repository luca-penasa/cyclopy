# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 23:57:30 2013

@author: luca
"""
from __future__ import division
import numpy as np

def RotateOnX(point, angle): #note is a clockwise rotation
    angle_rad = angle * np.pi /180
    #get the point in a 4d coord system
    point4 = np.matrix(np.concatenate( [point, [1]] )).T
    
    RotX = np.eye(4)

    RotX[1,1] = RotX[2,2] = np.cos(angle_rad)
    RotX[1,2] = np.sin(angle_rad)    
    RotX[2,1] = - np.sin(angle_rad)
    
    
    rot_point = RotX * point4
    
    return np.array(rot_point).T[0][0:3]
    

def RotateOnZ(point, angle): #note is a clockwise rotation
    angle_rad = angle * np.pi /180
    #get the point in a 4d coord system
    point4 = np.matrix(np.concatenate( [point, [1]] )).T
    
    RotX = np.eye(4)

    RotX[0,0] = RotX[1,1] = np.cos(angle_rad)
    RotX[0,1] = np.sin(angle_rad)    
    RotX[1,0] = - np.sin(angle_rad)
    
    
    rot_point = RotX * point4
    
    return np.array(rot_point).T[0][0:3]
    
def DipDirectionAndDipAngleToDipVector(dipdirection, dipangle):
    v = np.array([0,1,0]) #unit vector oriented with y   
    v = RotateOnX(v, dipangle) #now the vector have the same dip angle but is oriented the wring direction
    v = RotateOnZ(v, dipdirection) #now it also have right direction
    return v

def DipDirectionToStrikeVector(dipdirection):
    """
    strike may be also oriented the other way, change the sign if you want
    """
    strike_angle = dipdirection - 90
    v = np.array([0,1,0]) #unit vector oriented with y   
    
    return RotateOnZ(v, strike_angle)    
    
def DipDirectionAndDipAngleToPlaneNormal(dipdirection, dipangle):
    dip_vector = DipDirectionAndDipAngleToDipVector(dipdirection, dipangle)
    strike_vector = DipDirectionToStrikeVector(dipdirection)
    
    return np.cross(dip_vector, strike_vector)
    
def ProjectOnVector(v1, v2):
    """ Project V1 on V2"""
    v2norm =v2 / np.linalg.norm(v2)
    return v1.dot(v2norm) * v2norm
    
def VectorNorm(v):
    return np.linalg.norm(v)
    
def Project3DPointsOnStratigraphicNormal(points, normal):
    values = []    
    for point in points:
        projected = ProjectOnVector(point, normal)
        value = VectorNorm(projected)
        values.append(value)
        
    return np.array(values)

    
    
    