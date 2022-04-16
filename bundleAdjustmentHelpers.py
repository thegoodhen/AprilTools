from __future__ import print_function
import math
import numpy as np

import urllib
import urllib.request
import bz2
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
import time
from scipy.optimize import least_squares
# This is a sample Python script.


def prepareBundleAdjustmentParamsForSingleFrame(_pointCoords3D, _actualPointCoords2D,_rvec, _tvec, _cameraMatrix, _distCoeffs):
    '''

    :param _pointCoords3D: Individual 3D points, as a numpy array n x 3
    :param _actualPointCoords2D: Observed corresponding 2D points, as a numpy array n x 2
    :param _rvec: A Rodriguez rotation vector, determining the rotation of the camera, as returned by solvePnP
    :param _tvec: A translation vector, determining the position of the camera, as returned by solvePnP
    :param _cameraMatrix: The camera matrix of intrinsic parameters
    :param _distCoeffs: Distortion coefficients of the camera
    :return:
    '''
    _actualPointCoords2D=(_actualPointCoords2D).tolist()
    wx=_rvec[0][0]
    wy=_rvec[1][0]
    wz=_rvec[2][0]
    tx=_tvec[0][0]
    ty=_tvec[1][0]
    tz=_tvec[2][0]

    print("KVAK!")
    theta=math.sqrt((wx*wx)+(wy*wy)+(wz*wz))
    print(theta)
    omega=np.array([[0, -wz, wy], [wz, 0, -wx], [-wy, wx, 0]])
    print(omega)
    R=np.eye(3)+(np.sin(theta)/theta)*omega+((1-np.cos(theta))/(theta*theta))*(omega@omega)
    print(R)
    tmat=np.hstack((R,_tvec))


    #TODO: Also include the distortion in the calculations of the jacobians
    #https://learnopencv.com/understanding-lens-distortion/

    f=_cameraMatrix[0][0]#TODO: there are 2 focal lengths for x and y (non-square pixels!) Make sure to handle that!
    u0=_cameraMatrix[1][2]
    v0=_cameraMatrix[2][2]

    s1 = wx ** 2 + wy ** 2 + wz ** 2
    s2 = math.sin((s1) ** (1 / 2))
    s3 = math.cos((s1) ** (1 / 2))

    A_list=[]
    B_list=[]
    projectedPointCoords2D_list=[]


    for it in range(_pointCoords3D.shape[0]):
        print(it)
        X = _pointCoords3D[it, 0]
        Y = _pointCoords3D[it, 1]
        Z = _pointCoords3D[it, 2]

        uvs = _cameraMatrix @ tmat @ [[X], [Y], [Z], [1]]
        _projectedPointCoords2D = [0, 0]
        _projectedPointCoords2D[0] = uvs[0][0] / uvs[2][0]
        _projectedPointCoords2D[1] = uvs[1][0]/ uvs[2][0]
        projectedPointCoords2D_list.append(_projectedPointCoords2D)



        Jx_wx = -(X * (f * ( (2 * wx * (s3 - 1) * (wy ** 2 + wz ** 2)) / (s1) ** 2 + (wx * s2 * (wy ** 2 + wz ** 2)) / (s1) ** ( 3 / 2)) - u0 * ( (wx * wy * s2) / (s1) ** (3 / 2) - (wx * wy * s3) / (s1) - (wz * (s3 - 1)) / (s1) + ( wx ** 2 * wz * s2) / (s1) ** (3 / 2) + (2 * wx ** 2 * wz * (s3 - 1)) / ( s1) ** 2)) - Y * (f * ( (wx * wz * s2) / (s1) ** (3 / 2) - (wx * wz * s3) / (s1) - (wy * (s3 - 1)) / (s1) + ( wx ** 2 * wy * s2) / (s1) ** (3 / 2) + (2 * wx ** 2 * wy * (s3 - 1)) / (s1) ** 2) + u0 * ( s2 / (s1) ** (1 / 2) + (wx ** 2 * s3) / (s1) - ( wx ** 2 * s2) / (s1) ** (3 / 2) + ( 2 * wx * wy * wz * (s3 - 1)) / ( s1) ** 2 + (wx * wy * wz * s2) / (s1) ** ( 3 / 2))) + Z * (u0 * ( (2 * wx * (s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) ** 2 - (2 * wx * (s3 - 1)) / (s1) + ( wx * s2 * (wx ** 2 + wy ** 2)) / (s1) ** (3 / 2)) - f * ((wx * wy * s3) / (s1) - (wz * (s3 - 1)) / ( s1) - (wx * wy * s2) / (s1) ** (3 / 2) + (wx ** 2 * wz * s2) / (s1) ** (3 / 2) + (2 * wx ** 2 * wz * ( s3 - 1)) / (s1) ** 2))) / (tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) - ((X * ( (wx * wy * s2) / (s1) ** (3 / 2) - (wx * wy * s3) / (s1) - (wz * (s3 - 1)) / (s1) + ( wx ** 2 * wz * s2) / (s1) ** (3 / 2) + (2 * wx ** 2 * wz * (s3 - 1)) / (s1) ** 2) - Z * ((2 * wx * ( s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) ** 2 - (2 * wx * (s3 - 1)) / (s1) + (wx * s2 * ( wx ** 2 + wy ** 2)) / (s1) ** (3 / 2)) + Y * (s2 / (s1) ** (1 / 2) + (wx ** 2 * s3) / (s1) - ( wx ** 2 * s2) / (s1) ** (3 / 2) + (2 * wx * wy * wz * (s3 - 1)) / (s1) ** 2 + (wx * wy * wz * s2) / ( s1) ** (3 / 2))) * (f * tx - Y * ( f * ((wz * s2) / (s1) ** (1 / 2) + (wx * wy * (s3 - 1)) / (s1)) - u0 * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1))) + tz * u0 - X * (u0 * ( (wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) - f * (((s3 - 1) * (wy ** 2 + wz ** 2)) / ( s1) + 1)) + Z * (f * ((wy * s2) / (s1) ** (1 / 2) - (wx * wz * (s3 - 1)) / (s1)) + u0 * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) ** 2

        Jx_wy = -(X * (f * ((2 * wy * (s3 - 1) * (wy ** 2 + wz ** 2)) / (s1) ** 2 - (2 * wy * (s3 - 1)) / (s1) + ( wy * s2 * (wy ** 2 + wz ** 2)) / (s1) ** (3 / 2)) - u0 * ( (wy ** 2 * s2) / (s1) ** (3 / 2) - (wy ** 2 * s3) / (s1) - s2 / (s1) ** (1 / 2) + ( 2 * wx * wy * wz * (s3 - 1)) / (s1) ** 2 + (wx * wy * wz * s2) / (s1) ** ( 3 / 2))) + Z * (u0 * ( (2 * wy * (s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) ** 2 - (2 * wy * (s3 - 1)) / (s1) + ( wy * s2 * (wx ** 2 + wy ** 2)) / (s1) ** (3 / 2)) - f * ( s2 / (s1) ** (1 / 2) + (wy ** 2 * s3) / (s1) - ( wy ** 2 * s2) / (s1) ** (3 / 2) + ( 2 * wx * wy * wz * (s3 - 1)) / ( s1) ** 2 + (wx * wy * wz * s2) / (s1) ** ( 3 / 2))) - Y * (f * ( (wy * wz * s2) / (s1) ** (3 / 2) - (wy * wz * s3) / (s1) - (wx * (s3 - 1)) / (s1) + ( wx * wy ** 2 * s2) / (s1) ** (3 / 2) + (2 * wx * wy ** 2 * (s3 - 1)) / (s1) ** 2) + u0 * (( wx * wy * s3) / ( s1) - ( wz * ( s3 - 1)) / ( s1) - ( wx * wy * s2) / ( s1) ** ( 3 / 2) + ( wy ** 2 * wz * s2) / ( s1) ** ( 3 / 2) + ( 2 * wy ** 2 * wz * ( s3 - 1)) / ( s1) ** 2))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) - ((Y * ( (wx * wy * s3) / (s1) - (wz * (s3 - 1)) / (s1) - (wx * wy * s2) / (s1) ** (3 / 2) + ( wy ** 2 * wz * s2) / (s1) ** (3 / 2) + (2 * wy ** 2 * wz * (s3 - 1)) / (s1) ** 2) - Z * ((2 * wy * ( s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) ** 2 - (2 * wy * (s3 - 1)) / (s1) + (wy * s2 * ( wx ** 2 + wy ** 2)) / (s1) ** (3 / 2)) + X * ((wy ** 2 * s2) / (s1) ** (3 / 2) - (wy ** 2 * s3) / ( s1) - s2 / (s1) ** (1 / 2) + (2 * wx * wy * wz * (s3 - 1)) / (s1) ** 2 + (wx * wy * wz * s2) / (s1) ** ( 3 / 2))) * (f * tx - Y * ( f * ((wz * s2) / (s1) ** (1 / 2) + (wx * wy * (s3 - 1)) / (s1)) - u0 * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1))) + tz * u0 - X * (u0 * ( (wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) - f * (((s3 - 1) * (wy ** 2 + wz ** 2)) / ( s1) + 1)) + Z * (f * ((wy * s2) / (s1) ** (1 / 2) - (wx * wz * (s3 - 1)) / (s1)) + u0 * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) ** 2

        Jx_wz = -(Z * (u0 * ( (2 * wz * (s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) ** 2 + (wz * s2 * (wx ** 2 + wy ** 2)) / (s1) ** ( 3 / 2)) - f * ((wy * wz * s3) / (s1) - (wx * (s3 - 1)) / (s1) - (wy * wz * s2) / (s1) ** (3 / 2) + ( wx * wz ** 2 * s2) / (s1) ** (3 / 2) + (2 * wx * wz ** 2 * (s3 - 1)) / (s1) ** 2)) - Y * (u0 * ( (wx * wz * s3) / (s1) - (wy * (s3 - 1)) / (s1) - (wx * wz * s2) / (s1) ** (3 / 2) + ( wy * wz ** 2 * s2) / (s1) ** (3 / 2) + (2 * wy * wz ** 2 * (s3 - 1)) / (s1) ** 2) + f * (( wz ** 2 * s2) / ( s1) ** ( 3 / 2) - ( wz ** 2 * s3) / ( s1) - s2 / ( s1) ** ( 1 / 2) + ( 2 * wx * wy * wz * ( s3 - 1)) / ( s1) ** 2 + ( wx * wy * wz * s2) / ( s1) ** ( 3 / 2))) + X * ( f * ((2 * wz * (s3 - 1) * (wy ** 2 + wz ** 2)) / (s1) ** 2 - (2 * wz * (s3 - 1)) / (s1) + ( wz * s2 * (wy ** 2 + wz ** 2)) / (s1) ** (3 / 2)) - u0 * ( (wy * wz * s2) / (s1) ** (3 / 2) - (wy * wz * s3) / (s1) - (wx * (s3 - 1)) / ( s1) + (wx * wz ** 2 * s2) / (s1) ** (3 / 2) + (2 * wx * wz ** 2 * (s3 - 1)) / ( s1) ** 2))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) - ((X * ( (wy * wz * s2) / (s1) ** (3 / 2) - (wy * wz * s3) / (s1) - (wx * (s3 - 1)) / (s1) + ( wx * wz ** 2 * s2) / (s1) ** (3 / 2) + (2 * wx * wz ** 2 * (s3 - 1)) / (s1) ** 2) - Z * ((2 * wz * ( s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) ** 2 + (wz * s2 * (wx ** 2 + wy ** 2)) / (s1) ** (3 / 2)) + Y * ( (wx * wz * s3) / ( s1) - (wy * (s3 - 1)) / ( s1) - ( wx * wz * s2) / ( s1) ** (3 / 2) + ( wy * wz ** 2 * s2) / ( s1) ** (3 / 2) + ( 2 * wy * wz ** 2 * ( s3 - 1)) / ( s1) ** 2)) * ( f * tx - Y * (f * ( (wz * s2) / ( s1) ** (1 / 2) + ( wx * wy * ( s3 - 1)) / ( s1)) - u0 * (( wx * s2) / ( s1) ** ( 1 / 2) - ( wy * wz * ( s3 - 1)) / ( s1))) + tz * u0 - X * ( u0 * (( wy * s2) / ( s1) ** ( 1 / 2) + ( wx * wz * ( s3 - 1)) / ( s1)) - f * ( ( ( s3 - 1) * ( wy ** 2 + wz ** 2)) / ( s1) + 1)) + Z * ( f * (( wy * s2) / ( s1) ** ( 1 / 2) - ( wx * wz * ( s3 - 1)) / ( s1)) + u0 * ( ( ( s3 - 1) * ( wx ** 2 + wy ** 2)) / ( s1) + 1)))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) ** 2

        Jx_tx = f / (tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1))

        Jx_ty = 0

        Jx_tz = u0 / (tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) - (f * tx - Y * ( f * ((wz * s2) / (s1) ** (1 / 2) + (wx * wy * (s3 - 1)) / (s1)) - u0 * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1))) + tz * u0 - X * (u0 * ( (wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) - f * (((s3 - 1) * (wy ** 2 + wz ** 2)) / ( s1) + 1)) + Z * (f * ((wy * s2) / (s1) ** (1 / 2) - (wx * wz * (s3 - 1)) / (s1)) + u0 * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) ** 2

        Jx_X = (((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) * (f * tx - Y * ( f * ((wz * s2) / (s1) ** (1 / 2) + (wx * wy * (s3 - 1)) / (s1)) - u0 * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1))) + tz * u0 - X * (u0 * ( (wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) - f * (((s3 - 1) * (wy ** 2 + wz ** 2)) / ( s1) + 1)) + Z * (f * ((wy * s2) / (s1) ** (1 / 2) - (wx * wz * (s3 - 1)) / (s1)) + u0 * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) ** 2 - ( u0 * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) - f * ( ((s3 - 1) * (wy ** 2 + wz ** 2)) / (s1) + 1)) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1))

        Jx_Y = -(f * ((wz * s2) / (s1) ** (1 / 2) + (wx * wy * (s3 - 1)) / (s1)) - u0 * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) - ( ((wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) * (f * tx - Y * ( f * ((wz * s2) / (s1) ** (1 / 2) + (wx * wy * (s3 - 1)) / (s1)) - u0 * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1))) + tz * u0 - X * (u0 * ( (wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) - f * (((s3 - 1) * ( wy ** 2 + wz ** 2)) / (s1) + 1)) + Z * (f * ( (wy * s2) / (s1) ** (1 / 2) - (wx * wz * (s3 - 1)) / (s1)) + u0 * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / ( s1) + 1)))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) ** 2

        Jx_Z = (f * ((wy * s2) / (s1) ** (1 / 2) - (wx * wz * (s3 - 1)) / (s1)) + u0 * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) - ( (((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1) * (f * tx - Y * ( f * ((wz * s2) / (s1) ** (1 / 2) + (wx * wy * (s3 - 1)) / (s1)) - u0 * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1))) + tz * u0 - X * (u0 * ( (wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) - f * (((s3 - 1) * ( wy ** 2 + wz ** 2)) / (s1) + 1)) + Z * (f * ( (wy * s2) / (s1) ** (1 / 2) - (wx * wz * (s3 - 1)) / (s1)) + u0 * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / ( s1) + 1)))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) ** 2

        Jy_wx = -(Y * (f * ((2 * wx * (s3 - 1) * (wx ** 2 + wz ** 2)) / (s1) ** 2 - (2 * wx * (s3 - 1)) / (s1) + ( wx * s2 * (wx ** 2 + wz ** 2)) / (s1) ** (3 / 2)) - v0 * ( s2 / (s1) ** (1 / 2) + (wx ** 2 * s3) / (s1) - (wx ** 2 * s2) / (s1) ** (3 / 2) + ( 2 * wx * wy * wz * (s3 - 1)) / (s1) ** 2 + (wx * wy * wz * s2) / (s1) ** ( 3 / 2))) + Z * (v0 * ( (2 * wx * (s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) ** 2 - (2 * wx * (s3 - 1)) / (s1) + ( wx * s2 * (wx ** 2 + wy ** 2)) / (s1) ** (3 / 2)) - f * ((wx ** 2 * s2) / (s1) ** (3 / 2) - ( wx ** 2 * s3) / (s1) - s2 / (s1) ** (1 / 2) + (2 * wx * wy * wz * (s3 - 1)) / (s1) ** 2 + ( wx * wy * wz * s2) / (s1) ** ( 3 / 2))) - X * (f * ( (wx * wz * s3) / (s1) - (wy * (s3 - 1)) / (s1) - (wx * wz * s2) / (s1) ** (3 / 2) + ( wx ** 2 * wy * s2) / (s1) ** (3 / 2) + (2 * wx ** 2 * wy * (s3 - 1)) / (s1) ** 2) + v0 * (( wx * wy * s2) / ( s1) ** ( 3 / 2) - ( wx * wy * s3) / ( s1) - ( wz * ( s3 - 1)) / ( s1) + ( wx ** 2 * wz * s2) / ( s1) ** ( 3 / 2) + ( 2 * wx ** 2 * wz * ( s3 - 1)) / ( s1) ** 2))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) - ((X * ( (wx * wy * s2) / (s1) ** (3 / 2) - (wx * wy * s3) / (s1) - (wz * (s3 - 1)) / (s1) + ( wx ** 2 * wz * s2) / (s1) ** (3 / 2) + (2 * wx ** 2 * wz * (s3 - 1)) / (s1) ** 2) - Z * ((2 * wx * ( s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) ** 2 - (2 * wx * (s3 - 1)) / (s1) + (wx * s2 * ( wx ** 2 + wy ** 2)) / (s1) ** (3 / 2)) + Y * (s2 / (s1) ** (1 / 2) + (wx ** 2 * s3) / (s1) - ( wx ** 2 * s2) / (s1) ** (3 / 2) + (2 * wx * wy * wz * (s3 - 1)) / (s1) ** 2 + (wx * wy * wz * s2) / ( s1) ** (3 / 2))) * (f * ty + X * ( f * ((wz * s2) / (s1) ** (1 / 2) - (wx * wy * (s3 - 1)) / (s1)) - v0 * ( (wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1))) + tz * v0 + Y * (v0 * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + f * (((s3 - 1) * (wx ** 2 + wz ** 2)) / ( s1) + 1)) - Z * (f * ((wx * s2) / (s1) ** (1 / 2) + (wy * wz * (s3 - 1)) / (s1)) - v0 * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) ** 2

        Jy_wy = -(Y * (f * ( (2 * wy * (s3 - 1) * (wx ** 2 + wz ** 2)) / (s1) ** 2 + (wy * s2 * (wx ** 2 + wz ** 2)) / (s1) ** ( 3 / 2)) - v0 * ( (wx * wy * s3) / (s1) - (wz * (s3 - 1)) / (s1) - (wx * wy * s2) / (s1) ** (3 / 2) + ( wy ** 2 * wz * s2) / (s1) ** (3 / 2) + (2 * wy ** 2 * wz * (s3 - 1)) / ( s1) ** 2)) - X * (f * ( (wy * wz * s3) / (s1) - (wx * (s3 - 1)) / (s1) - (wy * wz * s2) / (s1) ** (3 / 2) + ( wx * wy ** 2 * s2) / (s1) ** (3 / 2) + (2 * wx * wy ** 2 * (s3 - 1)) / (s1) ** 2) + v0 * ( (wy ** 2 * s2) / (s1) ** (3 / 2) - (wy ** 2 * s3) / ( s1) - s2 / (s1) ** (1 / 2) + ( 2 * wx * wy * wz * (s3 - 1)) / ( s1) ** 2 + (wx * wy * wz * s2) / (s1) ** ( 3 / 2))) + Z * (v0 * ( (2 * wy * (s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) ** 2 - (2 * wy * (s3 - 1)) / (s1) + ( wy * s2 * (wx ** 2 + wy ** 2)) / (s1) ** (3 / 2)) - f * ((wx * wy * s2) / (s1) ** (3 / 2) - ( wx * wy * s3) / (s1) - (wz * (s3 - 1)) / (s1) + (wy ** 2 * wz * s2) / (s1) ** (3 / 2) + ( 2 * wy ** 2 * wz * ( s3 - 1)) / (s1) ** 2))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) - ((Y * ( (wx * wy * s3) / (s1) - (wz * (s3 - 1)) / (s1) - (wx * wy * s2) / (s1) ** (3 / 2) + ( wy ** 2 * wz * s2) / (s1) ** (3 / 2) + (2 * wy ** 2 * wz * (s3 - 1)) / (s1) ** 2) - Z * ((2 * wy * ( s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) ** 2 - (2 * wy * (s3 - 1)) / (s1) + (wy * s2 * ( wx ** 2 + wy ** 2)) / (s1) ** (3 / 2)) + X * ((wy ** 2 * s2) / (s1) ** (3 / 2) - (wy ** 2 * s3) / ( s1) - s2 / (s1) ** (1 / 2) + (2 * wx * wy * wz * (s3 - 1)) / (s1) ** 2 + (wx * wy * wz * s2) / (s1) ** ( 3 / 2))) * (f * ty + X * ( f * ((wz * s2) / (s1) ** (1 / 2) - (wx * wy * (s3 - 1)) / (s1)) - v0 * ( (wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1))) + tz * v0 + Y * (v0 * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + f * (((s3 - 1) * (wx ** 2 + wz ** 2)) / ( s1) + 1)) - Z * (f * ((wx * s2) / (s1) ** (1 / 2) + (wy * wz * (s3 - 1)) / (s1)) - v0 * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) ** 2

        Jy_wz = -(Z * (v0 * ( (2 * wz * (s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) ** 2 + (wz * s2 * (wx ** 2 + wy ** 2)) / (s1) ** ( 3 / 2)) - f * ((wx * wz * s2) / (s1) ** (3 / 2) - (wx * wz * s3) / (s1) - (wy * (s3 - 1)) / (s1) + ( wy * wz ** 2 * s2) / (s1) ** (3 / 2) + (2 * wy * wz ** 2 * (s3 - 1)) / (s1) ** 2)) - X * (v0 * ( (wy * wz * s2) / (s1) ** (3 / 2) - (wy * wz * s3) / (s1) - (wx * (s3 - 1)) / (s1) + ( wx * wz ** 2 * s2) / (s1) ** (3 / 2) + (2 * wx * wz ** 2 * (s3 - 1)) / (s1) ** 2) + f * (s2 / ( s1) ** (1 / 2) + (wz ** 2 * s3) / (s1) - (wz ** 2 * s2) / (s1) ** (3 / 2) + (2 * wx * wy * wz * (s3 - 1)) / ( s1) ** 2 + ( wx * wy * wz * s2) / ( s1) ** ( 3 / 2))) + Y * ( f * ((2 * wz * (s3 - 1) * (wx ** 2 + wz ** 2)) / (s1) ** 2 - (2 * wz * (s3 - 1)) / (s1) + ( wz * s2 * (wx ** 2 + wz ** 2)) / (s1) ** (3 / 2)) - v0 * ( (wx * wz * s3) / (s1) - (wy * (s3 - 1)) / (s1) - (wx * wz * s2) / (s1) ** ( 3 / 2) + (wy * wz ** 2 * s2) / (s1) ** (3 / 2) + ( 2 * wy * wz ** 2 * (s3 - 1)) / (s1) ** 2))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) - ((X * ( (wy * wz * s2) / (s1) ** (3 / 2) - (wy * wz * s3) / (s1) - (wx * (s3 - 1)) / (s1) + ( wx * wz ** 2 * s2) / (s1) ** (3 / 2) + (2 * wx * wz ** 2 * (s3 - 1)) / (s1) ** 2) - Z * ((2 * wz * ( s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) ** 2 + (wz * s2 * (wx ** 2 + wy ** 2)) / (s1) ** (3 / 2)) + Y * ( (wx * wz * s3) / ( s1) - (wy * (s3 - 1)) / ( s1) - ( wx * wz * s2) / ( s1) ** (3 / 2) + ( wy * wz ** 2 * s2) / ( s1) ** (3 / 2) + ( 2 * wy * wz ** 2 * ( s3 - 1)) / ( s1) ** 2)) * ( f * ty + X * (f * ( (wz * s2) / ( s1) ** (1 / 2) - ( wx * wy * ( s3 - 1)) / ( s1)) - v0 * (( wy * s2) / ( s1) ** ( 1 / 2) + ( wx * wz * ( s3 - 1)) / ( s1))) + tz * v0 + Y * ( v0 * (( wx * s2) / ( s1) ** ( 1 / 2) - ( wy * wz * ( s3 - 1)) / ( s1)) + f * ( ( ( s3 - 1) * ( wx ** 2 + wz ** 2)) / ( s1) + 1)) - Z * ( f * (( wx * s2) / ( s1) ** ( 1 / 2) + ( wy * wz * ( s3 - 1)) / ( s1)) - v0 * ( ( ( s3 - 1) * ( wx ** 2 + wy ** 2)) / ( s1) + 1)))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) ** 2

        Jy_tx = 0

        Jy_ty = f / (tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1))

        Jy_tz = v0 / (tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) - (f * ty + X * ( f * ((wz * s2) / (s1) ** (1 / 2) - (wx * wy * (s3 - 1)) / (s1)) - v0 * ( (wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1))) + tz * v0 + Y * (v0 * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + f * (((s3 - 1) * (wx ** 2 + wz ** 2)) / ( s1) + 1)) - Z * (f * ((wx * s2) / (s1) ** (1 / 2) + (wy * wz * (s3 - 1)) / (s1)) - v0 * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) ** 2

        Jy_X = (f * ((wz * s2) / (s1) ** (1 / 2) - (wx * wy * (s3 - 1)) / (s1)) - v0 * ( (wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) + ( ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) * (f * ty + X * ( f * ((wz * s2) / (s1) ** (1 / 2) - (wx * wy * (s3 - 1)) / (s1)) - v0 * ( (wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1))) + tz * v0 + Y * (v0 * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + f * (((s3 - 1) * ( wx ** 2 + wz ** 2)) / (s1) + 1)) - Z * (f * ( (wx * s2) / (s1) ** (1 / 2) + (wy * wz * (s3 - 1)) / (s1)) - v0 * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / ( s1) + 1)))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) ** 2

        Jy_Y = (v0 * ((wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + f * ( ((s3 - 1) * (wx ** 2 + wz ** 2)) / (s1) + 1)) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) - ( ((wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) * (f * ty + X * ( f * ((wz * s2) / (s1) ** (1 / 2) - (wx * wy * (s3 - 1)) / (s1)) - v0 * ( (wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1))) + tz * v0 + Y * (v0 * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + f * (((s3 - 1) * ( wx ** 2 + wz ** 2)) / (s1) + 1)) - Z * (f * ( (wx * s2) / (s1) ** (1 / 2) + (wy * wz * (s3 - 1)) / (s1)) - v0 * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / ( s1) + 1)))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) ** 2

        Jy_Z = -(f * ((wx * s2) / (s1) ** (1 / 2) + (wy * wz * (s3 - 1)) / (s1)) - v0 * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) - ( (((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1) * (f * ty + X * ( f * ((wz * s2) / (s1) ** (1 / 2) - (wx * wy * (s3 - 1)) / (s1)) - v0 * ( (wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1))) + tz * v0 + Y * (v0 * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + f * (((s3 - 1) * ( wx ** 2 + wz ** 2)) / (s1) + 1)) - Z * (f * ( (wx * s2) / (s1) ** (1 / 2) + (wy * wz * (s3 - 1)) / (s1)) - v0 * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / ( s1) + 1)))) / ( tz - X * ((wy * s2) / (s1) ** (1 / 2) + (wx * wz * (s3 - 1)) / (s1)) + Y * ( (wx * s2) / (s1) ** (1 / 2) - (wy * wz * (s3 - 1)) / (s1)) + Z * ( ((s3 - 1) * (wx ** 2 + wy ** 2)) / (s1) + 1)) ** 2

        A=([[Jx_wx, Jx_wy, Jx_wz, Jx_tx, Jx_ty, Jx_tz],[Jy_wx,Jy_wy,Jy_wz, Jy_tx, Jy_ty, Jy_tz]])
        B=([[Jx_X, Jx_Y, Jx_Z],[Jy_X, Jy_Y, Jy_Z]])
        A_list.append(A)
        B_list.append(B)

    retVal=[_actualPointCoords2D, (projectedPointCoords2D_list), (A_list), (B_list)]
    return retVal
    print("KOKOKOKO KLOF!!!!")
    print(A)
    print(B)
    kokon=_projectedPointCoords2D-_actualPointCoords2D




def read_bal_data(file_name):
    with bz2.open(file_name, "rt") as file:
        n_cameras, n_points, n_observations = map(
            int, file.readline().split())

        camera_indices = np.empty(n_observations, dtype=int)
        point_indices = np.empty(n_observations, dtype=int)
        points_2d = np.empty((n_observations, 2))

        for i in range(n_observations):
            camera_index, point_index, x, y = file.readline().split()
            camera_indices[i] = int(camera_index)
            point_indices[i] = int(point_index)
            points_2d[i] = [float(x), float(y)]

        camera_params = np.empty(n_cameras * 9)
        for i in range(n_cameras * 9):
            camera_params[i] = float(file.readline())
        camera_params = camera_params.reshape((n_cameras, -1))

        points_3d = np.empty(n_points * 3)
        for i in range(n_points * 3):
            points_3d[i] = float(file.readline())
        points_3d = points_3d.reshape((n_points, -1))

    return camera_params, points_3d, camera_indices, point_indices, points_2d


def rotate(points, rot_vecs):
    """Rotate points by given rotation vectors.

    Rodrigues' rotation formula is used.
    """
    theta = np.linalg.norm(rot_vecs, axis=1)[:, np.newaxis]
    with np.errstate(invalid='ignore'):
        v = rot_vecs / theta
        v = np.nan_to_num(v)
    dot = np.sum(points * v, axis=1)[:, np.newaxis]
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    return cos_theta * points + sin_theta * np.cross(v, points) + dot * (1 - cos_theta) * v

def project(points, camera_params):
    """Convert 3-D points to 2-D by projecting onto images."""
    points_proj = rotate(points, camera_params[:, :3])
    points_proj += camera_params[:, 3:6]
    points_proj = points_proj[:, :2] / points_proj[:, 2, np.newaxis]
    f = camera_params[:, 6]
    u0 = camera_params[:, 7]
    v0 = camera_params[:, 8]
    k1 = camera_params[:, 9]
    k2 = camera_params[:, 10]
    n = np.sum(points_proj**2, axis=1)
    r = 1 + k1 * n + k2 * n**2
    points_proj *= (r * f)[:, np.newaxis]
    points_proj+=(np.array([u0, v0]).transpose())
    return points_proj


def fun(params, n_cameras, n_points, camera_indices, point_indices, points_2d):
    """Compute residuals.

    `params` contains camera parameters and 3-D coordinates.
    """
    camera_params = params[:n_cameras * 11].reshape((n_cameras, 11))
    points_3d = params[n_cameras * 11:-5].reshape((n_points, 3))
    common_camera_params=params[-5:]
    common_camera_params=np.hstack(([0,0,0,0,0,0],common_camera_params))
    points_proj = project(points_3d[point_indices], common_camera_params+camera_params[camera_indices])
    return (points_proj - points_2d).ravel()

def bundle_adjustment_sparsity(n_cameras, n_points, camera_indices, point_indices, camera_intrinsics_base_mask,camera_intrinsics_frame_mask):
    '''
    Generate a sparsity jacobian for the bundle adjustment.
    :param n_cameras: Number of positions, from which the points were observed. In an optimal case, it's the total number of frames of the footage.
    :param n_points: Number of observed 3D points. Equal to the total number of detected unique markers times 4
    :param camera_indices: A
    :param point_indices:
    :return:
    '''
    camera_intrinsics_count=5#focal length, principal x, principal y, dist 1, dist 2
    camera_extrinsics_count=6#rotation x 3, translation x 3
    m = camera_indices.size * 2
    n = n_cameras * (camera_extrinsics_count+camera_intrinsics_count) + n_points * 3 + camera_intrinsics_count
    A = lil_matrix((m, n), dtype=int)

    i = np.arange(camera_indices.size)
    camera_params_frame_mask=np.hstack(([1,1,1,1,1,1],camera_intrinsics_frame_mask))

    for s in range(len(camera_params_frame_mask)):
        A[2 * i, camera_indices * (camera_extrinsics_count+camera_intrinsics_count) + s] = camera_params_frame_mask[s]#even rows, camera parameters
        A[2 * i + 1, camera_indices * (camera_extrinsics_count+camera_intrinsics_count) + s] = camera_params_frame_mask[s]#odd rows, camera parameters

    for s in range(3):
        A[2 * i, n_cameras * (camera_extrinsics_count+camera_intrinsics_count) + point_indices * 3 + s] = 1#even rows, point indices
        A[2 * i + 1, n_cameras * (camera_extrinsics_count+camera_intrinsics_count) + point_indices * 3 + s] = 1#odd rows, point indices

    for s in range(len(camera_intrinsics_base_mask)):
        A[2*i,n_cameras * (camera_extrinsics_count+camera_intrinsics_count)+(n_points*3)+s ]=camera_intrinsics_base_mask[s]#even rows, base camera intrinsics
        A[2*i+1,n_cameras * (camera_extrinsics_count+camera_intrinsics_count)+(n_points*3)+s ]=camera_intrinsics_base_mask[s]#odd rows, base camera intrinsics

    return A

def mainFunc(camera_params,points_3d,camera_indices,point_indices,points_2d):
    n_cameras = camera_params.shape[0]
    n_points = points_3d.shape[0]

    n = 11 * n_cameras + 3 * n_points
    m = 2 * points_2d.shape[0]

    print("n_cameras: {}".format(n_cameras))
    print("n_points: {}".format(n_points))
    print("Total number of parameters: {}".format(n))
    print("Total number of residuals: {}".format(m))
    x0 = np.hstack((camera_params.ravel(), points_3d.ravel(),[0,0,0,0,0]))

    f0 = fun(x0, n_cameras, n_points, camera_indices, point_indices, points_2d)
    plt.plot(f0)
    A = bundle_adjustment_sparsity(n_cameras, n_points, camera_indices, point_indices,[0,0,0,0,0],[1,1,1,1,1])
    t0 = time.time()
    res = least_squares(fun, x0, jac_sparsity=A, verbose=2, x_scale='jac', ftol=1e-5, method='trf',
                       args=(n_cameras, n_points, camera_indices, point_indices, points_2d))
    kokodak=res.x
    t1 = time.time()
    print("Optimization took {0:.0f} seconds".format(t1 - t0))
    plt.plot(res.fun)
    plt.show()





# See PyCharm help at https://www.jetbrains.com/help/pycharm/
