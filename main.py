import numpy as np
import cv2#pip install opencv-contrib-python
import time
import math

import cv2.aruco as aruco#no need to handle separately

from mpl_toolkits.mplot3d import Axes3D  #pip install matplotlib # noqa: F401 unused import

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sparseba as sba
from bundleAdjustmentHelpers import *
import random



def get_camera_matrix(wpx, hpx, wmm, fmm):
    retMat=np.zeros(shape=(3,3))
    sizePx=wmm/wpx#size of a pixel in mm
    retMat[0,0]=fmm/sizePx
    retMat[1,1]=fmm/sizePx
    retMat[0,2]=wpx/2
    retMat[1,2]=hpx/2
    retMat[2,2]=1
    print("Kokon")
    #print(retMat)
    return retMat

def get_random_transformation_matrix(rvec,tvec):
    rvec=(np.multiply((np.random.random_sample(rvec.shape)-(np.ones(rvec.shape)/2))*2,rvec))
    tvec=(np.multiply((np.random.random_sample(tvec.shape)-(np.ones(tvec.shape)/2))*2,tvec))
    return vectors_to_matrix(rvec,tvec)


def vectors_to_matrix(rvec, tvec):
    dest, jac=cv2.Rodrigues(rvec)
    tvec=np.transpose(tvec)
    dest=np.hstack((dest,tvec))
    dest=np.vstack((dest,np.array([0,0,0,1])))
    return dest

def get_first_index_of_nparray(theArray, theValue):
    result=np.where(theArray==theValue)
    result=result[0]
    if(len(result)==0):
        return None
    result=result[0]
    return result

'''
@brief return the 4 coordinates of the corners of a marker, in the coordinate system of the marker itself, as column vectors
'''
def get_corner_base_coordinates(markerSize):
    tempMat=np.matrix([[-markerSize/2,markerSize/2,0],[markerSize/2,markerSize/2,0],[markerSize/2,-markerSize/2,0],[-markerSize/2,-markerSize/2,0]])
    return np.transpose(tempMat)

'''
@brief given the corners of a marker as a matrix of column vectors and a 4-by-4 affinity matrix, transform the corners using the matrix and return their
transformed version as a matrix of column vectors
'''
def transform_corners_with_matrix(corners, mat):
    rows,cols=mat.shape
    ones=np.ones(shape=(1,cols))
    corners=np.vstack((corners,ones))
    tempMat=mat.dot(corners)
    tempMat=tempMat[:-1,:]
    return tempMat

def get_transformation_matrix_from_base_to_marker(rvecs, tvecs, ids, marker_matrix_pairs, target_marker_id):
    reference_marker_id=None#id (determined by the marker graphics) of the marker that is currently visible together with the target marker, and position of which wrt. base is known
    base_to_ref_matrix=None#the matrix, which transforms the coordinate system from base into the coordinate system of the reference marker

    #find some marker that has a known location wrt. base and is present in the "ids" list; that is, a marker with known position, which is on the frame together with the target_marker_id marker
    for marker in marker_matrix_pairs:
        if marker in ids and marker != target_marker_id:
            reference_marker_id=marker
            base_to_ref_matrix=(marker_matrix_pairs[marker])
            break;

    target_marker_index=get_first_index_of_nparray(ids, target_marker_id)#array index in the input arrays of the target marker
    reference_marker_index=get_first_index_of_nparray(ids, reference_marker_id)#array index in the input arrays of the reference marker

    if(reference_marker_index is None or target_marker_index is None):
        return None

    camera_to_ref_matrix=vectors_to_matrix(rvecs[reference_marker_index],tvecs[reference_marker_index])
    camera_to_base_matrix=camera_to_ref_matrix.dot(np.linalg.inv(base_to_ref_matrix))

    camera_to_target_marker_matrix=vectors_to_matrix(rvecs[target_marker_index],tvecs[target_marker_index])
    base_to_target_marker_matrix=(np.linalg.inv(camera_to_base_matrix)).dot(camera_to_target_marker_matrix)
    return base_to_target_marker_matrix


def get_transformation_matrix_from_base_to_marker_for_frame(frame_list,frame_number,marker_matrix_pairs, target_marker_id):
    return get_transformation_matrix_from_base_to_marker(frame_list[frame_number]['rvecs'],frame_list[frame_number]['tvecs'],frame_list[frame_number]['ids'], marker_matrix_pairs, target_marker_id)

def get_corners_of_marker_in_base_marker_coordinate_system(rvecs,tvecs, ids, base_marker_id, target_marker_id, markerSize):
    tempMat=get_transformation_matrix_from_base_to_marker(rvecs,tvecs,ids,base_marker_id,target_marker_id)
    if tempMat is None:
        return None
    corner_base_coordinates=get_corner_base_coordinates(markerSize)
    return transform_corners_with_matrix(corner_base_coordinates,tempMat)


'''
@brief given the base to target transformation matrix and a target marker id, iterate through all the frames; on each frame,
determine how the target marker would look on the camera (2D image), using the already-determined absolute positions of other markers on the frame;
the difference between how it would look like and how it looks adds to the reprojection error, which is then returned

'''
def get_reprojection_error(base_to_target_marker_matrix, marker_matrix_pairs, target_marker_id, camera_matrix, dist_coeffs, marker_size, frame_list, max_error):
    currentErr=0

    for frame_number in range(1,len(frame_list)):#get the sum of reprojection error for all frames
        rvecs=frame_list[frame_number]['rvecs']
        tvecs=frame_list[frame_number]['tvecs']
        corners=frame_list[frame_number]['corners']
        ids=frame_list[frame_number]['ids']

        for marker in marker_matrix_pairs:#for all the markers that have a known position wrt. base...
            if ids is None:#no markers on the current frame
                continue
            if marker in ids and marker != target_marker_id:#..and are present on the current frame
                reference_marker_id=marker
                base_to_ref_matrix=marker_matrix_pairs[marker]

                ref_marker_index =get_first_index_of_nparray(ids, reference_marker_id)
                target_marker_index=get_first_index_of_nparray(ids, target_marker_id)
                if(ref_marker_index is None or target_marker_index is None):
                    continue
                    #return None

                camera_to_ref_marker_matrix=vectors_to_matrix(rvecs[ref_marker_index],tvecs[ref_marker_index])
                camera_to_target_marker_matrix=camera_to_ref_marker_matrix.dot(np.linalg.inv(base_to_ref_matrix)).dot(base_to_target_marker_matrix)
                corners3D=get_corner_base_coordinates(marker_size)
                corners3D=transform_corners_with_matrix(corners3D,camera_to_target_marker_matrix)
                corners2DNew,jac=cv2.projectPoints(corners3D,(0,0,0),(0,0,0),camera_matrix,dist_coeffs)
                #print(corners[target_marker_index])
                corners[target_marker_index]=np.reshape(corners[target_marker_index],(2,-1))
                corners2DNew=np.reshape(corners2DNew,(2,-1))
                currentErr=currentErr+(np.sum(np.square(corners[target_marker_index]-corners2DNew)))
                if(currentErr>max_error):#pruning in case the error is already too large; we don't need to know just how large it is!
                    return currentErr
    return currentErr

'''
@brief take a list of pairs between a marker id and its position wrt. base, list of frame data (containing the positions of marker corners, transformation vectors of the
markers wrt. camera, etc.), camera matrix and the distortion coefficents of the camera; use this information to fill in the unknown positions of the markers wrt. base . 
The function iterates over all frames and for each of them calculates the best guess for the position of the markers wrt. base, based on the currently available info.
It then uses reprojection error to determine which frame gives the best guess; the newly estimated marker is added to the marker-matrix pair list, together with the
estimated transformation matrix.
'''
def update_marker_matrix_pairs(marker_matrix_pairs, valid_markers_list,frame_list,camera_matrix,distortion_coefficients):
    for marker_id in valid_markers_list:#for all the potentially valid markers
        bestMat=0
        bestErr=10000000
        currentErr=0
        t = time.time()
        processedFrames=0

        if marker_id in marker_matrix_pairs:
            continue

        for i in range(len(frame_list)):
            if frame_list[i]['ids'] is None:
                continue
            if marker_id not in frame_list[i]['ids']:#the specified marker is not on the current frame
                continue

            transMat=get_transformation_matrix_from_base_to_marker_for_frame(frame_list,i,marker_matrix_pairs,marker_id)
            if not transMat is None:
                processedFrames=processedFrames+1
                guessMat=transMat
                for j in range(1):
                    currentErr=(get_reprojection_error(guessMat,marker_matrix_pairs,marker_id,camera_matrix,distortion_coefficients,0.1,frame_list,bestErr))
                    if currentErr<bestErr:
                        bestErr=currentErr
                        bestMat=guessMat
                    #guessMat=transMat.dot(get_random_transformation_matrix(np.array([[0.1,0.1,0.1]]),np.array([[0.1,0.1,0.1]])))

        #print(processedFrames)
        #print(bestErr)
        #print(bestMat)
        marker_matrix_pairs[marker_id]=bestMat
    print(time.time()-t)

def plot_markers(marker_matrix_pairs, marker_size,ax):
    for pair in marker_matrix_pairs:
        marker_matrix=pair['matrix']
        corners=get_corner_base_coordinates(marker_size)
        newCorners=transform_corners_with_matrix(corners,marker_matrix)
        ax.scatter(newCorners[0,:],newCorners[1,:],newCorners[2,:])

def markers_to_3d_points(marker_matrix_pairs, marker_size):
    points=[]
    newList=sorted(marker_matrix_pairs,key= lambda d:d['marker_id'])
    for pair in newList:
        marker_matrix=pair['matrix']
        corners=get_corner_base_coordinates(marker_size)
        newCorners=transform_corners_with_matrix(corners,marker_matrix)
        newCorners=newCorners.transpose()
        points.extend(newCorners.tolist())
    return points


def prepare_2d_3d_correspondances_for_pnp_solver(frame_list,frame_number,marker_matrix_pairs,marker_size):
    frame_data=frame_list[frame_number]
    ids=frame_data['ids'].ravel()
    corners=frame_data['corners']

    #get rid of all the markers that we are not interested in, but that were still detected
    valid_marker_ids=marker_matrix_pairs.keys()
    valid_indices=[i for i,elem in enumerate(ids)if elem in valid_marker_ids]
    marker_ids_to_point_ids_dict={}#a lookup dictionary to convert the marker indices to point indices
    for i in range(len(valid_indices)):
        marker_ids_to_point_ids_dict[valid_indices[i]]=i*4

    points2D=[]
    points3D=[]
    pointIDs=[]
    ids=ids[valid_indices]
    corners=np.array(corners)
    corners=corners[valid_indices]
    sorted_indices=np.argsort(ids.ravel())
    ids=ids[sorted_indices]
    corners=corners[sorted_indices]
    for i in range (len(ids)):
        #the points are stored in an odd way; sometimes, the shape is (2,4), other times (4,2,1), seemingly randomly...
        points2Dtemp=corners[i]
        points2Dtemp=np.squeeze(points2Dtemp)
        points2Dtemp=points2Dtemp.reshape((4,2))

        if(len(points2D)==0):
            points2D=points2Dtemp
        else:
            points2D=np.vstack((points2D,points2Dtemp))

        for j in range(4):
            temp=ids[i]
            pointIDs.append(marker_ids_to_point_ids_dict[temp]+j)

        for marker in marker_matrix_pairs:#todo: maybe exit the loop once we find it?
            if  marker == ids[i]:
                reference_marker_id=(marker)
                base_to_target_matrix=(marker_matrix_pairs[marker])
                points3Dtemp=get_corner_base_coordinates(marker_size)
                points3Dtemp=transform_corners_with_matrix(points3Dtemp,base_to_target_matrix)
                points3Dtemp=np.multiply(points3Dtemp,[[1,1,1,1],[-1,-1,-1,-1],[1,1,1,1]])#todo: check if this is correct
                if(len(points3D)==0):
                    points3D=points3Dtemp
                else:
                    points3D=np.hstack((points3D,points3Dtemp))
    return (pointIDs, points2D, points3D)


def matrix_to_xyz_euler(rotation_matrix):
    rm=rotation_matrix
    if rm[2,0] < 1:
        if rm[2,0] > -1:
            ty=math.asin(-rm[2,0])
            tz=math.atan2(rm[1,0],rm[0,0])
            tx=math.atan2(rm[2,1],rm[2,2])
        else:
            ty=math.pi/2
            tz=-math.atan2(-rm[1,2],rm[1,1])
            tx=0
    else:
        ty=-math.pi/2
        tz=math.atan2(-rm[1,2],rm[1,1])
        tx=0
    return np.array([tx,ty,tz])

def matrix_to_translation_vector(the_matrix):
    return the_matrix[0:-1,3]

def write_line_to_file(the_file,the_line):
    the_file.write((','.join('%0.8f'%x for x in the_line))+'\n')

def generate_tracking_line_for_frame(frame_data,frame_number,marker_matrix_pairs,marker_size,camera_matrix,distortion_coefficients):
    if(frame_number==174):
        print("kokodaaaak")
    points2D,points3D=prepare_2d_3d_correspondances_for_pnp_solver(frame_data,frame_number,marker_matrix_pairs,marker_size)
    points3D=np.array(np.transpose(points3D))
    retval, rvecCam,tvecCam=cv2.solvePnP(points3D,points2D,camera_matrix,distortion_coefficients)
    tmat=vectors_to_matrix(np.transpose(rvecCam),np.transpose(tvecCam))
    prepareBundleAdjustmentParamsForSingleFrame(points3D,points2D,rvecCam,tvecCam,camera_matrix,distortion_coefficients)

    # cor=get_corner_base_coordinates(0.095)
  #  cor=np.multiply(cor,[[1,1,1,1],[-1,-1,-1,-1],[1,1,1,1]])
    #cor=transform_corners_with_matrix(cor,tmat)
    #print(cor)
    returnVec=np.hstack((np.array(frame_number),matrix_to_xyz_euler(tmat),matrix_to_translation_vector(tmat)))
    if(returnVec[6]<0):#object behind camera, as opposed to in front of the camera
        returnVec[1:7]=-returnVec[1:7]
        returnVec[3]=-(returnVec[3]+math.pi)

    returnVec[2]=-returnVec[2]#flip y rot axis
    returnVec[3]=-returnVec[3]#flip z rot axis
    returnVec[5]=-returnVec[5]#flip y trans axis
    returnVec[6]=-returnVec[6]#flip z trans axis

    #print(retval)
    #print(rvecCam)
    #print(tvecCam)
    #print(returnVec)
    return returnVec

def prepare_data_for_bundle_adjustment(frame_data, marker_matrix_pairs, marker_size, camera_matrix, distortion_coefficients):
    retPoints2D=[]
    retPoints2DProjected=[]
    retAList=[]
    retBList=[]
    pointIDs=[]
    viewportIDs=[]

    global framesCount
    for i in range(framesCount):
        pointIDsFrame, points2D,points3D=prepare_2d_3d_correspondances_for_pnp_solver(frame_data,i,marker_matrix_pairs,marker_size)
        points3D=np.array(np.transpose(points3D))
        retval, rvecCam,tvecCam=cv2.solvePnP(points3D,points2D,camera_matrix,distortion_coefficients)
        #rvecCam[2]=rvecCam[2]+0.1
        [points2DActualFrame, points2DProjectedFrame,Aframe,Bframe]=prepareBundleAdjustmentParamsForSingleFrame(points3D,points2D,rvecCam,tvecCam,camera_matrix,distortion_coefficients)
        retPoints2D.extend(points2DActualFrame)
        retPoints2DProjected.extend(points2DProjectedFrame)
        retAList.extend(Aframe)
        retBList.extend(Bframe)
        viewportIDsFrame=[i]*len(pointIDsFrame)
        viewportIDs.extend(viewportIDsFrame)
        pointIDs.extend(pointIDsFrame)
    return [pointIDs, viewportIDs, retPoints2D, retPoints2DProjected,retAList, retBList]

def get_camera_vector(_rvec,_tvec, _camera_matrix):
    retVal=[]
    retVal.extend(np.ravel(_rvec))
    retVal.extend(np.ravel(_tvec))
    f=camera_matrix[0,0]
    u0=camera_matrix[0,2]
    v0=camera_matrix[1,2]
    #TODO: distortion coefficients
    retVal.extend([f])
    retVal.extend([u0])
    retVal.extend([v0])
    retVal.extend([0,0])
    return retVal

def prepare_data_for_bundle_adjustment2(frame_data, marker_matrix_pairs, marker_size, camera_matrix, distortion_coefficients):
    retPoints2D=[]
    retPoints3D=[]
    retPoints2DProjected=[]
    retAList=[]
    retBList=[]
    pointIDs=[]
    viewportIDs=[]
    camVectorList=[]

    global framesCount
    for i in range(framesCount):
        pointIDsFrame, points2D,points3D=prepare_2d_3d_correspondances_for_pnp_solver(frame_data,i,marker_matrix_pairs,marker_size)
        points3D=np.array(np.transpose(points3D))
        retPoints3D.extend(points3D)
        retval, rvecCam,tvecCam=cv2.solvePnP(points3D,points2D,camera_matrix,distortion_coefficients)
        camVector=get_camera_vector(rvecCam,tvecCam,camera_matrix)
        #rvecCam[2]=rvecCam[2]+0.1
        #[points2DActualFrame, points2DProjectedFrame,Aframe,Bframe]=prepareBundleAdjustmentParamsForSingleFrame(points3D,points2D,rvecCam,tvecCam,camera_matrix,distortion_coefficients)
        retPoints2D.extend(points2D)
        #retPoints2DProjected.extend(points2DProjectedFrame)
        #retAList.extend(Aframe)
        #retBList.extend(Bframe)
        viewportIDsFrame=[i]*len(pointIDsFrame)
        viewportIDs.extend(viewportIDsFrame)
        pointIDs.extend(pointIDsFrame)
        camVectorList.append(camVector)
    return [np.array(pointIDs), np.array(viewportIDs), np.array(retPoints2D), np.array(retPoints3D),np.array(camVectorList)]


def testVideoConversion():
    vidcap = cv2.VideoCapture("L:\\personal\\tracker\\newCode\\footage\\tracking_footage.mp4")
    #vidcap = cv2.VideoCapture("L:\\personal\\tracker\\CameraParameterEstimatorCharuko\\workdir\\realCharuko\\arucoboard.mp4")
    success,image = vidcap.read()
    count = 0
    success = True
    while success:
      success,image = vidcap.read()
      cv2.imwrite("L:\\personal\\tracker\\newCode\\footage\\tracking_footage\\frame%04d.jpg"%count,image)
      #cv2.imwrite("L:\\personal\\tracker\\CameraParameterEstimatorCharuko\\workdir\\realCharuko\\pngs\\frame%d.png"%count,image)
      #cv2.imwrite("L:\\personal\\tracker\\newCode\\footage\\calib_footage\\frame%d.jpg" % count, image)     # save frame as JPEG file
      if cv2.waitKey(10) == 27:                     # exit if Escape is hit
          break
      count += 1




camera_matrix=get_camera_matrix(1920,1080,36,35)
camera_vector=[get_camera_vector([0,0,0],[0,0,0],camera_matrix)]
p3D=[[1.0,2,3]]

corners2DNew, jac = cv2.projectPoints(np.array(p3D), (0, 0, 0), (0, 0, 0), camera_matrix, np.array([0.0,0,0,0]))
corners2DNew2=project(p3D,np.array(camera_vector))

#prepareBundleAdjustmentParamsForSinglePoint([1, 2, 3],[-1,1],[1,2,3],[1,2,3],camera_matrix,None)
#exit()



#print(get_camera_matrix(640,480,36,35))

#cap = cv2.VideoCapture(0)


#camera_matrix=get_camera_matrix(1920,1080,36,35)
#kokodak()
#testVideoConversion()
#kokodak()

#fig = plt.figure()
#fig = plt.figure(figsize=plt.figaspect(1)*1.5) #Adjusts the aspect ratio and enlarges the figure (text does not enlarge)
#ax = fig.add_subplot(111, projection='3d')
#ax = fig.gca(projection='3d')
#ax.auto_scale_xyz([0,2], [0,2], [0,2])
aruco_dict = aruco.Dictionary_get(aruco.DICT_4X4_50)
parameters =  aruco.DetectorParameters_create()
#parameters.cornerRefinementMethod=aruco.CORNER_REFINE_APRILTAG
parameters.cornerRefinementMethod=aruco.CORNER_REFINE_CONTOUR
#parameters.cornerRefinementMethod=aruco
#camera_matrix=get_camera_matrix(1920,1080,36,35)
camera_matrix=get_camera_matrix(1920,1080,36,27.2)
camera_matrix=np.array([[1449.2,0.0,1007.2],[0.0,1449.2,542.28],[0.0,0.0,1.0]])
#distCoeffs=np.array([[3.916,-112.5,0.0,0.0,460.7]])
#distCoeffs=np.array([[3.916],[-112.5],[-0.00363],[-0.001854]])
distCoeffs=np.zeros(shape=(1,4))

frameData={}
frameDataList=[]

framesCount=200

for i in range(framesCount):
    frameData=frameData.copy()
    #filename="L:\\personal\\tracker\\testAnimCube\\%04d.png"%(i+1)
    #filename="L:\\personal\\tracker\\testAnimDetermineAxes\\%04d.png"%(i+1)
    #filename="L:\\personal\\tracker\\testAnimCrazyMotion3D\\%04d.png"%(i+1)
    filename="D:\\personal\\AprilTools\\AprilTools\\footage\\tracking_footage\\frame%04d.jpg"%(i+1)
    #filename="L:\\personal\\tracker\\newCode\\footage\\tracking_footage\\frame%04d.jpg"%(i+1)
    #print(filename)
    frame = cv2.imread(filename)#cap.read()
    gray = frame#cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    corners, ids, rejectedImgPoints = aruco.detectMarkers(gray, aruco_dict, parameters=parameters)
    rvecs, tvecs, _objPoints = aruco.estimatePoseSingleMarkers(corners,0.1,camera_matrix,distCoeffs)
    #newCorners=get_corners_of_marker_in_base_marker_coordinate_system(rvecs,tvecs,ids,0,2,1)
    #print(newCorners)
    frameData['corners']=list(corners)
    frameData['ids']=ids
    frameData['rvecs']=rvecs
    frameData['tvecs']=tvecs
    frameDataList.append(frameData)
    print(i)
    #if not newCorners is None:
    #    ax.scatter(newCorners[0,:],newCorners[1,:],newCorners[2,:])

markerMatrixPairs={}
markerMatrixPairs[0]=np.identity(4);#TODO:make this user-adjustable
#base_matrix_pair={}#marker id vs. the transformation matrix from base to the marker; the base can be a base marker or world
valid_markers_list=[0,1,2,3]
#markerMatrixPairList=[]
#base_matrix_pair['marker_id']=0
#base_matrix_pair['matrix']=np.identity(4)
#markerMatrixPairList.append(base_matrix_pair)

update_marker_matrix_pairs(markerMatrixPairs, valid_markers_list,frameDataList,camera_matrix,distCoeffs)
update_marker_matrix_pairs(markerMatrixPairs, valid_markers_list,frameDataList,camera_matrix,distCoeffs)
update_marker_matrix_pairs(markerMatrixPairs, valid_markers_list,frameDataList,camera_matrix,distCoeffs)
update_marker_matrix_pairs(markerMatrixPairs, valid_markers_list,frameDataList,camera_matrix,distCoeffs)
update_marker_matrix_pairs(markerMatrixPairs, valid_markers_list,frameDataList,camera_matrix,distCoeffs)
update_marker_matrix_pairs(markerMatrixPairs, valid_markers_list,frameDataList,camera_matrix,distCoeffs)

#markerPoints=np.array(markers_to_3d_points(markerMatrixPairList,0.1))
#ax.scatter(markerPoints[:,0],markerPoints[:,1],markerPoints[:,2])
#markerMatrixPairList[0]['matrix']= markerMatrixPairList[0]['matrix'].dot(get_random_transformation_matrix(np.array([[0.5,0,0]]),np.array([[0.03,0.02,0.06]])))
#markerMatrixPairList[1]['matrix']= markerMatrixPairList[1]['matrix'].dot(get_random_transformation_matrix(np.array([[0.5,0,0]]),np.array([[0.03,0.02,0.06]])))
#markerMatrixPairList[2]['matrix']= markerMatrixPairList[2]['matrix'].dot(get_random_transformation_matrix(np.array([[0.5,0,0]]),np.array([[0.03,0.02,0.06]])))
#markerMatrixPairList[3]['matrix']= markerMatrixPairList[3]['matrix'].dot(get_random_transformation_matrix(np.array([[0.5,0,0]]),np.array([[0.03,0.02,0.06]])))
#markerPoints=np.array(markers_to_3d_points(markerMatrixPairList,0.1))
#ax.scatter(markerPoints[:,0],markerPoints[:,1],markerPoints[:,2])
#plot_markers(markerMatrixPairList, 0.1,ax)


#filename="L:\\personal\\tracker\\testNewTracking.txt"
filename = "D:\\personal\\AprilTools\\AprilTools\\testNewTracking.txt"
#file=open(filename,'w')

#for i in range(framesCount):
#    line=generate_tracking_line_for_frame(frameDataList,i,markerMatrixPairList,0.1,camera_matrix,distCoeffs)
#    write_line_to_file(file,line)
#    print(i)
#file.close()
#kokodak=get_random_transformation_matrix(np.array([0,0,0]),np.array([0.2,0.2,0.2]))
#markerMatrixPairList[0]['matrix']= markerMatrixPairList[0]['matrix'].dot(get_random_transformation_matrix(np.array([[0,0,0]]),np.array([[0.2,0.2,0.2]])))
#camera_matrix=np.array([[1000.2,0.0,107.2],[0.0,1000.2,1542.28],[0.0,0.0,1.0]])
[point_indices,viewport_indices,points2D,points3D,camVector]=prepare_data_for_bundle_adjustment2(frameDataList,markerMatrixPairs,0.1,camera_matrix,distCoeffs)
mainFunc(camVector,points3D,viewport_indices,point_indices,points2D)
#kokodak=sba.core.SBA(viewport_indices,point_indices)
#[delta_a,delta_b]=kokodak.compute(np.array(x_true),np.array(x_pred),np.array(A),np.array(B))
#markerPoints=markerPoints+delta_b
#ax.scatter(markerPoints[:,0],markerPoints[:,1],markerPoints[:,2])
plt.show()

'''
print(rvecs)
print(tvecs)
theMat=(vectors_to_matrix(rvecs[0], tvecs[0]))
print(theMat)
print("kokokokodaaaaak")
corners2=get_corner_base_coordinates(1)
print(corners2)
print("klof")
print(transform_corners_with_matrix(corners2,theMat))

origin=np.transpose(np.array([[0,0,0,1],[0,0,1,1]]))
newOrigin=theMat.dot(origin)
print(newOrigin)

'''

#gray = aruco.drawDetectedMarkers(gray, corners)



cv2.imshow('frame',gray)
cv2.waitKey()

#cv2.destroyAllWindows()
