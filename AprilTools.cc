/*
 * Coded by thegoodhen 2019, theegodhen@gmail.com
 * Released under GNU GPL 3, see file LICENCE for more info.
 *
 * Focal length estimation based on "Using Vanishing Points for Camera Calibration and Coarse 3D Reconstruction from a Single Image" by E. Guillou, D. Meneveaux, E. Maisel, K. Bouatouch.
 * Based around the AprilTags library by University of Michigan.
 * Copyright notice from AprilTags library below:
 * */


/* Copyright (C) 2013-2016, The Regents of The University of Michigan.
All rights reserved.

This software was developed in the APRIL Robotics Lab under the
direction of Edwin Olson, ebolson@umich.edu. This software may be
available under alternative licensing terms; contact the address above.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the Regents of The University of Michigan.
*/


#include <iostream>

#include <opencv2/opencv.hpp>

//#include <opencv2/core.hpp>
//#include <opencv2/imgcodecs.hpp>
//#include <opencv2/highgui.hpp>

#include <string.h>
#include <vector>
#include <regex>

extern "C"{
#include <dirent.h>
#include <stdio.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#include "apriltag.h"
#include "tag36h11.h"
#include "tag25h9.h"
#include "tag16h5.h"
#include "tagCircle21h7.h"
#include "tagCircle49h12.h"
#include "tagCustom48h12.h"
#include "tagStandard41h12.h"
#include "tagStandard52h13.h"

#include "common/getopt.h"
#include "common/image_u8.h"
#include "common/image_u8x4.h"
#include "common/pjpeg.h"
#include "common/zarray.h"
#include "common/g2d.h"

#include "apriltag_pose.h"
}

uint8_t compareFilenames(std::string, std::string);
bool isFileAcceptable(char*);
bool has_suffix(const std::string &, const std::string &);
int getFilenameNumber(std::string);
bool quick=false;

double getPixelFFromIntersections(double, double* , double*, double* );

double fmm=-1;//focal length in mm
double fpx=-1;//focal length in px
double sensorWidth=-1;//sensor width in mm
double sensorWidthPx=-1;//sensor width in px


/*
 * Fill the remaining intrinsic parameters of the camera, given the known ones*/
void fillIntrinsics()
{
  double pxSize=-1;//pixel size in mm
  if(fmm!=-1 && fpx!=-1)
  {
    pxSize=fmm/fpx;
    sensorWidth=sensorWidthPx*pxSize;
    return;
  }
  if(fmm!=-1 && sensorWidth!=-1)
  {
    pxSize=sensorWidth/sensorWidthPx;
    fpx=fmm/pxSize;
    return;
  }
  if(fpx!=-1 && sensorWidth!=-1)
  {
    pxSize=sensorWidth/sensorWidthPx;
    fmm=fpx*pxSize;
    return;
  }
}

/* Print the intrinsic camera parameters to stdout
 * */
void printIntrinsics()
{
  printf("Focal length: %.2f mm (%.2f px); sensor width: %.2f mm\r\n",fmm,fpx,sensorWidth);
}

/*
 * Calculate the pixel focal length, based on the detection object provided and the principal point.
 * */
double getPixelF(apriltag_detection_t* det, double* principal)
{
 g2d_line_t lAB;
 g2d_line_t lCD;
 g2d_line_t lAD;
 g2d_line_t lBC;
 double pA[2]={det->p[0][0],det->p[0][1]};
 double pB[2]={det->p[1][0],det->p[1][1]};
 double pC[2]={det->p[2][0],det->p[2][1]};
 double pD[2]={det->p[3][0],det->p[3][1]};

 g2d_line_init_from_points(&lAB, pA, pB);
 g2d_line_init_from_points(&lCD, pC, pD);
 g2d_line_init_from_points(&lAD, pA, pD);
 g2d_line_init_from_points(&lBC, pB, pC);

 double m[2];
 double n[2];
 int int1Found=g2d_line_intersect_line(&lAB, &lCD, m);
 int int2Found=g2d_line_intersect_line(&lAD, &lBC, n);

 //If the intersections are "too far", this suggest the lines being close to parallel. In such cases, the result is numerically unstable and should be discarded.
 double dist1=sqrt((det->c[0]-m[0])*(det->c[0]-m[0])+(det->c[1]-m[1])*(det->c[1]-m[1]));
 double dist2=sqrt((det->c[0]-n[0])*(det->c[0]-n[0])+(det->c[1]-n[1])*(det->c[1]-n[1]));
 double thres=5*principal[0];//since the principal point is in the middle of the image, we can use it to dynamically adjust the threshold distance between the tag center and principal points, beyond which the result gets discarded, based on the resolution


 if(int1Found && int2Found && dist1<thres && dist2<thres)
 {
   return getPixelFFromIntersections(1,principal,m,n);
 }
 return -1;
}


/**
 * Given alpha(pixel aspect ratio), principal point (usually the center of the image) and two intersections (m, n), calculate the focal length in pixels and return it. 
 * Focal length estimation based on "Using Vanishing Points for Camera Calibration and Coarse 3D Reconstruction from a Single Image" by E. Guillou, D. Meneveaux, E. Maisel, K. Bouatouch.
 * */
double getPixelFFromIntersections(double alpha, double* principal, double*m, double* n)
{
  double um=m[0];
  double vm=m[1];
  double un=n[0];
  double vn=n[1];
  double u0=principal[0];
  double v0=principal[1];
  double A=(1/alpha)*(um-u0)*(un-u0);
  double B=(v0-vm)*(v0-vn);
  return sqrt(-(A)-(B));
}

/**
 * Given a std::vector of double values, calculate the median value and return it.
 * */
double calculateMedian(std::vector<double>theVec)
{
  if(theVec.size()==0)
  {
    return -1;
  }
  sort(theVec.begin(),theVec.end());
  if(theVec.size()%2==0)//even number of elements
  {
    return (theVec[(theVec.size()/2)-1]+theVec[(theVec.size()/2)-1])/2;
  }
  else
  {
    return theVec[theVec.size()/2];
  }
}


/**
 * Given a path to some folder, find files that contain a number in their name and then sort them, outputting the sorted vector of string filenames.
 * */
std::vector<std::string> getSortedFilenames(const char* path)
{
  std::vector<std::string> returnVec;
  DIR *d;
  struct dirent *dir;
  d = opendir(path);
  printf("\r\nsearching for files...\r\n");
  if (d)
  {
    while ((dir = readdir(d)) != NULL)
    {
      if(isFileAcceptable(dir->d_name))
      {
	returnVec.push_back(dir->d_name);
      }
      //printf("%s\n", dir->d_name);
    }
	  closedir(d);
  }
  //printf("sorting...\r\n");
  sort(returnVec.begin(),returnVec.end(),compareFilenames);
  //printf("sorted\r\n");
  return returnVec;
}

/**
 * Given two filenames (as std::string), which contain a number in them, parse the two numbers and then return 1 if the number in the first one is smaller and 0 otherwise.
 * */
uint8_t compareFilenames(std::string str1, std::string str2)
{
  int number;
  int number2;

  std::string output = std::regex_replace(
      str1,
      std::regex("[^0-9]*([0-9]+).*"),
      std::string("$1")
      );

  std::string output2 = std::regex_replace(
      str2,
      std::regex("[^0-9]*([0-9]+).*"),
      std::string("$1")
      );

  number=atoi(output.c_str());
  number2=atoi(output2.c_str());

  if(number>=number2)
  {
    return 0;
  }

  return 1;
}

int getFilenameNumber(std::string fileName)
{
  std::string output = std::regex_replace(
      fileName,
      std::regex("[^0-9]*([0-9]+).*"),
      std::string("$1")
      );
 int number=atoi(output.c_str());
 return number;
}


/**
 *
 * Convert rotation matrix to Euler angles, returning them as a matd_t* vector.
 **/
matd_t* getEulers(matd_t* Mr)
{
  double rx,ry,rz;
 if (MATD_EL(Mr,2,0)<1)
 {
   if(MATD_EL(Mr,2,0)>-1)
   {
     ry=asin(-MATD_EL(Mr,2,0));
     rz=atan2(MATD_EL(Mr,1,0),MATD_EL(Mr,0,0));
     rx=atan2(MATD_EL(Mr,2,1),MATD_EL(Mr,2,2));
   }
   else
   {
     ry=M_PI/2;
     rz=-atan2(-MATD_EL(Mr,1,2),MATD_EL(Mr,1,1));
     rx=0;
   }
 }
 else
 {
    ry=-M_PI/2;
    rz=atan2(-MATD_EL(Mr,1,2),MATD_EL(Mr,1,1));
    rx=0;
 }
 double data[]={rx,ry,rz};
 return matd_create_data(1,3,data);
}


/**
 * Attempt to detect tags in the image at 100% scale. If this fails, scale the image down twice, attempt to do it again. 
 * If it still fails, scale down to quarter the original size and reattempt.
 * Detecting the tag on a lower resolution image is more robust to motion blur, but yields less accurate results.
 *
 * Update the output arguments, which indicate how many images had to be rescaled to get the detection and on how many images the detection failed completely.
 * */
zarray_t* detectTagsRobust(apriltag_detector_t* td, image_u8_t* im, int* blurryPicsCount, int* badPicsCount)
{
    if(quick)
    {
      td->quad_decimate=2;
    }
    else
    {
      td->quad_decimate=1;//TODO: fix the refining accordingly
    }
    zarray_t *detections = apriltag_detector_detect(td, im);
    bool picIsBlurry=false;
    if(zarray_size(detections)==0)//no tags were detected
    {
      apriltag_detections_destroy(detections);//not sure if necessary...
      //printf("Failed to detect the tag, lowering resolution to 1/2 to counter blur...\r\n");
      td->quad_decimate=2;
      detections = apriltag_detector_detect(td, im);
      picIsBlurry=true;
    }

    if(zarray_size(detections)==0)//no tags were detected
    {
      //printf("Failed to detect the tag, lowering resolution to 1/4 to counter blur...\r\n");
      apriltag_detections_destroy(detections);//not sure if necessary...
      td->quad_decimate=4;
      detections = apriltag_detector_detect(td, im);
      picIsBlurry=true;
    }

    if(zarray_size(detections)==0)//no tags were detected
    {
      (*badPicsCount)++;
    }
    else if (picIsBlurry)
    {
      (*blurryPicsCount)++;
    }
    return detections;

}

//from https://stackoverflow.com/questions/20446201/how-to-check-if-string-ends-with-txt/20446257
/*
 * Return whether the std::string provided ends with the provided suffix.
*/
bool has_suffix(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

/**
 * Given a filename as a char*, decide whether it is acceptable as a frame of a sequence. To pass the test, the file needs to have a valid image suffix and has to contain a number.
 * */
bool isFileAcceptable(char* fileName)
{
  std::regex theRegex(".*[0-9].*");
  std::string fileNameStr=std::string(fileName);
  bool doesFileHaveNumber=std::regex_match(std::string(fileNameStr),theRegex);
  return (doesFileHaveNumber)&&((has_suffix(fileNameStr,".jpg"))||(has_suffix(fileNameStr,".JPG"))||(has_suffix(fileNameStr,".bmp"))||(has_suffix(fileNameStr,".BMP"))||(has_suffix(fileNameStr,".png"))||(has_suffix(fileNameStr,".PNG"))||(has_suffix(fileNameStr,".jpeg"))||(has_suffix(fileNameStr,".JPEG"))||(has_suffix(fileNameStr,".dib"))||(has_suffix(fileNameStr,".DIB"))||(has_suffix(fileNameStr,".jp2")));
}

int main(int argc, char *argv[])
{
  printf("\r\n");
  printf("\r\n");
  //imageTest();
  getopt_t *getopt = getopt_create();

  getopt_add_bool(getopt, 'h', "help", 0, "Show this help");
  getopt_add_string(getopt, 'p', "path", "", "Path to the source image sequence");
  getopt_add_double(getopt, 'f', "focal-length-mm", "-1", "Focal length in mm");
  getopt_add_double(getopt, 'F', "focal-length-pixels", "-1", "Focal length in pixels");
  getopt_add_double(getopt, 'w', "sensor-width", "-1", "Camera sensor width in mm");
  getopt_add_double(getopt, 's', "tag-size", "-1", "Tag size (black border, side of the square) in mm");
  getopt_add_bool(getopt, 'e', "estimate-focal-length", 0, "Do not track the marker; instead, estimate the camera focal length in pixels from the provided footage.");
  getopt_add_bool(getopt, 'q', "quick", 0, "Speed up the process at the expense of reduced accuracy.");

  if (!getopt_parse(getopt, argc, argv, 1) || getopt_get_bool(getopt, "help")) {
    printf("Usage: %s [options] <input files>\n", argv[0]);
    getopt_do_usage(getopt);
    exit(0);
  }

  const zarray_t *inputs = getopt_get_extra_args(getopt);
  const char *pathTemp = getopt_get_string(getopt, "path");
  std::string pathString=std::string(pathTemp);
  if(pathString.back()!='/' && pathString.back()!='\\')
  {
    pathString+="/";
  }
  const char* path=pathString.c_str();
  

  //printf(path);
  //zarray_get(inputs, 0, &path);

  apriltag_family_t *tf = NULL;
  const char *famname = "tag36h11";//getopt_get_string(getopt, "family");
  if (!strcmp(famname, "tag36h11")) {
    tf = tag36h11_create();
  } else if (!strcmp(famname, "tag25h9")) {
    tf = tag25h9_create();
  } else if (!strcmp(famname, "tag16h5")) {
    tf = tag16h5_create();
  } else if (!strcmp(famname, "tagCircle21h7")) {
    tf = tagCircle21h7_create();
  } else if (!strcmp(famname, "tagCircle49h12")) {
    tf = tagCircle49h12_create();
  } else if (!strcmp(famname, "tagStandard41h12")) {
    tf = tagStandard41h12_create();
  } else if (!strcmp(famname, "tagStandard52h13")) {
    tf = tagStandard52h13_create();
  } else if (!strcmp(famname, "tagCustom48h12")) {
    tf = tagCustom48h12_create();
  } else {
    printf("Unrecognized tag family name. Use e.g. \"tag36h11\".\n");
    exit(-1);
  }

  apriltag_detector_t *td = apriltag_detector_create();
  apriltag_detector_add_family_bits(td, tf, 1);//getopt_get_int(getopt, "hamming"));
  td->quad_decimate = 1;//getopt_get_double(getopt, "decimate");
  td->quad_sigma = 0;//getopt_get_double(getopt, "blur");
  td->nthreads = 1;//getopt_get_int(getopt, "threads");
  td->debug = 0;//getopt_get_bool(getopt, "debug");
  td->refine_edges = 1;//getopt_get_bool(getopt, "refine-edges");
  double tagSize=getopt_get_double(getopt, "tag-size");
  fmm=getopt_get_double(getopt,"focal-length-mm");
  fpx=getopt_get_double(getopt,"focal-length-pixels");
  sensorWidth=getopt_get_double(getopt,"sensor-width");
  quick = getopt_get_bool(getopt, "quick");

  bool estimateF = getopt_get_bool(getopt, "estimate-focal-length");
  if(!estimateF)
  {
    if(fmm==-1 && fpx==-1)
    {
      printf("ERROR. You need to specify either a pair of focal length in mm (--focal-length-mm) and the sensor width (--sensor-width), or specify the focal length in pixels (--focal-length-pixels).\r\n");
      printf("To estimate the focal length from your footage, first run apriltools.exe --path PATH_TO_FOOTAGE --estimate-focal-length, then rerun it on your desired footage, providing the newly obtained estimate using the parameter --focal-length-pixels");
      return 1;
    }

    if(fmm!=-1 &&sensorWidth==-1)
    {
      printf("ERROR. When providing the focal length in millimeters, you also need to provide the sensor width (using the --sensor-width parameter). You can also use this tool to estimate the focal length in pixels and then use that instead.");
      return 1;
    }

    if(fmm!=-1 && fpx!=-1 &&sensorWidth!=-1)
    {
      printf("ERROR: task was overconstrained. You may only specify 2 of the following at the same time: {focal-length-mm, focal-length-px, sensor-width}.");
    }

    if(sensorWidth==-1)
    {
      sensorWidth=36;
      printf("INFO: sensor width was not specified, defaulting to %.1f mm.", sensorWidth);
    }

    printf("\r\n");
    if(tagSize<0)
    {
      printf("WARNING: TAG SIZE UNSPECIFIED. DEFAULTING TO 173mm, but this is only true if your tag was printed at 100%% !!\r\n");
      tagSize=173;
    }
  }


  int quiet = 0;//getopt_get_bool(getopt, "quiet");

  int maxiters = 1;//getopt_get_int(getopt, "iters");

  const int hamm_hist_max = 10;




  if(strcmp(path,"")==0)
  {
    printf("No path specified. Use apriltools.exe --path (PATH_TO_FILE_SEQUENCE) to specify the input path. For more options, see apriltools --help .");
    return 1;
  }
  //char* path=R"(L:\personal\tracker\testAnim\)";
  FILE *fptr;
  if(!estimateF)
  {
  	fptr = fopen(((std::string)(path+(std::string)"track.txt")).c_str(),"w");
  }

  std::vector<std::string> filesList=getSortedFilenames(path);
  if(filesList.size()==0)
  {
    printf("No files found in the specified directory. The file names are expected to contain a number and end with .jpg, .jp2, .jpeg, .png or .bmp.\r\n");
    return 1;
  }
  printf("Found %d valid files in the specified directory.\r\n",filesList.size());


  int total_quads = 0;
  int total_hamm_hist[hamm_hist_max];
  memset(total_hamm_hist, 0, sizeof(total_hamm_hist));
  double total_time = 0;

  std::vector<double> focalLengths;

  int blurryPicsCount=0;
  int badPicsCount=0;
  for (int i=0;i<filesList.size();i++) {
    printf("Processed %d/%d files, out of which %d blurry (?) and %d were unusable (no tag found).\r\n",i+1,filesList.size(),blurryPicsCount,badPicsCount);

    const char* fileName=filesList[i].c_str();
    int frameNo=getFilenameNumber(std::string(fileName));
    int hamm_hist[hamm_hist_max];
    memset(hamm_hist, 0, sizeof(hamm_hist));

    //char path[100];
    //sprintf(path,"../../../../realFootage4/New Folder/%05d.jpg",i);
    //sprintf(path,"test_image.png",i);
    //zarray_get(inputs, input, &path);
    //printf("loading %s %s\n", path, fileName);

    image_u8_t *im = NULL;
    //printf("loading");

    cv::Mat img = cv::imread(path+(std::string)fileName, 0);//TODO: deallocate?
    image_u8_t img_header = { .width = img.cols,
      .height = img.rows,
      .stride = img.cols,
      .buf = img.data
    };

    im=&img_header;


    //return 1;
    if (img.data == NULL) {
      printf("couldn't load %s%s\r\n", path,fileName);
      continue;
    }

    //return 1;
    //TADY to uz nefacha...
    zarray_t *detections = detectTagsRobust(td, im, &blurryPicsCount, &badPicsCount);

    if(zarray_size(detections)==0)
    {
      printf("No tags detected on frame %s%s\r\n", path, fileName);
      continue;
    }

    //apriltag_detection_t* det;
    //zarray_get(detections,0,&det);


    /*
    float fx=2000;
    float fy=2000;
    float cx=960;
    float cy=540;
    float tagsize=97/1000.0f;
    */

    if(sensorWidthPx==-1)//not initialized yet
    {
      sensorWidthPx=img.cols;
      fillIntrinsics();
      if(!estimateF)
      {
        fprintf(fptr,"%s%s\r\n",path,fileName);
        fprintf(fptr,"%.4f, %.4f, %.4f, %d, %d,%d,%d,%d,%d\r\n",fmm,sensorWidth,tagSize, 0, 0,0,0,0,0);//padding for easier import
      }
    }

    float fx=fpx;//622.22;
    float fy=fpx;//622.22;
    float cx=img.cols/2.0;//320;
    float cy=img.rows/2.0;//240;


    //printf("%.2f\r\n",info.tagsize);
    //printf("%.2f\r\n",info.fx);

    double principal[2]={cx,cy};

    if(estimateF)
    {
      apriltag_detection_t* det;
      // base focal length off of first detection
      // TODO: get better estimate by looking at all detections?
      zarray_get(detections,0,&det);
      double pixelF=getPixelF(det, principal);
      if(pixelF>0)
      {
        focalLengths.push_back(pixelF);
      }
      double med=calculateMedian(focalLengths);
      printf("Focal length estimates median: %.2f; last estimate: %.2f\r\n",med, pixelF);

    }
    else
    {
      for (int i = 0; i < zarray_size(detections); i++) {
        apriltag_detection_t* det;
        zarray_get(detections, i, &det);

        apriltag_detection_info_t info;
        info.det = det;
        info.tagsize = tagSize/1000.0;
        info.fx = fx;
        info.fy = fy;
        info.cx = cx;
        info.cy = cy;

        // Then call estimate_tag_pose.
        apriltag_pose_t pose;
        double err = estimate_tag_pose(&info, &pose);
        matd_t* Mr =pose.R;
        matd_t* Mt =pose.t;

        matd_t* euler=getEulers(Mr);

        char line[1000];
        // save frame number, tag id, decision margin, tag pose to file
        sprintf(line,"%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\r\n", frameNo, det->id, det->decision_margin, MATD_EL(euler,0,2),MATD_EL(euler,0,1),MATD_EL(euler,0,0),MATD_EL(Mt,0,0),MATD_EL(Mt,1,0),MATD_EL(Mt,2,0));
        fprintf(fptr,"%s",line);
        matd_destroy(euler);
      }
    }




    apriltag_detections_destroy(detections);

    //if (!quiet) {
      //timeprofile_display(td->tp);
    //}

    total_quads += td->nquads;

    /*
    if (!quiet)
      printf("hamm ");

    for (int i = 0; i < hamm_hist_max; i++)
      printf("%5d ", hamm_hist[i]);
      */

    double t =  timeprofile_total_utime(td->tp) / 1.0E3;
    total_time += t;
    //printf("%12.3f ", t);
    //printf("%5d", td->nquads);

    //printf("\n");

    //image_u8_destroy(im);//TODO: DO NOT DO THIS
  }

  printf("\r\n");
  if(!estimateF)
  {
    fclose(fptr);
    printIntrinsics();
    printf("Camera track saved to: %s\r\n",(path+(std::string)"track.txt").c_str());
    printf("Open Blender and go to file-> import -> Apriltools tracking data (install the plugin if you haven't already) to import the tracking data.\r\n");
  }
  else
  {
    //something else
     double med=calculateMedian(focalLengths);
     printf("Focal length estimation complete. Best guess: %.2f px. Use apriltools --path PATH_TO_FOOTAGE --focal-length-pixels %.2f --tag-size TAG_SQUARE_SIZE_IN_MILLIMETERS to track your footage now.\r\n", med, med);
  }

  /**
  printf("Summary\n");

  printf("hamm ");

  for (int i = 0; i < hamm_hist_max; i++)
    printf("%5d ", total_hamm_hist[i]);
  printf("%12.3f ", total_time);
  printf("%5d", total_quads);
  printf("\n");
  */


  // don't deallocate contents of inputs; those are the argv
  apriltag_detector_destroy(td);

  if (!strcmp(famname, "tag36h11")) {
    tag36h11_destroy(tf);
  } else if (!strcmp(famname, "tag25h9")) {
    tag25h9_destroy(tf);
  } else if (!strcmp(famname, "tag16h5")) {
    tag16h5_destroy(tf);
  } else if (!strcmp(famname, "tagCircle21h7")) {
    tagCircle21h7_destroy(tf);
  } else if (!strcmp(famname, "tagCircle49h12")) {
    tagCircle49h12_destroy(tf);
  } else if (!strcmp(famname, "tagStandard41h12")) {
    tagStandard41h12_destroy(tf);
  } else if (!strcmp(famname, "tagStandard52h13")) {
    tagStandard52h13_destroy(tf);
  } else if (!strcmp(famname, "tagCustom48h12")) {
    tagCustom48h12_destroy(tf);
  }

  getopt_destroy(getopt);

  return 0;
}
