# AprilTools
A marker-based camera tracker for Blender based around the AprilTags Library by University of Michigan

## Usage
Download everything (Clone or download -> Download ZIP); unpack.
Print out one of the pages of the pdf files in the ./bin/ directory;
Shoot the footage with this printed out page in frame (you can cut out just the tag along the dotted line)
Convert this footage to a sequence of .png or .jpg files (you can use Blender for this)

  Install the .py plugin located in ./bin/ into your Blender 2.80 (Inside Blender: Edit/Preferences/Addons/Install...) then go in the bin folder, right click on the folder (not on any file, just on empty space), select "Run command prompt here". Inside the command prompt, you can run (by typing it in and then pressing enter):

apriltools.exe --path "PATH_TO_THE_FOLDER_WITH_THE_IMAGE_SEQUENCE_INSIDE" --focal-length-mm FOCAL_LENGTH_IN_MILLIMETERS --sensor-width SENSOR_WIDTH_IN_MILLIMETERS

That is, for example

apriltools.exe --path "C:/myFootage/" --focal-length-mm 35.00 --sensor-width 36.00

If you don't know the focal length of your camera in mm and the sensor width, you can instead provide the focal length in pixels by running

apriltools.exe --path "PATH_TO_THE_FOLDER_WITH_THE_IMAGE_SEQUENCE_INSIDE" --focal-length-pixels FOCAL_LENGTH_IN_PIXELS

That is, for example

apriltools.exe --path "C:/myFootage/" --focal-length-pixels 2000

This application can determine the focal length in pixels automatically from the footage. You only have to do this once and then you can reuse the result. To estimate the focal length, run:


apriltools.exe --path "PATH_TO_THE_FOLDER_WITH_THE_IMAGE_SEQUENCE_INSIDE" --estimate-focal-length

That is, for example

apriltools.exe --path "C:/myFootage/"  --estimate-focal-length

This will NOT track the footage, but instead it will calculate the focal length in pixels, allowing you to use this info to track your footage.

Having installed the plugin into Blender, you can go to File->Import->Import AprilTools tracking data (inside Blender).
You then need to choose whether you want the tracking data to be interpreted as a tag being static and the camera orbiting around it (default), or whether you want the camera to be static and the tag to move.

The tracking data will be located in the folder your footage is in. 
Optionally, you can now right-click the camera and on the right Blender panel go to Context->Background Images; there you click on Add Image and you select your footage.

The importer will generate a sequence of keyframes for the camera or the tag, based on what you selected as a target. In some cases, this sequence can be off by a single frame with respect to the actual footage. To fix this, select the keyframes and hit "g" to move them around, until you see the marker Blender object match up with the tag on the video.

## Donations

This software is and always will be free, but you can support the development here:

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=QHXE57NKFC3QL)

ETH: 0x2E4d8B0319D3e1E6428aC9B7dDf862D376A99e21

BTC: 1DTzg2yxkgudJNwYmGrbg6T9DPLNbDMkbq
