import bpy
import numpy as np
from mathutils import Euler
from mathutils import Matrix
bl_info = {"name": "Import AprilTools tracking data", "category": "Import-Export", "blender":(2,80,0)}



def loadSequence(theCam, thePath, frameCount):
    seq = bpy.data.images.load(thePath)
    seq.source = 'SEQUENCE'
    theCam.data.show_background_images = True
    bpy.context.scene.frame_end=frameCount
    bg = theCam.data.background_images.new()
    bg.image = seq
    setTrackLen(theCam,frameCount)
    
    
def setTrackLen(theCam, frameCount):
    bg = theCam.data.background_images[0]
    #seq.frame_end=frameCount
    bg.image_user.frame_duration=frameCount
    bpy.context.scene.frame_end=frameCount


def read_some_data(context, filepath, fixed_tag):
    print("Importing Apriltools data..")
    bpy.context.scene.unit_settings.system='METRIC'
   
    f=open(filepath)
    firstImage=f.readline()
    firstImage=firstImage.rstrip()
    print(firstImage)
    f.close();
    
    track=np.genfromtxt(filepath, delimiter=',',skip_header=1)
        
        
    cam=bpy.data.objects.get("Camera")

    if cam is not None:#camera already exists
        bpy.ops.object.delete({"selected_objects":[cam]})

    # delete objects starting with "apriltag" to get rid of old tags
    # hopefully this isn't too slow with a lot of objects...
    for o in bpy.data.objects:
        if o.name[0:8] == "apriltag":
            bpy.ops.object.delete({"selected_objects": [o]})

    cf=0
    firstLine=True
    fmm=-1
    sensor_width=-1
    tagSize=-1

    header = track[0]
    fmm = header[0]
    sensor_width = header[1]
    tagSize = header[2]

    bpy.ops.object.camera_add(enter_editmode=False,location=(0,0,0),rotation=(0,0,0))
    bpy.context.active_object.name = 'Camera'
    cam=bpy.data.objects['Camera']
    cam.data.lens=fmm
    cam.data.sensor_width=sensor_width

    loadSequence(cam,firstImage,int(track[len(track)-1][0]))
    # get rid of header line
    track = track[1:]

    # store the transforms of each tag in the frame so you can adjust for the camera location
    tags_in_frame = {}
    # keep track of the frame number from last loop
    last_frame = int(track[0][0])

    for line in track:
        frame = line[0]
        id = int(line[1]) # the tag id

        # it's the next frame, so actually set the positions
        # we need to wait for the next frame before setting the positions since the position of the tags relies on the camera
        if frame != last_frame:
            for k in tags_in_frame.keys():
                # the object to move. I start it as the camera since it's possible the marker doesn't exist yet
                obj = cam
                if k == fixed_tag:
                    # If this is the fixed tag, invert the transformation and move the camera instead of the tag
                    tags_in_frame[k].invert()
                else:
                    # create a tag plane if it doesn't exist yet
                    if 'apriltag_' + str(k) not in bpy.data.objects:
                        bpy.ops.mesh.primitive_plane_add(size=tagSize/1000,enter_editmode=False,location=(0,0,0))
                        bpy.context.active_object.name = 'apriltag_' + str(k)
                        marker = bpy.data.objects['apriltag_' + str(k)]
                    # set the object to the current tag
                    obj = bpy.data.objects['apriltag_' + str(k)]

                    # if there is a fixed tag, the camera must be moving
                    # apply the camera transform to the tags so they are in the correct position relative to the camera
                    if fixed_tag != -1:
                        tags_in_frame[k] = tags_in_frame[fixed_tag] @ tags_in_frame[k]

                obj.matrix_world = tags_in_frame[k]

                obj.keyframe_insert(data_path='rotation_euler',frame=frame-1)
                obj.keyframe_insert(data_path='location',frame=frame-1)

            last_frame = frame
            tag_in_frame = {}

        raw_rotation = line[3:6]
        rotation = np.copy(raw_rotation)

        rotation[0] = raw_rotation[2]
        rotation[2] = raw_rotation[0]
        rotation[1] *= -1
        rotation[2] *= -1

        location = line[6:9]
        location[1] *= -1
        location[2] *= -1

        Rx = Matrix.Rotation(rotation[0], 4, 'X')
        Ry = Matrix.Rotation(rotation[1], 4, 'Y')
        Rz = Matrix.Rotation(rotation[2], 4, 'Z')
        translation = Matrix.Translation(tuple(location)) # the location transformation
        tags_in_frame[id] = (translation@Rz@Ry@Rx) # add the rotation to make the full transform


    return {'FINISHED'}


# ImportHelper is a helper class, defines filename and
# invoke() function which calls the file selector.
from bpy_extras.io_utils import ImportHelper
from bpy.props import StringProperty, BoolProperty, EnumProperty, IntProperty
from bpy.types import Operator


class ImportMovement(Operator, ImportHelper):
    """This appears in the tooltip of the operator and in the generated docs"""
    bl_idname = "import_test.some_data"  # important since its how bpy.ops.import_test.some_data is constructed
    bl_label = "Import AprilTools camera tracking data"

    # ImportHelper mixin class uses this
    filename_ext = ".txt"

    filter_glob: StringProperty(
        default="*.txt",
        options={'HIDDEN'},
        maxlen=255,  # Max internal buffer length, longer would be clamped.
    )

    # List of operator properties, the attributes will be assigned
    # to the class instance from the operator settings before calling.
    #use_setting: BoolProperty(
    #    name="Example Boolean",
    #    description="Example Tooltip",
    #    default=True,
    #)

    # It should be possible to have multiple fixed tags, but there would need to be a way to combine them into one camera position
    # In order to do this you need to know the difference between the positions of the two (or more) tags
    # This could be determined from a frame of the video (the first one where both tags are visible?)
    # Or it could be measured and inputted manually by the user
    # For now, I'm going to leave it with a single fixed tag since that is a likely use case
    fixed_tag: IntProperty(
        name="Fixed tag",
        description="""The id of a tag that is in a fixed position, used to find the camera position
Set to -1 to assume the camera is fixed""",
        default=-1,
    )

    def execute(self, context):
        return read_some_data(context, self.filepath, self.fixed_tag)


# Only needed if you want to add into a dynamic menu
def menu_func_import(self, context):
    self.layout.operator(ImportMovement.bl_idname, text="Apriltools tracking data")


def register():
    bpy.utils.register_class(ImportMovement)
    bpy.types.TOPBAR_MT_file_import.append(menu_func_import)


def unregister():
    bpy.utils.unregister_class(ImportMovement)
    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import)


if __name__ == "__main__":
    register()

    # test call
    bpy.ops.import_test.some_data('INVOKE_DEFAULT')
