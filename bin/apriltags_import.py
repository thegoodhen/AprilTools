import bpy
import numpy as np
from mathutils import Euler
from mathutils import Matrix
bl_info = {"name": "Import AprilTools tracking data", "category": "Import-Export", "blender":(2,80,0)}

def read_some_data(context, filepath, type):
    print("Importing Apriltools data..")
    bpy.context.scene.unit_settings.system='METRIC'
    #f = open(filepath, 'r', encoding='utf-8')
    #data = f.read()
    #f.close()
    #track=np.genfromtxt(filepath, delimiter=',')
    print(type)
    track=np.genfromtxt(filepath, delimiter=',')
        
        
    cam=bpy.data.objects.get("Camera")
    marker=bpy.data.objects.get("marker")

    if marker is not None:#marker already exists
        bpy.ops.object.delete({"selected_objects":[marker]})
        
    if cam is not None:#camera already exists
        bpy.ops.object.delete({"selected_objects":[cam]})


    bpy.ops.mesh.primitive_plane_add(size=1,enter_editmode=False,location=(0,0,0))
    bpy.context.active_object.name = 'marker'
    marker=bpy.data.objects['marker']

    bpy.ops.object.camera_add(enter_editmode=False,location=(0,0,0),rotation=(0,0,0))
    bpy.context.active_object.name = 'Camera'
    cam=bpy.data.objects['Camera']
        

  

    cf=0
    firstLine=True
    fmm=-1
    sensor_width=-1
    tagSize=-1

    for line in track:
        print(line)
        if firstLine:
            fmm=line[0]
            sensor_width=line[1]
            tagSize=line[2]
            firstLine=False
            continue
                
        
        cf=line[0]
        the_rot_wrong=line[1:4]
        the_rot=np.copy(the_rot_wrong)
        #print("slepice")
        #print(the_rot)
        the_rot[0]=the_rot_wrong[2]
        the_rot[2]=the_rot_wrong[0]
        
        print(the_rot)
        #the_rot[0]=-the_rot[0]
        the_rot[1]=-the_rot[1]
        the_rot[2]=-the_rot[2]
        
        the_loc=line[4:7]
        the_loc[1]=-the_loc[1]
        the_loc[2]=-the_loc[2]
        
        Rx=Matrix.Rotation(the_rot[0],4,'X')
        Ry=Matrix.Rotation(the_rot[1],4,'Y')
        Rz=Matrix.Rotation(the_rot[2],4,'Z')
        Txyz=Matrix.Translation(tuple(the_loc))
        
        Mt=(Txyz@Rz@Ry@Rx)
        if (type == 'MOVING_CAM'):
            obj=cam
            Mt.invert()
        else: 
            obj=marker
            
        obj.matrix_world=Mt
            
        #obj.rotation_euler=Euler(tuple(the_rot),'XYZ')
        #obj.location=tuple(the_loc)
        obj.keyframe_insert(data_path='rotation_euler',frame=cf)
        obj.keyframe_insert(data_path='location',frame=cf)
        
    cam=bpy.data.objects['Camera']
    marker=bpy.data.objects['marker']
    cam.data.lens=fmm
    cam.data.sensor_width=sensor_width
    marker.dimensions=[tagSize/1000,tagSize/1000,0.01]

    return {'FINISHED'}


# ImportHelper is a helper class, defines filename and
# invoke() function which calls the file selector.
from bpy_extras.io_utils import ImportHelper
from bpy.props import StringProperty, BoolProperty, EnumProperty
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

    type: EnumProperty(
        name="Movement target",
        description="Whether the camera should be static and the tag moving, or vice versa.",
        items=(
            ('MOVING_CAM', "Moving camera", "Tag is static, camera is moving around it."),
            ('MOVING_TAG', "Moving tag", "Camera is static, tag is moving."),
        ),
        default='MOVING_CAM',
    )

    def execute(self, context):
        return read_some_data(context, self.filepath, self.type)


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
