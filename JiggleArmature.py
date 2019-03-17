#Copyright (c) 2019 Simón Flores (https://github.com/cheece)

#Permission is hereby granted, free of charge, 
#to any person obtaining a copy of this software 
#and associated documentation files (the "Software"), 
#to deal in the Software without restriction, 
#including without limitation the rights to use, 
#copy, modify, merge, publish, distribute, sublicense,
#and/or sell copies of the Software, and to permit 
#persons to whom the Software is furnished to do so,
#subject to the following conditions:The above copyright to_quaternion
#notice and this permission notice shall be included 
#in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY 
#OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
#LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
#FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
#EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
#FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
#AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
#OTHER DEALINGS IN THE SOFTWARE.

#this a jiggle bone animation tool

#to enable jiggle physics first enable "jiggle scene" in the scene properties and then enable jiggle bone on the bones
 
bl_info = {
    "name": "Jiggle Armature",
    "author": "Simón Flores",
    "version": (2, 1, 0),
    "blender": (2, 77, 0),
    "description": "Jiggle bone animation tool",
    "warning": "",
    "wiki_url": "",
    "category": "Animation",
}


import bpy
from bpy.types import Menu, Panel, UIList
from rna_prop_ui import PropertyPanel
import os
import bmesh
from bpy_extras.io_utils import ExportHelper
from bpy_extras.io_utils import ImportHelper
import mathutils
import math
import os
from mathutils import *
from bpy.app.handlers import persistent
from bpy_extras.io_utils import ExportHelper
from bpy.props import StringProperty, BoolProperty, EnumProperty
from bpy.types import Operator


FIX_DISPLACEMENT=False
SHOW_DISPLACEMENT = True

    
class JiggleArmature(bpy.types.PropertyGroup):
#    enabled = bpy.props.BoolProperty(default=False, update = funp("enabled"))
    fps = bpy.props.FloatProperty(name = "simulation fps",default = 24)
    time_acc = bpy.props.FloatProperty(default= 0.0)
    #Kd = bpy.props.FloatProperty(name = "angular damping",min=0.0, max=1.0,default = 0.01, update = funp("Kd"))
   # Ks = bpy.props.FloatProperty(name = "stiffness",min=0.0 , max = 1.0, default = 0.8, update = funp("Ks"))
    
class JiggleScene(bpy.types.PropertyGroup):
    test_mode = bpy.props.BoolProperty(default=False)
    sub_steps = bpy.props.IntProperty(min=1, default = 2)
    iterations = bpy.props.IntProperty(min=1, default = 4)
    fix_iterations = bpy.props.IntProperty(min=1, default = 0)
    last_frame = bpy.props.IntProperty()
    length_fix_iters = bpy.props.IntProperty(min=0, default = 2)
    show_displacement = bpy.props.BoolProperty( default = False)
    fix_displacement= bpy.props.BoolProperty( default = True)
def writemat(f, M):
    M = M.transposed()
    for i in M:
        for j in i:
            f.write(str(j)+" ")
    f.write("\n")
 
class JiggleArmaturePanel(bpy.types.Panel):
    bl_idname = "Armature_PT_jiggle"
    bl_label = "Jiggle Armature"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "data"
    bl_options = {'DEFAULT_CLOSED'}
    @classmethod
    def poll(cls, context):
        return (context.object is not None and context.object.type == 'ARMATURE')

    def draw(self, context):
        layout = self.layout
        col = layout.column() 
        col.prop(context.object.data.jiggle,"fps")
        

class JiggleScenePanel(bpy.types.Panel):
    bl_idname = "Scene_PT_jiggle"
    bl_label = "Jiggle Scene"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "scene"
    bl_options = {'DEFAULT_CLOSED'}

    def draw_header(self, context):
        layout = self.layout
        bon = context.bone
        layout.prop(context.scene.jiggle,"test_mode", text="")


    def draw(self, context):
        layout = self.layout
        col = layout.column()
        #col.prop(context.scene.jiggle,"test_mode")
        #col.prop(context.scene.jiggle,"sub_steps")
        col.prop(context.scene.jiggle,"iterations")
        col.prop(context.scene.jiggle,"fix_iterations")
       # col.prop(context.scene.jiggle,"show_displacement")
       # col.prop(context.scene.jiggle,"fix_displacement") 
        col.operator("jiggle.bake", text="Bake Selected").a = False
        col.operator("jiggle.bake", text="Bake All").a = True



def funKd(self,context):
    b = context.bone
    o = context.object
    arm = o.data
    for b2 in arm.bones:
        b2.jiggle.Kd = b.jiggle.Kd

def funKs(self,context):
    b = context.bone
    o = context.object
    arm = o.data
    for b2 in arm.bones:
        b2.jiggle.Ks = b.jiggle.Ks
        
def funmass(self,context):
    b = context.bone
    o = context.object
    arm = o.data
    for b2 in arm.bones:
        b2.jiggle.mass = b.jiggle.mass
        
inop = False        
def funp(prop):
    def f(self,context):
        global inop
        if(inop):
            return 
        inop = True
        b = context.bone
        o = context.object
        arm = o.data
        for b2 in arm.bones:
            if(b2.select):
                setattr(b2.jiggle, prop, getattr(b.jiggle,prop)) 
        inop = False
    return f
class JiggleBone(bpy.types.PropertyGroup):
    enabled = bpy.props.BoolProperty(default=False, update = funp("enabled"))
    Kld = bpy.props.FloatProperty(name = "linear damping",min=0.0, max=1.0,default = 0.01, update = funp("Kld"))
    Kd = bpy.props.FloatProperty(name = "angular damping",min=0.0, max=1.0,default = 0.01, update = funp("Kd"))
    Ks = bpy.props.FloatProperty(name = "stiffness",min=0.0 , max = 1.0, default = 0.8, update = funp("Ks"))
    mass = bpy.props.FloatProperty(min=0.0001, default = 1.0, update = funp("mass"))  
    #M = bpy.props.FloatVectorProperty(size=9,subtype='MATRIX')    
    R = bpy.props.FloatVectorProperty(name="rotation", size=4,subtype='QUATERNION')
    W = bpy.props.FloatVectorProperty(size=3,subtype='XYZ') #angular velocity
    P = bpy.props.FloatVectorProperty(size=3,subtype='XYZ')
    V = bpy.props.FloatVectorProperty(size=3,subtype='XYZ')	#linear velocity, ok? 
    
    use_custom_rest = bpy.props.BoolProperty(default=False, update = funp("use_custom_rest"))
    
    rest = bpy.props.FloatVectorProperty(name="rotation", size=4,subtype='QUATERNION')
    
    control = bpy.props.FloatProperty(name = "control",min=0.0, max=1.0,default = 1, update = funp("control"))
    control_bone = bpy.props.StringProperty(name = "control bone")
    debug = bpy.props.StringProperty()
    
def skew(v):
    m = Matrix.Identity(3)
    m[0][0], m[1][0], m[2][0] =     0, -v[2],  v[1]
    m[0][1], m[1][1], m[2][1] =  v[2],     0, -v[0]
    m[0][2], m[1][2], m[2][2] = -v[1], v[0],     0
    return m #m.transposed()
def iskew(m):
    v = Vector((m[1][2] - m[2][1], m[2][0]- m[0][2], m[0][1]-m[1][0]))*0.5
    return v #m.transposed()
def setm(om, m):
    for i in range(3):
        for j in range(3):
            om[i][j] = m[i][j]
           #posed()
def setq(om, m):
    for i in range(4):
        om[i]= m[i]
def getm(m, om):
    for i in range(3):
        for j in range(3):
            m[i][j]=om[i][j]
           
class ResetJigglePropsOperator(bpy.types.Operator):
    bl_idname = "jiggle.reset"
    bl_label = "Reset State"
    def execute(self, context):
        scene = context.scene
        for o in scene.objects:
            if(o.select and o.type == 'ARMATURE' ):
                arm = o.data
                ow = o.matrix_world
                scale = maxis(ow,0).length
                iow = ow.inverted()
                i=0
                for b in o.pose.bones:
                    if(b.bone.select):                    
                        M = ow*b.matrix #ow*Sbp.wmat* Sb.rmat #im
                        #l ,r,s = M.decompose()
                        
                        Jb = b.bone.jiggle
                        setq(Jb.R, M.to_quaternion().normalized())
                        Jb.V = Vector((0,0,0))
                        Jb.P =  mpos(M)+maxis(M,1)*b.bone.length*0.5
                        Jb.W = Vector((0,0,0))
                        #Jb.M = Matrix(M)
        return {'FINISHED'}
class SetRestJigglePropsOperator(bpy.types.Operator):
    bl_idname = "jiggle.set_rest"
    bl_label = "Set Rest"
    def execute(self, context):
        scene = context.scene
        for o in scene.objects:
            if(o.select and o.type == 'ARMATURE' ):
                arm = o.data
                ow = o.matrix_world
                scale = maxis(ow,0).length
                iow = ow.inverted()
                i=0
                for b in o.pose.bones:
                    if(b.bone.select):                    
                        M = b.parent.matrix.inverted()*b.matrix #ow*Sbp.wmat* Sb.rmat #im 
                        Jb = b.bone.jiggle
                        setq(Jb.rest, M.to_quaternion().normalized()) 
                        Jb.use_custom_rest = True
        return {'FINISHED'}

class JiggleBonePanel(bpy.types.Panel):
    bl_idname = "Bone_PT_jiggle_bone"
    bl_label = "Jiggle Bone"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "bone"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return ( context.bone is not None and context.object is not None and context.object.type == 'ARMATURE')

    def draw_header(self, context):
        layout = self.layout
        bon = context.bone
        layout.prop(bon.jiggle, "enabled", text="")

    def draw(self, context):
        layout = self.layout

        bon = context.bone
        col = layout.column()
        layout.enabled = context.scene.jiggle.test_mode
        if(not context.scene.jiggle.test_mode):
            col.label(text= "jigglescene disabled, see scene properties")   
            
        if(bon.jiggle.enabled):
            col.prop(bon.jiggle,"Ks")   
            col.prop(bon.jiggle,"Kd")  
            col.prop(bon.jiggle,"Kld")             
            col.prop(bon.jiggle,"mass")  
            col.prop(bon.jiggle,"debug")
            col.prop(bon.jiggle,"control")
            col.prop(bon.jiggle,"control_bone")
            
            col.prop(bon,"jiggle_tp")    
            
            col.operator("jiggle.reset")
            col.prop(bon.jiggle,"use_custom_rest") 
            if(bon.jiggle.use_custom_rest):
                col.prop(bon.jiggle,"rest")
            col.operator("jiggle.set_rest")
            if(bon.parent==None):                
                col.label(text= "warning: jibblebones without parent will fall",icon='COLOR_RED')   

def centerM(wb,l):
    ax = maxis(wb,1).normalized()
    wb[0][3] += ax[0]*l*0.5
    wb[1][3] += ax[1]*l*0.5
    wb[2][3] += ax[2]*l*0.5


#taken from https://github.com/InteractiveComputerGraphics/PositionBasedDynamics/blob/master/PositionBasedDynamics/PositionBasedRigidBodyDynamics.cpp
def computeMatrixK(connector,invMass,x,inertiaInverseW,K):
    if (invMass != 0.0):
        v = connector - x
        a = v[0]
        b = v[1]
        c = v[2]
        if(True):
            j11 = inertiaInverseW[0][0]
            j12 = inertiaInverseW[1][0]
            j13 = inertiaInverseW[2][0]
            j22 = inertiaInverseW[1][1]
            j23 = inertiaInverseW[2][1] 
            j33 = inertiaInverseW[2][2]

            K[0][0] = c*c*j22 - b*c*(j23 + j23) + b*b*j33 + invMass
            K[1][0] = -(c*c*j12) + a*c*j23 + b*c*j13 - a*b*j33
            K[2][0] = b*c*j12 - a*c*j22 - b*b*j13 + a*b*j23
            K[0][1] = K[1][0]
            K[1][1] = c*c*j11 - a*c*(j13 + j13) + a*a*j33 + invMass
            K[2][1] = -(b*c*j11) + a*c*j12 + a*b*j13 - a*a*j23
            K[0][2] = K[2][0]
            K[1][2] = K[2][1]
            K[2][2] = b*b*j11 - a*b*(j12 + j12) + a*a*j22 + invMass
        else:
            j11 = inertiaInverseW[0][0]
            j12 = inertiaInverseW[0][1]
            j13 = inertiaInverseW[0][2]
            j22 = inertiaInverseW[1][1]
            j23 = inertiaInverseW[1][2]
            j33 = inertiaInverseW[2][2]

            K[0][0] = c*c*j22 - b*c*(j23 + j23) + b*b*j33 + invMass
            K[0][1] = -(c*c*j12) + a*c*j23 + b*c*j13 - a*b*j33
            K[0][2] = b*c*j12 - a*c*j22 - b*b*j13 + a*b*j23
            K[1][0] = K[0][1]
            K[1][1] = c*c*j11 - a*c*(j13 + j13) + a*a*j33 + invMass
            K[1][2] = -(b*c*j11) + a*c*j12 + a*b*j13 - a*a*j23
            K[2][0] = K[0][2]
            K[2][1] = K[1][2]
            K[2][2] = b*b*j11 - a*b*(j12 + j12) + a*a*j22 + invMass
    else:
        K.zero()

class JB:
    def __init__(self, b,M,p):
        self.M  = M.copy()
        self.length = b.bone.length*maxis(M,0).length 
        self.b = b
        self.parent = p
        self.rest = None
        self.rest_w = None
        self.w = 0
        self.X = None
        self.P = None
        self.R = None
        self.Q = None
        self.cQ = None
        self.Kc = 0
        self.iI = Matrix.Identity(3) #first naive approach
        self.iIw = self.iI
    def computeI(self):
        self.iI = Matrix.Identity(3)*(self.w/(self.l*self.l)*5.0/2.0)
    def updateIW(self):
        rot = self.Q.to_matrix()
        self.iIw = rot*self.iI*rot.transposed()
 
def propB(ow,b, l, p):
    j = JB(b, ow*b.matrix, p)
    l.append(j)
    for c in b.children:
        propB(ow,c,l,j)
def maxis(M,i):
    return Vector((M[0][i],M[1][i],M[2][i]))    
def saxis(M,i,v):
    M[0][i] = v[0]
    M[1][i] = v[1]
    M[2][i] = v[2]        
def mpos(M):
    return Vector((M[0][3],M[1][3],M[2][3]))   
def ort(M):
    a = M[0]
    b = M[1]
    c = M[2]
    a = a.normalized()
    b = (b- a*a.dot(b)).normalized()
    c = (c - a*a.dot(c) - b*b.dot(c)).normalized()
    M = Matrix.Identity(3)
    M[0] = a
    M[1] = b
    M[2] = c
    return M 
def qadd(a,b):
    return Quaternion((a[0]+b[0],a[1]+b[1],a[2]+b[2],a[3]+b[3]))
def qadd2(a,b):
    a.x+=b.x
    a.y+=b.y
    a.z+=b.z
    a.w+=b.w
    
def normR(m):
    for i in range(3):
        saxis(m,i, maxis(m,i).normalized())

def locSpringGradient1(Q0,Q1,P0,P1,v0,lf):
    Q0x = Q0.x
    Q0y = Q0.y
    Q0z = Q0.z
    Q0w = Q0.w
    Q1x = Q1.x
    Q1y = Q1.y
    Q1z = Q1.z
    Q1w = Q1.w  
    v0x = v0.x
    v0y = v0.y
    v0z = v0.z 
    P0x = P0.x
    P0y = P0.y
    P0z = P0.z
    P1x = P1.x
    P1y = P1.y
    P1z = P1.z    

        
    tmp0 = -16.0*Q0w
    tmp1 = tmp0*Q0x
    tmp2 = tmp1*Q0z
    tmp3 = tmp1*Q0y
    tmp4 = tmp3*v0x
    tmp5 = tmp2*v0x
    tmp6 = tmp0*Q0y*Q0z
    tmp7 = -16.0*Q0x*Q0y*Q0z
    tmp8 = tmp6*v0x
    tmp9 = -8.0*Q0x*Q1w
    tmp10 = -8.0*Q0z
    tmp11 = -8.0*Q0w*Q1w
    tmp12 = tmp7*v0x
    tmp13 = tmp11*Q1z*lf
    tmp14 = -8.0*Q0y*Q1w*Q1x*lf
    tmp15 = tmp9*Q1z*lf
    tmp16 = 16.0*Q0w*Q0x*Q0z
    tmp17 = -8.0*Q0y*Q1w
    tmp18 = tmp17*Q1z*lf
    tmp19 = -8.0*Q0y*Q1y*Q1z*lf
    tmp20 = 16.0*Q0w*Q0y*Q0z*v0x
    tmp21 = 16.0*Q0w*Q0x
    tmp22 = 8.0*Q0w*Q1w
    tmp23 = tmp10*Q1x*Q1y*lf
    tmp24 = 16.0*Q0x*Q0y*Q0z*v0x
    tmp25 = tmp21*Q0y*v0x
    tmp26 = 16.0*Q0w
    tmp27 = tmp26*Q0y*Q0z
    tmp28 = tmp21*Q0y
    tmp29 = tmp16*v0x
    tmp30 = 16.0*Q0x*Q0y*Q0z
    tmp31 = 8.0*Q0x*Q1x*Q1y*lf
    tmp32 = tmp22*Q1x*lf
    tmp33 = 8.0*Q0y*Q1w
    tmp34 = 8.0*Q0w*Q1x*Q1y*lf
    tmp35 = 8.0*Q0w*Q1y*Q1z*lf
    tmp36 = 8.0*Q0z
    tmp37 = tmp36*Q1w
    tmp38 = 8.0*Q0x*Q1w*Q1x*lf
    tmp39 = 8.0*Q0x*Q1y*Q1z*lf
    tmp40 = tmp37*Q1x*lf
    tmp41 = 8.0*Q0y*Q1x*Q1y*lf
    tmp42 = tmp37*Q1z*lf
    tmp43 = tmp36*Q1y*Q1z*lf
    tmp44 = -4.0*Q0w*Q0w*Q0x
    tmp45 = -4.0*Q0w*Q0w*Q0y
    tmp46 = -4.0*Q0x*Q0y*Q0y
    tmp47 = -4.0*Q0w*Q0w*Q0z
    tmp48 = -4.0*Q0x*Q0z*Q0z
    tmp49 = -4.0*Q0w*Q1x*Q1x*lf
    tmp50 = -4.0*Q0y*Q0z*Q0z
    tmp51 = -4.0*Q0w*Q1z*Q1z*lf
    tmp52 = -4.0*Q0x*Q1w*Q1w*lf
    tmp53 = -4.0*Q0y*Q1x*Q1x*lf
    tmp54 = -4.0*Q0y*Q1z*Q1z*lf
    tmp55 = -4.0*Q0z
    tmp56 = -4.0*Q0x*Q1y*Q1y*lf
    tmp57 = -8.0*Q0w*Q0w
    tmp58 = tmp55*Q1x*Q1x*lf
    tmp59 = tmp57*Q0x
    tmp60 = tmp55*Q1z*Q1z*lf
    tmp61 = tmp57*Q0y
    tmp62 = tmp57*Q0z
    tmp63 = tmp59*v0x
    tmp64 = -8.0*Q0x*Q0x*Q0y
    tmp65 = -8.0*Q0x*Q0x*Q0z
    tmp66 = tmp61*v0x
    tmp67 = -8.0*Q0x*Q0x
    tmp68 = tmp62*v0x
    tmp69 = -8.0*Q0y*Q0y*Q0z
    tmp70 = -8.0*Q0w*Q0x*Q0x*v0x
    tmp71 = -8.0*Q0w*Q0x*Q0x
    tmp72 = tmp64*v0x
    tmp73 = -8.0*Q0w*Q0y*Q0y
    tmp74 = tmp73*v0x
    tmp75 = tmp65*v0x
    tmp76 = -8.0*Q0w*Q0z*Q0z
    tmp77 = tmp76*v0x
    tmp78 = -8.0*Q0x*Q0y*Q0y
    tmp79 = tmp69*v0x
    tmp80 = -8.0*Q0x*Q0z*Q0z
    tmp81 = -8.0*Q0w*Q0x*Q0y
    tmp82 = -8.0*Q0w*Q0x*Q0z
    tmp83 = tmp78*v0x
    tmp84 = -8.0*Q0y*Q0z*Q0z
    tmp85 = -8.0*Q0w*Q0y*Q0z
    tmp86 = tmp80*v0x
    tmp87 = -8.0*Q0x*Q0y*Q0z
    tmp88 = 4.0*Q0w*Q0w*Q0x
    tmp89 = tmp84*v0x
    tmp90 = 4.0*Q0w*Q0w*Q0y
    tmp91 = 4.0*Q0w*Q0w*Q0z
    tmp92 = 4.0*Q0x*Q0y*Q0y
    tmp93 = 4.0*Q0x*Q0z*Q0z
    tmp94 = 4.0*Q0w*Q1w*Q1w*lf
    tmp95 = 4.0*Q0y*Q0z*Q0z
    tmp96 = 4.0*Q0w*Q1y*Q1y*lf
    tmp97 = 4.0*Q0y*Q1w*Q1w*lf
    tmp98 = 4.0*Q0x*Q1x*Q1x*lf
    tmp99 = 4.0*Q0z
    tmp100 = 4.0*Q0x*Q1z*Q1z*lf
    tmp101 = 4.0*Q0y*Q1y*Q1y*lf
    tmp102 = tmp99*Q1w*Q1w*lf
    tmp103 = 8.0*Q0w*Q0w*Q0x
    tmp104 = tmp99*Q1y*Q1y*lf
    tmp105 = tmp103*v0x
    tmp106 = 8.0*Q0w*Q0w*Q0y
    tmp107 = 8.0*Q0x*Q0x*Q0y
    tmp108 = tmp106*v0x
    tmp109 = 8.0*Q0x*Q0x
    tmp110 = 8.0*Q0w*Q0w*Q0z*v0x
    tmp111 = 8.0*Q0w*Q0x*Q0x*v0x
    tmp112 = tmp107*v0x
    tmp113 = 8.0*Q0w*Q0y*Q0y*v0x
    tmp114 = tmp109*Q0z*v0x
    tmp115 = 8.0*Q0w*Q0z*Q0z*v0x
    tmp116 = 8.0*Q0y*Q0y*Q0z*v0x
    tmp117 = 8.0*Q0x*Q0y*Q0y
    tmp118 = 8.0*Q0x*Q0z*Q0z
    tmp119 = 8.0*Q0w*Q0x*Q0y
    tmp120 = 8.0*Q0w*Q0x*Q0z
    tmp121 = tmp117*v0x
    tmp122 = tmp118*v0x
    tmp123 = 8.0*Q0y*Q0z*Q0z
    tmp124 = 8.0*Q0w*Q0y*Q0z
    tmp125 = 8.0*Q0x*Q0y*Q0z
    tmp126 = -4.0*Q0x*Q0x*Q0x
    tmp127 = tmp123*v0x
    tmp128 = -4.0*Q0w*Q0w*Q0w
    tmp129 = -4.0*Q0y*Q0y*Q0y*v0x
    tmp130 = -4.0*Q0z*Q0z*Q0z*v0x
    tmp131 = -8.0*Q0w
    tmp132 = tmp128*v0x
    tmp133 = -4.0*P0x*Q0z
    tmp134 = Q0x*Q0x*Q0x
    tmp135 = tmp131*Q0x
    tmp136 = tmp131*Q0y
    tmp137 = -4.0*P0y*Q0x
    tmp138 = -4.0*Q0x*Q0x*Q0y
    tmp139 = tmp131*Q0z
    tmp140 = -4.0*Q0x*Q0x*Q0z
    tmp141 = -4.0*P0z*Q0y
    tmp142 = -4.0*Q0z*Q0z*Q0z
    tmp143 = tmp136*Q1w*lf
    tmp144 = tmp126*v0x
    tmp145 = -8.0*Q0y
    tmp146 = -8.0*Q0x*Q0y
    tmp147 = -4.0*P1x*Q0w
    tmp148 = 8.0*Q0w
    tmp149 = tmp148*Q0y
    tmp150 = tmp136*Q1z*lf
    tmp151 = tmp139*Q1x*lf
    tmp152 = -8.0*Q0x
    tmp153 = tmp152*Q0z
    tmp154 = -4.0*P1x*Q0y
    tmp155 = -4.0*Q0y*Q0y*Q0z
    tmp156 = tmp146*Q1z*lf
    tmp157 = -4.0*Q0y*Q0y*Q0y
    tmp158 = -4.0*P0z
    tmp159 = tmp145*Q0z
    tmp160 = -4.0*P1y*Q0w
    tmp161 = tmp148*Q0x
    tmp162 = tmp161*Q1x*lf
    tmp163 = tmp148*Q0z
    tmp164 = -4.0*P1y*Q0z
    tmp165 = 8.0*Q0x*Q0y
    tmp166 = 8.0*Q0x
    tmp167 = tmp163*Q1w*lf
    tmp168 = -4.0*P1z*Q0w
    tmp169 = tmp166*Q0z
    tmp170 = tmp165*Q1y*lf
    tmp171 = -4.0*P1z*Q0x
    tmp172 = tmp161*Q1z*lf
    tmp173 = 8.0*Q0y
    tmp174 = -4.0*P1x*Q0x
    tmp175 = tmp169*Q1x*lf
    tmp176 = -4.0*Q0w*Q0x*Q0x
    tmp177 = tmp173*Q0z
    tmp178 = tmp177*Q1w*lf
    tmp179 = -4.0*P1y*Q0y
    tmp180 = tmp169*Q1y*lf
    tmp181 = tmp177*Q1y*lf
    tmp182 = -4.0*Q0w*Q0w
    tmp183 = -4.0*Q0w*Q0y*Q0y
    tmp184 = -4.0*Q0x*Q0x
    tmp185 = tmp182*Q1z*lf
    tmp186 = -4.0*Q0w*Q0z*Q0z
    tmp187 = -4.0*P1z*Q0z
    tmp188 = 4.0*tmp134
    tmp189 = tmp184*Q1w*lf
    tmp190 = Q0w*Q0w*Q0w
    tmp191 = -4.0*Q0y*Q0y
    tmp192 = Q0y*Q0y*Q0y
    tmp193 = Q0z*Q0z*Q0z
    tmp194 = 4.0*tmp192*v0x
    tmp195 = tmp184*Q1y*lf
    tmp196 = -4.0*Q0z*Q0z
    tmp197 = tmp184*Q1z*lf
    tmp198 = 4.0*tmp190*v0x
    tmp199 = 4.0*Q0x*Q0x*Q0z
    tmp200 = tmp191*Q1y*lf
    tmp201 = 4.0*Q0x*Q0x*Q0y
    tmp202 = -8.0*Q1w
    tmp203 = tmp191*Q1x*lf
    tmp204 = 4.0*tmp193*v0x
    tmp205 = 4.0*P0x*Q0w
    tmp206 = tmp202*Q1x
    tmp207 = tmp191*Q1z*lf
    tmp208 = tmp196*Q1y*lf
    tmp209 = tmp188*v0x
    tmp210 = 4.0*P0x*Q0y
    tmp211 = 4.0*Q0y*Q0y*Q0z
    tmp212 = 4.0*tmp192
    tmp213 = 4.0*Q0x*Q0x
    tmp214 = 4.0*P0y*Q0w
    tmp215 = 4.0*Q0w*Q0w
    tmp216 = tmp215*Q1y*lf
    tmp217 = tmp215*Q1w*lf
    tmp218 = 4.0*Q0y*Q0y
    tmp219 = tmp215*Q1x*lf
    tmp220 = 4.0*P0y*Q0z
    tmp221 = 4.0*P0x*Q0x
    tmp222 = tmp218*Q1w*lf
    tmp223 = 4.0*P0z*Q0w
    tmp224 = 4.0*Q0z*Q0z
    tmp225 = tmp224*Q1w*lf
    tmp226 = 4.0*P0z*Q0x
    tmp227 = 4.0*P0y*Q0y
    tmp228 = tmp213*Q1x*lf
    tmp229 = 4.0*P1x*Q0z
    tmp230 = tmp224*Q1x*lf
    tmp231 = -2.0*Q0x
    tmp232 = 8.0*Q1x
    tmp233 = 8.0*Q1w*Q1x
    tmp234 = tmp224*Q1z*lf
    tmp235 = 4.0*P1y*Q0x
    tmp236 = -2.0*Q0w
    tmp237 = -4.0*Q0w*Q0z
    tmp238 = -4.0*Q1w*Q1w
    tmp239 = tmp231*Q0y
    tmp240 = -4.0*Q0x*Q0y
    tmp241 = -4.0*Q0y
    tmp242 = -2.0*Q0y*Q0z
    tmp243 = 4.0*P0z*Q0z
    tmp244 = -4.0*P0x
    tmp245 = 4.0*P1z*Q0y
    tmp246 = 4.0*Q0w*Q0x*Q0x
    tmp247 = tmp231*Q0z
    tmp248 = 4.0*Q0x*Q0y
    tmp249 = 4.0*P1z
    tmp250 = 4.0*Q0w
    tmp251 = 4.0*Q0y
    tmp252 = -4.0*Q0w
    tmp253 = -4.0*P1x
    tmp254 = -4.0*P1y
    tmp255 = Q0x*Q0x
    tmp256 = 2.0*Q0w
    tmp257 = -4.0*Q0x
    tmp258 = tmp250*Q0y*Q0y
    tmp259 = tmp252*Q0y
    tmp260 = tmp250*Q0x
    tmp261 = tmp252*Q0x
    tmp262 = tmp257*Q0z
    tmp263 = -4.0*P1z
    tmp264 = -4.0*Q1x
    tmp265 = tmp241*Q0z
    tmp266 = 8.0*Q0w*Q0w
    tmp267 = -4.0*Q1w
    tmp268 = -4.0*P0y
    tmp269 = 4.0*Q0x*Q0z
    tmp270 = tmp250*Q0z*Q0z
    tmp271 = -4.0*Q1y
    tmp272 = Q1w*Q1w
    tmp273 = 4.0*Q0x
    tmp274 = -2.0*Q0y*Q0y
    tmp275 = tmp251*Q0z
    tmp276 = tmp250*Q0z
    tmp277 = 4.0*Q1w
    tmp278 = tmp266*Q0z
    tmp279 = 4.0*P0y
    tmp280 = tmp148*tmp255
    tmp281 = Q1y*Q1y
    tmp282 = 4.0*Q1x
    tmp283 = Q0w*Q0w
    tmp284 = tmp250*Q0y
    tmp285 = 4.0*P0x
    tmp286 = 4.0*P0z
    tmp287 = -2.0*tmp283
    tmp288 = Q0y*Q0y
    tmp289 = 4.0*P1x
    tmp290 = tmp148*tmp288
    tmp291 = tmp109*Q0z
    tmp292 = -2.0*tmp255
    tmp293 = 2.0*Q0z*Q0z
    tmp294 = 4.0*P1y
    tmp295 = -2.0*Q0z*Q0z
    tmp296 = Q1x*Q1x
    tmp297 = Q0z*Q0z
    tmp298 = 8.0*tmp272
    tmp299 = Q1z*Q1z
    tmp300 = 2.0*tmp255
    tmp301 = tmp148*tmp297
    tmp302 = 2.0*tmp283
    tmp303 = 8.0*tmp288*Q0z
    tmp304 = 2.0*tmp288
    tmp305 = 8.0*Q1w
    tmp306 = v0x*v0x
    tmp307 = tmp236*Q0x*v0y+tmp247*v0x+tmp242*v0y+-2.0*Q1w*Q1x*lf+-2.0*Q1y*Q1z*lf+tmp256*Q0y*v0x-tmp283*v0z-tmp297*v0z+tmp255*v0z+tmp288*v0z-P0z+P1z
    tmp308 = tmp236*Q0z*v0x+tmp239*v0x+tmp242*v0z+tmp256*Q0x*v0z-tmp283*v0y-tmp288*v0y-tmp272*lf-tmp281*lf+tmp255*v0y+tmp297*v0y+tmp296*lf+tmp299*lf-P0y+P1y
    tmp309 = v0y*v0y
    tmp310 = v0z*v0z
    tmp311 = 4.0*tmp193
    tmp312 = tmp236*Q0y*v0z+tmp239*v0y+tmp247*v0z+-2.0*Q1x*Q1y*lf+tmp256*Q0z*v0y+2.0*Q1w*Q1z*lf-tmp283*v0x-tmp255*v0x+tmp288*v0x+tmp297*v0x-P0x+P1x
    tmp313 = 4.0*tmp190
    tmp314 = lf*lf
    c2 = tmp308*tmp308+tmp307*tmp307+tmp312*tmp312
    dpQ0x = tmp4*v0z+tmp5*v0y+tmp7*v0y*v0z+tmp9*Q1x*lf*v0z+tmp15*v0x+tmp152*Q1y*Q1z*lf*v0z+tmp18*v0y+tmp10*Q1w*Q1z*lf*v0z+tmp25*v0z+tmp29*v0y+tmp30*v0y*v0z+tmp32*v0y+tmp35*v0y+tmp31*v0x+tmp41*v0y+tmp40*v0x+tmp36*Q1x*Q1y*lf*v0z+tmp43*v0x+tmp138*v0x*v0y+tmp140*v0x*v0z+tmp176*v0y*v0z+tmp252*tmp272*lf*v0z+tmp252*tmp281*lf*v0z+tmp52*v0y+tmp56*v0y+tmp53*v0x+tmp54*v0x+tmp66*v0y+tmp68*v0z+tmp72*v0y+tmp75*v0z+tmp79*v0z+tmp71*v0y*v0z+tmp73*v0y*v0z+tmp76*v0y*v0z+tmp85*tmp306+tmp85*tmp309+tmp85*tmp310+tmp89*v0y+tmp201*v0x*v0y+tmp199*v0x*v0z+tmp246*v0y*v0z+tmp250*tmp296*lf*v0z+tmp250*tmp299*lf*v0z+tmp98*v0y+tmp100*v0y+tmp97*v0x+tmp101*v0x+tmp108*v0y+tmp110*v0z+tmp112*v0y+tmp114*v0z+tmp116*v0z+tmp280*v0y*v0z+tmp290*v0y*v0z+tmp301*v0y*v0z+tmp124*tmp306+tmp124*tmp309+tmp124*tmp310+tmp127*v0y+tmp44*tmp309+tmp44*tmp310+tmp128*v0y*v0z+tmp129*v0y+tmp130*v0z+tmp268*Q0w*v0z+tmp137*v0y+tmp158*Q0x*v0z+tmp174*v0x+tmp154*v0y+tmp253*Q0z*v0z+tmp179*v0x+tmp168*v0y+tmp187*v0x+tmp46*tmp306+tmp46*tmp309+tmp48*tmp306+tmp48*tmp310+tmp88*tmp306+tmp313*v0y*v0z+tmp194*v0y+tmp204*v0z+tmp221*v0x+tmp210*v0y+tmp285*Q0z*v0z+tmp227*v0x+tmp223*v0y+tmp243*v0x+tmp294*Q0w*v0z+tmp235*v0y+tmp249*Q0x*v0z+tmp92*tmp310+tmp93*tmp309+tmp103*tmp309+tmp103*tmp310+tmp117*tmp306+tmp117*tmp309+tmp118*tmp306+tmp118*tmp310+tmp188*tmp306+tmp188*tmp309+tmp188*tmp310
    dpQ0y = tmp3*v0y*v0z+tmp8*v0y+tmp12*v0z+tmp11*Q1x*lf*v0x+tmp13*v0z+tmp131*Q1y*Q1z*lf*v0x+tmp15*v0y+tmp14*v0z+tmp145*Q1x*Q1y*lf*v0x+tmp19*v0z+tmp28*v0y*v0z+tmp20*v0y+tmp24*v0z+tmp34*v0z+tmp31*v0y+tmp33*Q1z*lf*v0x+tmp40*v0y+tmp43*v0y+tmp155*v0y*v0z+tmp183*v0x*v0z+tmp46*v0x*v0y+tmp257*tmp296*lf*v0x+tmp257*tmp299*lf*v0x+tmp53*v0y+tmp54*v0y+tmp58*v0z+tmp60*v0z+tmp63*v0y+tmp62*v0y*v0z+tmp65*v0y*v0z+tmp69*v0y*v0z+tmp70*v0z+tmp74*v0z+tmp77*v0z+tmp82*tmp306+tmp82*tmp309+tmp82*tmp310+tmp83*v0y+tmp86*v0y+tmp211*v0y*v0z+tmp258*v0x*v0z+tmp92*v0x*v0y+tmp273*tmp272*lf*v0x+tmp273*tmp281*lf*v0x+tmp97*v0y+tmp101*v0y+tmp102*v0z+tmp104*v0z+tmp105*v0y+tmp278*v0y*v0z+tmp291*v0y*v0z+tmp303*v0y*v0z+tmp111*v0z+tmp113*v0z+tmp115*v0z+tmp120*tmp306+tmp120*tmp309+tmp120*tmp310+tmp121*v0y+tmp122*v0y+tmp45*tmp306+tmp45*tmp310+tmp132*v0z+tmp138*tmp306+tmp138*tmp309+tmp144*v0y+tmp142*v0y*v0z+tmp244*Q0y*v0x+tmp158*Q0w*v0x+tmp141*v0z+tmp147*v0z+tmp174*v0y+tmp254*Q0x*v0x+tmp179*v0y+tmp164*v0z+tmp187*v0y+tmp50*tmp309+tmp50*tmp310+tmp90*tmp309+tmp198*v0z+tmp201*tmp310+tmp209*v0y+tmp311*v0y*v0z+tmp205*v0z+tmp221*v0y+tmp279*Q0x*v0x+tmp227*v0y+tmp220*v0z+tmp243*v0y+tmp289*Q0y*v0x+tmp249*Q0w*v0x+tmp245*v0z+tmp95*tmp306+tmp106*tmp306+tmp106*tmp310+tmp107*tmp306+tmp107*tmp309+tmp123*tmp309+tmp123*tmp310+tmp212*tmp306+tmp212*tmp309+tmp212*tmp310
    dpQ0z = tmp2*v0y*v0z+tmp8*v0z+tmp12*v0y+tmp131*Q1x*Q1y*lf*v0y+tmp15*v0z+tmp23*v0x+tmp16*v0y*v0z+tmp20*v0z+tmp24*v0y+tmp22*Q1z*lf*v0y+tmp38*v0x+tmp31*v0z+tmp39*v0x+tmp33*Q1x*lf*v0y+tmp173*Q1y*Q1z*lf*v0y+tmp40*v0z+tmp42*v0x+tmp43*v0z+tmp186*v0x*v0y+tmp49*v0x+tmp51*v0x+tmp48*v0x*v0z+tmp50*v0y*v0z+tmp53*v0z+tmp54*v0z+tmp55*tmp272*lf*v0y+tmp55*tmp281*lf*v0y+tmp63*v0z+tmp61*v0y*v0z+tmp64*v0y*v0z+tmp70*v0y+tmp74*v0y+tmp77*v0y+tmp81*tmp306+tmp81*tmp309+tmp81*tmp310+tmp83*v0z+tmp86*v0z+tmp84*v0y*v0z+tmp270*v0x*v0y+tmp94*v0x+tmp96*v0x+tmp93*v0x*v0z+tmp95*v0y*v0z+tmp97*v0z+tmp101*v0z+tmp99*tmp296*lf*v0y+tmp99*tmp299*lf*v0y+tmp105*v0z+tmp106*v0y*v0z+tmp107*v0y*v0z+tmp111*v0y+tmp113*v0y+tmp115*v0y+tmp119*tmp306+tmp119*tmp309+tmp119*tmp310+tmp121*v0z+tmp122*v0z+tmp123*v0y*v0z+tmp47*tmp306+tmp47*tmp309+tmp132*v0y+tmp140*tmp306+tmp140*tmp310+tmp144*v0z+tmp155*tmp309+tmp155*tmp310+tmp157*v0y*v0z+tmp244*Q0w*v0y+tmp133*v0x+tmp268*Q0z*v0y+tmp174*v0z+tmp160*v0x+tmp179*v0z+tmp171*v0x+tmp263*Q0y*v0y+tmp187*v0z+tmp91*tmp310+tmp198*v0y+tmp199*tmp309+tmp209*v0z+tmp211*tmp306+tmp212*v0y*v0z+tmp221*v0z+tmp214*v0x+tmp227*v0z+tmp226*v0x+tmp286*Q0y*v0y+tmp243*v0z+tmp289*Q0w*v0y+tmp229*v0x+tmp294*Q0z*v0y+tmp278*tmp306+tmp278*tmp309+tmp291*tmp306+tmp291*tmp310+tmp303*tmp309+tmp303*tmp310+tmp311*tmp306+tmp311*tmp309+tmp311*tmp310
    dpQ0w = tmp4*v0y+tmp5*v0z+tmp6*v0y*v0z+tmp13*v0x+tmp14*v0x+tmp18*v0z+tmp19*v0x+tmp23*v0y+tmp25*v0y+tmp29*v0z+tmp27*v0y*v0z+tmp32*v0z+tmp34*v0x+tmp35*v0z+tmp38*v0y+tmp39*v0y+tmp41*v0z+tmp42*v0y+tmp44*v0y*v0z+tmp45*v0x*v0z+tmp47*v0x*v0y+tmp49*v0y+tmp51*v0y+tmp52*v0z+tmp56*v0z+tmp58*v0x+tmp60*v0x+tmp59*v0y*v0z+tmp66*v0z+tmp68*v0y+tmp72*v0z+tmp75*v0y+tmp79*v0y+tmp78*v0y*v0z+tmp80*v0y*v0z+tmp87*tmp306+tmp87*tmp309+tmp87*tmp310+tmp89*v0z+tmp88*v0y*v0z+tmp90*v0x*v0z+tmp91*v0x*v0y+tmp94*v0y+tmp96*v0y+tmp98*v0z+tmp100*v0z+tmp102*v0x+tmp104*v0x+tmp103*v0y*v0z+tmp108*v0z+tmp110*v0y+tmp112*v0z+tmp114*v0y+tmp116*v0y+tmp117*v0y*v0z+tmp118*v0y*v0z+tmp125*tmp306+tmp125*tmp309+tmp125*tmp310+tmp127*v0z+tmp126*v0y*v0z+tmp129*v0z+tmp130*v0y+tmp133*v0y+tmp137*v0z+tmp141*v0x+tmp147*v0x+tmp154*v0z+tmp160*v0y+tmp164*v0x+tmp168*v0z+tmp171*v0y+tmp176*tmp309+tmp176*tmp310+tmp183*tmp306+tmp183*tmp310+tmp186*tmp306+tmp186*tmp309+tmp188*v0y*v0z+tmp194*v0z+tmp204*v0y+tmp205*v0x+tmp210*v0z+tmp214*v0y+tmp220*v0x+tmp223*v0z+tmp226*v0y+tmp229*v0y+tmp235*v0z+tmp245*v0x+tmp246*tmp306+tmp258*tmp309+tmp270*tmp310+tmp280*tmp309+tmp280*tmp310+tmp290*tmp306+tmp290*tmp310+tmp301*tmp306+tmp301*tmp309+tmp313*tmp306+tmp313*tmp309+tmp313*tmp310
    dpQ1x = tmp143*v0x+tmp151*v0x+tmp139*Q1y*lf*v0y+tmp146*Q1x*lf*v0x+tmp159*Q1x*lf*v0z+tmp161*Q1w*lf*v0y+tmp162*v0z+tmp149*Q1y*lf*v0z+tmp170*v0y+tmp169*Q1w*lf*v0x+tmp180*v0z+tmp178*v0y+tmp182*Q1x*lf*v0y+tmp189*v0z+tmp191*Q1w*lf*v0z+tmp203*v0y+tmp200*v0x+tmp208*v0x+tmp202*Q1y*Q1z*tmp314+tmp217*v0z+tmp216*v0x+tmp228*v0y+tmp213*Q1y*lf*v0x+tmp225*v0z+tmp230*v0y+tmp305*Q1y*Q1z*tmp314+tmp238*Q1x*tmp314+tmp268*Q1x*lf+tmp253*Q1y*lf+tmp263*Q1w*lf+tmp264*tmp281*tmp314+tmp285*Q1y*lf+tmp286*Q1w*lf+tmp294*Q1x*lf+tmp282*tmp299*tmp314+tmp298*Q1x*tmp314+tmp232*tmp281*tmp314+4.0*tmp296*Q1x*tmp314
    dpQ1y = tmp135*Q1y*lf*v0z+tmp150*v0x+tmp151*v0y+tmp172*v0y+tmp149*Q1x*lf*v0z+tmp163*Q1y*lf*v0x+tmp165*Q1x*lf*v0y+tmp170*v0x+tmp175*v0z+tmp169*Q1z*lf*v0x+tmp181*v0z+tmp177*Q1z*lf*v0y+tmp195*v0y+tmp197*v0z+tmp203*v0x+tmp207*v0z+tmp196*Q1x*lf*v0x+tmp208*v0y+tmp206*Q1z*tmp314+tmp219*v0x+tmp216*v0y+tmp215*Q1z*lf*v0z+tmp228*v0x+tmp218*Q1y*lf*v0y+tmp234*v0z+tmp233*Q1z*tmp314+-4.0*tmp296*Q1y*tmp314+tmp253*Q1x*lf+tmp254*Q1y*lf+tmp263*Q1z*lf+tmp271*tmp299*tmp314+4.0*tmp272*Q1y*tmp314+tmp285*Q1x*lf+tmp279*Q1y*lf+tmp286*Q1z*lf+8.0*tmp296*Q1y*tmp314+8.0*Q1y*tmp299*tmp314+4.0*tmp281*Q1y*tmp314
    dpQ1z = tmp143*v0z+tmp136*Q1y*lf*v0x+tmp139*Q1z*lf*v0x+tmp146*Q1w*lf*v0y+tmp156*v0x+tmp153*Q1w*lf*v0z+tmp159*Q1z*lf*v0z+tmp161*Q1y*lf*v0y+tmp172*v0z+tmp167*v0y+tmp180*v0x+tmp181*v0y+tmp182*Q1w*lf*v0x+tmp185*v0y+tmp189*v0x+tmp195*v0z+tmp200*v0z+tmp207*v0y+tmp206*Q1y*tmp314+tmp216*v0z+tmp213*Q1z*lf*v0y+tmp222*v0x+tmp225*v0x+tmp224*Q1y*lf*v0z+tmp234*v0y+tmp233*Q1y*tmp314+tmp238*Q1z*tmp314+-4.0*tmp281*Q1z*tmp314+tmp244*Q1w*lf+tmp268*Q1z*lf+tmp263*Q1y*lf+4.0*tmp296*Q1z*tmp314+tmp286*Q1y*lf+tmp289*Q1w*lf+tmp294*Q1z*lf+tmp298*Q1z*tmp314+8.0*tmp281*Q1z*tmp314+4.0*tmp299*Q1z*tmp314
    dpQ1w = tmp135*Q1w*lf*v0z+tmp136*Q1x*lf*v0x+tmp150*v0z+tmp156*v0y+tmp153*Q1z*lf*v0z+tmp162*v0y+tmp167*v0x+tmp163*Q1z*lf*v0y+tmp165*Q1w*lf*v0x+tmp175*v0x+tmp178*v0z+tmp177*Q1x*lf*v0y+tmp185*v0x+tmp189*v0y+tmp184*Q1x*lf*v0z+tmp197*v0x+tmp203*v0z+tmp196*Q1w*lf*v0y+-8.0*Q1x*Q1y*Q1z*tmp314+tmp217*v0y+tmp219*v0z+tmp222*v0y+tmp218*Q1z*lf*v0x+tmp230*v0z+tmp234*v0x+tmp232*Q1y*Q1z*tmp314+tmp244*Q1z*lf+tmp254*Q1w*lf+tmp263*Q1x*lf+tmp267*tmp296*tmp314+tmp267*tmp299*tmp314+tmp279*Q1w*lf+tmp286*Q1x*lf+tmp289*Q1z*lf+tmp277*tmp281*tmp314+tmp305*tmp296*tmp314+tmp305*tmp299*tmp314+4.0*tmp272*Q1w*tmp314
    dpP0x = tmp237*v0y+tmp267*Q1z*lf+tmp284*v0z+tmp248*v0y+tmp269*v0z+tmp282*Q1y*lf+tmp274*v0x+tmp295*v0x+tmp302*v0x+tmp300*v0x+2.0*P0x+-2.0*P1x
    dpP0y = tmp261*v0z+tmp276*v0x+tmp248*v0x+tmp275*v0z+tmp292*v0y+tmp295*v0y+-2.0*tmp296*lf+-2.0*tmp299*lf+tmp302*v0y+tmp304*v0y+2.0*tmp272*lf+2.0*tmp281*lf+2.0*P0y+-2.0*P1y
    dpP0z = tmp259*v0x+tmp260*v0y+tmp269*v0x+tmp275*v0y+tmp277*Q1x*lf+4.0*Q1y*Q1z*lf+tmp292*v0z+tmp274*v0z+tmp302*v0z+tmp293*v0z+2.0*P0z+-2.0*P1z
    dpP1x = tmp259*v0z+tmp240*v0y+tmp262*v0z+tmp264*Q1y*lf+tmp276*v0y+tmp277*Q1z*lf+tmp287*v0x+tmp292*v0x+tmp304*v0x+tmp293*v0x+-2.0*P0x+2.0*P1x
    dpP1y = tmp237*v0x+tmp240*v0x+tmp265*v0z+tmp260*v0z+tmp287*v0y+tmp274*v0y+-2.0*tmp272*lf+-2.0*tmp281*lf+tmp300*v0y+tmp293*v0y+2.0*tmp296*lf+2.0*tmp299*lf+-2.0*P0y+2.0*P1y
    dpP1z = tmp261*v0y+tmp262*v0x+tmp265*v0y+tmp267*Q1x*lf+tmp271*Q1z*lf+tmp284*v0x+tmp287*v0z+tmp295*v0z+tmp300*v0z+tmp304*v0z+-2.0*P0z+2.0*P1z


    return c2,dpQ0x,dpQ0y,dpQ0z,dpQ0w,dpQ1x,dpQ1y,dpQ1z,dpQ1w,dpP0x,dpP0y,dpP0z,dpP1x,dpP1y,dpP1z
    
def locSpringGradient2(Q0,Q1,P0,P1,v0,lf):
    Q0x = Q0.x
    Q0y = Q0.y
    Q0z = Q0.z
    Q0w = Q0.w
    Q1x = Q1.x
    Q1y = Q1.y
    Q1z = Q1.z
    Q1w = Q1.w  
    v0x = v0.x
    v0y = v0.y
    v0z = v0.z 
    P0x = P0.x
    P0y = P0.y
    P0z = P0.z
    P1x = P1.x
    P1y = P1.y
    P1z = P1.z    
    
    tmp0 = (Q0x*Q0x*v0y+Q0z*Q0z*v0y+Q1x*Q1x*lf+Q1z*Q1z*lf+-1.0*Q0w*Q0w*v0y+-1.0*Q0y*Q0y*v0y+-1.0*Q1w*Q1w*lf+-1.0*Q1y*Q1y*lf+-P0y+-2.0*Q0w*Q0z*v0x+-2.0*Q0x*Q0y*v0x+-2.0*Q0y*Q0z*v0z+2.0*Q0w*Q0x*v0z+P1y)*(Q0x*Q0x*v0y+Q0z*Q0z*v0y+Q1x*Q1x*lf+Q1z*Q1z*lf+-1.0*Q0w*Q0w*v0y+-1.0*Q0y*Q0y*v0y+-1.0*Q1w*Q1w*lf+-1.0*Q1y*Q1y*lf+-P0y+-2.0*Q0w*Q0z*v0x+-2.0*Q0x*Q0y*v0x+-2.0*Q0y*Q0z*v0z+2.0*Q0w*Q0x*v0z+P1y)+(Q0x*Q0x*v0z+Q0y*Q0y*v0z+-1.0*Q0w*Q0w*v0z+-1.0*Q0z*Q0z*v0z+-P0z+-2.0*Q0w*Q0x*v0y+-2.0*Q0x*Q0z*v0x+-2.0*Q0y*Q0z*v0y+-2.0*Q1w*Q1x*lf+-2.0*Q1y*Q1z*lf+2.0*Q0w*Q0y*v0x+P1z)*(Q0x*Q0x*v0z+Q0y*Q0y*v0z+-1.0*Q0w*Q0w*v0z+-1.0*Q0z*Q0z*v0z+-P0z+-2.0*Q0w*Q0x*v0y+-2.0*Q0x*Q0z*v0x+-2.0*Q0y*Q0z*v0y+-2.0*Q1w*Q1x*lf+-2.0*Q1y*Q1z*lf+2.0*Q0w*Q0y*v0x+P1z)+(Q0y*Q0y*v0x+Q0z*Q0z*v0x+-1.0*Q0w*Q0w*v0x+-1.0*Q0x*Q0x*v0x+-P0x+-2.0*Q0w*Q0y*v0z+-2.0*Q0x*Q0y*v0y+-2.0*Q0x*Q0z*v0z+-2.0*Q1x*Q1y*lf+2.0*Q0w*Q0z*v0y+2.0*Q1w*Q1z*lf+P1x)*(Q0y*Q0y*v0x+Q0z*Q0z*v0x+-1.0*Q0w*Q0w*v0x+-1.0*Q0x*Q0x*v0x+-P0x+-2.0*Q0w*Q0y*v0z+-2.0*Q0x*Q0y*v0y+-2.0*Q0x*Q0z*v0z+-2.0*Q1x*Q1y*lf+2.0*Q0w*Q0z*v0y+2.0*Q1w*Q1z*lf+P1x)
    tmp1 = Q0x*Q0x
    tmp2 = Q1z*Q1z
    tmp3 = Q0w*Q0w*Q0w
    tmp4 = Q0z*Q0z*Q0z
    tmp5 = v0x*v0x
    tmp6 = Q0x*Q0x*Q0x
    tmp7 = Q0w*Q0w
    tmp8 = v0z*v0z
    tmp9 = v0y*v0y
    tmp10 = Q1w*Q1w
    tmp11 = 1.0/sqrt(tmp0)
    tmp12 = Q0y*Q0y*Q0y
    tmp13 = lf*lf
    tmp14 = Q0z*Q0z
    tmp15 = Q1x*Q1x
    tmp16 = Q0y*Q0y
    tmp17 = Q1y*Q1y
    c2 = sqrt(tmp0)
    dpQ0x = -2.0*tmp11*tmp7*tmp9*Q0x+-2.0*tmp11*tmp7*tmp8*Q0x+-2.0*tmp11*tmp3*v0y*v0z+-2.0*tmp11*tmp1*Q0w*v0y*v0z+-2.0*tmp11*tmp1*Q0y*v0x*v0y+-2.0*tmp11*tmp1*Q0z*v0x*v0z+-2.0*tmp11*tmp16*tmp5*Q0x+-2.0*tmp11*tmp16*tmp9*Q0x+-2.0*tmp11*tmp12*v0x*v0y+-2.0*tmp11*tmp14*tmp5*Q0x+-2.0*tmp11*tmp14*tmp8*Q0x+-2.0*tmp11*tmp4*v0x*v0z+-2.0*tmp11*tmp10*Q0w*lf*v0z+-2.0*tmp11*tmp10*Q0x*lf*v0y+-2.0*tmp11*tmp15*Q0y*lf*v0x+-2.0*tmp11*tmp17*Q0w*lf*v0z+-2.0*tmp11*tmp17*Q0x*lf*v0y+-2.0*tmp11*tmp2*Q0y*lf*v0x+-2.0*tmp11*P0y*Q0w*v0z+-2.0*tmp11*P0y*Q0x*v0y+-2.0*tmp11*P0z*Q0x*v0z+-2.0*tmp11*P1x*Q0x*v0x+-2.0*tmp11*P1x*Q0y*v0y+-2.0*tmp11*P1x*Q0z*v0z+-2.0*tmp11*P1y*Q0y*v0x+-2.0*tmp11*P1z*Q0w*v0y+-2.0*tmp11*P1z*Q0z*v0x+-4.0*tmp11*tmp7*Q0y*v0x*v0y+-4.0*tmp11*tmp7*Q0z*v0x*v0z+-4.0*tmp11*tmp1*Q0w*v0y*v0z+-4.0*tmp11*tmp1*Q0y*v0x*v0y+-4.0*tmp11*tmp1*Q0z*v0x*v0z+-4.0*tmp11*tmp16*Q0w*v0y*v0z+-4.0*tmp11*tmp16*Q0z*v0x*v0z+-4.0*tmp11*tmp14*Q0w*v0y*v0z+-4.0*tmp11*tmp14*Q0y*v0x*v0y+-4.0*tmp11*tmp5*Q0w*Q0y*Q0z+-4.0*tmp11*tmp9*Q0w*Q0y*Q0z+-4.0*tmp11*tmp8*Q0w*Q0y*Q0z+-4.0*tmp11*Q0x*Q1w*Q1x*lf*v0z+-4.0*tmp11*Q0x*Q1w*Q1z*lf*v0x+-4.0*tmp11*Q0x*Q1y*Q1z*lf*v0z+-4.0*tmp11*Q0y*Q1w*Q1z*lf*v0y+-4.0*tmp11*Q0z*Q1w*Q1z*lf*v0z+-8.0*tmp11*Q0w*Q0x*Q0y*v0x*v0z+-8.0*tmp11*Q0w*Q0x*Q0z*v0x*v0y+-8.0*tmp11*Q0x*Q0y*Q0z*v0y*v0z+2.0*tmp11*tmp7*tmp5*Q0x+2.0*tmp11*tmp3*v0y*v0z+2.0*tmp11*tmp1*Q0w*v0y*v0z+2.0*tmp11*tmp1*Q0y*v0x*v0y+2.0*tmp11*tmp1*Q0z*v0x*v0z+2.0*tmp11*tmp6*tmp5+2.0*tmp11*tmp6*tmp9+2.0*tmp11*tmp6*tmp8+2.0*tmp11*tmp16*tmp8*Q0x+2.0*tmp11*tmp12*v0x*v0y+2.0*tmp11*tmp14*tmp9*Q0x+2.0*tmp11*tmp4*v0x*v0z+2.0*tmp11*tmp10*Q0y*lf*v0x+2.0*tmp11*tmp15*Q0w*lf*v0z+2.0*tmp11*tmp15*Q0x*lf*v0y+2.0*tmp11*tmp17*Q0y*lf*v0x+2.0*tmp11*tmp2*Q0w*lf*v0z+2.0*tmp11*tmp2*Q0x*lf*v0y+2.0*tmp11*P0x*Q0x*v0x+2.0*tmp11*P0x*Q0y*v0y+2.0*tmp11*P0x*Q0z*v0z+2.0*tmp11*P0y*Q0y*v0x+2.0*tmp11*P0z*Q0w*v0y+2.0*tmp11*P0z*Q0z*v0x+2.0*tmp11*P1y*Q0w*v0z+2.0*tmp11*P1y*Q0x*v0y+2.0*tmp11*P1z*Q0x*v0z+4.0*tmp11*tmp7*tmp9*Q0x+4.0*tmp11*tmp7*tmp8*Q0x+4.0*tmp11*tmp7*Q0y*v0x*v0y+4.0*tmp11*tmp7*Q0z*v0x*v0z+4.0*tmp11*tmp1*Q0w*v0y*v0z+4.0*tmp11*tmp1*Q0y*v0x*v0y+4.0*tmp11*tmp1*Q0z*v0x*v0z+4.0*tmp11*tmp16*tmp5*Q0x+4.0*tmp11*tmp16*tmp9*Q0x+4.0*tmp11*tmp16*Q0w*v0y*v0z+4.0*tmp11*tmp16*Q0z*v0x*v0z+4.0*tmp11*tmp14*tmp5*Q0x+4.0*tmp11*tmp14*tmp8*Q0x+4.0*tmp11*tmp14*Q0w*v0y*v0z+4.0*tmp11*tmp14*Q0y*v0x*v0y+4.0*tmp11*tmp5*Q0w*Q0y*Q0z+4.0*tmp11*tmp9*Q0w*Q0y*Q0z+4.0*tmp11*tmp8*Q0w*Q0y*Q0z+4.0*tmp11*Q0w*Q1w*Q1x*lf*v0y+4.0*tmp11*Q0w*Q1y*Q1z*lf*v0y+4.0*tmp11*Q0x*Q1x*Q1y*lf*v0x+4.0*tmp11*Q0y*Q1x*Q1y*lf*v0y+4.0*tmp11*Q0z*Q1w*Q1x*lf*v0x+4.0*tmp11*Q0z*Q1x*Q1y*lf*v0z+4.0*tmp11*Q0z*Q1y*Q1z*lf*v0x+8.0*tmp11*Q0w*Q0x*Q0y*v0x*v0z+8.0*tmp11*Q0w*Q0x*Q0z*v0x*v0y+8.0*tmp11*Q0x*Q0y*Q0z*v0y*v0z
    dpQ0y = -2.0*tmp11*tmp7*tmp5*Q0y+-2.0*tmp11*tmp7*tmp8*Q0y+-2.0*tmp11*tmp3*v0x*v0z+-2.0*tmp11*tmp1*tmp5*Q0y+-2.0*tmp11*tmp1*tmp9*Q0y+-2.0*tmp11*tmp6*v0x*v0y+-2.0*tmp11*tmp16*Q0w*v0x*v0z+-2.0*tmp11*tmp16*Q0x*v0x*v0y+-2.0*tmp11*tmp16*Q0z*v0y*v0z+-2.0*tmp11*tmp14*tmp9*Q0y+-2.0*tmp11*tmp14*tmp8*Q0y+-2.0*tmp11*tmp4*v0y*v0z+-2.0*tmp11*tmp15*Q0x*lf*v0x+-2.0*tmp11*tmp15*Q0y*lf*v0y+-2.0*tmp11*tmp15*Q0z*lf*v0z+-2.0*tmp11*tmp2*Q0x*lf*v0x+-2.0*tmp11*tmp2*Q0y*lf*v0y+-2.0*tmp11*tmp2*Q0z*lf*v0z+-2.0*tmp11*P0x*Q0y*v0x+-2.0*tmp11*P0z*Q0w*v0x+-2.0*tmp11*P0z*Q0y*v0z+-2.0*tmp11*P1x*Q0w*v0z+-2.0*tmp11*P1x*Q0x*v0y+-2.0*tmp11*P1y*Q0x*v0x+-2.0*tmp11*P1y*Q0y*v0y+-2.0*tmp11*P1y*Q0z*v0z+-2.0*tmp11*P1z*Q0z*v0y+-4.0*tmp11*tmp7*Q0x*v0x*v0y+-4.0*tmp11*tmp7*Q0z*v0y*v0z+-4.0*tmp11*tmp1*Q0w*v0x*v0z+-4.0*tmp11*tmp1*Q0z*v0y*v0z+-4.0*tmp11*tmp16*Q0w*v0x*v0z+-4.0*tmp11*tmp16*Q0x*v0x*v0y+-4.0*tmp11*tmp16*Q0z*v0y*v0z+-4.0*tmp11*tmp14*Q0w*v0x*v0z+-4.0*tmp11*tmp14*Q0x*v0x*v0y+-4.0*tmp11*tmp5*Q0w*Q0x*Q0z+-4.0*tmp11*tmp9*Q0w*Q0x*Q0z+-4.0*tmp11*tmp8*Q0w*Q0x*Q0z+-4.0*tmp11*Q0w*Q1w*Q1x*lf*v0x+-4.0*tmp11*Q0w*Q1w*Q1z*lf*v0z+-4.0*tmp11*Q0w*Q1y*Q1z*lf*v0x+-4.0*tmp11*Q0x*Q1w*Q1z*lf*v0y+-4.0*tmp11*Q0y*Q1w*Q1x*lf*v0z+-4.0*tmp11*Q0y*Q1x*Q1y*lf*v0x+-4.0*tmp11*Q0y*Q1y*Q1z*lf*v0z+-8.0*tmp11*Q0w*Q0x*Q0y*v0y*v0z+-8.0*tmp11*Q0w*Q0y*Q0z*v0x*v0y+-8.0*tmp11*Q0x*Q0y*Q0z*v0x*v0z+2.0*tmp11*tmp7*tmp9*Q0y+2.0*tmp11*tmp3*v0x*v0z+2.0*tmp11*tmp1*tmp8*Q0y+2.0*tmp11*tmp6*v0x*v0y+2.0*tmp11*tmp16*Q0w*v0x*v0z+2.0*tmp11*tmp16*Q0x*v0x*v0y+2.0*tmp11*tmp16*Q0z*v0y*v0z+2.0*tmp11*tmp12*tmp5+2.0*tmp11*tmp12*tmp9+2.0*tmp11*tmp12*tmp8+2.0*tmp11*tmp14*tmp5*Q0y+2.0*tmp11*tmp4*v0y*v0z+2.0*tmp11*tmp10*Q0x*lf*v0x+2.0*tmp11*tmp10*Q0y*lf*v0y+2.0*tmp11*tmp10*Q0z*lf*v0z+2.0*tmp11*tmp17*Q0x*lf*v0x+2.0*tmp11*tmp17*Q0y*lf*v0y+2.0*tmp11*tmp17*Q0z*lf*v0z+2.0*tmp11*P0x*Q0w*v0z+2.0*tmp11*P0x*Q0x*v0y+2.0*tmp11*P0y*Q0x*v0x+2.0*tmp11*P0y*Q0y*v0y+2.0*tmp11*P0y*Q0z*v0z+2.0*tmp11*P0z*Q0z*v0y+2.0*tmp11*P1x*Q0y*v0x+2.0*tmp11*P1z*Q0w*v0x+2.0*tmp11*P1z*Q0y*v0z+4.0*tmp11*tmp7*tmp5*Q0y+4.0*tmp11*tmp7*tmp8*Q0y+4.0*tmp11*tmp7*Q0x*v0x*v0y+4.0*tmp11*tmp7*Q0z*v0y*v0z+4.0*tmp11*tmp1*tmp5*Q0y+4.0*tmp11*tmp1*tmp9*Q0y+4.0*tmp11*tmp1*Q0w*v0x*v0z+4.0*tmp11*tmp1*Q0z*v0y*v0z+4.0*tmp11*tmp16*Q0w*v0x*v0z+4.0*tmp11*tmp16*Q0x*v0x*v0y+4.0*tmp11*tmp16*Q0z*v0y*v0z+4.0*tmp11*tmp14*tmp9*Q0y+4.0*tmp11*tmp14*tmp8*Q0y+4.0*tmp11*tmp14*Q0w*v0x*v0z+4.0*tmp11*tmp14*Q0x*v0x*v0y+4.0*tmp11*tmp5*Q0w*Q0x*Q0z+4.0*tmp11*tmp9*Q0w*Q0x*Q0z+4.0*tmp11*tmp8*Q0w*Q0x*Q0z+4.0*tmp11*Q0w*Q1x*Q1y*lf*v0z+4.0*tmp11*Q0x*Q1x*Q1y*lf*v0y+4.0*tmp11*Q0y*Q1w*Q1z*lf*v0x+4.0*tmp11*Q0z*Q1w*Q1x*lf*v0y+4.0*tmp11*Q0z*Q1y*Q1z*lf*v0y+8.0*tmp11*Q0w*Q0x*Q0y*v0y*v0z+8.0*tmp11*Q0w*Q0y*Q0z*v0x*v0y+8.0*tmp11*Q0x*Q0y*Q0z*v0x*v0z
    dpQ0z = -2.0*tmp11*tmp7*tmp5*Q0z+-2.0*tmp11*tmp7*tmp9*Q0z+-2.0*tmp11*tmp3*v0x*v0y+-2.0*tmp11*tmp1*tmp5*Q0z+-2.0*tmp11*tmp1*tmp8*Q0z+-2.0*tmp11*tmp6*v0x*v0z+-2.0*tmp11*tmp16*tmp9*Q0z+-2.0*tmp11*tmp16*tmp8*Q0z+-2.0*tmp11*tmp12*v0y*v0z+-2.0*tmp11*tmp14*Q0w*v0x*v0y+-2.0*tmp11*tmp14*Q0x*v0x*v0z+-2.0*tmp11*tmp14*Q0y*v0y*v0z+-2.0*tmp11*tmp10*Q0z*lf*v0y+-2.0*tmp11*tmp15*Q0w*lf*v0x+-2.0*tmp11*tmp15*Q0y*lf*v0z+-2.0*tmp11*tmp17*Q0z*lf*v0y+-2.0*tmp11*tmp2*Q0w*lf*v0x+-2.0*tmp11*tmp2*Q0y*lf*v0z+-2.0*tmp11*P0x*Q0w*v0y+-2.0*tmp11*P0x*Q0z*v0x+-2.0*tmp11*P0y*Q0z*v0y+-2.0*tmp11*P1x*Q0x*v0z+-2.0*tmp11*P1y*Q0w*v0x+-2.0*tmp11*P1y*Q0y*v0z+-2.0*tmp11*P1z*Q0x*v0x+-2.0*tmp11*P1z*Q0y*v0y+-2.0*tmp11*P1z*Q0z*v0z+-4.0*tmp11*tmp7*Q0x*v0x*v0z+-4.0*tmp11*tmp7*Q0y*v0y*v0z+-4.0*tmp11*tmp1*Q0w*v0x*v0y+-4.0*tmp11*tmp1*Q0y*v0y*v0z+-4.0*tmp11*tmp16*Q0w*v0x*v0y+-4.0*tmp11*tmp16*Q0x*v0x*v0z+-4.0*tmp11*tmp14*Q0w*v0x*v0y+-4.0*tmp11*tmp14*Q0x*v0x*v0z+-4.0*tmp11*tmp14*Q0y*v0y*v0z+-4.0*tmp11*tmp5*Q0w*Q0x*Q0y+-4.0*tmp11*tmp9*Q0w*Q0x*Q0y+-4.0*tmp11*tmp8*Q0w*Q0x*Q0y+-4.0*tmp11*Q0w*Q1x*Q1y*lf*v0y+-4.0*tmp11*Q0x*Q1w*Q1z*lf*v0z+-4.0*tmp11*Q0z*Q1x*Q1y*lf*v0x+-8.0*tmp11*Q0w*Q0x*Q0z*v0y*v0z+-8.0*tmp11*Q0w*Q0y*Q0z*v0x*v0z+-8.0*tmp11*Q0x*Q0y*Q0z*v0x*v0y+2.0*tmp11*tmp7*tmp8*Q0z+2.0*tmp11*tmp3*v0x*v0y+2.0*tmp11*tmp1*tmp9*Q0z+2.0*tmp11*tmp6*v0x*v0z+2.0*tmp11*tmp16*tmp5*Q0z+2.0*tmp11*tmp12*v0y*v0z+2.0*tmp11*tmp14*Q0w*v0x*v0y+2.0*tmp11*tmp14*Q0x*v0x*v0z+2.0*tmp11*tmp14*Q0y*v0y*v0z+2.0*tmp11*tmp4*tmp5+2.0*tmp11*tmp4*tmp9+2.0*tmp11*tmp4*tmp8+2.0*tmp11*tmp10*Q0w*lf*v0x+2.0*tmp11*tmp10*Q0y*lf*v0z+2.0*tmp11*tmp15*Q0z*lf*v0y+2.0*tmp11*tmp17*Q0w*lf*v0x+2.0*tmp11*tmp17*Q0y*lf*v0z+2.0*tmp11*tmp2*Q0z*lf*v0y+2.0*tmp11*P0x*Q0x*v0z+2.0*tmp11*P0y*Q0w*v0x+2.0*tmp11*P0y*Q0y*v0z+2.0*tmp11*P0z*Q0x*v0x+2.0*tmp11*P0z*Q0y*v0y+2.0*tmp11*P0z*Q0z*v0z+2.0*tmp11*P1x*Q0w*v0y+2.0*tmp11*P1x*Q0z*v0x+2.0*tmp11*P1y*Q0z*v0y+4.0*tmp11*tmp7*tmp5*Q0z+4.0*tmp11*tmp7*tmp9*Q0z+4.0*tmp11*tmp7*Q0x*v0x*v0z+4.0*tmp11*tmp7*Q0y*v0y*v0z+4.0*tmp11*tmp1*tmp5*Q0z+4.0*tmp11*tmp1*tmp8*Q0z+4.0*tmp11*tmp1*Q0w*v0x*v0y+4.0*tmp11*tmp1*Q0y*v0y*v0z+4.0*tmp11*tmp16*tmp9*Q0z+4.0*tmp11*tmp16*tmp8*Q0z+4.0*tmp11*tmp16*Q0w*v0x*v0y+4.0*tmp11*tmp16*Q0x*v0x*v0z+4.0*tmp11*tmp14*Q0w*v0x*v0y+4.0*tmp11*tmp14*Q0x*v0x*v0z+4.0*tmp11*tmp14*Q0y*v0y*v0z+4.0*tmp11*tmp5*Q0w*Q0x*Q0y+4.0*tmp11*tmp9*Q0w*Q0x*Q0y+4.0*tmp11*tmp8*Q0w*Q0x*Q0y+4.0*tmp11*Q0w*Q1w*Q1z*lf*v0y+4.0*tmp11*Q0x*Q1w*Q1x*lf*v0x+4.0*tmp11*Q0x*Q1x*Q1y*lf*v0z+4.0*tmp11*Q0x*Q1y*Q1z*lf*v0x+4.0*tmp11*Q0y*Q1w*Q1x*lf*v0y+4.0*tmp11*Q0y*Q1y*Q1z*lf*v0y+4.0*tmp11*Q0z*Q1w*Q1x*lf*v0z+4.0*tmp11*Q0z*Q1w*Q1z*lf*v0x+4.0*tmp11*Q0z*Q1y*Q1z*lf*v0z+8.0*tmp11*Q0w*Q0x*Q0z*v0y*v0z+8.0*tmp11*Q0w*Q0y*Q0z*v0x*v0z+8.0*tmp11*Q0x*Q0y*Q0z*v0x*v0y
    dpQ0w = -2.0*tmp11*tmp7*Q0x*v0y*v0z+-2.0*tmp11*tmp7*Q0y*v0x*v0z+-2.0*tmp11*tmp7*Q0z*v0x*v0y+-2.0*tmp11*tmp1*tmp9*Q0w+-2.0*tmp11*tmp1*tmp8*Q0w+-2.0*tmp11*tmp6*v0y*v0z+-2.0*tmp11*tmp16*tmp5*Q0w+-2.0*tmp11*tmp16*tmp8*Q0w+-2.0*tmp11*tmp12*v0x*v0z+-2.0*tmp11*tmp14*tmp5*Q0w+-2.0*tmp11*tmp14*tmp9*Q0w+-2.0*tmp11*tmp4*v0x*v0y+-2.0*tmp11*tmp10*Q0x*lf*v0z+-2.0*tmp11*tmp15*Q0w*lf*v0y+-2.0*tmp11*tmp15*Q0z*lf*v0x+-2.0*tmp11*tmp17*Q0x*lf*v0z+-2.0*tmp11*tmp2*Q0w*lf*v0y+-2.0*tmp11*tmp2*Q0z*lf*v0x+-2.0*tmp11*P0x*Q0z*v0y+-2.0*tmp11*P0y*Q0x*v0z+-2.0*tmp11*P0z*Q0y*v0x+-2.0*tmp11*P1x*Q0w*v0x+-2.0*tmp11*P1x*Q0y*v0z+-2.0*tmp11*P1y*Q0w*v0y+-2.0*tmp11*P1y*Q0z*v0x+-2.0*tmp11*P1z*Q0w*v0z+-2.0*tmp11*P1z*Q0x*v0y+-4.0*tmp11*tmp7*Q0x*v0y*v0z+-4.0*tmp11*tmp7*Q0y*v0x*v0z+-4.0*tmp11*tmp7*Q0z*v0x*v0y+-4.0*tmp11*tmp1*Q0y*v0x*v0z+-4.0*tmp11*tmp1*Q0z*v0x*v0y+-4.0*tmp11*tmp16*Q0x*v0y*v0z+-4.0*tmp11*tmp16*Q0z*v0x*v0y+-4.0*tmp11*tmp14*Q0x*v0y*v0z+-4.0*tmp11*tmp14*Q0y*v0x*v0z+-4.0*tmp11*tmp5*Q0x*Q0y*Q0z+-4.0*tmp11*tmp9*Q0x*Q0y*Q0z+-4.0*tmp11*tmp8*Q0x*Q0y*Q0z+-4.0*tmp11*Q0w*Q1w*Q1z*lf*v0x+-4.0*tmp11*Q0y*Q1w*Q1x*lf*v0x+-4.0*tmp11*Q0y*Q1w*Q1z*lf*v0z+-4.0*tmp11*Q0y*Q1y*Q1z*lf*v0x+-4.0*tmp11*Q0z*Q1x*Q1y*lf*v0y+-8.0*tmp11*Q0w*Q0x*Q0y*v0x*v0y+-8.0*tmp11*Q0w*Q0x*Q0z*v0x*v0z+-8.0*tmp11*Q0w*Q0y*Q0z*v0y*v0z+2.0*tmp11*tmp7*Q0x*v0y*v0z+2.0*tmp11*tmp7*Q0y*v0x*v0z+2.0*tmp11*tmp7*Q0z*v0x*v0y+2.0*tmp11*tmp3*tmp5+2.0*tmp11*tmp3*tmp9+2.0*tmp11*tmp3*tmp8+2.0*tmp11*tmp1*tmp5*Q0w+2.0*tmp11*tmp6*v0y*v0z+2.0*tmp11*tmp16*tmp9*Q0w+2.0*tmp11*tmp12*v0x*v0z+2.0*tmp11*tmp14*tmp8*Q0w+2.0*tmp11*tmp4*v0x*v0y+2.0*tmp11*tmp10*Q0w*lf*v0y+2.0*tmp11*tmp10*Q0z*lf*v0x+2.0*tmp11*tmp15*Q0x*lf*v0z+2.0*tmp11*tmp17*Q0w*lf*v0y+2.0*tmp11*tmp17*Q0z*lf*v0x+2.0*tmp11*tmp2*Q0x*lf*v0z+2.0*tmp11*P0x*Q0w*v0x+2.0*tmp11*P0x*Q0y*v0z+2.0*tmp11*P0y*Q0w*v0y+2.0*tmp11*P0y*Q0z*v0x+2.0*tmp11*P0z*Q0w*v0z+2.0*tmp11*P0z*Q0x*v0y+2.0*tmp11*P1x*Q0z*v0y+2.0*tmp11*P1y*Q0x*v0z+2.0*tmp11*P1z*Q0y*v0x+4.0*tmp11*tmp7*Q0x*v0y*v0z+4.0*tmp11*tmp7*Q0y*v0x*v0z+4.0*tmp11*tmp7*Q0z*v0x*v0y+4.0*tmp11*tmp1*tmp9*Q0w+4.0*tmp11*tmp1*tmp8*Q0w+4.0*tmp11*tmp1*Q0y*v0x*v0z+4.0*tmp11*tmp1*Q0z*v0x*v0y+4.0*tmp11*tmp16*tmp5*Q0w+4.0*tmp11*tmp16*tmp8*Q0w+4.0*tmp11*tmp16*Q0x*v0y*v0z+4.0*tmp11*tmp16*Q0z*v0x*v0y+4.0*tmp11*tmp14*tmp5*Q0w+4.0*tmp11*tmp14*tmp9*Q0w+4.0*tmp11*tmp14*Q0x*v0y*v0z+4.0*tmp11*tmp14*Q0y*v0x*v0z+4.0*tmp11*tmp5*Q0x*Q0y*Q0z+4.0*tmp11*tmp9*Q0x*Q0y*Q0z+4.0*tmp11*tmp8*Q0x*Q0y*Q0z+4.0*tmp11*Q0w*Q1w*Q1x*lf*v0z+4.0*tmp11*Q0w*Q1x*Q1y*lf*v0x+4.0*tmp11*Q0w*Q1y*Q1z*lf*v0z+4.0*tmp11*Q0x*Q1w*Q1x*lf*v0y+4.0*tmp11*Q0x*Q1y*Q1z*lf*v0y+4.0*tmp11*Q0y*Q1x*Q1y*lf*v0z+4.0*tmp11*Q0z*Q1w*Q1z*lf*v0y+8.0*tmp11*Q0w*Q0x*Q0y*v0x*v0y+8.0*tmp11*Q0w*Q0x*Q0z*v0x*v0z+8.0*tmp11*Q0w*Q0y*Q0z*v0y*v0z
    dpQ1x = -2.0*tmp11*tmp7*Q1x*lf*v0y+-2.0*tmp11*tmp1*Q1w*lf*v0z+-2.0*tmp11*tmp16*Q1w*lf*v0z+-2.0*tmp11*tmp16*Q1x*lf*v0y+-2.0*tmp11*tmp16*Q1y*lf*v0x+-2.0*tmp11*tmp14*Q1y*lf*v0x+-2.0*tmp11*tmp10*tmp13*Q1x+-2.0*tmp11*tmp17*tmp13*Q1x+-2.0*tmp11*P0y*Q1x*lf+-2.0*tmp11*P1x*Q1y*lf+-2.0*tmp11*P1z*Q1w*lf+-4.0*tmp11*tmp13*Q1w*Q1y*Q1z+-4.0*tmp11*Q0w*Q0y*Q1w*lf*v0x+-4.0*tmp11*Q0w*Q0z*Q1x*lf*v0x+-4.0*tmp11*Q0w*Q0z*Q1y*lf*v0y+-4.0*tmp11*Q0x*Q0y*Q1x*lf*v0x+-4.0*tmp11*Q0y*Q0z*Q1x*lf*v0z+2.0*tmp11*tmp7*Q1w*lf*v0z+2.0*tmp11*tmp7*Q1y*lf*v0x+2.0*tmp11*tmp1*Q1x*lf*v0y+2.0*tmp11*tmp1*Q1y*lf*v0x+2.0*tmp11*tmp14*Q1w*lf*v0z+2.0*tmp11*tmp14*Q1x*lf*v0y+2.0*tmp11*Q1x*Q1x*Q1x*tmp13+2.0*tmp11*tmp2*tmp13*Q1x+2.0*tmp11*P0x*Q1y*lf+2.0*tmp11*P0z*Q1w*lf+2.0*tmp11*P1y*Q1x*lf+4.0*tmp11*tmp10*tmp13*Q1x+4.0*tmp11*tmp17*tmp13*Q1x+4.0*tmp11*tmp13*Q1w*Q1y*Q1z+4.0*tmp11*Q0w*Q0x*Q1w*lf*v0y+4.0*tmp11*Q0w*Q0x*Q1x*lf*v0z+4.0*tmp11*Q0w*Q0y*Q1y*lf*v0z+4.0*tmp11*Q0x*Q0y*Q1y*lf*v0y+4.0*tmp11*Q0x*Q0z*Q1w*lf*v0x+4.0*tmp11*Q0x*Q0z*Q1y*lf*v0z+4.0*tmp11*Q0y*Q0z*Q1w*lf*v0y
    dpQ1y = -2.0*tmp11*tmp1*Q1y*lf*v0y+-2.0*tmp11*tmp1*Q1z*lf*v0z+-2.0*tmp11*tmp16*Q1x*lf*v0x+-2.0*tmp11*tmp16*Q1z*lf*v0z+-2.0*tmp11*tmp14*Q1x*lf*v0x+-2.0*tmp11*tmp14*Q1y*lf*v0y+-2.0*tmp11*tmp15*tmp13*Q1y+-2.0*tmp11*tmp2*tmp13*Q1y+-2.0*tmp11*P1x*Q1x*lf+-2.0*tmp11*P1y*Q1y*lf+-2.0*tmp11*P1z*Q1z*lf+-4.0*tmp11*tmp13*Q1w*Q1x*Q1z+-4.0*tmp11*Q0w*Q0x*Q1y*lf*v0z+-4.0*tmp11*Q0w*Q0y*Q1z*lf*v0x+-4.0*tmp11*Q0w*Q0z*Q1x*lf*v0y+2.0*tmp11*tmp7*Q1x*lf*v0x+2.0*tmp11*tmp7*Q1y*lf*v0y+2.0*tmp11*tmp7*Q1z*lf*v0z+2.0*tmp11*tmp1*Q1x*lf*v0x+2.0*tmp11*tmp16*Q1y*lf*v0y+2.0*tmp11*tmp14*Q1z*lf*v0z+2.0*tmp11*tmp10*tmp13*Q1y+2.0*tmp11*Q1y*Q1y*Q1y*tmp13+2.0*tmp11*P0x*Q1x*lf+2.0*tmp11*P0y*Q1y*lf+2.0*tmp11*P0z*Q1z*lf+4.0*tmp11*tmp15*tmp13*Q1y+4.0*tmp11*tmp2*tmp13*Q1y+4.0*tmp11*tmp13*Q1w*Q1x*Q1z+4.0*tmp11*Q0w*Q0x*Q1z*lf*v0y+4.0*tmp11*Q0w*Q0y*Q1x*lf*v0z+4.0*tmp11*Q0w*Q0z*Q1y*lf*v0x+4.0*tmp11*Q0x*Q0y*Q1x*lf*v0y+4.0*tmp11*Q0x*Q0y*Q1y*lf*v0x+4.0*tmp11*Q0x*Q0z*Q1x*lf*v0z+4.0*tmp11*Q0x*Q0z*Q1z*lf*v0x+4.0*tmp11*Q0y*Q0z*Q1y*lf*v0z+4.0*tmp11*Q0y*Q0z*Q1z*lf*v0y
    dpQ1z = -2.0*tmp11*tmp7*Q1w*lf*v0x+-2.0*tmp11*tmp7*Q1z*lf*v0y+-2.0*tmp11*tmp1*Q1w*lf*v0x+-2.0*tmp11*tmp1*Q1y*lf*v0z+-2.0*tmp11*tmp16*Q1y*lf*v0z+-2.0*tmp11*tmp16*Q1z*lf*v0y+-2.0*tmp11*tmp10*tmp13*Q1z+-2.0*tmp11*tmp17*tmp13*Q1z+-2.0*tmp11*P0x*Q1w*lf+-2.0*tmp11*P0y*Q1z*lf+-2.0*tmp11*P1z*Q1y*lf+-4.0*tmp11*tmp13*Q1w*Q1x*Q1y+-4.0*tmp11*Q0w*Q0y*Q1w*lf*v0z+-4.0*tmp11*Q0w*Q0y*Q1y*lf*v0x+-4.0*tmp11*Q0w*Q0z*Q1z*lf*v0x+-4.0*tmp11*Q0x*Q0y*Q1w*lf*v0y+-4.0*tmp11*Q0x*Q0y*Q1z*lf*v0x+-4.0*tmp11*Q0x*Q0z*Q1w*lf*v0z+-4.0*tmp11*Q0y*Q0z*Q1z*lf*v0z+2.0*tmp11*tmp7*Q1y*lf*v0z+2.0*tmp11*tmp1*Q1z*lf*v0y+2.0*tmp11*tmp16*Q1w*lf*v0x+2.0*tmp11*tmp14*Q1w*lf*v0x+2.0*tmp11*tmp14*Q1y*lf*v0z+2.0*tmp11*tmp14*Q1z*lf*v0y+2.0*tmp11*tmp15*tmp13*Q1z+2.0*tmp11*Q1z*Q1z*Q1z*tmp13+2.0*tmp11*P0z*Q1y*lf+2.0*tmp11*P1x*Q1w*lf+2.0*tmp11*P1y*Q1z*lf+4.0*tmp11*tmp10*tmp13*Q1z+4.0*tmp11*tmp17*tmp13*Q1z+4.0*tmp11*tmp13*Q1w*Q1x*Q1y+4.0*tmp11*Q0w*Q0x*Q1y*lf*v0y+4.0*tmp11*Q0w*Q0x*Q1z*lf*v0z+4.0*tmp11*Q0w*Q0z*Q1w*lf*v0y+4.0*tmp11*Q0x*Q0z*Q1y*lf*v0x+4.0*tmp11*Q0y*Q0z*Q1y*lf*v0y
    dpQ1w = -2.0*tmp11*tmp7*Q1z*lf*v0x+-2.0*tmp11*tmp1*Q1w*lf*v0y+-2.0*tmp11*tmp1*Q1x*lf*v0z+-2.0*tmp11*tmp1*Q1z*lf*v0x+-2.0*tmp11*tmp16*Q1x*lf*v0z+-2.0*tmp11*tmp14*Q1w*lf*v0y+-2.0*tmp11*tmp15*tmp13*Q1w+-2.0*tmp11*tmp2*tmp13*Q1w+-2.0*tmp11*P0x*Q1z*lf+-2.0*tmp11*P1y*Q1w*lf+-2.0*tmp11*P1z*Q1x*lf+-4.0*tmp11*tmp13*Q1x*Q1y*Q1z+-4.0*tmp11*Q0w*Q0x*Q1w*lf*v0z+-4.0*tmp11*Q0w*Q0y*Q1x*lf*v0x+-4.0*tmp11*Q0w*Q0y*Q1z*lf*v0z+-4.0*tmp11*Q0x*Q0y*Q1z*lf*v0y+-4.0*tmp11*Q0x*Q0z*Q1z*lf*v0z+2.0*tmp11*tmp7*Q1w*lf*v0y+2.0*tmp11*tmp7*Q1x*lf*v0z+2.0*tmp11*tmp16*Q1w*lf*v0y+2.0*tmp11*tmp16*Q1z*lf*v0x+2.0*tmp11*tmp14*Q1x*lf*v0z+2.0*tmp11*tmp14*Q1z*lf*v0x+2.0*tmp11*Q1w*Q1w*Q1w*tmp13+2.0*tmp11*tmp17*tmp13*Q1w+2.0*tmp11*P0y*Q1w*lf+2.0*tmp11*P0z*Q1x*lf+2.0*tmp11*P1x*Q1z*lf+4.0*tmp11*tmp15*tmp13*Q1w+4.0*tmp11*tmp2*tmp13*Q1w+4.0*tmp11*tmp13*Q1x*Q1y*Q1z+4.0*tmp11*Q0w*Q0x*Q1x*lf*v0y+4.0*tmp11*Q0w*Q0z*Q1w*lf*v0x+4.0*tmp11*Q0w*Q0z*Q1z*lf*v0y+4.0*tmp11*Q0x*Q0y*Q1w*lf*v0x+4.0*tmp11*Q0x*Q0z*Q1x*lf*v0x+4.0*tmp11*Q0y*Q0z*Q1w*lf*v0z+4.0*tmp11*Q0y*Q0z*Q1x*lf*v0y
    dpP0x = tmp11*tmp7*v0x+tmp11*tmp1*v0x+tmp11*P0x+-1.0*tmp11*tmp16*v0x+-1.0*tmp11*tmp14*v0x+-1.0*tmp11*P1x+-2.0*tmp11*Q0w*Q0z*v0y+-2.0*tmp11*Q1w*Q1z*lf+2.0*tmp11*Q0w*Q0y*v0z+2.0*tmp11*Q0x*Q0y*v0y+2.0*tmp11*Q0x*Q0z*v0z+2.0*tmp11*Q1x*Q1y*lf
    dpP0y = tmp11*tmp7*v0y+tmp11*tmp16*v0y+tmp11*tmp10*lf+tmp11*tmp17*lf+tmp11*P0y+-1.0*tmp11*tmp1*v0y+-1.0*tmp11*tmp14*v0y+-1.0*tmp11*tmp15*lf+-1.0*tmp11*tmp2*lf+-1.0*tmp11*P1y+-2.0*tmp11*Q0w*Q0x*v0z+2.0*tmp11*Q0w*Q0z*v0x+2.0*tmp11*Q0x*Q0y*v0x+2.0*tmp11*Q0y*Q0z*v0z
    dpP0z = tmp11*tmp7*v0z+tmp11*tmp14*v0z+tmp11*P0z+-1.0*tmp11*tmp1*v0z+-1.0*tmp11*tmp16*v0z+-1.0*tmp11*P1z+-2.0*tmp11*Q0w*Q0y*v0x+2.0*tmp11*Q0w*Q0x*v0y+2.0*tmp11*Q0x*Q0z*v0x+2.0*tmp11*Q0y*Q0z*v0y+2.0*tmp11*Q1w*Q1x*lf+2.0*tmp11*Q1y*Q1z*lf
    dpP1x = tmp11*tmp16*v0x+tmp11*tmp14*v0x+tmp11*P1x+-1.0*tmp11*tmp7*v0x+-1.0*tmp11*tmp1*v0x+-1.0*tmp11*P0x+-2.0*tmp11*Q0w*Q0y*v0z+-2.0*tmp11*Q0x*Q0y*v0y+-2.0*tmp11*Q0x*Q0z*v0z+-2.0*tmp11*Q1x*Q1y*lf+2.0*tmp11*Q0w*Q0z*v0y+2.0*tmp11*Q1w*Q1z*lf
    dpP1y = tmp11*tmp1*v0y+tmp11*tmp14*v0y+tmp11*tmp15*lf+tmp11*tmp2*lf+tmp11*P1y+-1.0*tmp11*tmp7*v0y+-1.0*tmp11*tmp16*v0y+-1.0*tmp11*tmp10*lf+-1.0*tmp11*tmp17*lf+-1.0*tmp11*P0y+-2.0*tmp11*Q0w*Q0z*v0x+-2.0*tmp11*Q0x*Q0y*v0x+-2.0*tmp11*Q0y*Q0z*v0z+2.0*tmp11*Q0w*Q0x*v0z
    dpP1z = tmp11*tmp1*v0z+tmp11*tmp16*v0z+tmp11*P1z+-1.0*tmp11*tmp7*v0z+-1.0*tmp11*tmp14*v0z+-1.0*tmp11*P0z+-2.0*tmp11*Q0w*Q0x*v0y+-2.0*tmp11*Q0x*Q0z*v0x+-2.0*tmp11*Q0y*Q0z*v0y+-2.0*tmp11*Q1w*Q1x*lf+-2.0*tmp11*Q1y*Q1z*lf+2.0*tmp11*Q0w*Q0y*v0x


    return c2,dpQ0x,dpQ0y,dpQ0z,dpQ0w,dpQ1x,dpQ1y,dpQ1z,dpQ1w,dpP0x,dpP0y,dpP0z,dpP1x,dpP1y,dpP1z

def locSpringSimple(Jb):
    
    Q0 = Jb.parent.Q
    Q1 = Jb.Q
    w0 = Jb.parent.w
    w1 = Jb.w  
    
    v0 = Jb.rest_p 
    
    P0 = Jb.parent.P
    P1 = Jb.P 
    lf = Jb.l*0.5

   
    ld = Jb.parent.P+Jb.parent.Q*v0 - Jb.P-Jb.Q*Vector((0,-lf,0))
    ww = 1.0/(w0+w1)
    if w0>0:
        Jb.parent.P+= -ld*w0*ww 
    Jb.P+= ld*w1*ww   

K1 = Matrix().to_3x3()
K2 = Matrix().to_3x3()

def locSpring(Jb):
    global K1
    global K2
    
    Q0 = Jb.parent.Q
    Q1 = Jb.Q
    w0 = Jb.parent.w
    w1 = Jb.w  
    
    v0 = Jb.rest_p 
    
    P0 = Jb.parent.P
    P1 = Jb.P 
    lf = Jb.l*0.5
    
    Jb.updateIW()
    Jb.parent.updateIW()
    
    connector0 = Jb.parent.P+Jb.parent.Q*v0
    connector1 = Jb.P+Jb.Q*Vector((0,-lf,0))
 
    computeMatrixK(connector0, w0, P0, Jb.parent.iIw, K1)
    computeMatrixK(connector1, w1, P1, Jb.iIw, K2)
    
    Kinv = (K1 + K2).inverted() 
   # print(connector0,connector1)
   # print(Jb.iIw)
    
    pt = Kinv * (connector1 - connector0)
    if (w0 != 0.0):
        r0 = connector0 - P0
        Jb.parent.P += w0*pt
        ot = (Jb.parent.iIw * (r0.cross(pt)))
        
        otQ = Quaternion()
        otQ.x =ot[0]
        otQ.y = ot[1]
        otQ.z =  ot[2]
        otQ.w = 0         
        Jb.parent.Q = qadd(Jb.parent.Q, otQ*Jb.parent.Q*0.5).normalized() 
        
    if (w1 != 0.0):
        r1 = connector1 - P1
        Jb.P += -w1*pt
        ot = (Jb.iIw * (r1.cross(-pt)))
        
        otQ = Quaternion()
        otQ.x =ot[0]
        otQ.y = ot[1]
        otQ.z =  ot[2]
        otQ.w = 0         
        Jb.Q = qadd(Jb.Q, otQ*Jb.Q*0.5).normalized() 
         
    
def quatSpringGradient1(Q0,Q1,r):

    Q0x = Q0.x
    Q0y = Q0.y
    Q0z = Q0.z
    Q0w = Q0.w
    Q1x = Q1.x
    Q1y = Q1.y
    Q1z = Q1.z
    Q1w = Q1.w
    rx  =  r.x
    ry  =  r.y
    rz  =  r.z
    rw  =  r.w
    
    
    tmp0 = -2.0*Q0w*Q0y
    tmp1 = -2.0*Q0w*Q0x
    tmp2 = -2.0*Q0w*Q1x
    tmp3 = -2.0*Q0w*Q1w
    tmp4 = -2.0*Q0x*Q1w
    tmp5 = -2.0*Q0y*Q1w
    tmp6 = -2.0*Q0w*Q0z
    tmp7 = -2.0*Q0x*Q1x
    tmp8 = -2.0*Q0x*Q0y
    tmp9 = -2.0*Q0y*Q1x
    tmp10 = -2.0*Q0x*Q0z
    tmp11 = -2.0*Q0z*Q1w
    tmp12 = -2.0*Q0z*Q1x
    tmp13 = -2.0*Q0y*Q0z
    tmp14 = 2.0*Q0w*Q1w
    tmp15 = 2.0*Q0x*Q1w
    tmp16 = 2.0*Q0w*Q0x
    tmp17 = 2.0*Q0w
    tmp18 = tmp17*Q0y
    tmp19 = tmp17*Q1x
    tmp20 = tmp17*Q0z
    tmp21 = 2.0*Q0y*Q1w
    tmp22 = 2.0*Q0y*Q1x
    tmp23 = 2.0*Q0x*Q0y
    tmp24 = 2.0*Q0x*Q1x
    tmp25 = 2.0*Q0z*Q1w
    tmp26 = 2.0*Q0x*Q0z
    tmp27 = 2.0*Q0y*Q0z
    tmp28 = 2.0*Q0z*Q1x
    tmp29 = -2.0*Q0w
    tmp30 = -2.0*Q0x
    tmp31 = -2.0*Q1x
    tmp32 = 2.0*Q0w*Q0w
    tmp33 = -2.0*Q0y
    tmp34 = -2.0*Q0z
    tmp35 = 2.0*Q0x*Q0x
    tmp36 = -2.0*Q1y
    tmp37 = -2.0*Q1z
    tmp38 = 2.0*Q0y*Q0y
    tmp39 = 2.0*Q0z*Q0z
    tmp40 = Q1w*Q1w
    tmp41 = -Q0x*Q1w-Q0y*Q1z+Q0w*Q1x+Q0z*Q1y-rx
    tmp42 = Q1x*Q1x
    tmp43 = 2.0*Q0x
    tmp44 = -Q0x*Q1y-Q0z*Q1w+Q0w*Q1z+Q0y*Q1x-rz
    tmp45 = 2.0*Q1w
    tmp46 = Q1y*Q1y
    tmp47 = -Q0y*Q1w-Q0z*Q1x+Q0w*Q1y+Q0x*Q1z-ry
    tmp48 = 2.0*Q0y
    tmp49 = Q0w*Q1w+Q0x*Q1x+Q0y*Q1y+Q0z*Q1z-rw
    tmp50 = 2.0*Q0z
    tmp51 = Q1z*Q1z
    c = tmp41*tmp41+tmp44*tmp44+tmp47*tmp47+tmp49*tmp49
    dQ0x = tmp3*Q1x+tmp29*Q1y*Q1z+tmp5*Q1z+tmp9*Q1y+tmp11*Q1y+tmp12*Q1z+tmp14*Q1x+tmp17*Q1y*Q1z+tmp21*Q1z+tmp22*Q1y+tmp25*Q1y+tmp28*Q1z+tmp31*rw+tmp37*ry+tmp43*tmp40+tmp43*tmp42+tmp43*tmp46+tmp43*tmp51+tmp45*rx+2.0*Q1y*rz
    dQ0y = tmp3*Q1y+tmp2*Q1z+tmp4*Q1z+tmp7*Q1y+tmp11*Q1x+tmp34*Q1y*Q1z+tmp14*Q1y+tmp19*Q1z+tmp15*Q1z+tmp24*Q1y+tmp25*Q1x+tmp50*Q1y*Q1z+tmp31*rz+tmp36*rw+tmp48*tmp40+tmp48*tmp42+tmp48*tmp46+tmp48*tmp51+tmp45*ry+2.0*Q1z*rx
    dQ0z = tmp3*Q1z+tmp2*Q1y+tmp4*Q1y+tmp7*Q1z+tmp5*Q1x+tmp33*Q1y*Q1z+tmp14*Q1z+tmp19*Q1y+tmp15*Q1y+tmp24*Q1z+tmp21*Q1x+tmp48*Q1y*Q1z+tmp36*rx+tmp37*rw+tmp50*tmp40+tmp50*tmp42+tmp50*tmp46+tmp50*tmp51+tmp45*rz+2.0*Q1x*ry
    dQ0w = tmp4*Q1x+tmp30*Q1y*Q1z+tmp5*Q1y+tmp9*Q1z+tmp11*Q1z+tmp12*Q1y+tmp15*Q1x+tmp43*Q1y*Q1z+tmp21*Q1y+tmp22*Q1z+tmp25*Q1z+tmp28*Q1y+-2.0*Q1w*rw+tmp31*rx+tmp36*ry+tmp37*rz+tmp17*tmp40+tmp17*tmp42+tmp17*tmp46+tmp17*tmp51
    dQ1x = tmp1*Q1w+tmp0*Q1z+tmp6*Q1y+tmp8*Q1y+tmp10*Q1z+tmp13*Q1w+tmp16*Q1w+tmp18*Q1z+tmp20*Q1y+tmp23*Q1y+tmp26*Q1z+tmp27*Q1w+tmp29*rx+tmp30*rw+tmp33*rz+tmp32*Q1x+tmp35*Q1x+tmp38*Q1x+tmp39*Q1x+tmp50*ry
    dQ1y = tmp1*Q1z+tmp0*Q1w+tmp6*Q1x+tmp8*Q1x+tmp10*Q1w+tmp13*Q1z+tmp16*Q1z+tmp18*Q1w+tmp20*Q1x+tmp23*Q1x+tmp26*Q1w+tmp27*Q1z+tmp29*ry+tmp33*rw+tmp34*rx+tmp32*Q1y+tmp35*Q1y+tmp38*Q1y+tmp39*Q1y+tmp43*rz
    dQ1z = tmp1*Q1y+tmp0*Q1x+tmp6*Q1w+tmp8*Q1w+tmp10*Q1x+tmp13*Q1y+tmp16*Q1y+tmp18*Q1x+tmp20*Q1w+tmp23*Q1w+tmp26*Q1x+tmp27*Q1y+tmp29*rz+tmp30*ry+tmp34*rw+tmp32*Q1z+tmp35*Q1z+tmp38*Q1z+tmp39*Q1z+tmp48*rx
    dQ1w = tmp1*Q1x+tmp0*Q1y+tmp6*Q1z+tmp8*Q1z+tmp10*Q1y+tmp13*Q1x+tmp16*Q1x+tmp18*Q1y+tmp20*Q1z+tmp23*Q1z+tmp26*Q1y+tmp27*Q1x+tmp29*rw+tmp32*Q1w+tmp35*Q1w+tmp38*Q1w+tmp39*Q1w+tmp43*rx+tmp48*ry+tmp50*rz

    return c, dQ0x,dQ0y,dQ0z,dQ0w,dQ1x,dQ1y,dQ1z,dQ1w
sqrt = math.sqrt

def quatSpringGradient2(Q0,Q1,r):
    
    Q0x = Q0.x
    Q0y = Q0.y
    Q0z = Q0.z
    Q0w = Q0.w
    Q1x = Q1.x
    Q1y = Q1.y
    Q1z = Q1.z
    Q1w = Q1.w
    rx  =  r.x
    ry  =  r.y
    rz  =  r.z
    rw  =  r.w
    
    
    tmp0 = sqrt(((((((((-(Q0x*Q1w)-(Q0y*Q1z))+(Q0w*Q1x))+(Q0z*Q1y))-rx)*((((-(Q0x*Q1w)-(Q0y*Q1z))+(Q0w*Q1x))+(Q0z*Q1y))-rx))+(((((-(Q0x*Q1y)-(Q0z*Q1w))+(Q0w*Q1z))+(Q0y*Q1x))-rz)*((((-(Q0x*Q1y)-(Q0z*Q1w))+(Q0w*Q1z))+(Q0y*Q1x))-rz)))+(((((-(Q0y*Q1w)-(Q0z*Q1x))+(Q0w*Q1y))+(Q0x*Q1z))-ry)*((((-(Q0y*Q1w)-(Q0z*Q1x))+(Q0w*Q1y))+(Q0x*Q1z))-ry)))+((((((Q0w*Q1w)+(Q0x*Q1x))+(Q0y*Q1y))+(Q0z*Q1z))-rw)*(((((Q0w*Q1w)+(Q0x*Q1x))+(Q0y*Q1y))+(Q0z*Q1z))-rw))))
    tmp1 = 1.0/tmp0*Q0w*Q0y
    tmp2 = 1.0/tmp0*Q0w*Q1x
    tmp3 = 1.0/tmp0*Q0w*Q0x
    tmp4 = 1.0/tmp0*Q0x*Q1w
    tmp5 = 1.0/tmp0*Q0w*Q1w
    tmp6 = 1.0/tmp0*Q0y*Q1w
    tmp7 = 1.0/tmp0*Q0w*Q0z
    tmp8 = 1.0/tmp0*Q0x*Q1x
    tmp9 = 1.0/tmp0*Q0y*Q1x
    tmp10 = 1.0/tmp0*Q0x*Q0y
    tmp11 = 1.0/tmp0*Q0x*Q0z
    tmp12 = 1.0/tmp0*Q0z*Q1w
    tmp13 = 1.0/tmp0*Q0z*Q1x
    tmp14 = 1.0/tmp0*Q0y*Q0z
    tmp15 = 1.0/tmp0*Q0w*Q0w
    tmp16 = Q1w*Q1w
    tmp17 = 1.0/tmp0*Q0x*Q0x
    tmp18 = Q1x*Q1x
    tmp19 = 1.0/tmp0*Q0y*Q0y
    tmp20 = 1.0/tmp0
    tmp21 = Q1y*Q1y
    tmp22 = tmp20*Q0z*Q0z
    tmp23 = Q1z*Q1z
    tmp24 = tmp20*Q0x
    tmp25 = tmp20*Q0y
    tmp26 = tmp4*Q1x
    tmp27 = tmp24*Q1y*Q1z
    tmp28 = tmp3*Q1y
    tmp29 = tmp20*Q0z
    tmp30 = tmp3*Q1z
    tmp31 = tmp3*Q1w
    tmp32 = tmp5*Q1y
    tmp33 = tmp5*Q1z
    tmp34 = tmp1*Q1z
    tmp35 = tmp5*Q1x
    tmp36 = tmp1*Q1x
    tmp37 = tmp1*Q1w
    tmp38 = tmp6*Q1y
    tmp39 = tmp7*Q1y
    tmp40 = tmp2*Q1z
    tmp41 = tmp7*Q1x
    tmp42 = tmp9*Q1z
    tmp43 = tmp2*Q1y
    tmp44 = tmp3*Q1x
    tmp45 = tmp7*Q1w
    tmp46 = tmp20*Q0w*Q1y*Q1z
    tmp47 = tmp10*Q1x
    tmp48 = tmp4*Q1z
    tmp49 = tmp10*Q1y
    tmp50 = tmp10*Q1w
    tmp51 = tmp6*Q1z
    tmp52 = tmp4*Q1y
    tmp53 = tmp1*Q1y
    tmp54 = tmp12*Q1z
    tmp55 = -Q0x*Q1w-Q0y*Q1z+Q0w*Q1x+Q0z*Q1y-rx
    tmp56 = tmp20*Q1w
    tmp57 = tmp11*Q1z
    tmp58 = tmp9*Q1y
    tmp59 = tmp7*Q1z
    tmp60 = tmp11*Q1x
    tmp61 = tmp8*Q1y
    tmp62 = tmp13*Q1y
    tmp63 = tmp11*Q1w
    tmp64 = tmp8*Q1z
    tmp65 = -Q0x*Q1y-Q0z*Q1w+Q0w*Q1z+Q0y*Q1x-rz
    tmp66 = tmp14*Q1y
    tmp67 = tmp14*Q1w
    tmp68 = tmp12*Q1y
    tmp69 = tmp14*Q1z
    tmp70 = tmp20*Q1x
    tmp71 = -Q0y*Q1w-Q0z*Q1x+Q0w*Q1y+Q0x*Q1z-ry
    tmp72 = tmp6*Q1x
    tmp73 = tmp10*Q1z
    tmp74 = tmp12*Q1x
    tmp75 = tmp20*Q1y
    tmp76 = Q0w*Q1w+Q0x*Q1x+Q0y*Q1y+Q0z*Q1z-rw
    tmp77 = tmp29*Q1y*Q1z
    tmp78 = tmp25*Q1y*Q1z
    tmp79 = tmp13*Q1z
    tmp80 = tmp11*Q1y
    tmp81 = tmp20*Q0w
    tmp82 = tmp20*Q1z
    tmp83 = tmp14*Q1x
    c = tmp0
    dQ0x = tmp35+tmp46+tmp51+tmp58+tmp68+tmp79+tmp24*tmp16+tmp24*tmp18+tmp24*tmp21+tmp24*tmp23+tmp56*rx+tmp75*rz-tmp35-tmp46-tmp51-tmp58-tmp68-tmp79-tmp70*rw-tmp82*ry
    dQ0y = tmp32+tmp40+tmp48+tmp61+tmp74+tmp77+tmp25*tmp16+tmp25*tmp18+tmp25*tmp21+tmp25*tmp23+tmp56*ry+tmp82*rx-tmp32-tmp40-tmp48-tmp61-tmp74-tmp77-tmp70*rz-tmp75*rw
    dQ0z = tmp33+tmp43+tmp52+tmp64+tmp72+tmp78+tmp29*tmp16+tmp29*tmp18+tmp29*tmp21+tmp29*tmp23+tmp56*rz+tmp70*ry-tmp33-tmp43-tmp52-tmp64-tmp72-tmp78-tmp75*rx-tmp82*rw
    dQ0w = tmp26+tmp27+tmp38+tmp42+tmp54+tmp62+tmp81*tmp16+tmp81*tmp18+tmp81*tmp21+tmp81*tmp23-tmp26-tmp27-tmp38-tmp42-tmp54-tmp62-tmp56*rw-tmp70*rx-tmp75*ry-tmp82*rz
    dQ1x = tmp31+tmp34+tmp39+tmp49+tmp57+tmp67+tmp15*Q1x+tmp17*Q1x+tmp19*Q1x+tmp22*Q1x+tmp29*ry-tmp31-tmp34-tmp39-tmp49-tmp57-tmp67-tmp81*rx-tmp24*rw-tmp25*rz
    dQ1y = tmp30+tmp37+tmp41+tmp47+tmp63+tmp69+tmp15*Q1y+tmp17*Q1y+tmp19*Q1y+tmp22*Q1y+tmp24*rz-tmp30-tmp37-tmp41-tmp47-tmp63-tmp69-tmp81*ry-tmp25*rw-tmp29*rx
    dQ1z = tmp28+tmp36+tmp45+tmp50+tmp60+tmp66+tmp15*Q1z+tmp17*Q1z+tmp19*Q1z+tmp22*Q1z+tmp25*rx-tmp28-tmp36-tmp45-tmp50-tmp60-tmp66-tmp81*rz-tmp24*ry-tmp29*rw
    dQ1w = tmp44+tmp53+tmp59+tmp73+tmp80+tmp83+tmp15*Q1w+tmp17*Q1w+tmp19*Q1w+tmp22*Q1w+tmp24*rx+tmp25*ry+tmp29*rz-tmp44-tmp53-tmp59-tmp73-tmp80-tmp83-tmp81*rw

    return c, dQ0x,dQ0y,dQ0z,dQ0w,dQ1x,dQ1y,dQ1z,dQ1w
    

def quatSpring(Jb,r=None,k=None):
    
    #Jb.Q = Jb.parent.Q*r
   # return 
    
    Q0 = Jb.parent.Q
    Q1 = Jb.Q
    w0 = Jb.parent.w
    w1 = Jb.w
    if(r==None):
        r = Jb.rest.to_quaternion()
    if(k==None):
        k = Jb.k 
        
    ra = Q0.inverted()*Q1
    if ra.dot(r) < 0:
        r = -r
        
    c, dQ0x,dQ0y,dQ0z,dQ0w,dQ1x,dQ1y,dQ1z,dQ1w = quatSpringGradient2(Q0,Q1,r)
    
    div = dQ0x*dQ0x*w0 + \
          dQ0y*dQ0y*w0 + \
          dQ0z*dQ0z*w0 + \
          dQ0w*dQ0w*w0 + \
          dQ1x*dQ1x*w1 + \
          dQ1y*dQ1y*w1 + \
          dQ1z*dQ1z*w1 + \
          dQ1w*dQ1w*w1 
    
    if(div> 1e-8):
        s = -c/div
        if(w0>0.0):
            #qadd2(Q0, Quaternion((dQ0w*s*w0*k,dQ0x*s*w0*k,dQ0y*s*w0*k,dQ0z*s*w0*k)))
            #
            Q0.x+=dQ0x*s*w0*k
            Q0.y+=dQ0y*s*w0*k
            Q0.z+=dQ0z*s*w0*k 
            Q0.w+=dQ0w*s*w0*k    
            Jb.parent.Q = Q0.normalized()    
        #qadd2(Q1, Quaternion((dQ1w*s*w1*k,dQ1x*s*w1*k,dQ1y*s*w1*k,dQ1z*s*w1*k)))
        Q1.x+=dQ1x*s*w1*k
        Q1.y+=dQ1y*s*w1*k
        Q1.z+=dQ1z*s*w1*k
        Q1.w+=dQ1w*s*w1*k
        Jb.Q = Q1.normalized()
        
def step(scene):
    global iters 
    global dt
    global ctx
    global cc
    #sub_steps = 1#*scene.jiggle.sub_steps)
    dt = 1.0/(scene.render.fps) 
  
    for o in scene.objects:
        if(o.type == 'ARMATURE'):
            
            arm = o.data
            ow = o.matrix_world.copy()
            scale =maxis(ow,0).length
            #normR(ow)
            iow = ow.inverted()
            iow3 = ow.to_3x3().inverted()
            
            i=0
            arm.jiggle.time_acc+= dt* arm.jiggle.fps 
            while arm.jiggle.time_acc > 1:
                arm.jiggle.time_acc-=1
            #for j in range( sub_steps):  
                bl = []
                wt = []
                        
                for b in o.pose.bones:
                    if(b.parent==None):
                        propB(ow,b,bl,None)
                hooks= []
                
                bl2 = []
                for wb in bl:
                    b = wb.b
                    wb.rest_w = b.bone.matrix_local.copy()*scale
                    saxis(wb.rest_w,3, maxis(wb.rest_w,3)+maxis(wb.rest_w,1)*b.bone.length*0.5)
                    
                for wb in bl:
                    b = wb.b
                    crest = b
                    # o------ -> ---o---
                    wb.restW = b.bone.matrix_local.copy() * scale #??
                  #  wb.Q = wb.Q.normalized()
                    
                    M = wb.M
                    if(b.bone.jiggle.enabled):
                        Jb = b.bone.jiggle
                        wb.X = wb.P = Jb.P
                        wb.R = wb.Q = Jb.R
                        wb.rest =  wb.rest_w #b.bone.matrix_local #
                        if(b.parent!=None):
                            wb.rest = wb.parent.rest_w.inverted()*wb.rest_w
                       # saxis(wb.rest,3, mpos(wb.rest)*scale)
                
                        wb.rest_base = b.bone.matrix_local
                        if(b.parent!=None):
                            wb.rest_base = b.parent.bone.matrix_local.inverted()*wb.rest_base
                            
                        wb.rest_p = wb.parent.rest_w.inverted()* (maxis(wb.rest_w,3)- maxis(wb.rest_w,1)*b.bone.length*0.5)# mpos(wb.rest)
                        wb.l = b.bone.length*scale
                        wb.w = 1.0/Jb.mass
                        wb.k = 1- pow(1-Jb.Ks, 1/scene.jiggle.iterations)
                        Jb.V*= 1.0-Jb.Kld
                        Jb.V+= scene.gravity*dt
                        Jb.W*= 1.0-Jb.Kd
                        qv = Quaternion()
                        qv.x =Jb.W[0]
                        qv.y = Jb.W[1]
                        qv.z = Jb.W[2]
                        qv.w = 0         
                                
                        wb.Q = qadd(wb.Q, qv*wb.Q*dt*0.5).normalized() 
                        
                        wb.P = wb.X + Jb.V*dt
                        wb.computeI()
                        
                        
                        if(Jb.control_bone!=""):
                            if(Jb.control_bone in o.pose.bones):
                                cb = o.pose.bones[Jb.control_bone]
                                clm = cb.matrix
                                if(cb.parent!=None):
                                    clm = cb.parent.matrix.inverted()*clm
                                wb.cQ = clm.to_quaternion().normalized()
                                wb.Kc = 1- pow(1-Jb.control, 1/scene.jiggle.iterations)
                            
                            
                            
                            
                        bl2.append(wb)
                      #  print(wb.Q)
                    else:
                        wb.w = 0
                        wb.X = wb.P = mpos(M)+maxis(M,1)*b.bone.length*0.5
                        wb.R = wb.Q = M.to_quaternion().normalized()
                        
                         
                                   
                for i in range(scene.jiggle.iterations):             
                    for wb in bl2:
                        b = wb.b        
                        if(b.parent==None):
                            continue                
                        Jb = b.bone.jiggle                        
                        locSpring(wb)                             
                    for wb in bl2:
                        b = wb.b        
                        if(b.parent==None):
                            continue                
                        Jb = b.bone.jiggle                     
                        quatSpring(wb,Jb.rest if Jb.use_custom_rest else wb.rest.to_quaternion().normalized())  
                        if(wb.cQ!=None):
                            quatSpring(wb, wb.cQ, wb.Kc)    
                for i in range(scene.jiggle.fix_iterations):               
                    for wb in bl2:
                        b = wb.b        
                        if(b.parent==None):
                            continue                
                        Jb = b.bone.jiggle                        
                       # locSpring(wb)      
                         
               # if FIX_DISPLACEMENT:
               #     for wb in bl2:
               #         if(wb.parent!=None):
               #             wb.P = wb.parent.P+wb.parent.Q.normalized()*wb.rest_p + maxis(wb.Q.normalized().to_matrix(),1)*b.bone.length*0.5*scale
               #     
                for wb in bl2:
                    b = wb.b
                    Jb = b.bone.jiggle
                    
                    wb.Q = wb.Q.normalized()
                    m = wb.Q.to_matrix()
                    for i in range(3):
                        for j in range(3):
                            wb.M[i][j] = m[i][j]*scale      
                    
                    #if FIX_DISPLACEMENT:
                     #   wb.P = wb.parent.P+wb.parent.Q*wb.rest_p + maxis(wb.M,1)*b.bone.length*0.5
                  #  R = Jb.R.to_matrix()       
                   #+  wb.P-wb.Q*Vector((0,1,0))*(b.bone.length*0.5*scale) - wb.parent.P+wb.parent.Q*wb.rest_p #tail
                    
                    Jb.V = (wb.P - wb.X)/dt
                    Jb.P = wb.P.copy()
                    qv = wb.Q*Jb.R.conjugated() #qadd(wb.Q,-Jb.R)*Jb.R.conjugated()#
                    Jb.W = Vector((qv.x,qv.y,qv.z))*(2/dt)
                    Jb.R = wb.Q
                    
                  
                                 
                    #
                    cp = Jb.P - maxis(wb.M,1)*b.bone.length*0.5
                    #if scene.jiggle.show_displacement:
                    #    cp = Jb.P - maxis(wb.M,1)*b.bone.length*0.5
                    #else:
                    #    cp =  wb.parent.P+wb.parent.Q.normalized()*wb.rest_p
                    #                 
                    wb.M[0][3]= cp[0]
                    wb.M[1][3]= cp[1]
                    wb.M[2][3]= cp[2]
               # for wb in bl:
                #   wb.R*=scale   
                
               # if scene.jiggle.fix_displacement:
                #    for wb in bl2:
               #         b = wb.b
                 #       Jb = b.bone.jiggle
                        
                                
                #if( FIX_DISPLACEMENT):
                #    for wb in bl2:
                #        b = wb.b
                #        Jb = b.bone.jiggle
                #        if(wb.parent!=None):
                #            Jb.P =  wb.P = wb.parent.P+wb.parent.Q*wb.rest_p + wb.Q*Vector((0,1,0))*(b.bone.length*0.5*scale)
                #            cp = Jb.P - maxis(wb.M,1)*b.bone.length*0.5
                #            wb.M[0][3]= cp[0]
                #            wb.M[1][3]= cp[1]
                #            wb.M[2][3]= cp[2]
                for wb in bl2:
                    b = wb.b
                    pM = ow
                    if(b.parent!=None):
                        pM = wb.parent.M
                    mb =  (pM*wb.rest_base).inverted()*wb.M
                    
                    #if SHOW_DISPLACEMENT:
                    #saxis(mb,3,Vector((0,0,0)))
                   # if not ENABLE_DISPLACEMENT:
                   #     saxis(mb,3,Vector((0,0,0)))
                    b.matrix_basis =mb
                   # print("im activ")
                   # b.matrix = iow*wb.M
                    
    scene.jiggle.last_frame+= 1
   #print("updated")
    
@persistent
def update_post(scene, tm = False):
    global iters 
    global dt
    global ctx
    global cc

#backing = False
@persistent
def update(scene, tm = False):
    global iters 
    global dt
    global ctx
    global cc 
    dt = 1.0/(scene.render.fps*scene.jiggle.sub_steps)
    if(not (scene.jiggle.test_mode or tm)):# or (backing and not tm)):
        return 
    step(scene)        

def bake(aa):
    print("bake all:",aa)
    global ctx
   # global backing 
    
    scene = bpy.context.scene
    scene.frame_set(scene.frame_start)
    
    for o in scene.objects:
        if(o.type == 'ARMATURE' and (o.select or aa)):
            
            
            arm = o.data
            ow = o.matrix_world
            scale = maxis(ow,0).length
            iow = ow.inverted()
            i=0
            for b in o.pose.bones:
                b.bone.select = (b.bone.select or aa) and b.bone.jiggle.enabled
                if(b.bone.select or aa and b.bone.jiggle.enabled):                    
                    M = ow*b.matrix #ow*Sbp.wmat* Sb.rmat #im
                    #l ,r,s = M.decompose()
                    
                    Jb = b.bone.jiggle
                    setq(Jb.R, M.to_quaternion().normalized())
                    Jb.V = Vector((0,0,0))
                    Jb.P = mpos(M)
                    Jb.W = Vector((0,0,0))
    
    ltm = scene.jiggle.test_mode                
    scene.jiggle.test_mode = False
   # backing = True
    for i in range(scene.frame_start, scene.frame_end):
        scene.frame_set(i)
        update(scene,tm=True)
        print("frame: ",i)
        for o in scene.objects:
            if( (o.select or aa) and o.type == 'ARMATURE' ):
                scene.objects.active = o
                m = o.mode == 'POSE'
                
                if(not m):
                    bpy.ops.object.posemode_toggle()
                
                bpy.ops.anim.keyframe_insert_menu(type='LocRotScale')
        
                if(not m):
                    bpy.ops.object.posemode_toggle()
    #backing= False
    scene.jiggle.test_mode = ltm
                   
class BakeOperator(bpy.types.Operator):
    a = bpy.props.BoolProperty()
    bl_idname = "jiggle.bake"
    bl_label = "Bake Animation"
    def execute(self, context):
        bake(self.a)           
        return {'FINISHED'}    

def register():
    bpy.app.handlers.frame_change_pre.append(update) 
    
    bpy.utils.register_class(JiggleArmaturePanel)
    bpy.utils.register_class(JiggleScene)
    bpy.types.Scene.jiggle = bpy.props.PointerProperty(type = JiggleScene)
    bpy.utils.register_class(JiggleScenePanel)
    bpy.utils.register_class(JiggleArmature)
    bpy.utils.register_class(JiggleBone) 
    bpy.utils.register_class(BakeOperator)

    bpy.types.Armature.jiggle = bpy.props.PointerProperty(type = JiggleArmature)
    bpy.types.Bone.jiggle = bpy.props.PointerProperty(type = JiggleBone)
    bpy.types.Bone.jiggle_tp = bpy.props.FloatProperty()
    bpy.utils.register_class(ResetJigglePropsOperator)

    bpy.utils.register_class(SetRestJigglePropsOperator)
    bpy.utils.register_class(JiggleBonePanel)
def unregister():    

    bpy.utils.unregister_class(JiggleScene)
    bpy.utils.unregister_class(JiggleScenePanel)
    bpy.utils.unregister_class(JiggleBone)
    bpy.utils.unregister_class(JiggleArmature)
    bpy.utils.unregister_class(SetRestJigglePropsOperator)
    bpy.utils.unregister_class(ResetJigglePropsOperator)
    bpy.utils.unregister_class(JiggleBonePanel) 
    bpy.utils.unregister_class(BakeOperator)
    bpy.utils.unregister_class(JiggleArmaturePanel)
    bpy.app.handlers.frame_change_pre.remove(update) 
if __name__ == '__main__':
	register()
