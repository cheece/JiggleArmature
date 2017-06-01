#Copyright (c) 2017 Simón Flores (https://github.com/cheece)

#Permission is hereby granted, free of charge, 
#to any person obtaining a copy of this software 
#and associated documentation files (the "Software"), 
#to deal in the Software without restriction, 
#including without limitation the rights to use, 
#copy, modify, merge, publish, distribute, sublicense,
#and/or sell copies of the Software, and to permit 
#persons to whom the Software is furnished to do so,
#subject to the following conditions:The above copyright 
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

#to enable jiggle physics on a bone just click on jigglebone 
 
bl_info = {
    "name": "Jiggle Armature",
    "author": "Simón Flores",
    "version": (2, 0, 0),
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


    
    
class JiggleScene(bpy.types.PropertyGroup):
    test_mode = bpy.props.BoolProperty(default=True)
    sub_steps = bpy.props.IntProperty(min=1, default = 2)
    iterations = bpy.props.IntProperty(min=1, default = 2)
    last_frame = bpy.props.IntProperty()
    length_fix_iters = bpy.props.IntProperty(min=0, default = 2)


bpy.utils.register_class(JiggleScene)
bpy.types.Scene.jiggle = bpy.props.PointerProperty(type = JiggleScene)

class JiggleScenePanel(bpy.types.Panel):
    bl_idname = "Scene_PT_jiggle"
    bl_label = "Jiggle Scene"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "scene"
    bl_options = {'DEFAULT_CLOSED'}


    def draw(self, context):
        layout = self.layout
        col = layout.column()
        col.prop(context.scene.jiggle,"test_mode")
        col.prop(context.scene.jiggle,"sub_steps")
        col.prop(context.scene.jiggle,"iterations")
        #col.prop(context.scene.jiggle,"length_fix_iters")
        col.operator("jiggle.bake")



bpy.utils.register_class(JiggleScenePanel)

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
    enabled = bpy.props.BoolProperty(default=False)
    Kd = bpy.props.FloatProperty(name = "damping",min=0.0, default = 0.1, update = funp("Kd"))
    Ks = bpy.props.FloatProperty(name = "linear spring force",min=0.0 ,default = 100, update = funp("Ks"))
    Kr = bpy.props.FloatProperty(name = "angular spring force",min=0.0, default = 100)
    mass = bpy.props.FloatProperty(min=0.0001, default = 1.0, update = funp("mass"))  
    #M = bpy.props.FloatVectorProperty(size=9,subtype='MATRIX')    
    R = bpy.props.FloatVectorProperty(size=4,subtype='QUATERNION')
    V = bpy.props.FloatVectorProperty(size=3,subtype='XYZ')
    P = bpy.props.FloatVectorProperty(size=3,subtype='XYZ')
  
bpy.utils.register_class(JiggleBone)

bpy.types.Bone.jiggle = bpy.props.PointerProperty(type = JiggleBone)

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
                iow = ow.inverted()
                i=0
                for b in o.pose.bones:
                    if(b.bone.select):                    
                        M = ow*b.matrix #ow*Sbp.wmat* Sb.rmat #im
                        #l ,r,s = M.decompose()
                        
                        Jb = b.bone.jiggle
                        setq(Jb.R, M.to_quaternion())
                        l = b.bone.length
                        tg = Vector((M[0][1]*l+M[0][3],M[1][1]*l+M[1][3],M[2][1]*l+M[2][3]))#(M[0][1]+M[0][3],M[1][1]+M[1][3],M[2][1]+M[2][3]))
                        
                        Jb.V = Vector((0,0,0))
                        Jb.P = tg 
                        #Jb.M = Matrix(M)
        return {'FINISHED'}
    
bpy.utils.register_class(ResetJigglePropsOperator)

class JiggleBonePanel(bpy.types.Panel):
    bl_idname = "Bone_PT_jiggle_bone"
    bl_label = "Jiggle Bone"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "bone"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return (context.bone is not None and context.object is not None and context.object.type == 'ARMATURE')

    def draw_header(self, context):
        layout = self.layout
        bon = context.bone
        layout.prop(bon.jiggle, "enabled", text="")

    def draw(self, context):
        layout = self.layout

        bon = context.bone
        col = layout.column()
        if(bon.jiggle.enabled):
            col.prop(bon.jiggle,"Ks")          
            col.prop(bon.jiggle,"Kd")
            col.prop(bon.jiggle,"mass")
        col.operator("jiggle.reset")
            
bpy.utils.register_class(JiggleBonePanel)


class JB:
    def __init__(self, b,w,p):
        self.w  = w
        self.b = b
        self.parent = p
        
def propB(ow,b, l, p):
    j = JB(b, ow*b.matrix, p)
    l.append(j)
    for c in b.children:
        propB(ow,c,l,j)
def mpos(M):
    return Vector((M[0][3],M[1][3],M[2][3]))      
def step(scene):
    global iters 
    global dt
    global ctx
    global cc
    dt = 1.0/(scene.render.fps*scene.jiggle.sub_steps)
  
    for o in scene.objects:
        if(o.type == 'ARMATURE'):
            bl = []
            wt = []
            
            arm = o.data
            ow = o.matrix_world
            iow = ow.inverted()
            iow3 = ow.to_3x3().inverted()
            i=0
            for b in o.pose.bones:
                if(b.parent==None):
                    propB(ow,b,bl,None)
            for j in range( scene.jiggle.sub_steps):  
                    
                for i in range(scene.jiggle.length_fix_iters):            
                    for wb in bl:
                        b = wb.b
                        if(b.bone.jiggle.enabled):
                            Jb = b.bone.jiggle
                            O = mpos(wb.w)
                            if(b.parent.bone.jiggle.enabled):
                                O = b.parent.bone.jiggle.P
                            OP = Jb.P -O    
                            l = b.bone.length
                            dP = -OP*(OP.length - l)/OP.length
                            w0 = 0
                            w1 = 1
                            if(b.parent.bone.jiggle.enabled):
                               w0 = 1.0/b.parent.bone.jiggle.mass
                               w1 = Jb.mass
                            Jb.P += dP*w1/(w0 + w1)
                            if(b.parent.bone.jiggle.enabled):
                                b.parent.bone.jiggle.P += -dP*w0/(w0 + w1)    
                                                 
                for wb in bl:
                    b = wb.b
                    if(b.bone.jiggle.enabled):
                        tgr = (wb.parent.w * (b.bone.parent.matrix_local.inverted() * b.bone.matrix_local)).to_3x3()
                        Jb = b.bone.jiggle
                        #print("m: " + str(Jb.M) )
                        #Jb.V = Vector((1,0,0))
                        R = Jb.R.to_matrix()
                        R2 = iow3* R                  
                        M = wb.w 
                        l = b.bone.length
                        P2 =  Vector((R[0][1]*l+M[0][3],R[1][1]*l+M[1][3],R[2][1]*l+M[2][3]))#(M[0][1]+M[0][3],M[1][1]+M[1][3],M[2][1]+M[2][3]))
                        op = Vector((M[0][3],M[1][3],M[2][3]))
                        ax = P2 - op
                        tx = Jb.P - op
                        ax.normalize()
                        tx.normalize()
                        cosa = min(1,max(-1,ax.dot(tx)))
                        rax = ax.cross(tx)
                        #tq =  (P2- Jb.P).cross(Jb.P - op))
                        print(cosa)
                        rr = mathutils.Matrix.Rotation(math.acos(cosa),3, rax)
                        R = rr*R
                        #R = R + skew(-tq)*R*dt
                        #tv = (P2-Jb.P)/dt    
                        
                        skm = tgr*R.inverted() - Matrix.Identity(3)
                        chv = iskew(skm)
                        print("skm: " + str(skm) +  ":V:" + str(chv))
                        Jb.V+= -min(1.0,Jb.Kd*dt)*Jb.V
                        
                        #gtq = #gravity torque
                        
                        Jbp = b.parent.bone.jiggle
                        F = chv*Jb.Ks
                        if(Jbp.enabled):
                            Jbp.V += -F* (1.0/Jbp.mass) * dt 
                        Jb.V += F * (1.0/Jb.mass) *dt
                        Jb.V += -(Jb.P - op).cross(scene.gravity) * dt 
                            
                         #Jb.V += -tq/dt
                        #Jb.V += -tv*dt
                        #print("dm:" + str(dM))
                        #R = R + dM
                        
                        setq(Jb.R,R.to_quaternion())
                        Jb.R.normalize()
                        R = Jb.R.to_matrix()
                        setm(wb.w,R)
                        #m = b.matrix.copy()
                        #setm(m, iow3* R)
                        #b.matrix = m #, Jb.M.to_4x4()*iow)
                        #
                        print("matrix changed")
                for wb in bl:
                    b = wb.b
                    if(b.bone.jiggle.enabled):
                        Jb = b.bone.jiggle
                        R = Jb.R.to_matrix()
                        dM = skew(Jb.V)*R*dt
                        R = R + dM
                        setq(Jb.R,R.to_quaternion())
                        Jb.R.normalize()
                        R = Jb.R.to_matrix()
                        setm(wb.w,R)
                        M = wb.w 
                        Jb.P = Vector((R[0][1]*l+M[0][3],R[1][1]*l+M[1][3],R[2][1]*l+M[2][3]))#(M[0][1]+M[0][3],M[1][1]+M[1][3],M[2][1]+M[2][3]))
                                
            for wb in bl:
                b = wb.b
                if(b.bone.jiggle.enabled):  
                    b.matrix = iow*wb.w  
                    
    scene.jiggle.last_frame+= 1
    print("updated")
@persistent
def update(scene, tm = False):
    global iters 
    global dt
    global ctx
    global cc
    dt = 1.0/(scene.render.fps*scene.jiggle.sub_steps)
    print(scene.jiggle.test_mode)
    if(not (scene.jiggle.test_mode or tm)):
        return
   # if(scene.frame_current == scene.jiggle.last_frame):
   #     return        
    print("beg2 " + str(scene.frame_current)+ " " +  str(scene.jiggle.last_frame))
    if(scene.frame_current <  scene.jiggle.last_frame or scene.frame_current == scene.frame_start): #frame break
        scene.jiggle.last_frame = scene.frame_current
        for o in scene.objects:
            if( o.type == 'ARMATURE'):
                arm = o.data
                ow = o.matrix_world
                iow = ow.inverted()
                i=0
                for b in o.pose.bones:
                    if(b.bone.jiggle.enabled):
                        M = ow*b.matrix #ow*Sbp.wmat* Sb.rmat #im

                        Jb = b.bone.jiggle
                        setq(Jb.R, M.to_quaternion())
                        l = b.bone.length
                        tg = Vector((M[0][1]*l+M[0][3],M[1][1]*l+M[1][3],M[2][1]*l+M[2][3]))#(M[0][1]+M[0][3],M[1][1]+M[1][3],M[2][1]+M[2][3]))
                        Jb.V = Vector((0,0,0))
                        Jb.P = tg 
        if(scene.frame_current <= bpy.context.scene.frame_start):
            return
    nframes = scene.frame_current - scene.jiggle.last_frame
    for i in range(nframes):     
        #
        step(scene)        

bpy.app.handlers.frame_change_post.append(update) 
