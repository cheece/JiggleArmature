#Copyright (c) 2016 Simón Flores (https://github.com/cheece)

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

#to enable jiggle physics on a bone first enable Jiggle on the
#armature it belongs to and then enable Jiggle Bone on the bone 
 
bl_info = {
    "name": "Jiggle Armature",
    "author": "Simón Flores",
    "version": (0, 1, 0),
    "blender": (2, 73, 0),
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
class JiggleArmature(bpy.types.PropertyGroup):
    enabled = bpy.props.BoolProperty(default=False)
    iterations = bpy.props.IntProperty(min=1, default = 2)
    test_mode = bpy.props.BoolProperty(default=True)
    
class JiggleBone(bpy.types.PropertyGroup):
    enabled = bpy.props.BoolProperty(default=False)
    Kd = bpy.props.FloatProperty(name = "damping",min=0.0, max=1.0, default = 0.1)
    Kl = bpy.props.FloatProperty(name = "length conservation",min=0.0, max=1.0, default = 0.5)
    Ks = bpy.props.FloatProperty(name = "shape conservation",min=0.0, max=1.0, default = 0.5)
    Kv = bpy.props.FloatProperty(name = "volume conservation",min=0.0, max=1.0, default = 0.5)
    mass = bpy.props.FloatProperty(min=0.0001, default = 1.0)
    
    X = bpy.props.FloatVectorProperty(size=3,subtype='XYZ')
    V = bpy.props.FloatVectorProperty(size=3,subtype='XYZ')
    P = bpy.props.FloatVectorProperty(size=3,subtype='XYZ')   
    
    index = bpy.props.IntProperty()
    max_deformation = bpy.props.FloatProperty(default = 10.2)
    debug = bpy.props.StringProperty(name = "debug object")

dt = 1.0/24.0
iters = 1

def getS(b):
    return ctx.states[b.bone.jiggle.index]

def applyConstraint(b,C,dC0, dC1, K):
    K1 =  1- pow(1-max(min(1,K),0),1.0/iters)   
    Jb = b.bone.jiggle
    w1  = 1.0/Jb.mass
    w0 = 0        
    Jbp = None
    if(b.parent != None and b.parent.bone.jiggle.enabled):
        w0 = 1.0/b.parent.bone.jiggle.mass
        Jbp = b.parent.bone.jiggle 
    div =     (dC0.dot(dC0)*w0 +dC1.dot(dC1)*w1 )  
    if(div > 0.0001):
        s = C/div
        Jb.P+= -dC1*s*w1*K1
        if(Jbp!=None):
            Jbp.P+= -dC0*s*w0*K1
def updateMat(ow, iow,b):

    Jb = b.bone.jiggle  
    Jbp = b.parent.bone.jiggle
    Sb = getS(b)
    Sbp = getS(b.parent)
        
    l = b.bone.length
    
    par = b.parent
    aM = ow*Sbp.wmat* Sb.rmat
    
    O = Vector((aM[0][3] ,aM[1][3] ,aM[2][3] ))
    V = Jb.P- O #Vector((Jb.P.x-O.x ,Jb.P.y-O.y ,Jb.P.z-O.z ))
    lv = V.length
    #if(lv > (Jb.max_deformation+1)*l):
        #V = V*(Jb.max_deformation+1)*l/lv
        #Jb.P = O+ V
    nV = V/lv
    cur = Vector((aM[0][1],aM[1][1],aM[2][1]))
    ncur = cur.normalized()
    cos = nV.dot(ncur)
    axis = nV.cross(ncur)
    if(cos>1):
        cos=1
    if(cos<-1):
        cos = -1    
    
    la = axis.length
    if(la>0):
        axis/=la#.normalized()
        #print(nV, ncur,cos)
        aM = Matrix.Rotation(-math.acos(cos),4, axis)*aM
        
        
        aM[0][3] = O.x #Jb.X.x-aM[0][3]    
        aM[1][3] = O.y #Jb.X.y-aM[1][3]
        aM[2][3] = O.z #Jb.X.z-aM[2][3]
        cur = getAxis(aM,1)
        cur *= (lv/l)/cur.length
        
        aM[0][1] = cur.x #Jb.X.x-aM[0][3]    
        aM[1][1] = cur.y #Jb.X.y-aM[1][3]
        aM[2][1] = cur.z #Jb.X.z-aM[2][3]
        
        
        Sb.wmat = iow*aM
        #scene = bpy.context.scene
       # b.scale.y = V.length/l    
    
  
def getAxis(M,i):
    return   Vector((M[0][i],M[1][i],M[2][i]))


def resetBone(ow, iow,b):
    par = b.parent
    im = par.bone.matrix_local.inverted()* b.bone.matrix_local
    M = ow*par.matrix* im    
    l = b.bone.length    
    tg = Vector((M[0][1]*l+M[0][3],M[1][1]*l+M[1][3],M[2][1]*l+M[2][3]))#(M[0][1]+M[0][3],M[1][1]+M[1][3],M[2][1]+M[2][3]))
    
    Jb = b.bone.jiggle
    Jb.X=tg
    Jb.P = tg
    Jb.V= Vector((0,0,0))
def updateBone(ow, iow,b):
    par = b.parent
    Jb = b.bone.jiggle
    Jbp = b.parent.bone.jiggle
    Sb = getS(b)
    Sbp = getS(b.parent)
    
    #im = par.bone.matrix_local.inverted()* b.bone.matrix_local
    M = ow*Sbp.wmat* Sb.rmat #im
    
    l = b.bone.length
    
    tg = Vector((M[0][1]*l+M[0][3],M[1][1]*l+M[1][3],M[2][1]*l+M[2][3]))#(M[0][1]+M[0][3],M[1][1]+M[1][3],M[2][1]+M[2][3]))
    
    aM = M.copy() #ow*b.matrix.copy()
    


    N = getAxis(aM, 1)
    Nl = N*l
    X0 = getAxis(aM, 3)
    X1 = Jb.P
    X01 = X1-X0
    Xl1 = X1-X0-Nl
    applyConstraint(b,X01.dot(N)-l,-N,N,Jb.Kv)
    if(par.bone.jiggle.enabled):
        updateMat(ow,iow,par)           
    updateMat(ow,iow,b)
    X1 = Jb.P
    X01 = X1-X0
    Xl1 = X1-X0-Nl
    applyConstraint(b,X01.dot(X01)-l*l,-X01*2,X01*2,Jb.Kl) 
    if(par.bone.jiggle.enabled):
        updateMat(ow,iow,par)           
    updateMat(ow,iow,b)
    X1 = Jb.P
    X01 = X1-X0
    Xl1 = X1-X0-Nl
    applyConstraint(b,(Xl1).dot(Xl1),-Xl1*2,Xl1*2,Jb.Ks) 
    if(par.bone.jiggle.enabled):
        updateMat(ow,iow,par)           
    updateMat(ow,iow,b)

    
    
    
def sortbb(b, ol):
    if(b.bone.jiggle.enabled):
        ol.append(b)
    for c in b.children:
        sortbb(c,ol)
def sortb(o, ol):
    for b in o.pose.bones:
        if(b.parent==None):
            sortbb(b,ol)
class BoneState:
    def __init__(self):
        self.wmat = None
        self.rmat = None       
        
class Ctx:    
    def __init__(self):
        self.states = [] 
ctx = None    

def bake():
    scene = bpy.context.scene
    scene.frame_set(scene.frame_start)
    for o in scene.objects:
        if(o.type == 'ARMATURE' and o.data.jiggle.enabled):
            arm = o.data
            ow = o.matrix_world
            iow = ow.inverted() 
            for b in  o.pose.bones:
                b.bone.select = b.bone.jiggle.enabled
                if(b.bone.jiggle.enabled):  
                    arm.bones.active = b.bone        
                    Jb = b.bone.jiggle             
                    resetBone(ow,iow,b)
                
    for i in range(scene.frame_start, scene.frame_end):
        scene.frame_set(i)
        update(scene,tm=True)
        
        for o in scene.objects:
            if(o.type == 'ARMATURE' and o.data.jiggle.enabled):
                scene.objects.active = o
                m = o.mode == 'POSE'
                
                if(not m):
                    bpy.ops.object.posemode_toggle()
                
                bpy.ops.anim.keyframe_insert_menu(type='LocRotScale')
        
                if(not m):
                    bpy.ops.object.posemode_toggle()
                

                
        
        
class BakeOperator(bpy.types.Operator):
    bl_idname = "jiggle.bake"
    bl_label = "Bake Animation"
    def execute(self, context):
        bake()           
        return {'FINISHED'}    
bpy.utils.register_class(BakeOperator)
def update(scene, tm = False):
    global iters 
    global dt
    global ctx
    dt = 1.0/scene.render.fps
    #print("beg")
    for o in scene.objects:
        if(o.type == 'ARMATURE' and o.data.jiggle.enabled and (o.data.jiggle.test_mode or tm)):
            arm = o.data
            ow = o.matrix_world
            iow = ow.inverted()
            iters = o.data.jiggle.iterations
            ol = []            
            sortb(o,ol)
            ctx = Ctx()
            i=0
            for b in o.pose.bones:
                s = BoneState()
                s.wmat = b.matrix
                M = b.bone.matrix_local
                if(b.parent!=None):
                    M = b.parent.bone.matrix_local.inverted()*M                    
                s.rmat = M     
                ctx.states.append(s)  
                b.bone.jiggle.index = i   
                i+=1      
                
            for b in  ol:            
                Jb = b.bone.jiggle
                Jb.V+= scene.gravity*dt
                Jb.V*= (1.0-Jb.Kd)
                Jb.P = Jb.X+Jb.V*dt                
                if(bpy.context.scene.frame_current == bpy.context.scene.frame_start):
                    resetBone(ow,iow,b)
                    
            for i in range(iters):
                for b in ol:
                    #print(b.name)    
                    updateBone(ow,iow,b)
 
            for b in ol:
                Jb = b.bone.jiggle
                Jb.V= (Jb.P-Jb.X)/dt
                Jb.X = Jb.P      
                Sb = getS(b)
                Sbp = getS(b.parent)
                b.matrix_basis = Sb.rmat.inverted()*Sbp.wmat.inverted()*Sb.wmat
                               
                if(Jb.debug in scene.objects):
                    scene.objects[Jb.debug].location = Jb.X
class JiggleBonePanel(bpy.types.Panel):
    bl_idname = "Bone_PT_jiggle_bone"
    bl_label = "Jiggle Bone"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "bone"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return (context.object is not None and context.object.type == 'ARMATURE' and context.object.data.jiggle.enabled)

    def draw_header(self, context):
        layout = self.layout
        bon = context.bone
        layout.prop(bon.jiggle, "enabled", text="")

    def draw(self, context):
        layout = self.layout

        bon = context.bone
        if(bon.jiggle.enabled):
            col = layout.column()
            col.prop(bon.jiggle,"Kv")   
            col.prop(bon.jiggle,"Kl")   
            col.prop(bon.jiggle,"Ks")          
            col.prop(bon.jiggle,"Kd")
            col.prop(bon.jiggle,"mass")
            col.prop(bon.jiggle,"debug")
            col.prop(bon.jiggle,"max_deformation")


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

    def draw_header(self, context):
        layout = self.layout
        arm = context.object.data
        layout.prop(arm.jiggle, "enabled", text="")

    def draw(self, context):
        layout = self.layout
        arm = context.object.data
        if(arm.jiggle.enabled):
            col = layout.column()
            col.prop(arm.jiggle,"iterations")
            col.prop(arm.jiggle,"test_mode")
            col.operator("jiggle.bake")




def register():
    bpy.utils.register_class(JiggleArmature)
    bpy.types.Armature.jiggle = bpy.props.PointerProperty(type = JiggleArmature)
    bpy.utils.register_class(JiggleBone)
    bpy.types.Bone.jiggle = bpy.props.PointerProperty(type = JiggleBone)
    bpy.utils.register_class(JiggleBonePanel)
    bpy.utils.register_class(JiggleArmaturePanel)
    bpy.app.handlers.frame_change_post.append(update) 
def unregister():    
    bpy.utils.unregister_class(JiggleArmature)
    bpy.utils.unregister_class(JiggleBone)
    bpy.utils.unregister_class(JiggleBonePanel)
    bpy.utils.unregister_class(JiggleArmaturePanel)
   
if __name__ == '__main__':
	register()

