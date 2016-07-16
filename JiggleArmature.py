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
    "version": (0, 2, 1),
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
from bpy.app.handlers import persistent

class JiggleScene(bpy.types.PropertyGroup):
    test_mode = bpy.props.BoolProperty(default=True)
    sub_steps = bpy.props.IntProperty(min=1, default = 2)

class JiggleArmature(bpy.types.PropertyGroup):
    enabled = bpy.props.BoolProperty(default=False)
    iterations = bpy.props.IntProperty(min=1, default = 2)
    
class JiggleBone(bpy.types.PropertyGroup):
    enabled = bpy.props.BoolProperty(default=False)
    Kd = bpy.props.FloatProperty(name = "damping",min=0.0, max=1.0, default = 0.1)
    Kl = bpy.props.FloatProperty(name = "length conservation",min=0.0, max=1.0, default = 0.5)
    Ks = bpy.props.FloatProperty(name = "shape conservation",min=0.0, max=1.0, default = 0.5)
    Kv = bpy.props.FloatProperty(name = "volume conservation",min=0.0, max=1.0, default = 0.5)
    mass = bpy.props.FloatProperty(min=0.0001, default = 1.0)  
   # Sv = bpy.props.FloatProperty(name = "bone volume conservation",min=0.0, max=1.0, default = 0.0)
    
       
    X = bpy.props.FloatVectorProperty(size=3,subtype='XYZ')
    V = bpy.props.FloatVectorProperty(size=3,subtype='XYZ')
    P = bpy.props.FloatVectorProperty(size=3,subtype='XYZ')   
    
    M = bpy.props.FloatVectorProperty(size=16,subtype='MATRIX')   
    rotV = bpy.props.FloatProperty(default=0.0)
    Kr = bpy.props.FloatProperty(name = "rotation conservation",min=0.0, max=1.0, default = 1.0)
    
    
    index = bpy.props.IntProperty()
    max_deformation = bpy.props.FloatProperty(default = 10.2)
    debug = bpy.props.StringProperty(name = "debug object")

dt = 1.0/24.0
iters = 1
def copyObject(A, B):
    A.enabled = B.enabled
    A.Kd   = B.Kd
    A.Kl   = B.Kl  
    A.Ks   = B.Ks  
    A.Kv   = B.Kv  
    A.mass = B.mass           
    A.X    = B.X   
    A.V    = B.V   
    A.P    = B.P   
    A.max_deformation = B.max_deformation
            
def getS(b):
    return ctx.states[b.bone.jiggle.index]

def applyConstraint(b,C,dC0, dC1, K):
    #return
    K1 =  1- pow(1-max(min(1,K),0),1.0/iters)   
    Jb = b.bone.jiggle
    w1  = 1.0/Jb.mass
    w0 = 0        
    Jbp = None
    if(b.parent != None and b.parent.bone.jiggle.enabled):
        w0 = 1.0/b.parent.bone.jiggle.mass
        Jbp = b.parent.bone.jiggle 
    div =     (dC0.dot(dC0)*w0 +dC1.dot(dC1)*w1 )  
    if(div > 0.000001):
        
        s = C/div
        #~ l0 = (-dC1*s*w1*K1).length
        #~ if(l0 > 10):
            #~ raise Exception('value error 1 '+ str(l0) + ' ' + str(div) + ñ)
        Jb.P+= -dC1*s*w1*K1
        if(Jbp!=None):
            Jbp.P+= -dC0*s*w0*K1                
            #~ l1 = ( -dC0*s*w0*K1).length
            #~ if(l1 > 10):
                #~ raise Exception('value error 2 '+ str(l1) + ' ' + str(div))
def nomu(rm, pm,ipm, sc):
    gmu =pm*rm
    for i in range(3):
        v = getAxis(gmu,i)
        v.length = sc[i]
        setAxis(gmu,i,v)
    return ipm*gmu    
def updateMat(b):

    Jb = b.bone.jiggle  
    Jbp = b.parent.bone.jiggle
    Sb = getS(b)
    Sbp = getS(b.parent)
        
    l = b.bone.length
    
    par = b.parent
    aM = nomu(Sb.rmat,Sbp.wmat,Sbp.iwmat,Sb.scale) #Sb.rmat #Sbp.wmat* Sb.rmat#ow*Sbp.wmat
    
    lP = Sbp.iwmat*Jb.P
    O = Vector((aM[0][3] ,aM[1][3] ,aM[2][3] ))
    V = lP- O #Vector((Jb.P.x-O.x ,Jb.P.y-O.y ,Jb.P.z-O.z ))
    lv = V.length
    ll = (Jb.P - Sbp.wmat*O).length
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
    if(cur.length > 100.001):
        print(cur.length)
    if(la>0):
        axis/=la#.normalized()
        #print(nV, ncur,cos)
        aM = Matrix.Rotation(-math.acos(cos),4, axis)*aM
        
        
        aM[0][3] = O.x #Jb.X.x-aM[0][3]    
        aM[1][3] = O.y #Jb.X.y-aM[1][3]
        aM[2][3] = O.z #Jb.X.z-aM[2][3]
        #deform = (lv/l)
        cur = getAxis(aM,1)
        cur *= (ll/l)/cur.length
        if((lv/l)/cur.length >2.0):
            print((lv/l)/cur.length)
        
        aM[0][1] = cur.x #Jb.X.x-aM[0][3]    
        aM[1][1] = cur.y #Jb.X.y-aM[1][3]
        aM[2][1] = cur.z #Jb.X.z-aM[2][3]
        
        #if(Jb.Sv>0.0 and lv>0.0):            
            #rv = math.sqrt(l/lv)#/cur.length
            #ax = getAxis(aM,0)
            #az = getAxis(aM,2)
            #ax.length = rv*(Jb.Sv) + ax.length*(1-Jb.Sv)
            #az.length = rv*(Jb.Sv) + az.length*(1-Jb.Sv)
            #setAxis(aM,0,ax)
            #setAxis(aM,2,az)
            
    lmat = normalizeM(aM) #Sbp.wmat.inverted()*aM)
    Sb.lmat = lmat
    Sb.setW(Sbp.wmat*lmat)
    
        #iow*aM
        #scene = bpy.context.scene
       # b.scale.y = V.length/l    
    
  
def getAxis(M,i):
    return   Vector((M[0][i],M[1][i],M[2][i]))

def setAxis(M,i,v):
    M[0][i] = v[0]
    M[1][i] = v[1]
    M[2][i] = v[2]
def resetBone(ow, iow,b):
    par = b.parent
    Sb = getS(b)
    Sbp = getS(b.parent)
    M = Sbp.wmat* Sb.rmat #im  
    l = b.bone.length    
    tg = Vector((M[0][1]*l+M[0][3],M[1][1]*l+M[1][3],M[2][1]*l+M[2][3]))#(M[0][1]+M[0][3],M[1][1]+M[1][3],M[2][1]+M[2][3]))
    
    Jb = b.bone.jiggle
    Jb.X=tg
    Jb.P = tg
    Jb.V= Vector((0,0,0))
def updateBone(b):
    par = b.parent
    Jb = b.bone.jiggle
    Jbp = b.parent.bone.jiggle
    Sb = getS(b)
    Sbp = getS(b.parent)
    
    #im = par.bone.matrix_local.inverted()* b.bone.matrix_local
    M = Sbp.wmat* Sb.rmat #im
    
    
     
    aM = M.copy() #ow*b.matrix.copy()
    


    N = getAxis(aM, 1)
    l = b.bone.length*Sb.scale[1]
    
    N/= N.length
    Nl = N*l
    X0 = getAxis(aM, 3)
    X1 = Jb.P
    X01 = X1-X0
    Xl1 = X1-X0-Nl
    applyConstraint(b,X01.dot(N)-l,-N,N,Jb.Kv)
    if(par.bone.jiggle.enabled):
        updateMat(par)           
    updateMat(b)
    X1 = Jb.P
    X01 = X1-X0
    Xl1 = X1-X0-Nl
    applyConstraint(b,X01.dot(X01)-l*l,-X01*2,X01*2,Jb.Kl) #length conservation
    if(par.bone.jiggle.enabled):
        updateMat(par)           
    updateMat(b)
    X1 = Jb.P
    X01 = X1-X0
    Xl1 = X1-X0-Nl
    applyConstraint(b,(Xl1).dot(Xl1),-Xl1*2,Xl1*2,Jb.Ks) 
    if(par.bone.jiggle.enabled):
        updateMat(par)     
    if(b.bone.jiggle_control>0.0):
        tm = Sbp.wmat*Sb.lm
        N = getAxis(tm, 1)            
        N/= N.length
        tp = X0 + N*l
        Jb.P += (tp - Jb.P)*b.bone.jiggle_control
    updateMat(b)
def controlBone(b):
    par = b.parent

    
    
    
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
        self.iwmat = None
        self.rmat = None       
        self.lmat = None
        self.scale = mathutils.Vector((1,1,1))
        self.control_P = None
        self.lm = None
    def setW(self,w):
        self.wmat = w
        self.iwmat = w.inverted()
        
class Ctx:    
    def __init__(self):
        self.states = [] 
ctx = None    

def bake():
    global ctx
    scene = bpy.context.scene
    scene.frame_set(scene.frame_start)
    for o in scene.objects:
        if(o.type == 'ARMATURE' and o.data.jiggle.enabled):
            arm = o.data
            ow = o.matrix_world
            iow = ow.inverted() 
            
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
        
        
class CopyJigglePropsOperator(bpy.types.Operator):
    bl_idname = "jiggle.copy_bone"
    bl_label = "Copy Jiggle Bone Properties"
    def execute(self, context):
        bone = context.bone
        for b in context.object.data.bones:
            if(b.select and b!=bone):
                copyObject(b.jiggle, bone.jiggle)
        return {'FINISHED'}    
def normalizeM(m):
    loc, rot, sca = m.decompose()
    mat_loc = mathutils.Matrix.Translation(loc)

    # create an identitiy matrix
    mat_sca = mathutils.Matrix()#.Scale(sca)
    mat_sca[0][0] = sca[0]
    mat_sca[1][1] = sca[1]
    mat_sca[2][2] = sca[2]
    mat_sca[3][3] = 1
    
    # create a rotation matrix
    mat_rot = rot.to_matrix()
    mat_rot.resize_4x4()
    # combine transformations
    return  mat_loc * mat_rot * mat_sca
@persistent
def update(scene, tm = False):
    global iters 
    global dt
    global ctx
    dt = 1.0/(scene.render.fps*scene.jiggle.sub_steps)
    #print("beg")
    for si in range(scene.jiggle.sub_steps):
            
        for o in scene.objects:
            if(o.type == 'ARMATURE' and o.data.jiggle.enabled and (scene.jiggle.test_mode or tm)):
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
                    s.setW(ow*b.matrix.copy())
                    M = b.bone.matrix_local
                    if(b.parent!=None):
                        M = b.parent.bone.matrix_local.inverted()*M                    
                    s.rmat = M     
                    M2 = ow*b.bone.matrix_local                                     
                    loc, rot, sca = M2.decompose()
                    s.scale = sca
                    if(b.bone.jiggle.enabled and b.parent!=None):
                        lM = b.parent.matrix.inverted() * b.matrix
                        s.control_P = getAxis(lM,1)
                        s.lm = lM
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
                        updateBone(b)
     
                for b in ol:
                    Jb = b.bone.jiggle
                    Jb.V= (Jb.P-Jb.X)/dt
                    Jb.X = Jb.P      
                    updateMat(b)
                    Sb = getS(b)
                    Sbp = getS(b.parent)
                    b.matrix_basis = Sb.rmat.inverted()*Sb.lmat#Sbp.wmat.inverted()*Sb.wmat
                    #b.matrix = iow*Sb.wmat               
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
        return (context.bone is not None and context.object is not None and context.object.type == 'ARMATURE' and context.object.data.jiggle.enabled)

    def draw_header(self, context):
        layout = self.layout
        bon = context.bone
        layout.prop(bon.jiggle, "enabled", text="")

    def draw(self, context):
        layout = self.layout

        bon = context.bone
        col = layout.column()
        if(bon.jiggle.enabled):
            col.prop(bon.jiggle,"Kv")   
            col.prop(bon.jiggle,"Kl")   
            col.prop(bon.jiggle,"Ks")          
            col.prop(bon.jiggle,"Kd")
            col.prop(bon.jiggle,"mass")
            col.prop(bon.jiggle,"debug")
            col.prop(bon.jiggle,"max_deformation")
            col.separator()
            col.prop(bon,"jiggle_control")
        col.operator("jiggle.copy_bone")
        col.operator("jiggle.reset")
            
           # col.prop(bon.jiggle,"Sv")


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
        col.operator("jiggle.bake")




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




class ResetJigglePropsOperator(bpy.types.Operator):
    bl_idname = "jiggle.reset"
    bl_label = "Reset State"
    def execute(self, context):
        scene = context.scene
        for o in scene.objects:
            if(o.select and o.type == 'ARMATURE' and o.data.jiggle.enabled and (scene.jiggle.test_mode or tm)):
                arm = o.data
                ow = o.matrix_world
                iow = ow.inverted()
                i=0
                for b in o.pose.bones:
                    if(b.bone.select):                    
                        M = ow *b.matrix #ow*Sbp.wmat* Sb.rmat #im
                        l = b.bone.length
                        tg = Vector((M[0][1]*l+M[0][3],M[1][1]*l+M[1][3],M[2][1]*l+M[2][3]))#(M[0][1]+M[0][3],M[1][1]+M[1][3],M[2][1]+M[2][3]))
                        
                        Jb = b.bone.jiggle
                        Jb.X=tg
                        Jb.P = tg
                        Jb.V= Vector((0,0,0))
        return {'FINISHED'}
def register():
    bpy.utils.register_class(JiggleScene)
    bpy.types.Scene.jiggle = bpy.props.PointerProperty(type = JiggleScene)
    bpy.utils.register_class(JiggleArmature)
    bpy.types.Armature.jiggle = bpy.props.PointerProperty(type = JiggleArmature)
    bpy.utils.register_class(JiggleBone)
    bpy.types.Bone.jiggle = bpy.props.PointerProperty(type = JiggleBone)
    bpy.utils.register_class(JiggleBonePanel)
    bpy.utils.register_class(JiggleArmaturePanel)
    bpy.utils.register_class(JiggleScenePanel)
    bpy.utils.register_class(BakeOperator)
    bpy.utils.register_class(CopyJigglePropsOperator)
    bpy.utils.register_class(ResetJigglePropsOperator)
    bpy.types.Bone.jiggle_control = bpy.props.FloatProperty(name = "control",min=0.0, max=1.0, default = 0.0)

    
    bpy.app.handlers.frame_change_post.append(update) 
    #print("Jiggle Armature Registered")
def unregister():    
    bpy.utils.unregister_class(JiggleArmature)
    bpy.utils.unregister_class(JiggleBone)
    bpy.utils.unregister_class(JiggleBonePanel)
    bpy.utils.unregister_class(JiggleArmaturePanel)
    bpy.utils.unregister_class(JiggleScene)
    bpy.utils.unregister_class(BakeOperator)
    bpy.utils.unregister_class(CopyJigglePropsOperator)
    bpy.utils.unregister_class(ResetJigglePropsOperator)
    
   
    bpy.app.handlers.frame_change_post.remove(update) 
if __name__ == '__main__':
	register()


