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


#based on the Position Based Dynamics paper by Müller et al. http://matthias-mueller-fischer.ch/publications/posBasedDyn.pdf

bl_info = {
    "name": "Jiggle Armature",
    "author": "Simón Flores",
    "version": (2, 2, 1),
    "blender": (2, 79, 0),
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



    
class JiggleArmature(bpy.types.PropertyGroup): 
    enabled= bpy.props.BoolProperty(default=True)
    fps= bpy.props.FloatProperty(name = "simulation fps",default = 24)
    time_acc= bpy.props.FloatProperty(default= 0.0)
    
    
class JiggleScene(bpy.types.PropertyGroup):
    test_mode= bpy.props.BoolProperty(default=False)
    sub_steps= bpy.props.IntProperty(min=1, default = 2)
    iterations= bpy.props.IntProperty(min=1, default = 4) 
    last_frame= bpy.props.IntProperty()
    



class JARM_PT_armature(bpy.types.Panel):
    bl_idname = "ARMATURE_PT_jiggle"
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
        col.prop(context.object.data.jiggle,"enabled")
        col.prop(context.object.data.jiggle,"fps")
        

class JARM_PT_scene(bpy.types.Panel):
    bl_idname = "SCENE_PT_jiggle"
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
        
        col.prop(context.scene.jiggle,"iterations")
        
        
        col.operator("jiggle.bake", text="Bake Selected").a = False
        col.operator("jiggle.bake", text="Bake All").a = True



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
                setattr(b2, prop, getattr(b,prop)) 
        inop = False
    return f 


def setq(om, m):
    for i in range(4):
        om[i]= m[i]
        
class JARM_OT_reset(bpy.types.Operator):
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
                        M = ow*b.matrix
                        
                        
                        Jb = b.bone 
                        setq(Jb.jiggle_R, M.to_quaternion().normalized())
                        Jb.jiggle_V = Vector((0,0,0))
                        Jb.jiggle_P =  mpos(M)+maxis(M,1)*b.bone.length*0.5
                        Jb.jiggle_W = Vector((0,0,0))
                        
        return {'FINISHED'}
        
class JARM_OT_set_rest(bpy.types.Operator):
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
                        Jb = b.bone 
                        setq(Jb.jiggle_rest, M.to_quaternion().normalized()) 
                        Jb.jiggle_use_custom_rest = True
        return {'FINISHED'}

class JARM_PT_bone(bpy.types.Panel):
    bl_idname = "BONE_PT_jiggle_bone"
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
        layout.prop(bon , "jiggle_enabled", text="")

    def draw(self, context):
        layout = self.layout
        armature = context.object.data
        
        bon = context.bone
        col = layout.column()
        layout.enabled = context.scene.jiggle.test_mode and armature.jiggle.enabled
        
        if(not context.scene.jiggle.test_mode):
            col.label(text= "JiggleArmature is disabled for the scene, see scene properties")    
            
        if(not armature.jiggle.enabled):
            col.label(text= "JiggleArmature is disabled for the armature, see the armature properties")   
            
        if(bon.jiggle_enabled):  

            col.prop(context.bone,"jiggle_Ks")    
            col.prop(bon,"jiggle_Kd")  
            col.prop(bon,"jiggle_Kld")             
            col.prop(bon,"jiggle_mass")  
            
            col.prop_search(bon,"jiggle_control_object",bpy.data,"objects")
            if(bon.jiggle_control_object in bpy.data.objects): 
                o = bpy.data.objects[bon.jiggle_control_object]
                if(o.type == 'ARMATURE'):
                     col.prop_search(bon,"jiggle_control_bone",o.data,"bones")
                col.prop(bon,"jiggle_control") 
            
            col.operator("jiggle.reset")
            col.prop(bon,"jiggle_use_custom_rest") 
            if(bon.jiggle_use_custom_rest):
                col.prop(bon,"jiggle_rest")
            col.operator("jiggle.set_rest")
            if(bon.parent==None):                
                col.label(text= "warning: jibblebones without parent will fall",icon='COLOR_RED')   

def centerM(wb,l):
    ax = maxis(wb,1).normalized()
    wb[0][3] += ax[0]*l*0.5
    wb[1][3] += ax[1]*l*0.5
    wb[2][3] += ax[2]*l*0.5


#adapted from https://github.com/InteractiveComputerGraphics/PositionBasedDynamics/blob/master/PositionBasedDynamics/PositionBasedRigidBodyDynamics.cpp
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
        self.Kc = 0
        self.cQ = None
        self.X = None
        self.P = None
        self.R = None
        self.Q = None
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
         

sqrt = math.sqrt

#NOTE: the following gradient computation implementation was automatically generated, if possible, it should be change for a clearer implementation 
def quatSpringGradient2(Q0,Q1,r):
    """
    returns the gradient of C = |Q0*r - Q1|^2 wrt Q0 and Q1  
    """
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
    
    
    tmp0 = math.sqrt(((((((((-(Q0x*Q1w)-(Q0y*Q1z))+(Q0w*Q1x))+(Q0z*Q1y))-rx)*((((-(Q0x*Q1w)-(Q0y*Q1z))+(Q0w*Q1x))+(Q0z*Q1y))-rx))+(((((-(Q0x*Q1y)-(Q0z*Q1w))+(Q0w*Q1z))+(Q0y*Q1x))-rz)*((((-(Q0x*Q1y)-(Q0z*Q1w))+(Q0w*Q1z))+(Q0y*Q1x))-rz)))+(((((-(Q0y*Q1w)-(Q0z*Q1x))+(Q0w*Q1y))+(Q0x*Q1z))-ry)*((((-(Q0y*Q1w)-(Q0z*Q1x))+(Q0w*Q1y))+(Q0x*Q1z))-ry)))+((((((Q0w*Q1w)+(Q0x*Q1x))+(Q0y*Q1y))+(Q0z*Q1z))-rw)*(((((Q0w*Q1w)+(Q0x*Q1x))+(Q0y*Q1y))+(Q0z*Q1z))-rw))))
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
            
            Q0.x+=dQ0x*s*w0*k
            Q0.y+=dQ0y*s*w0*k
            Q0.z+=dQ0z*s*w0*k 
            Q0.w+=dQ0w*s*w0*k    
            Jb.parent.Q = Q0.normalized()    
            
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
    
    dt = 1.0/(scene.render.fps) 
  
    for o in scene.objects:
        if(o.type == 'ARMATURE' and o.data.jiggle.enabled): 
            
            arm = o.data
            ow = o.matrix_world.copy()
            scale =maxis(ow,0).length 
            
            iow = ow.inverted()
            iow3 = ow.to_3x3().inverted()
            
            i=0
            arm.jiggle.time_acc+= dt* arm.jiggle.fps 
            while arm.jiggle.time_acc > 1:
                arm.jiggle.time_acc-=1 
                
                bl = []
                wt = []
                        
                for b in o.pose.bones:
                    if(b.parent==None):
                        propB(ow,b,bl,None)
                hooks= []
                
                bl2 = []
                for wb in bl:
                    b = wb.b
                    wb.rest_w = b.bone.matrix_local.copy()
                    saxis(wb.rest_w,3, maxis(wb.rest_w,3)*scale)  
                    saxis(wb.rest_w,3, maxis(wb.rest_w,3)+maxis(wb.rest_w,1)*b.bone.length*0.5*scale)
                    
                for wb in bl:
                    b = wb.b
                    crest = b 
                    
                    wb.restW = b.bone.matrix_local.copy() * scale 
                    saxis(wb.restW,3, maxis(wb.restW,3)*scale)    
                    
                    M = wb.M 
                    if(b.bone.jiggle_enabled):
                        Jb = b.bone 
                        wb.X = wb.P = Jb.jiggle_P
                        wb.R = wb.Q = Jb.jiggle_R
                        wb.rest =  wb.rest_w  
                        if(b.parent!=None):
                            wb.rest = wb.parent.rest_w.inverted()*wb.rest_w 
                
                        wb.rest_base = b.bone.matrix_local
                        if(b.parent!=None):
                            wb.rest_base = b.parent.bone.matrix_local.inverted()*wb.rest_base
                            
                        wb.rest_p = wb.parent.rest_w.inverted()* (maxis(wb.rest_w,3)- maxis(wb.rest_w,1)*b.bone.length*0.5*scale)# mpos(wb.rest)
                        
                        wb.l = b.bone.length*scale
                        wb.w = 1.0/Jb.jiggle_mass
                        wb.k = 1- pow(1-Jb.jiggle_Ks, 1/scene.jiggle.iterations)
                        Jb.jiggle_V*= 1.0-Jb.jiggle_Kld
                        Jb.jiggle_V+= scene.gravity*dt
                        Jb.jiggle_W*= 1.0-Jb.jiggle_Kd
                        qv = Quaternion()
                        qv.x =Jb.jiggle_W[0]
                        qv.y = Jb.jiggle_W[1]
                        qv.z = Jb.jiggle_W[2]
                        qv.w = 0         
                                
                        wb.Q = qadd(wb.Q, qv*wb.Q*dt*0.5).normalized() 
                        
                        wb.P = wb.X + Jb.jiggle_V*dt
                        wb.computeI()
                        
                        #control object/bone constraint
                        if(Jb.jiggle_control_object in bpy.data.objects): 
                            
                            target_object =  bpy.data.objects[Jb.jiggle_control_object]
                            
                            target_matrix =  target_object.matrix_local
                            
                            if(target_object.type =='ARMATURE' and Jb.jiggle_control_bone in target_object.pose.bones):
                                cb = target_object.pose.bones[Jb.jiggle_control_bone]
                                target_matrix = cb.matrix
                                if(cb.parent!=None):
                                    target_matrix = cb.parent.matrix.inverted()*target_matrix
                                    
                            wb.cQ = target_matrix.to_quaternion().normalized()
                            wb.Kc = 1- pow(1-Jb.jiggle_control, 1.0/scene.jiggle.iterations)
                        
                        
                            
                        bl2.append(wb) 
                    else:
                        wb.w = 0
                        wb.X = wb.P = mpos(M)+maxis(M,1)*b.bone.length*0.5
                        wb.R = wb.Q = M.to_quaternion().normalized()
                               
                        M = ow*b.matrix  
                        
                        Jb = b.bone 
                        setq(Jb.jiggle_R, M.to_quaternion().normalized())
                        Jb.jiggle_V = Vector((0,0,0))
                        Jb.jiggle_P =  mpos(M)+maxis(M,1)*b.bone.length*0.5
                        Jb.jiggle_W = Vector((0,0,0))
                         
                                   
                for i in range(scene.jiggle.iterations):      
                    #parent constraint       
                    for wb in bl2:
                        b = wb.b        
                        if(b.parent==None):
                            continue                
                        Jb = b.bone                
                        locSpring(wb)                 
                    #spring constraint                    
                    for wb in bl2:
                        b = wb.b        
                        if(b.parent==None):
                            continue                
                        Jb = b.bone                   
                        quatSpring(wb,  Jb.jiggle_rest if Jb.jiggle_use_custom_rest else wb.rest.to_quaternion().normalized())   
                        if(wb.cQ!=None):
                            quatSpring(wb, wb.cQ, wb.Kc)    


                for wb in bl2:
                    b = wb.b
                    Jb = b.bone 
                    
                    wb.Q = wb.Q.normalized()
                    m = wb.Q.to_matrix()
                    for i in range(3):
                        for j in range(3):
                            wb.M[i][j] = m[i][j]*scale      
                    wb.M[3][3]=1
 
                    Jb.jiggle_V = (wb.P - wb.X)/dt
                    Jb.jiggle_P = wb.P.copy()
                    qv = wb.Q*Jb.jiggle_R.conjugated() 
                    Jb.jiggle_W = Vector((qv.x,qv.y,qv.z))*(2/dt)
                    Jb.jiggle_R = wb.Q
 
                    cp = Jb.jiggle_P - maxis(wb.M,1)*b.bone.length*0.5
             
                    wb.M[0][3]= cp[0]
                    wb.M[1][3]= cp[1]
                    wb.M[2][3]= cp[2]
 
                for wb in bl2:
                    b = wb.b
                    pM = ow
                    if(b.parent!=None):
                        pM = wb.parent.M
                    mb =  (pM*wb.rest_base).inverted()*wb.M
 
                    b.matrix_basis = mb
 
                    
    scene.jiggle.last_frame+= 1
 
    
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

def bake(bake_all):
    print("baking " + ("all" if(bake_all) else "selected") + "...")
    global ctx
   # global backing 
    
    scene = bpy.context.scene
    scene.frame_set(scene.frame_start)
    
    for o in scene.objects:
        if(o.type == 'ARMATURE' and (o.select or bake_all)):
            
            
            arm = o.data
            ow = o.matrix_world
            scale = maxis(ow,0).length
            iow = ow.inverted()
            i=0
            for b in o.pose.bones:
                b.bone.select = (b.bone.select or bake_all) and b.bone.jiggle_enabled
                if(b.bone.select or bake_all and b.bone.jiggle_enabled):                    
                    M = ow*b.matrix  
                    
                    Jb = b.bone 
                    setq(Jb.jiggle_R, M.to_quaternion().normalized())
                    Jb.jiggle_V = Vector((0,0,0))
                    Jb.jiggle_P = mpos(M)
                    Jb.jiggle_W = Vector((0,0,0))
    
    ltm = scene.jiggle.test_mode                
    scene.jiggle.test_mode = False 
    for i in range(scene.frame_start, scene.frame_end):
        scene.frame_set(i)
        update(scene,tm=True)
        print("frame: ",i)
        for o in scene.objects:
            if( (o.select or bake_all) and o.type == 'ARMATURE' ):
                bpy.context.scene.objects.active = o 
                m = o.mode == 'POSE'
                
                if(not m):
                    bpy.ops.object.posemode_toggle()
                
                bpy.ops.anim.keyframe_insert_menu(type='LocRotScale')
        
                if(not m):
                    bpy.ops.object.posemode_toggle() 
    scene.jiggle.test_mode = ltm
                   
class JARM_OT_bake(bpy.types.Operator):
    a = bpy.props.BoolProperty()
    bl_idname = "jiggle.bake"
    bl_label = "Bake Animation"
    def execute(self, context):
        bake(self.a)           
        return {'FINISHED'}    

classes = ( 
    JARM_PT_armature,
    JiggleScene,
    JARM_PT_scene,
    JiggleArmature,
    JARM_OT_bake,
    JARM_OT_set_rest,
    JARM_OT_reset,
    JARM_PT_bone
)

def register():
    
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls) 
        
    bpy.app.handlers.frame_change_post.append(update) 
     
    bpy.types.Scene.jiggle = bpy.props.PointerProperty(type = JiggleScene) 

    bpy.types.Armature.jiggle = bpy.props.PointerProperty(type = JiggleArmature,options={'ANIMATABLE'}) 
    
    bpy.types.Bone.jiggle_enabled = bpy.props.BoolProperty(default=False, update = funp("jiggle_enabled"))
    bpy.types.Bone.jiggle_Kld=bpy.props.FloatProperty(name = "linear damping",min=0.0, max=1.0,default = 0.01, update = funp("jiggle_Kld"))
    bpy.types.Bone.jiggle_Kd =bpy.props.FloatProperty(name = "angular damping",min=0.0, max=1.0,default = 0.01, update = funp("jiggle_Kd"))
    bpy.types.Bone.jiggle_Ks =bpy.props.FloatProperty(name = "stiffness",min=0.0 , max = 1.0, default = 0.8, update = funp("jiggle_Ks"))
    bpy.types.Bone.jiggle_mass =bpy.props.FloatProperty(name = "mass",min=0.0001, default = 1.0, update = funp("jiggle_mass"))   
    bpy.types.Bone.jiggle_R = bpy.props.FloatVectorProperty(name="rotation", size=4,subtype='QUATERNION')
    bpy.types.Bone.jiggle_W = bpy.props.FloatVectorProperty(size=3,subtype='XYZ') #angular velocity
    bpy.types.Bone.jiggle_P = bpy.props.FloatVectorProperty(size=3,subtype='XYZ')
    bpy.types.Bone.jiggle_V = bpy.props.FloatVectorProperty(size=3,subtype='XYZ')	#linear velocity, ok? 
    
    bpy.types.Bone.jiggle_use_custom_rest =bpy.props.BoolProperty(default=False, name="Use Custom Rest Pose",  update = funp("jiggle_use_custom_rest"))
    
    bpy.types.Bone.jiggle_rest =bpy.props.FloatVectorProperty(name="rotation", size=4,subtype='QUATERNION')
    
    bpy.types.Bone.jiggle_control =bpy.props.FloatProperty(name = "control",min=0.0, max=1.0,default = 1, update = funp("jiggle_control"))
    bpy.types.Bone.jiggle_control_object =bpy.props.StringProperty(name = "control object")
    bpy.types.Bone.jiggle_control_bone =bpy.props.StringProperty(name = "control bone")
    
      
def unregister():    

    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
         
    bpy.app.handlers.frame_change_post.remove(update) 
if __name__ == '__main__':
	register()
