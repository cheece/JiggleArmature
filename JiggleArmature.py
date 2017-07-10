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

#to enable jiggle physics first enable "jiggle scene" in the scene properties and then enable jiggle bone on the bones
 
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
    test_mode = bpy.props.BoolProperty(default=False)
    sub_steps = bpy.props.IntProperty(min=1, default = 2)
    iterations = bpy.props.IntProperty(min=1, default = 4)
    last_frame = bpy.props.IntProperty()
    length_fix_iters = bpy.props.IntProperty(min=0, default = 2)



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
        #col.prop(context.scene.jiggle,"length_fix_iters")
        #col.operator("jiggle.bake")



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
    W = bpy.props.FloatVectorProperty(size=3,subtype='XYZ')
    P = bpy.props.FloatVectorProperty(size=3,subtype='XYZ')
    V = bpy.props.FloatVectorProperty(size=3,subtype='XYZ')
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
                        Jb.P = mpos(M)
                        Jb.W = Vector((0,0,0))
                        #Jb.M = Matrix(M)
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
            col.operator("jiggle.reset")
            if(bon.parent==None):                
                col.label(text= "warning: jibblebones without parent will fall",icon='COLOR_RED')   

def centerM(wb,l):
    ax = maxis(wb,1).normalized()
    wb[0][3] += ax[0]*l*0.5
    wb[1][3] += ax[1]*l*0.5
    wb[2][3] += ax[2]*l*0.5

class JB:
    def __init__(self, b,M,p):
        self.M  = M.copy()
        self.length = b.bone.length*maxis(M,0).length
       # self.P += maxis(M,1)*0.5*b.bone.length
        self.b = b
        self.parent = p
        self.rest = None
        self.w = 0
        self.Cx = None
    def sample(self, t):
        return self.P + maxis(self.R,1)*(t*self.length)
    def applyImpulse(self, p, I,t=0.5):
        if(self.w>0.0):
            I = I*self.w
            #ip = self.M.inverted()*p
            #p = self.sample(t)
            C = self.sample(t)
            r = p-C
            tg = p+I
            r2 = tg-C
            ax = r.cross(r2)
            lax = ax.length
            if(lax > 0.000001):
                r.normalize()
                r2.normalize()
                cos = r.dot(r2)
                ax.normalize()
                
                if(cos < -0.9999999):
                    cos = -0.9999999
                if(cos >  0.9999999):
                    cos =  0.9999999
                ag = math.acos(cos)
                mr = Matrix.Rotation(ag , 3, ax)#.to_quaternion()
                self.R = (mr*self.R).normalized()
                #self.Q = (mr*self.Q).normalized()
                self.P+= C- self.sample(t)
                p = mr*(p-C) + C#self.sample(t)
                 #p = self.M*ip
            #print((tg-p).length)
            self.P+= tg-p
    @property
    def P(self):
        return Vector((self.M[0][3],self.M[1][3],self.M[2][3]))
    @P.setter
    def P(self, x):
        self.M[0][3]= x[0]
        self.M[1][3]= x[1]
        self.M[2][3]= x[2]
    @property
    def Q(self):
        return self.M.to_quaternion()    
    @Q.setter
    def Q(self, x):
        m = x.to_matrix()
        for i in range(3):
            for j in range(3):
                self.M[i][j] = m[i][j]         
    @property
    def R(self):
        return self.M.to_3x3()    
    @R.setter
    def R(self, x):
        for i in range(3):
            for j in range(3):
                self.M[i][j] = x[i][j]     
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
        
def quatSpring(Jb,r=None,k=None):
    Q0 = Jb.parent.Q
    Q1 = Jb.Q
    w0 = Jb.parent.w
    w1 = Jb.w
    if(r==None):
        r = Jb.rest.to_quaternion()
        k = Jb.k
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
    
    tmp0=Q0w*rx+Q0x*rw+Q0y*rz+Q0z*ry*(-1)
    tmp1=Q0w*rw+Q0x*rx*(-1)+Q0y*ry*(-1)+Q0z*rz*(-1)
    tmp2=Q0w*rz+Q0x*ry+Q0y*rx*(-1)+Q0z*rw
    tmp3=Q0w*ry+Q0x*rz*(-1)+Q0y*rw+Q0z*rx
    tmp4=Q1w*tmp0+Q1x*(-1)*tmp1+Q1y*(-1)*tmp2+tmp3*Q1z
    tmp5=Q1w*tmp3+tmp2*Q1x+Q1y*(-1)*tmp1+Q1z*(-1)*tmp0
    tmp6=Q1w*tmp2+Q1x*(-1)*tmp3+tmp0*Q1y+Q1z*(-1)*tmp1
    tmp7=Q1x*rz
    tmp8=Q1y*rx
    tmp9=rw*Q1w
    tmp10=Q1z*ry
    c = pow(tmp4,2)+pow(tmp5,2)+pow(tmp6,2) 
    dQ0x = tmp6*2*(tmp7+rw*Q1y+Q1z*rx+ry*Q1w)+tmp5*2*(ry*Q1x+tmp8+Q1z*rw*(-1)+rz*(-1)*Q1w)+tmp4*2*(Q1x*rx+Q1y*ry*(-1)+rz*(-1)*Q1z+tmp9) 
    dQ0y = tmp6*2*(Q1x*rw*(-1)+rz*Q1y+tmp10+rx*(-1)*Q1w)+tmp5*2*(rx*(-1)*Q1x+Q1y*ry+Q1z*rz*(-1)+tmp9)+tmp4*2*(Q1x*ry+tmp8+rw*Q1z+rz*Q1w) 
    dQ0z = tmp6*2*(Q1x*rx*(-1)+ry*(-1)*Q1y+Q1z*rz+tmp9)+tmp5*2*(rw*Q1x+Q1y*rz+tmp10+rx*Q1w)+tmp4*2*(tmp7+Q1y*rw*(-1)+rx*Q1z+ry*(-1)*Q1w) 
    dQ0w = tmp6*2*(tmp8+(-1)*rw*Q1z+(-1)*ry*Q1x+Q1w*rz)+tmp5*2*((-1)*rx*Q1z+(-1)*rw*Q1y+tmp7+Q1w*ry)+tmp4*2*((-1)*rz*Q1y+tmp10+(-1)*rw*Q1x+Q1w*rx) 
    dQ1x = tmp6*(-2)*tmp3+2*tmp5*tmp2+tmp4*(-2)*tmp1 
    dQ1y = 2*tmp6*tmp0+tmp5*(-2)*tmp1+tmp4*(-2)*tmp2 
    dQ1z = tmp6*(-2)*tmp1+tmp5*(-2)*tmp0+2*tmp4*tmp3 
    dQ1w = 2*tmp5*tmp3+2*tmp6*tmp2+2*tmp4*tmp0 

    
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
    sub_steps = 1#*scene.jiggle.sub_steps)
    dt = 1.0/(scene.render.fps*sub_steps) 
  
    for o in scene.objects:
        if(o.type == 'ARMATURE'):
            
            arm = o.data
            ow = o.matrix_world.copy()
            scale = maxis(ow,0).length
            #normR(ow)
            iow = ow.inverted()
            iow3 = ow.to_3x3().inverted()
            
            i=0
            for j in range( sub_steps):  
                bl = []
                wt = []
                        
                for b in o.pose.bones:
                    if(b.parent==None):
                        propB(ow,b,bl,None)
                hooks= []
                
                bl2 = []
                for wb in bl:
                    b = wb.b
                    crest = b
                    # o------ -> ---o---
                    wb.restW = b.bone.matrix_local.copy() * scale
                    wb.Q = wb.Q.normalized()
                    
                    if(b.bone.jiggle.enabled):
                    #(wb.parent.restW.inverted() * wb.restW) #
                        Jb = b.bone.jiggle
                        wb.rest =  b.bone.matrix_local #
                        if(b.parent!=None):
                            wb.rest = b.bone.parent.matrix_local.inverted()*wb.rest
                        wb.Kc = 0
                        if(Jb.control_bone!=""):
                            if(Jb.control_bone in o.pose.bones):
                                cb = o.pose.bones[Jb.control_bone]
                                clm = cb.matrix
                                if(cb.parent!=None):
                                    clm = cb.parent.matrix.inverted()*clm
                                wb.cQ = clm.to_quaternion().normalized()
                                wb.Kc = 1- pow(1-Jb.control, 1/scene.jiggle.iterations)
                            
                        
                        wb.rest_base = wb.rest.copy()
                        saxis(wb.rest,3, mpos(wb.rest)*scale)
                        wb.length = b.bone.length*scale
                        wb.irest = wb.rest.inverted()
                        wb.w = 1.0/Jb.mass
                        wb.k = 1- pow(1-Jb.Ks, 1/scene.jiggle.iterations)
                        Jb.V*= 1.0-Jb.Kld
                        Jb.V+= scene.gravity*dt
                        Jb.W*= 1.0-Jb.Kd
                        R = Jb.R.to_matrix()
                        wb.R = R.normalized()
                        wb.P = Jb.P.copy()
                        wb.Cx = wb.sample(0.5)
                        qv = Quaternion()
                        qv.x = Jb.W[0]
                        qv.y = Jb.W[1]
                        qv.z = Jb.W[2]
                        qv.w = 0         
                        cv = wb.Cx + Jb.V*dt                
                        wb.Q = qadd(wb.Q, qv*wb.Q*dt*0.5).normalized()    #newton's first law   
                        wb.P += cv - wb.sample(0.5)#the same
                        
                        bl2.append(wb)
                        
                                   
                for i in range(scene.jiggle.iterations):                        
                    for wb in bl2:
                        b = wb.b
                        if(b.parent==None):
                            continue
                        Jb = b.bone.jiggle                        
                        Pc =  wb.P
                        target_m = wb.parent.M*wb.rest
                        Pt = mpos(target_m)
                        if(Jb.debug in scene.objects):
                            scene.objects[Jb.debug].location = Pc
                        W = wb.w + wb.parent.w
                        I = (Pc-Pt)/W                       
                        
                        wb.applyImpulse(Pc,-I)
                        wb.parent.applyImpulse(Pt,I)
                        
                #for i in range(scene.jiggle.iterations):                        
                    for wb in bl2:
                        b = wb.b        
                        if(b.parent==None):
                            continue                
                        Jb = b.bone.jiggle                        
                        quatSpring(wb)                       
                        if(wb.Kc>0.0):
                            quatSpring(wb,wb.cQ, wb.Kc)
                for wb in bl2:
                    b = wb.b
                    if(b.parent==None):
                        continue
                    Jb = b.bone.jiggle                    
                    wb.P = wb.parent.M*mpos(wb.rest)
                    wb.R = wb.Q.normalized().to_matrix() 
               
                                     
                for wb in bl2:
                    b = wb.b
                    Jb = b.bone.jiggle
                    R = Jb.R.to_matrix()              
                    Jb.V = (wb.sample(0.5) - wb.Cx)/dt
                    Jb.P = wb.P.copy()
                    wb.Q = wb.Q.normalized()
                    qv = wb.Q*Jb.R.conjugated() #qadd(wb.Q,-Jb.R)*Jb.R.conjugated()#
                    Jb.W = Vector((qv.x,qv.y,qv.z))*(2/dt)
                    Jb.R = wb.Q
                                              
                for wb in bl:
                    wb.R*=scale                
                for wb in bl2:
                    b = wb.b
                    pM = ow
                    if(b.parent!=None):
                        pM = wb.parent.M
                    b.matrix_basis = (pM*wb.rest_base).inverted()*wb.M
                   # b.matrix = iow*wb.M
                    
    scene.jiggle.last_frame+= 1
    print("updated")
    
@persistent
def update_post(scene, tm = False):
    global iters 
    global dt
    global ctx
    global cc
    
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
                        Jb.V = Vector((0,0,0))
                        Jb.P = mpos(M)          
                        Jb.W = Vector((0,0,0))
                    
        if(scene.frame_current <= bpy.context.scene.frame_start):
            return
    nframes = scene.frame_current - scene.jiggle.last_frame
    for i in range(nframes):     
        #
        step(scene)        

def register():
    bpy.app.handlers.frame_change_pre.append(update) 
    
                
    bpy.utils.register_class(JiggleScene)
    bpy.types.Scene.jiggle = bpy.props.PointerProperty(type = JiggleScene)
    bpy.utils.register_class(JiggleScenePanel)
    bpy.utils.register_class(JiggleBone)

    bpy.types.Bone.jiggle = bpy.props.PointerProperty(type = JiggleBone)
    bpy.utils.register_class(ResetJigglePropsOperator)

    bpy.utils.register_class(JiggleBonePanel)
def unregister():    

    bpy.utils.unregister_class(JiggleScene)
    bpy.utils.unregister_class(JiggleScenePanel)
    bpy.utils.unregister_class(JiggleBone)
    bpy.utils.unregister_class(ResetJigglePropsOperator)
    bpy.utils.unregister_class(JiggleBonePanel)
    bpy.app.handlers.frame_change_pre.remove(update) 
if __name__ == '__main__':
	register()
