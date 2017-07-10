# JiggleArmature 2
 Jiggle bone animation tool for blender
 
this is a jiggle bone animation add-on for blender,
to enable jiggle physics first enable Jiggle scene in the 
scene settings and then enable Jiggle Bone on the bones

features:
 
 - easy animation of softbody effects with only bones
 - real time simulation of the wiggle physics
 - multiple edition of parameters for selected jiggle bones
 - jiggle bone control through control bones
 - unconditionally stable simulation
limitations:
 
 - the phisical behaviour depends on fps of the scene and the selected number of iterations
 
changelog from version 1
 
 - volume, length and shape parameters reduced to one global **stiffness** parameter
 - jiggle bones now can't stretch (because it's not possible to set non ortogonal matrices to the blender bones)
 - removed sub steps parameter provisionally (because it introduced annoying discontinuities) 
 - removed bake button, if you need to record the animation just enable automatic keyframe insertion (the red button) and play the animation
 - **rotation based physics**, this version allows rotational wiggle effects about the bone axes  
 - **control bones**, you can control now the rotation of the jiggle bones using control bones 

more info at (soon): 
 http://cheece.github.io/JiggleArmature/
