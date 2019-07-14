# JiggleArmature
 Jiggle bone animation tool for blender
 
this is a jiggle bone animation add-on for blender,
to enable jiggle physics first enable Jiggle scene in the 
scene settings and then enable Jiggle Bone on the bones

features:
 
 - easy animation of softbody effects using just bones
 - real time simulation  
 - jiggle bone control through control bones and animatable custom rest poses
 - unconditionally stable simulation based on Position Based Dynamics
 
changelog  

* **Added checkbox to disable the simulation per armature**
* **Enabled animation of all JiggleBone properties** even the rest poses 
* **Improved control bone constraint** now supports objects and bones in other armatures
* **Simulation frame rate independent from scene frame rate**  for consistent simulation , the fps of the scene no longer affects the simulation behavior,  **time sub-stepping**  was removed since it’s no longer needed, you can change the simulation frame rate for an armature in the JiggleArmature settings in the Armature tab
* **Custom rest pose**  you can set custom rest poses for the bones
* **Bake simulation button is back**  for selected bones and for all bones (in the scene tab)
* **Improved solver**  it’s more stable than previous versions even on large chains of bones

[![Usage demo video](https://img.youtube.com/vi/Irv5P7PEc1k/0.jpg)](https://www.youtube.com/watch?v=Irv5P7PEc1k "Usage demonstration")
 
