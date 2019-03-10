# JiggleArmature
 Jiggle bone animation tool for blender
 
this is a jiggle bone animation add-on for blender,
to enable jiggle physics first enable Jiggle scene in the 
scene settings and then enable Jiggle Bone on the bones

features:
 
 - easy animation of softbody effects using just bones
 - real time simulation  
 - jiggle bone control through control bones
 - unconditionally stable simulation 
 
changelog  
 
* **Simulation frame rate independent from scene frame rate**  for consistent simulation , the fps of the scene no longer affects the simulation behavior,  **time sub-stepping**  was removed since it’s no longer needed, you can change the simulation frame rate for an armature in the JiggleArmature settings in the Armature tab
* **Custom rest pose**  you can set custom rest poses for the bones
* **Bake simulation button is back**  for selected bones and for all bones (in the scene tab)
* **Improved solver**  it’s more stable than previous versions even on large chains of boneses
