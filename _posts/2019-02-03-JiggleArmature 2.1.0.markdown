---
layout: post
title:  "1.2.0 released!"
date:   2019-02-03 00:50:42 -0430
categories: version release
---
The new version [0.2.1](https://github.com/cheece/JiggleArmature/) is ready!.

the features implemented in this version are:


* **Simulation frame rate independent from scene frame rate**  for consistent simulation , the fps of the scene no longer affects the simulation behavior,  **time sub-stepping**  was removed since it’s no longer needed, you can change the simulation frame rate for an armature in the JiggleArmature settings in the Armature tab
* **Custom rest pose**  you can set custom rest poses for the bones
* **Bake simulation button is back**  for selected bones and for all bones (in the scene tab)
* **Improved solver**  it’s more stable than previous versions even on large chains of boneses
