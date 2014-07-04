1. There are 5 different scenes in "worlds" folder. Using main.cpp @line 18 to change different scene.
2. The current releasing height is 1.3. if user needs to releas jello in a higher height, it may need to reduce the time step in Euler integration.  
3. There is a video in the output folder.  
4. I tried to use the same ks kd in different integration methods. RK4 definitely could accept the ks kd which are currently used on Euler integration. 
   However I am not so satisfied with the performance, so I decided to give RK4 its specified parameters.
Extra:
1. I create the function to colliside with other shapes.
2. User could press "Ctrl" and "mouse middle button" ,and meanwhile drag the mouse to generate a corresponding forces to push the jello
3. Create Leapfrog integration.

1.	What is the effect of the Ks and Kd parameters on the jello? 
Ks is the spring coefficient, when the spring is elongated in an identical length, the higher of the Ks, the grater force the spring will generate. 
In the jello simulation, the Ks provides the rebounding force and is relative to the distance between two points. 
On the other hand, Kd is the damping coefficient and is relative with the particle velocity. 
The higher velocity, the greater damping force will be generated, the force direction is opposite to the particle velocity. 
Damping force provides the cushion effect to mitigate the violent velocity change.

2.	What are the benefits and the drawbacks of the collision system used here? What are some different ways in which it could be improved?
In our collision system, we push back the particles which penetrate the surface back to the surface in one step and then provide a spring force to let it rebound. 
The benefit is that this algorithm is pretty simple and fast to calculate. 
The drawback is that if the velocity of the particle is too fast, the particle will penetrate too deep and cause either visual or computing problems. 
One way to improve this problem is using adaptive time step. 
When program finds the velocity is too fast and cause the particle to penetrate too deep, the program could go back to last step and try a smaller time step to integrate. 
The new time step could depend on the original penetration depth. 

3.	From lecture, what is an example of a stiff constrained system?
The stiff constrained system's response relies on object stiffness.
An infinitely hard collision would cause to a impulsive force response and dramatical velocity change, just like two rocks collide. 


4.	From lecture, what is the difference between and explicit and implicit integration scheme?
Explicit integration uses the current state to calculate future state situation. 
On the contrary, implicit integration uses the future state to compute the current state situation. 
The benefit of using implicit integration is because it is more stable. 
Implicit integration could guarantee the result convergence when providing a positive time step.

5.	Does the jello behave realistically?
Not exactly, the jello lack the Poisson effect. 
That means the volume of this jello in the simulation is not incompressible. 
Realistic jello should be volume conservative.

