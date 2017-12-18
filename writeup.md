# Project Write-Up

## Co-ordinate Transformation  
The controller and the motion model are defined in the car's frame of reference. The transformation between the global co-ordinates to the car's frame of reference is   
```
x_c = x * cos(psi) + y * sin(psi)  
y_c = -x * sin(psi) + y * cos(psi)
psi_c = 0
```

## Motion Model
The motion model is describred by the following equations:

Given the vehicle state - {x, y, psi, v, cte, epsi, delta, a}  
where:  
*x*      - X position of the car in the car's frame of reference  
*y*      - Y position of the car in the car's frame of reference  
*psi*    - Oriention of the car with respect to its longitudinal axis  
*v*      - Velocity of the car along its longitudinal axis  
*cte*    - Cross-Track Error  
*epsi*   - Error in the orientation of the car  
*delta*  - Steering angle (control input)  
*a*      - Acceleration/Throttle (control input)  

*Lf* - a constant obtained by measuring the radius formed by running the vehicle in the simulator around in a circle with a constant steering angle and velocity on a flat terrain.  
*y = f(x)* - A polynomial describing the desired trajectory. This is calculated by fitting a polynomial to the desired waypoints in the car's frame of reference.    

The motion model is defining the transition of the car from t to t + 1 is  
```
dt = (t + 1) - t  
psides[t] = arctan(f'(x[t])) - The desired psi at time t
x[t+1] = x[t] + v * cos(psi)  
y[t+1] = y[t] + v * sin(psi)  
psi[t+1] = psi[t] +  v[t] / Lf * delta[t] * dt  
v[t+1] = v[t] + a[t] * dt  
cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt  
epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt  
```

## The MPC Controller  
The MPC controller is an optimization based controller. 
The optimization (minimization) problem is solved with the following parameters:

### Cost Function
The cost function is composed of the following components  
1. Cross track error (cte) - This measures the trajectory error in x and y
2. Car orientation error (epsi) - This measures the orientation error of the car
3. Velocity error - This measures the error between actual and desired velocities)
4. Control inputs - These are included so that the solution is achieved with the lowest possible control inputs
5. Change in control inputs - These are included so that the solution is achieved with the lowest possible rate of change of control inputs resulting in a smooth trajectory

Individual components are weighted differently to achieve the desired result of following the desired trajectory with smooth control inputs.  

### Constraints
The constraints of the optimization problem are defined by the motion model described above.

The motion model constraint is  
```State[t+1] - State[t] = 0```

### Control variables bounds
Upper and lower bounds are specified for the control variables involved.  
For steering control input (delta), this is limited to -25 deg to +25 deg  
For acceleration input, this is limited to -1 m/s2 to 1 m/s2  

### Initialization
The initial state of the car is provided by the simulator. The initial state of car is then modified to take into account the latency of the controls. The motion model is applied to predict the new initial state of the car with a timestep equal to the latency.  

### Model Hyper-Parameters Selection
dt - Size of each timestep in seconds  
N - Number of timesteps  

The selection of the above parameters determines the time horizon considered by the controller in predicting the car's trajectory. 

**N** - Increases the number of variables to optimize and hence increases the computation cost  
**dt** - A lower *dt* better approximates the continuous trajectory function and hence reduces discretization error. 

The selection of *N* and *dt* depends on the time horizon of the predicted trajectory that is desired. In the case of a car, anything more than a second is not useful since the conditions would have changed sufficiently during this duration. 

#### Tuning
The parameters *N* and *dt* were tuned keeping the time horizon less than a second and arriving at the highest possible *dt* giving the desired accuracy.