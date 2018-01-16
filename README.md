# Cornell Modelling Competition in Mathematics (CMCM) 2017
This manuscript was our submission to the CMCM 2017 competition, a Cornell wide competition to solve an open ended problem in the real world using modelling techniques.

## Abstract

In this manuscript, we present a model that simulates the diesel fuel usage of hybrid and
diesel buses on certain TCAT bus routes in Ithaca. Our model makes use of elevation data
to calculate the physical power exerted by the engine at a series of discrete points along the
route. When modeling the hybrid bus, we take into consideration the proper conditions for
battery power to be used, and use this information to determine the total energy consumed by
hybrid and diesel buses along each route. We then solve the allocation problem using linear
programming to find the optimal distribution of hybrid buses in order to maximize fuel savings
and maintain the current number of buses along each route. We determined that the optimal
allocation is two hybrid buses on route 10, two buses on route 11, one bus on route 15, one
bus on route 17, and two buses on route 81. TCAT can save $28,630 when compared to an
allocation of no hybrid buses to these routes and can save $12,411 when compared to the least
optimal allocation of hybrid buses. Finally, we discuss further extensions and improvements
to our model to improve its accuracy.

The final manuscript is under the title _CMCM Manuscript - Optimum Allocation of TCAT Hybrid Buses to Maximize Fuel Efficiency_.
Other files are supplementary.
<p align="center">
<img src="https://github.com/RaghavBatra/CMCM_2017/blob/master/map_routes_terrain.png" width="750" height="500">
</p>

## More information on the below is contained in the manuscript
* *route_10_0* and *route_10_5* are folders that contain graphs of battery consumption, fuel consumption and electric consumption
* *routes_csv_data* is a folder that contains the elevation data for each route in a CSV format
* *slope* is a folder that computes the elevation slope for each route, stored as a CSV
* *smoothing_elev_vs_dist* is a folder that contains the smoothened slope: elev_profile_smooth_10/17 are images of such plots
* *find_and_smoothen_slope.m* is a MATLAB script that smoothens the initial data
* *find_optimal.m* is a MATLAB script to solve the bus allocation problem using linear programming

 
 ## Contributors
* Raghav Batra (rb698)
* Michael Galbato (mag386)
* Red Guiliano (rg586) 
