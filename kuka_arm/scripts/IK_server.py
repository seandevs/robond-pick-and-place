#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *
from sympy.matrices import Matrix
import numpy as np
from numpy import array
import pprint

pp = pprint.PrettyPrinter(indent=4)


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
		
        ### Your FK code here
        # Create symbols
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')
        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
        theta_z, theta_y = symbols('theta_z theta_y') # for correction
	#
	#   
	# Create Modified DH parameters
	#
	#            
	# Define Modified DH Transformation matrix
	#
	#
        dh = {alpha0:      0,       a0: 0,   d1: 0.75,        q1: q1,
              alpha1: -pi/2.,    a1: 0.35,      d2: 0,  q2: q2-pi/2., 
              alpha2:      0,    a2: 1.25,      d3: 0,        q3: q3,
              alpha3: -pi/2.,  a3: -0.054,   d4: 1.50,        q4: q4, 
              alpha4:  pi/2.,       a4: 0,      d5: 0,        q5: q5,
              alpha5: -pi/2.,       a5: 0,      d6: 0,        q6: q6,
              alpha6:      0,       a6: 0,      d7: 0.303,    q7: 0}
	# Create individual transformation matrices
	#
	#
        def transform(alpha, a, d, q):
            return Matrix([[           cos(q),           -sin(q),           0,             a],
                           [sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                           [sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
                           [                0,                 0,           0,             1]]) 
	# Extract rotation matrices from the transformation matrices
	#
	#
        ###
        T0_1 = transform(alpha0, a0, d1, q1).subs(dh)
        T1_2 = transform(alpha1, a1, d2, q2).subs(dh)
        T2_3 = transform(alpha2, a2, d3, q3).subs(dh)
        T3_4 = transform(alpha3, a3, d4, q4).subs(dh)
        T4_5 = transform(alpha4, a4, d5, q5).subs(dh)
        T5_6 = transform(alpha5, a5, d6, q6).subs(dh)
        T6_G = transform(alpha6, a6, d7, q7).subs(dh)

        T0_2 = T0_1 * T1_2
        T0_3 = T0_2 * T2_3
        T0_4 = T0_3 * T3_4
        T0_5 = T0_4 * T4_5
        T0_6 = T0_5 * T5_6
        T0_G = T0_6 * T6_G


        # Initialize service response
        '''
        correctional transform
        use actual values instead of symbols for performance
        '''

        def R_x(roll): 
            return Matrix([[          1,          0,           0,    0],
                           [          0,  cos(roll),  -sin(roll),    0],
                           [          0,  sin(roll),   cos(roll),    0],
                           [          0,          0,           0,    1]])

        # use -90 degrees for q -np.pi/2
        def R_y(pitch):
            return Matrix([[ cos(pitch),          0,  sin(pitch),    0],
                           [          0,          1,           0,    0],
                           [-sin(pitch),          0,  cos(pitch),    0],
                           [          0,          0,           0,    1]])

        # use 180 degrees for q np.pi
        def R_z(yaw): 
            return Matrix([[   cos(yaw),  -sin(yaw),           0,    0],
                           [   sin(yaw),   cos(yaw),           0,    0],
                           [          0,          0,           1,    0],
                           [          0,          0,           0,    1]])

        R_fix = R_z(theta_z) * R_y(theta_y)

        Total_transform = T0_G * R_fix
        # pp.pprint(Total_transform.evalf(subs={q1: 1.81, q2: 0.53, q3: -0.22, q4:
        #    .37, q5: .36, q6:4.98, theta_y: -np.pi/2., theta_z: np.pi }))
        

        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

	    # Extract end-effector position and orientation from request
	    # px,py,pz = end-effector position
	    # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])
     
            ### Your IK code here 
	    # Compensate for rotation discrepancy between DH parameters and Gazebo
	    #
	    #
            r, p, y = symbols('r p y') 
            EE = Matrix([[px], [py], [pz]])
            d = dh[d7]
            Rrpy = R_z(y)[0:-1, 0:-1] * R_y(p)[0:-1, 0:-1] * R_x(r)[0:-1, 0:-1] * R_fix[0:-1, 0:-1]
            Rrpy = Rrpy.evalf(subs={r: roll, p: pitch, y: yaw, theta_z: np.pi, 
                    theta_y: -np.pi/2.})
            WC = EE - d * Rrpy * Matrix([0, 0, 1])


	    # Calculate joint angles using Geometric IK method
	    #
	    #
            ###
            # find thetas
            x = WC[0]
            y = WC[1] 
            z = WC[2] - dh[d1]

            r = sqrt((x*x + y*y) - dh[a1])
            A = 1.501 # from measurement
            B = sqrt(r*r + z*z)
            C = 1.25 # from measurement

            angle_a = acos((B*B + C*C - A*A) / (2.*B*C))
            angle_b = acos((C*C + A*A - B*B) / (2.*C*A))

            theta_1 = atan2(WC[1], WC[2])
            theta_2 = np.pi/2. - (angle_a + atan2(z, r))
            theta_3 = np.pi/2. - (angle_b + (1 + dh[a3])) # to account for offset in original diagram

            #calculate R3_6
            R0_6 = Rrpy
            R0_3 = T0_1[0:3, 0:3] * T1_2[0:3, 0:3] * T2_3[0:3, 0:3]
            R0_3 = R0_3.evalf(subs={q1: theta_1, q2: theta_2, q3: theta_3})
            R3_6 = R0_3.inv('LU') * R0_6
            pp.pprint("R3_6 {}".format(R3_6))
            sys.stdout.flush()

            #https://gamedev.stackexchange.com/questions/50963/how-to-extract-euler-angles-from-transformation-matrix
            theta_4 = atan2(R3_6[2,2], -R3_6[0,2])
            theta_5 = atan2(sqrt(R3_6[0,2]*R3_6[0,2] + R3_6[2,2]*R3_6[2,2]), R3_6[1,2])
            theta_6 = atan2(-R3_6[1,1], R3_6[1,0])
            #theta_5 = atan2(sqrt(R3_6[1,1]*R3_6[1,1] + R3_6[1,0]*R3_6[1,0]), R3_6[1,2])
		
            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
	    joint_trajectory_point.positions = [theta_1, theta_2, theta_3,
                    theta_4, theta_5, theta_6]
	    joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
