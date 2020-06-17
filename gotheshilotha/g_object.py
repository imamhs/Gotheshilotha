# Copyright (c) 2019, Md Imam Hossain (emamhd at gmail dot com)
# see LICENSE.txt for details

"""
Racing subject object
"""

from math import hypot, radians, sqrt
from obosthan import OVector2D
from numpy import polyfit, polyval
from scipy import spatial
from shuddo import S_moving_average_data

class GTS_object:

    def __init__(self, _stride_duration, _max_speed, _max_acceleration, _max_yaw_rate):
        self.time = []  # time list
        self.time_interval = 0.0    # time period between samples
        self.stride_duration = _stride_duration
        self.stride_cycle_steps = 0 # number of samples per stride
        self.sampling_rate = 0  # number of samples per second
        self.max_displacement = 0.0 # max displacement in a sample
        self.max_speed = _max_speed
        self.max_acceleration = _max_acceleration
        self.max_angular_displacement = 0.0
        self.max_yaw_rate = _max_yaw_rate
        self.coord = [] # X and Y coordinates
        self.displacement = []
        self.heading = []
        self.distance = []
        self.speed = []
        self.average_speed = []
        self.speed_model = [] # only available for adjusted data
        self.acceleration = []
        self.acceleration_model = []
        self.yaw_rate = []
        self.angular_displacement = []
        self.distance_to_lure = []
        self.lane_position = 0
        self.sampling_size = 0  # total number of samples
        self.curvature = [] # curvature of path
        self.radius_of_curvature = []   # radius of curvature of path
        self.adjusted_data = False  # wheather data is adjusted to match stride frequency or not
        self.distance_to_lure_path = []

    def clean_coord(self, _factor=2):
        x_coord, y_coord = list(zip(*self.coord))

        xnew_coord = S_moving_average_data(x_coord, _smoothing=_factor)
        ynew_coord = S_moving_average_data(y_coord, _smoothing=_factor)

        self.coord.clear()
        self.coord = list(zip(xnew_coord,ynew_coord))

    def calculate_dynamics_factors(self, _lure_points=[], _cal_lpd=False):

        total_distance = 0.0

        previous_heading = OVector2D(0, 0)
        current_heading = OVector2D(0, 0)

        for i in range(self.sampling_size):

            if (i == 0):
                self.displacement.append(0.0)
                self.distance.append(0.0)
            else:
                distance = hypot(self.coord[i][0]-self.coord[i-1][0], self.coord[i][1]-self.coord[i-1][1])
                if (distance >= self.max_displacement):
                    total_distance = total_distance + (self.distance[i-1]-self.distance[i-2])
                    self.displacement.append(self.displacement[i-1])
                    self.distance.append(total_distance)
                else:
                    total_distance = total_distance + distance
                    self.displacement.append(distance)
                    self.distance.append(total_distance)

            if (i == 0):
                self.speed.append(0.0)
                self.average_speed.append(0)
            else:
                speed = self.displacement[i] / self.time_interval
                if speed < 100.0:
                    self.speed.append(speed)
                else:
                    self.speed.append(self.speed[i - i])
                self.average_speed.append(self.distance[i] / self.time[i])

            if (i == 0 or i == 1):
                self.acceleration.append(0.0)
            else:
                acceleration = (self.speed[i] - self.speed[i - 1]) / self.time_interval
                if abs(acceleration) < self.max_acceleration:
                    self.acceleration.append(acceleration)
                else:
                    self.acceleration.append(self.acceleration[i - i])

            if (i == 0 or i == 1):
                self.angular_displacement.append(0.0)
                self.yaw_rate.append(0.0)
                if i == 1:
                    previous_heading.define_line(self.coord[i-1][0], self.coord[i-1][1], self.coord[i][0], self.coord[i][1])
                    self.heading.append(0.0)
                    self.heading.append(previous_heading.angle)
            else:
                current_heading.define_line(self.coord[i-1][0], self.coord[i-1][1], self.coord[i][0], self.coord[i][1])
                self.heading.append(current_heading.angle)
                angular_displacement = current_heading.angle - previous_heading.angle
                if abs(angular_displacement) > 300:
                    if current_heading.angle > 270:
                        angular_displacement = (360-current_heading.angle) + previous_heading.angle
                    elif previous_heading.angle > 270:
                        angular_displacement = -((360 - previous_heading.angle) + current_heading.angle)
                previous_heading[0] = current_heading[0]
                previous_heading[1] = current_heading[1]

                if abs(angular_displacement) < self.max_angular_displacement:
                    self.angular_displacement.append(angular_displacement)
                else:
                    self.angular_displacement.append(self.angular_displacement[i - i])

                yaw_rate = radians(angular_displacement) / self.time_interval
                if abs(yaw_rate) < self.max_yaw_rate:
                    self.yaw_rate.append(yaw_rate)
                else:
                    self.yaw_rate.append(self.yaw_rate[i - i])

            self.heading[0] = self.heading[1]
            
            if (i != 0 and i != (self.sampling_size - 1)):
                point1_x = self.coord[i-1][0]
                point1_y = self.coord[i-1][1]
                point2_x = self.coord[i][0]
                point2_y = self.coord[i][1]
                point3_x = self.coord[i+1][0]
                point3_y = self.coord[i+1][1]
                side_a = sqrt(((point2_x - point1_x) ** 2) + ((point2_y - point1_y) ** 2))
                side_b = sqrt(((point3_x - point2_x) ** 2) + ((point3_y - point2_y) ** 2))
                side_c = sqrt(((point3_x - point1_x) ** 2) + ((point3_y - point1_y) ** 2))
                semiperimeter = (0.5) * (side_a + side_b + side_c)
                delta = semiperimeter * (semiperimeter - side_a) * (semiperimeter - side_b) * (semiperimeter - side_c)
                triangle_area = 0
                if (delta > 0):
                    triangle_area = sqrt(delta)
                else:
                    triangle_area = 0
                if (triangle_area == 0):
                    self.radius_of_curvature.append(float('inf'))
                    self.curvature.append(0.0)
                else:
                    radius_of_curvature = (side_a * side_b * side_c) / (4 * triangle_area)
                    curvature = 1 / radius_of_curvature
                    if radius_of_curvature < 1000:
                        self.radius_of_curvature.append(radius_of_curvature)
                        self.curvature.append(curvature)
                    else:
                        self.radius_of_curvature.append(float('inf'))
                        self.curvature.append(0.0)
            else:
                self.radius_of_curvature.append(float('inf'))
                self.curvature.append(0.0)

            if _cal_lpd == True:
                nearby_lure_points_distances, nearby_lure_points_indices = spatial.KDTree(_lure_points).query((self.coord[i][0],self.coord[i][1]), 2)
                nearby_lure_point1 = (_lure_points[nearby_lure_points_indices[0]][0], _lure_points[nearby_lure_points_indices[0]][1])
                nearby_lure_point2 = (_lure_points[nearby_lure_points_indices[1]][0], _lure_points[nearby_lure_points_indices[1]][1])
                self.distance_to_lure_path.append(((((nearby_lure_point2[1]-nearby_lure_point1[1])*self.coord[i][0])-((nearby_lure_point2[0]-nearby_lure_point1[0])*self.coord[i][1])+(nearby_lure_point2[0]*nearby_lure_point1[1])-(nearby_lure_point2[1]*nearby_lure_point1[0]))/sqrt(((nearby_lure_point2[1]-nearby_lure_point1[1])**2)+((nearby_lure_point2[0]-nearby_lure_point1[0])**2))))

        if self.adjusted_data == True:

            p9 = polyfit(self.time, self.speed, 9)

            y_p = polyval(p9, self.time)
            self.speed_model = y_p
            self.speed_model[0] = 0.0

            for i in range(self.sampling_size):
                if (i == 0):
                    self.acceleration_model.append(0.0)
                else:
                    acceleration_model = (self.speed_model[i] - self.speed_model[i - 1]) / self.time_interval
                    if acceleration < self.max_acceleration:
                        self.acceleration_model.append(acceleration_model)
                    else:
                        self.acceleration_model.append(self.acceleration_model[i - i])