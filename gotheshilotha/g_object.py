# Copyright (c) 2019-2021, Md Imam Hossain (emamhd at gmail dot com)
# see LICENSE.txt for details

"""
Racing subject object
"""

from math import hypot, radians, sqrt, ceil, degrees
from obosthan import OVector2D
from shuddo import S_moving_average_filter

class GTS_object:

    def __init__(self, _stride_duration, _max_speed, _max_acceleration, _max_yaw_rate, _max_radius_of_curvature=1000):
        self.stride_duration = _stride_duration
        self.max_speed = _max_speed
        self.max_acceleration = _max_acceleration
        self.max_yaw_rate = _max_yaw_rate
        self.max_radius_of_curvature = _max_radius_of_curvature  # considered as straight
        self.sampling_size = -1  # total number of data points
        self.time_interval = -1 # time period between data points
        self.max_displacement = -1  # max displacement between data points
        self.max_angular_displacement = -1
        self.min_radius_of_curvature = -1
        self.sampling_rate = -1  # number of data points per second
        self.stride_cycle_steps = -1  # number of data points per stride
        self.stripped_data = False  # flag for whether data points are stripped down to match stride frequency
        self.time = []  # time
        self.coord = [] # X and Y coordinates tuple
        self.displacement = []
        self.heading = []
        self.distance = []
        self.speed = []
        self.average_speed = []
        self.acceleration = []
        self.yaw_rate = []
        self.angular_displacement = []
        self.curvature = [] # path curvature
        self.radius_of_curvature = []   # path radius of curvature
        self.distance_to_lure = []
        self.average_distance_to_others = []
        self.offset_to_track = []
        self.start_position_identifier = 0

    def adjust_data_sampling(self):

        adjusted_indices = []

        i = 0
        while i < self.sampling_size:
            adjusted_indices.append(i)
            i += self.stride_cycle_steps

        time = []
        coord = []

        for i in range(len(adjusted_indices)):
            time.append(self.time[adjusted_indices[i]])
            coord.append(self.coord[adjusted_indices[i]])

        self.time = time
        self.coord = coord
        self.sampling_size = len(self.time)
        self.time_interval = max(self.time) / self.sampling_size
        self.max_displacement = self.time_interval * self.max_speed
        self.max_angular_displacement = self.time_interval * degrees(self.max_yaw_rate)
        self.min_radius_of_curvature = self.max_speed / self.max_yaw_rate
        self.sampling_rate = int(1 / self.time_interval)
        self.stride_cycle_steps = ceil(self.stride_duration / self.time_interval)
        self.stripped_data = True

    def clean_coord(self, _factor=2):

        x_coord, y_coord = list(zip(*self.coord))

        xnew_coord = S_moving_average_filter(x_coord, _smoothing=_factor)
        ynew_coord = S_moving_average_filter(y_coord, _smoothing=_factor)

        self.coord.clear()
        self.coord = list(zip(xnew_coord, ynew_coord))

    def clean_dynamics_results(self, _factor=2):

        xnew_displacement = S_moving_average_filter(self.displacement, _smoothing=_factor)
        xnew_heading = S_moving_average_filter(self.heading, _smoothing=_factor)
        xnew_speed = S_moving_average_filter(self.speed, _smoothing=_factor)
        xnew_acceleration = S_moving_average_filter(self.acceleration, _smoothing=_factor)
        xnew_yaw_rate = S_moving_average_filter(self.yaw_rate, _smoothing=_factor)
        xnew_angular_displacement = S_moving_average_filter(self.angular_displacement, _smoothing=_factor)
        xnew_curvature = S_moving_average_filter(self.curvature, _smoothing=_factor)
        xnew_radius_of_curvature = S_moving_average_filter(self.radius_of_curvature, _smoothing=_factor)

        self.displacement = xnew_displacement
        self.heading = xnew_heading
        self.speed = xnew_speed
        self.acceleration = xnew_acceleration
        self.yaw_rate = xnew_yaw_rate
        self.angular_displacement = xnew_angular_displacement
        self.curvature = xnew_curvature
        self.radius_of_curvature = xnew_radius_of_curvature

    def calculate_dynamics(self):

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
                    total_distance += self.max_displacement
                    self.displacement.append(self.max_displacement)
                    self.distance.append(total_distance)
                else:
                    total_distance += distance
                    self.displacement.append(distance)
                    self.distance.append(total_distance)

            if (i == 0):
                self.speed.append(0.0)
                self.average_speed.append(0.0)
            else:
                speed = self.displacement[i] / self.time_interval
                if speed < self.max_speed:
                    self.speed.append(speed)
                else:
                    self.speed.append(self.max_speed)
                self.average_speed.append(self.distance[i] / self.time[i])

            if (i == 0 or i == 1):
                self.acceleration.append(0.0)
            else:
                acceleration = (self.speed[i] - self.speed[i - 1]) / self.time_interval
                if abs(acceleration) < self.max_acceleration:
                    self.acceleration.append(acceleration)
                else:
                    self.acceleration.append((abs(acceleration)/acceleration)*self.max_acceleration)

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
                    self.angular_displacement.append((abs(angular_displacement)/angular_displacement)*self.max_angular_displacement)

                yaw_rate = radians(angular_displacement) / self.time_interval
                if abs(yaw_rate) < self.max_yaw_rate:
                    self.yaw_rate.append(yaw_rate)
                else:
                    self.yaw_rate.append((abs(yaw_rate)/yaw_rate)*self.max_yaw_rate)

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
                    if radius_of_curvature > self.min_radius_of_curvature and radius_of_curvature < self.max_radius_of_curvature:
                        self.radius_of_curvature.append(radius_of_curvature)
                        self.curvature.append(curvature)
                    else:
                        if radius_of_curvature < self.min_radius_of_curvature:
                            self.radius_of_curvature.append(self.min_radius_of_curvature)
                            self.curvature.append(1/self.min_radius_of_curvature)
                        elif radius_of_curvature > self.max_radius_of_curvature:
                            self.radius_of_curvature.append(self.max_radius_of_curvature)
                            self.curvature.append(1 / self.max_radius_of_curvature)
            else:
                self.radius_of_curvature.append(0.0)
                self.curvature.append(0.0)

        self.heading[0] = self.heading[1]

        if self.stripped_data == True:
            pass