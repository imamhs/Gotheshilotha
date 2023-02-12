# Copyright (c) 2019-2023, Md Imam Hossain (emamhd at gmail dot com)
# see LICENSE.txt for details

"""
Racing subject object
"""

from math import hypot, radians, sqrt, ceil, degrees
from numpy import polyfit, polyval
from obosthan import OVector2D
from shuddo import S_moving_average_filter, S_uniform_spread_data

from .g_constants import UNDEFINED
from .g_constants import STRAIGHT_RADIUS

class GTS_object:

    def __init__(self, _stride_duration, _speed_limit, _acceleration_limit, _yaw_rate_limit, _radius_of_curvature_limit=STRAIGHT_RADIUS):

        self.stride_duration = _stride_duration
        self.speed_limit = _speed_limit
        self.acceleration_limit = _acceleration_limit
        self.yaw_rate_limit = _yaw_rate_limit
        self.radius_of_curvature_limit = _radius_of_curvature_limit  # this value is considered as straight
        self.sampling_size = UNDEFINED  # total number of data points
        self.time_interval = []  # time period between data points
        self.time_interval_avg = UNDEFINED  # average time period between data points
        self.displacement_limit = UNDEFINED
        self.angular_displacement_limit = UNDEFINED
        self.min_radius_of_curvature = UNDEFINED
        self.sampling_rate = UNDEFINED  # number of data points per second
        self.stride_cycle_steps = UNDEFINED  # number of data points per stride
        self.stripped_data = False  # flag for whether data points are stripped down to match stride frequency
        self.corrupted_data = False
        self.time = []  # time
        self.coord = []  # X and Y coordinates tuple
        self.displacement = []
        self.heading = []
        self.distance = []
        self.yaw = []
        self.speed = []
        self.average_speed = []
        self.acceleration = []
        self.yaw_rate = []
        self.angular_displacement = []
        self.curvature = []  # path curvature
        self.radius_of_curvature = []   # path radius of curvature
        self.centrifugal_acceleration = []
        self.distance_to_target = []
        self.average_distance_to_others = []
        self.min_distance_to_others = []
        self.offset_to_track = []
        self.start_position_identifier = 0

    def reset_data(self):  # clears dynamics calculation data

        self.sampling_size = UNDEFINED
        self.time_interval.clear()
        self.time_interval_avg = UNDEFINED
        self.displacement_limit = UNDEFINED
        self.angular_displacement_limit = UNDEFINED
        self.min_radius_of_curvature = UNDEFINED
        self.sampling_rate = UNDEFINED  # number of data points per second
        self.stride_cycle_steps = UNDEFINED  # number of data points per stride
        self.displacement.clear()
        self.heading.clear()
        self.distance.clear()
        self.yaw.clear()
        self.speed.clear()
        self.average_speed.clear()
        self.acceleration.clear()
        self.yaw_rate.clear()
        self.angular_displacement.clear()
        self.curvature.clear()
        self.radius_of_curvature.clear()
        self.centrifugal_acceleration.clear()
        self.distance_to_target.clear()
        self.average_distance_to_others.clear()
        self.min_distance_to_others.clear()
        self.offset_to_track.clear()

    def calculate_data_sampling(self):

        self.sampling_size = len(self.time)

        self.time_interval.clear()

        for i in range(self.sampling_size):
            if i == 0:
                self.time_interval.append(0)
            else:
                self.time_interval.append(self.time[i]-self.time[i-1])

        self.time_interval_avg = sum(self.time_interval[:-1]) / (self.sampling_size - 1)
        self.displacement_limit = self.time_interval_avg * self.speed_limit
        self.angular_displacement_limit = self.time_interval_avg * degrees(self.yaw_rate_limit)
        self.min_radius_of_curvature = self.speed_limit / self.yaw_rate_limit
        self.sampling_rate = int(1 / self.time_interval_avg)
        self.stride_cycle_steps = ceil(self.stride_duration / self.time_interval_avg)

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

        self.calculate_data_sampling()

        self.stripped_data = True

    def uniform_results_sampling(self):

        data_time = []
        data_speed = []
        data_curvature = []

        for i in range(self.sampling_size):

            data_time.append((self.distance[i], self.time[i]))
            data_speed.append((self.distance[i], self.speed[i]))
            data_curvature.append((self.distance[i], self.curvature[i]))

        udata_time = S_uniform_spread_data(data_time, self.sampling_size)
        udata_speed = S_uniform_spread_data(data_speed, self.sampling_size)
        udata_curvature = S_uniform_spread_data(data_curvature, self.sampling_size)

        self.sampling_size = len(udata_time)

        self.time = []
        self.distance = []
        self.speed = []
        self.curvature = []

        for i in range(self.sampling_size):

            self.time.append(udata_time[i][1])
            self.distance.append(udata_time[i][0])
            self.speed.append(udata_speed[i][1])
            self.curvature.append(udata_curvature[i][1])

    def clean_coord(self, _factor=2):

        x_coord, y_coord = list(zip(*self.coord))

        xnew_coord = S_moving_average_filter(x_coord, _smoothing=_factor)
        ynew_coord = S_moving_average_filter(y_coord, _smoothing=_factor)

        self.coord.clear()
        self.coord = list(zip(xnew_coord, ynew_coord))

    def clean_dynamics_results(self, _factor=2):

        xnew_displacement = S_moving_average_filter(self.displacement, _smoothing=_factor)
        xnew_speed = S_moving_average_filter(self.speed, _smoothing=_factor)
        xnew_heading = S_moving_average_filter(self.heading, _smoothing=_factor)
        xnew_acceleration = S_moving_average_filter(self.acceleration, _smoothing=_factor)
        xnew_angular_displacement = S_moving_average_filter(self.angular_displacement, _smoothing=_factor)
        xnew_yaw_rate = S_moving_average_filter(self.yaw_rate, _smoothing=_factor)
        xnew_curvature = S_moving_average_filter(self.curvature, _smoothing=_factor)
        xnew_radius_of_curvature = S_moving_average_filter(self.radius_of_curvature, _smoothing=_factor)
        xnew_centrifugal_acceleration = S_moving_average_filter(self.centrifugal_acceleration, _smoothing=_factor)

        self.displacement = xnew_displacement
        self.speed = xnew_speed
        self.heading = xnew_heading
        self.acceleration = xnew_acceleration
        self.angular_displacement = xnew_angular_displacement
        self.yaw_rate = xnew_yaw_rate
        self.curvature = xnew_curvature
        self.radius_of_curvature = xnew_radius_of_curvature
        self.centrifugal_acceleration = xnew_centrifugal_acceleration

    def calculate_dynamics(self):

        total_distance = 0.0
        total_yaw = 0.0
        previous_heading = OVector2D(0, 0)
        current_heading = OVector2D(0, 0)

        for i in range(self.sampling_size):

            if i == 0:

                self.displacement.append(0.0)
                self.distance.append(0.0)
                self.speed.append(0.0)
                self.average_speed.append(0.0)

            else:
                displacement = hypot(self.coord[i][0]-self.coord[i-1][0], self.coord[i][1]-self.coord[i-1][1])

                if displacement >= self.displacement_limit:
                    if len(self.displacement) > 2:
                        displacement_slope = (self.displacement[-1]-self.displacement[-2])/self.time_interval[i]
                        displacement = (displacement_slope*self.time_interval[i])+self.displacement[-1]
                    else:
                        displacement = 0.0

                total_distance += displacement

                self.displacement.append(displacement)
                self.distance.append(total_distance)

                speed = displacement / self.time_interval[i]

                if speed > self.speed_limit:
                    if len(self.speed) > 2:
                        speed_slope = (self.speed[-1]-self.speed[-2])/self.time_interval[i]
                        speed = (speed_slope*self.time_interval[i])+self.speed[-1]
                    else:
                        speed = 0.0

                self.speed.append(speed)
                self.average_speed.append(total_distance / self.time[i])

            if (i == 0) or (i == 1):

                self.acceleration.append(0.0)
                self.angular_displacement.append(0.0)
                self.yaw_rate.append(0.0)
                self.yaw.append(0.0)

                if i == 0:
                    self.heading.append(0.0)
                else:
                    previous_heading.define_line(self.coord[0][0], self.coord[0][1], self.coord[1][0], self.coord[1][1])
                    self.heading.append(previous_heading.angle)

            else:

                current_heading.define_line(self.coord[i-1][0], self.coord[i-1][1], self.coord[i][0], self.coord[i][1])
                self.heading.append(current_heading.angle)

                acceleration = (self.speed[-1] - self.speed[-2])/self.time_interval[i]

                if abs(acceleration) > self.acceleration_limit:
                    acceleration_slope = (self.acceleration[-1]-self.acceleration[-2])/self.time_interval[i]
                    acceleration = (acceleration_slope*self.time_interval[i])+self.acceleration[-1]

                self.acceleration.append(acceleration)

                angular_displacement = current_heading.angle - previous_heading.angle
                if abs(angular_displacement) > 300:
                    if current_heading.angle > 270:
                        angular_displacement = (360-current_heading.angle) + previous_heading.angle
                    elif previous_heading.angle > 270:
                        angular_displacement = -((360 - previous_heading.angle) + current_heading.angle)
                previous_heading[0] = current_heading[0]
                previous_heading[1] = current_heading[1]

                if abs(angular_displacement) > self.angular_displacement_limit:
                    angular_displacement_slope = (self.angular_displacement[-1]-self.angular_displacement[-2])/self.time_interval[i]
                    angular_displacement = (angular_displacement_slope*self.time_interval[i])+self.angular_displacement[-1]

                self.angular_displacement.append(angular_displacement)

                yaw_rate = radians(angular_displacement)/self.time_interval[i]
                if abs(yaw_rate) > self.yaw_rate_limit:
                    yaw_rate_slope = (self.yaw_rate[-1]-self.yaw_rate[-2])/self.time_interval[i]
                    yaw_rate = (yaw_rate_slope*self.time_interval[i])+self.yaw_rate[-1]

                self.yaw_rate.append(yaw_rate)
                total_yaw += abs(yaw_rate * self.time_interval[i])
                self.yaw.append(total_yaw)

            if (i != 0) and (i != (self.sampling_size - 1)):
                point1_x = self.coord[i-1][0]
                point1_y = self.coord[i-1][1]
                point2_x = self.coord[i][0]
                point2_y = self.coord[i][1]
                point3_x = self.coord[i+1][0]
                point3_y = self.coord[i+1][1]
                side_a = sqrt(((point2_x - point1_x)**2) + ((point2_y - point1_y)**2))
                side_b = sqrt(((point3_x - point2_x)**2) + ((point3_y - point2_y)**2))
                side_c = sqrt(((point3_x - point1_x)**2) + ((point3_y - point1_y)**2))
                semiperimeter = 0.5 * (side_a + side_b + side_c)
                delta = semiperimeter * (semiperimeter - side_a) * (semiperimeter - side_b) * (semiperimeter - side_c)
                triangle_area = 0.0
                if delta > 0:
                    triangle_area = sqrt(delta)
                else:
                    triangle_area = 0.0
                if triangle_area == 0.0:
                    self.radius_of_curvature.append(float('inf'))
                    self.curvature.append(0.0)
                    self.centrifugal_acceleration.append(0.0)
                else:
                    radius_of_curvature = (side_a * side_b * side_c) / (4 * triangle_area)
                    curvature = 1 / radius_of_curvature
                    if (radius_of_curvature >= self.min_radius_of_curvature) and (radius_of_curvature <= self.radius_of_curvature_limit):
                        self.radius_of_curvature.append(radius_of_curvature)
                        self.curvature.append(curvature)
                        self.centrifugal_acceleration.append((self.speed[i]**2)*curvature)
                    else:
                        if radius_of_curvature < self.min_radius_of_curvature:
                            self.radius_of_curvature.append(self.min_radius_of_curvature)
                            self.curvature.append(1 / self.min_radius_of_curvature)
                            self.centrifugal_acceleration.append((self.speed[i]**2) * self.curvature[-1])
                        elif radius_of_curvature > self.radius_of_curvature_limit:
                            self.radius_of_curvature.append(self.radius_of_curvature_limit)
                            self.curvature.append(1 / self.radius_of_curvature_limit)
                            self.centrifugal_acceleration.append((self.speed[i]**2) * self.curvature[-1])
            else:
                self.radius_of_curvature.append(0.0)
                self.curvature.append(0.0)
                self.centrifugal_acceleration.append(0.0)

        self.heading[0] = self.heading[1]

        radius_of_curvature_slope = (self.radius_of_curvature[1] - self.radius_of_curvature[2]) / self.time_interval[1]
        radius_of_curvature = (radius_of_curvature_slope * self.time_interval[1]) + self.radius_of_curvature[1]

        self.radius_of_curvature[0] = radius_of_curvature
        self.curvature[0] = (1 / self.radius_of_curvature[0])
        self.centrifugal_acceleration[0] = ((self.speed[0] ** 2) * self.curvature[0])

        radius_of_curvature_slope = (self.radius_of_curvature[-2] - self.radius_of_curvature[-3]) / self.time_interval[-2]
        radius_of_curvature = (radius_of_curvature_slope * self.time_interval[-2]) + self.radius_of_curvature[-2]

        self.radius_of_curvature[-1] = radius_of_curvature
        self.curvature[-1] = (1 / self.radius_of_curvature[-1])
        self.centrifugal_acceleration[-1] = ((self.speed[-1] ** 2) * self.curvature[-1])

        if self.stripped_data is True:

            p9 = polyfit(self.time, self.speed, 9)

            self.speed_model = polyval(p9, self.time)
            self.speed_model[0] = 0

            p9 = polyfit(self.time, self.acceleration, 9)

            self.acceleration_model = polyval(p9, self.time)
            self.acceleration_model[0] = 0

