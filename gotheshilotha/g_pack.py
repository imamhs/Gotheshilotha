# Copyright (c) 2019-2021, Md Imam Hossain (emamhd at gmail dot com)
# see LICENSE.txt for details

"""
Racing subjects pack object
"""

from csv import reader
from math import hypot, ceil, degrees
from obosthan import OPoint2D
from .g_object import GTS_object
from shuddo import S_get_cluster_centroid_data, S_find_sample_distance_values

class GTS_pack:

    def __init__(self, _max_distance):
        self.SPEED_OF_LIGHT = 299792458
        self.max_distance = _max_distance  # max separation distance between objects
        self.num_of_objects = -1
        self.sample_size = -1
        self.stripped_data = False  # flag for whether data points are stripped down to match stride frequency
        self.racing_objects = [] # list of racers
        self.centroid_coord = []
        self.centroid_distance = [] # distance travelled by the objects' centroid
        self.lure_separation_distance = []  # distance to lure from the leading object
        self.average_lure_distance = []  # average of objects' distances to the lure
        self.average_centroid_distance = [] # average of objects' distances to the objects' centroid
        self.average_objects_distance = []  # average distance of objects' distances
        self.average_objects_maximum_distance = []  # average distance of objects' maximum distances
        self.average_objects_speed = [] # average speed of objects' speeds
        self.average_objects_average_speed = []
        self.average_objects_heading = []
        self.average_objects_yaw_rate = []  # average yaw rate of objects' yaw rates
        self.average_objects_curvature = []  # average curvature of objects' paths
        self.minmax_objects_speed = []  # minimum and maximum speed of objects
        self.minmax_objects_yaw_rate = [] # minimum and maximum yaw of objects
        self.minmax_objects_curvature = []  # minimum and maximum curvature of objects
        self.minmax_centroid_distance = []  # minimum and maximum curvature of objects

    def is_number(self, s):
        try:
            float(s)  # for int, long and float
        except ValueError:
            try:
                complex(s)  # for complex
            except ValueError:
                return False
        return True

    def clean_data(self):

        for i in range(self.num_of_objects):
            self.racing_objects[i].clean_coord()

    def clean_objects_dynamics(self):

        for i in range(self.num_of_objects):
            self.racing_objects[i].clean_dynamics_results()

    def calculate_objects_dynamics(self):

        for i in range(self.num_of_objects):
            self.racing_objects[i].calculate_dynamics()

    def adjust_data_sampling(self):

        for i in range(self.num_of_objects):
            self.racing_objects[i].adjust_data_sampling()

        self.sample_size = len(self.racing_objects[0].time)
        self.stripped_data = True

    def calculate_dynamics(self):

        nob = self.num_of_objects - 1

        centroid_distance = 0

        max_displacement = 0

        for ii in range(nob):
            if self.racing_objects[ii+1].max_displacement > max_displacement:
                max_displacement = self.racing_objects[ii+1].max_displacement

        for i in range(self.sample_size):

            racing_objects_coords = []

            distances = []
            total_speed = 0
            total_speed_average = 0
            total_heading = 0
            total_yaw_rate = 0
            total_curvature = 0

            min_speed = self.racing_objects[1].max_speed
            max_speed = -self.SPEED_OF_LIGHT

            min_curvature = 1000
            max_curvature = 0

            min_yaw_rate = self.racing_objects[1].max_yaw_rate
            max_yaw_rate = -self.SPEED_OF_LIGHT

            for ii in range(nob):
                racing_objects_coords.append((self.racing_objects[ii+1].coord[i][0], self.racing_objects[ii+1].coord[i][1]))

            average_lure_separation_distance, min_lure_separation_distance, max_lure_separation_distance, lure_separation_distances = S_find_sample_distance_values(racing_objects_coords, (self.racing_objects[0].coord[i][0], self.racing_objects[0].coord[i][1]))

            self.average_lure_distance.append(average_lure_separation_distance)
            self.lure_separation_distance.append(min_lure_separation_distance)

            centroid_position_x, centroid_position_y = S_get_cluster_centroid_data(racing_objects_coords)

            self.centroid_coord.append(OPoint2D(centroid_position_x, centroid_position_y))

            average_centroid_distance, min_centroid_distance, max_centroid_distance, centroid_distances = S_find_sample_distance_values(racing_objects_coords, (centroid_position_x, centroid_position_y))

            if average_centroid_distance < self.max_distance:
                self.average_centroid_distance.append(average_centroid_distance)
            else:
                self.average_centroid_distance.append(0.0)

            self.minmax_centroid_distance.append((min_centroid_distance, max_centroid_distance))

            if i == 0:
                self.centroid_distance.append(0.0)
            else:
                distance = hypot(self.centroid_coord[i][0]-self.centroid_coord[i-1][0], self.centroid_coord[i][1]-self.centroid_coord[i-1][1])
                if distance >= max_displacement:
                    centroid_distance += max_displacement
                    self.centroid_distance.append(centroid_distance)
                else:
                    centroid_distance += distance
                    self.centroid_distance.append(centroid_distance)

            for ii in range(nob):

                if lure_separation_distances[ii] < self.max_distance:
                    self.racing_objects[ii+1].distance_to_lure.append(lure_separation_distances[ii])
                else:
                    self.racing_objects[ii+1].distance_to_lure.append(0.0)

                distances.append(self.racing_objects[ii+1].distance[i])
                total_speed += self.racing_objects[ii+1].speed[i]
                total_speed_average += self.racing_objects[ii+1].average_speed[i]
                total_yaw_rate += self.racing_objects[ii+1].yaw_rate[i]
                total_curvature += self.racing_objects[ii+1].curvature[i]
                total_heading += self.racing_objects[ii+1].heading[i]

                if self.racing_objects[ii+1].speed[i] < min_speed:
                    min_speed = self.racing_objects[ii+1].speed[i]

                if self.racing_objects[ii+1].speed[i] > max_speed:
                    max_speed = self.racing_objects[ii+1].speed[i]

                if self.racing_objects[ii+1].curvature[i] < min_curvature:
                    min_curvature = self.racing_objects[ii+1].curvature[i]

                if self.racing_objects[ii+1].curvature[i] > max_curvature:
                    max_curvature = self.racing_objects[ii+1].curvature[i]

                if self.racing_objects[ii+1].yaw_rate[i] < min_yaw_rate:
                    min_yaw_rate = self.racing_objects[ii+1].yaw_rate[i]

                if self.racing_objects[ii+1].yaw_rate[i] > max_yaw_rate:
                    max_yaw_rate = self.racing_objects[ii+1].yaw_rate[i]

            self.average_objects_distance.append(sum(distances) / nob)
            self.average_objects_maximum_distance.append(max(distances))
            self.average_objects_speed.append(total_speed / nob)
            self.average_objects_average_speed.append(total_speed_average / nob)

            self.average_objects_yaw_rate.append(total_yaw_rate / nob)
            self.average_objects_curvature.append(total_curvature/nob)
            self.average_objects_heading.append(total_heading / nob)

            self.minmax_objects_speed.append((min_speed, max_speed))
            self.minmax_objects_curvature.append((min_curvature, max_curvature))
            self.minmax_objects_yaw_rate.append((min_yaw_rate, max_yaw_rate))

        if self.stripped_data:
            pass

    def load_race_data(self, filename, _stride_duration=0.282, _max_speed=24, _max_acceleration=20, _max_yaw_rate=1, _min_time_period=2.0, _skip_samples_head=0, _skip_samples_tail=0):

        print("Trying loading race data ...")

        csv_in = open(filename, 'r')
        csv_reader = reader(csv_in)
        buffer_list = list(csv_reader)
        total_rows = len(buffer_list)
        row_counter = 1

        self.racing_objects.clear()
        self.centroid_coord.clear()
        self.centroid_distance.clear()
        self.lure_separation_distance.clear()
        self.average_lure_distance.clear()
        self.average_centroid_distance.clear()
        self.average_objects_distance.clear()
        self.average_objects_maximum_distance.clear()
        self.average_objects_speed.clear()
        self.average_objects_average_speed.clear()
        self.average_objects_heading.clear()
        self.average_objects_yaw_rate.clear()
        self.average_objects_curvature.clear()
        self.minmax_objects_speed.clear()
        self.minmax_objects_yaw_rate.clear()
        self.minmax_objects_curvature.clear()
        self.minmax_centroid_distance.clear()

        try:

            for row in buffer_list:
                if len(self.racing_objects) == 0:

                    for i in range(int((len(row)-1)/2)):
                        self.racing_objects.append(GTS_object(_stride_duration, _max_speed, _max_acceleration, _max_yaw_rate))
                        if i > 0:
                            self.racing_objects[i].start_position_identifier = int(row[i+i+1][1:2])

                    self.num_of_objects = len(self.racing_objects)

                else:

                    if (row_counter >= _skip_samples_head and row_counter <= (total_rows - _skip_samples_tail)):

                        for i in range(len(self.racing_objects)):
                            self.racing_objects[i].time.append(float(row[0]))
                            self.racing_objects[i].coord.append((float(row[i+i+1]), float(row[i+i+2])))

                row_counter = row_counter + 1

        except:
            csv_in.close()
            return False

        csv_in.close()

        if self.racing_objects[0].time[-1] < _min_time_period:
            self.racing_objects.clear()
            self.num_of_objects = 0
            return False

        self.num_of_objects = len(self.racing_objects)
        self.sample_size = len(self.racing_objects[0].time)

        for i in range(self.num_of_objects):
            self.racing_objects[i].sampling_size = len(self.racing_objects[i].time)
            self.racing_objects[i].time_interval = max(self.racing_objects[i].time) / self.racing_objects[i].sampling_size
            self.racing_objects[i].max_displacement = self.racing_objects[i].time_interval * self.racing_objects[i].max_speed
            self.racing_objects[i].max_angular_displacement = self.racing_objects[i].time_interval * degrees(self.racing_objects[i].max_yaw_rate)
            self.racing_objects[i].min_radius_of_curvature = self.racing_objects[i].max_speed / self.racing_objects[i].max_yaw_rate
            self.racing_objects[i].sampling_rate = int(1 / self.racing_objects[i].time_interval)
            self.racing_objects[i].stride_cycle_steps = ceil(self.racing_objects[i].stride_duration / self.racing_objects[i].time_interval)

        return True
