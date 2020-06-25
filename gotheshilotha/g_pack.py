# Copyright (c) 2019, Md Imam Hossain (emamhd at gmail dot com)
# see LICENSE.txt for details

"""
Racing subjects pack object
"""

from csv import reader
from math import hypot, ceil, degrees
from obosthan import OPoint2D
from numpy import polyfit, polyval, empty
from .g_object import GTS_object

class GTS_pack:

    def __init__(self, _max_distance):
        self.num_of_objects = 0
        self.sample_size = 0
        self.racing_objects = [] # index zero is lure object
        self.max_distance = _max_distance # max distance between objects
        self.centroid_coord = []
        self.average_centroid_distance = [] # average of objects' distances to the objects' centroid
        self.centroid_distance = [] # distance travelled by the objects' centroid
        self.lure_separation_distance = [] # distance to lure from the leading object
        self.average_lure_distance = [] # average of objects' distances to the lure
        self.average_objects_yaw_rate = []  # average yaw rate of objects' yaw rates
        self.average_objects_speed = [] # average speed of objects' speeds
        self.average_objects_curvature = [] # average curvature of objects' paths
        self.average_objects_heading = []
        self.average_objects_average_speed = []
        self.average_objects_speed_model = [] # only available for adjusted data

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

    def readjust_data_sampling(self):

        adjusted_indices = []

        i = 0
        while i < self.racing_objects[1].sampling_size:
            adjusted_indices.append(i)
            i += self.racing_objects[1].stride_cycle_steps

        for i in range(self.num_of_objects):
            time = []
            coord = []
            for ii in range(len(adjusted_indices)):
                time.append(self.racing_objects[i].time[adjusted_indices[ii]])
                coord.append(self.racing_objects[i].coord[adjusted_indices[ii]])
            self.racing_objects[i].time = time
            self.racing_objects[i].coord = coord
            self.racing_objects[i].sampling_size = len(self.racing_objects[i].time)
            self.racing_objects[i].time_interval = max(self.racing_objects[i].time) / self.racing_objects[i].sampling_size
            self.racing_objects[i].max_displacement = self.racing_objects[i].time_interval * self.racing_objects[i].max_speed
            self.racing_objects[i].max_angular_displacement = self.racing_objects[i].time_interval * degrees(self.racing_objects[i].max_yaw_rate)
            self.racing_objects[i].min_radius_of_curvature = self.racing_objects[i].max_speed / self.racing_objects[i].max_yaw_rate
            self.racing_objects[i].sampling_rate = int(1 / self.racing_objects[i].time_interval)
            self.racing_objects[i].stride_cycle_steps = ceil(self.racing_objects[i].stride_duration / self.racing_objects[i].time_interval)
            self.racing_objects[i].adjusted_data = True

        self.sample_size = len(self.racing_objects[0].time)

    def calculate_dynamics(self):

        nob = self.num_of_objects - 1

        centroid_distance = 0

        for i in range(self.sample_size):
            centroid_position_x = 0.0
            centroid_position_y = 0.0
            lure_separation_distances = []

            average_centroid_distance = 0.0

            average_lure_distance = 0.0

            total_speed = 0
            total_curvature = 0
            total_heading = 0
            total_speed_average = 0

            total_yaw_rate = 0

            for ii in range(nob):

                centroid_position_x = centroid_position_x + self.racing_objects[ii+1].coord[i][0]
                centroid_position_y = centroid_position_y + self.racing_objects[ii+1].coord[i][1]
                lure_separation_distance = hypot(self.racing_objects[0].coord[i][0] - self.racing_objects[ii + 1].coord[i][0], self.racing_objects[0].coord[i][1] - self.racing_objects[ii + 1].coord[i][1])
                if (lure_separation_distance < self.max_distance):
                    lure_separation_distances.append(lure_separation_distance)
                    self.racing_objects[ii + 1].distance_to_lure.append(lure_separation_distance)
                else:
                    if ii == 0:
                        lure_separation_distances.append(0.0)
                    else:
                        lure_separation_distances.append(lure_separation_distances[ii-1])
                    self.racing_objects[ii + 1].distance_to_lure.append(0.0)

                average_lure_distance = average_lure_distance + self.racing_objects[ii + 1].distance_to_lure[i]

                total_speed += self.racing_objects[ii + 1].speed[i]
                total_speed_average += self.racing_objects[ii + 1].average_speed[i]
                total_yaw_rate += self.racing_objects[ii + 1].yaw_rate[i]
                total_curvature += self.racing_objects[ii + 1].curvature[i]
                total_heading += self.racing_objects[ii + 1].heading[i]

            self.lure_separation_distance.append(min(lure_separation_distances))

            centroid_position_x = centroid_position_x / nob
            centroid_position_y = centroid_position_y / nob

            for ii in range(nob):
                average_centroid_distance += hypot(centroid_position_x - self.racing_objects[ii + 1].coord[i][0], centroid_position_y - self.racing_objects[ii + 1].coord[i][1])

            self.centroid_coord.append(OPoint2D(centroid_position_x, centroid_position_y))

            if (i == 0):
                self.centroid_distance.append(0.0)
            else:
                distance = hypot(self.centroid_coord[i][0]-self.centroid_coord[i-1][0], self.centroid_coord[i][1]-self.centroid_coord[i-1][1])
                if (distance >= self.racing_objects[1].max_displacement):
                    centroid_distance += self.centroid_distance[i-1]
                    self.centroid_distance.append(centroid_distance)
                else:
                    centroid_distance += distance
                    self.centroid_distance.append(centroid_distance)

            average_centroid_distance = average_centroid_distance / nob

            if average_centroid_distance < self.max_distance:
                self.average_centroid_distance.append(average_centroid_distance)
            else:
                self.average_centroid_distance.append(0.0)

            average_lure_distance = average_lure_distance / nob

            self.average_lure_distance.append(average_lure_distance)

            self.average_objects_speed.append(total_speed / nob)
            self.average_objects_average_speed.append(total_speed_average / nob)

            self.average_objects_yaw_rate.append(total_yaw_rate / nob)
            self.average_objects_curvature.append(total_curvature/nob)
            self.average_objects_heading.append(total_heading / nob)

        if self.racing_objects[1].adjusted_data == True:

            p9 = polyfit(self.racing_objects[1].time, self.average_objects_speed, 9)

            y_p = polyval(p9, self.racing_objects[1].time)
            self.average_objects_speed_model = y_p
            self.average_objects_speed_model[0] = 0.0

    def load_data(self, filename, _stride_duration=0.282, _max_speed=24, _max_acceleration=20, _max_yaw_rate=1, _min_time_period=2.0, _skip_samples_head=0, _skip_samples_tail=0):

        print("Trying loading race data ...")

        csv_in = open(filename, 'r')
        csv_reader = reader(csv_in)
        buffer_list = list(csv_reader)
        total_rows = len(buffer_list)
        row_counter = 1

        self.racing_objects.clear()
        self.centroid_coord.clear()
        self.average_centroid_distance.clear()
        self.centroid_distance.clear()
        self.lure_separation_distance.clear()
        self.average_lure_distance.clear()
        self.average_objects_speed.clear()
        self.average_objects_curvature.clear()
        self.average_objects_heading.clear()
        self.average_objects_speed_model.clear()

        try:

            for row in buffer_list:
                if len(self.racing_objects) == 0:

                    for i in range(int((len(row)-1)/2)):
                        self.racing_objects.append(GTS_object(_stride_duration, _max_speed, _max_acceleration, _max_yaw_rate))
                        if i > 0:
                            self.racing_objects[i].lane_position = int(row[i+i+1][1:2])

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
