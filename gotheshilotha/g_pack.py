# Copyright (c) 2019-2022, Md Imam Hossain (emamhd at gmail dot com)
# see LICENSE.txt for details

"""
Racing subjects pack object
"""

from csv import reader, writer
from math import hypot, ceil, degrees
from numpy import empty
from scipy.spatial import KDTree
from obosthan import OPoint2D, OLine2D
from shuddo import S_get_cluster_centroid_data, S_find_sample_distance_data, S_standard_deviation_values
from .g_object import GTS_object

class GTS_pack:

    def __init__(self, _max_distance):
        self.SPEED_OF_LIGHT = 299792458.0
        self.max_distance = _max_distance  # max separation distance between objects
        self.num_of_objects = -1
        self.sample_size = -1
        self.stripped_data = False  # flag for whether data points are stripped down to match stride frequency
        self.racing_objects = [] # list of racers
        self.centroid_coord = []
        self.track_coord = []
        self.__track_coord = []
        self.__track_coord_tree = None
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
        self.average_objects_centrifugal_acceleration = []
        self.average_objects_offset_to_track = []  # average offset of objects' from the track path
        self.minmax_objects_speed = []  # minimum and maximum speed of objects
        self.minmax_objects_yaw_rate = [] # minimum and maximum yaw of objects
        self.minmax_objects_curvature = []  # minimum and maximum curvature of objects
        self.minmax_centroid_distance = []  # minimum and maximum distances to the objects' centroid
        self.minmax_objects_offset_to_track = []  # minimum and maximum offsets of objects' from the track path
        self.std_objects_speed = []  # std speed of objects
        self.std_objects_yaw_rate = []  # std yaw of objects
        self.std_objects_curvature = []  # std curvature of objects
        self.std_objects_heading = []

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

    def uniform_results_sampling(self):

        for i in range(self.num_of_objects):
            self.racing_objects[i].uniform_results_sampling()

    def calculate_dynamics_fundamental(self):

        nob = self.num_of_objects - 1

        for i in range(self.sample_size):

            racing_objects_speeds = []
            racing_objects_curvatures = []
            distances = []

            for ii in range(nob):
                racing_objects_speeds.append(self.racing_objects[ii+1].speed[i])
                racing_objects_curvatures.append(self.racing_objects[ii+1].curvature[i])
                distances.append(self.racing_objects[ii+1].distance[i])

            self.average_objects_distance.append(sum(distances) / nob)
            self.average_objects_speed.append(sum(racing_objects_speeds) / nob)
            self.average_objects_curvature.append(sum(racing_objects_curvatures) / nob)

    def calculate_dynamics(self):

        nob = self.num_of_objects - 1

        centroid_distance = 0

        max_displacement = 0

        for ii in range(nob):
            if self.racing_objects[ii+1].max_displacement > max_displacement:
                max_displacement = self.racing_objects[ii+1].max_displacement

        for i in range(self.sample_size):

            racing_objects_coords = []
            racing_objects_speeds = []
            racing_objects_yaw_rates = []
            racing_objects_curvatures = []
            racing_objects_centrifugal_accelerations = []
            racing_objects_headings = []

            distances = []
            total_speed_average = 0

            min_speed = self.racing_objects[1].speed_limit
            speed_limit = -self.SPEED_OF_LIGHT

            min_curvature = 1000
            max_curvature = 0

            min_yaw_rate = self.racing_objects[1].yaw_rate_limit
            yaw_rate_limit = -self.SPEED_OF_LIGHT

            for ii in range(nob):
                racing_objects_coords.append((self.racing_objects[ii+1].coord[i][0], self.racing_objects[ii+1].coord[i][1]))
                racing_objects_speeds.append(self.racing_objects[ii+1].speed[i])
                racing_objects_yaw_rates.append(self.racing_objects[ii+1].yaw_rate[i])
                racing_objects_curvatures.append(self.racing_objects[ii+1].curvature[i])
                racing_objects_centrifugal_accelerations.append(self.racing_objects[ii+1].centrifugal_acceleration[i])
                racing_objects_headings.append(self.racing_objects[ii+1].heading[i])

            self.std_objects_speed.append(S_standard_deviation_values(racing_objects_speeds))
            self.std_objects_yaw_rate.append(S_standard_deviation_values(racing_objects_yaw_rates))
            self.std_objects_curvature.append(S_standard_deviation_values(racing_objects_curvatures))
            self.std_objects_heading.append(S_standard_deviation_values(racing_objects_headings))

            average_lure_separation_distance, min_lure_separation_distance, max_lure_separation_distance, lure_separation_distances = S_find_sample_distance_data(racing_objects_coords, (self.racing_objects[0].coord[i][0], self.racing_objects[0].coord[i][1]))

            self.average_lure_distance.append(average_lure_separation_distance)
            self.lure_separation_distance.append(min_lure_separation_distance)

            centroid_position_x, centroid_position_y = S_get_cluster_centroid_data(racing_objects_coords)

            self.centroid_coord.append(OPoint2D(centroid_position_x, centroid_position_y))

            average_centroid_distance, min_centroid_distance, max_centroid_distance, centroid_distances = S_find_sample_distance_data(racing_objects_coords, (centroid_position_x, centroid_position_y))

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
                total_speed_average += self.racing_objects[ii+1].average_speed[i]

                if self.racing_objects[ii+1].speed[i] < min_speed:
                    min_speed = self.racing_objects[ii+1].speed[i]

                if self.racing_objects[ii+1].speed[i] > speed_limit:
                    speed_limit = self.racing_objects[ii+1].speed[i]

                if self.racing_objects[ii+1].curvature[i] < min_curvature:
                    min_curvature = self.racing_objects[ii+1].curvature[i]

                if self.racing_objects[ii+1].curvature[i] > max_curvature:
                    max_curvature = self.racing_objects[ii+1].curvature[i]

                if self.racing_objects[ii+1].yaw_rate[i] < min_yaw_rate:
                    min_yaw_rate = self.racing_objects[ii+1].yaw_rate[i]

                if self.racing_objects[ii+1].yaw_rate[i] > yaw_rate_limit:
                    yaw_rate_limit = self.racing_objects[ii+1].yaw_rate[i]

            self.average_objects_distance.append(sum(distances) / nob)
            self.average_objects_maximum_distance.append(max(distances))
            self.average_objects_speed.append(sum(racing_objects_speeds) / nob)
            self.average_objects_average_speed.append(total_speed_average / nob)

            self.average_objects_yaw_rate.append(sum(racing_objects_yaw_rates) / nob)
            self.average_objects_curvature.append(sum(racing_objects_curvatures) /nob)
            self.average_objects_centrifugal_acceleration.append(sum(racing_objects_centrifugal_accelerations) / nob)
            self.average_objects_heading.append(sum(racing_objects_headings) / nob)

            self.minmax_objects_speed.append((min_speed, speed_limit))
            self.minmax_objects_curvature.append((min_curvature, max_curvature))
            self.minmax_objects_yaw_rate.append((min_yaw_rate, yaw_rate_limit))

        if self.stripped_data and (self.__track_coord_tree != None):

            for i in range(self.sample_size):
                for ii in range(nob):

                    projected_track_segment_end_distance, projected_track_segment_point_index = self.__track_coord_tree.query(self.racing_objects[ii+1].coord[i], 2)
                    self.racing_objects[ii+1].offset_to_track.append(OLine2D(float(self.__track_coord[projected_track_segment_point_index[0], 0]), float(self.__track_coord[projected_track_segment_point_index[0], 1]), float(self.__track_coord[projected_track_segment_point_index[1], 0]), float(self.__track_coord[projected_track_segment_point_index[1], 1])).distance_to_point(self.racing_objects[ii+1].coord[i]))

                racing_objects_offset_to_track = []

                for ii in range(nob):
                    racing_objects_offset_to_track.append(self.racing_objects[ii + 1].offset_to_track[i])

                self.average_objects_offset_to_track.append(sum(racing_objects_offset_to_track) / nob)
                self.minmax_objects_offset_to_track.append((min(racing_objects_offset_to_track), max(racing_objects_offset_to_track)))

    def load_track_coord(self, filename):
        print("Trying loading track coord data ...")

        csv_in = open(filename, 'r')
        csv_reader = reader(csv_in)
        buffer_list = list(csv_reader)
        total_rows = len(buffer_list)
        row_counter = 1

        if total_rows == 0:
            return False

        try:

            for row in buffer_list:

                self.track_coord.append((float(row[0]), float(row[1])))
                row_counter = row_counter + 1

        except:
            csv_in.close()
            return False

        csv_in.close()

        self.__track_coord = empty((total_rows, 2))

        for i in range(total_rows):
            self.__track_coord[i, 0] = self.track_coord[i][0]
            self.__track_coord[i, 1] = self.track_coord[i][1]

        if len(self.__track_coord) > 0:
            self.__track_coord_tree = KDTree(self.__track_coord)

        return True

    def dump_race_data(self, filename):

       csv_out = open(filename, 'w', newline='')
       csv_writer = writer(csv_out)

       header = []

       header.append("Time")
       header.append("Lure X")
       header.append("Lure Y")

       for i in range(self.num_of_objects - 1):
           header.append("Lane " + str(i+1) + " X")
           header.append("Lane " + str(i+1) + " Y")

       csv_writer.writerow(header)

       for i in range(self.sample_size):
           row = []
           row.append(self.racing_objects[0].time[i])
           for ii in range(self.num_of_objects):
               row.append(self.racing_objects[ii].coord[i][0])
               row.append(self.racing_objects[ii].coord[i][1])
           csv_writer.writerow(row)

       csv_out.close()

    def dump_coord(self, gr, filename):

       csv_out = open(filename, 'w', newline='')
       csv_writer = writer(csv_out)

       header = []

       header.append('Time(s)')
       header.append('Distance(m)')
       header.append('Coord X(m)')
       header.append('Coord Y(m)')

       csv_writer.writerow(header)

       for i in range(self.sample_size):
           row = []
           row.append(self.racing_objects[gr].time[i])
           row.append(self.racing_objects[gr].distance[i])
           row.append(self.racing_objects[gr].coord[i][0])
           row.append(self.racing_objects[gr].coord[i][1])

           csv_writer.writerow(row)

       csv_out.close()

    def dump_race_dynamics(self, gr, filename):

       csv_out = open(filename, 'w', newline='')
       csv_writer = writer(csv_out)

       header = []

       header.append('Time(s)')
       header.append('Distance(m)')
       header.append('Stride length(m)')
       header.append('Speed(m/s)')
       header.append('Acceleration(m/s2)')
       header.append('Yaw(rad)')

       csv_writer.writerow(header)

       for i in range(self.sample_size):
           row = []
           row.append(self.racing_objects[gr].time[i])
           row.append(self.racing_objects[gr].distance[i])
           row.append(self.racing_objects[gr].displacement[i])
           row.append(self.racing_objects[gr].speed[i])
           row.append(self.racing_objects[gr].acceleration[i])
           row.append(self.racing_objects[gr].yaw[i])

           csv_writer.writerow(row)

       csv_out.close()

    def load_race_data(self, filename, _stride_duration=0.282, _speed_limit=25, _max_acceleration=20, _yaw_rate_limit=1, _min_time_period=2.0, _skip_samples_head=0, _skip_samples_tail=0):

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
        self.average_objects_centrifugal_acceleration.clear()
        self.minmax_objects_speed.clear()
        self.minmax_objects_yaw_rate.clear()
        self.minmax_objects_curvature.clear()
        self.minmax_centroid_distance.clear()
        self.minmax_objects_offset_to_track.clear()
        self.std_objects_speed.clear()
        self.std_objects_yaw_rate.clear()
        self.std_objects_curvature.clear()
        self.std_objects_heading.clear()

        try:

            for row in buffer_list:
                if len(self.racing_objects) == 0:

                    for i in range(int((len(row)-1)/2)):
                        self.racing_objects.append(GTS_object(_stride_duration, _speed_limit, _max_acceleration, _yaw_rate_limit))
                        if i > 0:
                            self.racing_objects[i].start_position_identifier = int(row[i+i+1][1:2])

                    self.num_of_objects = len(self.racing_objects)

                else:

                    if (row_counter >= (_skip_samples_head+2)) and (row_counter <= (total_rows - _skip_samples_tail)):

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
            self.racing_objects[i].displacement_limit = self.racing_objects[i].time_interval * self.racing_objects[i].speed_limit
            self.racing_objects[i].angular_displacement_limit = self.racing_objects[i].time_interval * degrees(self.racing_objects[i].yaw_rate_limit)
            self.racing_objects[i].min_radius_of_curvature = self.racing_objects[i].speed_limit / self.racing_objects[i].yaw_rate_limit
            self.racing_objects[i].sampling_rate = int(1 / self.racing_objects[i].time_interval)
            self.racing_objects[i].stride_cycle_steps = ceil(self.racing_objects[i].stride_duration / self.racing_objects[i].time_interval)

        return True
