"""Provider for a set of snapshots from a MD trajectory."""
from .provider import Provider
import os
from shutil import copyfile
import ase.io
import ase
import ase.io
import numpy as np
from scipy.spatial import distance
from asap3.analysis.rdf import RadialDistributionFunction


class SnapshotsProvider(Provider):
    """
    Filters snapshots from a given MD trajectory, with a user specified metric.

    Parameters
    ----------
    parameters : malada.utils.parametes.Parameters
        Parameters used to create this object.

    external_snapshots : string
        Path to a trajectory file containing snapshots. If not None,
        no parsing will be done.
    """

    def __init__(self, parameters, external_snapshots=None):
        super(SnapshotsProvider, self).__init__(parameters)
        self.external_snapshots = external_snapshots
        self.snapshot_file = None
        self.fitted_trajectory = None
        self.__equilibrated_rdf = None

    def provide(self, provider_path, trajectoryfile, temperaturefile):
        """
        Provide a trajectory file containing atomic snapshots.

        Parameters
        ----------
        provider_path : string
            Path in which to operate in.

        trajectoryfile : string
            Path to file containing the MD trajectory as ASE trajectory.

        temperaturefile : string
            File containing the temperatures from the MD run as numpy array.
        """
        file_name = self.parameters.element + \
                    str(self.parameters.number_of_atoms) + \
                    "_" + self.parameters.crystal_structure +\
                    "_" + str(self.parameters.temperature) +\
                    "_possible_snapshots"
        self.snapshot_file = os.path.join(provider_path, file_name+".traj")
        iteration_numbers = os.path.join(provider_path, file_name+".md_iterations.npy")
        md_trajectory = ase.io.trajectory.Trajectory(trajectoryfile)
        temperatures = np.load(temperaturefile)
        if self.external_snapshots is None:
            # Find out where to begin with the parsing.
            first_snapshot = self.__get_first_snapshot(md_trajectory)

            # Now determine the value for the distance metric.
            distance_metric = self.__determine_distance_metric(md_trajectory)

            # Now parse the entire trajectory.
            self.__parse_trajectory(md_trajectory, first_snapshot,
                                    temperatures, distance_metric,
                                    self.snapshot_file, iteration_numbers)
        else:
            copyfile(self.external_snapshots, self.snapshot_file)
            print("Getting <<snapshots>>.npy"
                  " files from disc.")

    def __get_first_snapshot(self, trajectory):
        if self.parameters.snapshot_parsing_beginning < 0:
            # In this case we automatically fit a function to the trajectory
            # and use that to detect from where to begin parsing.
            self.__fit_trajectory(trajectory)

        else:
            return self.parameters.snapshot_parsing_beginning

    def __determine_distance_metric(self, trajectoryfile):
        if self.parameters.distance_metric_snapshots_cutoff < 0 and \
                self.parameters.snapshot_parsing_criterion != "random":
            raise Exception(
                "Automatic detection of distance cutoff not supported")
        else:
            return self.parameters.distance_metric_snapshots_cutoff

    def __parse_trajectory(self, trajectory, beginning_snapshot, temperatures,
                           distance_metric, filename_traj,
                           filename_numbers):
        allowed_temp_diff_K = (self.parameters.
                               snapshot_parsing_temperature_tolerance_percent
                               / 100) * self.parameters.temperature
        current_snapshot = beginning_snapshot
        begin_snapshot = beginning_snapshot+1
        end_snapshot = len(trajectory)
        j = 0
        md_iteration = []
        for i in range(begin_snapshot, end_snapshot):
            if self.__check_if_snapshot_is_valid(trajectory[i],
                                                 temperatures[i],
                                                 trajectory[current_snapshot],
                                                 temperatures[current_snapshot],
                                                 distance_metric,
                                                 allowed_temp_diff_K):
                current_snapshot = i
                md_iteration.append(current_snapshot)
                j+=1
        np.random.shuffle(md_iteration)
        for i in range(0, len(md_iteration)):
            if i == 0:
                traj_writer = ase.io.trajectory.TrajectoryWriter(filename_traj, mode='w')
            else:
                traj_writer = ase.io.trajectory.TrajectoryWriter(filename_traj, mode='a')
            atoms_to_write = self._enforce_pbc(trajectory[md_iteration[i]])
            traj_writer.write(atoms=atoms_to_write)
        np.save(filename_numbers, md_iteration)
        print(j, "possible snapshots found in MD trajectory.")
        if j < self.parameters.number_of_snapshots:
            raise Exception("Not enough snapshots found in MD trajectory. "
                            "Please run a longer MD calculation.")

    def __check_if_snapshot_is_valid(self, snapshot_to_test, temp_to_test,
                                     reference_snapshot, reference_temp,
                                     distance_metric,
                                     allowed_temp_diff):
        distance = self.__calculate_distance_between_snapshots(snapshot_to_test, reference_snapshot)
        temp_diff = np.abs(temp_to_test-reference_temp)
        if distance > distance_metric and temp_diff < allowed_temp_diff:
            return True
        else:
            return False

    def __calculate_distance_between_snapshots(self, snapshot1, snapshot2,
                                               dist_type="simple_geometric",
                                               save_rdf1 = False):
        if self.parameters.distance_metric_snapshots == "realspace":
            positions1 = snapshot1.get_positions()
            positions2 = snapshot2.get_positions()
            if self.parameters.distance_metric_snapshots_reduction == "minimal_distance":
                result = np.amin(distance.cdist(positions1, positions2), axis=0)
                result = np.mean(result)

            elif self.parameters.distance_metric_snapshots_reduction == "cosine_distance":
                number_of_atoms = snapshot1.get_number_of_atoms()
                result = distance.cosine(np.reshape(positions1, [number_of_atoms*3]),
                                         np.reshape(positions2, [number_of_atoms*3]))

            else:
                raise Exception("Unknown distance metric reduction.")
        elif self.parameters.distance_metric_snapshots == "rdf":
            rng = np.min(
                np.linalg.norm(snapshot1.get_cell(), axis=0)) - self.\
                parameters.distance_metric_snapshots_rdf_tolerance
            if save_rdf1 is True:
                if self.__equilibrated_rdf is None:
                    self.__equilibrated_rdf = RadialDistributionFunction(snapshot1, \
                                                  rng, \
                                                  self.parameters.distance_metric_snapshots_rdf_bins).get_rdf()
                rdf1 = self.__equilibrated_rdf
            else:
                rdf1 = RadialDistributionFunction(snapshot1, \
                                                  rng, \
                                                  self.parameters.distance_metric_snapshots_rdf_bins).get_rdf()
            rdf2 = RadialDistributionFunction(snapshot2, \
                                              rng, \
                                              self.parameters.distance_metric_snapshots_rdf_bins).get_rdf()

            if self.parameters.distance_metric_snapshots_reduction == "minimal_distance":
                raise Exception("Combination of distance metric and reduction "
                                "not supported.")

            elif self.parameters.distance_metric_snapshots_reduction == "cosine_distance":
                result = distance.cosine(rdf1, rdf2)

            else:
                raise Exception("Unknown distance metric reduction.")
        else:
            raise Exception("Unknown distance metric selected.")

        return result

    def __denoise(self, signal):
        denoised_signal = np.convolve(signal, np.ones(
            self.parameters.distance_metrics_denoising_width) / self.parameters.distance_metrics_denoising_width, mode='same')
        return denoised_signal

    def fit_trajectory(self, trajectory):
        # First, we ned to calculate the reduced metrics for the trajectory.
        # For this, we calculate the distance between all the snapshots
        # and the last one.
        distance_to_last = []
        for idx, step in enumerate(trajectory):
            distance_to_last.append(self.__calculate_distance_between_snapshots(trajectory[-1], step, save_rdf1=True))

        # Now, we denoise the distance metrics.
        distance_to_last = self.__denoise(distance_to_last)

        return distance_to_last




