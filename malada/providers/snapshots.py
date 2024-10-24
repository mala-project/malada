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
import pickle


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
        self.distance_metrics_denoised = None
        self.distance_metrics = None
        self.distances_realspace = None
        self.__saved_rdf = None
        self.first_snapshot = None
        self.distance_metric_cutoff = None
        self.average_distance_equilibrated = None
        self.first_considered_snapshot = None
        self.last_considered_snapshot = None

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
        file_name = (
            self.parameters.element
            + str(self.parameters.number_of_atoms)
            + "_"
            + self.parameters.crystal_structure
            + "_"
            + str(self.parameters.temperature)
            + "_possible_snapshots"
        )
        self.snapshot_file = os.path.join(provider_path, file_name + ".traj")
        iteration_numbers = os.path.join(
            provider_path, file_name + ".md_iterations.npy"
        )
        md_trajectory = ase.io.trajectory.Trajectory(trajectoryfile)
        temperatures = np.load(temperaturefile)
        if self.external_snapshots is None:
            # Find out where to begin with the parsing.
            self.first_snapshot = self.__get_first_snapshot(md_trajectory)

            # Now determine the value for the distance metric.
            self.distance_metric_cutoff = self.__determine_distance_metric(
                md_trajectory
            )

            # Now parse the entire trajectory.
            self.__parse_trajectory(
                md_trajectory,
                temperatures,
                self.snapshot_file,
                iteration_numbers,
            )
        else:
            copyfile(self.external_snapshots, self.snapshot_file)
            print("Getting <<snapshots>>.npy" " files from disc.")

    def __get_first_snapshot(self, trajectory):
        if self.parameters.snapshot_parsing_beginning < 0:
            # In this case we automatically fit a function to the trajectory
            # and use that to detect from where to begin parsing.
            return self.analyze_trajectory(trajectory)

        else:
            return self.parameters.snapshot_parsing_beginning

    def __determine_distance_metric(self, trajectory):
        if self.parameters.distance_metric_snapshots_cutoff < 0:
            return self.analyze_distance_metric(trajectory)
        else:
            return self.parameters.distance_metric_snapshots_cutoff

    def __parse_trajectory(
        self, trajectory, temperatures, filename_traj, filename_numbers
    ):
        allowed_temp_diff_K = (
            self.parameters.snapshot_parsing_temperature_tolerance_percent
            / 100
        ) * self.parameters.temperature
        current_snapshot = self.first_snapshot
        begin_snapshot = self.first_snapshot + 1
        end_snapshot = len(trajectory)
        j = 0
        md_iteration = []
        for i in range(begin_snapshot, end_snapshot):
            if self.__check_if_snapshot_is_valid(
                trajectory[i],
                temperatures[i],
                trajectory[current_snapshot],
                temperatures[current_snapshot],
                self.distance_metric_cutoff,
                allowed_temp_diff_K,
            ):
                current_snapshot = i
                md_iteration.append(current_snapshot)
                j += 1
        np.random.shuffle(md_iteration)
        for i in range(0, len(md_iteration)):
            if i == 0:
                traj_writer = ase.io.trajectory.TrajectoryWriter(
                    filename_traj, mode="w"
                )
            else:
                traj_writer = ase.io.trajectory.TrajectoryWriter(
                    filename_traj, mode="a"
                )
            atoms_to_write = self.enforce_pbc(trajectory[md_iteration[i]])
            traj_writer.write(atoms=atoms_to_write)
        np.save(filename_numbers, md_iteration)
        print(j, "possible snapshots found in MD trajectory.")
        if j < self.parameters.number_of_snapshots:
            raise Exception(
                "Not enough snapshots found in MD trajectory. "
                "Please run a longer MD calculation."
            )

    def __check_if_snapshot_is_valid(
        self,
        snapshot_to_test,
        temp_to_test,
        reference_snapshot,
        reference_temp,
        distance_metric,
        allowed_temp_diff,
    ):
        distance = self._calculate_distance_between_snapshots(
            snapshot_to_test,
            reference_snapshot,
            "realspace",
            "minimal_distance",
        )
        temp_diff = np.abs(temp_to_test - reference_temp)
        if distance > distance_metric and temp_diff < allowed_temp_diff:
            return True
        else:
            return False

    def _calculate_distance_between_snapshots(
        self, snapshot1, snapshot2, distance_metric, reduction, save_rdf1=False
    ):
        if distance_metric == "realspace":
            positions1 = snapshot1.get_positions()
            positions2 = snapshot2.get_positions()
            if reduction == "minimal_distance":
                result = np.amin(
                    distance.cdist(positions1, positions2), axis=0
                )
                result = np.mean(result)

            elif reduction == "cosine_distance":
                number_of_atoms = snapshot1.get_number_of_atoms()
                result = distance.cosine(
                    np.reshape(positions1, [number_of_atoms * 3]),
                    np.reshape(positions2, [number_of_atoms * 3]),
                )

            else:
                raise Exception("Unknown distance metric reduction.")
        elif distance_metric == "rdf":
            rng = (
                np.min(np.linalg.norm(snapshot1.get_cell(), axis=0))
                - self.parameters.distance_metric_snapshots_rdf_tolerance
            )
            if save_rdf1 is True:
                if self.__saved_rdf is None:
                    self.__saved_rdf = RadialDistributionFunction(
                        snapshot1,
                        rng,
                        self.parameters.distance_metric_snapshots_rdf_bins,
                    ).get_rdf()
                rdf1 = self.__saved_rdf
            else:
                rdf1 = RadialDistributionFunction(
                    snapshot1,
                    rng,
                    self.parameters.distance_metric_snapshots_rdf_bins,
                ).get_rdf()
            rdf2 = RadialDistributionFunction(
                snapshot2,
                rng,
                self.parameters.distance_metric_snapshots_rdf_bins,
            ).get_rdf()

            if reduction == "minimal_distance":
                raise Exception(
                    "Combination of distance metric and reduction "
                    "not supported."
                )

            elif reduction == "cosine_distance":
                result = distance.cosine(rdf1, rdf2)

            else:
                raise Exception("Unknown distance metric reduction.")
        else:
            raise Exception("Unknown distance metric selected.")

        return result

    def __denoise(self, signal):
        denoised_signal = np.convolve(
            signal,
            np.ones(self.parameters.distance_metrics_denoising_width)
            / self.parameters.distance_metrics_denoising_width,
            mode="same",
        )
        return denoised_signal

    def analyze_trajectory(
        self, trajectory, equilibrated_snapshot=None, distance_threshold=None
    ):
        """
        Calculate distance metrics/first equilibrated timestep on a trajectory.

        For this step, the RDF+Cosine distance will be used as a distance
        metric. Only the first snapshot is return, all the other quantities
        can be accessed as member variables of the object calling this
        function.

        Parameters
        ----------
        trajectory : ase.io.Trajectory
            Trajectory to be analyzed.

        equilibrated_snapshot : ase.Atoms
            An equilibrated snapshot. Will usually be read from the trajectory
            itself, but may be provided by the user if desired.

        distance_threshold : float
            Distance threshold to be used. Usually determined by analyzing
            the end of the trajectory, but may be provided by user if
            e.g. multiple trajectories are to be compared.

        Returns
        -------
        first_snapshot : int
            First snapshot for which the trajectory is equilibrated.
        """
        # First, we ned to calculate the reduced metrics for the trajectory.
        # For this, we calculate the distance between all the snapshots
        # and the last one.
        self.distance_metrics = []
        if equilibrated_snapshot is None:
            equilibrated_snapshot = trajectory[-1]
        for idx, step in enumerate(trajectory):
            self.distance_metrics.append(
                self._calculate_distance_between_snapshots(
                    equilibrated_snapshot,
                    step,
                    "rdf",
                    "cosine_distance",
                    save_rdf1=True,
                )
            )

        # Now, we denoise the distance metrics.
        self.distance_metrics_denoised = self.__denoise(self.distance_metrics)

        # Which snapshots are considered depends on how we denoise the
        # distance metrics.
        self.first_considered_snapshot = (
            self.parameters.distance_metrics_denoising_width
        )
        self.last_considered_snapshot = (
            np.shape(self.distance_metrics_denoised)[0]
            - self.parameters.distance_metrics_denoising_width
        )
        considered_length = (
            self.last_considered_snapshot - self.first_considered_snapshot
        )

        # Next, the average of the presumed equilibrated part is calculated,
        # and then the first N number of times teps which are below this
        # average is calculated.
        self.average_distance_equilibrated = distance_threshold
        if self.average_distance_equilibrated is None:
            self.average_distance_equilibrated = np.mean(
                self.distance_metrics_denoised[
                    considered_length
                    - int(
                        self.parameters.distance_metrics_estimated_equilibrium
                        * considered_length
                    ) : self.last_considered_snapshot
                ]
            )
        is_below = True
        counter = 0
        first_snapshot = None
        for idx, dist in enumerate(self.distance_metrics_denoised):
            if (
                idx >= self.first_considered_snapshot
                and idx <= self.last_considered_snapshot
            ):
                if is_below:
                    counter += 1
                if dist < self.average_distance_equilibrated:
                    is_below = True
                if dist >= self.average_distance_equilibrated:
                    counter = 0
                    is_below = False
                if (
                    counter
                    == self.parameters.distance_metrics_below_average_counter
                ):
                    first_snapshot = idx
                    break

        print("First equilibrated timestep of trajectory is", first_snapshot)
        return first_snapshot

    def analyze_distance_metric(self, trajectory):
        """
        Calculate the cutoff for the distance metric.

        The distance metric used here is realspace (i.e. the smallest
        displacement of an atom between two snapshots). The cutoff gives
        a lower estimate for the oscillations of the trajectory. Any distance
        above this cutoff can be attributed to the oscillations in the
        trajectory. Any cutoff below is the consquence of temporal
        neighborhood of these snapshots.

        Parameters
        ----------
        trajectory : ase.io.Trajectory
            Trajectory to be analyzed.

        Returns
        -------
        cutoff : float
            Cutoff below which two snapshots can be assumed to be similar
            to each other to a degree that suggests temporal neighborhood.

        """
        # distance metric usef for the snapshot parsing (realspace similarity
        # of the snapshot), we first find the center of the equilibrated part
        # of the trajectory and calculate the differences w.r.t to to it.
        center = (
            int(
                (
                    np.shape(self.distance_metrics_denoised)[0]
                    - self.first_snapshot
                )
                / 2
            )
            + self.first_snapshot
        )
        width = int(
            self.parameters.distance_metrics_estimated_equilibrium
            * np.shape(self.distance_metrics_denoised)[0]
        )
        self.distances_realspace = []
        self.__saved_rdf = None
        for i in range(center - width, center + width):
            self.distances_realspace.append(
                self._calculate_distance_between_snapshots(
                    trajectory[center],
                    trajectory[i],
                    "realspace",
                    "minimal_distance",
                    save_rdf1=True,
                )
            )

        # From these metrics, we assume mean - 2.576 std as limit.
        # This translates to a confidence interval of ~99%, which should
        # make any coincidental similarites unlikely.
        cutoff = np.mean(self.distances_realspace) - 2.576 * np.std(
            self.distances_realspace
        )
        print("Distance metric cutoff is", cutoff)
        return cutoff
