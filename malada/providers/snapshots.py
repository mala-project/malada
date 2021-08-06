from .provider import Provider
import os
from shutil import copyfile
import ase.io
import ase
import ase.io
import numpy as np
from scipy.spatial.distance import cdist


class SnapshotsProvider(Provider):
    """Performs a DFT-MD calculation and provides an ASE trjactory.."""
    def __init__(self, parameters, external_snapshots=None):
        super(SnapshotsProvider, self).__init__(parameters)
        self.external_snapshots = external_snapshots
        self.snapshot_file = None

    def provide(self, provider_path, trajectoryfile, temperaturefile):
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


    def __get_first_snapshot(self, trajectoryfile):
        if self.parameters.snapshot_parsing_beginning < 0:
            raise Exception(
                "Automatic detection of first equilibrated snapshot"
                " not supported")
        else:
            return self.parameters.snapshot_parsing_beginning

    def __determine_distance_metric(self, trajectoryfile):
        if self.parameters.distance_metric_snapshots_cutoff < 0:
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
            atoms_to_write = self.__enforce_pbc(trajectory[md_iteration[i]])
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


    def __calculate_distance_between_snapshots(self, snapshot1, snapshot2, dist_type="simple_geometric"):
        positions1 = snapshot1.get_positions()
        positions2 = snapshot2.get_positions()
        if self.parameters.distance_metric_snapshots=="realspace":
            result = np.amin(cdist(positions1, positions2), axis=0)
            result = np.mean(result)
        return result

    @staticmethod
    def __enforce_pbc(atoms):
        """
        Explictly enforeces the PBC on an ASE atoms object.

        QE (and potentially other codes?) do that internally. Meaning that the
        raw positions of atoms (in Angstrom) can lie outside of the unit cell.
        When setting up the DFT calculation, these atoms get shifted into
        the unit cell. Since we directly use these raw positions for the
        descriptor calculation, we need to enforce that in the ASE atoms
        objects, the atoms are explicitly in the unit cell.

        Parameters
        ----------
        atoms : ase.atoms
            The ASE atoms object for which the PBC need to be enforced.

        Returns
        -------
        new_atoms : ase.atoms
            The ASE atoms object for which the PBC have been enforced.
        """
        new_atoms = atoms.copy()
        new_atoms.set_scaled_positions(new_atoms.get_scaled_positions())

        # This might be unecessary, but I think it is nice to have some sort of
        # metric here.
        rescaled_atoms = 0
        for i in range(0, len(atoms)):
            if False in (np.isclose(new_atoms[i].position,
                          atoms[i].position, atol=0.001)):
                rescaled_atoms += 1
        return new_atoms


