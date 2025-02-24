import os
import subprocess


def run_command(command):
    """Run a shell command and handle errors."""
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")


def main():
    base_dir = "/media/laurin/T7/axon_radius_mapping_v2/data/histology/processed"
    for i in range(1, 3):  # Loop from 1 to 2
        formatted_index = f"sub-{i:02d}"
        dir_path = os.path.join(base_dir, formatted_index)

        if os.path.isdir(dir_path):
            print(f"Processing directory {i}: {dir_path}")

            run_command(
                f"antsRegistrationSyN.sh -d 2 -t s -f cc_mask.nii.gz "
                f"-m {dir_path}/images/{formatted_index}_mask.nii.gz "
                f"-o {dir_path}/images/{formatted_index}_ants "
                f"-i {dir_path}/images/{formatted_index}_ants0GenericAffineOld.mat"
            )

            run_command(
                f"antsApplyTransforms -d 2 -r probmap_2d.nii.gz -n linear "
                f"-t {dir_path}/images/{formatted_index}_ants1Warp.nii.gz "
                f"-t {dir_path}/images/{formatted_index}_ants0GenericAffine.mat "
                f"-i {dir_path}/images/{formatted_index}_mask.nii.gz "
                f"-o {dir_path}/images/{formatted_index}_maskAligned.nii.gz"
            )

            run_command(
                f"antsApplyTransformsToPoints -d 2 "
                f"-t [{dir_path}/images/{formatted_index}_ants0GenericAffine.mat,1] "
                f"-t {dir_path}/images/{formatted_index}_ants1InverseWarp.nii.gz "
                f"-i {dir_path}/images/{formatted_index}_roiLocations.csv "
                f"-o {dir_path}/images/{formatted_index}_roiLocationsAligned.csv"
            )

        else:
            print(f"Directory {dir_path} does not exist.")


if __name__ == "__main__":
    main()