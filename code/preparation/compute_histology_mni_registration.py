import glob
import os
import tempfile

import ants
import numpy as np
import pandas as pd

from dwi_processing.utils import run_command


def compute_histology_mni_registration():
    """
    Updates histology ROI information with MNI coordinates (both physical
    and voxel coordinates).

    Registration uses a two-dimensional tissue mask per histology donor and a
    two-dimensional fractional anisotropy (FA) template in MNI space (extracted
    from the FSL HCP-1065 FA atlas, thresholded at FA >= 0.3)

    Returns
    -------
    None
    """

    # load data
    arv_data_path = os.getenv("MRV_DATA_PATH")
    histo_raw_data_path = os.path.join(arv_data_path, "histology/rawdata")
    cc_atlas_mask_file = os.path.join(
        arv_data_path,
        "mri_in_vivo/templates/masks",
        "space-MNI152NLin6Asym_desc-FAThresholded2D_mask.nii.gz")
    cc_atlas_mask_ants = ants.image_read(cc_atlas_mask_file)
    t1_atlas_image_file = os.path.join(
        arv_data_path,
        "mri_in_vivo/templates/anat",
        "space-MNI152NLin6Asym_T1w.nii.gz")
    t1_atlas_image_ants = ants.image_read(t1_atlas_image_file)
    subject_paths = glob.glob(os.path.join(histo_raw_data_path, "sub-*"))
    mni_roi_coordinates_file = os.path.realpath(
        f"{histo_raw_data_path}/desc-noMNI_roiinfo.tsv")
    roi_info_df = pd.read_csv(
        mni_roi_coordinates_file, sep="\t")

    # allocate array for MNI voxel coordinates
    coords_voxel_3d = np.zeros((roi_info_df.shape[0], 3), dtype=float)

    # fix x coordinate to mid-sagittal slice
    coords_voxel_3d[:, 0] = 90

    for subject_path in subject_paths:
        subject_str = os.path.basename(subject_path)
        subject_id = subject_str.split("-")[1]  # Get the second part

        with tempfile.TemporaryDirectory() as temp_dir:

            # estimate registration: histology <-> MNI space
            run_command(
                f"antsRegistrationSyN.sh "
                f"-y 1 "
                f"-d 2 "
                f"-t s "
                f"-f {cc_atlas_mask_file} "
                f"-m {subject_path}/{subject_str}_label-cc_mask.nii.gz "
                f"-o {temp_dir}/ants "
                f"-i {subject_path}/{subject_str}"
                f"_from-native_to-mni_desc-initial_affine.mat"
            )

            # transform histology tissue mask to MNI space
            run_command(
                f"antsApplyTransforms "
                f"-d 2 "
                f"-r {cc_atlas_mask_file} "
                f"-n linear "
                f"-t {temp_dir}/ants1Warp.nii.gz "
                f"-t {temp_dir}/ants0GenericAffine.mat "
                f"-i {subject_path}/{subject_str}_label-cc_mask.nii.gz "
                f"-o {subject_path}/{subject_str}"
                f"_label-cc_space-mni_mask.nii.gz"
            )

            # filter data for the current subject
            subject_coords_df = roi_info_df[
                roi_info_df["subject_id"] == subject_id]
            mni_roi_coordinates_df = subject_coords_df[
                ["histo_voxel_x", "histo_voxel_y"]].rename(
                columns={"histo_voxel_x": "x", "histo_voxel_y": "y"}
            )

            # transform ROI centroids from MNI to native space
            coords_phyiscal = ants.apply_transforms_to_points(
                dim=2,
                points=mni_roi_coordinates_df,
                transformlist=[
                    f"{temp_dir}/ants0GenericAffine.mat",
                    f"{temp_dir}/ants1InverseWarp.nii.gz"],
                whichtoinvert=[True, False]
            )

            # compute MNI voxel coordinates
            for index, row in coords_phyiscal.iterrows():
                coords_voxel_3d[index, 1:3] = \
                    ants.transform_physical_point_to_index(
                    cc_atlas_mask_ants,
                    row.tolist())

        # compute physical coordinates in MNI space
        coords_physical_3d = ((coords_voxel_3d
                               @ t1_atlas_image_ants.direction.T)
                              + t1_atlas_image_ants.origin)
        coords_physical_3d[:, 0:2] = coords_physical_3d[:, 0:2]

        # add MNI coordinates to ROI info table
        roi_info_df["mni_voxel_x"] = coords_voxel_3d[:, 0]
        roi_info_df["mni_voxel_y"] = coords_voxel_3d[:, 1]
        roi_info_df["mni_voxel_z"] = coords_voxel_3d[:, 2]
        roi_info_df["mni_physical_x"] = coords_physical_3d[:, 0]
        roi_info_df["mni_physical_y"] = coords_physical_3d[:, 1]
        roi_info_df["mni_physical_z"] = coords_physical_3d[:, 2]

        # write updated ROI info table to .tsv file
        roi_info_df.to_csv(os.path.join(f"{histo_raw_data_path}/roiinfo.tsv"),
                             sep="\t", index=False)


if __name__ == "__main__":
    compute_histology_mni_registration()
