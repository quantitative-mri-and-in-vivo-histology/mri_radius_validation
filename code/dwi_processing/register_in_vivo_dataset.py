import os
import subprocess
import glob
import pandas as pd
import shutil
import ants
import antspynet
import json
import numpy as np
import argparse
from dwi_processing.utils import run_command, smooth_masked


def register_in_vivo_dataset(output_dir, subject=None):
    """
    Registers in-vivo diffusion MRI data to OMM and histology space.

    Images from native in-vivo dMRI space are mapped to OMM space. Histological
    ROI locations are mapped to native in-vivo dMRI space, relying on OMM space
    as an intermediate spatial reference.

    Parameters
    ----------
    output_dir : str
        Path to the output directory where processed data will be stored.
    subject : str, optional
        If specified, only this subject (e.g., 'sub-01') will be processed.

    Returns
    -------
    None
    """
    # set input directory paths
    arv_data_path = os.getenv("MRV_DATA_PATH")
    if arv_data_path is None:
        raise ValueError("MRV_DATA_PATH environment variable is not set.")
    mri_data_path = os.path.realpath(
        os.path.join(arv_data_path, "mri_in_vivo"))
    mri_raw_data_path = os.path.join(mri_data_path, "rawdata")
    processed_data_path = os.path.realpath(output_dir)
    histology_data_path = os.path.realpath(
        os.path.join(arv_data_path, "histology"))

    # Dynamically find subject folders
    if args.subject is not None:
        subjects = [f"sub-{subject}"]
    else:
        subjects = glob.glob(os.path.join(mri_raw_data_path, "sub-*"))

    # define OMM template files
    t1_omm_template_file = os.path.join(
        mri_data_path, "templates/anat",
        "space-omm_T1w.nii.gz")
    t1_omm_template_brain_mask_file = os.path.join(
        mri_data_path, "templates/masks",
        "space-omm_label-brain_mask.nii.gz")
    fa_omm_template_file = os.path.join(
        mri_data_path, "templates/dwi",
        "space-omm_FA.nii.gz")
    corpus_callosum_mask_file = os.path.join(
        mri_data_path, "templates/masks",
        "space-omm_label-ccMid_mask.nii.gz")
    bids_omm_id = "omm"

    # load threshold definitions
    cc_eval_params_file = os.path.join(
        "../parameters/cc_evaluation_parameters.json")
    with open(cc_eval_params_file, 'r') as file:
        cc_eval_params = json.load(file)
    fa_min_thres = cc_eval_params["fa_min"]
    r_eff_min_thres = cc_eval_params["r_eff_min"]

    # Loop through samples
    for subject_path in subjects:
        subject = os.path.basename(subject_path)
        root_output_dir = os.path.join(processed_data_path, subject)
        dwi_output_dir = os.path.join(processed_data_path, subject, "dwi")
        temp_dir = os.path.join(processed_data_path, subject, "temp")
        dwi_out_file_prefix = os.path.join(dwi_output_dir, f"{subject}")
        root_out_file_prefix = os.path.join(root_output_dir, f"{subject}")
        temp_file_prefix = os.path.join(temp_dir, f"{subject}")

        # create output directory if it doesn't exist
        os.makedirs(temp_dir, exist_ok=True)

        # load T1w images
        t1_image_file = os.path.join(processed_data_path, f"{subject}", "anat",
                                     f"{subject}_T1w.nii.gz")
        t1_mask_file = os.path.join(processed_data_path, f"{subject}", "anat",
                                    f"{subject}_label-brain_mask.nii.gz")

        # mask FA image for registration
        run_command(
            f"mrcalc -force {dwi_out_file_prefix}_FA.nii.gz "
            f"{dwi_out_file_prefix}_label-brain_mask.nii.gz "
            f"-mul {temp_file_prefix}_desc-brainMasked_FA.nii.gz"
        )

        # register to OMM
        run_command(
            f"antsRegistration --verbose 1 "
            f"--random-seed 0 "
            "--dimensionality 3 --float 0 --collapse-output-transforms 1 "
            f"--output {temp_file_prefix}_from-native_to-{bids_omm_id}_ants "
            "--interpolation BSpline[3] --use-histogram-matching 0 "
            "--winsorize-image-intensities [ 0.005,0.995 ] "
            # Rigid stage
            "--transform Rigid[ 0.1 ] "
            f"-x [ NULL,NULL ] "
            f"--metric MI[ "
            f"{t1_omm_template_file}, "
            f"{t1_image_file}, "
            f"0.75, 32, Regular, 0.25 ] "
            f"--metric MI[ "
            f"{fa_omm_template_file}, "
            f"{temp_file_prefix}_desc-brainMasked_FA.nii.gz, "
            f"0.25, 32, Regular, 0.25 ] "
            "--convergence [ 1000x500x250x100,1e-6,10 ] "
            "--shrink-factors 8x4x2x1 "
            "--smoothing-sigmas 3x2x1x0vox "
            # Affine stage
            "--transform Affine[ 0.1 ] "
            f"-x [ NULL,NULL ] "
            f"--metric MI[ "
            f"{t1_omm_template_file}, "
            f"{t1_image_file}, "
            f"0.3, 32, Regular, 0.25 ] "
            f"--metric MI[ "
            f"{fa_omm_template_file}, "
            f"{temp_file_prefix}_desc-brainMasked_FA.nii.gz, "
            f"0.7, 32, Regular, 0.25 ] "
            "--convergence [ 1000x500x250x100,1e-6,10 ] "
            "--shrink-factors 8x4x2x1 "
            "--smoothing-sigmas 3x2x1x0vox "
            # SyN (nonlinear) stage
            "--transform SyN[ 0.1,3,0 ] "
            f"-x [ {t1_omm_template_brain_mask_file}, {t1_mask_file} ] "
            f"--metric CC[ "
            f"{t1_omm_template_file}, "
            f"{t1_image_file}, "
            f"0.2, 4 ] "
            f"--metric CC[ "
            f"{fa_omm_template_file}, "
            f"{temp_file_prefix}_desc-brainMasked_FA.nii.gz, "
            f"0.8, 4 ] "
            "--convergence [ 100x70x50x20,1e-6,10 ] "
            "--shrink-factors 8x4x2x1 "
            "--smoothing-sigmas 3x2x1x0vox"
        )
        shutil.move(
            f"{temp_file_prefix}"
            f"_from-native_to-{bids_omm_id}_ants0GenericAffine.mat",
            f"{root_out_file_prefix}"
            f"_from-native_to-{bids_omm_id}_affine.mat")
        shutil.move(
            f"{temp_file_prefix}"
            f"_from-native_to-{bids_omm_id}_ants1Warp.nii.gz",
            f"{root_out_file_prefix}"
            f"_from-native_to-{bids_omm_id}_warp.nii.gz")
        shutil.move(
            f"{temp_file_prefix}"
            f"_from-native_to-{bids_omm_id}_ants1InverseWarp.nii.gz",
            f"{root_out_file_prefix}"
            f"_from-{bids_omm_id}_to-native_warp.nii.gz")

        # transform FA to OMM space
        fa_template_ants = ants.image_read(fa_omm_template_file)
        fa_image_ants = ants.image_read(f"{dwi_out_file_prefix}_FA.nii.gz")
        fa_image_ants_reg_to_omm = ants.apply_transforms(
            fixed=fa_template_ants,
            moving=fa_image_ants,
            transformlist=[
                f"{root_out_file_prefix}"
                f"_from-native_to-{bids_omm_id}_warp.nii.gz",
                f"{root_out_file_prefix}"
                f"_from-native_to-{bids_omm_id}_affine.mat"],
            whichtoinvert=[False, False],
            interpolator="bSpline")
        ants.image_write(
            fa_image_ants_reg_to_omm,
            f"{dwi_out_file_prefix}_space-{bids_omm_id}_FA.nii.gz")

        # # transform effective radius to OMM space
        effective_radius_image_ants = ants.image_read(
            f"{dwi_out_file_prefix}_effectiveRadius.nii.gz"
        )
        effective_radius_image_ants_reg_to_omm = ants.apply_transforms(
            fixed=fa_template_ants,
            moving=effective_radius_image_ants,
            transformlist=[
                f"{root_out_file_prefix}"
                f"_from-native_to-{bids_omm_id}_warp.nii.gz",
                f"{root_out_file_prefix}"
                f"_from-native_to-{bids_omm_id}_affine.mat"],
            whichtoinvert=[False, False],
            interpolator="bSpline")
        ants.image_write(
            effective_radius_image_ants_reg_to_omm,
            f"{dwi_out_file_prefix}"
            f"_space-{bids_omm_id}_effectiveRadius.nii.gz")

        # transform coarse corpus callosum mask to native space
        cc_mask_image_ants = ants.image_read(corpus_callosum_mask_file)
        cc_mask_image_ants_reg_to_native = ants.apply_transforms(
            fixed=effective_radius_image_ants,
            moving=cc_mask_image_ants,
            transformlist=[
                f"{root_out_file_prefix}"
                f"_from-native_to-{bids_omm_id}_affine.mat",
                f"{root_out_file_prefix}"
                f"_from-{bids_omm_id}_to-native_warp.nii.gz"],
            whichtoinvert=[True, False],
            interpolator="nearestNeighbor")
        ants.image_write(
            cc_mask_image_ants_reg_to_native,
            f"{dwi_out_file_prefix}_label-ccMid_mask.nii.gz")

        # transform histology ROI locations to native space
        omm_roi_coordinates_file = os.path.join(histology_data_path, "rawdata",
                                                "roiinfo.tsv")
        roi_info_df = pd.read_csv(
            omm_roi_coordinates_file, sep="\t")
        omm_roi_coordinates_df = roi_info_df[
            ["omm_physical_x", "omm_physical_y", "omm_physical_z"]].rename(
            columns={"omm_physical_x": "x", "omm_physical_y": "y",
                     "omm_physical_z": "z"}
        )
        coords_phyiscal = ants.apply_transforms_to_points(
            dim=3,
            points=omm_roi_coordinates_df,
            transformlist=[
                f"{root_out_file_prefix}"
                f"_from-native_to-{bids_omm_id}_warp.nii.gz",
                f"{root_out_file_prefix}"
                f"_from-native_to-{bids_omm_id}_affine.mat"],
            whichtoinvert=[False, False]
        )
        coords_voxel = np.zeros(
            (coords_phyiscal.shape[0], coords_phyiscal.shape[1]),
            dtype=float)
        for index, row in coords_phyiscal.iterrows():
            coords_voxel[index, :] = ants.transform_physical_point_to_index(
                effective_radius_image_ants,
                row.tolist())

        # add native coordinates to data frame and save
        output_df = roi_info_df.copy()
        output_df["native_voxel_x"] = np.round(coords_voxel[:, 0]).astype(int)
        output_df["native_voxel_y"] = np.round(coords_voxel[:, 1]).astype(int)
        output_df["native_voxel_z"] = np.round(coords_voxel[:, 2]).astype(int)
        output_df["native_physical_x"] = coords_phyiscal["x"]
        output_df["native_physical_y"] = coords_phyiscal["y"]
        output_df["native_physical_z"] = coords_phyiscal["z"]
        output_df.to_csv(os.path.join(f"{dwi_out_file_prefix}_roiinfo.tsv"),
                         sep="\t", index=False)

        # create one label image per histology subject, showing ROI locations
        # in native space (each ROI identified by unique_id)
        unique_subjects = roi_info_df["subject_id"].unique()
        for histo_subject in unique_subjects:

            # filter data for the current subject
            subject_df = roi_info_df[
                roi_info_df["subject_id"] == histo_subject]

            # create all-zero image with correct dimensions
            label_image = effective_radius_image_ants.clone()
            label_image[:] = 0
            label_image = label_image.astype("uint8")

            # assign unique_id labels to ROIs in the image
            for idx, row in subject_df.iterrows():
                transformed_index = coords_voxel[idx, :]
                transformed_index = np.round(transformed_index)
                label_image[tuple(map(int, transformed_index))] = row[
                    "unique_id"]

            # save the labeled image
            label_image_filename = (f"{dwi_out_file_prefix}"
                                    f"_histosub-{histo_subject}"
                                    f"_label-rois_dseg.nii.gz")
            ants.image_write(label_image, label_image_filename)

        # create spatially smoothed effective radius map using mask-constrained
        # smoothing (see https://doi.org/10.1016/j.neuroimage.2008.09.041)
        run_command(
            "mrcalc -force "
            f"{dwi_out_file_prefix}_FA.nii.gz {fa_min_thres} -abs -gt "
            f"{dwi_out_file_prefix}_label-ccMid_mask.nii.gz -mult "
            f"{dwi_out_file_prefix}_effectiveRadius.nii.gz -finite -mult "
            f"{dwi_out_file_prefix}_effectiveRadius.nii.gz {r_eff_min_thres}"
            f" -abs -gt -mult "
            f"{temp_file_prefix}_label-cc_desc-strict_mask.nii.gz"
        )
        smooth_masked(
            f"{dwi_out_file_prefix}_effectiveRadius.nii.gz",
            f"{dwi_out_file_prefix}"
            f"_roi-cc_desc-smoothed_effectiveRadius.nii.gz",
            f"{temp_file_prefix}_label-cc_desc-strict_mask.nii.gz",
            fwhm=3.75,
            temp_dir=temp_dir
        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Register in-vivo MRI dataset to OMM space.")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Path to the output directory.")
    parser.add_argument("--subject", type=str, default=None,
                        help="Optional: Process a specific subject.")
    args = parser.parse_args()
    register_in_vivo_dataset(args.output_dir, args.subject)
