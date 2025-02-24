import argparse
import glob
import json
import os
import re

import matlab.engine
import pandas as pd

from dwi_processing.utils import run_command, strip_extension


def get_shell_images_and_bvals(subject_path, sample):
    """
    Retrieve all shell images and corresponding b-values for a given sample in
    a subject's directory.

    Parameters:
    subject_path (str): Path to the subject's directory containing DWI images.
    sample (str): Sample identifier (e.g., "sample-01").

    Returns:
    tuple: A list of file paths corresponding to images per b-value and a list
    of extracted b-values.
    """
    image_per_bval_files = glob.glob(f"{subject_path}/*{sample}_bval-*.nii.gz")

    if not image_per_bval_files:
        raise ValueError("No files matching the pattern were found.")

    bvals = []
    for file in image_per_bval_files:
        match = re.search(r"_bval-(\d+)", file)
        if not match:
            raise ValueError(
                f"Filename does not contain a valid b-value: {file}")
        bvals.append(int(match.group(1)))

    return image_per_bval_files, bvals


def get_unique_samples(data_path):
    """
    Identify unique samples within a dataset based on filenames matching the
    'sample-*' pattern.

    Parameters:
    data_path (str): Path to the directory containing MRI data.

    Returns:
    list: Sorted list of unique sample identifiers found in the dataset.
    """
    unique_samples = set()
    for filename in os.listdir(data_path):
        if filename.endswith(".nii.gz"):
            match = re.search(r"sample-\d+", filename)
            if match:
                unique_samples.add(match.group(0))
    return list(sorted(unique_samples))


def process_ex_vivo_dataset(output_base_dir, subject_id=None, sample_id=None):
    """
    Process ex-vivo MRI dataset, including, preprocessing and estimation of
    powder average signal, fiber orientation, dot fraction, and effective
    radius.

    Parameters:
    output_dir (str): Path to the output directory.
    """

    # set input directory/file paths
    arv_data_path = os.getenv("MRV_DATA_PATH")
    if arv_data_path is None:
        raise ValueError("MRV_DATA_PATH environment variable is not set.")
    raw_data_path = os.path.join(arv_data_path,
                                 "mri_ex_vivo",
                                 "rawdata")
    masks_path = os.path.join(arv_data_path,
                                 "mri_ex_vivo",
                                 "masks")
    processed_data_path = os.path.realpath(output_base_dir)
    mr_protocol_params_file = os.path.join(
        "../parameters/mr_protocol_ex_vivo_experimental.json")
    tissue_params_file = os.path.join(
        "../parameters/tissue_params_ex_vivo.json")

    # tissue and sequence parameters
    with open(mr_protocol_params_file, 'r') as file:
        mr_protocol_params = json.load(file)
    diffusion_gradient_separation = (
        mr_protocol_params)["diffusion_gradient_separation"]
    diffusion_gradient_duration = (
        mr_protocol_params)["diffusion_gradient_duration"]
    if len(set(diffusion_gradient_separation)) > 1:
        raise ValueError(
            "Diffusion gradient separation contains multiple differing values.")
    else:
        diffusion_gradient_separation = diffusion_gradient_separation[0]
    if len(set(diffusion_gradient_duration)) > 1:
        raise ValueError(
            "Diffusion gradient duration contains multiple differing values.")
    else:
        diffusion_gradient_duration = diffusion_gradient_duration[0]

    with open(tissue_params_file, 'r') as file:
        tissue_params = json.load(file)
    diffusivity_axoplasm = tissue_params["d_0"]

    # define b-values for noddi, effective radius and dot fraction estimation
    b_min_noddi = 500
    b_max_noddi = 10000
    b_min_r_eff = 20000
    b_max_r_eff = 100000
    b_min_dot_fraction = 90000

    # define maximum angle with mean fiber direction for dot fraction estimation
    angle_max_dot_fraction = 20

    # Dynamically find all subject dwi folders (sub-*/dwi)
    subject_paths = glob.glob(os.path.join(raw_data_path, "sub-*/dwi"))

    # Filter subjects if specified
    if subject_id:
        subject_paths = [p for p in subject_paths if f"sub-{subject_id}" in p]
        if not subject_paths:
            raise ValueError(f"No data found for subject {subject_id}.")

    # prepare matlab engine
    matlab_engine = matlab.engine.start_matlab()
    matlab_engine.addpath(matlab_engine.genpath(os.getcwd()))

    # Loop through subjects
    for subject_path in subject_paths:
        subject = os.path.basename(os.path.dirname(subject_path))
        output_dir = os.path.join(processed_data_path, subject, "dwi")
        temp_dir = os.path.join(processed_data_path, subject, "temp")
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(temp_dir, exist_ok=True)

        # Get samples
        samples = get_unique_samples(subject_path)
        if sample_id:
            samples = [s for s in samples if s == f"sample-{sample_id}"]
            if not samples:
                raise ValueError(f"No data found for sample {sample_id} "
                                 f"in subject {subject_id}.")

        # Loop through samples
        for sample in samples:
            image_per_bval_files, bvals = get_shell_images_and_bvals(
                os.path.join(raw_data_path, subject, "dwi"), sample)
            tissue_mask_file = os.path.join(
                masks_path, subject, "dwi",
                f"{subject}_{sample}_label-tissue_mask.nii.gz")

            # Preprocess images per shell
            preprocessed_dwi_per_bval_file = []
            powder_average_per_bval_files = []
            noisemap_per_bval_files = []
            for bval, image_per_bval_file in zip(bvals, image_per_bval_files):

                in_file_prefix = strip_extension(image_per_bval_file)
                temp_file_prefix = os.path.join(
                    temp_dir,
                    f"{subject}_{sample}_bval-{bval}")

                # convert to .mif
                run_command(
                    f"mrconvert -force "
                    f"{image_per_bval_file} "
                    f"{temp_file_prefix}_dwi.mif "
                    f"-fslgrad "
                    f"{in_file_prefix}.bvec "
                    f"{in_file_prefix}.bval")

                # estimate noise level
                run_command(
                    f"dwidenoise -force "
                    f"{temp_file_prefix}_dwi.mif "
                    f"{temp_file_prefix}_desc-denoised_dwi.mif "
                    f"-noise {temp_file_prefix}_noisemap.mif")

                # correct for Gibbs ringing
                run_command(
                    f"mrdegibbs -force "
                    f"{temp_file_prefix}_dwi.mif "
                    f"{temp_file_prefix}_desc-degibbsed_dwi.mif")

                # compute mean b=0 image
                run_command(
                    f"dwiextract -bzero "
                    f"{temp_file_prefix}_desc-degibbsed_dwi.mif - "
                    f"-fslgrad {in_file_prefix}.bvec {in_file_prefix}.bval "
                    f"| mrmath -force "
                    f"- "
                    f"mean -axis 3 "
                    f"{temp_file_prefix}_desc-degibbsedMeanB0_dwi.mif")

                # normalize dwi to b=0
                run_command(
                    f"mrcalc -force "
                    f"{temp_file_prefix}_desc-degibbsed_dwi.mif "
                    f"{temp_file_prefix}_desc-degibbsedMeanB0_dwi.mif "
                    f"-divide {temp_file_prefix}_desc-preprocessed_dwi.mif")

                # normalize noisemap to b=0
                run_command(
                    f"mrcalc -force "
                    f"{temp_file_prefix}_noisemap.mif "
                    f"{temp_file_prefix}_desc-degibbsedMeanB0_dwi.mif "
                    f"-divide {temp_file_prefix}_noisemap.nii.gz")

                # convert preprocessed dwi no nifti
                run_command(
                    f"mrconvert -force "
                    f"{temp_file_prefix}_desc-preprocessed_dwi.mif "
                    f"{temp_file_prefix}_desc-preprocessed_dwi.nii.gz "
                    f"-export_grad_fsl "
                    f"{temp_file_prefix}_desc-preprocessed_dwi.bvec "
                    f"{temp_file_prefix}_desc-preprocessed_dwi.bval")

                # compute powder average
                matlab_engine.compute_powder_average_dwi(
                    f"{temp_file_prefix}_desc-preprocessed_dwi.nii.gz",
                    f"{temp_file_prefix}_desc-preprocessed_dwi.bvec",
                    f"{temp_file_prefix}_desc-preprocessed_dwi.bval",
                    f"{temp_file_prefix}_powderAverage.nii.gz",
                    "noisemap_file", f"{temp_file_prefix}_noisemap.nii.gz",
                    nargout=0)

                # store per-shell images
                preprocessed_dwi_per_bval_file.append(
                    f"{temp_file_prefix}_desc-preprocessed_dwi.mif")
                noisemap_per_bval_files.append(
                    f"{temp_file_prefix}_noisemap.nii.gz")
                powder_average_per_bval_files.append(
                    f"{temp_file_prefix}_powderAverage.nii.gz")

            out_file_prefix = os.path.join(output_dir,
                                           f"{subject}_{sample}")
            temp_file_prefix = os.path.join(temp_dir,
                                            f"{subject}_{sample}")

            # merge dwi into single file
            run_command(
                f"mrcat -force "
                f"{' '.join(preprocessed_dwi_per_bval_file)} "
                f"{temp_file_prefix}_desc-preprocessed_dwi.mif")
            run_command(
                f"mrconvert -force "
                f"{temp_file_prefix}_desc-preprocessed_dwi.mif "
                f"{out_file_prefix}_dwi.nii.gz "
                f"-export_grad_fsl {out_file_prefix}_dwi.bvec "
                f"{out_file_prefix}_dwi.bval")

            # merge powder averages into single file
            run_command(
                    f"mrcat -force {' '.join(powder_average_per_bval_files)} "
                    f"{out_file_prefix}_powderAverage.nii.gz")
            pd.DataFrame(bvals).transpose().to_csv(
                  f"{out_file_prefix}_powderAverage.bval",
                sep=" ", index=False, header=None)

            # merge noisemaps into single file
            run_command(
                f"mrcat -force "
                f"{' '.join(noisemap_per_bval_files)} "
                f"{out_file_prefix}_noisemap.nii.gz")
            pd.DataFrame(bvals).transpose().to_csv(
                f"{out_file_prefix}_noisemap.bval",
                index=False, sep=" ", header=None)

            # estimate main fibre direction using noddi
            b_noddi = [0] + [b for b in bvals if b_min_noddi <= b <= b_max_noddi]
            b_min = min(bvals)
            index_b_min = bvals.index(b_min)
            run_command(
                f"fslroi "
                f"{out_file_prefix}_noisemap.nii.gz "
                f"{temp_file_prefix}_noddiNoisemap.nii.gz "
                f"{index_b_min} "
                f"1")
            run_command(
                f"dwiextract -force "
                f"{temp_file_prefix}_desc-preprocessed_dwi.mif "
                f"{temp_file_prefix}_noddiShells.nii.gz "
                f"-shells {','.join(map(str, b_noddi))} "
                f"-export_grad_fsl {temp_file_prefix}_noddiShells.bvec "
                f"{temp_file_prefix}_noddiShells.bval")
            noddi_file_prefix = os.path.join(temp_dir,
                                             f"{subject}_{sample}_noddi")
            matlab_engine.addpath(matlab_engine.genpath(
                "external/NODDI"),
                nargout=0)
            matlab_engine.addpath(matlab_engine.genpath(
                "external/nifti_matlab"),
                nargout=0)
            matlab_engine.compute_noddi_ex_vivo_dwi(
                f"{temp_file_prefix}_noddiShells.nii.gz",
                f"{temp_file_prefix}_noddiShells.bvec",
                f"{temp_file_prefix}_noddiShells.bval",
                noddi_file_prefix,
                f"{temp_file_prefix}_noddiNoisemap.nii.gz",
                1e-9 * diffusivity_axoplasm,
                "mask_file", tissue_mask_file,
                nargout=0)
            run_command(f"mrcat -force "
                        f"{temp_file_prefix}_noddi_fibredirs_xvec.nii "
                        f"{temp_file_prefix}_noddi_fibredirs_yvec.nii "
                        f"{temp_file_prefix}_noddi_fibredirs_zvec.nii "
                        f"{out_file_prefix}_fibreDirs.nii.gz")

            # estimate dot fraction
            b_max = max(bvals)
            index_b_max = bvals.index(b_max)
            run_command(
                f"fslroi "
                f"{out_file_prefix}_noisemap.nii.gz "
                f"{temp_file_prefix}_desc-dotFractionShell_noisemap.nii.gz "
                f"{index_b_max} "
                f"1")
            run_command(
                f"dwiextract -force"
                f" {temp_file_prefix}_desc-preprocessed_dwi.mif "
                f"{temp_file_prefix}_desc-dotFractionShell_dwi.nii.gz "
                f"-shells {b_max} "
                f"-export_grad_fsl"
                f" {temp_file_prefix}_desc-dotFractionShell_dwi.bvec "
                f"{temp_file_prefix}_desc-dotFractionShell_dwi.bval")
            matlab_engine.addpath(matlab_engine.genpath("external/NODDI"),
                                  nargout=0)
            matlab_engine.compute_dot_fraction_dwi(
                f"{temp_file_prefix}_desc-dotFractionShell_dwi.nii.gz",
                f"{temp_file_prefix}_desc-dotFractionShell_dwi.bvec",
                f"{temp_file_prefix}_desc-dotFractionShell_dwi.bval",
                f"{out_file_prefix}_fibreDirs.nii.gz",
                f"{out_file_prefix}_dotFraction.nii.gz",
                f"{temp_file_prefix}_desc-dotFractionShell_noisemap.nii.gz",
                angle_max_dot_fraction,
                b_min_dot_fraction,
                "mask_file", tissue_mask_file,
                nargout=0)

            # estimate dot-free powder average via subtraction
            run_command(
                f"mrcalc "
                f"{out_file_prefix}_powderAverage.nii.gz "
                f"{out_file_prefix}_dotFraction.nii.gz "
                f"-subtract "
                f"- "
                f"| "
                f"mrcalc -force "
                f"- 0 -gt - 0 -if "
                f"{out_file_prefix}_desc-dotFree_powderAverage.nii.gz")

            # estimate effective radius from dot-free powder average
            matlab_engine.compute_effective_radius_multi_shell_dwi(
                f"{out_file_prefix}_desc-dotFree_powderAverage.nii.gz",
                f"{out_file_prefix}_powderAverage.bval",
                f"{out_file_prefix}_effectiveRadius.nii.gz",
                diffusion_gradient_separation,
                diffusion_gradient_duration,
                diffusivity_axoplasm,
                "Neuman",
                "mask_file", tissue_mask_file,
                "min_bval", b_min_r_eff,
                "max_bval", b_max_r_eff,
                nargout=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process ex-vivo MRI dataset and compute effective radius.")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Path to the output directory.")
    parser.add_argument("--subject", type=str, default=None,
                        help="Optional: Process only a specific subject "
                             "(e.g., 'ev01').")
    parser.add_argument("--sample", type=str, default=None,
                        help="Optional: Process only a specific sample "
                             "(e.g., '01').")
    args = parser.parse_args()

    process_ex_vivo_dataset(args.output_dir, args.subject, args.sample)
