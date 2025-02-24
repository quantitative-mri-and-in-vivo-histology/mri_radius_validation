import argparse
import glob
import json
import os
import shutil

import ants
import antspynet
import matlab.engine

from dwi_processing.utils import run_command, combine_warps


def process_in_vivo_dataset(output_dir, subject=None):
    """
    Processes in-vivo diffusion MRI data, including
    effective radius estimation.

    The preprocessing pipeline follows a two-pass approach:
    1. A forward run applies all preprocessing steps:
       - Gibbs artifact correction
       - susceptibility, eddy current, and motion correction
       - gradient distortion correction
       - shell alignment
       - registration of all images to anatomical native space
    2. The pipeline then revisits pre-eddy data and applies all transformations
       in a single interpolation step, minimizing data degradation.

    After preprocessing, the script estimates the effective radius.

    Parameters
    ----------
    output_dir : str
        Path to the output directory where processed data will be stored.
    subject : str, optional
        Process only a specific subject (e.g., 'iv01'). If None, all subjects
        are processed.

    Returns
    -------
    None
        The processed dataset is written to `output_dir`.
    """
    # set input directory/file paths
    arv_data_path = os.getenv("MRV_DATA_PATH")
    if arv_data_path is None:
        raise ValueError("MRV_DATA_PATH environment variable is not set.")
    mri_data_path = os.path.realpath(
        os.path.join(arv_data_path, "mri_in_vivo"))
    mri_raw_data_path = os.path.join(mri_data_path, "rawdata")
    processed_data_path = os.path.realpath(output_dir)
    sl_spec_file = os.path.join(mri_data_path, "resources/slspec.txt")
    connectom_coeff_file = os.path.join(
        mri_data_path,
        "resources/connectom_coeff.grad")
    mr_protocol_params_file = os.path.join(
        "../parameters/mr_protocol_in_vivo_experimental.json")
    tissue_params_file = os.path.join(
        "../parameters/tissue_params_in_vivo.json")

    # Dynamically find all subject folders (sub-*)
    if subject is not None:
        subjects = [f"sub-{args.subject}"]
    else:
        subjects = glob.glob(os.path.join(mri_raw_data_path, "sub-*"))

    # parse sequence parameters for effective radius estimation
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

    # parse tissue parameters for effective radius estimation
    with open(tissue_params_file, 'r') as file:
        tissue_params = json.load(file)
    diffusivity_axoplasm = tissue_params["d_0"]

    # prepare matlab engine
    matlab_engine = matlab.engine.start_matlab()
    matlab_engine.addpath(matlab_engine.genpath(os.getcwd()))

    # Loop through subjects
    for subject_path in subjects:
        subject = os.path.basename(subject_path)
        dwi_output_dir = os.path.join(processed_data_path, subject, "dwi")
        anat_output_dir = os.path.join(processed_data_path, subject, "anat")
        temp_dir = os.path.join(processed_data_path, subject, "temp")

        # Create output directory if it doesn't exist
        os.makedirs(dwi_output_dir, exist_ok=True)
        os.makedirs(anat_output_dir, exist_ok=True)
        os.makedirs(temp_dir, exist_ok=True)

        # define input files and output file prefixes
        axon_radius_in_file_prefix = os.path.join(
            mri_raw_data_path,
            f"{subject}",
            "dwi",
            f"{subject}_acq-axonRadius")
        axon_radius_temp_file_prefix = os.path.join(
            temp_dir,
            f"{subject}_acq-axonRadius")
        fast_dki_in_file_prefix = os.path.join(
            mri_raw_data_path,
            f"{subject}",
            "dwi",
            f"{subject}_acq-fastDki")
        fast_dki_temp_file_prefix = os.path.join(
            temp_dir,
            f"{subject}_acq-fastDki")
        t1_image_file = os.path.join(
            mri_raw_data_path,
            f"{subject}",
            "anat",
            f"{subject}_T1w.nii.gz")
        rpe_image_file = os.path.join(
            mri_raw_data_path,
            f"{subject}",
            "fmap",
            f"{subject}_dir-pa_epi.nii.gz")
        in_file_prefixes = [fast_dki_in_file_prefix, axon_radius_in_file_prefix]

        # Run initial preprocessing pass
        for in_file_prefix in in_file_prefixes:

            temp_file_prefix = os.path.join(
                temp_dir,
                os.path.basename(in_file_prefix))

            # get readout time from .json file
            with open(f"{in_file_prefix}_dwi.json", 'r') as file:
                in_file_json = json.load(file)
            readout_time = in_file_json["TotalReadoutTime"]

            # convert to mif for preprocessing
            run_command(
                f"mrconvert -force "
                f"{in_file_prefix}_dwi.nii.gz "
                f"{temp_file_prefix}_dwi.mif "
                f"-fslgrad "
                f"{in_file_prefix}_dwi.bvec "
                f"{in_file_prefix}_dwi.bval")

            # correct gibbs ringing
            run_command(
                f"mrdegibbs -force "
                f"{temp_file_prefix}_dwi.mif "
                f"{temp_file_prefix}_desc-degibbsed_dwi.mif"
            )

            # create reverse phase encoding B0 pair for susceptibility
            # correction
            run_command(
                f"dwiextract -force "
                f"-bzero "
                f"{temp_file_prefix}_dwi.mif "
                f"- "
                f"| "
                f"mrmath -force "
                f"- "
                f"mean "
                f"{temp_file_prefix}_dir-ap_desc-b0_dwi.mif "
                f"-axis 3")
            run_command(
                f"mrmath -force "
                f"{rpe_image_file} "
                f"mean "
                f"{temp_file_prefix}_dir-pa_desc-b0_dwi.mif "
                f"-axis 3")
            run_command(
                f"mrcat -force -axis 3 "
                f"{temp_file_prefix}_dir-ap_desc-b0_dwi.mif "
                f"{temp_file_prefix}_dir-pa_desc-b0_dwi.mif "
                f"{temp_file_prefix}_desc-rpePair_dwi.mif")

            # correct for susceptibility, eddy current and motion artifacts
            if os.path.isdir(f"{temp_file_prefix}_dwifslpreproc"):
                shutil.rmtree(f"{temp_file_prefix}_dwifslpreproc")
            eddy_options = (
                "--repol "
                "--dont_sep_offs_move "
                "--data_is_shelled "
                "--dfields "
                "--flm=cubic "
                "--ol_type=both "
                "--mporder=13 "
                "--cnr_maps "
                "--residuals "
                "--nvoxhp=4000 "
                "--niter=8 "
                "--fwhm=10,4,2,0,0,0,0,0 "
                "--s2v_niter=8 "
            )
            run_command(
                f"dwifslpreproc -force "
                f"{temp_file_prefix}_desc-degibbsed_dwi.mif "
                f"{temp_file_prefix}_desc-dwifslpreproc_dwi.mif "
                f"-nocleanup "
                f"-rpe_pair "
                f"-align_seepi "
                f"-pe_dir ap "
                f"-se_epi {temp_file_prefix}_desc-rpePair_dwi.mif "
                f"-eddy_slspec {sl_spec_file} "
                f"-readout_time {readout_time} "
                f"-scratch {temp_file_prefix}_dwifslpreproc "
                f"-eddyqc_all {temp_file_prefix}_eddyqc "
                f"-eddy_options ' {eddy_options}'"
            )

            # Copy displacement fields and delete dwifslpreproc folder
            dfields_dir = os.path.join(f"{temp_file_prefix}_eddyqc", "dfields")
            os.makedirs(dfields_dir, exist_ok=True)
            dwi_files = glob.glob(os.path.join(
                f"{temp_file_prefix}_dwifslpreproc",
                "dwi_post_eddy.eddy_displacement_fields.*.nii.gz"))
            if os.path.isdir(f"{temp_file_prefix}_dwifslpreproc"):
                for file in dwi_files:
                    shutil.copy(file, dfields_dir)
                shutil.copy(
                    f"{temp_file_prefix}_dwifslpreproc/"
                    f"dwi_post_eddy.nii.gz",
                    f"{temp_file_prefix}_eddyqc")

                shutil.rmtree(f"{temp_file_prefix}_dwifslpreproc")

            # convert back to nifti
            run_command(
                f"mrconvert -force "
                f"{temp_file_prefix}_desc-dwifslpreproc_dwi.mif "
                f"{temp_file_prefix}_desc-dwifslpreproc_dwi.nii.gz "
                f"-export_grad_fsl "
                f"{temp_file_prefix}_desc-preprocessed_dwi.bvec "
                f"{temp_file_prefix}_desc-preprocessed_dwi.bval")

            # estimate gradient distortion correction
            run_command(
                f"fslroi "
                f"{temp_file_prefix}_eddyqc/dwi_post_eddy.nii.gz "
                f"{temp_file_prefix}_eddyqc/dwi_post_eddy_first_vol.nii.gz "
                f"0 1")
            run_command(
                f"gradient_unwarp.py "
                f"--numpoints 128 "
                f"--interp_order 2 "
                f"-n "
                f"{temp_file_prefix}_eddyqc/dwi_post_eddy_first_vol.nii.gz "
                f"{temp_file_prefix}_eddyqc/"
                f"dwi_post_eddy_first_vol_unwarped.nii.gz "
                f"siemens -g {connectom_coeff_file} ")
            shutil.move("fullWarp_abs.nii.gz",
                        f"{temp_file_prefix}_desc-gnlcAbs_warp.nii.gz")
            run_command(
                f"convertwarp --abs --relout "
                f"--ref={temp_file_prefix}_eddyqc/"
                f"dwi_post_eddy_first_vol.nii.gz "
                f"--warp1={temp_file_prefix}_desc-gnlcAbs_warp.nii.gz "
                f"--out={temp_file_prefix}_desc-gnlcRel_warp.nii.gz")

            # apply gradient distortion correction
            run_command(
                f"applywarp --rel --interp=spline "
                f"-i {temp_file_prefix}_eddyqc/dwi_post_eddy.nii.gz "
                f"-o {temp_file_prefix}_desc-preprocessed_dwi.nii.gz "
                f"-r {temp_file_prefix}_eddyqc/dwi_post_eddy_first_vol.nii.gz "
                f"-w {temp_file_prefix}_desc-gnlcRel_warp.nii.gz")

            # convert to mrtrix order
            run_command(
                "mrconvert -force -strides -1,-2,3,4 "
                f"{temp_file_prefix}_desc-preprocessed_dwi.nii.gz "
                f"{temp_file_prefix}_desc-preprocessed_dwi.nii.gz -force"
            )

            # create brain mask
            run_command(
                f"fslroi {temp_file_prefix}_desc-preprocessed_dwi.nii.gz "
                f"{temp_file_prefix}_desc-preprocessedFirstVol_dwi.nii.gz "
                f"0 1")
            run_command(
                f"bet {temp_file_prefix}_desc-preprocessedFirstVol_dwi.nii.gz "
                f"{temp_file_prefix}"
                f"_desc-preprocessedFirstVolMasked_dwi.nii.gz "
                f"-R -m")
            shutil.move(
                f"{temp_file_prefix}"
                f"_desc-preprocessedFirstVolMasked_dwi_mask.nii.gz",
                f"{temp_file_prefix}_label-brain_mask.nii.gz")

            # extract mean b0 for registration of fastDki and axonRadius frames
            run_command(
                f"dwiextract -force "
                f"-bzero "
                f"{temp_file_prefix}_desc-preprocessed_dwi.nii.gz "
                f"- "
                f"-fslgrad "
                f"{temp_file_prefix}_desc-preprocessed_dwi.bvec "
                f"{temp_file_prefix}_desc-preprocessed_dwi.bval "
                f"| "
                f"mrmath -force "
                f"- "
                f"mean "
                f"{temp_file_prefix}_desc-preprocessedMeanB0_dwi.nii.gz "
                f"-axis 3")

            # estimate noisemap
            run_command(
                f"dwidenoise -force {temp_file_prefix}_dwi.mif "
                f"{temp_file_prefix}_desc-denoised_dwi.mif "
                f"-noise {temp_file_prefix}_noisemap.nii.gz")

            # compute spherical average
            matlab_engine.compute_powder_average_dwi(
                f"{temp_file_prefix}_desc-preprocessed_dwi.nii.gz",
                f"{temp_file_prefix}_desc-preprocessed_dwi.bvec",
                f"{temp_file_prefix}_desc-preprocessed_dwi.bval",
                f"{temp_file_prefix}_powderAverage.nii.gz",
                "noisemap_file",
                f"{temp_file_prefix}_noisemap.nii.gz",
                "mask_file",
                f"{temp_file_prefix}_label-brain_mask.nii.gz",
                nargout=0
            )

        temp_file_prefix = os.path.join(temp_dir, f"{subject}")
        dwi_out_file_prefix = os.path.join(dwi_output_dir, f"{subject}")
        anat_out_file_prefix = os.path.join(anat_output_dir, f"{subject}")

        # extract powder averages and first b=0 volumes for registration
        run_command(
            f"fslroi {fast_dki_temp_file_prefix}_powderAverage.nii.gz "
            f"{fast_dki_temp_file_prefix}_desc-b500_powderAverage.nii.gz "
            f"0 1")
        run_command(
            f"fslroi {fast_dki_temp_file_prefix}_powderAverage.nii.gz "
            f"{fast_dki_temp_file_prefix}_desc-b1000_powderAverage.nii.gz "
            f"1 1")
        run_command(
            f"fslroi {fast_dki_temp_file_prefix}_powderAverage.nii.gz "
            f"{fast_dki_temp_file_prefix}_desc-b2500_powderAverage.nii.gz "
            f"2 1")
        run_command(
            f"fslroi {axon_radius_temp_file_prefix}_powderAverage.nii.gz "
            f"{axon_radius_temp_file_prefix}_desc-b6000_powderAverage.nii.gz "
            f"0 1")
        run_command(
            f"fslroi {axon_radius_temp_file_prefix}_powderAverage.nii.gz "
            f"{axon_radius_temp_file_prefix}_desc-b30450_powderAverage.nii.gz "
            f"1 1")
        run_command(
            f"fslroi {axon_radius_temp_file_prefix}"
            f"_desc-preprocessed_dwi.nii.gz "
            f"{axon_radius_temp_file_prefix}"
            f"_desc-preprocessedFirstVol_dwi.nii.gz "
            f"0 1")
        run_command(
            f"fslroi {fast_dki_temp_file_prefix}_desc-preprocessed_dwi.nii.gz "
            f"{fast_dki_temp_file_prefix}_desc-preprocessedFirstVol_dwi.nii.gz "
            f"0 1")

        # transform brain mask for T1w image
        t1_image_ants = ants.image_read(t1_image_file)
        t1_brain_mask_ants = antspynet.brain_extraction(
            t1_image_ants,
            "t1")
        t1_brain_mask_ants = ants.threshold_image(t1_brain_mask_ants,
                                                  low_thresh=0.01)
        ants.image_write(
            t1_brain_mask_ants,
            f"{anat_out_file_prefix}_label-brain_mask.nii.gz"
        )
        t1_image_ants_denoised = ants.denoise_image(
            t1_image_ants)
        t1_image_n4 = ants.n4_bias_field_correction(
            t1_image_ants_denoised)
        ants.image_write(t1_image_n4,
                         f"{anat_out_file_prefix}_T1w.nii.gz")

        # register b=30450 to b=6000 shell using powder averages
        run_command(
            f"antsRegistration --verbose 1 "
            f"--random-seed 0 "
            "--dimensionality 3 --float 0 --collapse-output-transforms 1 "
            f"--output {temp_file_prefix}_from-b30450_to-b6000_ants "
            "--interpolation Linear --use-histogram-matching 0 "
            "--winsorize-image-intensities [ 0.005,0.995 ] "
            # Rigid stage
            "--transform Rigid[ 0.1 ] "
            f"-x [ NULL,NULL ] "
            f"--metric MI[ "
            f"{axon_radius_temp_file_prefix}_desc-b6000_powderAverage.nii.gz, "
            f"{axon_radius_temp_file_prefix}_desc-b30450_powderAverage.nii.gz, "
            f"1.0, 32, Regular, 0.25 ] "
            "--convergence [ 3000x1500x750x375,1e-8,10 ] "
            "--shrink-factors 8x4x2x1 "
            "--smoothing-sigmas 3x2x1x0vox"
        )
        shutil.move(
            f"{temp_file_prefix}_from-b30450_to-b6000_ants0GenericAffine.mat",
            f"{temp_file_prefix}_from-b30450_to-b6000_affine.mat")

        # register fastDki to T1w image
        run_command(
            f"antsRegistration --verbose 1 "
            f"--random-seed 0 "
            "--dimensionality 3 --float 0 --collapse-output-transforms 1 "
            f"--output {temp_file_prefix}_from-fastDki_to-anat_ants "
            "--interpolation Linear --use-histogram-matching 0 "
            "--winsorize-image-intensities [ 0.005,0.995 ] "
            # Rigid stage (unmasked; emphasis on b=0 images with skull)
            "--transform Rigid[ 0.1 ] "
            f"-x [ NULL,NULL ] "
            f"--metric MI[ "
            f"{anat_out_file_prefix}_T1w.nii.gz, "
            f"{fast_dki_temp_file_prefix}"
            f"_desc-preprocessedFirstVol_dwi.nii.gz, "
            f"1.0, 32, Regular, 0.25 ] "
            "--convergence [ 3000x1500x750x375,1e-8,10 ] "
            "--shrink-factors 8x4x2x1 "
            "--smoothing-sigmas 3x2x1x0vox "
            # Rigid stage (masked; weight b=0 and powder average equally)
            "--transform Rigid[ 0.1 ] "
            f"-x [ {anat_out_file_prefix}_label-brain_mask.nii.gz, "
            f"{fast_dki_temp_file_prefix}_label-brain_mask.nii.gz ] "
            f"--metric MI[ "
            f"{anat_out_file_prefix}_T1w.nii.gz, "
            f"{fast_dki_temp_file_prefix}"
            f"_desc-preprocessedFirstVol_dwi.nii.gz, "
            f"0.25, 32, Regular, 0.25 ] "
            f"--metric MI[ {anat_out_file_prefix}_T1w.nii.gz, "
            f"{fast_dki_temp_file_prefix}_desc-b500_powderAverage.nii.gz, "
            f"0.25, 32, Regular, 0.25 ] "
            f"--metric MI[ {anat_out_file_prefix}_T1w.nii.gz, "
            f"{fast_dki_temp_file_prefix}_desc-b1000_powderAverage.nii.gz, "
            f"0.25, 32, Regular, 0.25 ] "
            f"--metric MI[ {anat_out_file_prefix}_T1w.nii.gz, "
            f"{fast_dki_temp_file_prefix}_desc-b2500_powderAverage.nii.gz,"
            f" 0.25, 32, Regular, 0.25 ] "
            "--convergence [ 3000x1500x750x375,1e-8,10 ] "
            "--shrink-factors 8x4x2x1 "
            "--smoothing-sigmas 3x2x1x0vox"
        )
        shutil.move(
            f"{temp_file_prefix}"
            f"_from-fastDki_to-anat_ants0GenericAffine.mat",
            f"{temp_file_prefix}_from-fastDki_to-anat_affine.mat")

        # register b=6000 to T1w image
        run_command(
            f"antsRegistration --verbose 1 "
            f"--random-seed 0 "
            "--dimensionality 3 --float 0 --collapse-output-transforms 1 "
            f"--output {temp_file_prefix}_from-b6000_to-anat_ants "
            "--interpolation Linear --use-histogram-matching 0 "
            "--winsorize-image-intensities [ 0.005,0.995 ] "
            # Rigid stage (unmasked; emphasis on b=0 images with skull)
            "--transform Rigid[ 0.1 ] "
            f"-x [ NULL,NULL ] "
            f"--metric MI[ "
            f"{anat_out_file_prefix}_T1w.nii.gz, "
            f"{axon_radius_temp_file_prefix}"
            f"_desc-preprocessedFirstVol_dwi.nii.gz, "
            f"0.75, 32, Regular, 0.25 ] "
            f"--metric MI["
            f" {anat_out_file_prefix}_T1w.nii.gz, "
            f"{axon_radius_temp_file_prefix}_desc-b6000_powderAverage.nii.gz, "
            f"0.25, 32, Regular, 0.25 ] "
            "--convergence [ 3000x1500x750x375,1e-8,10 ] "
            "--shrink-factors 8x4x2x1 "
            "--smoothing-sigmas 3x2x1x0vox "
            # Rigid stage (masked; weight b=0 and powder average equally)
            "--transform Rigid[ 0.1 ] "
            f"-x [ {anat_out_file_prefix}_label-brain_mask.nii.gz, "
            f"{axon_radius_temp_file_prefix}_label-brain_mask.nii.gz ] "
            f"--metric MI[ "
            f"{anat_out_file_prefix}_T1w.nii.gz, "
            f"{axon_radius_temp_file_prefix}"
            f"_desc-preprocessedFirstVol_dwi.nii.gz, "
            f"0.5, 32, Regular, 0.25 ] "
            f"--metric MI[ "
            f"{anat_out_file_prefix}_T1w.nii.gz, "
            f"{axon_radius_temp_file_prefix}_desc-b6000_powderAverage.nii.gz, "
            f"0.5, 32, Regular, 0.25 ] "
            "--convergence [ 3000x1500x750x375,1e-8,10 ] "
            "--shrink-factors 8x4x2x1 "
            "--smoothing-sigmas 3x2x1x0vox"
        )
        shutil.move(
            f"{temp_file_prefix}_from-b6000_to-anat_ants0GenericAffine.mat",
            f"{temp_file_prefix}_from-b6000_to-anat_affine.mat")

        # compute b=30450 -> T1w transform by combining
        # b=30450 -> b=6000 and b=6000 -> T1w
        run_command(
            "antsApplyTransforms -d 3 "
            f"-r {anat_out_file_prefix}_T1w.nii.gz "
            f"-i {axon_radius_temp_file_prefix}"
            f"_desc-preprocessedFirstVol_dwi.nii.gz "
            f"-o Linear[{temp_file_prefix}_from-b30450_to-anat_affine.mat,0] "
            f"-t {temp_file_prefix}_from-b6000_to-anat_affine.mat "
            f"-t {temp_file_prefix}_from-b30450_to-b6000_affine.mat"
        )

        # prepare rigid transformations
        # (fastDki -> T1w, b=30450 -> T1w and b=6000 -> T1w)
        # for application to eddy output
        run_command(
            "c3d_affine_tool "
            f"-ref {axon_radius_temp_file_prefix}_eddyqc/"
            f"dwi_post_eddy_first_vol.nii.gz "
            f"-src {fast_dki_temp_file_prefix}_eddyqc/"
            f"dwi_post_eddy_first_vol.nii.gz "
            f"-itk {temp_file_prefix}_from-fastDki_to-anat_affine.mat "
            f"-ras2fsl "
            f"-o {temp_file_prefix}_from-fastDki_to-anat_affine.txt"
        )
        run_command(
            "c3d_affine_tool "
            f"-ref {axon_radius_temp_file_prefix}_eddyqc/"
            f"dwi_post_eddy_first_vol.nii.gz "
            f"-src {axon_radius_temp_file_prefix}_eddyqc/"
            f"dwi_post_eddy_first_vol.nii.gz "
            f"-itk {temp_file_prefix}_from-b6000_to-anat_affine.mat "
            f"-ras2fsl "
            f"-o {temp_file_prefix}_from-b6000_to-anat_affine.txt"
        )
        run_command(
            "c3d_affine_tool "
            f"-ref {axon_radius_temp_file_prefix}_eddyqc/"
            f"dwi_post_eddy_first_vol.nii.gz "
            f"-src {fast_dki_temp_file_prefix}_eddyqc/"
            f"dwi_post_eddy_first_vol.nii.gz "
            f"-itk {temp_file_prefix}_from-b30450_to-anat_affine.mat "
            f"-ras2fsl "
            f"-o {temp_file_prefix}_from-b30450_to-anat_affine.txt"
        )

        # apply transformations up to here in a single step (eddy,
        # gradient distortions, alignment of different images in native space),
        # thereby transforming each image to native anatomical space.
        combine_warps(
            f"{fast_dki_temp_file_prefix}_eddyqc",
            f"{fast_dki_temp_file_prefix}_desc-preprocessed_dwi.nii.gz",
            extra_warp_file=f"{fast_dki_temp_file_prefix}"
                            f"_desc-gnlcAbs_warp.nii.gz",
            global_flirt=f"{temp_file_prefix}_from-fastDki_to-anat_affine.txt"
        )
        combine_warps(
            f"{axon_radius_temp_file_prefix}_eddyqc",
            f"{axon_radius_temp_file_prefix}"
            f"_desc-preprocessed_dwi.nii.gz",
            extra_warp_file=f"{axon_radius_temp_file_prefix}"
                            f"_desc-gnlcAbs_warp.nii.gz",
            flirt_mats={
                0: f"{temp_file_prefix}_from-b6000_to-anat_affine.txt",
                6000: f"{temp_file_prefix}_from-b6000_to-anat_affine.txt",
                30450: f"{temp_file_prefix}_from-b30450_to-anat_affine.txt"},
            bval_file=f"{axon_radius_temp_file_prefix}"
                      f"_desc-preprocessed_dwi.bval"
        )

        # convert to mrtrix format
        run_command(
            "mrconvert -force -strides -1,-2,3,4 "
            f"{fast_dki_temp_file_prefix}_desc-preprocessed_dwi.nii.gz "
            f"{fast_dki_temp_file_prefix}_desc-preprocessed_dwi.nii.gz"
        )
        run_command(
            "mrconvert -force -strides -1,-2,3,4 "
            f"{axon_radius_temp_file_prefix}_desc-preprocessed_dwi.nii.gz "
            f"{axon_radius_temp_file_prefix}_desc-preprocessed_dwi.nii.gz"
        )

        # merge preprocessed dwi images into single file
        run_command(
            "mrconvert -force "
            f"{fast_dki_temp_file_prefix}_desc-preprocessed_dwi.nii.gz "
            f"{fast_dki_temp_file_prefix}_desc-preprocessed_dwi.mif "
            f"-fslgrad "
            f"{fast_dki_temp_file_prefix}_desc-preprocessed_dwi.bvec "
            f"{fast_dki_temp_file_prefix}_desc-preprocessed_dwi.bval"
        )
        run_command(
            "mrconvert -force "
            f"{axon_radius_temp_file_prefix}_desc-preprocessed_dwi.nii.gz "
            f"{axon_radius_temp_file_prefix}_desc-preprocessed_dwi.mif "
            f"-fslgrad "
            f"{axon_radius_temp_file_prefix}_desc-preprocessed_dwi.bvec "
            f"{axon_radius_temp_file_prefix}_desc-preprocessed_dwi.bval"
        )
        run_command(
            "mrcat -force "
            f"{fast_dki_temp_file_prefix}_desc-preprocessed_dwi.mif "
            f"{axon_radius_temp_file_prefix}_desc-preprocessed_dwi.mif "
            f"{temp_file_prefix}_dwi.mif"
        )
        run_command(
            "mrconvert -force "
            f"{temp_file_prefix}_dwi.mif "
            f"{dwi_out_file_prefix}_dwi.nii.gz "
            f"-export_grad_fsl "
            f"{dwi_out_file_prefix}_dwi.bvec "
            f"{dwi_out_file_prefix}_dwi.bval"
        )

        # transform fastDki noisemap to b=6000 frame, assuming that the
        # lower b-values provide more accurate noise estimates than b>=6000
        # # create brain mask
        run_command(
            f"fslroi {dwi_out_file_prefix}_dwi.nii.gz "
            f"{temp_file_prefix}_desc-firstVol_dwi.nii.gz "
            f"0 1")
        run_command(
            "antsApplyTransforms -d 3 "
            f"-i {fast_dki_temp_file_prefix}_noisemap.nii.gz "
            f"-o {dwi_out_file_prefix}_noisemap.nii.gz "
            f"-r {temp_file_prefix}_desc-firstVol_dwi.nii.gz "
            f"-t {temp_file_prefix}_from-fastDki_to-anat_affine.mat"
        )

        # resample T1w brain mask to dwi resolution
        run_command(
            "mrtransform -force "
            f"{anat_out_file_prefix}_label-brain_mask.nii.gz "
            f"{dwi_out_file_prefix}_label-brain_mask.nii.gz "
            f"-template {dwi_out_file_prefix}_dwi.nii.gz "
            f"-interp nearest"
        )

        # compute spherical average
        matlab_engine.compute_powder_average_dwi(
            f"{dwi_out_file_prefix}_dwi.nii.gz",
            f"{dwi_out_file_prefix}_dwi.bvec",
            f"{dwi_out_file_prefix}_dwi.bval",
            f"{dwi_out_file_prefix}_powderAverage.nii.gz",
            "noisemap_file",
            f"{dwi_out_file_prefix}_noisemap.nii.gz",
            "mask_file",
            f"{dwi_out_file_prefix}_label-brain_mask.nii.gz",
            nargout=0
        )

        # compute effective radius
        matlab_engine.compute_effective_radius_power_law_ratio_dwi(
            f"{dwi_out_file_prefix}_powderAverage.nii.gz",
            f"{dwi_out_file_prefix}_powderAverage.bval",
            f"{dwi_out_file_prefix}_effectiveRadius.nii.gz",
            diffusion_gradient_separation,
            diffusion_gradient_duration,
            diffusivity_axoplasm,
            "mask_file",
            f"{dwi_out_file_prefix}_label-brain_mask.nii.gz",
            "min_bval", 6000,
            "max_bval", 30450,
            nargout=0
        )

        # extract DKI shells
        run_command(
            f"dwiextract -force "
            f"-shells 0,1000,2500 "
            f"{dwi_out_file_prefix}_dwi.nii.gz "
            f"{temp_file_prefix}_desc-dkiShells_dwi.nii.gz "
            f"-fslgrad "
            f"{dwi_out_file_prefix}_dwi.bvec "
            f"{dwi_out_file_prefix}_dwi.bval "
            f"-export_grad_fsl "
            f"{temp_file_prefix}_desc-dkiShells_dwi.bvec "
            f"{temp_file_prefix}_desc-dkiShells_dwi.bval")

        # compute DKI-based fractional anisotropy
        run_command(
            f"dwi2tensor -force "
            f"{temp_file_prefix}_desc-dkiShells_dwi.nii.gz "
            f"{temp_file_prefix}_D.nii.gz "
            f"-dkt {temp_file_prefix}_K.nii.gz "
            f"-fslgrad "
            f"{temp_file_prefix}_desc-dkiShells_dwi.bvec "
            f"{temp_file_prefix}_desc-dkiShells_dwi.bval")
        run_command(
            f"tensor2metric -force "
            f"{temp_file_prefix}_D.nii.gz "
            f"-fa {dwi_out_file_prefix}_FA.nii.gz")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process in-vivo MRI dataset and compute effective radius."
    )
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Path to output directory.")
    parser.add_argument("--subject", type=str, default=None,
                        help="Process only a specific subject.")

    args = parser.parse_args()
    process_in_vivo_dataset(args.output_dir, args.subject)
