"""
Utility functions for processing diffusion MRI data.

This module provides helper functions for diffusion MRI preprocessing,
including running shell commands, handling file operations, and applying
image processing steps such as smoothing and warping.
"""
import os
import re
import subprocess
from multiprocessing import Pool

import numpy as np
import psutil


def run_command(command_str, read_output=False, print_cmd=True):
    """
    Execute a shell command with optional output capture and logging.

    Parameters
    ----------
    command_str : str
        The command to execute as a string.
    read_output : bool, optional
        If True, captures and returns the command output. Defaults to False.
    print_cmd : bool, optional
        If True, prints the command before execution. Defaults to True.

    Returns
    -------
    str or None
        If `read_output` is True, returns the command output as a string.
        Otherwise, returns None.
    """

    if print_cmd:
        print("Command: {}".format(command_str))
    if read_output:
        return os.popen(command_str).read()
    else:
        subprocess.call(command_str, shell=True)


def strip_extension(filename):
    """
    Remove the file extension from a given filename.

    Parameters
    ----------
    filename : str
        The input filename.

    Returns
    -------
    str
        The filename without any extensions.
    """
    return filename.split(".", 1)[0]


def smooth_masked(map_in_file, map_out_file, mask_file, fwhm, temp_dir):
    """
        Smooth an input map within a masked region using Gaussian filtering.

        Parameters
        ----------
        map_in_file : str
            Path to the input map file.
        map_out_file : str
            Path to the output smoothed map file.
        mask_file : str
            Path to the binary mask file defining the smoothing region.
        fwhm : float
            Full width at half maximum (FWHM) for Gaussian smoothing.
        temp_dir : str
            Directory for storing temporary intermediate files.

        Returns
        -------
        None
            The smoothed and masked output map is written to `map_out_file`.
        """
    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir, exist_ok=True)

    # apply mask on effective radius map
    run_command(
        "mrcalc -force "
        f"{map_in_file} "
        f"{mask_file} "
        "-mul "
        f"{temp_dir}/map_masked.nii.gz"
    )

    # smooth masked effective radius map
    run_command(
        f"mrfilter -force "
        f"{temp_dir}/map_masked.nii.gz "
        f"smooth -fwhm {fwhm} "
        f"{temp_dir}/map_smoothed.nii.gz"
    )

    # create weights as smoothed map
    run_command(
        "mrfilter -force "
        f"{mask_file} "
        f"smooth -fwhm {fwhm} "
        f"- "
        f"| "
        f"mrcalc -force "
        f"- 0.000001 -max "
        f"{temp_dir}/weights.nii.gz"
    )

    # Normalize the smoothed values
    run_command(
        f"mrcalc -force "
        f"{temp_dir}/map_smoothed.nii.gz "
        f"{temp_dir}/weights.nii.gz "
        f"-div "
        f"{temp_dir}/map_weighted.nii.gz"
    )

    # Apply mask to region
    run_command(
        f"mrcalc -force "
        f"{temp_dir}/map_weighted.nii.gz "
        f"{mask_file} "
        f"-mul "
        f"{map_out_file}"
    )


def get_sorted_eddy_displacement_fields(folder_path):
    """
    Retrieve and return a list of sorted eddy displacement field file paths.

    Parameters
    ----------
    folder_path : str
        Path to the folder containing eddy displacement field files.

    Returns
    -------
    list of str
        Sorted list of file paths for eddy displacement fields.
    """
    pattern = re.compile(
        r"dwi_post_eddy\.eddy_displacement_fields\.(\d+)\.nii\.gz")
    matched_files = [
        (int(match.group(1)), os.path.join(folder_path, file))
        for file in os.listdir(folder_path) if (match := pattern.match(file))
    ]
    return [file for _, file in sorted(matched_files)]


def combine_warps(eddy_dir, out_file, extra_warp_file=None, bval_file=None,
                  flirt_mats=None, global_flirt=None):
    """
    Concatenates and applies warp fields from eddy correction, gradient
    nonlinearity correction (optional), and optionally a global FLIRT
    transformation or per-b-value FLIRT transforms.

    Source:
        Adapted from https://github.com/mrphysics-bonn/AxonDiameter
        (function: combine_warps)

    Copyright (c) 2024, Marten Veldmann <marten.veldmann@dzne.de>

    Licensed under the MIT License (see LICENSES/MIT.txt for full text).

    Parameters
    ----------
    eddy_dir : str
        Path to the directory containing eddy-corrected DWI data.
    out_file : str
        Output file path for the merged, corrected DWI data.
    extra_warp_file : str, optional
        Additional warp field to apply.
    bval_file : str, optional
        Path to the b-values file.
    flirt_mats : dict, optional
        Dictionary mapping b-values to FLIRT transformation matrices.
    global_flirt : str, optional
        Path to a single global FLIRT transformation matrix.

    Returns
    -------
    None
        The function writes the corrected DWI data to the specified output file.
    """
    # Ensure mutual exclusivity of FLIRT matrices
    if global_flirt and flirt_mats:
        raise ValueError(
            "Error: You cannot specify both a global matrix and per-b-value "
            "matrices. Choose one.")

    # Load b-values once (if provided)
    bvals = np.loadtxt(bval_file, dtype=int) if bval_file else None

    # Validate FLIRT matrices if provided
    if flirt_mats:
        if not isinstance(flirt_mats, dict):
            raise TypeError(
                "Error: flirt_mats must be a dictionary with b-values as keys "
                "and matrix filenames as values.")

        if bvals is None:
            raise ValueError(
                "Error: If providing per-b-value FLIRT transforms, a b-value "
                "file (-b) must also be provided.")

        # Ensure all provided b-values in flirt_mats exist in bval_file
        for b in flirt_mats.keys():
            if b not in bvals:
                raise ValueError(
                    f"Error: b-value {b} specified in flirt_mats is missing "
                    f"from the bval file.")

    # Paths
    dwi = os.path.join(eddy_dir, "eddy_outlier_free_data.nii.gz")
    dfields_dir = os.path.join(eddy_dir, "dfields")
    split_dir = os.path.join(eddy_dir, "split_data")
    comb_warp_dir = os.path.join(eddy_dir, "comb_warp")
    warp_dir = os.path.join(eddy_dir, "warped_data")

    # Ensure directories exist
    os.makedirs(split_dir, exist_ok=True)
    os.makedirs(comb_warp_dir, exist_ok=True)
    os.makedirs(warp_dir, exist_ok=True)

    # Split DWI
    run_command(" ".join(["fslsplit", dwi, split_dir + "/data", "-t"]))
    dwi_files = sorted(
        [os.path.join(split_dir, item) for item in os.listdir(split_dir)]
    )

    # Read displacement fields
    dfield_files = get_sorted_eddy_displacement_fields(dfields_dir)

    # Initialize b-values array (use zero if no bval file is given)
    bval = bvals if bvals is not None else np.zeros(len(dwi_files), dtype=int)

    # Initialize command lists

    cmd_list_convert = []
    cmd_list_apply = []
    cmd_list_avg_jac = []
    cmd_list_apply_jac = []

    for i, file in enumerate(dwi_files):
        comb_warp_file = os.path.join(comb_warp_dir, f"comb_warp.{i}")
        jac_file = os.path.join(comb_warp_dir, f"jac.{i}")
        warped_file = os.path.join(warp_dir, f"warped_vol.{i}")

        # Construct convertwarp command
        cmd = [
            "convertwarp",
            "-o", comb_warp_file,
            "-r", file,
            "-j", jac_file,
            f"--warp1={dfield_files[i]}"
        ]

        # Append extra warp file if provided
        if extra_warp_file:
            cmd.append(f"--warp2={extra_warp_file}")

        # Apply FLIRT transformation if specified
        if global_flirt:
            cmd.append(f"--postmat={global_flirt}")
        elif bval[i] in flirt_mats:
            cmd.append(f"--postmat={flirt_mats[bval[i]]}")

        cmd_list_convert.append(cmd)

        # Apply warp to image
        cmd_list_apply.append([
            "applywarp",
            "-i", file,
            "-r", file,
            "-o", warped_file,
            "-w", comb_warp_file,
            "--interp=spline"
        ])

        # Compute mean Jacobian determinant
        cmd_list_avg_jac.append(["fslmaths", jac_file, "-Tmean", jac_file])

        # Apply Jacobian determinant correction
        cmd_list_apply_jac.append(
            ["fslmaths", warped_file, "-mul", jac_file, warped_file])

    # Run commands in parallel
    cores = psutil.cpu_count(logical=False)
    pool = Pool(processes=cores)

    for cmd_list in [cmd_list_convert, cmd_list_apply, cmd_list_avg_jac,
                     cmd_list_apply_jac]:
        res = [pool.apply_async(run_command, (" ".join(cmd), ))
               for cmd in cmd_list]
        [val.get() for val in res]

    pool.close()

    # Merge DWI volumes
    cmd = ["fslmerge", "-t", out_file] + [
        os.path.join(warp_dir, f"warped_vol.{i}") for i in
        range(len(dwi_files))]
    run_command(" ".join(cmd))
