# RELION Subtomogram Extraction Bug: Tomogram Centering Discrepancy (Integer vs. Float Division)

If the user does not specify tomogram dimensions in the tomograms.star starfile, there can be a discrepancy in the centering of the tomogram during subtomogram extraction, causing a 0.5px shift, which can cause possibly noticable negative effects on downstream processing. 

## Bug Description

When extracting subtomograms, as the first transformation as part of the set of 3D affine transformations to generate the projection matrix (to project 3d pixel coordinates to 2D pixel coordinates to then crop), RELION does a centering of the tomogram, which is reliant on the tomogram dimensions. https://github.com/3dem/relion/blob/master/src/jaz/tomography/tomogram.cpp#L17-64 :

```cpp
void Tomogram::setProjectionMatrix(int f, RFLOAT xtilt, RFLOAT ytilt, RFLOAT zrot, RFLOAT xshift_angst, RFLOAT yshift_angst)
{

 /* From Alister Burt
 *
 * tilt_image_center = tilt_image_dimensions / 2
 * specimen_center = tomogram_dimensions / 2
 *
 * # Transformations, defined in order of application
 * s0 = S(-specimen_center)  # put specimen center-of-rotation at the origin
 * r0 = Rx(euler_angles['rlnTomoXTilt'])  # rotate specimen around X-axis
 * r1 = Ry(euler_angles['rlnTomoYTilt'])  # rotate specimen around Y-axis
 * r2 = Rz(euler_angles['rlnTomoZRot'])  # rotate specimen around Z-axis
 * s1 = S(specimen_shifts)  # shift projected specimen in xy (camera) plane
 * s2 = S(tilt_image_center)  # move specimen back into tilt-image coordinate system
 *
 * # compose matrices
 * transformations = s2 @ s1 @ r2 @ r1 @ r0 @ s0
 *
 */

    if (optics.pixelSize < 0.001)
        REPORT_ERROR("BUG: Tomogram::getProjectionMatrix encountered pixel size of zero");

    d4Matrix s0, s1, s2, r0, r1, r2;

    // Get specimen center
    t3Vector<double> specimen_center((double)int(w0/2), (double)int(h0/2), (double)int(d0/2) );
    s0 = s0.translation(-specimen_center);

    // Get specimen shifts (in pixels)
    t3Vector<double> specimen_shifts(xshift_angst / optics.pixelSize, yshift_angst / optics.pixelSize, 0.);
    s1 = s1.translation(specimen_shifts);

    // Get tilt image center
    std::vector<long int> tilt_image_center_int = stack.getSizeVector();
    t3Vector<double> tilt_image_center((double)int(tilt_image_center_int[0]/2), (double)int(tilt_image_center_int[1]/2), 0.);
    s2 = s2.translation(tilt_image_center);

    // Get rotation matrices
    t3Vector<double> xaxis(1., 0., 0.), yaxis(0., 1., 0.), zaxis(0., 0., 1.);
    r0 = r0.rotation(xaxis, xtilt);
    r1 = r1.rotation(yaxis, ytilt);
    r2 = r2.rotation(zaxis, zrot);

    projectionMatrices[f] = s1 * s2 * r2 * r1 * r0 * s0;

}
```

Specifically:

```cpp
    // Get specimen center
    t3Vector<double> specimen_center((double)int(w0/2), (double)int(h0/2), (double)int(d0/2) );
    s0 = s0.translation(-specimen_center);
```

Here, an **integer** divide is done, as opposed to other occurrences in the code, that have **float** divisions. In particular, the particle decentering done in `ParticleSet::getParticleCoordDecenteredPixel` utilizes a float-divided tomogram center `tomo_centre`:
https://github.com/3dem/relion/blob/master/src/jaz/tomography/particle_set.cpp#L429-L619

```cpp
d3Vector ParticleSet::getParticleCoordDecenteredPixel(ParticleIndex particle_id, const gravis::d3Vector &tomo_centre, RFLOAT tiltSeriesPixelSize) const
{
	d3Vector out;

    if (partTable.containsLabel(EMDL_IMAGE_CENT_COORD_X_ANGST) &&
        partTable.containsLabel(EMDL_IMAGE_CENT_COORD_Y_ANGST) &&
        partTable.containsLabel(EMDL_IMAGE_CENT_COORD_Z_ANGST))
    {
        partTable.getValue(EMDL_IMAGE_CENT_COORD_X_ANGST, out.x, particle_id.value);
        partTable.getValue(EMDL_IMAGE_CENT_COORD_Y_ANGST, out.y, particle_id.value);
        partTable.getValue(EMDL_IMAGE_CENT_COORD_Z_ANGST, out.z, particle_id.value);

        // Inside Jasenko's code, all coordinates are in decentered pixels, convert now from centered Angstroms (centre_tomo is in pixels)
        out /= tiltSeriesPixelSize;
        out += tomo_centre;
    }

...

	return out;
}
```

In the subtomogram extraction code, this `tomo_centre` value can be traced back to `SubtomoProgram::processTomograms`, where this function provides `tomogram.centre` ultimately as the `tomo_centre` argument (`SubtomoProgram::processTomograms` -> `ParticleSet::getTrajectoryInPixels` -> `ParticleSet::getPosition` -> `ParticleSet::getParticleCoordDecenteredPixel`):
https://github.com/3dem/relion/blob/master/src/jaz/tomography/programs/subtomo.cpp#L722

```cpp
            const std::vector<d3Vector> traj = particleSet.getTrajectoryInPixels(
                    part_id, fc, tomogram.centre, tomogram.optics.pixelSize, !apply_offsets);
```


`tomogram.centre` is derived from the tomogram dimensions, which is computed in the `TomogramSet::LoadTomogram` function (also called in `SubtomoProgram::processTomograms` some lines before `ParticleSet::getTrajectoryInPixels`) where the tomogram dimensions are set to the default values of (-999, -999, -999) if not provided in the tomograms.star file:

```cpp
Tomogram TomogramSet::loadTomogram(int index, bool loadImageData, bool loadEvenFramesOnly, bool loadOddFramesOnly,
                                   int _w0, int _h0, int _d0) const //Set loadEven/OddFramesOnly to True to loadImageData from rlnTomoMicrographNameEven/Odd rather than rlnMicrographName
{
	Tomogram out;

    std::string tomoName;
    i3Vector stackSize;

    globalTable.getValueSafely(EMDL_TOMO_NAME, tomoName, index);
    const MetaDataTable& m = tomogramTables[index];

    if (_w0 != -999 && _h0 != -999 && _d0 != -999)
    {
        out.w0 = _w0;
        out.h0 = _h0;
        out.d0 = _d0;
    }
    else if (globalTable.containsLabel(EMDL_TOMO_SIZE_X) &&
        globalTable.containsLabel(EMDL_TOMO_SIZE_Y) &&
        globalTable.containsLabel(EMDL_TOMO_SIZE_Z))
    {
        globalTable.getValueSafely(EMDL_TOMO_SIZE_X, out.w0, index);
        globalTable.getValueSafely(EMDL_TOMO_SIZE_Y, out.h0, index);
        globalTable.getValueSafely(EMDL_TOMO_SIZE_Z, out.d0, index);
    }
    else
    {
        out.w0 = -999;
        out.h0 = -999;
        out.d0 = -999;
    }
```

and later on in that same function:

```cpp
    out.centre = d3Vector(out.w0/2.0, out.h0/2.0, out.d0/2.0);
```

### And so here we can find a discrepancy in tomogram centre values (by 0.5px) between the `Tomogram::setProjectionMatrix` function and `ParticleSet::getParticleCoordDecenteredPixel` function when the tomogram dimensions are not provided, because the former uses an integer divide, while the latter uses a float divide on the default dimensions (-999, -999, -999).

## Some other integer divides that might be worth noting (but not directly related to the bug):

- Similar tilt image center calculations using integer division is done in `Tomogram::setProjectionMatrix` and `Tomogram::getProjectionAnglesFromMatrix`:
  - https://github.com/3dem/relion/blob/master/src/jaz/tomography/tomogram.cpp#L53
  - https://github.com/3dem/relion/blob/master/src/jaz/tomography/tomogram.cpp#L103
- `ParticleSet::getMatrix4x4` does a similar integer divide when generating a projection matrix
  - https://github.com/3dem/relion/blob/master/src/jaz/tomography/particle_set.cpp#L429-L457

```cpp
d4Matrix ParticleSet::getMatrix4x4(ParticleIndex particle_id, const gravis::d3Vector &tomo_centre, double w, double h, double d) const
{
        d3Matrix A = getMatrix3x3(particle_id);
        d3Vector pos = getPosition(particle_id, tomo_centre, true);

        int cx = ((int)w) / 2;
        int cy = ((int)h) / 2;
        int cz = ((int)d) / 2;

        gravis::d4Matrix Tc(
                1, 0, 0, -cx,
                0, 1, 0, -cy,
                0, 0, 1, -cz,
                0, 0, 0, 1);
    ...
}
```

## Bug Analysis & Repository Breakdown

### This repository's contents is intended to provide a minimal reproducible example to illustrate the bug in RELION subtomogram extraction.

It may be worth noting that the reason that the dimensions only matter in their parity (even/odd-ness) is because the provided particle picks are already in **centered** Angstrom coordinates (which later get converted to **centered** pixel coordinates), and so a translation to any arbitrary center (derived from the dimensions) will not affect the extraction, as long as the center values are consistent (which is not the case when the dimensions are odd and the integer and float divides yield different results). **So tomogram dimensions can be provided incorrectly, and as long as they are even, the extraction will work fine. And so in this repository we demonstrate this.**

### File breakdown

This repository contains a minimal example to illustrate this bug, with the tilt series and pick generated with synthetic data using [polnet](https://github.com/anmartinezs/polnet) with the following files:

- relion_project/input/particle.star: The particle star file containing the coordinates of two particles to be extracted.
- relion_project/input/tilt_series_filled.star: The tilt series star file containing the properly provided tomogram dimensions (630, 630, 200).
- relion_project/input/tilt_series_odd.star: The tilt series star file containing incorrect tomogram dimensions (odd - (313, 207, 115)) (used to illustrate that the dimensions only matter in their parity).
- relion_project/input/tilt_series_even.star: The tilt series star file containing incorrect tomogram dimensions (even - (1002, -402, 400)) (used to illustrate that the dimensions only matter in their parity).
- relion_project/input/tilt_series_blank.star: The tilt series star file for the blank particles (used to illustrate the buggy behavior of the default -999 dimensions).
- relion_project/input/tilt_series_-999.star: The tilt series star file for the -999 particles (used to illustrate the buggy behavior of the default -999 dimensions).
- relion_project/input/TS_1.star: The tilt series file containing per-tilt data.
- relion_project/input/TS_1.mrcs: The actual tilt series MRC data (2d_stack).


**Note: Scaling the tomogram (using `scipy.ndimage.zoom`) might induce some artifacts, so this has not been completely explored and tested in this repository yet.**
- relion_project/input/tilt_series_filled_scaled.star: The tilt series star file for the filled particles with the tiltseries scaled by 4x (to 2.5A) to simulate a binned tomogram.
- relion_project/input/tilt_series_blank_scaled.star: The tilt series star file for the blank particles with the tiltseries scaled by 4x (to 2.5A) to simulate a binned tomogram.
- relion_project/input/particle_scaled.star: The corresponding particle star file for the scaled tilt series, with an updated optics group.
- relion_project/input/TS_1_scaled.star: The tilt series file containing per-tilt data for the scaled tilt series.
- relion_project/input/TS_1_scaled.mrcs: The actual tilt series MRC data (2d_stack) for the scaled tilt series. **Note that this file is not included, as it is too large, it must be generated using: `python scale_mrcs.py --input-file relion_project/input/TS_1.mrcs --output-file relion_project/input/TS_1_scaled.mrcs --scale-factor 4 --input-voxel-size 10`**

The script `extract.sh` is used to extract the subtomograms from the tilt series using RELION's `relion_tomo_subtomo` command, which will generate the subtomograms in the `Extract/job001_*` directories.

The script `compare_subtomos.py` is used to compare the extracted subtomograms from the different starfiles (which only differ in tomogram dimensions), and it generates heatmaps and summary statistics of the differences.

### Viewing the differences in the subtomogram extraction RELION jobs:

The folders labeled `mrc_comparison_results_*` are the results of the comparisons between the different subtomogram extractions. The `mrc_comparison_results_1_*` files are for the first particle and the `mrc_comparison_results_2_*` files are for the second particle. Each folder contains a comparison at the specified slices/sections (1, 16, 31) of the subtomograms. The `*_absolute_difference.png` files show a the absolute difference at a per-pixel scale between the two subtomograms. For the absolute comparisons, the full orange heatmaps indicate that the two subtomograms are identical, while the heatmaps with a range of colors indicate the differences between the two subtomogram sections.

- "blank" vs "filled -999": Not including any tomogram dimensions vs providing dimensions (-999, -999, -999) are identical.
- "blank" vs "odd": Not including any tomogram dimensions vs providing odd dimensions (313, 207, 115) are identical.
- "even" vs "odd": Incorrect even dimensions (1002, -402, 400) vs incorrect odd dimensions (313, 207, 115) are different.
- "filled" vs "filled -999": Correctly provided dimensions (630, 630, 200) vs providing dimensions (-999, -999, -999) are different.
- "filled" vs "blank": Correctly provided dimensions (630, 630, 200) vs not including any dimensions is different.
- "filled" vs "even": Correctly provided dimensions (630, 630, 200) vs incorrect even dimensions (1002, -402, 400) are identical.
- "filled" vs "odd": Correctly provided dimensions (630, 630, 200) vs incorrect odd dimensions (313, 207, 115) are different.

### Commands used to generate the visualizations

```bash
rm -rf mrc_comparison_results*

# Compare filled (correct dimensions) vs blank (default -999 dimensions)
python compare_subtomos.py \
--mrc-file-1 relion_project/Extract/job001_filled/Subtomograms/session1_TS_1/1_stack2d.mrcs \
--mrc-file-2 relion_project/Extract/job001_blank/Subtomograms/session1_TS_1/1_stack2d.mrcs \
--mrc-file-1-name filled \
--mrc-file-2-name blank \
--output-dir mrc_comparison_results_1_filled_vs_blank \
--sections 1 16 31

python compare_subtomos.py \
--mrc-file-1 relion_project/Extract/job001_filled/Subtomograms/session1_TS_1/2_stack2d.mrcs \
--mrc-file-2 relion_project/Extract/job001_blank/Subtomograms/session1_TS_1/2_stack2d.mrcs \
--mrc-file-1-name filled \
--mrc-file-2-name blank \
--output-dir mrc_comparison_results_2_filled_vs_blank \
--sections 1 16 31

# Compare filled (correct dimensions) vs filled in -999 (incorrect dimensions, what default dimensions are set to in RELION)
python compare_subtomos.py \
--mrc-file-1 relion_project/Extract/job001_filled/Subtomograms/session1_TS_1/1_stack2d.mrcs \
--mrc-file-2 relion_project/Extract/job001_-999/Subtomograms/session1_TS_1/1_stack2d.mrcs \
--mrc-file-1-name filled \
--mrc-file-2-name "filled -999" \
--output-dir mrc_comparison_results_1_filled_vs_-999 \
--sections 1 16 31

python compare_subtomos.py \
--mrc-file-1 relion_project/Extract/job001_filled/Subtomograms/session1_TS_1/2_stack2d.mrcs \
--mrc-file-2 relion_project/Extract/job001_-999/Subtomograms/session1_TS_1/2_stack2d.mrcs \
--mrc-file-1-name filled \
--mrc-file-2-name "filled -999" \
--output-dir mrc_comparison_results_2_filled_vs_-999 \
--sections 1 16 31

# Compare filled (correct dimensions) vs even dimensions (incorrect dimensions, correct parity)
python compare_subtomos.py \
--mrc-file-1 relion_project/Extract/job001_filled/Subtomograms/session1_TS_1/1_stack2d.mrcs \
--mrc-file-2 relion_project/Extract/job001_even/Subtomograms/session1_TS_1/1_stack2d.mrcs \
--mrc-file-1-name filled \
--mrc-file-2-name "even" \
--output-dir mrc_comparison_results_1_filled_vs_even \
--sections 1 16 31

python compare_subtomos.py \
--mrc-file-1 relion_project/Extract/job001_filled/Subtomograms/session1_TS_1/2_stack2d.mrcs \
--mrc-file-2 relion_project/Extract/job001_even/Subtomograms/session1_TS_1/2_stack2d.mrcs \
--mrc-file-1-name filled \
--mrc-file-2-name "even" \
--output-dir mrc_comparison_results_2_filled_vs_even \
--sections 1 16 31

# Compare filled (correct dimensions) vs even dimensions (incorrect dimensions, odd parity)
python compare_subtomos.py \
--mrc-file-1 relion_project/Extract/job001_filled/Subtomograms/session1_TS_1/1_stack2d.mrcs \
--mrc-file-2 relion_project/Extract/job001_odd/Subtomograms/session1_TS_1/1_stack2d.mrcs \
--mrc-file-1-name filled \
--mrc-file-2-name "odd" \
--output-dir mrc_comparison_results_1_filled_vs_odd \
--sections 1 16 31

python compare_subtomos.py \
--mrc-file-1 relion_project/Extract/job001_filled/Subtomograms/session1_TS_1/2_stack2d.mrcs \
--mrc-file-2 relion_project/Extract/job001_odd/Subtomograms/session1_TS_1/2_stack2d.mrcs \
--mrc-file-1-name filled \
--mrc-file-2-name "odd" \
--output-dir mrc_comparison_results_2_filled_vs_odd \
--sections 1 16 31

# Compare blank (default -999 dimensions) vs filled in -999 (both are identical, to illustrate the bug)
python compare_subtomos.py \
--mrc-file-1 relion_project/Extract/job001_blank/Subtomograms/session1_TS_1/1_stack2d.mrcs \
--mrc-file-2 relion_project/Extract/job001_-999/Subtomograms/session1_TS_1/1_stack2d.mrcs \
--mrc-file-1-name "blank" \
--mrc-file-2-name "filled -999" \
--output-dir mrc_comparison_results_1_blank_vs_filled_-999 \
--sections 1 16 31

python compare_subtomos.py \
--mrc-file-1 relion_project/Extract/job001_blank/Subtomograms/session1_TS_1/2_stack2d.mrcs \
--mrc-file-2 relion_project/Extract/job001_-999/Subtomograms/session1_TS_1/2_stack2d.mrcs \
--mrc-file-1-name blank \
--mrc-file-2-name "filled -999" \
--output-dir mrc_comparison_results_2_blank_vs_filled_-999 \
--sections 1 16 31

# Compare blank (default -999 dimensions) vs odd (incorrect dimensions, odd parity) (both are identical, to illustrate the bug)
python compare_subtomos.py \
--mrc-file-1 relion_project/Extract/job001_blank/Subtomograms/session1_TS_1/1_stack2d.mrcs \
--mrc-file-2 relion_project/Extract/job001_odd/Subtomograms/session1_TS_1/1_stack2d.mrcs \
--mrc-file-1-name blank \
--mrc-file-2-name "odd" \
--output-dir mrc_comparison_results_1_blank_vs_odd \
--sections 1 16 31

python compare_subtomos.py \
--mrc-file-1 relion_project/Extract/job001_blank/Subtomograms/session1_TS_1/2_stack2d.mrcs \
--mrc-file-2 relion_project/Extract/job001_odd/Subtomograms/session1_TS_1/2_stack2d.mrcs \
--mrc-file-1-name blank \
--mrc-file-2-name "odd" \
--output-dir mrc_comparison_results_2_blank_vs_odd \
--sections 1 16 31

# Compare even (incorrect dimensions, even parity) vs odd (incorrect dimensions, odd parity)
python compare_subtomos.py \
--mrc-file-1 relion_project/Extract/job001_even/Subtomograms/session1_TS_1/1_stack2d.mrcs \
--mrc-file-2 relion_project/Extract/job001_odd/Subtomograms/session1_TS_1/1_stack2d.mrcs \
--mrc-file-1-name "even" \
--mrc-file-2-name "odd" \
--output-dir mrc_comparison_results_1_even_vs_odd \
--sections 1 16 31

python compare_subtomos.py \
--mrc-file-1 relion_project/Extract/job001_even/Subtomograms/session1_TS_1/2_stack2d.mrcs \
--mrc-file-2 relion_project/Extract/job001_odd/Subtomograms/session1_TS_1/2_stack2d.mrcs \
--mrc-file-1-name "even" \
--mrc-file-2-name "odd" \
--output-dir mrc_comparison_results_2_even_vs_odd \
--sections 1 16 31

# See above (File Breakdown) for why testing binned data is inconclusive and not fully explored yet.
# Compare filled (correct dimensions) vs blank (default -999 dimensions) for scaled data at 256px
# python compare_subtomos.py \
# --mrc-file-1 relion_project/Extract/job001_filled_scaled_256/Subtomograms/session1_TS_1/1_stack2d.mrcs \
# --mrc-file-2 relion_project/Extract/job001_blank_scaled_256/Subtomograms/session1_TS_1/1_stack2d.mrcs \
# --mrc-file-1-name "filled (scaled by 4)" \
# --mrc-file-2-name "blank (scaled by 4)" \
# --output-dir mrc_comparison_results_1_filled_vs_blank_scaled_256 \
# --sections 1 16 31

# python compare_subtomos.py \
# --mrc-file-1 relion_project/Extract/job001_filled_scaled_256/Subtomograms/session1_TS_1/2_stack2d.mrcs \
# --mrc-file-2 relion_project/Extract/job001_blank_scaled_256/Subtomograms/session1_TS_1/2_stack2d.mrcs \
# --mrc-file-1-name "filled (scaled by 4)" \
# --mrc-file-2-name "blank (scaled by 4)" \
# --output-dir mrc_comparison_results_2_filled_vs_blank_scaled_256 \
# --sections 1 16 31

# # Compare filled (correct dimensions) vs blank (default -999 dimensions) for scaled data at 64px (bin 1)
# python compare_subtomos.py \
# --mrc-file-1 relion_project/Extract/job001_filled_scaled_64/Subtomograms/session1_TS_1/1_stack2d.mrcs \
# --mrc-file-2 relion_project/Extract/job001_blank_scaled_64/Subtomograms/session1_TS_1/1_stack2d.mrcs \
# --mrc-file-1-name "filled (scaled by 4)" \
# --mrc-file-2-name "blank (scaled by 4)" \
# --output-dir mrc_comparison_results_1_filled_vs_blank_scaled_64 \
# --sections 1 16 31

# python compare_subtomos.py \
# --mrc-file-1 relion_project/Extract/job001_filled_scaled_64/Subtomograms/session1_TS_1/2_stack2d.mrcs \
# --mrc-file-2 relion_project/Extract/job001_blank_scaled_64/Subtomograms/session1_TS_1/2_stack2d.mrcs \
# --mrc-file-1-name "filled (scaled by 4)" \
# --mrc-file-2-name "blank (scaled by 4)" \
# --output-dir mrc_comparison_results_2_filled_vs_blank_scaled_64 \
# --sections 1 16 31

# # Compare filled (correct dimensions) vs blank (default -999 dimensions) for scaled data at 64px (bin 4)
# python compare_subtomos.py \
# --mrc-file-1 relion_project/Extract/job001_filled_scaled_64_bin_4/Subtomograms/session1_TS_1/1_stack2d.mrcs \
# --mrc-file-2 relion_project/Extract/job001_blank_scaled_64_bin_4/Subtomograms/session1_TS_1/1_stack2d.mrcs \
# --mrc-file-1-name "filled (scaled by 4, bin 4)" \
# --mrc-file-2-name "blank (scaled by 4, bin 4)" \
# --output-dir mrc_comparison_results_1_filled_vs_blank_scaled_64_bin_4 \
# --sections 1 16 31

# python compare_subtomos.py \
# --mrc-file-1 relion_project/Extract/job001_filled_scaled_64_bin_4/Subtomograms/session1_TS_1/2_stack2d.mrcs \
# --mrc-file-2 relion_project/Extract/job001_blank_scaled_64_bin_4/Subtomograms/session1_TS_1/2_stack2d.mrcs \
# --mrc-file-1-name "filled (scaled by 4, bin 4)" \
# --mrc-file-2-name "blank (scaled by 4, bin 4)" \
# --output-dir mrc_comparison_results_2_filled_vs_blank_scaled_64_bin_4 \
# --sections 1 16 31
```
