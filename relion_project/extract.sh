`which relion_tomo_subtomo` --p input/particle.star --t input/tilt_series_filled.star --theme classic --o Extract/job001_filled/ --b 64 --bin 1 --min_frames 1 --float16  --stack2d  --j 1  --pipeline_control Extract/job001_filled/
`which relion_tomo_subtomo` --p input/particle.star --t input/tilt_series_odd.star --theme classic --o Extract/job001_odd/ --b 64 --bin 1 --min_frames 1 --float16  --stack2d  --j 1  --pipeline_control Extract/job001_odd/
`which relion_tomo_subtomo` --p input/particle.star --t input/tilt_series_even.star --theme classic --o Extract/job001_even/ --b 64 --bin 1 --min_frames 1 --float16  --stack2d  --j 1  --pipeline_control Extract/job001_even/
`which relion_tomo_subtomo` --p input/particle.star --t input/tilt_series_blank.star --theme classic --o Extract/job001_blank/ --b 64 --bin 1 --min_frames 1 --float16  --stack2d  --j 1  --pipeline_control Extract/job001_blank/
`which relion_tomo_subtomo` --p input/particle.star --t input/tilt_series_-999.star --theme classic --o Extract/job001_-999/ --b 64 --bin 1 --min_frames 1 --float16  --stack2d  --j 1  --pipeline_control Extract/job001_-999/


`which relion_tomo_subtomo` --p input/particle_scaled.star --t input/tilt_series_filled_scaled.star --theme classic --o Extract/job001_filled_scaled_256/ --b 256 --bin 1 --min_frames 1 --float16  --stack2d  --j 1  --pipeline_control Extract/job001_filled_scaled_256/
`which relion_tomo_subtomo` --p input/particle_scaled.star --t input/tilt_series_blank_scaled.star --theme classic --o Extract/job001_blank_scaled_256/ --b 256 --bin 1 --min_frames 1 --float16  --stack2d  --j 1  --pipeline_control Extract/job001_blank_scaled_256/
`which relion_tomo_subtomo` --p input/particle_scaled.star --t input/tilt_series_filled_scaled.star --theme classic --o Extract/job001_filled_scaled_64/ --b 64 --bin 1 --min_frames 1 --float16  --stack2d  --j 1  --pipeline_control Extract/job001_filled_scaled_64/
`which relion_tomo_subtomo` --p input/particle_scaled.star --t input/tilt_series_blank_scaled.star --theme classic --o Extract/job001_blank_scaled_64/ --b 64 --bin 1 --min_frames 1 --float16  --stack2d  --j 1  --pipeline_control Extract/job001_blank_scaled_64/

# See note in README.md about why binning has not been tested as part of analyzing the bug.
# `which relion_tomo_subtomo` --p input/particle_scaled.star --t input/tilt_series_filled_scaled.star --theme classic --o Extract/job001_filled_scaled_64_bin_4/ --b 64 --bin 4 --min_frames 1 --float16  --stack2d  --j 1  --pipeline_control Extract/job001_filled_scaled_64_bin_4/
# `which relion_tomo_subtomo` --p input/particle_scaled.star --t input/tilt_series_blank_scaled.star --theme classic --o Extract/job001_blank_scaled_64_bin_4/ --b 64 --bin 4 --min_frames 1 --float16  --stack2d  --j 1  --pipeline_control Extract/job001_blank_scaled_64_bin_4/
