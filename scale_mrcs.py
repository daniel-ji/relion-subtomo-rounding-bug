import argparse
import mrcfile
import numpy as np
from scipy.ndimage import zoom

def scale_mrcs(input_file: str, input_voxel_size: float, output_file: str, scale_factor: float) -> None:
    # since an mrcs file is a stack of 2D images
    scale_factor = (1, scale_factor, scale_factor) 
    with mrcfile.open(input_file, mode='r') as mrc:
        data = mrc.data
        scaled_data = zoom(data, scale_factor, order=1)

    with mrcfile.new(output_file, overwrite=True) as mrc_out:
        mrc_out.set_data(scaled_data)
        mrc_out.voxel_size = input_voxel_size / scale_factor[1]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scale MRC files by a given factor.")
    parser.add_argument("--input-file", type=str, help="Path to the input MRC file.")
    parser.add_argument("--input-voxel-size", type=float, help="Voxel size of the input MRC file.")
    parser.add_argument("--output-file", type=str, help="Path to save the scaled MRC file.")
    parser.add_argument("--scale-factor", type=float, help="Scaling factor for the MRC data.")
    
    args = parser.parse_args()
    
    scale_mrcs(args.input_file, args.input_voxel_size, args.output_file, args.scale_factor)
    print(f"Scaled MRC file saved to {args.output_file}")