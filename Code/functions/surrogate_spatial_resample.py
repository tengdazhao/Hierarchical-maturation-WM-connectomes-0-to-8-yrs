from brainsmash.workbench.geo import volume
from brainsmash.mapgen.sampled import Sampled
import scipy.io as sio
coord_file = "voxel_coordinatesageType10NgE.nii.txt"
output_dir = "../0~8"

filenames = volume(coord_file, output_dir)

brain_map = "valueforsp.nii.txt"
gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'], resample=True)
surrogate_maps = gen(n=10000)
sio.savemat('surrogate_maps_type_resample.mat',{'surrogate_maps':surrogate_maps})
