# spPipline
spatial data (dartfish) processing and analysis

## Installation
### Requirements
* [SimpleElastix](https://github.com/SuperElastix/SimpleElastix)
* [FIJI](https://downloads.imagej.net/fiji/latest/fiji-linux64.zip)
* [Cellpose](https://github.com/mouseland/cellpose)

- for current version, you can clone the github repository

```bash
git clone https://github.com/huqiwen0313/spPipline.git
cd spPipline
```  

## Usage
```bash
cd spPipline
python -m spPipline [arguments]
```
Also copy "_codebook" to the output path 

Arguments | Description | default
-----------|----------- | -------
`-r or --raw` | input path of Raw files. | "../0_Raw"
`-o or --output` | output path contains processed files. | "./"
`-s or --sigma` | sigma value for gaussian filter. | 0.7
`-rd or --rnd_list` | list of decoding cycle rounds. | ["0_anchor", "1_dc0", "2_dc1", "3_dc2", "4_dc3", "5_dc4", "6_dc5", "7_DRAQ5"]
`-nf or --nfovs` | number of fov (filed of view) | 80
`-cr or --channel_DIC_reference` | DIC channel for reference cycle | 'ch03'
`-cd or --channel_DIC` | DIC channel for (non-reference) decoding cycles | 'ch03'
`-co or ----cycle_other` | other data-containing folders which need to be aligned but are not names 'CycleXX' | ['0_anchor', '7_DRAQ5']
`-cdo or --channel_DIC_other` | DIC channel for other data-containing folders | {'0_anchor': 'ch01', '7_DRAQ5': 'ch01'}
`-gx or --grid_size_x` | grid size (x-axis) for each tile | 8
`-gy or --grid_size_y` | grid size (y-axis) for each tile | 10
`-to or --tile_overlap` | overlap for each tile | 15
`-ij or --ij_path` | The path to imagej’s executables | "/home/qiwenhu/software/Fiji.app/ImageJ-linux64" (need to set with your own path)
`-sr or --stitchRef` | The round to be used as the reference for stitching | "dc3"


## output file structure (processed data)
Example:

<pre>
SAMPLE (BUKMAP_20200303F_DFv1_ROI1)
    ├── _codebook
    ├── 1_Projected
    ├── 2_Registered
    │   ├── stitched
    │   └── FOVs
    ├── 3_Decoded
    │   ├── output_Starfish
    │   └── data_Starfish
    └── 4_CellAssignment
</pre>
