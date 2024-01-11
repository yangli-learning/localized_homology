# Localized Homology

A python implementation of the localized homology algorithm in the paper:

A. Zomorodian, G. Carlsson, Localized Homology, Computational Geometry 41 (2008) 126–148


## File structure
### Modules
- `blowup_complex`: main data structure for blowup complex
- `persistence`: localized persistence solver based on cechmate and phat library
- `complex` : construct simplicial complex from point cloud
- `simplicial_cover` : data structure for simplicial covers
- `datasets` : dataset tools for converting input format 
- `util` : misc helper functions

### Demos
- `main.py`: main demo file
- `test_semi_local.py`: test file for semi-localized persistence (in-progress)
- `notebooks/` : some jupytor notebooks for education purposes


-----
*contributors*: Yang Li, Zixi Zhao