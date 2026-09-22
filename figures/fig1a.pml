# Panel A: keseluruhan lipatan hAChE, lokasi gorge ditandai donepezil
reinitialize
load docking/4ey7_rec_clean.pdb, ache
load docking/4ey7_ligand_ref.pdb, dpz
remove hydrogens
hide everything
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_shadows, 0
set specular, 0.12
set ambient, 0.22
set direct, 0.6
set orthoscopic, on
set cartoon_fancy_helices, 1
set cartoon_transparency, 0.0

show cartoon, ache
spectrum count, deepblue cyan green yellow red, ache and name CA
show spheres, dpz
set sphere_scale, 0.45, dpz
color grey20, dpz

orient ache
turn x, 90
zoom ache, 5
ray 1500, 1700
png figures/_p1a.png, dpi=300
