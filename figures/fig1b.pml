# Panel B: arsitektur gorge - CAS, PAS, gerbang Tyr337, dan triad katalitik
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
set ambient, 0.24
set direct, 0.58
set orthoscopic, on
set cartoon_transparency, 0.62
set cartoon_side_chain_helper, 1

show cartoon, ache
color grey70, ache

select cas,  ache and resi 86  and sidechain
select pas,  ache and resi 286 and sidechain
select gate, ache and resi 337 and sidechain
select tri,  ache and resi 203+447+334 and sidechain
show sticks, cas or pas or gate or tri
set stick_radius, 0.19
color marine,    cas  and elem C
color firebrick, pas  and elem C
color orange,    gate and elem C
color forest,    tri  and elem C
util.cnc("cas or pas or gate or tri")

show sticks, dpz
set stick_radius, 0.26, dpz
color yellow, dpz and elem C
util.cnc("dpz")

# jarak CAS-PAS: rentang yang harus dijembatani bis-ligan
distance span, ache and resi 86 and name CZ3, ache and resi 286 and name CZ3
set dash_color, black
set dash_gap, 0.35
set dash_width, 3.0
set label_size, 26
set label_color, black
set label_font_id, 7

orient (cas or pas or dpz)
turn z, -90
turn y, -12
zoom (cas or pas or dpz or gate), 4.5
ray 1500, 1700
png figures/_p1b.png, dpi=300
