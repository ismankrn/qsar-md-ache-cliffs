# Gambar 1: struktur hAChE (4EY7) dengan gorge katalitik
reinitialize
load docking/4ey7_rec_clean.pdb, ache
load docking/4ey7_ligand_ref.pdb, dpz_xtal

hide everything
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_shadows, 0
set specular, 0.15
set cartoon_transparency, 0.15
set transparency, 0.55
set surface_quality, 1
set two_sided_lighting, on
set ambient, 0.18
set direct, 0.55
set reflect, 0.35
set orthoscopic, on

# rangka protein
show cartoon, ache
spectrum count, skyblue white salmon, ache and name CA
set cartoon_fancy_helices, 1

# permukaan molekuler semi-transparan untuk memperlihatkan kedalaman gorge
create surf, ache
show surface, surf
color grey70, surf

# residu penting pelapis gorge
select cas,   ache and resi 86            # Trp86  - dasar gorge
select pas,   ache and resi 286           # Trp286 - mulut gorge
select gate,  ache and resi 337           # Tyr337 - gerbang ayun
select triad, ache and resi 203+447+334   # triad katalitik
show sticks, (cas or pas or gate or triad) and sidechain
set stick_radius, 0.22
color marine,    cas   and elem C
color firebrick, pas   and elem C
color orange,    gate  and elem C
color forest,    triad and elem C
util.cnc("(cas or pas or gate or triad) and sidechain")

# donepezil kristal menandai lokasi gorge
show spheres, dpz_xtal
set sphere_scale, 0.32, dpz_xtal
color yellow, dpz_xtal and elem C
util.cnc("dpz_xtal")

# orientasi: sumbu gorge menghadap pembaca
orient (cas or pas or dpz_xtal)
turn x, -12
turn y, 8
zoom ache, 3

set label_size, 20
set label_color, black
set label_font_id, 7
set label_position, (0, 1.6, 0)
label cas and name CZ3,   "Trp86 (CAS)"
label pas and name CZ3,   "Trp286 (PAS)"
label gate and name OH,   "Tyr337"

ray 3000, 2250
png figures/fig1_protein_structure.png, dpi=300
