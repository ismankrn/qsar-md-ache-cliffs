# Gambar 2: validasi redocking donepezil + rongga pengikatan
reinitialize
load docking/4ey7_rec_clean.pdb, ache
load docking/4ey7_ligand_ref.pdb, xtal
load docking/E20_docked_pose1.sdf, docked
remove hydrogens
hide everything
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_shadows, 0
set specular, 0.10
set ambient, 0.28
set direct, 0.52
set orthoscopic, on
set two_sided_lighting, on
set surface_quality, 1

# objek permukaan TERPISAH supaya pewarnaan residu tidak merembes ke permukaan
create cavity, ache
show surface, cavity
color grey72, cavity
set transparency, 0.38, cavity

# residu penanda CAS / PAS
select cas,  ache and resi 86  and sidechain
select gate, ache and resi 337 and sidechain
select pas, ache and resi 286 and sidechain
show sticks, cas or pas or gate
set stick_radius, 0.17, (cas or pas or gate)
color skyblue,  cas and elem C
color salmon,   pas and elem C
color orange,   gate and elem C
util.cnc("cas or pas or gate")

# dua pose donepezil - inti gambar
show sticks, xtal
show sticks, docked
set stick_radius, 0.26, xtal
set stick_radius, 0.26, docked
color green,   xtal   and elem C
color magenta, docked and elem C
util.cnc("xtal or docked")

orient xtal
turn x, -15
turn y, 18
zoom xtal, 6.5
ray 1900, 1700
png figures/_p2.png, dpi=300
