# Panel C: huperzine A CROSS-docking ke 4EY7 - Tyr337 menghalangi pose sejati
reinitialize
load docking/4ey7_rec_clean.pdb, rec
load docking/panel_HUP_xtal_in4ey7.pdb, xtal
load docking/panel_HUP_cross.sdf, dock
load docking/panel_Y337_4ey5_in4ey7.pdb, y337_open
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

create cavity, rec
show surface, cavity
color grey72, cavity
set transparency, 0.38, cavity

select cas,  rec and resi 86  and sidechain
select gate, rec and resi 337 and sidechain
show sticks, cas or gate
set stick_radius, 0.17, (cas or gate)
color skyblue, cas  and elem C
color orange,  gate and elem C
util.cnc("cas or gate")

# rotamer Tyr337 milik 4EY5 (posisi yang dibutuhkan huperzine A)
show sticks, y337_open and sidechain
set stick_radius, 0.13, y337_open
color wheat, y337_open
set stick_transparency, 0.45, y337_open

show sticks, xtal
show sticks, dock
set stick_radius, 0.26, xtal
set stick_radius, 0.26, dock
color green,   xtal and elem C
color magenta, dock and elem C
util.cnc("xtal or dock")

orient xtal
turn x, -15
turn y, 18
zoom xtal, 6.5
ray 1900, 1700
png figures/_p2c.png, dpi=300
