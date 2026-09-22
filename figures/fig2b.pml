# Panel B: huperzine A SELF-docking ke reseptornya sendiri (4EY5)
reinitialize
load docking/4ey5_rec_clean.pdb, rec
load docking/panel_HUP_xtal_4ey5.pdb, xtal
load docking/panel_HUP_self.sdf, dock
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

show sticks, xtal
show sticks, dock
set stick_radius, 0.26
color green,     xtal and elem C
color magenta, dock and elem C
util.cnc("xtal or dock")

orient xtal
turn x, -15
turn y, 18
zoom xtal, 6.5
ray 1900, 1700
png figures/_p2b.png, dpi=300
