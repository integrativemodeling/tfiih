"""Simple script to upgrade the RMF files in ../outputs/ to modern RMF/PMI
   style, so we can read them with modern IMP.

   Note:
    - Modern IMP checks to make sure that every floppy body has mass, and so
      will fail to read the old RMF files (since they don't have mass). You will
      need to use a modified version of IMP with the test in
      modules/atom/src/Hierarchy.cpp disabled.
    - IMP.rmf will try to score the restraints read from the old RMF files,
      which won't work since no internal coordinates are present. Work around
      this by disabling the rsf_->evaluate() call in
      modules/rmf/src/restraint_io.cpp.
"""

from __future__ import print_function
import sys
import os
import IMP.rmf
import RMF
import tempfile
import shutil

def fix_residue_names(res, comp_name):
    for c in res.get_children():
        name = c.get_name()
        if name.startswith('%s_' % comp_name):
            name = name[len(comp_name) + 1:]
        if name.endswith('_pdb'):
            name = name[:-4]
        c.set_name(name)

def fix_fragment_names(res, comp_name):
    for c in res.get_children():
        c.set_name(c.get_name().replace('-', '_'))

def handle_ssl2_res1(m, res1, beads, c):
    """Special handling for the 298-538 resolution=1 part of ssl2"""
    # Old RMF files represent this region with a single 298-538 fragment,
    # even though the input file has gaps at 308-311, 318-322 and 340-344.
    # Split into multiple fragments and create an (approximated) bead for
    # each gap.
    children = {}
    for child in c.get_children():
        c.remove_child(child)
        children[int(child.get_name())] = child
    for fragment_range in [(298, 307), (312, 317), (323, 339), (345, 538)]:
        name = 'ssl2_%d-%d_pdb' % fragment_range
        f = IMP.atom.Hierarchy.setup_particle(IMP.Particle(m))
        f.set_name(name)
        res1.add_child(f)
        for i in range(fragment_range[0], fragment_range[1] + 1):
            f.add_child(children[i])
    for bead_range in [(308, 311), (318, 322), (340, 344)]:
        name = 'ssl2_%d-%d_bead_floppy_body_rigid_body_member' % bead_range
        f = IMP.atom.Hierarchy.setup_particle(IMP.Particle(m))
        f.set_name(name)
        beads.add_child(f)
        # Approximate the coordinates as the average of the residues on
        # either side of the gap
        xyz1 = IMP.core.XYZ(children[bead_range[0] - 1]).get_coordinates()
        xyz2 = IMP.core.XYZ(children[bead_range[1] + 1]).get_coordinates()
        s = IMP.algebra.Sphere3D((xyz1 + xyz2) / 2., 1.)
        IMP.core.XYZR.setup_particle(f, s)
    IMP.atom.destroy(c)

def upgrade_component(comp):
    beads = IMP.atom.Hierarchy.setup_particle(IMP.Particle(m))
    beads.set_name("Beads")

    res1 = IMP.atom.Hierarchy.setup_particle(IMP.Particle(m))
    res1.set_name("%s_Res:1" % comp.get_name())
    res30 = IMP.atom.Hierarchy.setup_particle(IMP.Particle(m))
    res30.set_name("%s_Res:30" % comp.get_name())

    for c in comp.get_children():
        comp.remove_child(c)
        name = c.get_name()
        if 'bead' in name:
            beads.add_child(c)
        elif name.endswith('_Res:1'):
            fix_residue_names(c, comp.get_name())
            c.set_name(name[:-6])
            if name == 'ssl2_298-538_pdb_Res:1':
                handle_ssl2_res1(m, res1, beads, c)
            else:
                res1.add_child(c)
        elif name.endswith('_Res:30'):
            fix_fragment_names(c, comp.get_name())
            c.set_name(name[:-7])
            res30.add_child(c)
        else:
            raise ValueError("Don't know what to do with %s" % name)

    if comp.get_name() in ('tfb1', 'tfb5'):
        comp.add_child(res1)
        comp.add_child(res30)
        comp.add_child(beads)
    else:
        comp.add_child(beads)
        comp.add_child(res1)
        comp.add_child(res30)

m = IMP.Model()

rh_in = RMF.open_rmf_file_read_only(sys.argv[1])
num_frames = rh_in.get_number_of_frames()
del rh_in
tempdir = tempfile.mkdtemp()

for i in range(num_frames):
    rh_in = RMF.open_rmf_file_read_only(sys.argv[1])
    h = IMP.rmf.create_hierarchies(rh_in, m)
    rs = IMP.rmf.create_restraints(rh_in, m)
    IMP.rmf.load_frame(rh_in, RMF.FrameID(i))
    del rh_in
    for comp in h[0].get_children():
        upgrade_component(comp)
    rh_out = RMF.create_rmf_file(os.path.join(tempdir, '%d.rmf' % i))
    IMP.rmf.add_hierarchies(rh_out, h)
    IMP.rmf.add_restraints(rh_out, rs)
    IMP.rmf.save_frame(rh_out)
    del rh_out
    for x in h:
        IMP.atom.destroy(x)

rh_out = RMF.create_rmf_file(sys.argv[2])
rh_in = RMF.open_rmf_file_read_only(os.path.join(tempdir, '0.rmf'))
h = IMP.rmf.create_hierarchies(rh_in, m)
rs = IMP.rmf.create_restraints(rh_in, m)
IMP.rmf.add_hierarchies(rh_out, h)
IMP.rmf.add_restraints(rh_out, rs)
IMP.rmf.save_frame(rh_out)

for i in range(1, num_frames):
    rh_in = RMF.open_rmf_file_read_only(os.path.join(tempdir, '%d.rmf' % i))
    IMP.rmf.link_hierarchies(rh_in, h)
    IMP.rmf.link_restraints(rh_in, rs)
    IMP.rmf.load_frame(rh_in, RMF.FrameID(0))
    IMP.rmf.save_frame(rh_out)
    del rh_in
del rh_out

shutil.rmtree(tempdir)
