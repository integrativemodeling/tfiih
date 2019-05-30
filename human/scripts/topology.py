from __future__ import print_function, absolute_import
import IMP
import IMP.algebra
import IMP.core
import IMP.pmi1.representation as representation
from numpy.random import rand as nrrand
from numpy import array
import ihm
import IMP.pmi1.mmcif
import sys

def make_topology():
    rbmaxtrans=0.7
    fbmaxtrans=0.7

    #####################################################
    #create hierarchies and rigid bodies and flexible parts
    #####################################################

    m=IMP.Model()
    simo=representation.Representation(m, upperharmonic=False, disorderedlength=True)

    if '--mmcif' in sys.argv:
        # Record the modeling protocol to an mmCIF file
        po = IMP.pmi1.mmcif.ProtocolOutput(open('human-tfiih.cif', 'w'))
        simo.add_protocol_output(po)
        po.system.title = ('Architecture of the human general transcription '
                           'and DNA repair factor TFIIH')
    else:
        po = None
    simo.dry_run = '--dry-run' in sys.argv

    # --- Get atomic models available
    # Sub component 1: Kinase 
    simo.create_component("ccl1")
    simo.add_component_sequence("ccl1", "../inputs/human_tfiih.fasta")
    simo.autobuild_model("ccl1",'../inputs/KIN28_CCL1.pdb',"B",resolutions=[1,30],resrange=(1,323),missingbeadsize=30,color=60./360,attachbeads=True)
    simo.setup_component_sequence_connectivity("ccl1", resolution=30)

    simo.create_component("kin28")
    simo.add_component_sequence("kin28", "../inputs/human_tfiih.fasta")
    simo.autobuild_model("kin28",'../inputs/KIN28_CCL1.pdb',"A",resolutions=[1,30],resrange=(1,346), missingbeadsize=30,color=185./360,attachbeads=True)
    simo.setup_component_sequence_connectivity("kin28", resolution=30)

    simo.create_component("tfb3")
    simo.add_component_sequence("tfb3", "../inputs/human_tfiih.fasta")
    simo.autobuild_model("tfb3",'../inputs/TFB3.pdb',"A",resolutions=[1,30],resrange=(1,309), missingbeadsize=30,color=285./360,attachbeads=True)
    simo.setup_component_sequence_connectivity("tfb3", resolution=30)

    # Sub component 2: Rad3
    simo.create_component("rad3")
    simo.add_component_sequence("rad3", "../inputs/human_tfiih.fasta")
    # The original modeling ignored gaps in this PDB.
    # To reproduce the modeling as closely as possible, force modern PMI
    # to ignore these gaps too:
    old_rg_hier = IMP.pmi1.tools.get_residue_gaps_in_hierarchy
    def dummy_rg_hier(hierarchy, start, end):
        # Return only the leading and trailing gap, plus the PDB range
        return [[start, 15, 'gap'], [16, 702, 'cont'], [703, end, 'gap']]
    try:
        IMP.pmi1.tools.get_residue_gaps_in_hierarchy = dummy_rg_hier
        simo.autobuild_model("rad3",'../inputs/RAD3.pdb'," ",resolutions=[1,30],resrange=(1,760), missingbeadsize=30,color=0.4,attachbeads=True)
    finally:
        IMP.pmi1.tools.get_residue_gaps_in_hierarchy = old_rg_hier

    simo.setup_component_sequence_connectivity("rad3", resolution=30)

    # Sub component 3: Ssl2 and Tfiicore
    simo.create_component("ssl2")
    simo.add_component_sequence("ssl2", "../inputs/human_tfiih.fasta")
    # The original modeling ignored gaps in this PDB.
    # To reproduce the modeling as closely as possible, force modern PMI
    # to ignore these gaps too:
    old_rg_hier = IMP.pmi1.tools.get_residue_gaps_in_hierarchy
    def dummy_rg_hier(hierarchy, start, end):
        # Return only the leading and trailing gap, plus the PDB range
        return [[start, 273, 'gap'], [274, 668, 'cont'], [669, end, 'gap']]
    try:
        IMP.pmi1.tools.get_residue_gaps_in_hierarchy = dummy_rg_hier
        simo.autobuild_model("ssl2",'../inputs/SSL2.pdb'," ",resolutions=[1,30],resrange=(1,782), missingbeadsize=30,color=0.0,attachbeads=True)
    finally:
        IMP.pmi1.tools.get_residue_gaps_in_hierarchy = old_rg_hier
    simo.setup_component_sequence_connectivity("ssl2", resolution=30)
      
    simo.create_component("tfb1")
    simo.add_component_sequence("tfb1", "../inputs/human_tfiih.fasta")
    simo.autobuild_model("tfb1",'../inputs/TFB1.pdb'," ",resolutions=[1,30],resrange=(1,548), missingbeadsize=30,color=1.0,attachbeads=True)
    simo.setup_component_sequence_connectivity("tfb1", resolution=30)

    simo.create_component("tfb2")
    simo.add_component_sequence("tfb2", "../inputs/human_tfiih.fasta")
    simo.autobuild_model("tfb2",'../inputs/TFB2_TFB5.pdb',"A",resolutions=[1,30],resrange=(1,462), missingbeadsize=30,color=0.9,attachbeads=True)
    simo.setup_component_sequence_connectivity("tfb2", resolution=30)

    simo.create_component("tfb4")
    simo.add_component_sequence("tfb4", "../inputs/human_tfiih.fasta")
    # The original modeling ignored gaps in this PDB.
    # To reproduce the modeling as closely as possible, force modern PMI
    # to ignore these gaps too:
    old_rg_hier = IMP.pmi1.tools.get_residue_gaps_in_hierarchy
    def dummy_rg_hier(hierarchy, start, end):
        # Return only the leading and trailing gap, plus the PDB range
        return [[start, 2, 'gap'], [3, 231, 'cont'], [232, end, 'gap']]
    try:
        IMP.pmi1.tools.get_residue_gaps_in_hierarchy = dummy_rg_hier
        simo.autobuild_model("tfb4",'../inputs/TFB4.pdb'," ",resolutions=[1,30],resrange=(1,308), missingbeadsize=30,color=0.8,attachbeads=True)
    finally:
        IMP.pmi1.tools.get_residue_gaps_in_hierarchy = old_rg_hier
    simo.setup_component_sequence_connectivity("tfb4", resolution=30)

    simo.create_component("tfb5")
    simo.add_component_sequence("tfb5", "../inputs/human_tfiih.fasta")
    # The original modeling ignored gaps in this PDB.
    # To reproduce the modeling as closely as possible, force modern PMI
    # to ignore these gaps too:
    old_rg_hier = IMP.pmi1.tools.get_residue_gaps_in_hierarchy
    def dummy_rg_hier(hierarchy, start, end):
        # Return only the trailing gap plus the PDB range
        return [[start, 60, 'cont'], [61, end, 'gap']]
    try:
        IMP.pmi1.tools.get_residue_gaps_in_hierarchy = dummy_rg_hier
        simo.autobuild_model("tfb5",'../inputs/TFB2_TFB5.pdb',"B",resolutions=[1,30],resrange=(2,71), missingbeadsize=30,color=0.,attachbeads=True)
    finally:
        IMP.pmi1.tools.get_residue_gaps_in_hierarchy = old_rg_hier
    simo.setup_component_sequence_connectivity("tfb5", resolution=30)

    simo.create_component("ssl1")
    simo.add_component_sequence("ssl1", "../inputs/human_tfiih.fasta")
    simo.autobuild_model("ssl1",'../inputs/SSL1.pdb'," ",resolutions=[1,30],resrange=(1,395), missingbeadsize=30,color=0.6,attachbeads=True)
    simo.setup_component_sequence_connectivity("ssl1", resolution=30)

    # Make a set of subunits into a single rigid body, and set its coordinates
    def make_rigid_with_coord(simo, subunits, coords):
        randomize_coords = lambda c: tuple(1.*(nrrand(3)-0.5)+array(c))
        rb = simo.set_rigid_bodies(subunits)
        rb.set_coordinates(randomize_coords(coords))

    # --- set rigid bodies
    # Sub component 1: Kinase
    emxk,emyk,emzk = 80.35, -23.66, 73.14   
    make_rigid_with_coord(simo, ["ccl1","kin28"],(emxk,emyk,emzk))

    # Sub component 2: Tfb3
    emxk,emyk,emzk = 17.25, -41.81, -21.93   
    make_rigid_with_coord(simo, ["tfb3"],(emxk,emyk,emzk))

    # Sub component 2: Rad3
    emxr,emyr,emzr = 37.25, 14.39, -1.38
    make_rigid_with_coord(simo, ["rad3"],(emxr,emyr,emzr))

    # Sub component 3: Ssl2 + Tfiihcore
    emxt,emyt,emzt = -29.65, -1.27, -0.32
    make_rigid_with_coord(simo, [("tfb2"),("tfb5")], (emxt,emyt,emzt))
    make_rigid_with_coord(simo, ["tfb1"], (emxt,emyt,emzt))
    make_rigid_with_coord(simo, ["tfb4"],(emxt,emyt,emzt))
    make_rigid_with_coord(simo, ["ssl1"],(emxt,emyt,emzt))
    make_rigid_with_coord(simo, ["ssl2"],(emxt,emyt,emzt))

    # --- set simulation params
    simo.set_floppy_bodies()
    d=simo.get_particles_to_sample()
    simo.set_rigid_bodies_max_trans(rbmaxtrans)
    simo.set_floppy_bodies_max_trans(fbmaxtrans)

    def shuffle_configuration_no_translation(simo, bounding_box_length=300.):
        "shuffle configuration, used to restart the optimization"
        "it only works if rigid bodies were initialized"
        if len(simo.rigid_bodies)==0:
            print("MultipleStates: rigid bodies were not intialized")
        hbbl=bounding_box_length/2
        ub = IMP.algebra.Vector3D(-hbbl,-hbbl,-hbbl)
        lb = IMP.algebra.Vector3D( hbbl, hbbl, hbbl)
        bb = IMP.algebra.BoundingBox3D(ub, lb)
        for rb in simo.rigid_bodies:
            translation = (rb.get_x(), rb.get_y(), rb.get_z())
            rotation = IMP.algebra.get_random_rotation_3d()
            transformation = IMP.algebra.Transformation3D(rotation, translation)
            rb.set_reference_frame(IMP.algebra.ReferenceFrame3D(transformation))
        for fb in simo.floppy_bodies:
            translation = IMP.algebra.get_random_vector_in(bb)
            IMP.core.XYZ(fb).set_coordinates(translation)

    #--- re-orient initial orientation only
    shuffle_configuration_no_translation(simo)
    return m, simo, po
