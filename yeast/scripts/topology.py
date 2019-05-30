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

activated_components=[1,2,3,4]

def make_topology():
    rbmaxtrans=0.5
    fbmaxtrans=0.5

    #####################################################
    #create hierarchies and rigid bodies and flexible parts
    #####################################################

    m=IMP.Model()
    simo=representation.Representation(m, upperharmonic=False, disorderedlength=True)

    if '--mmcif' in sys.argv:
        # Record the modeling protocol to an mmCIF file
        po = IMP.pmi1.mmcif.ProtocolOutput(open('yeast-tfiih.cif', 'w'))
        simo.add_protocol_output(po)
        po.system.title = ('Architecture of the yeast general transcription '
                           'and DNA repair factor TFIIH')
        # Add publication
        po.system.citations.append(ihm.Citation.from_pubmed_id(26340423))
    else:
        po = None
    simo.dry_run = '--dry-run' in sys.argv

    if 1 in activated_components:
      
      simo.create_component("ccl1")
      simo.add_component_sequence("ccl1", "../inputs/yeast_tfiih.fasta")
      simo.autobuild_model("ccl1",'../inputs/CCL1_48-393.pdb',"A",resolutions=[1,30],resrange=(1,393),missingbeadsize=30,color=1.0,attachbeads=True)
      simo.setup_component_sequence_connectivity("ccl1", resolution=30)

      simo.create_component("kin28")
      simo.add_component_sequence("kin28", "../inputs/yeast_tfiih.fasta")
      simo.autobuild_model("kin28",'../inputs/KIN28_5-299.pdb',"A",resolutions=[1,30],resrange=(1,306), missingbeadsize=30,color=0.9,attachbeads=True)
      simo.setup_component_sequence_connectivity("kin28", resolution=30)

      simo.create_component("tfb3")
      simo.add_component_sequence("tfb3", "../inputs/yeast_tfiih.fasta")
      simo.autobuild_model("tfb3",'../inputs/TFB3_8-142.pdb',"A",resolutions=[1,30],resrange=(1,321), missingbeadsize=30,color=0.8,attachbeads=True)
      simo.setup_component_sequence_connectivity("tfb3", resolution=30)


    if 2 in activated_components:

      simo.create_component("rad3")
      simo.add_component_sequence("rad3", "../inputs/yeast_tfiih.fasta")
      simo.autobuild_model("rad3",'../inputs/RAD3_14-725.pdb',"A",resolutions=[1,30],resrange=(1,778), missingbeadsize=30,color=150./360,attachbeads=True)
      simo.setup_component_sequence_connectivity("rad3", resolution=30)


    if 4 in activated_components:

      simo.create_component("tfb1")
      simo.add_component_sequence("tfb1", "../inputs/yeast_tfiih.fasta")
      simo.autobuild_model("tfb1",'../inputs/TFB1_2-115.pdb',"A",resolutions=[1,30],resrange=(2,642), missingbeadsize=30,color=60./360,attachbeads=True)
      simo.setup_component_sequence_connectivity("tfb1", resolution=30)

      simo.create_component("tfb2")
      simo.add_component_sequence("tfb2", "../inputs/yeast_tfiih.fasta")
      #simo.autobuild_model("tfb2",'../inputs/TFB2_1-168.pdb',"A",resolutions=[1,30],resrange=(1,170), missingbeadsize=30,color=185./360,attachbeads=True)
      #simo.autobuild_model("tfb2",'../inputs/TFB2_186-417.pdb',"A",resolutions=[1,30],resrange=(171,417), missingbeadsize=30,color=185./360,attachbeads=True)
      simo.autobuild_model("tfb2",'../inputs/TFB2_392-513.pdb',"A",resolutions=[1,30],resrange=(1,513), missingbeadsize=30,color=185./360,attachbeads=True)
      simo.setup_component_sequence_connectivity("tfb2", resolution=30)

      simo.create_component("tfb4")
      simo.add_component_sequence("tfb4", "../inputs/yeast_tfiih.fasta")
      simo.autobuild_model("tfb4",'../inputs/TFB4_24-250.pdb',"A",resolutions=[1,30],resrange=(1,338), missingbeadsize=30,color=210./360,attachbeads=True)
      simo.setup_component_sequence_connectivity("tfb4", resolution=30)

      simo.create_component("tfb5")
      simo.add_component_sequence("tfb5", "../inputs/yeast_tfiih.fasta")
      simo.autobuild_model("tfb5",'../inputs/TFB5.pdb',"B",resolutions=[1,30],resrange=(2,72), missingbeadsize=30,color=0.,attachbeads=True)
      simo.setup_component_sequence_connectivity("tfb5", resolution=30)

      simo.create_component("ssl1")
      simo.add_component_sequence("ssl1", "../inputs/yeast_tfiih.fasta")
      simo.autobuild_model("ssl1",'../inputs/SSL1_123-302.pdb',"A",resolutions=[1,30],resrange=(1,302), missingbeadsize=30,color=285./360,attachbeads=True)
      simo.autobuild_model("ssl1",'../inputs/SSL1_386-455.pdb',"A",resolutions=[1,30],resrange=(303,461), missingbeadsize=30,color=285./360,attachbeads=True)
      simo.setup_component_sequence_connectivity("ssl1", resolution=30)


    if 1 in activated_components:
      simo.create_component("ssl2")
      simo.add_component_sequence("ssl2", "../inputs/yeast_tfiih.fasta")
      # The original modeling ignored the gaps 308-311, 318-322 and 340-344.
      # To reproduce the modeling as closely as possible, force modern PMI
      # to ignore these gaps too:
      old_rg_hier = IMP.pmi1.tools.get_residue_gaps_in_hierarchy
      def dummy_rg_hier(hierarchy, start, end):
          # Return only the first two entries (leading gap, the missing
          # residues at the N terminus; plus the pdb range)
          return [[start, 297, 'gap'], [298, end, 'cont']]
      try:
          IMP.pmi1.tools.get_residue_gaps_in_hierarchy = dummy_rg_hier
          simo.autobuild_model("ssl2",'../inputs/fit_gtfs_3.pdb',"C",resolutions=[1,30],resrange=(1,538), missingbeadsize=30,color=0.0,attachbeads=True)
      finally:
          IMP.pmi1.tools.get_residue_gaps_in_hierarchy = old_rg_hier
      simo.autobuild_model("ssl2",'../inputs/fit_gtfs_3.pdb',"D",resolutions=[1,30],resrange=(539,843), missingbeadsize=30,color=0.0,attachbeads=True)
      simo.setup_component_sequence_connectivity("ssl2", resolution=30)


    # Make a set of subunits into a single rigid body, and set its coordinates
    def make_rigid_with_coord(simo, subunits, coords):
        randomize_coords = lambda c: tuple(1.*(nrrand(3)-0.5)+array(c))
        rb = simo.set_rigid_bodies(subunits)
        rb.set_coordinates(randomize_coords(coords))

    # --- set rigid bodies
    if 1 in activated_components:
      emxk,emyk,emzk = 215.249,220.731,167.764 #21
      make_rigid_with_coord(simo, ["tfb3"],(emxk,emyk,emzk))
      make_rigid_with_coord(simo, ["ccl1","kin28"],(emxk,emyk,emzk))

    if 2 in activated_components:
      emxr,emyr,emzr = 147.133,219.683,188.062 #0
      make_rigid_with_coord(simo, ["rad3"],(emxr,emyr,emzr))

    if 4 in activated_components:
      emxt,emyt,emzt = 126.646,130.588,192.099 #90
      make_rigid_with_coord(simo, [("tfb2",(437,513)),("tfb5",(1,72))], (emxt,emyt,emzt))
      make_rigid_with_coord(simo, [("tfb2",(1,436))], (emxt,emyt,emzt))
      #make_rigid_with_coord(simo, [("tfb2",(1,168))], (emxt,emyt,emzt))
      #make_rigid_with_coord(simo, [("tfb2",(169,436))], (emxt,emyt,emzt))
      make_rigid_with_coord(simo, ["tfb1"], (emxt,emyt,emzt))
      make_rigid_with_coord(simo, ["tfb4"],(emxt,emyt,emzt))
      make_rigid_with_coord(simo, [("ssl1",(1,385))],(emxt,emyt,emzt))
      make_rigid_with_coord(simo, [("ssl1",(386,461))],(emxt,emyt,emzt))


    if 3 in activated_components:
      emxk,emyk,emzk = 212.272,129.032,178.206 #80
      make_rigid_with_coord(simo, ["ssl2"],(emxk,emyk,emzk))



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

    #re-orient initial positions
    shuffle_configuration_no_translation(simo)
    return m, simo, po
