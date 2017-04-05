from __future__ import print_function, absolute_import
import IMP
import IMP.algebra
import IMP.core
import IMP.pmi.representation as representation
from numpy.random import rand as nrrand
from numpy import array

def make_topology():
    rbmaxtrans=0.7
    fbmaxtrans=0.7

    #####################################################
    #create hierarchies and rigid bodies and flexible parts
    #####################################################

    m=IMP.Model()
    simo=representation.SimplifiedModel(m, upperharmonic=False, disorderedlength=True)

    # --- Get atomic models available
    # Sub component 1: Kinase 
    simo.add_component_name("ccl1")
    simo.autobuild_model("ccl1",'../inputs/KIN28_CCL1.pdb',"B",resolutions=[1,30],resrange=(1,323),missingbeadsize=30,color=60./360,attachbeads=True)
    simo.setup_component_sequence_connectivity("ccl1", resolution=30)

    simo.add_component_name("kin28")
    simo.autobuild_model("kin28",'../inputs/KIN28_CCL1.pdb',"A",resolutions=[1,30],resrange=(1,346), missingbeadsize=30,color=185./360,attachbeads=True)
    simo.setup_component_sequence_connectivity("kin28", resolution=30)

    simo.add_component_name("tfb3")
    simo.autobuild_model("tfb3",'../inputs/TFB3.pdb',"A",resolutions=[1,30],resrange=(1,309), missingbeadsize=30,color=285./360,attachbeads=True)
    simo.setup_component_sequence_connectivity("tfb3", resolution=30)

    # Sub component 2: Rad3
    simo.add_component_name("rad3")
    simo.autobuild_model("rad3",'../inputs/RAD3.pdb'," ",resolutions=[1,30],resrange=(1,760), missingbeadsize=30,color=0.4,attachbeads=True)
    simo.setup_component_sequence_connectivity("rad3", resolution=30)

    # Sub component 3: Ssl2 and Tfiicore
    simo.add_component_name("ssl2")
    simo.autobuild_model("ssl2",'../inputs/SSL2.pdb'," ",resolutions=[1,30],resrange=(1,782), missingbeadsize=30,color=0.0,attachbeads=True)
    simo.setup_component_sequence_connectivity("ssl2", resolution=30)
      
    simo.add_component_name("tfb1")
    simo.autobuild_model("tfb1",'../inputs/TFB1.pdb'," ",resolutions=[1,30],resrange=(1,548), missingbeadsize=30,color=1.0,attachbeads=True)
    simo.setup_component_sequence_connectivity("tfb1", resolution=30)

    simo.add_component_name("tfb2")
    simo.autobuild_model("tfb2",'../inputs/TFB2_TFB5.pdb',"A",resolutions=[1,30],resrange=(1,462), missingbeadsize=30,color=0.9,attachbeads=True)
    simo.setup_component_sequence_connectivity("tfb2", resolution=30)

    simo.add_component_name("tfb4")
    simo.autobuild_model("tfb4",'../inputs/TFB4.pdb'," ",resolutions=[1,30],resrange=(1,308), missingbeadsize=30,color=0.8,attachbeads=True)
    simo.setup_component_sequence_connectivity("tfb4", resolution=30)

    simo.add_component_name("tfb5")
    simo.autobuild_model("tfb5",'../inputs/TFB2_TFB5.pdb',"B",resolutions=[1,30],resrange=(1,71), missingbeadsize=30,color=0.,attachbeads=True)
    simo.setup_component_sequence_connectivity("tfb5", resolution=30)

    simo.add_component_name("ssl1")
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
    return m, simo
