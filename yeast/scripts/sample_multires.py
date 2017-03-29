#!/usr/bin/env python
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import os
import IMP.rmf
import RMF
from numpy import array
from numpy.random import rand as nrrand
import IMP.pmi.restraints as restraints
import IMP.pmi.representation as representation
import IMP.pmi.tools as tools
import IMP.pmi.samplers as samplers
import IMP.pmi.output as output

ncycl=100  #number of montecarlo steps cycles
rbmaxtrans=0.5
fbmaxtrans=0.5
outputobjects=[]
sampleobjects=[]




#####################################################
#create hierarchies and rigid bodies and flexible parts
#####################################################

m=IMP.Model()
simo=representation.SimplifiedModel(m, upperharmonic=False, disorderedlength=True)



activated_components=[1,2,3,4]


if 1 in activated_components:
  
  simo.add_component_name("ccl1")
  simo.autobuild_model("ccl1",'../inputs/CCL1_48-393.pdb',"A",resolutions=[1,30],resrange=(1,393),missingbeadsize=30,color=1.0,attachbeads=True)
  simo.setup_component_sequence_connectivity("ccl1", resolution=30)

  simo.add_component_name("kin28")
  simo.autobuild_model("kin28",'../inputs/KIN28_5-299.pdb',"A",resolutions=[1,30],resrange=(1,306), missingbeadsize=30,color=0.9,attachbeads=True)
  simo.setup_component_sequence_connectivity("kin28", resolution=30)

  simo.add_component_name("tfb3")
  simo.autobuild_model("tfb3",'../inputs/TFB3_8-142.pdb',"A",resolutions=[1,30],resrange=(1,321), missingbeadsize=30,color=0.8,attachbeads=True)
  simo.setup_component_sequence_connectivity("tfb3", resolution=30)


if 2 in activated_components:

  simo.add_component_name("rad3")
  simo.autobuild_model("rad3",'../inputs/RAD3_14-725.pdb',"A",resolutions=[1,30],resrange=(1,778), missingbeadsize=30,color=150./360,attachbeads=True)
  simo.setup_component_sequence_connectivity("rad3", resolution=30)


if 4 in activated_components:

  simo.add_component_name("tfb1")
  simo.autobuild_model("tfb1",'../inputs/TFB1_2-115.pdb',"A",resolutions=[1,30],resrange=(1,642), missingbeadsize=30,color=60./360,attachbeads=True)
  simo.setup_component_sequence_connectivity("tfb1", resolution=30)

  simo.add_component_name("tfb2")
  #simo.autobuild_model("tfb2",'../inputs/TFB2_1-168.pdb',"A",resolutions=[1,30],resrange=(1,170), missingbeadsize=30,color=185./360,attachbeads=True)
  #simo.autobuild_model("tfb2",'../inputs/TFB2_186-417.pdb',"A",resolutions=[1,30],resrange=(171,417), missingbeadsize=30,color=185./360,attachbeads=True)
  simo.autobuild_model("tfb2",'../inputs/TFB2_392-513.pdb',"A",resolutions=[1,30],resrange=(418,513), missingbeadsize=30,color=185./360,attachbeads=True)
  simo.setup_component_sequence_connectivity("tfb2", resolution=30)

  simo.add_component_name("tfb4")
  simo.autobuild_model("tfb4",'../inputs/TFB4_24-250.pdb',"A",resolutions=[1,30],resrange=(1,338), missingbeadsize=30,color=210./360,attachbeads=True)
  simo.setup_component_sequence_connectivity("tfb4", resolution=30)

  simo.add_component_name("tfb5")
  simo.autobuild_model("tfb5",'../inputs/TFB5.pdb',"B",resolutions=[1,30],resrange=(1,72), missingbeadsize=30,color=0.,attachbeads=True)
  simo.setup_component_sequence_connectivity("tfb5", resolution=30)

  simo.add_component_name("ssl1")
  simo.autobuild_model("ssl1",'../inputs/SSL1_123-302.pdb',"A",resolutions=[1,30],resrange=(1,302), missingbeadsize=30,color=285./360,attachbeads=True)
  simo.autobuild_model("ssl1",'../inputs/SSL1_386-455.pdb',"A",resolutions=[1,30],resrange=(303,461), missingbeadsize=30,color=285./360,attachbeads=True)
  simo.setup_component_sequence_connectivity("ssl1", resolution=30)


if 1 in activated_components:
  simo.add_component_name("ssl2")
  simo.autobuild_model("ssl2",'../inputs/fit_gtfs_3.pdb',"C",resolutions=[1,30],resrange=(1,538), missingbeadsize=30,color=0.0,attachbeads=True)
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
        print "MultipleStates: rigid bodies were not intialized"
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


prot=simo.prot
outputobjects.append(simo)
sampleobjects.append(simo)


##simo.draw_hierarchy_composition()

#####################################################
#restraint setup
#####################################################

#EV restraint
ev=restraints.ExcludedVolumeSphere(prot, resolution=30)
ev.add_to_model()
outputobjects.append(ev)


#cross-link restraint
xl=restraints.SigmoidCrossLinkMS(prot,\
                             '../inputs/Ranish_Kornberg_thiih_xlinks.txt', \
                             inflection=25.0,slope=5.0,amplitude=25.0,linear_slope=0.05,resolution=1)
xl.add_to_model()
outputobjects.append(xl)
#xl.plot_restraint(5.,5.,maxdist=100,npoints=20)


#GaussianEMRestraint
Rad3 = [0,3,5,12,26,28,35,38,42,48,58,65,67,68,70,79,108,113,118]
Kinase = [2,4,11,16,21,22,23,24,32,33,40,43,47,49,54,64,66,73,74,75,81,82,86,85,89,92,95,98,100,106,109,116]
Ssl2 = [1,7,10,15,39,50,52,57,60,63,76,77,78,80,83,91,103,104,110,119]
Tfiihcore = list(set(range(120))-set(Rad3)-set(Kinase)-set(Ssl2))

if 1 in activated_components:
  em1=restraints.GaussianEMRestraint(prot,'../inputs/fixed_holoIIH.mrc.txt',segment_anchors=Kinase,
                                segment_parts=['tfb3','ccl1','kin28'],resolution=30)
  em1.set_label('kinase_em')
  em1.add_to_model()
  sampleobjects.append(em1)
  outputobjects.append(em1)

if 2 in activated_components:
  em2=restraints.GaussianEMRestraint(prot,'../inputs/fixed_holoIIH.mrc.txt',segment_anchors=Rad3,segment_parts=['rad3'],resolution=30)
  em2.set_label('rad3_em')
  em2.add_to_model()
  sampleobjects.append(em2)
  outputobjects.append(em2)

if 3 in activated_components:
  em3=restraints.GaussianEMRestraint(prot,'../inputs/fixed_holoIIH.mrc.txt',segment_anchors=Ssl2,segment_parts=['ssl2'],resolution=30)
  em3.set_label('ssl2_em')
  em3.add_to_model()
  sampleobjects.append(em3)
  outputobjects.append(em3)

if 4 in activated_components:
  em4=restraints.GaussianEMRestraint(prot,'../inputs/fixed_holoIIH.mrc.txt',segment_anchors=Tfiihcore,
                 segment_parts=['tfb1','tfb2','tfb4','tfb5','ssl1'],resolution=30)
  em4.set_label('tfiihcore_em')
  em4.add_to_model()
  sampleobjects.append(em4)
  outputobjects.append(em4)




####################################################
# Monte Carlo
####################################################

mc=samplers.MonteCarlo(m,sampleobjects,5.)
#mc.set_simulated_annealing(15.,45.,50,10)
mc.set_simulated_annealing(25.,75.,50,10)
#cg=samplers.ConjugateGradients(m,sampleobjects)
outputobjects.append(mc)
#outputobjects.append(cg)


import sys
sw = tools.Stopwatch()
outputobjects.append(sw)

output = output.Output()
output.init_stat2("stat.dat", outputobjects,extralabels=["rmf_file","rmf_frame_index"])


#####################################################
#running simulation
#####################################################
nrmffiles=1
nframes=20#000
bestscore, step = 1000000000000000, 0
for k in range(nrmffiles):
  rmffile="models.rmf"
  output.init_rmf(rmffile, prot)
  output.add_restraints_to_rmf(rmffile,[xl])

  for i in range(nframes):
    mc.run(ncycl)
    print mc.get_frame_number()
    score = m.evaluate(False)

    output.set_output_entry("rmf_file",rmffile)
    output.set_output_entry("rmf_frame_index",step)
    if score < bestscore:
      output.write_stats2()
      print '\tBest score: ',score
      output.write_rmf(rmffile,step)
      bestscore = score
      step += 1
  output.close_rmf(rmffile)



