#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container
#import IMP.isd2
import random
import os
import IMP.rmf
import RMF
import inspect

from numpy import random
#from numpy import *

from math import cos
from math import sqrt

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

'''
#####################################################
#read EM maps and centroid coordinates
#####################################################
map_kinase = IMP.em.read_map('kinase.mrc', IMP.em.MRCReaderWriter())
map_rad3 = IMP.em.read_map('rad3.mrc', IMP.em.MRCReaderWriter())
map_ssl2 = IMP.em.read_map('ssl2.mrc', IMP.em.MRCReaderWriter())
map_tfb = IMP.em.read_map('tfb.mrc', IMP.em.MRCReaderWriter())
'''



#####################################################
#create hierarchies and rigid bodies and flexible parts
#####################################################

m=IMP.Model()
simo=representation.Representation(m, upperharmonic=False, disorderedlength=True)



activated_components=[1,2,3,4]


if 1 in activated_components:
  
  simo.create_component("ccl1")
  simo.autobuild_model("ccl1",'../inputs/CCL1_48-393.pdb',"A",resolutions=[1,30],resrange=(1,393),missingbeadsize=30,color=1.0,attachbeads=True)
  simo.setup_component_sequence_connectivity("ccl1", resolution=30)

  simo.create_component("kin28")
  simo.autobuild_model("kin28",'../inputs/KIN28_5-299.pdb',"A",resolutions=[1,30],resrange=(1,306), missingbeadsize=30,color=0.9,attachbeads=True)
  simo.setup_component_sequence_connectivity("kin28", resolution=30)

  simo.create_component("tfb3")
  simo.autobuild_model("tfb3",'../inputs/TFB3_8-142.pdb',"A",resolutions=[1,30],resrange=(1,321), missingbeadsize=30,color=0.8,attachbeads=True)
  simo.setup_component_sequence_connectivity("tfb3", resolution=30)


if 2 in activated_components:

  simo.create_component("rad3")
  simo.autobuild_model("rad3",'../inputs/RAD3_14-725.pdb',"A",resolutions=[1,30],resrange=(1,778), missingbeadsize=30,color=150./360,attachbeads=True)
  simo.setup_component_sequence_connectivity("rad3", resolution=30)


if 4 in activated_components:

  simo.create_component("tfb1")
  simo.autobuild_model("tfb1",'../inputs/TFB1_2-115.pdb',"A",resolutions=[1,30],resrange=(1,642), missingbeadsize=30,color=60./360,attachbeads=True)
  simo.setup_component_sequence_connectivity("tfb1", resolution=30)

  simo.create_component("tfb2")
  simo.autobuild_model("tfb2",'../inputs/TFB2_1-168.pdb',"A",resolutions=[1,30],resrange=(1,170), missingbeadsize=30,color=185./360,attachbeads=True)
  simo.autobuild_model("tfb2",'../inputs/TFB2_186-417.pdb',"A",resolutions=[1,30],resrange=(171,417), missingbeadsize=30,color=185./360,attachbeads=True)
  simo.autobuild_model("tfb2",'../inputs/TFB2_392-513.pdb',"A",resolutions=[1,30],resrange=(418,513), missingbeadsize=30,color=185./360,attachbeads=True)
  simo.setup_component_sequence_connectivity("tfb2", resolution=30)

  simo.create_component("tfb4")
  simo.autobuild_model("tfb4",'../inputs/TFB4_24-250.pdb',"A",resolutions=[1,30],resrange=(1,338), missingbeadsize=30,color=210./360,attachbeads=True)
  simo.setup_component_sequence_connectivity("tfb4", resolution=30)

  simo.create_component("tfb5")
  simo.autobuild_model("tfb5",'../inputs/TFB5.pdb',"B",resolutions=[1,30],resrange=(1,72), missingbeadsize=30,color=0.,attachbeads=True)
  simo.setup_component_sequence_connectivity("tfb5", resolution=30)

  simo.create_component("ssl1")
  simo.autobuild_model("ssl1",'../inputs/SSL1_123-302.pdb',"A",resolutions=[1,30],resrange=(1,302), missingbeadsize=30,color=285./360,attachbeads=True)
  simo.autobuild_model("ssl1",'../inputs/SSL1_386-455.pdb',"A",resolutions=[1,30],resrange=(303,461), missingbeadsize=30,color=285./360,attachbeads=True)
  simo.setup_component_sequence_connectivity("ssl1", resolution=30)


if 3 in activated_components:
  simo.create_component("ssl2")
  #simo.autobuild_model("ssl2",'fit_gtfs_3.pdb',"C",resolutions=[1,30],resrange=(1,538), missingbeadsize=30,color=0.0,attachbeads=True)
  #simo.autobuild_model("ssl2",'fit_gtfs_3.pdb',"D",resolutions=[1,30],resrange=(539,843), missingbeadsize=30,color=0.0,attachbeads=True)
  simo.autobuild_model("ssl2",'../inputs/SSL2_316-723.pdb',"A",resolutions=[1,30],resrange=(1,843), missingbeadsize=30,color=0.0,attachbeads=True)
  simo.setup_component_sequence_connectivity("ssl2", resolution=30)






# --- set rigid bodies
if 1 in activated_components:
  emxk,emyk,emzk = -68.03536154, 4.11344154, 33.2057 #99
  simo.set_rigid_bodies(["tfb3"],(emxk,emyk,emzk))
  simo.set_rigid_bodies(["ccl1","kin28"],(emxk,emyk,emzk))

if 2 in activated_components:
  emxr,emyr,emzr = -36.77662667, -44.77829333,  12.750732 #73
  simo.set_rigid_bodies(["rad3"],(emxr,emyr,emzr))

if 4 in activated_components:
  emxt,emyt,emzt = 8.31117245, -63.22024737,  54.85138763 #40
  simo.set_rigid_bodies([("tfb2",(437,513)),("tfb5",(1,72))], (emxt,emyt,emzt))
  #simo.set_rigid_bodies([("tfb2",(1,436))], (emxt,emyt,emzt))
  simo.set_rigid_bodies([("tfb2",(1,168))], (emxt,emyt,emzt))
  simo.set_rigid_bodies([("tfb2",(169,436))], (emxt,emyt,emzt))
  simo.set_rigid_bodies(["tfb1"], (emxt,emyt,emzt))
  simo.set_rigid_bodies(["tfb4"],(emxt,emyt,emzt))
  emxt,emyt,emzt = 8.31117245, -63.22024737,  54.85138763 #40
  simo.set_rigid_bodies([("ssl1",(1,385))],(emxt,emyt,emzt))
  simo.set_rigid_bodies([("ssl1",(386,461))],(emxt,emyt,emzt))


if 3 in activated_components:
  emxk,emyk,emzk = -23.7230975, -13.095193, 69.7855035
  simo.set_rigid_bodies(["ssl2"],(emxk,emyk,emzk))


simo.set_floppy_bodies()
d=simo.get_particles_to_sample()
simo.set_rigid_bodies_max_trans(rbmaxtrans)
simo.set_floppy_bodies_max_trans(fbmaxtrans)


#re-orient initial positions
simo.shuffle_configuration(translate=False)


prot=simo.get_hierarchy()
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
Polii = [0, 4, 5, 11, 12, 16, 18, 21, 22, 23, 24, 26, 30, 34, 35, 37, 41, 45, \
        46, 47, 48, 49, 50, 57, 58, 60, 61, 63, 68, 74, 75, 77, 78, 79, 83, 84,\
        90, 91, 96, 103, 106, 108, 114, 115, 118, 119, 121, 123, 124, 125, 130,\
        131, 133, 139, 142, 143, 144, 145, 150, 155, 157, 158, 163, 166, 168, \
        169, 172, 173, 174, 176, 177, 179, 181, 185, 186, 187, 190, 192, 195, \
        196, 201, 202, 207, 208, 211, 213, 218, 220, 222, 225, 227, 228, 229, \
        231, 235, 240, 241, 242, 243, 244, 250, 251, 254, 256, 261, 262, 264, \
        267, 276, 279, 280, 285, 287, 288, 290, 293, 294, 298]
Rad3 = [14, 20, 29, 56, 62, 73, 153, 156, 162, 171, 214, 239, 291, 296, 297]
Tfiif3 = [2, 13, 59, 70, 82, 93, 102, 110, 122, 160, 161, 164, 165, 191, 193, \
         197, 199, 200, 212, 219, 223, 226, 232, 257, 268, 271, 274, 278, 282, 283, 284, 289]
Kinase = [32, 64, 67, 89, 92, 99, 100, 129, 137, 167, 182, 233, 277]
Ssl2 = [14, 17, 19, 25, 27, 39, 43, 51, 53, 69, 85, 94, 95, 98, 113, 138, 140, \
        149, 152, 154, 180, 206, 246, 252, 253, 269, 273, 286]
Rest = [8, 9, 10, 36, 54, 55, 65, 66, 76, 80, 81, 87, 107, 112, 116, 127, 136, \
        148, 175, 183, 189, 209, 221, 224, 245, 249, 255, 265, 266, 292]
Tfiihcore = [2, 6, 28, 33, 38, 40, 42, 43, 44, 72, 86, 88, 104, 109, 117, 120, \
            126, 134, 146, 147, 170, 203, 204, 205, 210, 216, 217, 234, 238, 247,\
            248, 253, 258, 260, 272, 275, 281, 299]
Tfiie = [1, 3, 15, 31, 52, 71, 97, 101, 105, 128, 132, 135, 141, 147, 151, 159, 178,\
         184, 188, 194, 197, 198, 215, 217, 230, 237, 259, 263, 270, 275, 295]

if 1 in activated_components:
  em1=restraints.GaussianEMRestraint(prot,'../inputs/themap.txt',segment_anchors=Kinase,
                                segment_parts=['tfb3','ccl1','kin28'],resolution=30)
  em1.set_label('kinase_em')
  em1.add_to_model()
  sampleobjects.append(em1)
  outputobjects.append(em1)

if 2 in activated_components:
  em2=restraints.GaussianEMRestraint(prot,'../inputs/themap.txt',segment_anchors=Rad3,segment_parts=['rad3'],resolution=30)
  em2.set_label('rad3_em')
  em2.add_to_model()
  sampleobjects.append(em2)
  outputobjects.append(em2)

if 3 in activated_components:
  em3=restraints.GaussianEMRestraint(prot,'../inputs/themap.txt',segment_anchors=Ssl2,segment_parts=['ssl2'],resolution=30)
  em3.set_label('ssl2_em')
  em3.add_to_model()
  sampleobjects.append(em3)
  outputobjects.append(em3)

if 4 in activated_components:
  em4=restraints.GaussianEMRestraint(prot,'../inputs/themap.txt',segment_anchors=Tfiihcore,
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



