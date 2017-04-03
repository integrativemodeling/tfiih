#!/usr/bin/env python
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import os
import IMP.rmf
import RMF
import restraints
import topology
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.representation as representation
import IMP.pmi.tools as tools
import IMP.pmi.samplers as samplers
import IMP.pmi.output as output

run_number = 0 # Change for each independent run
ncycl=100  #number of montecarlo steps cycles
outputobjects=[]
sampleobjects=[]

m, simo = topology.make_topology()

prot=simo.prot
outputobjects.append(simo)
sampleobjects.append(simo)


##simo.draw_hierarchy_composition()

#####################################################
#restraint setup
#####################################################

#EV restraint
ev=IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo, resolution=30)
ev.add_to_model()
outputobjects.append(ev)


#cross-link restraint
xl=IMP.pmi.restraints.crosslinking.SigmoidalCrossLinkMS(simo,
                             '../inputs/Ranish_Kornberg_thiih_xlinks.txt',
                             inflection=25.0,slope=5.0,amplitude=25.0,linear_slope=0.05,resolution=1)
xl.add_to_model()
outputobjects.append(xl)
#xl.plot_restraint(5.,5.,maxdist=100,npoints=20)


#GaussianEMRestraint
Rad3 = [0,3,5,12,26,28,35,38,42,48,58,65,67,68,70,79,108,113,118]
Kinase = [2,4,11,16,21,22,23,24,32,33,40,43,47,49,54,64,66,73,74,75,81,82,86,85,89,92,95,98,100,106,109,116]
Ssl2 = [1,7,10,15,39,50,52,57,60,63,76,77,78,80,83,91,103,104,110,119]
Tfiihcore = list(set(range(120))-set(Rad3)-set(Kinase)-set(Ssl2))

if 1 in topology.activated_components:
  em1=restraints.GaussianEMRestraint(prot,'../inputs/fixed_holoIIH.mrc.txt',segment_anchors=Kinase,
                                segment_parts=['tfb3','ccl1','kin28'],resolution=30)
  em1.set_label('kinase_em')
  em1.add_to_model()
  sampleobjects.append(em1)
  outputobjects.append(em1)

if 2 in topology.activated_components:
  em2=restraints.GaussianEMRestraint(prot,'../inputs/fixed_holoIIH.mrc.txt',segment_anchors=Rad3,segment_parts=['rad3'],resolution=30)
  em2.set_label('rad3_em')
  em2.add_to_model()
  sampleobjects.append(em2)
  outputobjects.append(em2)

if 3 in topology.activated_components:
  em3=restraints.GaussianEMRestraint(prot,'../inputs/fixed_holoIIH.mrc.txt',segment_anchors=Ssl2,segment_parts=['ssl2'],resolution=30)
  em3.set_label('ssl2_em')
  em3.add_to_model()
  sampleobjects.append(em3)
  outputobjects.append(em3)

if 4 in topology.activated_components:
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
output.init_stat2("../outputs/stat_%d.dat" % run_number,
                  outputobjects,extralabels=["rmf_file","rmf_frame_index"])


#####################################################
#running simulation
#####################################################
nrmffiles=1
nframes=20 if '--test' in sys.argv else 20000
bestscore, step = 1000000000000000, 0
for k in range(nrmffiles):
  rmffile="../outputs/models_%d.0.rmf" % run_number
  output.init_rmf(rmffile, [prot])
  output.add_restraints_to_rmf(rmffile,[xl])
  rset = IMP.pmi.tools.get_restraint_set(m)

  for i in range(nframes):
    mc.optimize(ncycl)
    print mc.get_frame_number()
    m.update()
    score = rset.evaluate(False)

    output.set_output_entry("rmf_file",rmffile)
    output.set_output_entry("rmf_frame_index",step)
    if score < bestscore:
      output.write_stats2()
      print '\tBest score: ',score
      output.write_rmf(rmffile)
      bestscore = score
      step += 1
  output.close_rmf(rmffile)



