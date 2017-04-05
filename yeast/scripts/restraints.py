from __future__ import print_function, absolute_import
import IMP.isd
import IMP.atom
import IMP.pmi.tools as tools

def get_particles_by_resolution(prot,resolution):
    #this function does not work with the root hierarchy, but
    #individual proteins
    #for hier in prot.get_children():
    particles = []
    residues = set()

    for p in IMP.atom.get_leaves(prot):
        residues.update(IMP.atom.Fragment(p).get_residue_indexes())

    firstresn=min(residues)
    lastresn=max(residues)
    for nres in range(firstresn,lastresn+1):
        s=IMP.atom.Selection(prot,residue_index=nres)
        resolutions=[]

        # calculate the closest resolution for each set of particles that represent a residue
        ps=s.get_selected_particles()

        if len(ps)>0:
            for p in ps:
                resolutions.append(IMP.pmi.Resolution(IMP.pmi.Resolution(p)).get_resolution())
            closestres=min(resolutions, key=lambda x:abs(float(x)-float(resolution)))

            # now we get the particle
            for p in ps:
                if closestres==IMP.pmi.Resolution.get_resolution(IMP.pmi.Resolution(p)):
                    if not p in particles:
                        particles.append(p)
        else:
            print("get_particles_by_resolution> WARNING residue %d in molecule %s is missing" % (nres,prot.get_name()))

    return list(particles)

class GaussianEMRestraint():
    """Old version of the GaussianEMRestraint (from salilab/pmi@c522e0abc),
       updated to work with current IMP."""

    def __init__(self,prot,map_anchors_fn,segment_anchors=None,segment_parts=None,rigid=True,resolution=None):
        #segment parts should be a list containing protein names or a tuple with protein names and residue ranges:
        # [("ABC",1,100),"CYT","CDR"]

        #import IMP.multifit

        if segment_anchors==None: segment_anchors=[]
        if segment_parts==None: segment_parts=[]

        #dcoords=IMP.multifit.read_anchors_data(map_anchors_fn).points_
        self.prot=prot
        self.m=self.prot.get_model()



        if resolution!=None:
            hierarchy_anchors=[]
            for prot in self.prot.get_children():
                hierarchy_anchors+=get_particles_by_resolution(prot,resolution)


        if segment_parts!=None:
            model_anchors=[]
            for seg in segment_parts:
                if type(seg)==str:
                    s=IMP.atom.Selection(self.prot,molecule=seg)
                    ps=s.get_selected_particles()
                    model_anchors+=ps
                elif type(seg)==tuple:
                    s=IMP.atom.Selection(self.prot,molecule=seg[0],residue_indexes=range(seg[1],seg[2]+1))
                    ps=s.get_selected_particles()
                    model_anchors+=ps
            if resolution!=None:
                #get the intersection to remove redundant particles
                self.model_anchors=(list(set(model_anchors) & set(hierarchy_anchors)))

        for p in self.model_anchors:
            print(p.get_name())



        data = open(map_anchors_fn)
        D = data.readlines()
        data.close()
        dcoords={}
        for d in D:
            d=d.strip().split('|')
            if len(d)==6: dcoords[int(d[1])] = IMP.algebra.Vector3D(float(d[2]),float(d[3]),float(d[4]))

        # parameters
        self.model_sigmas=[15.0]*len(self.model_anchors)

        self.model_sigmas=[]
        self.model_weights=[]
        for p in self.model_anchors: self.model_sigmas.append(IMP.core.XYZR(p).get_radius())
        for p in self.model_weights: self.model_weights.append(len(IMP.atom.Fragment(p).get_residue_indexes()))



        #self.model_sigmas=[float(anch.get_as_xyzr().get_radius()) for anch in self.segment_parts]
        self.model_weights=[1.0]*len(self.model_anchors)
        self.density_sigmas=[15.0]*len(segment_anchors)
        self.density_weights=[5.0]*len(segment_anchors)
        self.sigmamaxtrans=0.1
        self.sigmamin=1.
        self.sigmamax=100.0
        self.sigmainit=10.0
        self.cutoff_dist_for_container=10.0
        self.rigid=rigid
        self.segment_anchors=segment_anchors
        self.segment_parts=segment_parts
        self.tabexp=True


        self.density_anchors=[]
        for d in dcoords:
            if self.segment_anchors==[]:
                p=IMP.Particle(self.m)
                self.density_anchors.append(p)
                IMP.core.XYZR.setup_particle(p,\
                                         IMP.algebra.Sphere3D(d,\
                                         self.density_sigmas[nd]*1.5))
            else:
                if d in self.segment_anchors:
                    p=IMP.Particle(self.m)
                    self.density_anchors.append(p)
                    IMP.core.XYZR.setup_particle(p,\
                                         IMP.algebra.Sphere3D(dcoords[d],\
                                         self.density_sigmas[0]*1.5))

        for np,p in enumerate(self.model_anchors):
            self.model_sigmas[np]=IMP.core.XYZR(p).get_radius()/1.5

        self.sigmaglobal=tools.SetupNuisance(self.m,self.sigmainit,
                 self.sigmamin,self.sigmamax,True).get_particle()
        print('setting up restraint')

        self.gaussianEM_restraint=IMP.isd.GaussianAnchorEMRestraint(
            self.model_anchors,self.model_sigmas,self.model_weights,
            self.density_anchors,self.density_sigmas,self.density_weights,
            self.sigmaglobal.get_particle(),self.cutoff_dist_for_container,
            self.rigid,self.tabexp)
        print('done setup')
        self.rs = IMP.RestraintSet(self.m, 'GaussianEMRestraint')
        self.rs.add_restraint(self.gaussianEM_restraint)

    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        IMP.pmi.tools.add_restraint_to_model(self.m, self.rs)

    def get_particles_to_sample(self):
        ps={}
        ps["Nuisances_GaussianEMRestraint_sigma_"+self.label]=([self.sigmaglobal],self.sigmamaxtrans)
        return ps

    def get_hierarchy(self):
        return self.prot

    def get_restraint_set(self):
        return self.rs

    def get_output(self):
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)
        output["GaussianEMRestraint_"+self.label]=str(self.rs.unprotected_evaluate(None))
        output["GaussianEMRestraint_sigma_"+self.label]=str(self.sigmaglobal.get_scale())
        return output
