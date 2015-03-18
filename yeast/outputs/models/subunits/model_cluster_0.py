import cPickle, base64
try:
	from SimpleSession.versions.v62 import beginRestore,\
	    registerAfterModelsCB, reportRestoreError, checkVersion
except ImportError:
	from chimera import UserError
	raise UserError('Cannot open session that was saved in a'
	    ' newer version of Chimera; update your version')
checkVersion([1, 9, 39798])
import chimera
from chimera import replyobj
replyobj.status('Restoring session...', \
    blankAfter=0)
replyobj.status('Beginning session restore...', \
    blankAfter=0, secondary=True)
beginRestore()

def restoreCoreModels():
	from SimpleSession.versions.v62 import init, restoreViewer, \
	     restoreMolecules, restoreColors, restoreSurfaces, \
	     restoreVRML, restorePseudoBondGroups, restoreModelAssociations
	molInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVRFyaWJib25JbnNpZGVDb2xvcnECSwBOfYdVCWJhbGxTY2FsZXEDSwBOfYdVCXBvaW50U2l6ZXEESwBOfYdVBWNvbG9ycQVLAE59h1UKcmliYm9uVHlwZXEGSwBOfYdVCnN0aWNrU2NhbGVxB0sATn2HVQxhcm9tYXRpY01vZGVxCEsATn2HVQp2ZHdEZW5zaXR5cQlLAE59h1UGaGlkZGVucQpLAE59h1UNYXJvbWF0aWNDb2xvcnELSwBOfYdVD3JpYmJvblNtb290aGluZ3EMSwBOfYdVCWF1dG9jaGFpbnENSwBOfYdVCG9wdGlvbmFscQ59VQ9sb3dlckNhc2VDaGFpbnNxD0sATn2HVQlsaW5lV2lkdGhxEEsATn2HVQ9yZXNpZHVlTGFiZWxQb3NxEUsATn2HVQRuYW1lcRJLAE59h1UPYXJvbWF0aWNEaXNwbGF5cRNLAE59h1UPcmliYm9uU3RpZmZuZXNzcRRLAE59h1UKcGRiSGVhZGVyc3EVXVUDaWRzcRZLAE59h1UOc3VyZmFjZU9wYWNpdHlxF0sATn2HVRBhcm9tYXRpY0xpbmVUeXBlcRhLAE59h1UUcmliYm9uSGlkZXNNYWluY2hhaW5xGUsATn2HVQdkaXNwbGF5cRpLAE59h3Uu'))
	resInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQZpbnNlcnRxAksATn2HVQtmaWxsRGlzcGxheXEDSwBOfYdVBG5hbWVxBEsATn2HVQVjaGFpbnEFSwBOfYdVDnJpYmJvbkRyYXdNb2RlcQZLAE59h1UCc3NxB0sATn2HVQhtb2xlY3VsZXEISwBOfYdVC3JpYmJvbkNvbG9ycQlLAE59h1UFbGFiZWxxCksATn2HVQpsYWJlbENvbG9ycQtLAE59h1UIZmlsbE1vZGVxDEsATn2HVQVpc0hldHENSwBOfYdVC2xhYmVsT2Zmc2V0cQ5LAE59h1UIcG9zaXRpb25xD11VDXJpYmJvbkRpc3BsYXlxEEsATn2HVQhvcHRpb25hbHERfVUEc3NJZHESSwBOfYd1Lg=='))
	atomInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQdyZXNpZHVlcQJLAE59h1UIdmR3Q29sb3JxA0sATn2HVQRuYW1lcQRLAE59h1UDdmR3cQVLAE59h1UOc3VyZmFjZURpc3BsYXlxBksATn2HVQVjb2xvcnEHSwBOfYdVCWlkYXRtVHlwZXEISwBOfYdVBmFsdExvY3EJSwBOfYdVBWxhYmVscQpLAE59h1UOc3VyZmFjZU9wYWNpdHlxC0sATn2HVQdlbGVtZW50cQxLAE59h1UKbGFiZWxDb2xvcnENSwBOfYdVDHN1cmZhY2VDb2xvcnEOSwBOfYdVD3N1cmZhY2VDYXRlZ29yeXEPSwBOfYdVBnJhZGl1c3EQSwBOfYdVCmNvb3JkSW5kZXhxEV1VC2xhYmVsT2Zmc2V0cRJLAE59h1USbWluaW11bUxhYmVsUmFkaXVzcRNLAE59h1UIZHJhd01vZGVxFEsATn2HVQhvcHRpb25hbHEVfXEWVQxzZXJpYWxOdW1iZXJxF4iIXYdzVQdkaXNwbGF5cRhLAE59h3Uu'))
	bondInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQVjb2xvcnECSwBOfYdVBWF0b21zcQNdVQVsYWJlbHEESwBOfYdVCGhhbGZib25kcQVLAE59h1UGcmFkaXVzcQZLAE59h1ULbGFiZWxPZmZzZXRxB0sATn2HVQhkcmF3TW9kZXEISwBOfYdVCG9wdGlvbmFscQl9VQdkaXNwbGF5cQpLAE59h3Uu'))
	crdInfo = cPickle.loads(base64.b64decode('gAJ9cQEu'))
	surfInfo = {'category': (0, None, {}), 'probeRadius': (0, None, {}), 'pointSize': (0, None, {}), 'name': [], 'density': (0, None, {}), 'colorMode': (0, None, {}), 'useLighting': (0, None, {}), 'transparencyBlendMode': (0, None, {}), 'molecule': [], 'smoothLines': (0, None, {}), 'lineWidth': (0, None, {}), 'allComponents': (0, None, {}), 'twoSidedLighting': (0, None, {}), 'customVisibility': [], 'drawMode': (0, None, {}), 'display': (0, None, {}), 'customColors': []}
	vrmlInfo = {'subid': (0, None, {}), 'display': (0, None, {}), 'id': (0, None, {}), 'vrmlString': [], 'name': (0, None, {})}
	colors = {'Ru': ((0.141176, 0.560784, 0.560784), 1, u'default'), 'Re': ((0.14902, 0.490196, 0.670588), 1, u'default'), 'Rf': ((0.8, 0, 0.34902), 1, u'default'), 'Ra': ((0, 0.490196, 0), 1, u'default'), 'slate gray': ((0.439216, 0.501961, 0.564706), 1, u'default'), 'Rn': ((0.258824, 0.509804, 0.588235), 1, u'default'), 'Rh': ((0.0392157, 0.490196, 0.54902), 1, u'default'), 'Be': ((0.760784, 1, 0), 1, u'default'), 'Ba': ((0, 0.788235, 0), 1, u'default'), 'Bh': ((0.878431, 0, 0.219608), 1, u'default'), 'Bi': ((0.619608, 0.309804, 0.709804), 1, u'default'), 'Bk': ((0.541176, 0.309804, 0.890196), 1, u'default'), 'Br': ((0.65098, 0.160784, 0.160784), 1, u'default'), 'orange': ((1, 0.498039, 0), 1, u'default'), 'H': ((1, 1, 1), 1, u'default'), 'P': ((1, 0.501961, 0), 1, u'default'), 'Os': ((0.14902, 0.4, 0.588235), 1, u'default'), 'orange red': ((1, 0.270588, 0), 1, u'default'), 'Ge': ((0.4, 0.560784, 0.560784), 1, u'default'), 'Gd': ((0.270588, 1, 0.780392), 1, u'default'), 'Ga': ((0.760784, 0.560784, 0.560784), 1, u'default'), 'Pr': ((0.85098, 1, 0.780392), 1, u'default'),
'Rb': ((0.439216, 0.180392, 0.690196), 1, u'default'), 'Pt': ((0.815686, 0.815686, 0.878431), 1, u'default'), 'Pu': ((0, 0.419608, 1), 1, u'default'), 'C': ((0.564706, 0.564706, 0.564706), 1, u'default'), 'Pb': ((0.341176, 0.34902, 0.380392), 1, u'default'), 'forest green': ((0.133333, 0.545098, 0.133333), 1, u'default'), 'Pa': ((0, 0.631373, 1), 1, u'default'), 'Pd': ((0, 0.411765, 0.521569), 1, u'default'), 'Cd': ((1, 0.85098, 0.560784), 1, u'default'), 'Po': ((0.670588, 0.360784, 0), 1, u'default'), 'Pm': ((0.639216, 1, 0.780392), 1, u'default'), 'purple': ((0.627451, 0.12549, 0.941176), 1, u'default'), 'Hs': ((0.901961, 0, 0.180392), 1, u'default'), 'Ho': ((0, 1, 0.611765), 1, u'default'), 'Hf': ((0.301961, 0.760784, 1), 1, u'default'), 'Hg': ((0.721569, 0.721569, 0.815686), 1, u'default'), 'He': ((0.85098, 1, 1), 1, u'default'), 'Md': ((0.701961, 0.0509804, 0.65098), 1, u'default'), 'Mg': ((0.541176, 1, 0), 1, u'default'), 'dim gray': ((0.411765, 0.411765, 0.411765), 1, u'default'), 'K': ((0.560784, 0.25098, 0.831373), 1, u'default'), 'Mn': ((0.611765, 0.478431, 0.780392), 1, u'default'),
'O': ((1, 0.0509804, 0.0509804), 1, u'default'), 'Mt': ((0.921569, 0, 0.14902), 1, u'default'), 'S': ((1, 1, 0.188235), 1, u'default'), 'W': ((0.129412, 0.580392, 0.839216), 1, u'default'), 'cornflower blue': ((0.392157, 0.584314, 0.929412), 1, u'default'), 'Eu': ((0.380392, 1, 0.780392), 1, u'default'), 'Zr': ((0.580392, 0.878431, 0.878431), 1, u'default'), 'Er': ((0, 0.901961, 0.458824), 1, u'default'), 'Ni': ((0.313725, 0.815686, 0.313725), 1, u'default'), 'No': ((0.741176, 0.0509804, 0.529412), 1, u'default'), 'Na': ((0.670588, 0.360784, 0.94902), 1, u'default'), 'Nb': ((0.45098, 0.760784, 0.788235), 1, u'default'), 'Nd': ((0.780392, 1, 0.780392), 1, u'default'), 'Ne': ((0.701961, 0.890196, 0.960784), 1, u'default'), 'Np': ((0, 0.501961, 1), 1, u'default'), 'Fr': ((0.258824, 0, 0.4), 1, u'default'), 'Fe': ((0.878431, 0.4, 0.2), 1, u'default'), 'Fm': ((0.701961, 0.121569, 0.729412), 1, u'default'), 'B': ((1, 0.709804, 0.709804), 1, u'default'), 'F': ((0.564706, 0.878431, 0.313725), 1, u'default'), 'Sr': ((0, 1, 0), 1, u'default'), 'Zn': ((0.490196, 0.501961, 0.690196), 1, u'default'),
'N': ((0.188235, 0.313725, 0.972549), 1, u'default'), 'Kr': ((0.360784, 0.721569, 0.819608), 1, u'default'), 'Si': ((0.941176, 0.784314, 0.627451), 1, u'default'), 'Sn': ((0.4, 0.501961, 0.501961), 1, u'default'), 'Sm': ((0.560784, 1, 0.780392), 1, u'default'), 'V': ((0.65098, 0.65098, 0.670588), 1, u'default'), 'Sc': ((0.901961, 0.901961, 0.901961), 1, u'default'), 'Sb': ((0.619608, 0.388235, 0.709804), 1, u'default'), 'Sg': ((0.85098, 0, 0.270588), 1, u'default'), 'Se': ((1, 0.631373, 0), 1, u'default'), 'Co': ((0.941176, 0.564706, 0.627451), 1, u'default'), 'Cm': ((0.470588, 0.360784, 0.890196), 1, u'default'), 'Cl': ((0.121569, 0.941176, 0.121569), 1, u'default'), 'Ca': ((0.239216, 1, 0), 1, u'default'), 'Cf': ((0.631373, 0.211765, 0.831373), 1, u'default'), 'Ce': ((1, 1, 0.780392), 1, u'default'), 'Xe': ((0.258824, 0.619608, 0.690196), 1, u'default'), 'Tm': ((0, 0.831373, 0.321569), 1, u'default'), 'Cs': ((0.341176, 0.0901961, 0.560784), 1, u'default'), 'Cr': ((0.541176, 0.6, 0.780392), 1, u'default'), 'Cu': ((0.784314, 0.501961, 0.2), 1, u'default'),
'La': ((0.439216, 0.831373, 1), 1, u'default'), 'Li': ((0.8, 0.501961, 1), 1, u'default'), 'Tl': ((0.65098, 0.329412, 0.301961), 1, u'default'), 'Lu': ((0, 0.670588, 0.141176), 1, u'default'), 'Lr': ((0.780392, 0, 0.4), 1, u'default'), 'Th': ((0, 0.729412, 1), 1, u'default'), 'Ti': ((0.74902, 0.760784, 0.780392), 1, u'default'), 'Te': ((0.831373, 0.478431, 0), 1, u'default'), 'Tb': ((0.188235, 1, 0.780392), 1, u'default'), 'Tc': ((0.231373, 0.619608, 0.619608), 1, u'default'), 'Ta': ((0.301961, 0.65098, 1), 1, u'default'), 'pink': ((1, 0.752941, 0.796078), 1, u'default'), 'Yb': ((0, 0.74902, 0.219608), 1, u'default'), 'Db': ((0.819608, 0, 0.309804), 1, u'default'), 'Dy': ((0.121569, 1, 0.780392), 1, u'default'), 'At': ((0.458824, 0.309804, 0.270588), 1, u'default'), 'I': ((0.580392, 0, 0.580392), 1, u'default'), 'U': ((0, 0.560784, 1), 1, u'default'), 'Y': ((0.580392, 1, 1), 1, u'default'), 'Ac': ((0.439216, 0.670588, 0.980392), 1, u'default'), 'Ag': ((0.752941, 0.752941, 0.752941), 1, u'default'), 'Ir': ((0.0901961, 0.329412, 0.529412), 1, u'default'),
'Am': ((0.329412, 0.360784, 0.94902), 1, u'default'), 'Al': ((0.74902, 0.65098, 0.65098), 1, u'default'), 'As': ((0.741176, 0.501961, 0.890196), 1, u'default'), 'Ar': ((0.501961, 0.819608, 0.890196), 1, u'default'), 'Au': ((1, 0.819608, 0.137255), 1, u'default'), 'Es': ((0.701961, 0.121569, 0.831373), 1, u'default'), 'In': ((0.65098, 0.458824, 0.45098), 1, u'default'), 'Mo': ((0.329412, 0.709804, 0.709804), 1, u'default')}
	materials = {u'default': ((0.85, 0.85, 0.85), 30)}
	pbInfo = {'category': [u'distance monitor'], 'bondInfo': [{'color': (0, None, {}), 'atoms': [], 'label': (0, None, {}), 'halfbond': (0, None, {}), 'labelColor': (0, None, {}), 'drawMode': (0, None, {}), 'display': (0, None, {})}], 'lineType': (1, 2, {}), 'color': (1, 0, {}), 'optional': {'fixedLabels': (True, False, (1, False, {}))}, 'display': (1, True, {}), 'showStubBonds': (1, False, {}), 'lineWidth': (1, 1, {}), 'stickScale': (1, 1, {}), 'id': [-2]}
	modelAssociations = {}
	colorInfo = (3, (u'white', (1, 1, 1, 1)), {(u'green', (0, 1, 0, 1)): [2], (u'yellow', (1, 1, 0, 1)): [0]})
	viewerInfo = {'cameraAttrs': {'center': (-84.53450012207, 2.3449478149414, 33.621021311414), 'fieldOfView': 24.735315116823, 'nearFar': (127.85411500303, -60.612072380205), 'ortho': False, 'eyeSeparation': 50.8, 'focal': 36.022327423096}, 'viewerAttrs': {'silhouetteColor': None, 'clipping': False, 'showSilhouette': False, 'showShadows': False, 'viewSize': 18.834988498185, 'labelsOnTop': True, 'depthCueRange': (0.5, 1), 'silhouetteWidth': 2, 'singleLayerTransparency': True, 'shadowTextureSize': 2048, 'backgroundImage': [None, 1, 2, 1, 0, 0], 'backgroundGradient': [('Chimera default', [(1, 1, 1, 1), (0, 0, 1, 1)], 1), 1, 0, 0], 'depthCue': True, 'highlight': 0, 'scaleFactor': 0.19446355467659, 'angleDependentTransparency': True, 'backgroundMethod': 0}, 'viewerHL': 2, 'cameraMode': 'mono', 'detail': 1.5, 'viewerFog': None, 'viewerBG': 1}

	replyobj.status("Initializing session restore...", blankAfter=0,
		secondary=True)
	from SimpleSession.versions.v62 import expandSummary
	init(dict(enumerate(expandSummary(colorInfo))))
	replyobj.status("Restoring colors...", blankAfter=0,
		secondary=True)
	restoreColors(colors, materials)
	replyobj.status("Restoring molecules...", blankAfter=0,
		secondary=True)
	restoreMolecules(molInfo, resInfo, atomInfo, bondInfo, crdInfo)
	replyobj.status("Restoring surfaces...", blankAfter=0,
		secondary=True)
	restoreSurfaces(surfInfo)
	replyobj.status("Restoring VRML models...", blankAfter=0,
		secondary=True)
	restoreVRML(vrmlInfo)
	replyobj.status("Restoring pseudobond groups...", blankAfter=0,
		secondary=True)
	restorePseudoBondGroups(pbInfo)
	replyobj.status("Restoring model associations...", blankAfter=0,
		secondary=True)
	restoreModelAssociations(modelAssociations)
	replyobj.status("Restoring camera...", blankAfter=0,
		secondary=True)
	restoreViewer(viewerInfo)

try:
	restoreCoreModels()
except:
	reportRestoreError("Error restoring core models")

	replyobj.status("Restoring extension info...", blankAfter=0,
		secondary=True)


try:
	import StructMeasure
	from StructMeasure.DistMonitor import restoreDistances
	registerAfterModelsCB(restoreDistances, 1)
except:
	reportRestoreError("Error restoring distances in session")


def restoreMidasBase():
	formattedPositions = {}
	import Midas
	Midas.restoreMidasBase(formattedPositions)
try:
	restoreMidasBase()
except:
	reportRestoreError('Error restoring Midas base state')


def restoreMidasText():
	from Midas import midas_text
	midas_text.aliases = {}
	midas_text.userSurfCategories = {}

try:
	restoreMidasText()
except:
	reportRestoreError('Error restoring Midas text state')

geomData = {'AxisManager': {}, 'CentroidManager': {}, 'PlaneManager': {}}

try:
	from StructMeasure.Geometry import geomManager
	geomManager._restoreSession(geomData)
except:
	reportRestoreError("Error restoring geometry objects in session")


def restoreSession_RibbonStyleEditor():
	import SimpleSession
	import RibbonStyleEditor
	userScalings = [('TPPfront', [[0.5, 0.5], [0.9, 0.25], [0.9, 0.25], [1.8, 0.25, 0.25, 0.25], [0.9, 0.25]]), ('default', [[0.25, 0.25], [1.5, 0.25], [1.5, 0.25], [0.25, 0.25, 3, 0.25], [1.5, 0.25]]), ('licorice', [[0.35, 0.35], [0.35, 0.35], [0.35, 0.35], [0.35, 0.35, 0.35, 0.35], [0.35, 0.35]]), ('irina', [[0.25, 0.2], [1.4, 0.1], [1.3, 0.1], [2.5, 0.15, 0.1, 0.15], [1.5, 0.15]]), ('TPP_Front', [[0.5, 0.5], [1.25, 0.25], [1.25, 0.25], [1.8, 0.25, 0.25, 0.25], [0.9, 0.25]]), ('irina_3', [[0.2, 0.2], [1, 0.2], [1, 0.2], [2, 0.2, 0.1, 0.1], [0.8, 0.2]]), ('irina_2', [[0.2, 0.2], [1.3, 0.2], [1.2, 0.2], [2, 0.2, 0.2, 0.2], [0.9, 0.2]])]
	userXSections = []
	userResidueClasses = []
	residueData = []
	flags = RibbonStyleEditor.NucleicDefault1
	SimpleSession.registerAfterModelsCB(RibbonStyleEditor.restoreState,
				(userScalings, userXSections,
				userResidueClasses, residueData, flags))
try:
	restoreSession_RibbonStyleEditor()
except:
	reportRestoreError("Error restoring RibbonStyleEditor state")

trPickle = 'gAJjQW5pbWF0ZS5UcmFuc2l0aW9ucwpUcmFuc2l0aW9ucwpxASmBcQJ9cQMoVQxjdXN0b21fc2NlbmVxBGNBbmltYXRlLlRyYW5zaXRpb24KVHJhbnNpdGlvbgpxBSmBcQZ9cQcoVQZmcmFtZXNxCEsBVQ1kaXNjcmV0ZUZyYW1lcQlLAVUKcHJvcGVydGllc3EKXXELVQNhbGxxDGFVBG5hbWVxDWgEVQRtb2RlcQ5VBmxpbmVhcnEPdWJVCGtleWZyYW1lcRBoBSmBcRF9cRIoaAhLFGgJSwFoCl1xE2gMYWgNaBBoDmgPdWJVBXNjZW5lcRRoBSmBcRV9cRYoaAhLAWgJSwFoCl1xF2gMYWgNaBRoDmgPdWJ1Yi4='
scPickle = 'gAJjQW5pbWF0ZS5TY2VuZXMKU2NlbmVzCnEBKYFxAn1xA1UHbWFwX2lkc3EEfXNiLg=='
kfPickle = 'gAJjQW5pbWF0ZS5LZXlmcmFtZXMKS2V5ZnJhbWVzCnEBKYFxAn1xA1UHZW50cmllc3EEXXEFc2Iu'
def restoreAnimation():
	'A method to unpickle and restore animation objects'
	# Scenes must be unpickled after restoring transitions, because each
	# scene links to a 'scene' transition. Likewise, keyframes must be 
	# unpickled after restoring scenes, because each keyframe links to a scene.
	# The unpickle process is left to the restore* functions, it's 
	# important that it doesn't happen prior to calling those functions.
	import SimpleSession
	from Animate.Session import restoreTransitions
	from Animate.Session import restoreScenes
	from Animate.Session import restoreKeyframes
	SimpleSession.registerAfterModelsCB(restoreTransitions, trPickle)
	SimpleSession.registerAfterModelsCB(restoreScenes, scPickle)
	SimpleSession.registerAfterModelsCB(restoreKeyframes, kfPickle)
try:
	restoreAnimation()
except:
	reportRestoreError('Error in Animate.Session')

def restoreLightController():
	import Lighting
	Lighting._setFromParams({'ratio': 1.25, 'brightness': 1.16, 'material': [30.0, (0.85, 0.85, 0.85), 1.0], 'back': [(0.35740674433659325, 0.6604015517481454, -0.6604015517481455), (1.0, 1.0, 1.0), 0.0], 'mode': 'two-point', 'key': [(-0.35740674433659325, 0.6604015517481454, 0.6604015517481455), (1.0, 1.0, 1.0), 1.0], 'contrast': 0.83, 'fill': [(0.25056280708573153, 0.25056280708573153, 0.9351131265310293), (1.0, 1.0, 1.0), 0.0]})
try:
	restoreLightController()
except:
	reportRestoreError("Error restoring lighting parameters")


def restore_surface_color_mapping():
 try:
  surface_color_state = \
   {
    'class': 'Surface_Colorings_State',
    'coloring_table': {},
    'geometry': None,
    'is_visible': False,
    'version': 2,
   }
  import SurfaceColor.session
  SurfaceColor.session.restore_surface_color_state(surface_color_state)
 except:
  reportRestoreError('Error restoring surface color mapping')

registerAfterModelsCB(restore_surface_color_mapping)


def restore_volume_data():
 volume_data_state = \
  {
   'class': 'Volume_Manager_State',
   'data_and_regions_state': [
     (
      {
       'available_subsamplings': {},
       'cell_angles': ( 90.0, 90.0, 90.0, ),
       'class': 'Data_State',
       'file_type': 'mrc',
       'grid_id': '',
       'name': 'tfiih_cluster_1_split_kin28.mrc',
       'path': 'tfiih_cluster_1_split_kin28.mrc',
       'rotation': (
         ( 1, 0, 0, ),
         ( 0, 1, 0, ),
         ( 0, 0, 1, ),
        ),
       'symmetries': ( ),
       'version': 6,
       'xyz_origin': None,
       'xyz_step': None,
      },
      [
       {
        'class': 'Volume_State',
        'default_rgba': ( 1, 1, 0.7, 1, ),
        'region': (
          ( 0, 0, 0, ),
          ( 18, 19, 19, ),
          ( 1, 1, 1, ),
         ),
        'region_list': {
          'class': 'Region_List_State',
          'current_index': 0,
          'named_regions': [ ],
          'region_list': [
            (
             ( 0, 0, 0, ),
             ( 18, 19, 19, ),
            ),
           ],
          'version': 1,
         },
        'rendering_options': {
          'box_faces': False,
          'bt_correction': 0,
          'cap_faces': 1,
          'class': 'Rendering_Options_State',
          'color_mode': 'auto8',
          'dim_transparency': 1,
          'dim_transparent_voxels': 1,
          'flip_normals': 0,
          'limit_voxel_count': 1,
          'line_thickness': 1.0,
          'linear_interpolation': 1,
          'maximum_intensity_projection': 0,
          'mesh_lighting': 1,
          'minimal_texture_memory': 0,
          'orthoplane_positions': ( 0, 0, 0, ),
          'orthoplanes_shown': ( False, False, False, ),
          'outline_box_linewidth': 1.0,
          'outline_box_rgb': ( 1.0, 1.0, 1.0, ),
          'projection_mode': 'auto',
          'show_outline_box': 0,
          'smooth_lines': 0,
          'smoothing_factor': 0.3,
          'smoothing_iterations': 2,
          'square_mesh': 1,
          'subdivide_surface': 0,
          'subdivision_levels': 1,
          'surface_smoothing': 0,
          'two_sided_lighting': 1,
          'version': 1,
          'voxel_limit': 1.0,
         },
        'representation': 'surface',
        'session_volume_id': '_>Sc\x0bj!YbufW9[\\bjH4I!xn\tkt@T+;N4',
        'solid_brightness_factor': 1.0,
        'solid_colors': [
          ( 1, 1, 0.7, 1, ),
          ( 1, 1, 0.7, 1, ),
          ( 1, 1, 0.7, 1, ),
         ],
        'solid_levels': [
          ( 0.0, 0, ),
          ( 0.37164158368706701, 0.99, ),
          ( 0.501743733882904, 1, ),
         ],
        'solid_model': None,
        'surface_brightness_factor': 1.0,
        'surface_colors': [
          ( 1.0, 0.0, 0.0, 1.0, ),
         ],
        'surface_levels': [ 0.24026230174765337, ],
        'surface_model': {
          'active': True,
          'class': 'Model_State',
          'clip_plane_normal': ( 0.0, 0.0, 0.0, ),
          'clip_plane_origin': ( 0.0, 0.0, 0.0, ),
          'clip_thickness': 5.0,
          'display': True,
          'id': 1,
          'name': u'tfiih_cluster_1_split_kin28.mrc',
          'osl_identifier': u'#1',
          'silhouette': True,
          'subid': 0,
          'use_clip_plane': False,
          'use_clip_thickness': False,
          'version': 5,
          'xform': {
            'class': 'Xform_State',
            'rotation_angle': 75.24069897560051,
            'rotation_axis': ( -0.14764781516820408, 0.7261596647515306, -0.6714851181999497, ),
            'translation': ( -95.82427864093718, 11.785462468339727, -23.801339356350766, ),
            'version': 1,
           },
         },
        'transparency_depth': 0.5,
        'transparency_factor': 0.0,
        'version': 6,
       },
      ],
     ),
     (
      {
       'available_subsamplings': {},
       'cell_angles': ( 90.0, 90.0, 90.0, ),
       'class': 'Data_State',
       'file_type': 'mrc',
       'grid_id': '',
       'name': 'tfiih_cluster_1_split_ccl1.mrc',
       'path': 'tfiih_cluster_1_split_ccl1.mrc',
       'rotation': (
         ( 1, 0, 0, ),
         ( 0, 1, 0, ),
         ( 0, 0, 1, ),
        ),
       'symmetries': ( ),
       'version': 6,
       'xyz_origin': None,
       'xyz_step': None,
      },
      [
       {
        'class': 'Volume_State',
        'default_rgba': ( 0.7, 0.7, 0.7, 1, ),
        'region': (
          ( 0, 0, 0, ),
          ( 23, 24, 22, ),
          ( 1, 1, 1, ),
         ),
        'region_list': {
          'class': 'Region_List_State',
          'current_index': 0,
          'named_regions': [ ],
          'region_list': [
            (
             ( 0, 0, 0, ),
             ( 23, 24, 22, ),
            ),
           ],
          'version': 1,
         },
        'rendering_options': {
          'box_faces': False,
          'bt_correction': 0,
          'cap_faces': 1,
          'class': 'Rendering_Options_State',
          'color_mode': 'auto8',
          'dim_transparency': 1,
          'dim_transparent_voxels': 1,
          'flip_normals': 0,
          'limit_voxel_count': 1,
          'line_thickness': 1.0,
          'linear_interpolation': 1,
          'maximum_intensity_projection': 0,
          'mesh_lighting': 1,
          'minimal_texture_memory': 0,
          'orthoplane_positions': ( 0, 0, 0, ),
          'orthoplanes_shown': ( False, False, False, ),
          'outline_box_linewidth': 1.0,
          'outline_box_rgb': ( 1.0, 1.0, 1.0, ),
          'projection_mode': 'auto',
          'show_outline_box': 0,
          'smooth_lines': 0,
          'smoothing_factor': 0.3,
          'smoothing_iterations': 2,
          'square_mesh': 1,
          'subdivide_surface': 0,
          'subdivision_levels': 1,
          'surface_smoothing': 0,
          'two_sided_lighting': 1,
          'version': 1,
          'voxel_limit': 1.0,
         },
        'representation': 'surface',
        'session_volume_id': 'B/`65w*W?:`#8~OZ*IC8Q~q\t&\nsTy%FW',
        'solid_brightness_factor': 1.0,
        'solid_colors': [
          ( 1.0, 1.0, 1.0, 1, ),
          ( 1.0, 1.0, 1.0, 1, ),
          ( 1.0, 1.0, 1.0, 1, ),
         ],
        'solid_levels': [
          ( 0.0, 0, ),
          ( 0.39316123369932177, 0.99, ),
          ( 0.6780980229377747, 1, ),
         ],
        'solid_model': None,
        'surface_brightness_factor': 1.0,
        'surface_colors': [
          ( 1.0, 0.2705882489681244, 0.0, 1.0, ),
         ],
        'surface_levels': [ 0.3269953080216683, ],
        'surface_model': {
          'active': True,
          'class': 'Model_State',
          'clip_plane_normal': ( 0.0, 0.0, 0.0, ),
          'clip_plane_origin': ( 0.0, 0.0, 0.0, ),
          'clip_thickness': 5.0,
          'display': True,
          'id': 0,
          'name': u'tfiih_cluster_1_split_ccl1.mrc',
          'osl_identifier': u'#0',
          'silhouette': True,
          'subid': 0,
          'use_clip_plane': False,
          'use_clip_thickness': False,
          'version': 5,
          'xform': {
            'class': 'Xform_State',
            'rotation_angle': 75.24069897560051,
            'rotation_axis': ( -0.14764781516820408, 0.7261596647515306, -0.6714851181999497, ),
            'translation': ( -95.82427864093718, 11.785462468339727, -23.801339356350766, ),
            'version': 1,
           },
         },
        'transparency_depth': 0.5,
        'transparency_factor': 0.0,
        'version': 6,
       },
      ],
     ),
     (
      {
       'available_subsamplings': {},
       'cell_angles': ( 90.0, 90.0, 90.0, ),
       'class': 'Data_State',
       'file_type': 'mrc',
       'grid_id': '',
       'name': 'tfiih_cluster_1_split_tfb2.mrc',
       'path': 'tfiih_cluster_1_split_tfb2.mrc',
       'rotation': (
         ( 1, 0, 0, ),
         ( 0, 1, 0, ),
         ( 0, 0, 1, ),
        ),
       'symmetries': ( ),
       'version': 6,
       'xyz_origin': None,
       'xyz_step': None,
      },
      [
       {
        'class': 'Volume_State',
        'default_rgba': ( 0.7, 1, 0.7, 1, ),
        'region': (
          ( 0, 0, 0, ),
          ( 27, 23, 20, ),
          ( 1, 1, 1, ),
         ),
        'region_list': {
          'class': 'Region_List_State',
          'current_index': 0,
          'named_regions': [ ],
          'region_list': [
            (
             ( 0, 0, 0, ),
             ( 27, 23, 20, ),
            ),
           ],
          'version': 1,
         },
        'rendering_options': {
          'box_faces': False,
          'bt_correction': 0,
          'cap_faces': 1,
          'class': 'Rendering_Options_State',
          'color_mode': 'auto8',
          'dim_transparency': 1,
          'dim_transparent_voxels': 1,
          'flip_normals': 0,
          'limit_voxel_count': 1,
          'line_thickness': 1.0,
          'linear_interpolation': 1,
          'maximum_intensity_projection': 0,
          'mesh_lighting': 1,
          'minimal_texture_memory': 0,
          'orthoplane_positions': ( 0, 0, 0, ),
          'orthoplanes_shown': ( False, False, False, ),
          'outline_box_linewidth': 1.0,
          'outline_box_rgb': ( 1.0, 1.0, 1.0, ),
          'projection_mode': 'auto',
          'show_outline_box': 0,
          'smooth_lines': 0,
          'smoothing_factor': 0.3,
          'smoothing_iterations': 2,
          'square_mesh': 1,
          'subdivide_surface': 0,
          'subdivision_levels': 1,
          'surface_smoothing': 0,
          'two_sided_lighting': 1,
          'version': 1,
          'voxel_limit': 1.0,
         },
        'representation': 'surface',
        'session_volume_id': 'Ia^P&_TPhP87b<w8@!.Tsa.}d=My,#W]',
        'solid_brightness_factor': 1.0,
        'solid_colors': [
          ( 0.7, 1, 0.7, 1, ),
          ( 0.7, 1, 0.7, 1, ),
          ( 0.7, 1, 0.7, 1, ),
         ],
        'solid_levels': [
          ( 0.0, 0, ),
          ( 0.4180397817492485, 0.99, ),
          ( 0.6406739950180054, 1, ),
         ],
        'solid_model': None,
        'surface_brightness_factor': 1.0,
        'surface_colors': [
          ( 0.0, 1.0, 1.0, 1.0, ),
         ],
        'surface_levels': [ 0.3279193988183268, ],
        'surface_model': {
          'active': True,
          'class': 'Model_State',
          'clip_plane_normal': ( 0.0, 0.0, 0.0, ),
          'clip_plane_origin': ( 0.0, 0.0, 0.0, ),
          'clip_thickness': 5.0,
          'display': True,
          'id': 6,
          'name': u'tfiih_cluster_1_split_tfb2.mrc',
          'osl_identifier': u'#6',
          'silhouette': True,
          'subid': 0,
          'use_clip_plane': False,
          'use_clip_thickness': False,
          'version': 5,
          'xform': {
            'class': 'Xform_State',
            'rotation_angle': 75.24069897560051,
            'rotation_axis': ( -0.14764781516820408, 0.7261596647515306, -0.6714851181999497, ),
            'translation': ( -95.82427864093718, 11.785462468339727, -23.801339356350766, ),
            'version': 1,
           },
         },
        'transparency_depth': 0.5,
        'transparency_factor': 0.0,
        'version': 6,
       },
      ],
     ),
     (
      {
       'available_subsamplings': {},
       'cell_angles': ( 90.0, 90.0, 90.0, ),
       'class': 'Data_State',
       'file_type': 'mrc',
       'grid_id': '',
       'name': 'tfiih_cluster_1_split_tfb5.mrc',
       'path': 'tfiih_cluster_1_split_tfb5.mrc',
       'rotation': (
         ( 1, 0, 0, ),
         ( 0, 1, 0, ),
         ( 0, 0, 1, ),
        ),
       'symmetries': ( ),
       'version': 6,
       'xyz_origin': None,
       'xyz_step': None,
      },
      [
       {
        'class': 'Volume_State',
        'default_rgba': ( 0.8, 0.8, 0.6, 1, ),
        'region': (
          ( 0, 0, 0, ),
          ( 14, 15, 13, ),
          ( 1, 1, 1, ),
         ),
        'region_list': {
          'class': 'Region_List_State',
          'current_index': 0,
          'named_regions': [ ],
          'region_list': [
            (
             ( 0, 0, 0, ),
             ( 14, 15, 13, ),
            ),
           ],
          'version': 1,
         },
        'rendering_options': {
          'box_faces': False,
          'bt_correction': 0,
          'cap_faces': 1,
          'class': 'Rendering_Options_State',
          'color_mode': 'auto8',
          'dim_transparency': 1,
          'dim_transparent_voxels': 1,
          'flip_normals': 0,
          'limit_voxel_count': 1,
          'line_thickness': 1.0,
          'linear_interpolation': 1,
          'maximum_intensity_projection': 0,
          'mesh_lighting': 1,
          'minimal_texture_memory': 0,
          'orthoplane_positions': ( 0, 0, 0, ),
          'orthoplanes_shown': ( False, False, False, ),
          'outline_box_linewidth': 1.0,
          'outline_box_rgb': ( 1.0, 1.0, 1.0, ),
          'projection_mode': 'auto',
          'show_outline_box': 0,
          'smooth_lines': 0,
          'smoothing_factor': 0.3,
          'smoothing_iterations': 2,
          'square_mesh': 1,
          'subdivide_surface': 0,
          'subdivision_levels': 1,
          'surface_smoothing': 0,
          'two_sided_lighting': 1,
          'version': 1,
          'voxel_limit': 1.0,
         },
        'representation': 'surface',
        'session_volume_id': "`O70L=o_V+eS%$Ml'1em^L>;_5/{3dh)",
        'solid_brightness_factor': 1.0,
        'solid_colors': [
          ( 1.0, 1.0, 0.7499999999999999, 1, ),
          ( 1.0, 1.0, 0.7499999999999999, 1, ),
          ( 1.0, 1.0, 0.7499999999999999, 1, ),
         ],
        'solid_levels': [
          ( 0.0, 0, ),
          ( 0.1857056731581688, 0.99, ),
          ( 0.293281227350235, 1, ),
         ],
        'solid_model': None,
        'surface_brightness_factor': 1.0,
        'surface_colors': [
          ( 0.4117647111415863, 0.4117647111415863, 0.4117647111415863, 1.0, ),
         ],
        'surface_levels': [ 0.1144480859503109, ],
        'surface_model': {
          'active': True,
          'class': 'Model_State',
          'clip_plane_normal': ( 0.0, 0.0, 0.0, ),
          'clip_plane_origin': ( 0.0, 0.0, 0.0, ),
          'clip_thickness': 5.0,
          'display': True,
          'id': 9,
          'name': u'tfiih_cluster_1_split_tfb5.mrc',
          'osl_identifier': u'#9',
          'silhouette': True,
          'subid': 0,
          'use_clip_plane': False,
          'use_clip_thickness': False,
          'version': 5,
          'xform': {
            'class': 'Xform_State',
            'rotation_angle': 75.24069897560051,
            'rotation_axis': ( -0.14764781516820408, 0.7261596647515306, -0.6714851181999497, ),
            'translation': ( -95.82427864093718, 11.785462468339727, -23.801339356350766, ),
            'version': 1,
           },
         },
        'transparency_depth': 0.5,
        'transparency_factor': 0.0,
        'version': 6,
       },
      ],
     ),
     (
      {
       'available_subsamplings': {},
       'cell_angles': ( 90.0, 90.0, 90.0, ),
       'class': 'Data_State',
       'file_type': 'mrc',
       'grid_id': '',
       'name': 'tfiih_cluster_1_split_tfb4.mrc',
       'path': 'tfiih_cluster_1_split_tfb4.mrc',
       'rotation': (
         ( 1, 0, 0, ),
         ( 0, 1, 0, ),
         ( 0, 0, 1, ),
        ),
       'symmetries': ( ),
       'version': 6,
       'xyz_origin': None,
       'xyz_step': None,
      },
      [
       {
        'class': 'Volume_State',
        'default_rgba': ( 0.6, 0.75, 0.9, 1, ),
        'region': (
          ( 0, 0, 0, ),
          ( 20, 22, 20, ),
          ( 1, 1, 1, ),
         ),
        'region_list': {
          'class': 'Region_List_State',
          'current_index': 0,
          'named_regions': [ ],
          'region_list': [
            (
             ( 0, 0, 0, ),
             ( 20, 22, 20, ),
            ),
           ],
          'version': 1,
         },
        'rendering_options': {
          'box_faces': False,
          'bt_correction': 0,
          'cap_faces': 1,
          'class': 'Rendering_Options_State',
          'color_mode': 'auto8',
          'dim_transparency': 1,
          'dim_transparent_voxels': 1,
          'flip_normals': 0,
          'limit_voxel_count': 1,
          'line_thickness': 1.0,
          'linear_interpolation': 1,
          'maximum_intensity_projection': 0,
          'mesh_lighting': 1,
          'minimal_texture_memory': 0,
          'orthoplane_positions': ( 0, 0, 0, ),
          'orthoplanes_shown': ( False, False, False, ),
          'outline_box_linewidth': 1.0,
          'outline_box_rgb': ( 1.0, 1.0, 1.0, ),
          'projection_mode': 'auto',
          'show_outline_box': 0,
          'smooth_lines': 0,
          'smoothing_factor': 0.3,
          'smoothing_iterations': 2,
          'square_mesh': 1,
          'subdivide_surface': 0,
          'subdivision_levels': 1,
          'surface_smoothing': 0,
          'two_sided_lighting': 1,
          'version': 1,
          'voxel_limit': 1.0,
         },
        'representation': 'surface',
        'session_volume_id': '|{!mk ] p}`8OKS(gt*ttyiF\nOjvbny\x0b',
        'solid_brightness_factor': 1.0,
        'solid_colors': [
          ( 0.6666666666666666, 0.8333333333333333, 1.0, 1, ),
          ( 0.6666666666666666, 0.8333333333333333, 1.0, 1, ),
          ( 0.6666666666666666, 0.8333333333333333, 1.0, 1, ),
         ],
        'solid_levels': [
          ( 0.0, 0, ),
          ( 0.46813217141628266, 0.99, ),
          ( 0.7528661489486694, 1, ),
         ],
        'solid_model': None,
        'surface_brightness_factor': 1.0,
        'surface_colors': [
          ( 0.3921568691730499, 0.5843137502670288, 0.929411768913269, 1.0, ),
         ],
        'surface_levels': [ 0.40266554976857205, ],
        'surface_model': {
          'active': True,
          'class': 'Model_State',
          'clip_plane_normal': ( 0.0, 0.0, 0.0, ),
          'clip_plane_origin': ( 0.0, 0.0, 0.0, ),
          'clip_thickness': 5.0,
          'display': True,
          'id': 8,
          'name': u'tfiih_cluster_1_split_tfb4.mrc',
          'osl_identifier': u'#8',
          'silhouette': True,
          'subid': 0,
          'use_clip_plane': False,
          'use_clip_thickness': False,
          'version': 5,
          'xform': {
            'class': 'Xform_State',
            'rotation_angle': 75.24069897560051,
            'rotation_axis': ( -0.14764781516820408, 0.7261596647515306, -0.6714851181999497, ),
            'translation': ( -95.82427864093718, 11.785462468339727, -23.801339356350766, ),
            'version': 1,
           },
         },
        'transparency_depth': 0.5,
        'transparency_factor': 0.0,
        'version': 6,
       },
      ],
     ),
     (
      {
       'available_subsamplings': {},
       'cell_angles': ( 90.0, 90.0, 90.0, ),
       'class': 'Data_State',
       'file_type': 'mrc',
       'grid_id': '',
       'name': 'tfiih_cluster_1_split_ssl2.mrc',
       'path': 'tfiih_cluster_1_split_ssl2.mrc',
       'rotation': (
         ( 1, 0, 0, ),
         ( 0, 1, 0, ),
         ( 0, 0, 1, ),
        ),
       'symmetries': ( ),
       'version': 6,
       'xyz_origin': None,
       'xyz_step': None,
      },
      [
       {
        'class': 'Volume_State',
        'default_rgba': ( 1, 0.7, 1, 1, ),
        'region': (
          ( 0, 0, 0, ),
          ( 23, 25, 24, ),
          ( 1, 1, 1, ),
         ),
        'region_list': {
          'class': 'Region_List_State',
          'current_index': 0,
          'named_regions': [ ],
          'region_list': [
            (
             ( 0, 0, 0, ),
             ( 23, 25, 24, ),
            ),
           ],
          'version': 1,
         },
        'rendering_options': {
          'box_faces': False,
          'bt_correction': 0,
          'cap_faces': 1,
          'class': 'Rendering_Options_State',
          'color_mode': 'auto8',
          'dim_transparency': 1,
          'dim_transparent_voxels': 1,
          'flip_normals': 0,
          'limit_voxel_count': 1,
          'line_thickness': 1.0,
          'linear_interpolation': 1,
          'maximum_intensity_projection': 0,
          'mesh_lighting': 1,
          'minimal_texture_memory': 0,
          'orthoplane_positions': ( 0, 0, 0, ),
          'orthoplanes_shown': ( False, False, False, ),
          'outline_box_linewidth': 1.0,
          'outline_box_rgb': ( 1.0, 1.0, 1.0, ),
          'projection_mode': 'auto',
          'show_outline_box': 0,
          'smooth_lines': 0,
          'smoothing_factor': 0.3,
          'smoothing_iterations': 2,
          'square_mesh': 1,
          'subdivide_surface': 0,
          'subdivision_levels': 1,
          'surface_smoothing': 0,
          'two_sided_lighting': 1,
          'version': 1,
          'voxel_limit': 1.0,
         },
        'representation': 'surface',
        'session_volume_id': '+y([uuzB]Coe6w.Ik]C[kYTdIwaR&R=S',
        'solid_brightness_factor': 1.0,
        'solid_colors': [
          ( 1, 0.7, 1, 1, ),
          ( 1, 0.7, 1, 1, ),
          ( 1, 0.7, 1, 1, ),
         ],
        'solid_levels': [
          ( 0.0, 0, ),
          ( 0.69771231333017347, 0.99, ),
          ( 1.0544239282608032, 1, ),
         ],
        'solid_model': None,
        'surface_brightness_factor': 1.0,
        'surface_colors': [
          ( 0.0, 0.0, 1.0, 1.0, ),
         ],
        'surface_levels': [ 0.4969970985597165, ],
        'surface_model': {
          'active': True,
          'class': 'Model_State',
          'clip_plane_normal': ( 0.0, 0.0, 0.0, ),
          'clip_plane_origin': ( 0.0, 0.0, 0.0, ),
          'clip_thickness': 5.0,
          'display': True,
          'id': 4,
          'name': u'tfiih_cluster_1_split_ssl2.mrc',
          'osl_identifier': u'#4',
          'silhouette': True,
          'subid': 0,
          'use_clip_plane': False,
          'use_clip_thickness': False,
          'version': 5,
          'xform': {
            'class': 'Xform_State',
            'rotation_angle': 75.24069897560051,
            'rotation_axis': ( -0.14764781516820408, 0.7261596647515306, -0.6714851181999497, ),
            'translation': ( -95.82427864093718, 11.785462468339727, -23.801339356350766, ),
            'version': 1,
           },
         },
        'transparency_depth': 0.5,
        'transparency_factor': 0.0,
        'version': 6,
       },
      ],
     ),
     (
      {
       'available_subsamplings': {},
       'cell_angles': ( 90.0, 90.0, 90.0, ),
       'class': 'Data_State',
       'file_type': 'mrc',
       'grid_id': '',
       'name': 'tfiih_cluster_1_split_ssl1.mrc',
       'path': 'tfiih_cluster_1_split_ssl1.mrc',
       'rotation': (
         ( 1, 0, 0, ),
         ( 0, 1, 0, ),
         ( 0, 0, 1, ),
        ),
       'symmetries': ( ),
       'version': 6,
       'xyz_origin': None,
       'xyz_step': None,
      },
      [
       {
        'class': 'Volume_State',
        'default_rgba': ( 0.7, 0.7, 1, 1, ),
        'region': (
          ( 0, 0, 0, ),
          ( 24, 22, 19, ),
          ( 1, 1, 1, ),
         ),
        'region_list': {
          'class': 'Region_List_State',
          'current_index': 0,
          'named_regions': [ ],
          'region_list': [
            (
             ( 0, 0, 0, ),
             ( 24, 22, 19, ),
            ),
           ],
          'version': 1,
         },
        'rendering_options': {
          'box_faces': False,
          'bt_correction': 0,
          'cap_faces': 1,
          'class': 'Rendering_Options_State',
          'color_mode': 'auto8',
          'dim_transparency': 1,
          'dim_transparent_voxels': 1,
          'flip_normals': 0,
          'limit_voxel_count': 1,
          'line_thickness': 1.0,
          'linear_interpolation': 1,
          'maximum_intensity_projection': 0,
          'mesh_lighting': 1,
          'minimal_texture_memory': 0,
          'orthoplane_positions': ( 0, 0, 0, ),
          'orthoplanes_shown': ( False, False, False, ),
          'outline_box_linewidth': 1.0,
          'outline_box_rgb': ( 1.0, 1.0, 1.0, ),
          'projection_mode': 'auto',
          'show_outline_box': 0,
          'smooth_lines': 0,
          'smoothing_factor': 0.3,
          'smoothing_iterations': 2,
          'square_mesh': 1,
          'subdivide_surface': 0,
          'subdivision_levels': 1,
          'surface_smoothing': 0,
          'two_sided_lighting': 1,
          'version': 1,
          'voxel_limit': 1.0,
         },
        'representation': 'surface',
        'session_volume_id': 'C=sqSr)rl,>/0\t\t`<BN9|\nr\rI8Bl|, M',
        'solid_brightness_factor': 1.0,
        'solid_colors': [
          ( 0.7, 0.7, 1, 1, ),
          ( 0.7, 0.7, 1, 1, ),
          ( 0.7, 0.7, 1, 1, ),
         ],
        'solid_levels': [
          ( 0.0, 0, ),
          ( 0.53223569908142088, 0.99, ),
          ( 0.7841987609863281, 1, ),
         ],
        'solid_model': None,
        'surface_brightness_factor': 1.0,
        'surface_colors': [
          ( 0.6274510025978088, 0.125490203499794, 0.9411764740943909, 1.0, ),
         ],
        'surface_levels': [ 0.4246772655296804, ],
        'surface_model': {
          'active': True,
          'class': 'Model_State',
          'clip_plane_normal': ( 0.0, 0.0, 0.0, ),
          'clip_plane_origin': ( 0.0, 0.0, 0.0, ),
          'clip_thickness': 5.0,
          'display': True,
          'id': 3,
          'name': u'tfiih_cluster_1_split_ssl1.mrc',
          'osl_identifier': u'#3',
          'silhouette': True,
          'subid': 0,
          'use_clip_plane': False,
          'use_clip_thickness': False,
          'version': 5,
          'xform': {
            'class': 'Xform_State',
            'rotation_angle': 75.24069897560051,
            'rotation_axis': ( -0.14764781516820408, 0.7261596647515306, -0.6714851181999497, ),
            'translation': ( -95.82427864093718, 11.785462468339727, -23.801339356350766, ),
            'version': 1,
           },
         },
        'transparency_depth': 0.5,
        'transparency_factor': 0.0,
        'version': 6,
       },
      ],
     ),
     (
      {
       'available_subsamplings': {},
       'cell_angles': ( 90.0, 90.0, 90.0, ),
       'class': 'Data_State',
       'file_type': 'mrc',
       'grid_id': '',
       'name': 'tfiih_cluster_1_split_rad3.mrc',
       'path': 'tfiih_cluster_1_split_rad3.mrc',
       'rotation': (
         ( 1, 0, 0, ),
         ( 0, 1, 0, ),
         ( 0, 0, 1, ),
        ),
       'symmetries': ( ),
       'version': 6,
       'xyz_origin': None,
       'xyz_step': None,
      },
      [
       {
        'class': 'Volume_State',
        'default_rgba': ( 0.7, 1, 1, 1, ),
        'region': (
          ( 0, 0, 0, ),
          ( 26, 26, 25, ),
          ( 1, 1, 1, ),
         ),
        'region_list': {
          'class': 'Region_List_State',
          'current_index': 0,
          'named_regions': [ ],
          'region_list': [
            (
             ( 0, 0, 0, ),
             ( 26, 26, 25, ),
            ),
           ],
          'version': 1,
         },
        'rendering_options': {
          'box_faces': False,
          'bt_correction': 0,
          'cap_faces': 1,
          'class': 'Rendering_Options_State',
          'color_mode': 'auto8',
          'dim_transparency': 1,
          'dim_transparent_voxels': 1,
          'flip_normals': 0,
          'limit_voxel_count': 1,
          'line_thickness': 1.0,
          'linear_interpolation': 1,
          'maximum_intensity_projection': 0,
          'mesh_lighting': 1,
          'minimal_texture_memory': 0,
          'orthoplane_positions': ( 0, 0, 0, ),
          'orthoplanes_shown': ( False, False, False, ),
          'outline_box_linewidth': 1.0,
          'outline_box_rgb': ( 1.0, 1.0, 1.0, ),
          'projection_mode': 'auto',
          'show_outline_box': 0,
          'smooth_lines': 0,
          'smoothing_factor': 0.3,
          'smoothing_iterations': 2,
          'square_mesh': 1,
          'subdivide_surface': 0,
          'subdivision_levels': 1,
          'surface_smoothing': 0,
          'two_sided_lighting': 1,
          'version': 1,
          'voxel_limit': 1.0,
         },
        'representation': 'surface',
        'session_volume_id': "D*b-6GPTd>*0l5tC[Y@B)\\cZJz='4nel",
        'solid_brightness_factor': 1.0,
        'solid_colors': [
          ( 0.7, 1, 1, 1, ),
          ( 0.7, 1, 1, 1, ),
          ( 0.7, 1, 1, 1, ),
         ],
        'solid_levels': [
          ( 0.0, 0, ),
          ( 0.58423492163419721, 0.99, ),
          ( 0.9114429354667664, 1, ),
         ],
        'solid_model': None,
        'surface_brightness_factor': 1.0,
        'surface_colors': [
          ( 0.13333334028720856, 0.545098066329956, 0.13333334028720856, 1.0, ),
         ],
        'surface_levels': [ 0.42058431975487764, ],
        'surface_model': {
          'active': True,
          'class': 'Model_State',
          'clip_plane_normal': ( 0.0, 0.0, 0.0, ),
          'clip_plane_origin': ( 0.0, 0.0, 0.0, ),
          'clip_thickness': 5.0,
          'display': True,
          'id': 2,
          'name': u'tfiih_cluster_1_split_rad3.mrc',
          'osl_identifier': u'#2',
          'silhouette': True,
          'subid': 0,
          'use_clip_plane': False,
          'use_clip_thickness': False,
          'version': 5,
          'xform': {
            'class': 'Xform_State',
            'rotation_angle': 75.24069897560051,
            'rotation_axis': ( -0.14764781516820408, 0.7261596647515306, -0.6714851181999497, ),
            'translation': ( -95.82427864093718, 11.785462468339727, -23.801339356350766, ),
            'version': 1,
           },
         },
        'transparency_depth': 0.5,
        'transparency_factor': 0.0,
        'version': 6,
       },
      ],
     ),
     (
      {
       'available_subsamplings': {},
       'cell_angles': ( 90.0, 90.0, 90.0, ),
       'class': 'Data_State',
       'file_type': 'mrc',
       'grid_id': '',
       'name': 'tfiih_cluster_1_split_tfb1.mrc',
       'path': 'tfiih_cluster_1_split_tfb1.mrc',
       'rotation': (
         ( 1, 0, 0, ),
         ( 0, 1, 0, ),
         ( 0, 0, 1, ),
        ),
       'symmetries': ( ),
       'version': 6,
       'xyz_origin': None,
       'xyz_step': None,
      },
      [
       {
        'class': 'Volume_State',
        'default_rgba': ( 1, 0.7, 0.7, 1, ),
        'region': (
          ( 0, 0, 0, ),
          ( 20, 22, 22, ),
          ( 1, 1, 1, ),
         ),
        'region_list': {
          'class': 'Region_List_State',
          'current_index': 0,
          'named_regions': [ ],
          'region_list': [
            (
             ( 0, 0, 0, ),
             ( 20, 22, 22, ),
            ),
           ],
          'version': 1,
         },
        'rendering_options': {
          'box_faces': False,
          'bt_correction': 0,
          'cap_faces': 1,
          'class': 'Rendering_Options_State',
          'color_mode': 'auto8',
          'dim_transparency': 1,
          'dim_transparent_voxels': 1,
          'flip_normals': 0,
          'limit_voxel_count': 1,
          'line_thickness': 1.0,
          'linear_interpolation': 1,
          'maximum_intensity_projection': 0,
          'mesh_lighting': 1,
          'minimal_texture_memory': 0,
          'orthoplane_positions': ( 0, 0, 0, ),
          'orthoplanes_shown': ( False, False, False, ),
          'outline_box_linewidth': 1.0,
          'outline_box_rgb': ( 1.0, 1.0, 1.0, ),
          'projection_mode': 'auto',
          'show_outline_box': 0,
          'smooth_lines': 0,
          'smoothing_factor': 0.3,
          'smoothing_iterations': 2,
          'square_mesh': 1,
          'subdivide_surface': 0,
          'subdivision_levels': 1,
          'surface_smoothing': 0,
          'two_sided_lighting': 1,
          'version': 1,
          'voxel_limit': 1.0,
         },
        'representation': 'surface',
        'session_volume_id': 'k<\x0bATz/}F`u_?e\x0bC$3kJs1yeG2BWo-/R',
        'solid_brightness_factor': 1.0,
        'solid_colors': [
          ( 1, 0.7, 0.7, 1, ),
          ( 1, 0.7, 0.7, 1, ),
          ( 1, 0.7, 0.7, 1, ),
         ],
        'solid_levels': [
          ( 0.0, 0, ),
          ( 0.6018524921894074, 0.99, ),
          ( 0.8639857769012451, 1, ),
         ],
        'solid_model': None,
        'surface_brightness_factor': 1.0,
        'surface_colors': [
          ( 1.0, 1.0, 0.0, 1.0, ),
         ],
        'surface_levels': [ 0.38088426280289944, ],
        'surface_model': {
          'active': True,
          'class': 'Model_State',
          'clip_plane_normal': ( 0.0, 0.0, 0.0, ),
          'clip_plane_origin': ( 0.0, 0.0, 0.0, ),
          'clip_thickness': 5.0,
          'display': True,
          'id': 5,
          'name': u'tfiih_cluster_1_split_tfb1.mrc',
          'osl_identifier': u'#5',
          'silhouette': True,
          'subid': 0,
          'use_clip_plane': False,
          'use_clip_thickness': False,
          'version': 5,
          'xform': {
            'class': 'Xform_State',
            'rotation_angle': 75.24069897560051,
            'rotation_axis': ( -0.14764781516820408, 0.7261596647515306, -0.6714851181999497, ),
            'translation': ( -95.82427864093718, 11.785462468339727, -23.801339356350766, ),
            'version': 1,
           },
         },
        'transparency_depth': 0.5,
        'transparency_factor': 0.0,
        'version': 6,
       },
      ],
     ),
     (
      {
       'available_subsamplings': {},
       'cell_angles': ( 90.0, 90.0, 90.0, ),
       'class': 'Data_State',
       'file_type': 'mrc',
       'grid_id': '',
       'name': 'tfiih_cluster_1_split_tfb3.mrc',
       'path': 'tfiih_cluster_1_split_tfb3.mrc',
       'rotation': (
         ( 1, 0, 0, ),
         ( 0, 1, 0, ),
         ( 0, 0, 1, ),
        ),
       'symmetries': ( ),
       'version': 6,
       'xyz_origin': None,
       'xyz_step': None,
      },
      [
       {
        'class': 'Volume_State',
        'default_rgba': ( 0.9, 0.75, 0.6, 1, ),
        'region': (
          ( 0, 0, 0, ),
          ( 21, 23, 22, ),
          ( 1, 1, 1, ),
         ),
        'region_list': {
          'class': 'Region_List_State',
          'current_index': 0,
          'named_regions': [ ],
          'region_list': [
            (
             ( 0, 0, 0, ),
             ( 21, 23, 22, ),
            ),
           ],
          'version': 1,
         },
        'rendering_options': {
          'box_faces': False,
          'bt_correction': 0,
          'cap_faces': 1,
          'class': 'Rendering_Options_State',
          'color_mode': 'auto8',
          'dim_transparency': 1,
          'dim_transparent_voxels': 1,
          'flip_normals': 0,
          'limit_voxel_count': 1,
          'line_thickness': 1.0,
          'linear_interpolation': 1,
          'maximum_intensity_projection': 0,
          'mesh_lighting': 1,
          'minimal_texture_memory': 0,
          'orthoplane_positions': ( 0, 0, 0, ),
          'orthoplanes_shown': ( False, False, False, ),
          'outline_box_linewidth': 1.0,
          'outline_box_rgb': ( 1.0, 1.0, 1.0, ),
          'projection_mode': 'auto',
          'show_outline_box': 0,
          'smooth_lines': 0,
          'smoothing_factor': 0.3,
          'smoothing_iterations': 2,
          'square_mesh': 1,
          'subdivide_surface': 0,
          'subdivision_levels': 1,
          'surface_smoothing': 0,
          'two_sided_lighting': 1,
          'version': 1,
          'voxel_limit': 1.0,
         },
        'representation': 'surface',
        'session_volume_id': '&RHA,g[g+e,0b@DyQ?]Sb-*A\rUI\r\x0c\nd0',
        'solid_brightness_factor': 1.0,
        'solid_colors': [
          ( 1.0, 0.8333333333333333, 0.6666666666666666, 1, ),
          ( 1.0, 0.8333333333333333, 0.6666666666666666, 1, ),
          ( 1.0, 0.8333333333333333, 0.6666666666666666, 1, ),
         ],
        'solid_levels': [
          ( 0.0, 0, ),
          ( 0.34685057244300843, 0.99, ),
          ( 0.5999836921691895, 1, ),
         ],
        'solid_model': None,
        'surface_brightness_factor': 1.0,
        'surface_colors': [
          ( 1.0, 0.49803921580314636, 0.0, 1.0, ),
         ],
        'surface_levels': [ 0.3115574140801152, ],
        'surface_model': {
          'active': True,
          'class': 'Model_State',
          'clip_plane_normal': ( 0.0, 0.0, 0.0, ),
          'clip_plane_origin': ( 0.0, 0.0, 0.0, ),
          'clip_thickness': 5.0,
          'display': True,
          'id': 7,
          'name': u'tfiih_cluster_1_split_tfb3.mrc',
          'osl_identifier': u'#7',
          'silhouette': True,
          'subid': 0,
          'use_clip_plane': False,
          'use_clip_thickness': False,
          'version': 5,
          'xform': {
            'class': 'Xform_State',
            'rotation_angle': 75.24069897560051,
            'rotation_axis': ( -0.14764781516820408, 0.7261596647515306, -0.6714851181999497, ),
            'translation': ( -95.82427864093718, 11.785462468339727, -23.801339356350766, ),
            'version': 1,
           },
         },
        'transparency_depth': 0.5,
        'transparency_factor': 0.0,
        'version': 6,
       },
      ],
     ),
    ],
   'version': 2,
  }
 from VolumeViewer import session
 session.restore_volume_data_state(volume_data_state)

try:
  restore_volume_data()
except:
  reportRestoreError('Error restoring volume data')


def restore_volume_dialog():
 volume_dialog_state = \
  {
   'adjust_camera': 0,
   'auto_show_subregion': 0,
   'box_padding': '0',
   'class': 'Volume_Dialog_State',
   'data_cache_size': '512',
   'focus_volume': '&RHA,g[g+e,0b@DyQ?]Sb-*A\rUI\r\x0c\nd0',
   'geometry': u'417x405+1566+790',
   'histogram_active_order': [ 2, 0, 1, ],
   'histogram_volumes': [ 'B/`65w*W?:`#8~OZ*IC8Q~q\t&\nsTy%FW', '_>Sc\x0bj!YbufW9[\\bjH4I!xn\tkt@T+;N4', '&RHA,g[g+e,0b@DyQ?]Sb-*A\rUI\r\x0c\nd0', ],
   'immediate_update': 1,
   'initial_colors': (
     ( 0.7, 0.7, 0.7, 1, ),
     ( 1, 1, 0.7, 1, ),
     ( 0.7, 1, 1, 1, ),
     ( 0.7, 0.7, 1, 1, ),
     ( 1, 0.7, 1, 1, ),
     ( 1, 0.7, 0.7, 1, ),
     ( 0.7, 1, 0.7, 1, ),
     ( 0.9, 0.75, 0.6, 1, ),
     ( 0.6, 0.75, 0.9, 1, ),
     ( 0.8, 0.8, 0.6, 1, ),
    ),
   'is_visible': True,
   'max_histograms': '3',
   'representation': 'surface',
   'selectable_subregions': 0,
   'show_on_open': 1,
   'show_plane': 1,
   'shown_panels': [ 'Threshold and Color', 'Display style', ],
   'subregion_button': 'middle',
   'use_initial_colors': 1,
   'version': 12,
   'voxel_limit_for_open': '256',
   'voxel_limit_for_plane': '256',
   'zone_radius': 2.0,
  }
 from VolumeViewer import session
 session.restore_volume_dialog_state(volume_dialog_state)

try:
  restore_volume_dialog()
except:
  reportRestoreError('Error restoring volume dialog')


def restoreRemainder():
	from SimpleSession.versions.v62 import restoreWindowSize, \
	     restoreOpenStates, restoreSelections, restoreFontInfo, \
	     restoreOpenModelsAttrs, restoreModelClip, restoreSilhouettes

	curSelIds =  []
	savedSels = []
	openModelsAttrs = { 'cofrMethod': 4 }
	windowSize = (1203, 1189)
	xformMap = {}
	fontInfo = {'face': ('Sans Serif', 'Normal', 16)}
	clipPlaneInfo = {}
	silhouettes = {0: True}

	replyobj.status("Restoring window...", blankAfter=0,
		secondary=True)
	restoreWindowSize(windowSize)
	replyobj.status("Restoring open states...", blankAfter=0,
		secondary=True)
	restoreOpenStates(xformMap)
	replyobj.status("Restoring font info...", blankAfter=0,
		secondary=True)
	restoreFontInfo(fontInfo)
	replyobj.status("Restoring selections...", blankAfter=0,
		secondary=True)
	restoreSelections(curSelIds, savedSels)
	replyobj.status("Restoring openModel attributes...", blankAfter=0,
		secondary=True)
	restoreOpenModelsAttrs(openModelsAttrs)
	replyobj.status("Restoring model clipping...", blankAfter=0,
		secondary=True)
	restoreModelClip(clipPlaneInfo)
	replyobj.status("Restoring per-model silhouettes...", blankAfter=0,
		secondary=True)
	restoreSilhouettes(silhouettes)

	replyobj.status("Restoring remaining extension info...", blankAfter=0,
		secondary=True)
try:
	restoreRemainder()
except:
	reportRestoreError("Error restoring post-model state")
from SimpleSession.versions.v62 import makeAfterModelsCBs
makeAfterModelsCBs()

from SimpleSession.versions.v62 import endRestore
replyobj.status('Finishing restore...', blankAfter=0, secondary=True)
endRestore({})
replyobj.status('', secondary=True)
replyobj.status('Restore finished.')

