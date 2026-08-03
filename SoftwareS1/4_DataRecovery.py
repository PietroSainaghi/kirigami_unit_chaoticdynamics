###############################################################################################################
#                                                                                                             #
# This Python generates an Abaqus model for the Kirigami geometry proposed.
###############################################################################################################

from abaqus import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import numpy as np
import os
import csv
import math



# Set what rows from the CSV to analyze
# If zero is not present, the linear buckling analysis run (job submitted in line 465) will not be used by the nonlinear analysis, and the code will throw an error if 'Kirigami_Buckle_0' is not present
limInf = 600
limTop = 603


nLoadSteps = 2

# Model Construction for each set of parameters
#Read Results and Extract Desired Data
from odbAccess import *

for m in range(limInf,limTop):
    f = open('Kirigami_DispMidNodes'+str(m)+'.txt', "w")
    g = open('Kirigami_CoorMidNodes'+str(m)+'.txt', "w")
    h = open('Kirigami_DispLeftNodes'+str(m)+'.txt', "w")
    i = open('Kirigami_CoorLeftNodes'+str(m)+'.txt', "w")
    j = open('Kirigami_DispRightNodes'+str(m)+'.txt', "w")
    k = open('Kirigami_CoorRightNodes'+str(m)+'.txt', "w")
    l = open('Kirigami_ReactionForce'+str(m)+'.txt', "w")
    ll = open('Kirigami_Disp'+str(m)+'.txt', "w")
    mm = open('Kirigami_CenterLeftU'+str(m)+'.txt', "w")	
    nn = open('Kirigami_CenterRightU'+str(m)+'.txt', "w")

    JobNameNLGeom='Kirigami_NLGeom_'+str(m)+'.odb'
    myOdb = session.openOdb(name=JobNameNLGeom)
    session.viewports['Viewport: 1'].setValues(displayedObject=myOdb)

    # Get the frame repository for the step, find number of frames (starts at frame 0)
    firstFrame = myOdb.steps['Traction0'].frames[0]
    lastFrame = myOdb.steps['Traction'+str(nLoadSteps-1)].frames[-1]

    # Isolate the instance, get the number of nodes and elements
    myInstance = myOdb.rootAssembly.instances['Kirigami'.upper()]

    #Get Force and Displacement Data
    rfnodes = myOdb.rootAssembly.nodeSets['FIXEDEND']
    utracknode = myOdb.rootAssembly.nodeSets['TRACKNODE']
    for nsteps in range(0, nLoadSteps):
        frameRepository = myOdb.steps['Traction'+str(nsteps)].frames
        for frame in frameRepository:
            rf=frame.fieldOutputs['RF'].getSubset(region=rfnodes)
            for fr in rf.values:
                l.write(str(fr.nodeLabel))
                l.write(',')
                l.write(str(fr.dataDouble[0]))
                l.write('\n')
            u2=frame.fieldOutputs['U'].getSubset(region=utracknode)
            for un in u2.values:
                ##ll.write(str(un.nodeLabel))
                ll.write(str(frame.description[31:42]))
                ll.write(',')
                ll.write(str(un.dataDouble[0]))
                ll.write('\n')
    l.close()
    ll.close()

    #Get Displacement Data at Certain Nodes
    leftCentralNode = myOdb.rootAssembly.nodeSets['LEFTCENTERNODE']
    rightCentralNode = myOdb.rootAssembly.nodeSets['RIGHTCENTERNODE']
    for nsteps in range(0, nLoadSteps):
        frameRepository = myOdb.steps['Traction'+str(nsteps)].frames
        for frame in frameRepository:
            u3l=frame.fieldOutputs['U'].getSubset(region=leftCentralNode)
            for zdispl in u3l.values:
                mm.write(str(zdispl.nodeLabel))
                mm.write(',')
                mm.write(str(zdispl.dataDouble[2]))
                mm.write('\n')
            u3r=frame.fieldOutputs['U'].getSubset(region=rightCentralNode)
            for zdispr in u3r.values:
                nn.write(str(zdispr.nodeLabel))
                nn.write(',')
                nn.write(str(zdispr.dataDouble[2]))
                nn.write('\n')
    mm.close()
    nn.close()
    session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=-1 )
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
    ##session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
    session.graphicsOptions.setValues(backgroundStyle=SOLID, backgroundColor='#FFFFFF')
    session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(annotations=OFF)
    ##session.viewports['Viewport: 1'].view.setValues(nearPlane=19.801, farPlane=32.2696, width=17.4624, height=8.54503, viewOffsetX=2.44495, viewOffsetY=-1.97831)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=21.178, farPlane=32.1762, width=18.3706, height=8.80315, viewOffsetX=1.37452, viewOffsetY=-1.68468)
    session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(visibleEdges=FREE)
    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U3'), )
    session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(triadPosition=(8, 50))
    session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(title=OFF, state=OFF, compass=OFF)
    session.printOptions.setValues(vpDecorations=OFF)
    session.printToFile('Kirigami_Disp'+str(m), format=PNG, canvasObjects=(session.viewports['Viewport: 1'], ))

    #Create Displacement Variable
    coords = firstFrame.fieldOutputs['COORD']
    disp = lastFrame.fieldOutputs['U']
    mnodes = myOdb.rootAssembly.nodeSets['MIDNODES']
    lnodes = myOdb.rootAssembly.nodeSets['LEFTNODES']
    rnodes = myOdb.rootAssembly.nodeSets['RIGHTNODES']
    CoordsCenter=coords.getSubset(region=mnodes)
    CoordsLeft=coords.getSubset(region=lnodes)
    CoordsRight=coords.getSubset(region=rnodes)
    DisplacementCenter=disp.getSubset(region=mnodes)
    DisplacementLeft=disp.getSubset(region=lnodes)
    DisplacementRight=disp.getSubset(region=rnodes)
    for u in DisplacementCenter.values:
        f.write(str(u.nodeLabel))
        f.write(',')
        f.write(str(u.dataDouble[0]))
        f.write(',')
        f.write(str(u.dataDouble[1]))
        f.write(',')
        f.write(str(u.dataDouble[2]))
        f.write('\n')
    f.close()
    for v in CoordsCenter.values:
        g.write(str(v.nodeLabel))
        g.write(',')
        g.write(str(v.dataDouble[0]))
        g.write(',')
        g.write(str(v.dataDouble[1]))
        g.write(',')
        g.write(str(v.dataDouble[2]))
        g.write('\n')
    g.close()
    for w in DisplacementLeft.values:
        h.write(str(w.nodeLabel))
        h.write(',')
        h.write(str(w.dataDouble[0]))
        h.write(',')
        h.write(str(w.dataDouble[1]))
        h.write(',')
        h.write(str(w.dataDouble[2]))
        h.write('\n')
    h.close()
    for x in CoordsLeft.values:
        i.write(str(x.nodeLabel))
        i.write(',')
        i.write(str(x.dataDouble[0]))
        i.write(',')
        i.write(str(x.dataDouble[1]))
        i.write(',')
        i.write(str(x.dataDouble[2]))
        i.write('\n')
    i.close()
    for y in DisplacementRight.values:
        j.write(str(y.nodeLabel))
        j.write(',')
        j.write(str(y.dataDouble[0]))
        j.write(',')
        j.write(str(y.dataDouble[1]))
        j.write(',')
        j.write(str(y.dataDouble[2]))
        j.write('\n')
    j.close()
    for z in CoordsRight.values:
        k.write(str(z.nodeLabel))
        k.write(',')
        k.write(str(z.dataDouble[0]))
        k.write(',')
        k.write(str(z.dataDouble[1]))
        k.write(',')
        k.write(str(z.dataDouble[2]))
        k.write('\n')
    k.close()
    myOdb.close()
##Merging of files and deletion
    filenames = ['Kirigami_DispMidNodes'+str(m)+'.txt', 'Kirigami_CoorMidNodes'+str(m)+'.txt']
    with open('Kirigami_MidNodes'+str(m)+'.txt', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                outfile.write(infile.read())

    os.remove('Kirigami_DispMidNodes'+str(m)+'.txt')
    os.remove('Kirigami_CoorMidNodes'+str(m)+'.txt')

    filenames = ['Kirigami_DispLeftNodes'+str(m)+'.txt', 'Kirigami_CoorLeftNodes'+str(m)+'.txt']
    with open('Kirigami_LeftNodes'+str(m)+'.txt', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                outfile.write(infile.read())
    os.remove('Kirigami_DispLeftNodes'+str(m)+'.txt')
    os.remove('Kirigami_CoorLeftNodes'+str(m)+'.txt')

    filenames = ['Kirigami_DispRightNodes'+str(m)+'.txt', 'Kirigami_CoorRightNodes'+str(m)+'.txt']
    with open('Kirigami_RightNodes'+str(m)+'.txt', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                outfile.write(infile.read())
    os.remove('Kirigami_DispRightNodes'+str(m)+'.txt')
    os.remove('Kirigami_CoorRightNodes'+str(m)+'.txt')