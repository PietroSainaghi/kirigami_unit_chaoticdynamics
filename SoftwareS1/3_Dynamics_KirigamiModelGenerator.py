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

# parameters for computer solver
Ncores4Sim = 16

#Geometric Parameters all in [mm] - Should be read from a Matlab Script that produces the random numbers
#length=10.55 Defined abajo
lengtharm=6.25
width=10
thick=0.1524
EYoung = 2.76e3
Poisson = 0.34
Dens = 1.42e-9
DampAlpha = 216.1073
DampBeta = 1e-9
nLoadSteps = 2
TotDisp = 1
MeshSize=0.0625

#Initialize empty arrays
w1=[]
w2=[]
t1=[]
t2=[]
deltax1=[]
deltax2=[]
deltay1=[]
deltay2=[]
deltay3=[]
deltay4=[]
a1=[]
a2=[]
a3=[]
a4=[]
l1=[]
l2=[]
PCr=[]
alfa=[]
freqHz=[]

# Open CSV file that contains geometric parameters, and import them
## CHANGE NAME OF FILE
with open('KirigamiFrequencyData_Damped.csv') as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
        w1.append(float(row[0]))
        w2.append(float(row[1]))
        t1.append(float(row[2]))
        t2.append(float(row[3]))
        deltax1.append(float(row[4]))
        deltax2.append(float(row[5]))
        deltay1.append(float(row[6]))
        deltay2.append(float(row[7]))
        deltay3.append(float(row[8]))
        deltay4.append(float(row[9]))
        a1.append(float(row[10]))
        a2.append(float(row[11]))
        a3.append(float(row[12]))
        a4.append(float(row[13]))
        l1.append(float(row[14]))
        l2.append(float(row[15]))
        PCr.append(float(row[17]))
        alfa.append(float(row[18]))
        freqHz.append(float(row[19]))

# Set what rows from the CSV to analyze
# If zero is not present, the linear buckling analysis run (job submitted in line 465) will not be used by the nonlinear analysis, and the code will throw an error if 'Kirigami_Buckle_0' is not present
limInf = 0
limTop = 10
# Model Construction for each set of parameters
for m in range(limInf,limTop):


################################## Nonlinear dynamics frequency analysis
    namefile="Kirigami_NLGeom_"+str(m)
    ## Some shortcuts
    executeOnCaeStartup()
    mdb.Model(modelType=STANDARD_EXPLICIT, name=namefile)
    s = mdb.models[namefile].ConstrainedSketch(name='__profile__', sheetSize=50)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)

    ##Create Geometry
    length=deltax1[m]+lengtharm+a2[m]+w2[m]+a3[m]+lengtharm+deltax1[m]
    s.Line(point1=(0,0), point2=(0,width))
    s.Line(point1=(0,width), point2=(deltax1[m]+lengtharm+a2[m],width))
    s.Line(point1=(deltax1[m]+lengtharm+a2[m],width),point2=(deltax1[m]+lengtharm+a2[m],width-l2[m]))
    s.Line(point1=(deltax1[m]+lengtharm+a2[m],width-l2[m]),point2=(deltax1[m]+lengtharm+a2[m]+w2[m],width-l2[m]))
    s.Line(point1=(deltax1[m]+lengtharm+a2[m]+w2[m],width-l2[m]),point2=(deltax1[m]+lengtharm+a2[m]+w2[m],width))
    s.Line(point1=(deltax1[m]+lengtharm+a2[m]+w2[m],width),point2=(length,width))
    s.Line(point1=(length,width),point2=(length,0))
    s.Line(point1=(length,0),point2=(length-deltax2[m]-lengtharm-a3[m],0))
    s.Line(point1=(length-deltax2[m]-lengtharm-a3[m],0),point2=(length-deltax2[m]-lengtharm-a3[m],l1[m]))
    s.Line(point1=(length-deltax2[m]-lengtharm-a3[m],l1[m]),point2=(length-deltax2[m]-lengtharm-a3[m]-w1[m],l1[m]))
    s.Line(point1=(length-deltax2[m]-lengtharm-a3[m]-w1[m],l1[m]),point2=(length-deltax2[m]-lengtharm-a3[m]-w1[m],0))
    s.Line(point1=(length-deltax2[m]-lengtharm-a3[m]-w1[m],0),point2=(0,0))
    s.Line(point1=(deltax1[m]+lengtharm,deltay1[m]),point2=(deltax1[m],deltay1[m]))
    s.Line(point1=(deltax1[m],deltay1[m]),point2=(deltax1[m],width-deltay2[m]))
    s.Line(point1=(deltax1[m],width-deltay2[m]),point2=(deltax1[m]+lengtharm,width-deltay2[m]))
    s.Line(point1=(deltax1[m]+lengtharm,deltay1[m]),point2=(deltax1[m]+lengtharm,deltay1[m]+t1[m]))
    s.Line(point1=(deltax1[m]+lengtharm,deltay1[m]+t1[m]),point2=(deltax1[m]+t1[m],deltay1[m]+t1[m]))
    s.Line(point1=(deltax1[m]+t1[m],deltay1[m]+t1[m]),point2=(deltax1[m]+t1[m],width-deltay2[m]-t1[m]))
    s.Line(point1=(deltax1[m]+t1[m],width-deltay2[m]-t1[m]),point2=(deltax1[m]+lengtharm,width-deltay2[m]-t1[m]))
    s.Line(point1=(deltax1[m]+lengtharm,width-deltay2[m]-t1[m]),point2=(deltax1[m]+lengtharm,width-deltay2[m]))
    s.Line(point1=(deltax1[m]+lengtharm+a1[m]+w1[m]+a4[m],deltay4[m]),point2=(length-deltax2[m],deltay4[m]))
    s.Line(point1=(length-deltax2[m],deltay4[m]),point2=(length-deltax2[m],width-deltay3[m]))
    s.Line(point1=(length-deltax2[m],width-deltay3[m]),point2=(length-deltax2[m]-lengtharm,width-deltay3[m]))
    s.Line(point1=(length-deltax2[m]-lengtharm,width-deltay3[m]),point2=(length-deltax2[m]-lengtharm,width-deltay3[m]-t2[m]))
    s.Line(point1=(length-deltax2[m]-lengtharm,width-deltay3[m]-t2[m]),point2=(length-deltax2[m]-t2[m],width-deltay3[m]-t2[m]))
    s.Line(point1=(length-deltax2[m]-t2[m],width-deltay3[m]-t2[m]),point2=(length-deltax2[m]-t2[m],deltay4[m]+t2[m]))
    s.Line(point1=(length-deltax2[m]-t2[m],deltay4[m]+t2[m]),point2=(length-deltax2[m]-lengtharm,deltay4[m]+t2[m]))
    s.Line(point1=(length-deltax2[m]-lengtharm,deltay4[m]+t2[m]),point2=(deltax1[m]+lengtharm+a1[m]+w1[m]+a4[m],deltay4[m]))

    p = mdb.models[namefile].Part(name='Part-1', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models[namefile].parts['Part-1']
    p.BaseShell(sketch=s)
    s.unsetPrimaryObject()
    p = mdb.models[namefile].parts['Part-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models[namefile].sketches['__profile__']


    ##  Create material 'Material'
    mdb.models[namefile].Material('Material')
    mdb.models[namefile].materials['Material'].Elastic(table=((EYoung, Poisson), ))
    mdb.models[namefile].materials['Material'].Density(table=((Dens, ), ))
    mdb.models[namefile].materials['Material'].Damping(alpha=DampAlpha, beta=DampBeta)
    

    ##  Create shell section
    mdb.models[namefile].HomogeneousShellSection(name='Kirigami_Shell',
        preIntegrate=OFF, material='Material', thickness=thick,
        poissonDefinition=DEFAULT, temperature=GRADIENT,
        integrationRule=SIMPSON, numIntPts=5)
    p = mdb.models[namefile].parts['Part-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    f = p.faces
    faces = f
    region = regionToolset.Region(faces=faces)
    p.SectionAssignment(region=region, sectionName='Kirigami_Shell', offset=0.0)

    planes=[p.DatumPlaneByPrincipalPlane(XZPLANE,width/2).id, p.DatumPlaneByPrincipalPlane(YZPLANE,deltax1[m]+lengtharm+a1[m]+w1[m]/2).id, p.DatumPlaneByPrincipalPlane(YZPLANE,deltax1[m]).id, p.DatumPlaneByPrincipalPlane(YZPLANE,length-deltax2[m]).id,
            p.DatumPlaneByPrincipalPlane(YZPLANE,deltax1[m]+t1[m]).id,p.DatumPlaneByPrincipalPlane(YZPLANE,length-deltax2[m]-t2[m]).id,
            p.DatumPlaneByPrincipalPlane(YZPLANE,deltax1[m]+lengtharm+a1[m]).id, p.DatumPlaneByPrincipalPlane(YZPLANE,deltax1[m]+lengtharm+a1[m]+w1[m]).id,
            p.DatumPlaneByPrincipalPlane(XZPLANE,l1[m]).id, p.DatumPlaneByPrincipalPlane(XZPLANE,width-l2[m]).id, p.DatumPlaneByPrincipalPlane(YZPLANE,deltax1[m]+lengtharm).id, p.DatumPlaneByPrincipalPlane(YZPLANE,deltax1[m]+lengtharm+a1[m]+w1[m]+a4[m]).id, p.DatumPlaneByPrincipalPlane(XZPLANE,deltay1[m]).id, p.DatumPlaneByPrincipalPlane(XZPLANE,deltay1[m]+t1[m]).id, p.DatumPlaneByPrincipalPlane(XZPLANE,width-deltay2[m]).id, p.DatumPlaneByPrincipalPlane(XZPLANE,width-deltay2[m]-t1[m]).id]

    if (abs((deltax1[m]+lengtharm+a1[m]+w1[m]+a4[m])-(deltax1[m]+lengtharm+a2[m]+w2[m]+a3[m])))>1e-3:
       p.DatumPlaneByPrincipalPlane(YZPLANE,deltax1[m]+lengtharm+a2[m]+w2[m]+a3[m]).id
    
    ## Partitioning
    for i in range(0,len(planes)):
        pickedRegions = f.getByBoundingBox(-1., -1., -1., length+1, width+1, 1.) 
        datumid = planes[i]
        p.PartitionFaceByDatumPlane(pickedRegions,p.datums[datumid])

    ##Create Assembly
    myAssembly = mdb.models[namefile].rootAssembly
    myInstance = myAssembly.Instance(name='Kirigami',
        part=p, dependent=ON)

    ##  Create geometry set 'FixedEnd', 'ActuationEnd', 'All Faces'
    e1 = myAssembly.instances['Kirigami'].edges
    f1 = myAssembly.instances['Kirigami'].faces
    n1=myAssembly.instances['Kirigami'].nodes
    myAssembly.Set(edges=e1.getByBoundingBox(-0.01, -1., -1., 0.01, width+0.01, 1.), name='FixedEnd')
    myAssembly.Set(edges=e1.getByBoundingBox(length-0.01, -1., -1., length+0.01, width+0.01, 1.), name='ActuationEnd')
    myAssembly.Set(faces=f1.getByBoundingBox(-1, -1., -1., length+0.01, width+0.01, 1.), name='AllFaces')
    regionReac = myAssembly.sets['FixedEnd']

    ##Create RF Point and RF Set
    myAssembly.ReferencePoint(point=(length,width/2,0))
    myAssembly.Set(referencePoints=((myAssembly.referencePoints.findAt((length,width/2,0),),)), name='setrpmaster')
    
    ##Create Step
    for nstep in range(0, nLoadSteps):
        if nstep == 0:
            #mdb.models[namefile].StaticStep(name='Traction'+str(nstep),
            #    previous='Initial',
            #    description='Nonlinear analysis: TractionLoad',
            #    timePeriod=1, adiabatic=OFF, maxNumInc=1000000,
            #    stabilization=None, timeIncrementationMethod=AUTOMATIC, initialInc=1e-6,
            #    minInc=1e-20, maxInc=1e-3, matrixSolver=SOLVER_DEFAULT, amplitude=RAMP,
            #    extrapolation=LINEAR, fullyPlastic="", nlgeom=ON)
            mdb.models[namefile].StaticStep(name='Traction'+str(nstep),
                previous='Initial',
                description='Nonlinear analysis: TractionLoad',
                timePeriod=1,
                initialInc=1e-3,
                minInc=1e-20,
                maxInc=0.01,
                maxNumInc=10000,
                nlgeom=ON)
        else:
            tPer=(100*(1./freqHz[m]))	
            significant_digits = 2
            rounded_tPer =  round(tPer, significant_digits - int(math.floor(math.log10(abs(tPer)))) - 1)
            mdb.models[namefile].ImplicitDynamicsStep(name='Traction'+str(nstep),
                previous='Traction'+str(nstep-1),
                description='Dynamic Periodic Analysis TractionLoad',
                timePeriod=rounded_tPer, maxNumInc=10000, initialInc=1e-6,
                minInc=1e-10, maxInc=1e-4, nlgeom=ON)


    ##Create Interaction
    regionRF = myAssembly.sets['setrpmaster']
    regionSlave = myAssembly.Set(edges=e1.getByBoundingBox(length-0.01, -1., -1., length+0.01, width+0.01, 1.), name='ActuationEnd')
    mdb.models[namefile].RigidBody(name='RBE2',refPointRegion=regionRF,tieRegion=regionSlave)

    ##Set Output Request

    ##General
    ##del mdb.models[namefile].historyOutputRequests['H-Output-1']
    tInt=((1/30.)*(1./freqHz[m]))	
    significant_digits = 4
    rounded_tInt =  round(tInt, significant_digits - int(math.floor(math.log10(abs(tInt)))) - 1)
    mdb.models[namefile].fieldOutputRequests['F-Output-1'].setValuesInStep(stepName='Traction0', variables=(
        'RF', 'U', 'COORD'), frequency=LAST_INCREMENT)
    mdb.models[namefile].historyOutputRequests['H-Output-1'].setValuesInStep(stepName='Traction0', variables=(
        'ALLIE', 'ALLSE', 'ALLAE', 'ALLWK'), frequency=LAST_INCREMENT)
    mdb.models[namefile].fieldOutputRequests['F-Output-1'].setValuesInStep(stepName='Traction1', variables=(
        'RF', 'U', 'COORD'), timeInterval=rounded_tInt)
    mdb.models[namefile].historyOutputRequests['H-Output-1'].setValuesInStep(stepName='Traction1', variables=(
        'ALLIE', 'ALLSE', 'ALLAE', 'ALLWK'), frequency=LAST_INCREMENT)
   
    ##Define Amplitude
    PLoad=0.5*alfa[m]*PCr[m]
    mdb.models[namefile].PeriodicAmplitude(name='Amp-1', frequency=2*math.pi*freqHz[m], start=0.0, a_0=PLoad, data=((0.0, PLoad), ), timeSpan=TOTAL)

    ##Apply Boundary Conditions
    region1=myAssembly.sets['FixedEnd']
    region2=myAssembly.sets['ActuationEnd']
    mdb.models[namefile].EncastreBC(name='Fixed End', createStepName='Traction0', region=region1)
    mdb.models[namefile].DisplacementBC(name='Actuation End', createStepName='Traction0', region=region2, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0)
    ##Apply Load
    ##Prestep
    mdb.models[namefile].ConcentratedForce(name='CF', createStepName='Traction0', region=regionRF, cf1=PLoad, cf2=0.0, cf3=0.0, distributionType=UNIFORM, field='', localCsys=None)
    ##Second step
    mdb.models[namefile].loads['CF'].deactivate('Traction1')
    mdb.models[namefile].ConcentratedForce(name='CF_Time', createStepName='Traction1', region=regionRF, amplitude='Amp-1', cf1=1.0, cf2=0.0, cf3=0.0, distributionType=UNIFORM, field='', localCsys=None)

    ##Mesh
    elemType1 = mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=S3,  elemLibrary=STANDARD)
    p.seedPart(size=MeshSize)
    faces1 = p.faces
    pickedRegions =(faces1, )
    ##p.setMeshControls(regions=pickedRegions, technique=STRUCTURED)
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
    p.generateMesh()

    ##Edit Keyword to apply imperfections
    partname = "Part-1"
    model = mdb.models[namefile]
    modelkwb = model.keywordBlock
    modelkwb.synchVersions(storeNodesAndElements=False)
    line_num = 0
    for n, line in enumerate(modelkwb.sieBlocks):
        if line.replace(" ","").lower() == "*Restart,write,frequency=0".lower():
            line_num = n
            break
    if line_num:
	JobNameBuckle='Kirigami_Buckle_0'
        #kwds = '*IMPERFECTION, FILE='+ JobNameBuckle +', STEP=1 \n 1,0.001 \n 2,0.001 \n 3,0.001 \n 4,0.001 \n 5,0.001 \n'
        kwds = '*IMPERFECTION, FILE='+ JobNameBuckle +', STEP=1 \n 1,0.0001 \n'
        modelkwb.insert(position=line_num-11, text=kwds)
    else:
        e = ("Error: Part '{}' was not found".format(partname),
             "in the Model KeywordBlock.")
        raise Exception(" ".join(e))
    
    #Define Output NodeSets
    myAssembly.regenerate()
    n2=myAssembly.instances['Kirigami'].nodes
    myAssembly.Set(nodes=n2.getByBoundingBox(deltax1[m]+t1[m]/2, width/2-0.01, -1., length-deltax2[m]-t2[m]/2, width/2+0.01, 1.), name='MidNodes')
    myAssembly.Set(nodes=n2.getByBoundingBox(deltax1[m]+t1[m]-0.01, deltay1[m]+t1[m]-0.01, -1., deltax1[m]+t1[m]+0.01, width-deltay2[m]-t1[m]+0.01, 1.), name='LeftNodes')
    myAssembly.Set(nodes=n2.getByBoundingBox(length-deltax2[m]-t2[m]-0.01, deltay4[m]+t2[m]-0.01, -1., length-deltax2[m]-t2[m]+0.01, width-deltay4[m]-t2[m]+0.01, 1.), name='RightNodes')
    myAssembly.Set(nodes=n2.getByBoundingBox(length-0.01, -0.01, -1., length+0.01, 0.01, 1.), name='TrackNode')
    regionUtrack = myAssembly.sets['TrackNode']
    myAssembly.Set(nodes=n2.getByBoundingBox(deltax1[m]+t1[m]/2, width/2-0.01, -1., deltax1[m]+t1[m]+0.01, width/2+0.01, 1.), name='LeftCenterNode')
    regionUtrack = myAssembly.sets['LeftCenterNode']
    myAssembly.Set(nodes=n2.getByBoundingBox(length-deltax2[m]-t2[m]-0.01, width/2-0.01, -1., length-deltax2[m]-t2[m]/2, width/2+0.01, 1.), name='RightCenterNode')
    regionUtrack = myAssembly.sets['RightCenterNode']
    
    ##  Create job
    JobNameNLGeom='Kirigami_NLGeom_'+str(m)
    mdb.Job(name=JobNameNLGeom, model=namefile, numCpus=Ncores4Sim, numDomains=Ncores4Sim*2, nodalOutputPrecision=FULL)
    myAssembly.regenerate()






############################ jobs construction and post-processing

##End of Model Construction - Run
## Submit the job (Run Buckling Analysis Just Once)
#for m in range(0,1):
#    mdb.jobs['Kirigami_Buckle_'+str(m)].submit(consistencyChecking=OFF)
#    mdb.jobs['Kirigami_Buckle_'+str(m)].waitForCompletion()
# for m in range(0,1):
#     mdb.jobs['Kirigami_Resonance_'+str(m)].submit(consistencyChecking=OFF)
#     mdb.jobs['Kirigami_Resonance_'+str(m)].waitForCompletion()
for m in range(limInf,limTop):
    mdb.jobs['Kirigami_NLGeom_'+str(m)].submit(consistencyChecking=OFF)
    mdb.jobs['Kirigami_NLGeom_'+str(m)].waitForCompletion()

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