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
Ncores4Sim = 4

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
limTop = 1
# Model Construction for each set of parameters
for m in range(limInf,limTop):


################################ resonance eigenvalue analysis

    namefile="Kirigami_Resonance_"+str(m)
    ## Some shortcuts
    executeOnCaeStartup()
    mdb.Model(modelType=STANDARD_EXPLICIT, name=namefile)
    s = mdb.models[namefile].ConstrainedSketch(name='__profile__', sheetSize=50)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    print(m)

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
    # Other material properties for damping to be added here
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
        print(i)
        pickedRegions = f.getByBoundingBox(-1., -1., -1., length+1, width+1, 1.) 
        datumid = planes[i]
        p.PartitionFaceByDatumPlane(pickedRegions,p.datums[datumid])

    ##Create Assembly
    myAssembly = mdb.models[namefile].rootAssembly
    myInstance = myAssembly.Instance(name='Kirigami',
        part=p, dependent=ON)
    ##Create RF Point and RF Set
    myAssembly.ReferencePoint(point=(length,width/2,0))
    myAssembly.Set(referencePoints=((myAssembly.referencePoints.findAt((length,width/2,0),),)), name='setrpmaster')
    ##  Create geometry set 'FixedEnd', 'ActuationEnd', 'All Faces'
    e1 = myAssembly.instances['Kirigami'].edges
    f1 = myAssembly.instances['Kirigami'].faces
    myAssembly.Set(edges=e1.getByBoundingBox(-0.01, -1., -1., 0.01, width+0.01, 1.), name='FixedEnd')
    myAssembly.Set(edges=e1.getByBoundingBox(length-0.01, -1., -1., length+0.01, width+0.01, 1.), name='ActuationEnd')
    myAssembly.Set(faces=f1.getByBoundingBox(-1, -1., -1., length+0.01, width+0.01, 1.), name='AllFaces')

    ##Create Step
    mdb.models[namefile].StaticStep(name='ResonanceStaticLoad',
        previous='Initial',
        description='Static Loads For Resonance Anlaysis')
    mdb.models[namefile].FrequencyStep(name='ResonanceEig',
        previous='ResonanceStaticLoad',
        description='Frequency Eigenvalue Analysis',
        numEigen=10, eigensolver=LANCZOS, maxEigen=20000.0)

    ##Create Interaction
    regionRF = myAssembly.sets['setrpmaster']
    regionSlave = myAssembly.Set(edges=e1.getByBoundingBox(length-0.01, -1., -1., length+0.01, width+0.01, 1.), name='ActuationEnd')
    mdb.models[namefile].RigidBody(name='RBE2',refPointRegion=regionRF,tieRegion=regionSlave)

    ##Set Output Request
    # Outputs stress and displacement
    mdb.models[namefile].fieldOutputRequests['F-Output-1'].setValues(variables=(
        'S', 'U'))

    ##Apply Boundary Conditions
    region1=myAssembly.sets['FixedEnd']
    region2=myAssembly.sets['ActuationEnd']
    mdb.models[namefile].EncastreBC(name='Fixed End', createStepName='ResonanceEig', region=region1)
    mdb.models[namefile].DisplacementBC(name='Actuation End', createStepName='ResonanceEig', region=regionRF, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0)
    mdb.models[namefile].EncastreBC(name='Fixed End', createStepName='ResonanceStaticLoad', region=region1)
    mdb.models[namefile].DisplacementBC(name='Actuation End', createStepName='ResonanceStaticLoad', region=regionRF, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0)
    mdb.models[namefile].ConcentratedForce(name='CF', createStepName='ResonanceStaticLoad', region=regionRF, cf1=1.0, cf2=0.0, cf3=0.0, distributionType=UNIFORM, field='', localCsys=None)

    ##Mesh
    # quadrilateral (S4R) shell elements where possible, tringular shell elements wherever not possible
    elemType1 = mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=S3,  elemLibrary=STANDARD)
    p.seedPart(size=MeshSize)
    faces1 = p.faces
    pickedRegions =(faces1, )
    ##p.setMeshControls(regions=pickedRegions, technique=STRUCTURED)
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
    p.generateMesh()

    #Define Output NodeSets
    myAssembly.regenerate()
    n1=myAssembly.instances['Kirigami'].nodes
    myAssembly.Set(nodes=n1.getByBoundingBox(deltax1[m]+t1[m]/2, width/2-0.01, -1., length-deltax2[m]-t2[m]/2, width/2+0.01, 1.), name='MidNodes')

    ##  Create job
    JobNameBuckle='Kirigami_Resonance_'+str(m)
    mdb.Job(name=JobNameBuckle, model=namefile, numCpus=Ncores4Sim, numDomains=8, nodalOutputPrecision=FULL)
    myAssembly.regenerate()









############################ jobs construction and post-processing

##End of Model Construction - Run
## Submit the job (Run Buckling Analysis Just Once)
for m in range(0,1):
    mdb.jobs['Kirigami_Resonance_'+str(m)].submit(consistencyChecking=OFF)
    mdb.jobs['Kirigami_Resonance_'+str(m)].waitForCompletion()
