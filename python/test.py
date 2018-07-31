import jpype


jvmPath = jpype.getDefaultJVMPath() 
jpype.startJVM(jvmPath, "-ea")
#, "-Djava.class.path=/home/mendozasmith/src/javaplex/src/matlab/for_distribution/lib/javaplex.jar")

java.lang.System.out.println("Hello")
sys.exit()
pck = jpype.JPackage('edu').stanford.math.plex4.examples

pck = jpype.JPackage('edu').stanford.math.plex4.examples
pce = pck.PointCloudExamples
a = pce.getHouseExample
print(dir(a))

sys.exit()
pcClass = jpype.JClass("edu.stanford.math.plex4.examples.PointCloudExamples")
obj = pcClass()
print(obj)
h = obj.getHouseExample()
print(h)
sys.exit()

f = jpype.JPackage("edu").stanford.math.plex4.examples.PointCloudExamples
print(f)

#"edu/stanford/math/plex4/api/Plex4.class"


api = JClass("edu.stanford.math.plex4.api")
print(api)

#a = Plex4.createVietorisRipsStream(point_cloud, 2,3,1)

jvmPath = jpype.getDefaultJVMPath() 
#jpype.startJVM(jvmPath, "-Djava.class.path=/home/di/eclipse_plugins/plugins/*.jar")
jpype.startJVM(jvmPath, "-Djava.class.path=/home/mendozasmith/src/javaplex/src/matlab/for_distribution/lib/javaplex.jar")

pointCloudExamples = JPackage("edu").stanford.math.plex4.examples.PointCloudExamples
point_cloud = pointCloudExamples.getHouseExample

s = JClass('edu.stanford.math.plex4.api.Plex4')
print(s)
sys.exit()
a = s.createVietorisRipsStream(point_cloud, 2,3,1)
jpype.shutdownJVM() 

crv = plex4.createVietorisRipsStream
stream = crv(point_cloud, 2, 3, 2)
print(stream)

