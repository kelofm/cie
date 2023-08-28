import numpy
from scipy import linalg
from matplotlib import pyplot
from matplotlib.widgets import Slider


resolution = 15
dimension  = 3




localPoints = [
    [-1, -1],
    [1, -1],
    [-1, 1],
    [1, 1]
]

matrix = numpy.array([
    [1, -1, -1, 1],
    [1, 1, -1, -1],
    [1, -1, 1, -1],
    [1, 1, 1, 1]
])

factorization = linalg.lu_factor(matrix)



class Transform:
    def __init__(self, coefficients):
        self._basis = [
            self.makeBasis(coefficients[0]),
            self.makeBasis(coefficients[1])
        ]

    @staticmethod
    def makeBasis(coefficients):
        return [
            lambda u: coefficients[0],
            lambda u: coefficients[1] * u[0],
            lambda u: coefficients[2] * u[1],
            lambda u: coefficients[3] * u[0] * u[1]
        ]

    @staticmethod
    def getDefaultCoefficients():
        return [
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0]
        ]

#    @staticmethod
#    def makeBasis(c):
#        return [
#            lambda u: c[0] * c[2],
#            lambda u: c[1] * c[2] * u[0],
#            lambda u: c[0] * c[3] * u[1],
#            lambda u: c[1] * c[3] * u[0] * u[1]
#        ]
#
#    @staticmethod
#    def getDefaultCoefficients():
#        return [
#            [0.0, 1.0, 1.0, 0.0],
#            [1.0, 0.0, 0.0, 1.0]
#        ]

    def __call__(self, u):
        return numpy.array([
            sum(term(u) for term in basis) for basis in self._basis
        ])


def convertToSurfaceMesh(vertices, size0, size1):
    geometry = {'vertices':vertices, 'faces':[]}
    geometry['faces'] = numpy.zeros( ( 2*(size0-1)*(size1-1), 3 ), dtype=numpy.uint32 )
    k = 0
    for j in range(size0-1):
        for i in range(size1-1):
            jni = j*size1+i
            geometry['faces'][k] = [ 
                jni, 
                jni+size1,
                jni+size1+1
                ]
            geometry['faces'][k+1] = [ 
                jni, 
                jni+1+size1,
                jni+1
                ]
            k+=2
    return geometry

def getMesh(coefficients):
    localSpace = numpy.array([-1.0, 1.0])
    localGrid  = numpy.linspace(localSpace[0], localSpace[1], num=resolution)
    discretization = []
    for x in localGrid:
        for y in localGrid:
            discretization.append([x,y])

    
    transform = Transform(coefficients)
    transformed = [transform(point) for point in discretization]

    mesh = convertToSurfaceMesh(transformed, resolution, resolution)
    vertices = mesh["vertices"]

    coordinates = numpy.transpose(transformed)

    globalCoordinates = numpy.transpose([transform(point) for point in localPoints], (1,0))
#    xCoefficients = linalg.lu_solve(factorization, globalCoordinates[0])
#    yCoefficients = linalg.lu_solve(factorization, globalCoordinates[1])
#    
#    xPolynomial = numpy.polynomial.Polynomial(xCoefficients)
#    print(xPolynomial.roots())
    

    return coordinates, mesh["faces"]


coefficients = Transform.getDefaultCoefficients()

domain = [[-2.0, 2.0],[-2.0, 2.0]]

figure, axes = pyplot.subplots()
coordinates, faces = getMesh(coefficients)
lines = pyplot.triplot(coordinates[0], coordinates[1], faces)
axes.set_xlim(domain[0])
axes.set_ylim(domain[1])
pyplot.subplots_adjust(left=0.25, bottom=0.25)




coefficientNames = ["x", "y"]
coefDomain = [-2, 2]


widgets = {
    "sliders" : [[], []],
    "axes"    : [[], []],
    "orientation" : ["horizontal", "vertical"],
    "offset"      : [[0.2, 0.02], [0.02, 0.2]],
    "dOffset"     : [[0.0, 0.04], [0.04, 0.0]],
    "width"       : [0.6, 0.03],
    "height"      : [0.03, 0.6]
}

for i_component in range(len(coefficients)):
    width = widgets["width"][i_component]
    height = widgets["height"][i_component]
    dOffset = widgets["dOffset"][i_component]
    offset0 = widgets["offset"][i_component]

    for i_coef in range(len(coefficients[0])):
        localName = coefficientNames[i_component] + str(i_coef)

        offset = [offset0[0] + i_coef * dOffset[0], offset0[1] + i_coef * dOffset[1]]
        diagonal = [offset[0] + width, offset[1] + height]

        widgets["axes"][i_component].append(pyplot.axes([offset[0], offset[1], width, height]))
        widgets["sliders"][i_component].append(Slider(
            ax=widgets["axes"][i_component][-1],
            label=localName,
            valmin=coefDomain[0],
            valmax=coefDomain[1],
            valinit=coefficients[i_component][i_coef],
            orientation=widgets["orientation"][i_component]
        ))


def update(input):
    coefs = [[slider.val for slider in sliders] for sliders in widgets["sliders"]]
    coordinates, faces = getMesh(coefs)
    axes.clear()
    axes.triplot(coordinates[0], coordinates[1], faces)
    axes.set_xlim(domain[0])
    axes.set_ylim(domain[1])


for sliders in widgets["sliders"]:
    for slider in sliders:
        slider.on_changed(update)


pyplot.show()