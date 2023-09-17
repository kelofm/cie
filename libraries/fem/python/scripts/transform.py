""" @brief Visualize spatial transformation types implemented in CiE."""

# --- External Imports ---
from matplotlib import pyplot
import numpy

# --- CiE Imports ---
import cie.fem

# --- STD Imports ---
import typing


class ImageLoader:

    def __init__(self, imageURL: str):
        import PIL.Image
        import urllib.request
        import io

        response = urllib.request.urlopen(imageURL)
        self.__image = PIL.Image.open(io.BytesIO(response.read())).rotate(-90.0).convert("RGB")
        self.__image = self.__image.resize((max(self.__image.size),
                                            max(self.__image.size)))


    @property
    def image(self):
        return self.__image



class InteractiveTransform:

    def __init__(self,
                 transformType = cie.fem.maths.ProjectiveTransform2D,
                 valueSource: typing.Optional[str] = None):
        # Interactive state
        self.__i_selectedVertex: typing.Optional[int] = None

        # Constants
        self.__clickTolerance = 1e-1

        # Construct transformation
        self.__transformedBasis = numpy.array([[-1.0, -1.0],
                                               [ 1.0, -1.0],
                                               [-1.0,  1.0],
                                               [ 1.0,  1.0]]).transpose((1, 0))
        self.__transformType = transformType
        self.__transform = transformType()
        self.__inverseTransform = self.__transform.makeInverse()
        self.__updateTransform()

        # Render parameters
        self.__gridDensity = 50
        self.__gridRadius = 1

        # Local space
        self.__xGrid, self.__yGrid = numpy.meshgrid(
            numpy.linspace(0 - self.__gridRadius,
                           0 + self.__gridRadius,
                           int(2 * self.__gridRadius * self.__gridDensity)),
            numpy.linspace(0 - self.__gridRadius,
                           0 + self.__gridRadius,
                           int(2 * self.__gridRadius * self.__gridDensity)))

        self.__colors: numpy.ndarray
        if valueSource:
            try:
                import scipy.ndimage
                for source in (valueSource, "https://media.tenor.com/nfav07DxamYAAAAC/nicolas-cage-face.gif", None):
                    valueSource = source
                    try:
                        image = ImageLoader(valueSource).image
                        xSamples, ySamples, colorSamples = numpy.meshgrid(
                            numpy.linspace(0,
                                        image.size[1],
                                        int(2 * self.__gridRadius * self.__gridDensity - 1)),
                            numpy.linspace(image.size[1] / 4,
                                        3 * image.size[1] / 4,
                                        int(2 * self.__gridRadius * self.__gridDensity - 1)),
                            numpy.arange(3)
                        )
                        self.__colors = scipy.ndimage.map_coordinates(
                            image,
                            numpy.stack((xSamples, ySamples, colorSamples))
                        ) / 255.0
                        break
                    except:
                        continue
            except ImportError as exception: # <== missing packages, fall back to boring colormap
                valueSource = None

        if not valueSource:
            self.__colors = sum(numpy.meshgrid(
                numpy.linspace(0 - self.__gridRadius,
                               0 + self.__gridRadius,
                               int(2 * self.__gridRadius * self.__gridDensity - 1)),
                numpy.linspace(0 - self.__gridRadius,
                               0 + self.__gridRadius,
                               int(2 * self.__gridRadius * self.__gridDensity - 1))
            ))

        # Gridlines in local space
        self.__gridCoordinates = numpy.linspace(0 - self.__gridRadius, 0 + self.__gridRadius, 5)
        self.__localGrid = numpy.ndarray((2*len(self.__gridCoordinates), 2, self.__gridDensity))

        # Transformed grid
        self.__transformedXGrid = numpy.ndarray(self.__xGrid.shape)
        self.__transformedYGrid = numpy.ndarray(self.__xGrid.shape)
        self.__inverseXGrid = numpy.ndarray(self.__xGrid.shape)
        self.__inverseYGrid = numpy.ndarray(self.__xGrid.shape)
        self.__updateGrid()

        # Gridlines in transformed space
        self.__transformedGrid = numpy.ndarray((2*len(self.__gridCoordinates), 2, self.__gridDensity))
        self.__inverseGrid = numpy.ndarray((2*len(self.__gridCoordinates), 2, self.__gridDensity))
        self.__updateGridLines()

        # Static graphics
        self.__figure, (self.__localAxes, self.__transformedAxes, self.__inverseAxes) = pyplot.subplots(1, 3)
        self.__localColorMesh = self.__localAxes.pcolorfast(
            self.__xGrid,
            self.__yGrid,
            self.__colors)

        self.__localGridLines = []
        for i in range(0, self.__localGrid.shape[0]):
            self.__localGridLines.append(self.__localAxes.plot(
                self.__localGrid[i,0,:],
                self.__localGrid[i,1,:],
                color = [0.5, 0.5, 0.5, 1.0],
                linestyle = "-."
            ))

        self.__localAxes.set_xlim(2.5 * numpy.array([-self.__gridRadius, self.__gridRadius]))
        self.__localAxes.set_ylim(2.5 * numpy.array([-self.__gridRadius, self.__gridRadius]))
        self.__localAxes.set_title("Local Space")
        self.__localAxes.set_aspect("equal", "box")

        # Dynamic graphics
        self.__transformedColorMesh = self.__transformedAxes.pcolorfast(
            self.__transformedXGrid,
            self.__transformedYGrid,
            self.__colors)
        self.__inverseColorMesh = self.__inverseAxes.pcolorfast(
            self.__inverseXGrid,
            self.__inverseYGrid,
            self.__colors)

        self.__transformedGridLines = []
        for i in range(0, self.__transformedGrid.shape[0]):
            self.__transformedGridLines.append(self.__transformedAxes.plot(
                self.__transformedGrid[i,0,:],
                self.__transformedGrid[i,1,:],
                color = [0.5, 0.5, 0.5, 1.0],
                linestyle = "-."))
        self.__inverseGridLines = []
        for i in range(0, self.__inverseGrid.shape[0]):
            self.__inverseGridLines.append(self.__inverseAxes.plot(
                self.__inverseGrid[i,0,:],
                self.__inverseGrid[i,1,:],
                color = [0.5, 0.5, 0.5, 1.0],
                linestyle = "-."))

        # Interactive graphics
        controlPoints = self.__getControlPoints()
        self.__transformedBasisVertices = self.__transformedAxes.scatter(
            controlPoints[:,0],
            controlPoints[:,1],
            color = [1.0, 0.0, 0.0],
            zorder = 2
        )

        for axes in (self.__transformedAxes, self.__inverseAxes):
            axes.set_xlim(2.5 * numpy.array([-self.__gridRadius, self.__gridRadius]))
            axes.set_ylim(2.5 * numpy.array([-self.__gridRadius, self.__gridRadius]))
            axes.set_aspect("equal", "box")

        self.__transformedAxes.set_title("Transformed Space")
        self.__inverseAxes.set_title("Inverse Space")

        # Register callbacks
        canvas = self.__figure.canvas
        canvas.mpl_connect("button_press_event", self.__onButtonPress)
        canvas.mpl_connect("button_release_event", self.__onButtonRelease)
        canvas.mpl_connect("motion_notify_event", self.__onMouseMovement)


    def run(self) -> None:
        self.__figure.tight_layout()
        pyplot.show()


    def __draw(self) -> None:
        # Dynamic graphics
        self.__transformedColorMesh.remove() # <== I just can't believe the mesh is immutable @todo
        self.__transformedColorMesh = self.__transformedAxes.pcolorfast(
            self.__transformedXGrid,
            self.__transformedYGrid,
            self.__colors,
            zorder = 0)

        self.__inverseColorMesh.remove() # <== I just can't believe the mesh is immutable @todo
        self.__inverseColorMesh = self.__inverseAxes.pcolorfast(
            self.__inverseXGrid,
            self.__inverseYGrid,
            self.__colors,
            zorder = 0)

        for i in range(0, self.__transformedGrid.shape[0]):
            self.__transformedGridLines[i][0].set_xdata(self.__transformedGrid[i,0,:])
            self.__transformedGridLines[i][0].set_ydata(self.__transformedGrid[i,1,:])

        for i in range(0, self.__inverseGrid.shape[0]):
            self.__inverseGridLines[i][0].set_xdata(self.__inverseGrid[i,0,:])
            self.__inverseGridLines[i][0].set_ydata(self.__inverseGrid[i,1,:])

        # Interactive graphics
        self.__updateTransformedBasis()
        self.__figure.canvas.draw()


    def __onButtonPress(self, event) -> None:
        selectedIndex = self.__getClickedTransformedBasisIndex(event)
        if selectedIndex is not None:
            self.__i_selectedVertex = selectedIndex


    def __onButtonRelease(self, event) -> None:
        if self.__i_selectedVertex is not None:
            self.__i_selectedVertex = None


    def __onMouseMovement(self, event) -> None:
        if self.__i_selectedVertex is not None:
            self.__transformedBasis[:,self.__i_selectedVertex] = [event.xdata, event.ydata]
            self.__updateTransform()
            self.__updateGrid()
            self.__updateGridLines()
            self.__draw()


    def __updateTransform(self) -> None:
        """ @brief Update the transfromation to project the transformed basis."""
        controlPoints = self.__getControlPoints()
        self.__transform = self.__transformType(controlPoints)
        self.__inverseTransform = self.__transform.makeInverse()


    def __updateGrid(self) -> None:
        """ @brief Update mesh in transformed space."""
        for i_row in range(self.__xGrid.shape[0]):
            for i_column in range(self.__xGrid.shape[1]):
                transformedPoint = self.__transform.evaluate((self.__xGrid[i_row, i_column],
                                                              self.__yGrid[i_row, i_column]))
                self.__transformedXGrid[i_row, i_column] = transformedPoint[0]
                self.__transformedYGrid[i_row, i_column] = transformedPoint[1]

                inversePoint = self.__inverseTransform.evaluate((self.__xGrid[i_row, i_column],
                                                                 self.__yGrid[i_row, i_column]))
                self.__inverseXGrid[i_row, i_column] = inversePoint[0]
                self.__inverseYGrid[i_row, i_column] = inversePoint[1]


    def __updateGridLines(self) -> None:
        """ @brief Update gridlines in transformed space."""
        for i_main, maincoordinate in enumerate(self.__gridCoordinates):
            for i_sub, subcoordinate in enumerate(numpy.linspace(self.__gridCoordinates[0], self.__gridCoordinates[-1], self.__gridDensity)):
                self.__localGrid[2*i_main,:,i_sub] = (maincoordinate, subcoordinate)
                self.__localGrid[2*i_main+1,:,i_sub] = (subcoordinate, maincoordinate)
                self.__transformedGrid[2*i_main,:,i_sub]   = self.__transform.evaluate((maincoordinate, subcoordinate))
                self.__transformedGrid[2*i_main+1,:,i_sub] = self.__transform.evaluate((subcoordinate, maincoordinate))
                self.__inverseGrid[2*i_main,:,i_sub]   = self.__inverseTransform.evaluate((maincoordinate, subcoordinate))
                self.__inverseGrid[2*i_main+1,:,i_sub] = self.__inverseTransform.evaluate((subcoordinate, maincoordinate))


    def __mapBasisIndex(self, basisIndex: int) -> typing.Optional[int]:
        """ @brief Return a map as"""
        map: dict[int,int]
        if isinstance(self.__transform, cie.fem.maths.OrthogonalScaleTransform2D):
            map = {
                0 : None,
                1 : None,
                2 : None,
                3 : 0
            }
        elif isinstance(self.__transform, cie.fem.maths.ScaleTranslateTransform2D):
            map = {
                0 : 0,
                1 : None,
                2 : None,
                3 : 1
            }
        elif isinstance(self.__transform, cie.fem.maths.TranslateScaleTransform2D):
            map = {
                0 : 0,
                1 : None,
                2 : None,
                3 : 1
            }
        elif isinstance(self.__transform, cie.fem.maths.AffineTransform2D):
            map = {
                0 : 0,
                1 : 1,
                2 : 2,
                3 : None
            }
        elif isinstance(self.__transform, cie.fem.maths.ProjectiveTransform2D):
            map = {
                0 : 0,
                1 : 1,
                2 : 2,
                3 : 3
            }
        else:
            raise TypeError(f"Unsupported transform type: {type(self.__transform)}")
        return map[basisIndex]


    def __getControlPoints(self) -> numpy.ndarray:
        indexMap = [i for i, j in enumerate(self.__mapBasisIndex(i) for i in range(self.__transformedBasis.shape[1])) if j is not None]
        return self.__transformedBasis[:,indexMap].transpose()


    #def __getClickedTransformedBasisIndex(self, event) -> typing.Optional[int]:
    #    controlPoints = self.__getControlPoints()
    #    for i in range(controlPoints.shape[0]):
    #        distance = numpy.linalg.norm(controlPoints[i,:] - [event.xdata, event.ydata], ord=1)
    #        if distance < self.__clickTolerance:
    #            return i
    #    return None
    def __getClickedTransformedBasisIndex(self, event) -> typing.Optional[int]:
        for i in range(self.__transformedBasis.shape[1]):
            distance = abs(self.__transformedBasis[0,i] - event.xdata) + abs(self.__transformedBasis[1,i] - event.ydata)
            if distance < self.__clickTolerance:
                if self.__mapBasisIndex(i) is not None:
                    return i
        return None


    def __updateTransformedBasis(self) -> None:
        self.__transformedBasisVertices.set_offsets(self.__getControlPoints())



if __name__ == "__main__":
    # --- STD Imports ---
    import argparse
    import sys

    parser = argparse.ArgumentParser("transform",
                                     description = "Spatial transform demo")
    parser.add_argument("-t",
                        "--transform-type",
                        type = str,
                        dest = "transformType",
                        choices = ["scale", "scale-translate", "translate-scale", "affine", "projective"],
                        default = "projective")
    parser.add_argument("-i",
                        "--image-url",
                        type = str,
                        dest = "imageURL",
                        default = "")

    try:
        arguments = parser.parse_args(sys.argv[1:])
    except Exception as exception:
        print(exception)
        exit(1)

    InteractiveTransform(
        transformType = {
            "scale"             : cie.fem.maths.OrthogonalScaleTransform2D,
            "scale-translate"   : cie.fem.maths.ScaleTranslateTransform2D,
            "translate-scale"   : cie.fem.maths.TranslateScaleTransform2D,
            "affine"            : cie.fem.maths.AffineTransform2D,
            "projective"        : cie.fem.maths.ProjectiveTransform2D
        }[arguments.transformType],
        valueSource = arguments.imageURL
    ).run()
    exit(0)
