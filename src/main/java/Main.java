import controller.FresnelTransform;
import org.apache.commons.math3.complex.Complex;
import view.Graph;
import static controller.Superposition.*;

public class Main {

    public static void main(String[] args) {

        int width = 128, height = 128; // Размеры изображения.
        int N = 3; // Порядок полиномов Эрмита.
        int n = 3, m = 3; // Порядки моды для формулы 5.
        // u, v ∈ [-A; A]
        double uRange, vRange;
        uRange = vRange = 32;
        // x, y ∈ [-B; B]
        double xRange, yRange;
        xRange = yRange = 32;
        double z = 100; // Расстояние.
        double lambda = 0.00063; // Длина волны.
        double k = 4 /*2 * Math.PI / lambda*/; // Волновое число.
        double yMin, xMin; //Эти значения - границы для интеграла.
        yMin = xMin = -32;
        double stepY = -2 * yMin / height;
        double stepX = -2 * xMin / width;
        double gauss = 7; // Гауссовый параметр.
        Complex[][] function = new Complex[width][height];
        double stepU = 2 * uRange / width;
        double stepV = 2 * vRange / height;
        double xStep = 2 * xRange / width;
        double yStep = 2 * yRange / height;
        Complex[][] coefficient = new Complex[N + 1][N + 1]; //Набор комплексных коэффициентов.

        // Задаем им значения.
        for (int i = 0; i < coefficient.length; i++) {
            for (int j = 0; j < coefficient[0].length; j++) {
                coefficient[i][j] = new Complex(1, 0);
                // coefficient[i][j] = new Complex(i, j);
            }
        }

        transform2D(width, height, n, m, z, k, yMin, xMin, stepY, stepX, gauss, -uRange, -vRange, function, stepU, stepV);
        transform2DSuperposition(width, height, N, z, k, yMin, xMin, stepY, stepX, gauss, -uRange, -vRange, function, stepU, stepV, coefficient);
        superposition(width, height, N, gauss, coefficient, -xRange, -yRange, xStep, yStep);
        transform2DPhase(width, height, n, m, z, k, yMin, xMin, stepY, stepX, gauss, -uRange, -vRange, function, stepU, stepV);
        transform2DPhaseOnly(width, height, N, uRange, vRange, xRange, yRange, z, k, yMin, xMin, stepY, stepX, gauss, function, stepU, stepV, xStep, yStep, coefficient);
    }

    private static void transform2DPhaseOnly(int width, int height, int N, double uRange, double vRange, double xRange, double yRange, double z, double k, double yMin, double xMin, double stepY, double stepX, double gauss, Complex[][] function, double stepU, double stepV, double xStep, double yStep, Complex[][] coefficient) {
        Complex[][] superposition = new Complex[width][height];

        for (int i = 0; i < superposition.length; i++) {
            for (int j = 0; j < superposition[0].length; j++) {
                superposition[i][j] = superposition2D(N, -xRange + i * xStep, -yRange + j * yStep, gauss, coefficient);
            }
        }

        Complex[][] result = Graph.phaseOnlyEncode(superposition);

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                function[i][j] = FresnelTransform.
                        transform2DPhaseOnly(-uRange + i * stepU, -vRange + j * stepV, z, k,
                                yMin, -yMin,
                                xMin, -xMin,
                                stepX, stepY,
                                N, gauss, width, height, result, coefficient);
            }
        }

        Graph.draw2DIntensity(function, "pictures/intensity9PhaseOnly.bmp");
        Graph.draw2DPhase(function, "pictures/phase9PhaseOnly.bmp");
        Graph.draw2DPhaseSecondVersion(function, "pictures/phase9PhaseOnlyAnother.bmp");
    }

    private static void transform2D(int width, int height, int n, int m, double z, double k, double yMin, double xMin, double stepY, double stepX, double gauss, double u, double v, Complex[][] function, double stepU, double stepV) {

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                function[i][j] = FresnelTransform.transform2D(u + i * stepU , v + j * stepV, z, k,
                        yMin, xMin,
                        stepX, stepY,
                        n, m, gauss, width, height);
            }
        }

        Graph.draw2DIntensity(function, "pictures/intensity5.bmp");
        Graph.draw2DPhase(function, "pictures/phase5.bmp");
        Graph.draw2DPhaseSecondVersion(function, "pictures/phase5Another.bmp");
    }

    private static void transform2DPhase(int width, int height, int n, int m, double z, double k, double yMin, double xMin, double stepY, double stepX, double gauss, double u, double v, Complex[][] function, double stepU, double stepV) {

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                function[i][j] = FresnelTransform.transform2DPhase(u + i * stepU , v + j * stepV, z, k,
                        yMin, xMin,
                        stepX, stepY,
                        n, m, gauss, width, height);
            }
        }

        Graph.draw2DIntensity(function, "pictures/intensity5Phase.bmp");
        Graph.draw2DPhase(function, "pictures/phase5Phase.bmp");
        Graph.draw2DPhaseSecondVersion(function, "pictures/phase5PhaseAnother.bmp");

    }

    private static void superposition(int width, int height, int N, double gauss, Complex[][] coefficient, double x, double y, double xStep, double yStep) {
        Complex[][] superposition = new Complex[width][height];

        for (int i = 0; i < superposition.length; i++) {
            for (int j = 0; j < superposition[0].length; j++) {
                superposition[i][j] = superposition2D(N, x + i * xStep, y + j * yStep, gauss, coefficient);
            }
        }

        Graph.draw2DIntensity(superposition, "pictures/intensity7.bmp");
        Graph.draw2DPhase(superposition, "pictures/phase7.bmp");
        Graph.draw2DPhaseSecondVersion(superposition, "pictures/phase7Another.bmp");
    }

    private static void transform2DSuperposition(int width, int height, int N, double z, double k, double yMin, double xMin, double stepY, double stepX, double gauss, double u, double v, Complex[][] function, double stepU, double stepV, Complex[][] coefficient) {

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                function[i][j] = FresnelTransform.transform2DSuperposition(u + i * stepU, v + j * stepV, z, k,
                                yMin, xMin, stepX, stepY,
                                N, gauss, coefficient, width, height);
            }
        }

        Graph.draw2DIntensity(function, "pictures/intensity8.bmp");
        Graph.draw2DPhase(function, "pictures/phase8.bmp");
        Graph.draw2DPhaseSecondVersion(function, "pictures/phase8Another.bmp");
    }
}