import controller.FourierTransform2D;
import controller.FresnelTransform;
import model.HermiteGaussianModes;
import org.apache.commons.math3.complex.Complex;
import view.Graph;

import static controller.Superposition.*;

public class Main {

    public static void main(String[] args) {

        int width = 128, height = 128; // Размеры изображения.
        int dimension = 128;
        int N = 3; // Порядок полиномов Эрмита.
        int n = 3, m = 3; // Порядки моды для формулы 5.
        // u, v ∈ [-A; A]
        double uRange, vRange;
        uRange = vRange = 1;
        double uvRange = 1;
        // x, y ∈ [-B; B]
        double xRange, yRange;
        xRange = yRange = 1;
        double xyRange = 1;
        double yMin, xMin; //Эти значения - границы для интеграла.
        yMin = xMin = -1;
        double z = 100; // Расстояние.
        double lambda = 0.00063; // Длина волны.
        double k = 2 * Math.PI / lambda; // Волновое число.
        double stepY = -2 * yMin / height;
        double stepX = -2 * xMin / width;
        double gauss = 0.3; // Гауссовый параметр.
        Complex[][] function = new Complex[width][height];
        double stepU = 2 * uRange / width;
        double stepV = 2 * vRange / height;
        double xStep = 2 * xRange / width;
        double yStep = 2 * yRange / height;
        Complex[][] coefficient = new Complex[N + 1][N + 1]; //Набор комплексных коэффициентов.

        // Задаем им значения.
        for (int i = 0; i < coefficient.length; i++) {
            for (int j = 0; j < coefficient[0].length; j++) {
                coefficient[i][j] = new Complex(1, 1);
            }
        }

        /*uRange = vRange = 0.2;
        yMin = xMin = -0.2;
        stepY = -2 * yMin / height;
        stepX = -2 * xMin / width;
        gauss = 0.015;
        stepU = 2 * uRange / width;
        stepV = 2 * vRange / height;*/


        mode2D(dimension, n, m, xyRange, gauss);

        transform2DByFresnel(dimension, uvRange, xyRange, n, m, gauss, z, k);
        // transform2DPhase(width, height, n, m, z, k, yMin, xMin, stepY, stepX, gauss, -uRange, -vRange, function, stepU, stepV);
        // superposition(width, height, N, gauss, coefficient, -xRange, -yRange, xStep, yStep);
        // transform2DSuperpositionFresnel(width, height, N, z, k, yMin, xMin, stepY, stepX, gauss, -uRange, -vRange, function, stepU, stepV, coefficient);

        // transform2DPhaseOnlyFresnel(width, height, N, uRange, vRange, xRange, yRange, z, k, yMin, xMin, stepY, stepX, gauss, function, stepU, stepV, xStep, yStep, coefficient);

        /*uRange = vRange = 2.5;
        yMin = xMin = -1;
        stepY = -2 * yMin / height;
        stepX = -2 * xMin / width;
        gauss = 0.015;
        stepU = 2 * uRange / width;
        stepV = 2 * vRange / height;*/

        // transform2DSuperpositionFourier(width, height, N, z, k, yMin, xMin, stepY, stepX, gauss, -uRange, -vRange, function, stepU, stepV, coefficient);
        // transform2DPhaseOnlyFourier(width, height, uRange, vRange, z, k, yMin, xMin, stepY, stepX, gauss, function, stepU, stepV);


        // modesSum(width, height, yMin, xMin, stepY, stepX, gauss, function);
    }

    /**
     * Рисует интенсивность и фазу моды Гаусса-Эрмита.
     *
     * @param dimension Размерность выходного изображения dimension x dimension.
     * @param n         Порядок моды Гаусса-Эрмита.
     * @param m         Порядок моды Гаусса-Эрмита.
     * @param xyRange   Диапазон значений по осям x и y, по которым рисуются изображения.
     * @param gauss     Гауссовый параметр.
     */
    private static void mode2D(int dimension, int n, int m, double xyRange, double gauss) {
        Complex[][] result = new Complex[dimension][dimension];
        double step = 2 * xyRange / dimension;

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                result[i][j] = new Complex(HermiteGaussianModes.hermiteGauss2D(n, m, -xyRange + i * step, -xyRange + j * step, gauss));
            }
        }

        Graph.draw2DIntensity(result, "pictures/intensityGH" + n + m + " xy" + xyRange + " gauss " + gauss + ".bmp");
        Graph.draw2DPhase(result, "pictures/phaseGH" + n + m + " xy" + xyRange + " gauss " + gauss + ".bmp");
    }

    /**
     * Рисует интенсивность и фазу преобразования Френеля от моды 2D Гаусса-Эрмита.
     *
     * @param dimension Размерность выходного изображения dimension x dimension.
     * @param uvRange   Диапазон значений по осям u и v, по которым рисуются изображения.
     * @param xyRange   Диапазон значений по осям x и y, по которым рисуются изображения.
     * @param n         Порядок моды Гаусса-Эрмита.
     * @param m         Порядок моды Гаусса-Эрмита.
     * @param gauss     Гауссовый параметр.
     * @param z         Расстояние.
     * @param k         Волновое число.
     */
    private static void transform2DByFresnel(int dimension, double uvRange, double xyRange, int n, int m, double gauss, double z, double k) {

        Complex[][] result = new Complex[dimension][dimension];
        double stepUV = 2 * uvRange / dimension;
        double stepXY = 2 * xyRange / dimension;

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                result[i][j] = FresnelTransform.transform2D(-uvRange + i * stepUV, -uvRange + j * stepUV, z, k,
                        -xyRange, -xyRange,
                        stepXY, stepXY,
                        n, m, gauss, dimension, dimension);
            }
        }

        Graph.draw2DIntensity(result, "pictures/intensityFresnelGH" + n + m + " xy" + xyRange + " gauss " + gauss + " uv" + uvRange + " z" + z + " k" + k + ".bmp");
        Graph.draw2DPhase(result, "pictures/phaseFresnelGH" + n + m + " xy" + xyRange + " gauss " + gauss + " uv" + uvRange + " z" + z + " k" + k + ".bmp");
    }

    private static void modesSum(int width, int height, double yMin, double xMin, double stepY, double stepX, double gauss, Complex[][] function) {
        Complex[][] one = new Complex[width][height];
        Complex[][] two = new Complex[width][height];
        Complex[][] three = new Complex[width][height];

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                one[i][j] = Complex.valueOf(HermiteGaussianModes.hermiteGauss2D(6, 0, xMin + i * stepX, yMin + j * stepY, gauss)).multiply(new Complex(1, 5));
                two[i][j] = Complex.valueOf(HermiteGaussianModes.hermiteGauss2D(0, 6, xMin + i * stepX, yMin + j * stepY, gauss)).multiply(new Complex(5, 1));
                three[i][j] = Complex.valueOf(HermiteGaussianModes.hermiteGauss2D(3, 3, xMin + i * stepX, yMin + j * stepY, gauss)).multiply(new Complex(10, 10));
                function[i][j] = new Complex(one[i][j].getReal() + two[i][j].getReal() + three[i][j].getReal(), one[i][j].getImaginary() + two[i][j].getImaginary() + three[i][j].getImaginary());
            }
        }

        Graph.draw2DIntensity(function, "pictures/intensity7.bmp");
        Graph.draw2DPhase(function, "pictures/phase7.bmp");
    }

    private static void transform2DPhaseOnlyFresnel(int width, int height, int N, double uRange, double vRange, double xRange, double yRange, double z, double k, double yMin, double xMin, double stepY, double stepX, double gauss, Complex[][] function, double stepU, double stepV, double xStep, double yStep, Complex[][] coefficient) {
        Complex[][] superposition = new Complex[width][height];
        Complex[][] one = new Complex[width][height];
        Complex[][] two = new Complex[width][height];
        Complex[][] three = new Complex[width][height];

        for (int i = 0; i < superposition.length; i++) {
            for (int j = 0; j < superposition[0].length; j++) {
                one[i][j] = Complex.valueOf(HermiteGaussianModes.hermiteGauss2D(6, 0, xMin + i * stepX, yMin + j * stepY, gauss)).multiply(new Complex(1, 5));
                two[i][j] = Complex.valueOf(HermiteGaussianModes.hermiteGauss2D(0, 6, xMin + i * stepX, yMin + j * stepY, gauss)).multiply(new Complex(5, 1));
                three[i][j] = Complex.valueOf(HermiteGaussianModes.hermiteGauss2D(3, 3, xMin + i * stepX, yMin + j * stepY, gauss)).multiply(new Complex(10, 10));
                superposition[i][j] = new Complex(one[i][j].getReal() + two[i][j].getReal() + three[i][j].getReal(), one[i][j].getImaginary() + two[i][j].getImaginary() + three[i][j].getImaginary());
                // superposition[i][j] = superposition2D(N, -xRange + i * xStep, -yRange + j * yStep, gauss, coefficient);
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

        Graph.draw2DIntensity(function, "pictures/intensitySuperpositionPhaseOnlyOutputByFresnel.bmp");
        Graph.draw2DPhase(function, "pictures/phaseSuperpositionPhaseOnlyOutputByFresnel.bmp");
    }

    private static void transform2DPhaseOnlyFourier(int width, int height, double uRange, double vRange, double z, double k, double yMin, double xMin, double stepY, double stepX, double gauss, Complex[][] function, double stepU, double stepV) {
        Complex[][] superposition = new Complex[width][height];
        Complex[][] one = new Complex[width][height];
        Complex[][] two = new Complex[width][height];
        Complex[][] three = new Complex[width][height];

        for (int i = 0; i < superposition.length; i++) {
            for (int j = 0; j < superposition[0].length; j++) {
                one[i][j] = Complex.valueOf(HermiteGaussianModes.hermiteGauss2D(6, 0, xMin + i * stepX, yMin + j * stepY, gauss)).multiply(new Complex(1, 5));
                two[i][j] = Complex.valueOf(HermiteGaussianModes.hermiteGauss2D(0, 6, xMin + i * stepX, yMin + j * stepY, gauss)).multiply(new Complex(5, 1));
                three[i][j] = Complex.valueOf(HermiteGaussianModes.hermiteGauss2D(3, 3, xMin + i * stepX, yMin + j * stepY, gauss)).multiply(new Complex(10, 10));
                superposition[i][j] = new Complex(one[i][j].getReal() + two[i][j].getReal() + three[i][j].getReal(), one[i][j].getImaginary() + two[i][j].getImaginary() + three[i][j].getImaginary());
            }
        }

        Complex[][] phaseOnly = Graph.phaseOnlyEncode(superposition);

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                function[i][j] = FourierTransform2D.
                        transform2DPhaseOnly(-uRange + i * stepU, -vRange + j * stepV, z, k,
                                yMin, xMin,
                                stepX, stepY, width, height, phaseOnly);
            }
        }

        Graph.draw2DIntensity(function, "pictures/intensitySuperpositionPhaseOnlyOutputByFourier.bmp");
        Graph.draw2DPhase(function, "pictures/phaseSuperpositionPhaseOnlyOutputByFourier.bmp");
    }

    private static void transform2DPhase(int width, int height, int n, int m, double z, double k, double yMin, double xMin, double stepY, double stepX, double gauss, double u, double v, Complex[][] function, double stepU, double stepV) {

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                function[i][j] = FresnelTransform.transform2DPhase(u + i * stepU, v + j * stepV, z, k,
                        yMin, xMin,
                        stepX, stepY,
                        n, m, gauss,
                        width, height);
            }
        }

        Graph.draw2DIntensity(function, "pictures/intensityModeArgumentOutput.bmp");
        Graph.draw2DPhase(function, "pictures/phaseModeArgumentOutput.bmp");
    }

    private static void superposition(int width, int height, int N, double gauss, Complex[][] coefficient, double x, double y, double xStep, double yStep) {
        Complex[][] superposition = new Complex[width][height];

        for (int i = 0; i < superposition.length; i++) {
            for (int j = 0; j < superposition[0].length; j++) {
                // superposition[i][j] = new Complex(HermiteGaussianModes.hermiteGauss2D(3, 3, x, y, gauss));
                superposition[i][j] = superposition2D(N, x + i * xStep, y + j * yStep, gauss, coefficient);
            }
        }

        Graph.draw2DIntensity(superposition, "pictures/intensitySuperpositionInput.bmp");
        Graph.draw2DPhase(superposition, "pictures/phaseSuperpositionInput.bmp");
    }

    private static void transform2DSuperpositionFresnel(int width, int height, int N, double z, double k, double yMin, double xMin, double stepY, double stepX, double gauss, double u, double v, Complex[][] function, double stepU, double stepV, Complex[][] coefficient) {

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                function[i][j] = FresnelTransform.transform2DSuperposition(u + i * stepU, v + j * stepV, z, k,
                        yMin, xMin, stepX, stepY,
                        N, gauss, coefficient, width, height);
            }
        }

        Graph.draw2DIntensity(function, "pictures/intensitySuperpositionOutputFresnel.bmp");
        Graph.draw2DPhase(function, "pictures/phaseSuperpositionOutputFresnel.bmp");
    }

    private static void transform2DSuperpositionFourier(int width, int height, int N, double z, double k, double yMin, double xMin, double stepY, double stepX, double gauss, double u, double v, Complex[][] function, double stepU, double stepV, Complex[][] coefficient) {

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                function[i][j] = FourierTransform2D.transform2DSuperposition(u + i * stepU, v + j * stepV, z, k,
                        yMin, xMin, stepX, stepY,
                        gauss, width, height);
            }
        }

        Graph.draw2DIntensity(function, "pictures/intensitySuperpositionOutput.bmp");
        Graph.draw2DPhase(function, "pictures/phaseSuperpositionOutput.bmp");
    }
}