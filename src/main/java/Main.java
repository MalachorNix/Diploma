import controller.FourierTransform;
import controller.FresnelTransform;
import org.apache.commons.math3.complex.Complex;
import view.Graph;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import static model.GHMode.hermiteGauss2D;

public class Main {

    public static void main(String[] args) throws IOException {

        int dimension = 128; // Размеры изображения.
        int n = 3, m = 3; // Порядки моды.
        // u, v ∈ [-A; A]
        double uvRange = 1;
        // x, y ∈ [-B; B]
        double xyRange = 1;
        double z = 100; // Расстояние.
        double lambda = 0.00063; // Длина волны.
        double k = 2 * Math.PI / lambda; // Волновое число.
        double gauss = 0.3; // Гауссовый параметр.
        double threshold = 0.33;
        double sigma0 = 0.6;

        // mode2D(dimension, xyRange, n, m, gauss);
        // superposition(dimension, xyRange, gauss);
        // fresnelTransform2D(dimension, uvRange, xyRange, n, m, gauss, z, k);
        // fresnelTransform2DPhaseOnly(dimension, uvRange, xyRange, n, m, gauss, z, k, sigma0);
        // fourierTransform2D(dimension, uvRange, xyRange, n, m, gauss, z, k);
        // fourierTransform2DPhaseOnly(dimension, uvRange, xyRange, n, m, gauss, z, k, sigma0);
        // fresnelTransform2DSuperposition(dimension, uvRange, xyRange, gauss, z, k);
        // fresnelTransform2DSuperpositionPhaseOnly(dimension, uvRange, xyRange, gauss, z, k);
        // transform2DSuperpositionFourier(dimension, uvRange, xyRange, gauss, z, k);
        // transform2DPhaseOnlyFourier(dimension, uvRange, xyRange, gauss, z, k);
        // expAndLesem(dimension, xyRange, n, m, gauss, sigma0);
        // partialCoding(dimension, xyRange, n, m, gauss, threshold, sigma0);
        // pcSuper(dimension, xyRange, gauss, threshold, sigma0);


        File real = new File("real.bmp");
        File except = new File("except.bmp");
        rms(real, except);
    }

    public static void rms(File realFile,  File exceptFile) throws IOException {
        BufferedImage bufferReal = ImageIO.read(realFile);
        BufferedImage bufferExcept = ImageIO.read(exceptFile);
        double[][] real = new double[bufferReal.getHeight()][bufferReal.getWidth()];
        double[][] except = new double[bufferExcept.getHeight()][bufferExcept.getWidth()];
        double deviation;
        double sum1 = 0, sum2 = 0;
        for (int i = 0; i < real.length; i++) {
            for (int j = 0; j < real[0].length; j++) {
                real[i][j] = bufferReal.getRGB(i, j);
                except[i][j] = bufferExcept.getRGB(i, j);
            }
        }

        for (int i = 0; i < real.length; i++) {
            for (int j = 0; j < real[0].length; j++) {
                double one = Math.abs(except[i][j]) * Math.abs(except[i][j]) - Math.abs(real[i][j]) * Math.abs(real[i][j]);
                sum1 += one * one;
                double two = Math.abs(except[i][j]) * Math.abs(except[i][j]) * Math.abs(except[i][j]) * Math.abs(except[i][j]);
                sum2 += two;
            }
        }

        double reverse = 1 / sum2;
        deviation = Math.sqrt(sum1 * reverse) * 100;
        System.out.println(deviation);

    }

    private static double standardDeviation(Complex[][] reality, Complex[][] expect) {
        double deviation;

        int length = reality.length;
        double absReality, absExpect, firstMultiplier, secondMultiplier, firstSum = 0, secondSum = 0;

        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                absReality = reality[i][j].abs();
                absExpect = expect[i][j].abs();
                firstMultiplier = absReality * absReality - absExpect * absExpect;
                firstMultiplier *= firstMultiplier;
                firstSum += firstMultiplier;

                secondMultiplier = Math.pow(absExpect, 4);
                secondSum += secondMultiplier;
            }
        }

        double one = firstSum;
        double two = Math.pow(secondSum, -1);
        deviation = Math.sqrt(one * two) * 100;

        return deviation;
    }

    private static Complex[][] expAndLesem(int dimension, double xyRange, int n, int m, double gauss, double sigma0) {
        Complex[][] exponent = new Complex[dimension][dimension];
        double step = 2 * xyRange / dimension;
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                double x = -xyRange + i * step;
                double y = -xyRange + j * step;
                double exp = Math.exp(-(x * x + y * y) / (2 * sigma0 * sigma0));
                Complex mode = new Complex(hermiteGauss2D(n, m, x, y, gauss));
                Complex phaseOnly = Complex.I.multiply(mode.getArgument()).exp();
                exponent[i][j] = Complex.valueOf(exp).multiply(phaseOnly);
            }
        }
        Graph.draw2DIntensity(exponent, "pictures/intensityEXP sigma" + sigma0 + ".bmp");
        Graph.draw2DPhase(exponent, "pictures/phaseEXP sigma" + sigma0 + ".bmp");
        return exponent;
    }

    private static void partialCoding(int dimension, double xyRange, int n, int m, double gauss, double threshold, double sigma0) {
        double step = 2 * xyRange / dimension;
        double[][] amplitude = new double[dimension][dimension];
        Complex[][] picture = new Complex[2 * dimension][2 * dimension];
        Complex[][] mode = new Complex[dimension][dimension];
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                mode[i][j] = Complex.valueOf(hermiteGauss2D(n, m, -xyRange + i * step, -xyRange + j * step, gauss));
                amplitude[i][j] = mode[i][j].abs();
                if (amplitude[i][j] > max) {
                    max = amplitude[i][j];
                }
                if (amplitude[i][j] < min) {
                    min = amplitude[i][j];
                }
            }
        }
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                double normAmplitude = normFrom0to1(amplitude[i][j], min, max);
                double x = -xyRange + i * step;
                double y = -xyRange + j * step;
                double exp = Math.exp(-(x * x + y * y) / (2 * sigma0 * sigma0));
                if (normAmplitude >= threshold) {
                    picture[2 * i][2 * j] = Complex.I.multiply(mode[i][j].getArgument()).exp();
                    picture[2 * i][2 * j + 1] = Complex.I.multiply(mode[i][j].getArgument()).exp();
                    picture[2 * i + 1][2 * j] = Complex.I.multiply(mode[i][j].getArgument()).exp();
                    picture[2 * i + 1][2 * j + 1] = Complex.I.multiply(mode[i][j].getArgument()).exp();

                    picture[2 * i][2 * j] = picture[2 * i][2 * j].multiply(exp);
                    picture[2 * i][2 * j + 1] = picture[2 * i][2 * j + 1].multiply(exp);
                    picture[2 * i + 1][2 * j] = picture[2 * i + 1][2 * j].multiply(exp);
                    picture[2 * i + 1][2 * j + 1] = picture[2 * i + 1][2 * j + 1].multiply(exp);
                } else {
                    Complex phi1, phi2;
                    phi1 = new Complex(mode[i][j].getArgument() + Math.acos(amplitude[i][j]));
                    phi2 = new Complex(mode[i][j].getArgument() - Math.acos(amplitude[i][j]));
                    Complex first = Complex.I.multiply(phi1).exp();
                    Complex second = Complex.I.multiply(phi2).exp();
                    picture[2 * i][2 * j] = first;
                    picture[2 * i][2 * j + 1] = second;
                    picture[2 * i + 1][2 * j] = first;
                    picture[2 * i + 1][2 * j + 1] = second;

                    picture[2 * i][2 * j] = picture[2 * i][2 * j].multiply(exp);
                    picture[2 * i][2 * j + 1] = picture[2 * i][2 * j + 1].multiply(exp);
                    picture[2 * i + 1][2 * j] = picture[2 * i + 1][2 * j].multiply(exp);
                    picture[2 * i + 1][2 * j + 1] = picture[2 * i + 1][2 * j + 1].multiply(exp);
                }
            }
        }
        // Complex[][] exponent = expAndLesem(dimension, xyRange, n, m, gauss, sigma0);
        int[][] reality = Graph.draw2DIntensity(picture, "pictures/ВЫХОД partInt" + "threshold" + threshold + " sigma0 " + sigma0 + ".bmp");
        // int[][] exp = Graph.draw2DIntensity(exponent, "pictures/gk;nrtu.bmp");
        Graph.draw2DPhase(picture, "pictures/ВЫХОД partPhase" + "threshold" + threshold + " sigma0 " + sigma0 + ".bmp");
        // System.out.println(Graph.standardDeviation(reality, exp) + " значения на эксп");
    }

    private static void pcSuper(int dimension, double xyRange, double gauss, double threshold, double sigma0) {
        Complex[][] superposition = new Complex[dimension][dimension];
        Complex[][] picture = new Complex[dimension][dimension];
        Complex[][] superExp = new Complex[dimension][dimension];
        double[][] amplitude = new double[dimension][dimension];
        double[][] normAmplitude = new double[dimension][dimension];
        double stepXY = 2 * xyRange / dimension;
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                double x = -xyRange + i * stepXY;
                double y = -xyRange + j * stepXY;
                Complex one = Complex.valueOf(hermiteGauss2D(6, 0, x, y, gauss)).multiply(Complex.valueOf(1));
                Complex two = Complex.valueOf(hermiteGauss2D(0, 6, x, y, gauss)).multiply(Complex.valueOf(5));
                Complex three = Complex.valueOf(hermiteGauss2D(3, 3, x, y, gauss)).multiply(Complex.valueOf(10));
                superposition[i][j] = Complex.valueOf(one.getReal() + two.getReal() + three.getReal(), one.getImaginary() + two.getImaginary() + three.getImaginary());
                amplitude[i][j] = superposition[i][j].abs();
                if (amplitude[i][j] > max) {
                    max = amplitude[i][j];
                }
                if (amplitude[i][j] < min) {
                    min = amplitude[i][j];
                }
            }
        }
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                double x = -xyRange + i * stepXY;
                double y = -xyRange + j * stepXY;
                double exp = Math.exp(-(x * x + y * y) / (2 * sigma0 * sigma0));
                normAmplitude[i][j] = normFrom0to1(amplitude[i][j], min, max);
                if (normAmplitude[i][j] >= threshold) {
                    picture[i][j] = Complex.I.multiply(superposition[i][j].getArgument()).exp();
                    superExp[i][j] = picture[i][j].multiply(exp);
                } else {
                    picture[i][j] = anotherPartCoding(superposition[i][j], amplitude[i][j]);
                    superExp[i][j] = picture[i][j].multiply(exp);
                }
            }
        }
        Graph.draw2DIntensity(picture, "pictures/IntSuper" + "threshold" + threshold + ".bmp");
        Graph.draw2DPhase(picture, "pictures/PhaseSuper" + "threshold" + threshold + ".bmp");
        Graph.draw2DIntensity(superExp, "pictures/IntSuperExp" + "threshold" + threshold + ".bmp");
        Graph.draw2DPhase(superExp, "pictures/PhaseSuperExp" + "threshold" + threshold + ".bmp");
    }

    private static Complex anotherPartCoding(Complex function, double amplitude) {
        double phi1, phi2;
        Complex picture;
        phi1 = function.getArgument() + Math.acos(amplitude);
        phi2 = function.getArgument() - Math.acos(amplitude);
        Complex first = Complex.I.multiply(phi1).exp();
        Complex second = Complex.I.multiply(phi2).exp();
        picture = Complex.valueOf(first.getReal() + second.getReal(), first.getImaginary() + second.getImaginary());
        return picture;
    }

    private static double normFrom0to1(double x, double xMin, double xMax) {
        return (x - xMin) / (xMax - xMin);
    }

    /**
     * Рисует интенсивность и фазу 2D моды Гаусса-Эрмита.
     *
     * @param dimension Размерность выходного изображения dimension x dimension.
     * @param n         Порядок моды Гаусса-Эрмита.
     * @param m         Порядок моды Гаусса-Эрмита.
     * @param xyRange   Диапазон значений по осям x и y, по которым рисуются изображения.
     * @param gauss     Гауссовый параметр.
     */
    private static void mode2D(int dimension, double xyRange, int n, int m, double gauss) {
        Complex[][] mode = new Complex[dimension][dimension];
        Complex[][] phaseOnly = new Complex[dimension][dimension];
        double step = 2 * xyRange / dimension;
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                double x = -xyRange + i * step;
                double y = -xyRange + j * step;
                mode[i][j] = new Complex(hermiteGauss2D(n, m, x, y, gauss));
                phaseOnly[i][j] = Complex.I.multiply(mode[i][j].getArgument()).exp();
            }
        }
        int[][] reality = Graph.draw2DIntensity(phaseOnly, "pictures/intensityPhaseGH" + n + m + " xy" + xyRange + " sigma " + gauss + ".bmp");
        int[][] expect = Graph.draw2DIntensity(mode, "pictures/intensityGH" + n + m + " xy" + xyRange + " sigma " + gauss + ".bmp");
        System.out.println(Graph.standardDeviation(expect, reality) + " пиксели интенсивности");
        reality = Graph.draw2DPhase(phaseOnly, "pictures/phasePhaseGH" + n + m + " xy" + xyRange + " gauss " + gauss + ".bmp");
        expect = Graph.draw2DPhase(mode, "pictures/phaseGH" + n + m + " xy" + xyRange + " gauss " + gauss + ".bmp");
        System.out.println(Graph.standardDeviation(expect, reality) + " пиксели фазы");
        System.out.println(standardDeviation(phaseOnly, mode) + " значения");

    }

    /**
     * Рисует интенсивность и фазу преобразования Френеля от 2D моды Гаусса-Эрмита.
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
    private static int[][] fresnelTransform2D(int dimension, double uvRange, double xyRange, int n, int m, double gauss, double z, double k) {
        Complex[][] result = new Complex[dimension][dimension];
        Complex[][] mode = new Complex[dimension][dimension];
        double stepUV = 2 * uvRange / dimension;
        double stepXY = 2 * xyRange / dimension;
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                mode[i][j] = new Complex(hermiteGauss2D(n, m, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss));
            }
        }
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                result[i][j] = FresnelTransform.transform2D(-uvRange + i * stepUV, -uvRange + j * stepUV, z, k,
                        xyRange, stepXY, n, m, gauss, dimension, mode);
            }
        }
        // mode2D(dimension, xyRange, n, m, gauss);

        int[][] intens = Graph.draw2DIntensity(result, "pictures/intensityFresnelGH" + n + m + " xy" + xyRange + " gauss " + gauss + " uv" + uvRange + " z" + z + ".bmp");
        Graph.draw2DPhase(result, "pictures/phaseFresnelGH" + n + m + " xy" + xyRange + " gauss " + gauss + " uv" + uvRange + " z" + z + ".bmp");
        return intens;
    }

    /**
     * Рисует интенсивность и фазу чисто фазового поля преобразования Френеля от 2D моды Гаусса-Эрмита.
     * (от exp^(i * arg(mode)))
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
    private static int[][] fresnelTransform2DPhaseOnly(int dimension, double uvRange, double xyRange, int n, int m,
                                                       double gauss, double z, double k, double sigma0) {
        Complex[][] result = new Complex[dimension][dimension];
        Complex[][] mode = new Complex[dimension][dimension];
        double stepUV = 2 * uvRange / dimension;
        double stepXY = 2 * xyRange / dimension;

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                mode[i][j] = new Complex(hermiteGauss2D(n, m, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss));
            }
        }

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                result[i][j] = FresnelTransform.transform2DPhase(-uvRange + i * stepUV, -uvRange + j * stepUV, z, k,
                        xyRange, stepXY, n, m, gauss, dimension, mode, sigma0);
            }
        }
        // mode2D(dimension, xyRange, n, m, gauss);
        int[][] intens = Graph.draw2DIntensity(result, "pictures/intensityFresnelPhaseGH" + n + m + " xy" + xyRange + " gauss " + gauss + " uv" + uvRange + " z" + z + " k" + k + ".bmp");
        Graph.draw2DPhase(result, "pictures/phaseFresnelPhaseGH" + n + m + " xy" + xyRange + " gauss " + gauss + " uv" + uvRange + " z" + z + " k" + k + ".bmp");
        return intens;
    }

    /**
     * Рисует интенсивность и фазу преобразования Фурье от 2D моды Гаусса-Эрмита.
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
    private static void fourierTransform2D(int dimension, double uvRange, double xyRange, int n, int m, double gauss, double z, double k) {
        Complex[][] result = new Complex[dimension][dimension];
        double stepUV = 2 * uvRange / dimension;
        double stepXY = 2 * xyRange / dimension;

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                result[i][j] = FourierTransform.transform2D(-uvRange + i * stepUV, -uvRange + j * stepUV, z, k,
                        -xyRange, -xyRange,
                        stepXY, stepXY,
                        n, m, gauss, dimension, dimension);
            }
        }
        mode2D(dimension, xyRange, n, m, gauss);
        Graph.draw2DIntensity(result, "pictures/intensityFourierGH" + n + m + " xy" + xyRange + " gauss " + gauss + " uv" + uvRange + " z" + z + " k" + k + ".bmp");
        Graph.draw2DPhase(result, "pictures/phaseFourierGH" + n + m + " xy" + xyRange + " gauss " + gauss + " uv" + uvRange + " z" + z + " k" + k + ".bmp");
    }

    /**
     * Рисует интенсивность и фазу чисто фазового поля преобразования Фурье от 2D моды Гаусса-Эрмита.
     * (от exp^(i * arg(mode)))
     *
     * @param dimension Размерность выходного изображения dimension x dimension.
     * @param uvRange   Диапазон значений по осям u и v, по которым рисуются изображения.
     * @param xyRange   Размерность выходного изображения dimension x dimension.
     * @param n         Порядок моды Гаусса-Эрмита.
     * @param m         Порядок моды Гаусса-Эрмита.
     * @param gauss     Гауссовый параметр.
     * @param z         Расстояние.
     * @param k         Волновое число.
     */
    private static void fourierTransform2DPhaseOnly(int dimension, double uvRange, double xyRange, int n, int m, double gauss, double z, double k, double sigma0) {

        Complex[][] result = new Complex[dimension][dimension];
        double stepUV = 2 * uvRange / dimension;
        double stepXY = 2 * xyRange / dimension;

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                result[i][j] = FourierTransform.transform2DPhaseOnly(-uvRange + i * stepUV, -uvRange + j * stepUV, z, k,
                        -xyRange, -xyRange,
                        stepXY, stepXY,
                        n, m, gauss, dimension, dimension, sigma0);
            }
        }
        // mode2D(dimension, xyRange, n, m, gauss);
        Graph.draw2DIntensity(result, "pictures/intensityPhaseFourierGH" + n + m + " xy" + xyRange + " gauss " + gauss + " uv" + uvRange + " z" + z + " k" + k + ".bmp");
        Graph.draw2DPhase(result, "pictures/phasePhaseFourierGH" + n + m + " xy" + xyRange + " gauss " + gauss + " uv" + uvRange + " z" + z + " k" + k + ".bmp");
    }

    /**
     * Рисует интенсивность и фазу суперпозицию 2D мод Гаусса-Эрмита.
     *
     * @param dimension Размерность выходного изображения dimension x dimension.
     * @param xyRange   Размерность выходного изображения dimension x dimension.
     * @param gauss     Гауссовый параметр.
     */
    private static void superposition(int dimension, double xyRange, double gauss) {
        Complex[][] result = new Complex[dimension][dimension];
        Complex[][] phaseOnly = new Complex[dimension][dimension];
        double stepXY = 2 * xyRange / dimension;
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                Complex one = Complex.valueOf(hermiteGauss2D(6, 0, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(1));
                Complex two = Complex.valueOf(hermiteGauss2D(0, 6, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(5));
                Complex three = Complex.valueOf(hermiteGauss2D(3, 3, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(10));
                result[i][j] = Complex.valueOf(one.getReal() + two.getReal() + three.getReal(), one.getImaginary() + two.getImaginary() + three.getImaginary());
                phaseOnly[i][j] = Complex.I.multiply(result[i][j].getArgument()).exp();
            }
        }
        Graph.draw2DIntensity(phaseOnly, "pictures/intensitySuperpositionPhaseOnly.bmp");
        Graph.draw2DPhase(phaseOnly, "pictures/phaseSuperpositionPhaseOnly.bmp");
        Graph.draw2DIntensity(result, "pictures/intensitySuperposition.bmp");
        Graph.draw2DPhase(result, "pictures/phaseSuperposition.bmp");
    }

    private static int[][] fresnelTransform2DSuperposition(int dimension, double uvRange, double xyRange, double gauss,
                                                           double z, double k) {
        double stepXY = 2 * xyRange / dimension;
        double stepUV = 2 * uvRange / dimension;
        Complex[][] result = new Complex[dimension][dimension];
        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[0].length; j++) {
                result[i][j] = FresnelTransform.transform2DSuperposition(-uvRange + i * stepUV, -uvRange + j * stepUV,
                        z, k, xyRange, stepXY, gauss, dimension);
            }
        }
        superposition(dimension, xyRange, gauss);
        int[][] intens = Graph.draw2DIntensity(result, "pictures/intensitySuperpositionOutputFresnel.bmp");
        Graph.draw2DPhase(result, "pictures/phaseSuperpositionOutputFresnel.bmp");
        return intens;
    }

    private static int[][] fresnelTransform2DSuperpositionPhaseOnly(int dimension, double uvRange, double xyRange,
                                                                    double gauss, double z, double k) {
        Complex[][] superposition;
        Complex[][] function = new Complex[dimension][dimension];
        double stepXY = 2 * xyRange / dimension;
        double stepUV = 2 * uvRange / dimension;
        superposition = superpos(dimension, xyRange, gauss, stepXY);
        Complex[][] phaseOnly = Graph.phaseOnlyEncode(superposition);
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                function[i][j] = FresnelTransform.transform2DPhaseOnly(-uvRange + i * stepUV, -uvRange + j * stepUV,
                        z, k, xyRange, stepXY, dimension, phaseOnly);
            }
        }
        // superposition(dimension, xyRange, gauss);
        int[][] intens = Graph.draw2DIntensity(function, "pictures/intensitySuperpositionPhaseOnlyOutputByFresnel.bmp");
        Graph.draw2DPhase(function, "pictures/phaseSuperpositionPhaseOnlyOutputByFresnel.bmp");
        return intens;
    }

    private static Complex[][] superpos(int dimension, double xyRange, double gauss, double stepXY) {
        Complex[][] superposition = new Complex[dimension][dimension];
        Complex[][] one = new Complex[dimension][dimension];
        Complex[][] two = new Complex[dimension][dimension];
        Complex[][] three = new Complex[dimension][dimension];
        for (int i = 0; i < superposition.length; i++) {
            for (int j = 0; j < superposition[0].length; j++) {
                one[i][j] = Complex.valueOf(hermiteGauss2D(6, 0, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(1));
                two[i][j] = Complex.valueOf(hermiteGauss2D(0, 6, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(5));
                three[i][j] = Complex.valueOf(hermiteGauss2D(3, 3, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(10));
                superposition[i][j] = Complex.valueOf(one[i][j].getReal() + two[i][j].getReal() + three[i][j].getReal(), one[i][j].getImaginary() + two[i][j].getImaginary() + three[i][j].getImaginary());
            }
        }
        return superposition;
    }

    private static void transform2DSuperpositionFourier(int dimension, double uvRange, double xyRange, double gauss, double z, double k) {

        Complex[][] function = new Complex[dimension][dimension];
        double stepUV = 2 * uvRange / dimension;
        double stepXY = 2 * xyRange / dimension;

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                function[i][j] = FourierTransform.transform2DSuperposition(-uvRange + i * stepUV, -uvRange + j * stepUV, z, k,
                        -xyRange, -xyRange, stepXY, stepXY,
                        gauss, dimension, dimension);
            }
        }
        superposition(dimension, xyRange, gauss);
        Graph.draw2DIntensity(function, "pictures/intensitySuperpositionOutput.bmp");
        Graph.draw2DPhase(function, "pictures/phaseSuperpositionOutput.bmp");
    }

    private static void transform2DPhaseOnlyFourier(int dimension, double uvRange, double xyRange, double gauss, double z, double k) {
        Complex[][] superposition;
        Complex[][] function = new Complex[dimension][dimension];
        double stepUV = 2 * uvRange / dimension;
        double stepXY = 2 * xyRange / dimension;

        superposition = superpos(dimension, xyRange, gauss, stepXY);

        Complex[][] phaseOnly = Graph.phaseOnlyEncode(superposition);

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                function[i][j] = FourierTransform.
                        transform2DPhaseOnly(-uvRange + i * stepUV, -uvRange + j * stepUV, z, k,
                                -xyRange, -xyRange,
                                stepXY, stepXY, dimension, dimension, phaseOnly);
            }
        }
        superposition(dimension, xyRange, gauss);
        Graph.draw2DIntensity(function, "pictures/intensitySuperpositionPhaseOnlyOutputByFourier.bmp");
        Graph.draw2DPhase(function, "pictures/phaseSuperpositionPhaseOnlyOutputByFourier.bmp");
    }
}