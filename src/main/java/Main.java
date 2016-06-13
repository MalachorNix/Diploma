import controller.FourierTransform2D;
import controller.FresnelTransform;
import org.apache.commons.math3.complex.Complex;
import view.Graph;

import static model.GHMode.hermiteGauss2D;

public class Main {

    public static void main(String[] args) {

        int dimension = 128; // Размеры изображения.
        int n = 3, m = 3; // Порядки моды для формулы 5.
        // u, v ∈ [-A; A]
        double uvRange = 1;
        // x, y ∈ [-B; B]
        double xyRange = 1;
        double z = 100; // Расстояние.
        double lambda = 0.00063; // Длина волны.
        double k = 2 * Math.PI / lambda; // Волновое число.
        double gauss = 0.3; // Гауссовый параметр.
        double threshold = 1;

        fresnelTransform2D(dimension, uvRange, xyRange, n, m, gauss, z, k);
        fresnelTransform2DPhaseOnly(dimension, uvRange, xyRange, n, m, gauss, z, k);
        fourierTransform2D(dimension, uvRange, xyRange, n, m, gauss, z, k);
        fourierTransform2DPhaseOnly(dimension, uvRange, xyRange, n, m, gauss, z, k);
        fresnelTransform2DSuperposition(dimension, uvRange, xyRange, gauss, z, k);
        fresnelTransform2DSuperpositionPhaseOnly(dimension, uvRange, xyRange, gauss, z, k);
        transform2DSuperpositionFourier(dimension, uvRange, xyRange, gauss, z, k);
        transform2DPhaseOnlyFourier(dimension, uvRange, xyRange, gauss, z, k);

        partialCoding(dimension, xyRange, n, m, gauss, threshold);

        partialCodingRectangularJ(dimension, xyRange, n, m, gauss, threshold);

        partialCodingRectangularI(dimension, xyRange, n, m, gauss, threshold);

        pcSuper(dimension, xyRange, gauss, threshold);

        pcSuperRectangularI(dimension, xyRange, gauss, threshold);
    }

    private static void pcSuper(int dimension, double xyRange, double gauss, double threshold) {
        Complex[][] superposition = new Complex[dimension][dimension];
        Complex[][] picture = new Complex[dimension][dimension];
        double[][] amplitude = new double[dimension][dimension];
        double[][] normAmplitude = new double[dimension][dimension];
        double stepXY = 2 * xyRange / dimension;
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        double phi1, phi2;
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                Complex one = Complex.valueOf(hermiteGauss2D(6, 0, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(1, 5));
                Complex two = Complex.valueOf(hermiteGauss2D(0, 6, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(5, 1));
                Complex three = Complex.valueOf(hermiteGauss2D(3, 3, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(10, 10));
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
                normAmplitude[i][j] = normFrom0to1(amplitude[i][j], min, max);
                if (normAmplitude[i][j] >= threshold) {
                    picture[i][j] = Complex.I.multiply(superposition[i][j].getArgument()).exp();
                } else {
                    phi1 = superposition[i][j].getArgument() + Math.acos(amplitude[i][j]);
                    phi2 = superposition[i][j].getArgument() - Math.acos(amplitude[i][j]);
                    Complex first = Complex.I.multiply(phi1).exp();
                    Complex second = Complex.I.multiply(phi2).exp();
                    picture[i][j] = Complex.valueOf(first.getReal() + second.getReal(), first.getImaginary() + second.getImaginary());
                }
            }
        }
        superposition(dimension, xyRange, gauss);
        Graph.draw2DIntensity(picture, "pictures/partIntSuper" + "threshold" + threshold + ".bmp");
        Graph.draw2DPhase(picture, "pictures/partPhaseSuper" + "threshold" + threshold + ".bmp");
    }

    private static void pcSuperRectangularI(int dimension, double xyRange, double gauss, double threshold) {
        Complex[][] superposition = new Complex[dimension][dimension];
        Complex[][] picture = new Complex[2 * dimension][dimension];
        double[][] amplitude = new double[dimension][dimension];
        double[][] normAmplitude = new double[dimension][dimension];
        double stepXY = 2 * xyRange / dimension;
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        double phi1, phi2;
        Complex temp;
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                Complex one = Complex.valueOf(hermiteGauss2D(6, 0, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(1, 5));
                Complex two = Complex.valueOf(hermiteGauss2D(0, 6, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(5, 1));
                Complex three = Complex.valueOf(hermiteGauss2D(3, 3, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(10, 10));
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
                normAmplitude[i][j] = normFrom0to1(amplitude[i][j], min, max);
                if (normAmplitude[i][j] >= threshold) {
                    temp = Complex.I.multiply(superposition[i][j].getArgument()).exp();
                    picture[2 * i][j] = temp;
                    picture[2 * i + 1][j] = temp;
                } else {
                    phi1 = superposition[i][j].getArgument() + Math.acos(amplitude[i][j]);
                    phi2 = superposition[i][j].getArgument() - Math.acos(amplitude[i][j]);
                    Complex first = Complex.I.multiply(phi1).exp();
                    Complex second = Complex.I.multiply(phi2).exp();
                    // picture[i][j] = Complex.valueOf(first.getReal() + second.getReal(), first.getImaginary() + second.getImaginary());
                    picture[2 * i][j] = first;
                    // picture[i][2 * j] = temp;
                    picture[2 * i + 1][j] = second;
                    // picture[i][2 * j + 1] = temp;
                }
            }
        }
        superposition(dimension, xyRange, gauss);
        Graph.draw2DIntensity(picture, "pictures/partIntSuper" + "threshold" + threshold + ".bmp");
        Graph.draw2DPhase(picture, "pictures/partPhaseSuper" + "threshold" + threshold + ".bmp");
    }

    private static void partialCoding(int dimension, double xyRange, int n, int m, double gauss, double threshold) {
        double step = 2 * xyRange / dimension;
        double[][] amplitude = new double[dimension][dimension];
        double[][] normAmplitude = new double[dimension][dimension];
        Complex[][] picture = new Complex[dimension][dimension];
        Complex[][] mode = new Complex[dimension][dimension];
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        double phi1, phi2;
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
                normAmplitude[i][j] = normFrom0to1(amplitude[i][j], min, max);
                if (normAmplitude[i][j] >= threshold) {
                    picture[i][j] = Complex.I.multiply(mode[i][j].getArgument()).exp();
                } else {
                    phi1 = mode[i][j].getArgument() + Math.acos(amplitude[i][j]);
                    phi2 = mode[i][j].getArgument() - Math.acos(amplitude[i][j]);
                    Complex first = Complex.I.multiply(phi1).exp();
                    Complex second = Complex.I.multiply(phi2).exp();
                    picture[i][j] = Complex.valueOf(first.getReal() + second.getReal(), first.getImaginary() + second.getImaginary());
                }
            }
        }
        mode2D(dimension, xyRange, n, m, gauss);
        Graph.draw2DIntensity(picture, "pictures/partInt" + "threshold" + threshold + ".bmp");
        Graph.draw2DPhase(picture, "pictures/partPhase" + "threshold" + threshold + ".bmp");
    }

    private static void partialCodingRectangularJ(int dimension, double xyRange, int n, int m, double gauss, double threshold) {
        double step = 2 * xyRange / dimension;
        double[][] amplitude = new double[dimension][dimension];
        double[][] normAmplitude = new double[dimension][dimension];
        Complex[][] picture = new Complex[dimension][2 * dimension];
        Complex[][] mode = new Complex[dimension][dimension];
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        double phi1, phi2;
        Complex temp;
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
                normAmplitude[i][j] = normFrom0to1(amplitude[i][j], min, max);
                if (normAmplitude[i][j] >= threshold) {
                    temp = Complex.I.multiply(mode[i][j].getArgument()).exp();
                    picture[i][2 * j] = temp;
                    picture[i][2 * j + 1] = temp;
                } else {
                    phi1 = mode[i][j].getArgument() + Math.acos(amplitude[i][j]);
                    phi2 = mode[i][j].getArgument() - Math.acos(amplitude[i][j]);
                    Complex first = Complex.I.multiply(phi1).exp();
                    Complex second = Complex.I.multiply(phi2).exp();
                    // temp = Complex.valueOf(first.getReal() + second.getReal(), first.getImaginary() + second.getImaginary());
                    picture[i][2 * j] = first;
                    // picture[i][2 * j] = temp;
                    picture[i][2 * j + 1] = second;
                    // picture[i][2 * j + 1] = temp;
                }
            }
        }
        mode2D(dimension, xyRange, n, m, gauss);
        Graph.draw2DIntensity(picture, "pictures/partInt" + "threshold" + threshold + "Rect.bmp");
        Graph.draw2DPhase(picture, "pictures/partPhase" + "threshold" + threshold + "Rect.bmp");
    }

    private static void partialCodingRectangularI(int dimension, double xyRange, int n, int m, double gauss, double threshold) {
        double step = 2 * xyRange / dimension;
        double[][] amplitude = new double[dimension][dimension];
        double[][] normAmplitude = new double[dimension][dimension];
        Complex[][] picture = new Complex[2 * dimension][dimension];
        Complex[][] mode = new Complex[dimension][dimension];
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        double phi1, phi2;
        Complex temp;
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
                normAmplitude[i][j] = normFrom0to1(amplitude[i][j], min, max);
                if (normAmplitude[i][j] >= threshold) {
                    temp = Complex.I.multiply(mode[i][j].getArgument()).exp();
                    picture[2 * i][j] = temp;
                    picture[2 * i + 1][j] = temp;
                } else {
                    phi1 = mode[i][j].getArgument() + Math.acos(amplitude[i][j]);
                    phi2 = mode[i][j].getArgument() - Math.acos(amplitude[i][j]);
                    Complex first = Complex.I.multiply(phi1).exp();
                    Complex second = Complex.I.multiply(phi2).exp();
                    // temp = Complex.valueOf(first.getReal() + second.getReal(), first.getImaginary() + second.getImaginary());
                    picture[2 * i][j] = first;
                    // picture[i][2 * j] = temp;
                    picture[2 * i + 1][j] = second;
                    // picture[i][2 * j + 1] = temp;
                }
            }
        }
        mode2D(dimension, xyRange, n, m, gauss);
        Graph.draw2DIntensity(picture, "pictures/partInt" + "threshold" + threshold + "Rect.bmp");
        Graph.draw2DPhase(picture, "pictures/partPhase" + "threshold" + threshold + "Rect.bmp");
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
        Complex[][] result = new Complex[dimension][dimension];
        Complex[][] phaseOnly = new Complex[dimension][dimension];
        double step = 2 * xyRange / dimension;
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                result[i][j] = new Complex(hermiteGauss2D(n, m, -xyRange + i * step, -xyRange + j * step, gauss));
                phaseOnly[i][j] = Complex.I.multiply(result[i][j].getArgument()).exp();
            }
        }
        Graph.draw2DIntensity(phaseOnly, "pictures/intensityPhaseGH" + n + m + " xy" + xyRange + " sigma " + gauss + ".bmp");
        Graph.draw2DPhase(phaseOnly, "pictures/phasePhaseGH" + n + m + " xy" + xyRange + " gauss " + gauss + ".bmp");
        Graph.draw2DIntensity(result, "pictures/intensityGH" + n + m + " xy" + xyRange + " sigma " + gauss + ".bmp");
        Graph.draw2DPhase(result, "pictures/phaseGH" + n + m + " xy" + xyRange + " gauss " + gauss + ".bmp");
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
    private static void fresnelTransform2D(int dimension, double uvRange, double xyRange, int n, int m, double gauss, double z, double k) {
        Complex[][] result = new Complex[dimension][dimension];
        double stepUV = 2 * uvRange / dimension;
        double stepXY = 2 * xyRange / dimension;
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                result[i][j] = FresnelTransform.transform2D(-uvRange + i * stepUV, -uvRange + j * stepUV, z, k,
                        xyRange, stepXY, n, m, gauss, dimension);
            }
        }
        mode2D(dimension, xyRange, n, m, gauss);
        Graph.draw2DIntensity(result, "pictures/intensityFresnelGH" + n + m + " xy" + xyRange + " gauss " + gauss + " uv" + uvRange + " z" + z + ".bmp");
        Graph.draw2DPhase(result, "pictures/phaseFresnelGH" + n + m + " xy" + xyRange + " gauss " + gauss + " uv" + uvRange + " z" + z + ".bmp");
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
    private static void fresnelTransform2DPhaseOnly(int dimension, double uvRange, double xyRange, int n, int m,
                                                    double gauss, double z, double k) {
        Complex[][] result = new Complex[dimension][dimension];
        double stepUV = 2 * uvRange / dimension;
        double stepXY = 2 * xyRange / dimension;
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                result[i][j] = FresnelTransform.transform2DPhase(-uvRange + i * stepUV, -uvRange + j * stepUV, z, k,
                        xyRange, stepXY, n, m, gauss, dimension);
            }
        }
        mode2D(dimension, xyRange, n, m, gauss);
        Graph.draw2DIntensity(result, "pictures/intensityFresnelPhaseGH" + n + m + " xy" + xyRange + " gauss " + gauss + " uv" + uvRange + " z" + z + " k" + k + ".bmp");
        Graph.draw2DPhase(result, "pictures/phaseFresnelPhaseGH" + n + m + " xy" + xyRange + " gauss " + gauss + " uv" + uvRange + " z" + z + " k" + k + ".bmp");
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
                result[i][j] = FourierTransform2D.transform2D(-uvRange + i * stepUV, -uvRange + j * stepUV, z, k,
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
    private static void fourierTransform2DPhaseOnly(int dimension, double uvRange, double xyRange, int n, int m, double gauss, double z, double k) {

        Complex[][] result = new Complex[dimension][dimension];
        double stepUV = 2 * uvRange / dimension;
        double stepXY = 2 * xyRange / dimension;

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                result[i][j] = FourierTransform2D.transform2DPhaseOnly(-uvRange + i * stepUV, -uvRange + j * stepUV, z, k,
                        -xyRange, -xyRange,
                        stepXY, stepXY,
                        n, m, gauss, dimension, dimension);
            }
        }
        mode2D(dimension, xyRange, n, m, gauss);
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
                Complex one = Complex.valueOf(hermiteGauss2D(6, 0, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(1, 5));
                Complex two = Complex.valueOf(hermiteGauss2D(0, 6, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(5, 1));
                Complex three = Complex.valueOf(hermiteGauss2D(3, 3, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(10, 10));
                result[i][j] = Complex.valueOf(one.getReal() + two.getReal() + three.getReal(), one.getImaginary() + two.getImaginary() + three.getImaginary());
                phaseOnly[i][j] = Complex.I.multiply(result[i][j].getArgument()).exp();
            }
        }
        Graph.draw2DIntensity(phaseOnly, "pictures/intensitySuperpositionPhaseOnly.bmp");
        Graph.draw2DPhase(phaseOnly, "pictures/phaseSuperpositionPhaseOnly.bmp");
        Graph.draw2DIntensity(result, "pictures/intensitySuperposition.bmp");
        Graph.draw2DPhase(result, "pictures/phaseSuperposition.bmp");
    }

    private static void fresnelTransform2DSuperposition(int dimension, double uvRange, double xyRange, double gauss,
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
        Graph.draw2DIntensity(result, "pictures/intensitySuperpositionOutputFresnel.bmp");
        Graph.draw2DPhase(result, "pictures/phaseSuperpositionOutputFresnel.bmp");
    }

    private static void fresnelTransform2DSuperpositionPhaseOnly(int dimension, double uvRange, double xyRange,
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
        superposition(dimension, xyRange, gauss);
        Graph.draw2DIntensity(function, "pictures/intensitySuperpositionPhaseOnlyOutputByFresnel.bmp");
        Graph.draw2DPhase(function, "pictures/phaseSuperpositionPhaseOnlyOutputByFresnel.bmp");
    }

    private static Complex[][] superpos(int dimension, double xyRange, double gauss, double stepXY) {
        Complex[][] superposition = new Complex[dimension][dimension];
        Complex[][] one = new Complex[dimension][dimension];
        Complex[][] two = new Complex[dimension][dimension];
        Complex[][] three = new Complex[dimension][dimension];
        for (int i = 0; i < superposition.length; i++) {
            for (int j = 0; j < superposition[0].length; j++) {
                one[i][j] = Complex.valueOf(hermiteGauss2D(6, 0, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(1, 5));
                two[i][j] = Complex.valueOf(hermiteGauss2D(0, 6, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(5, 1));
                three[i][j] = Complex.valueOf(hermiteGauss2D(3, 3, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(10, 10));
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
                function[i][j] = FourierTransform2D.transform2DSuperposition(-uvRange + i * stepUV, -uvRange + j * stepUV, z, k,
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
                function[i][j] = FourierTransform2D.
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