package controller;

import model.HermiteGaussianModes;
import org.apache.commons.math3.complex.Complex;


public final class FourierTransform1D {

    private FourierTransform1D() {

    }

    public static Complex full1DTransform(double u, double z, double a, double b, double step,
                                          double ksi, int n, double gauss, double k) {
        return firstMultiplier1DTransform(k, z).
                multiply(firstExponent(k, z).
                        multiply(secondExponent(k, u, z)).
                        multiply(transform1D(a, b, step, ksi, n, gauss)));
    }

    private static Complex firstMultiplier1DTransform(double k, double z) {
        return Complex.I.multiply(-(k / (2 * Math.PI * z))).sqrt();
    }

    private static Complex firstExponent(double k, double z) {
        return Complex.I.multiply(k * z).exp();
    }

    private static Complex secondExponent(double k, double u, double z) {
        return Complex.I.multiply((k * u * u) / (2 * z)).exp();
    }

    /**
     * Считает аналоговое одномерное преобразование Фурье, интегрируя методом правых прямоугольников.
     *
     * @param a    Нижняя граница интегрирования.
     * @param b    Верхняя граница интегрирования.
     * @param step Шаг разбиения.
     * @param ksi  Аргумент ξ.
     * @return Значение функции аналогового преобразования Фурье.
     */
    private static Complex transform1D(double a, double b, double step, double ksi, int n, double gauss) {

        Complex sum = new Complex(0, 0);
        Complex multi;

        double[] x = new double[(int) ((b - a) / step)];

        for (int i = 0; i < x.length; i++) { // Если прижмет, то переделать без массива.
            x[i] = a + step * i;
        }

        for (int i = 1; i < x.length; i++) {
            Complex function = new Complex(HermiteGaussianModes.hermiteGauss1D(n, x[i], gauss));
            multi = integrand(x[i], ksi, function).multiply(step);
            // sum = new Complex(sum.getReal() + multi.getReal(), sum.getImaginary() + multi.getImaginary());
            sum.add(multi);
        }

        return sum;
    }

    /**
     * Считает выражение подынтегральное выражение в формуле преобразования.
     *
     * @param x        Аругмент x.
     * @param ksi      Аргумент ξ.
     * @param function Значение функции, которую будем преобразовывать.
     * @return Значение подынтегрального выражения формулы преобразования.
     */
    private static Complex integrand(double x, double ksi, Complex function) {
        return function.multiply(exponent(x, ksi));
    }

    /**
     * Считает экспоненту, находящуюся под знаком интгерала, в формуле преобразования.
     *
     * @param x   Аргумент x.
     * @param ksi Аргумент ξ.
     * @return Значение экспоненты формулы преобразования.
     */
    private static Complex exponent(double x, double ksi) {
        return Complex.I.multiply(-2 * Math.PI * x * ksi).exp();
    }
}