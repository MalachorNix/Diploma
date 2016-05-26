package controller;

import model.GHMode;
import org.apache.commons.math3.complex.Complex;

public final class FresnelTransform {

    private FresnelTransform() {

    }

    public static Complex transform1D(double u, double z, double k, int n, double a, double step, double gauss, int width) {
        Complex first = firstMultiplier1DTransform(k, z);
        Complex second = firstExponent(k, z);
        Complex third = integral1D(u, z, k, a, step, n, gauss, width);
        return first.multiply(second).multiply(third);
    }

    public static Complex transform2D(double u, double v, double z, double k,
                                      double yMin, double xMin,
                                      double stepX, double stepY,
                                      int n, int m, double gauss,
                                      int width, int height) {

        Complex first = firstMultiplier2DTransform(k, z);
        Complex second = firstExponent(k, z);
        Complex third = integral2D(u, v, z, k, yMin, xMin, stepX, stepY, n, m, gauss, width, height);
        return first.multiply(second).multiply(third);
    }

    public static Complex transform2DPhase(double u, double v, double z, double k,
                                      double yMin, double xMin,
                                      double stepX, double stepY,
                                      int n, int m, double gauss,
                                      int width, int height) {

        Complex first = firstMultiplier2DTransform(k, z);
        Complex second = firstExponent(k, z);
        Complex third = integral2DPhase(u, v, z, k, yMin, xMin, stepX, stepY, n, m, gauss, width, height);
        return first.multiply(second).multiply(third);
    }

    public static Complex transform2DSuperposition(double u, double v, double z, double k,
                                                   double yMin,
                                                   double xMin,
                                                   double stepX, double stepY,
                                                   Complex superposition, int width, int height) {

        Complex first = firstMultiplier2DTransform(k, z);
        Complex second = firstExponent(k, z);
        Complex third = integral2DSuperposition(u, v, z, k, yMin, xMin, stepX, stepY, superposition, width, height);

        return first.multiply(second).multiply(third);
    }

    public static Complex transform2DSuperposition(double u, double v, double z, double k,
                                                   double yMin, double xMin,
                                                   double stepX, double stepY,
                                                   int N, double gauss, Complex[][] coefficient,
                                                   int width, int height) {

        Complex first = firstMultiplier2DTransform(k, z);
        Complex second = firstExponent(k, z);
        Complex third = integral2DSuperposition(u, v, z, k, yMin, xMin, stepX, stepY, N, gauss, coefficient, width, height);

        return first.multiply(second).multiply(third);
    }

    public static Complex transform2DPhaseOnly(double u, double v, double z, double k,
                                                   double yMin, double yMax,
                                                   double xMin, double xMax,
                                                   double stepX, double stepY,
                                                   int N, double gauss,
                                                   int width, int height, Complex[][] phaseOnly, Complex[][] coefficient) {

        Complex first = firstMultiplier2DTransform(k, z);
        Complex second = firstExponent(k, z);
        Complex third = integral2DPhaseOnly(u, v, z, k, yMin, yMax, xMin, xMax, stepX, stepY, N, gauss, phaseOnly, width, height, coefficient);

        return first.multiply(second).multiply(third);
    }

    public static Complex firstMultiplier1DTransform(double k, double z) {
        return Complex.I.multiply(-(k / (2 * Math.PI * z))).sqrt();
    }

    public static Complex firstMultiplier2DTransform(double k, double z) {
        return Complex.I.multiply(-(k / (2 * Math.PI * z)));
    }

    public static Complex firstExponent(double k, double z) {
        return Complex.I.multiply(k * z).exp();
    }

    public static Complex integral1D(double u, double z, double k, double xMin, double step, int n, double gauss, int width) {

        Complex sum = new Complex(0, 0);
        Complex multi;

        for (int i = 1; i < width; i++) {
            Complex function = new Complex(GHMode.hermiteGauss1D(n, xMin + i * step, gauss));
            multi = integrand1D(function, u, z, k, xMin + i * step).multiply(step);
            sum.add(multi);
        }
        return sum;
    }

    public static Complex integral2D(double u, double v, double z, double k,
                                     double yMin, double xMin,
                                     double stepX, double stepY,
                                     int n, int m, double gauss,
                                     int width, int height) {

        Complex sum = new Complex(0, 0);
        Complex multi;

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                Complex function = new Complex(GHMode.hermiteGauss2D(n, m, xMin + i * stepX, yMin + j * stepY, gauss));

                multi = integrand2D(function, k, z, xMin + i * stepX, yMin + j * stepY, u, v).multiply(stepX).multiply(stepY);

                sum = new Complex(sum.getReal() + multi.getReal(), sum.getImaginary() + multi.getImaginary());
            }
        }

        return sum;
    }

    private static Complex integral2DPhase(double u, double v, double z, double k,
                                           double yMin, double xMin,
                                           double stepX, double stepY,
                                           int n, int m, double gauss,
                                           int width, int height) {

        Complex sum = new Complex(0, 0);
        Complex multi;

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                Complex mode = new Complex(GHMode.hermiteGauss2D(n, m, xMin + i * stepX, yMin + j * stepY, gauss));

                multi = integrand2DPhase(mode, k, z, xMin + i * stepX, yMin + j * stepY, u, v).multiply(stepX).multiply(stepY);

                sum = new Complex(sum.getReal() + multi.getReal(), sum.getImaginary() + multi.getImaginary());
            }
        }

        return sum;
    }

    private static Complex integral2DSuperposition(double u, double v, double z, double k,
                                                   double yMin, double xMin,
                                                   double stepX, double stepY,
                                                   Complex superposition, int width, int height) {

        Complex sum = new Complex(0, 0);
        Complex multi;

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                multi = integrand2D(superposition, k, z, xMin + i * stepX, yMin + j * stepY, u, v).multiply(stepX).multiply(stepY);

                sum = new Complex(sum.getReal() + multi.getReal(), sum.getImaginary() + multi.getImaginary());
            }
        }

        return sum;
    }

    public static Complex integral2DSuperposition(double u, double v, double z, double k,
                                                  double yMin, double xMin,
                                                  double stepX, double stepY,
                                                  int N, double gauss, Complex[][] coefficient,
                                                  int width, int height) {

        Complex sum = new Complex(0, 0);
        Complex multi;

        Complex[][] one = new Complex[width][height];
        Complex[][] two = new Complex[width][height];
        Complex[][] three = new Complex[width][height];
        Complex[][] function = new Complex[width][height];

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                one[i][j] = Complex.valueOf(GHMode.hermiteGauss2D(6, 0, xMin + i * stepX, yMin + j * stepY, gauss)).multiply(new Complex(1, 5));
                two[i][j] = Complex.valueOf(GHMode.hermiteGauss2D(0, 6, xMin + i * stepX, yMin + j * stepY, gauss)).multiply(new Complex(5, 1));
                three[i][j] = Complex.valueOf(GHMode.hermiteGauss2D(3, 3, xMin + i * stepX, yMin + j * stepY, gauss)).multiply(new Complex(10, 10));
                function[i][j] = new Complex(one[i][j].getReal() + two[i][j].getReal() + three[i][j].getReal(), one[i][j].getImaginary() + two[i][j].getImaginary() + three[i][j].getImaginary());
                multi = integrand2D(function[i][j],
                        k, z, xMin + i * stepX, yMin + j * stepY, u, v).
                        multiply(stepX).multiply(stepY);
                /*multi = integrand2D(superposition2D(N, xMin + i * stepX, yMin + j * stepY, gauss, coefficient),
                        k, z, xMin + i * stepX, yMin + j * stepY, u, v).
                        multiply(stepX).multiply(stepY);*/
                sum = new Complex(sum.getReal() + multi.getReal(),
                        sum.getImaginary() + multi.getImaginary());
            }
        }

        return sum;
    }

    public static Complex integral2DPhaseOnly(double u, double v, double z, double k,
                                                  double yMin, double yMax,
                                                  double xMin, double xMax,
                                                  double stepX, double stepY,
                                                  int N, double gauss, Complex[][] phaseOnly,
                                                  int width, int height, Complex[][] coefficient) {

        Complex sum = new Complex(0, 0);
        Complex multi;

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                multi = integrand2D(phaseOnly[i][j] /*new Complex(superposition2D(N, xMin + i * stepX, yMin + j * stepY, gauss, coefficient).getArgument())*/,
                        k, z, xMin + i * stepX, yMin + j * stepY, u, v).
                        multiply(stepX).multiply(stepY);
                sum = new Complex(sum.getReal() + multi.getReal(),
                        sum.getImaginary() + multi.getImaginary());
            }
        }

        return sum;
    }

    public static Complex integrand1D(Complex function, double k, double z, double x, double u) {
        return function.multiply(transformExponent1D(k, z, x, u));
    }

    public static Complex integrand2D(Complex function, double k, double z, double x, double y,
                                      double u, double v) {
        return function.multiply(transformExponent2D(k, z, x, y, u, v));
    }

    private static Complex integrand2DPhase(Complex mode, double k, double z, double x, double y,
                                            double u, double v) {
        Complex phaseOnly = Complex.I.multiply(mode.getArgument()).exp();

        return transformExponent2D(k, z, x, y, u, v).multiply(phaseOnly);
    }

    public static Complex transformExponent1D(double k, double z, double x, double u) {
        return Complex.I.multiply(k * (x - u) * (x - u) / (2 * z)).exp();
    }

    public static Complex transformExponent2D(double k, double z, double x, double y,
                                              double u, double v) {
        return Complex.I.multiply(k * ((x - u) * (x - u) + (y - v) * (y - v)) / (2 * z)).exp();
        // return Complex.I.multiply(-k / z * (x * u + y * v)).exp();
    }
}