package controller;

import model.GHMode;
import org.apache.commons.math3.complex.Complex;

public class FourierTransform {

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

    public static Complex transform2DPhaseOnly(double u, double v, double z, double k,
                                           double yMin, double xMin,
                                           double stepX, double stepY,
                                           int n, int m, double gauss,
                                           int width, int height) {

        Complex first = firstMultiplier2DTransform(k, z);
        Complex second = firstExponent(k, z);
        Complex third = integral2DPhase(u, v, z, k, yMin, xMin, stepX, stepY, n, m, gauss, width, height);
        return first.multiply(second).multiply(third);
    }

    private static Complex firstMultiplier2DTransform(double k, double z) {
        return Complex.I.multiply(-(k / (2 * Math.PI * z)));
    }

    private static Complex firstExponent(double k, double z) {
        return Complex.I.multiply(k * z).exp();
    }

    private static Complex integral2D(double u, double v, double z, double k,
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

    private static Complex integrand2DPhase(Complex mode, double k, double z, double x, double y,
                                            double u, double v) {
        Complex phaseOnly = Complex.I.multiply(mode.getArgument()).exp();
        return phaseOnly.multiply(Complex.I.multiply(-k / z * (x * u + y * v)).exp());
    }

    public static Complex transform2DSuperposition(double u, double v, double z, double k,
                                                   double yMin, double xMin,
                                                   double stepX, double stepY,
                                                   double gauss,
                                                   int width, int height) {

        Complex first = Complex.I.multiply(-(k / (2 * Math.PI * z)));
        Complex second = Complex.I.multiply(k * z).exp();
        Complex third = integral2DSuperposition(u, v, z, k, yMin, xMin, stepX, stepY, gauss, width, height);

        return first.multiply(second).multiply(third);
    }

    private static Complex integral2DSuperposition(double u, double v, double z, double k,
                                                   double yMin, double xMin,
                                                   double stepX, double stepY,
                                                   double gauss,
                                                   int width, int height) {

        Complex sum = new Complex(0, 0);
        Complex multi;

        Complex[][] one = new Complex[width][height];
        Complex[][] two = new Complex[width][height];
        Complex[][] three = new Complex[width][height];
        Complex[][] function = new Complex[width][height];

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                one[i][j] = Complex.valueOf(GHMode.hermiteGauss2D(6, 0, xMin + i * stepX, yMin + j * stepY, gauss)).multiply(new Complex(1));
                two[i][j] = Complex.valueOf(GHMode.hermiteGauss2D(0, 6, xMin + i * stepX, yMin + j * stepY, gauss)).multiply(new Complex(5));
                three[i][j] = Complex.valueOf(GHMode.hermiteGauss2D(3, 3, xMin + i * stepX, yMin + j * stepY, gauss)).multiply(new Complex(10));
                function[i][j] = new Complex(one[i][j].getReal() + two[i][j].getReal() + three[i][j].getReal(), one[i][j].getImaginary() + two[i][j].getImaginary() + three[i][j].getImaginary());
                multi = integrand2D(function[i][j],
                        k, z, xMin + i * stepX, yMin + j * stepY, u, v).
                        multiply(stepX).multiply(stepY);
                sum = new Complex(sum.getReal() + multi.getReal(),
                        sum.getImaginary() + multi.getImaginary());
            }
        }

        return sum;
    }

    public static Complex transform2DPhaseOnly(double u, double v, double z, double k,
                                               double yMin,
                                               double xMin,
                                               double stepX, double stepY,
                                               int width, int height, Complex[][] phaseOnly) {

        Complex first = Complex.I.multiply(-(k / (2 * Math.PI * z)));
        Complex second = Complex.I.multiply(k * z).exp();
        Complex third = integral2DPhaseOnly(u, v, z, k, yMin, xMin, stepX, stepY, phaseOnly, width, height);

        return first.multiply(second).multiply(third);
    }

    private static Complex integral2DPhaseOnly(double u, double v, double z, double k,
                                               double yMin,
                                               double xMin,
                                               double stepX, double stepY,
                                               Complex[][] phaseOnly,
                                               int width, int height) {

        Complex sum = new Complex(0, 0);
        Complex multi;

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                multi = integrand2D(phaseOnly[i][j],
                        k, z, xMin + i * stepX, yMin + j * stepY, u, v).
                        multiply(stepX).multiply(stepY);

                sum = new Complex(sum.getReal() + multi.getReal(),
                        sum.getImaginary() + multi.getImaginary());
            }
        }

        return sum;
    }

    private static Complex integrand2D(Complex function, double k, double z, double x, double y, double u, double v) {
        return function.multiply(Complex.I.multiply(-k / z * (x * u + y * v)).exp());
    }

}
