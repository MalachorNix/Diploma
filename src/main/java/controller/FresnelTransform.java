package controller;

import model.HermiteGaussianModes;
import org.apache.commons.math3.complex.Complex;

public final class FresnelTransform {

    private FresnelTransform() {

    }

    public static Complex transform1D(double u, double z, double k, int n, double a, double b,
                                      double step, double gauss) {
        Complex first = firstMultiplier1DTransform(k, z);
        Complex second = firstExponent(k, z);
        Complex third = integral1D(u, z, k, a, b, step, n, gauss);
        return first.multiply(second).multiply(third);
    }

    public static Complex transform2D(double u, double v, double z, double k, double a, double b, double c, double d,
                                      double stepX, double stepY, int n, int m, double gauss) {
        Complex first = firstMultiplier2DTransform(k, z);
        Complex second = firstExponent(k, z);
        Complex third = integral2D(u, v, z, k, a, b, c, d, stepX, stepY, n, m, gauss);
        return first.multiply(second).multiply(third);
    }

    public static Complex transform2DSuperposition(double u, double v, double z, double k, double a, double b, double c, double d,
                                                   double stepX, double stepY, Complex superposition) {
        Complex result = new Complex(0);
        Complex first = firstMultiplier2DTransform(k, z);
        Complex second = firstExponent(k, z);
        Complex third = integral2DSuperposition(u, v, z, k, a, b, c, d, stepX, stepY, superposition);

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

    public static Complex integral1D(double u, double z, double k, double a, double b,
                                      double step, int n, double gauss) {

        Complex sum = new Complex(0, 0);
        Complex multi;

        double[] x = new double[(int) ((b - a) / step)];

        for (int i = 0; i < x.length; i++) { // Если прижмет, то переделать без массива.
            x[i] = a + step * i;
        }

        for (int i = 1; i < x.length; i++) {
            Complex function = new Complex(model.HermiteGaussianModes.hermiteGauss1D(n, x[i], gauss));
            multi = integrand1D(function, u, z, k, x[i]).multiply(step);
            sum.add(multi);
        }
        return sum;
    }

    public static Complex integral2D(double u, double v, double z, double k, double yMin, double yMax, double xMin, double xMax,
                                      double stepX, double stepY, int n, int m, double gauss) {

        Complex sum = new Complex(0, 0);
        Complex multi;

        double[] x = new double[(int) ((xMax - xMin) / stepX) + 1];
        double[] y = new double[(int) ((yMax - yMin) / stepY) + 1];

        for (int i = 0; i < x.length; i++) { // Если прижмет, то переделать без массива.
            x[i] = xMin + stepX * i;
        }

        for (int i = 0; i < y.length; i++) {
            y[i] = yMin + stepY * i;
        }

        for (int i = 0; i < x.length; i++) {
            for (int j = 0; j < y.length; j++) {
                Complex function = new Complex(model.HermiteGaussianModes.hermiteGauss2D(n, m, x[i], y[j], gauss));
                multi = integrand2D(function, k, z, x[i], y[j], u, v).multiply(stepX).multiply(stepY);
                sum = new Complex(sum.getReal() + multi.getReal(), sum.getImaginary() + multi.getImaginary());
            }
        }

        return sum;
    }

    public static Complex integral2DSuperposition(double u, double v, double z, double k, double yMin, double yMax, double xMin, double xMax,
                                      double stepX, double stepY, Complex superposition) {

        Complex sum = new Complex(0, 0);
        Complex multi;

        double[] x = new double[(int) ((xMax - xMin) / stepX) + 1];
        double[] y = new double[(int) ((yMax - yMin) / stepY) + 1];

        for (int i = 0; i < x.length; i++) { // Если прижмет, то переделать без массива.
            x[i] = xMin + stepX * i;
        }

        for (int i = 0; i < y.length; i++) {
            y[i] = yMin + stepY * i;
        }

        for (int i = 0; i < x.length; i++) {
            for (int j = 0; j < y.length; j++) {
                multi = integrand2D(superposition, k, z, x[i], y[j], u, v).multiply(stepX).multiply(stepY);
                sum = new Complex(sum.getReal() + multi.getReal(), sum.getImaginary() + multi.getImaginary());
            }
        }

        return sum;
    }

    public static Complex integrand1D(Complex function, double k, double z, double x, double u) {
        return function.multiply(transformExponent1D(k, z, x, u));
    }

    public static Complex integrand2D(Complex function, double k, double z, double x, double y, double u, double v) {
        return function.multiply(transformExponent2D(k, z, x, y, u, v));
    }

    public static Complex transformExponent1D(double k, double z, double x, double u) {
        return Complex.I.multiply(k * (x - u) * (x - u) / (2 * z)).exp();
    }

    public static Complex transformExponent2D(double k, double z, double x, double y, double u, double v) {
        return Complex.I.multiply(k * ((x - u) * (x - u) + (y - v) * (y - v)) / (2 * z)).exp();
    }
}