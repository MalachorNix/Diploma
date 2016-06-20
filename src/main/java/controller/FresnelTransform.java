package controller;

import model.GHMode;
import org.apache.commons.math3.complex.Complex;

public final class FresnelTransform {

    private FresnelTransform() {

    }

    public static Complex transform2D(double u, double v, double z, double k, double xyRange, double stepXY,
                                      int n, int m, double gauss, int dimension, Complex[][] mode) {
        Complex first = Complex.I.multiply(-(k / (2 * Math.PI * z)));
        Complex second = Complex.I.multiply(k * z).exp();
        Complex third = integral2D(u, v, z, k, xyRange, stepXY, n, m, gauss, dimension, mode);
        return first.multiply(second).multiply(third);
    }

    public static Complex transform2DPhase(double u, double v, double z, double k, double xyRange, double stepXY,
                                           int n, int m, double gauss, int dimension, Complex[][] mode) {
        Complex first = Complex.I.multiply(-(k / (2 * Math.PI * z)));
        Complex second = Complex.I.multiply(k * z).exp();
        Complex third = integral2DPhase(u, v, z, k, xyRange, stepXY, n, m, gauss, dimension, mode);
        return first.multiply(second).multiply(third);
    }

    public static Complex transform2DSuperposition(double u, double v, double z, double k, double xyRange,
                                                   double stepXY, double gauss, int dimension) {
        Complex first = Complex.I.multiply(-(k / (2 * Math.PI * z)));
        Complex second = Complex.I.multiply(k * z).exp();
        Complex third = integral2DSuperposition(u, v, z, k, xyRange, stepXY, gauss, dimension);
        return first.multiply(second).multiply(third);
    }

    public static Complex transform2DPhaseOnly(double u, double v, double z, double k, double xyRange, double stepXY,
                                               int dimension, Complex[][] phaseOnly) {
        Complex first = Complex.I.multiply(-(k / (2 * Math.PI * z)));
        Complex second = Complex.I.multiply(k * z).exp();
        Complex third = integral2DPhaseOnly(u, v, z, k, xyRange, stepXY, dimension, phaseOnly);
        return first.multiply(second).multiply(third);
    }

    private static Complex integral2D(double u, double v, double z, double k, double xyRange, double stepXY,
                                      int n, int m, double gauss, int dimension, Complex[][] function) {
        Complex sum = new Complex(0, 0);
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                // Complex function = Complex.valueOf(GHMode.hermiteGauss2D(n, m, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss));
                Complex multi = integrand2D(function[i][j], k, z, -xyRange + i * stepXY, -xyRange + j * stepXY, u, v).multiply(stepXY * stepXY);
                // Complex multi = integrand2D(function, k, z, -xyRange + i * stepXY, -xyRange + j * stepXY, u, v).multiply(stepXY * stepXY);
                sum = Complex.valueOf(sum.getReal() + multi.getReal(), sum.getImaginary() + multi.getImaginary());
            }
        }
        return sum;
    }

    private static Complex integral2DPhase(double u, double v, double z, double k, double xyRange, double stepXY,
                                           int n, int m, double gauss, int dimension, Complex[][] mode) {
        Complex sum = new Complex(0, 0);
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                // Complex mode = Complex.valueOf(GHMode.hermiteGauss2D(n, m, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss));
                // Complex multi = integrand2DPhase(mode, k, z, -xyRange + i * stepXY, -xyRange + j * stepXY, u, v).multiply(stepXY * stepXY);
                Complex multi = integrand2DPhase(mode[i][j], k, z, -xyRange + i * stepXY, -xyRange + j * stepXY, u, v).multiply(stepXY * stepXY);
                sum = Complex.valueOf(sum.getReal() + multi.getReal(), sum.getImaginary() + multi.getImaginary());
            }
        }
        return sum;
    }

    private static Complex integral2DSuperposition(double u, double v, double z, double k, double xyRange,
                                                   double stepXY, double gauss, int dimension) {
        Complex sum = new Complex(0, 0);
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                Complex one = Complex.valueOf(GHMode.hermiteGauss2D(6, 0, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(1));
                Complex two = Complex.valueOf(GHMode.hermiteGauss2D(0, 6, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(5));
                Complex three = Complex.valueOf(GHMode.hermiteGauss2D(3, 3, -xyRange + i * stepXY, -xyRange + j * stepXY, gauss)).multiply(Complex.valueOf(10));
                Complex function = Complex.valueOf(one.getReal() + two.getReal() + three.getReal(), one.getImaginary() + two.getImaginary() + three.getImaginary());
                Complex multi = integrand2D(function, k, z, -xyRange + i * stepXY, -xyRange + j * stepXY, u, v).
                        multiply(stepXY * stepXY);
                sum = Complex.valueOf(sum.getReal() + multi.getReal(), sum.getImaginary() + multi.getImaginary());
            }
        }
        return sum;
    }

    private static Complex integral2DPhaseOnly(double u, double v, double z, double k, double xyRange, double stepXY,
                                               int dimension, Complex[][] phaseOnly) {
        Complex sum = new Complex(0, 0);
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                Complex multi = integrand2D(phaseOnly[i][j], k, z, xyRange + i * stepXY, xyRange + j * stepXY, u, v)
                        .multiply(stepXY * stepXY);
                sum = Complex.valueOf(sum.getReal() + multi.getReal(), sum.getImaginary() + multi.getImaginary());
            }
        }
        return sum;
    }

    public static Complex integrand2D(Complex function, double k, double z, double x, double y, double u, double v) {
        double exp = Math.exp(-(x * x + y * y) / (2 * 0.6 * 0.6));
        return function.multiply(transformExponent2D(k, z, x, y, u, v)).multiply(exp);
        // return function.multiply(transformExponent2D(k, z, x, y, u, v));
    }

    private static Complex integrand2DPhase(Complex mode, double k, double z, double x, double y, double u, double v) {
        double exp = Math.exp(-(x * x + y * y) / (2 * 0.6 * 0.6));
        Complex phaseOnly = Complex.I.multiply(mode.getArgument()).exp();
        // return phaseOnly.multiply(transformExponent2D(k, z, x, y, u, v));
        return phaseOnly.multiply(transformExponent2D(k, z, x, y, u, v)).multiply(exp);
    }

    public static Complex transformExponent2D(double k, double z, double x, double y, double u, double v) {
        return Complex.I.multiply(k * ((x - u) * (x - u) + (y - v) * (y - v)) / (2 * z)).exp();
    }
}