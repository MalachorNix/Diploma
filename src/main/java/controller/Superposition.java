package controller;

import model.GHMode;
import org.apache.commons.math3.complex.Complex;

public final class Superposition {

    private Superposition() {

    }

    public static Complex superposition1D(int N, double x, double gauss, Complex[] coefficient) {
        Complex result = new Complex(0);

        for (int n = 0; n <= N; n++) {
            result.add(coefficient[n].multiply(GHMode.hermiteGauss1D(n, x, gauss)));
        }

        return result;
    }


    public static Complex superposition2D(int N, double x, double y, double gauss, Complex[][] coefficient) {
        Complex result = new Complex(0);
        Complex multi;

        /*for (int m = 0; m <= N; m++) {
            for (int n = 0; n <= N; n++) {
                // multi = coefficient[m][n].multiply(GHMode.hermiteGauss2D(n, m, x, y, gauss));
                multi = coefficient[m][n].multiply(GHMode.hermiteGauss2D(m, n, x, y, gauss));
                result = new Complex(result.getReal() + multi.getReal(), result.getImaginary() + multi.getImaginary());
            }
        }*/

        for (int m = 0; m <= N; m++) {
            multi = coefficient[m][m].multiply(GHMode.hermiteGauss2D(m, m, x, y, gauss));
            result = new Complex(result.getReal() + multi.getReal(), result.getImaginary() + multi.getImaginary());
        }

        return result;
    }
}