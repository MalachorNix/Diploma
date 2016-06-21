package view;

import controller.FresnelTransform;
import org.apache.commons.math3.complex.Complex;
import org.junit.Test;

import static org.junit.Assert.*;

public class GraphTest {
    @Test
    public void draw2DIntensity() throws Exception {
        /*int N = 10;
        Complex[][] coefficient = new Complex[N][N];

        for (int i = 0; i < coefficient.length; i++) {
            for (int j = 0; j < coefficient[0].length; j++) {
                coefficient[i][j] = new Complex(i, j);
            }
        }*/

        /*Complex superposition = Superposition.superposition2D(N, 4, 7, 11, coefficient);

        Complex[][] function = new Complex[100][100];

        double u = -1, v = -1;

        for (int i = 0; i < function.length; i++, u+=0.02) {
            for (int j = 0; j < function[0].length; j++, v+=0.02) {
                function[i][j] = FresnelTransform.transform2DSuperposition(u, v, 100, 4, -1, 1, -1, 1, 0.02, 0.02, superposition);
            }
            v = -1;
        }

        new Graph().draw2DIntensity(function);*/


        /*Complex[][] superposition = new Complex[100][100];

        double x = -1, y = -1;

        for (int i = 0; i < superposition.length; i++, x+=0.02) {
            for (int j = 0; j < superposition[0].length; j++, y+=0.02) {
                superposition[i][j] = Superposition.superposition2D(N, x, y, 11, coefficient);
            }
            y = -1;
        }


        new Graph().draw2DIntensity(superposition);*/
    }

    @Test
    public void draw2DPhase() throws Exception {
        /*int N = 10;
        Complex[][] coefficient = new Complex[N][N];

        for (int i = 0; i < coefficient.length; i++) {
            for (int j = 0; j < coefficient[0].length; j++) {
                coefficient[i][j] = new Complex(i, j);
            }
        }

        Complex superposition = Superposition.superposition2D(N, 4, 7, 11, coefficient);

        Complex[][] function = new Complex[128][128];

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                function[i][j] = FresnelTransform.transform2DSuperposition(i, j, 100, 4, 5, 6, 7, 8, 0.01, 0.01, superposition);
            }
        }

        new Graph().draw2DPhase(function);*/

        /*Complex[][] superposition = new Complex[128][128];

        for (int i = 0; i < superposition.length; i++) {
            for (int j = 0; j < superposition[0].length; j++) {
                superposition[i][j] = Superposition.superposition2D(N, i, j, 11, coefficient);
            }
        }


        new Graph().draw2DIntensity(superposition);*/
    }

}