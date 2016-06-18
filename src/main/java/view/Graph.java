package view;

import org.apache.commons.math3.complex.Complex;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;


public class Graph {

    public static void draw2DIntensity(Complex[][] function, String filename) {
        double[][] intensity = new double[function.length][function[0].length];
        double minAmplitude = Double.MAX_VALUE;
        double maxAmplitude = Double.MIN_VALUE;

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                intensity[i][j] = function[i][j].abs() * function[i][j].abs();

                if (intensity[i][j] > maxAmplitude) {
                    maxAmplitude = intensity[i][j];
                }

                if (intensity[i][j] < minAmplitude) {
                    minAmplitude = intensity[i][j];
                }

            }
        }

        writeImage(intensity, maxAmplitude, minAmplitude, filename);

    }

    public static void draw2DPhase(Complex[][] function, String filename) {
        double[][] phase = new double[function.length][function[0].length];
        double minPhase = Double.MAX_VALUE;
        double maxPhase = Double.MIN_VALUE;

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                phase[i][j] = function[i][j].getArgument();

                if (phase[i][j] < 0) {
                    phase[i][j] += 2 * Math.PI;
                }

                if (phase[i][j] > maxPhase) {
                    maxPhase = phase[i][j];
                }

                if (phase[i][j] < minPhase) {
                    minPhase = phase[i][j];
                }

            }
        }

        writeImage(phase, 2 * Math.PI, 0, filename);
    }

    private static void writeImage(double[][] function, double max, double min, String filename) {
        double stepNorm = (max - min) / 255;


        BufferedImage image = new BufferedImage(function.length, function[0].length, BufferedImage.TYPE_INT_RGB);

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                for (int k = 1; k < 255; k++) {
                    if (function[i][j] > min + stepNorm * k) {
                        image.setRGB(i, j, new Color(k - 1, k - 1, k - 1).getRGB());
                    }
                    if (function[i][j] == max) {
                        image.setRGB(i, j, new Color(255, 255, 255).getRGB());
                    }
                }
            }
        }

        try {
            ImageIO.write(image, "bmp", new File(filename));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static Complex[][] phaseOnlyEncode(Complex[][] superposition) {

        Complex[][] result = new Complex[superposition.length][superposition[0].length];

        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[0].length; j++) {
                result[i][j] = Complex.I.multiply(superposition[i][j].getArgument()).exp();
            }
        }


        return result;
    }
}