package view;

import org.apache.commons.math3.complex.Complex;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;


public class Graph {

    public static int[][] draw2DIntensity(Complex[][] function, String filename) {
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

        return writeImage(intensity, maxAmplitude, minAmplitude, filename);

    }

    public static int[][] draw2DPhase(Complex[][] function, String filename) {
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

        return writeImage(phase, 2 * Math.PI, 0, filename);
    }

    private static int[][] writeImage(double[][] function, double max, double min, String filename) {
        double stepNorm = (max - min) / 255;
        int[][] color = new int[function.length][function[0].length];

        BufferedImage image = new BufferedImage(function.length, function[0].length, BufferedImage.TYPE_INT_RGB);

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                for (int k = 1; k < 255; k++) {
                    if (function[i][j] > min + stepNorm * k) {
                        color[i][j] = k - 1;
                        image.setRGB(i, j, new Color(k - 1, k - 1, k - 1).getRGB());
                    }
                    if (function[i][j] == max) {
                        color[i][j] = 255;
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

        return color;
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

    public static double standardDeviation(int[][] reality, int[][] expect) {
        double deviation;

        int length = reality.length;
        double absReality, absExpect, firstMultiplier, secondMultiplier, firstSum = 0, secondSum = 0;

        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                absReality = Math.abs(reality[i][j]);
                absExpect = Math.abs(expect[i][j]);
                firstMultiplier = absReality * absReality - absExpect * absExpect;
                firstMultiplier *= firstMultiplier;
                firstSum += firstMultiplier;

                secondMultiplier = Math.pow(absExpect, 4);
                secondSum += secondMultiplier;
            }
        }

        double one = firstSum;
        double two  = Math.pow(secondSum, -1);
        deviation = Math.sqrt(one * two) * 100;

        return deviation;
    }
}