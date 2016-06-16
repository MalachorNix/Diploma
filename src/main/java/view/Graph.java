package view;

import controller.FresnelTransform;
import org.apache.commons.math3.complex.Complex;
import javax.imageio.ImageIO;
import javax.imageio.stream.ImageOutputStream;
import java.awt.*;
import java.awt.image.*;
import java.io.ByteArrayOutputStream;
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

    public static void convertTo8BitImage() {
        try {
            BufferedImage source = ImageIO.read(new File("intensity.bmp"));
            BufferedImage gray = new BufferedImage(source.getWidth(), source.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
            Graphics2D g2d = gray.createGraphics();
            g2d.drawImage(source, 0, 0, null);
            g2d.dispose();

            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            ImageOutputStream ios = ImageIO.createImageOutputStream(baos);
            ImageIO.write(gray, "bmp", ios);
            ImageIO.write(gray, "bmp", new File("test.bmp"));
            ios.close();

        } catch (IOException ex) {
            ex.printStackTrace();
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

    /*public static void draw2DPhaseOnlyIntensity() {
        int width = 128, height = 128; // Размеры изображения.
        int N = 3; // Порядок полиномов Эрмита.
        int n = 3, m = 3;
        double u = -32, v = -32;
        double z = 100; // Расстояние.
        double k = 4; // Волновое число.
        double yMin = -32;
        double yMax = 32;
        double xMin = -32;
        double xMax = 32;
        double stepY = -2 * yMin / height;
        double stepX = -2 * xMin / width;
        double gauss = 7; // Гауссовый параметр.
        double uChange = u; // Переменные, которые будут меняться.
        double vChange = v;
        Complex[][] function = new Complex[width][height];
        double stepU = -2 * u / width;
        double stepV = -2 * v / height;
        Complex[][] coefficient = new Complex[N + 1][N + 1]; //Набор комплексных коэффициентов.

        for (int i = 0; i < coefficient.length; i++) {
            for (int j = 0; j < coefficient[0].length; j++) {
                coefficient[i][j] = new Complex(i, j);
            }
        }

        Complex[][] superposition = new Complex[width][height];

        double x = -32, y = -32;
        double xStep = -2 * x / width;
        double yStep = -2 * y / height;

        for (int i = 0; i < superposition.length; i++) {
            for (int j = 0; j < superposition[0].length; j++) {
                superposition[i][j] = superposition2D(N, x + i * xStep, y + j * yStep, gauss, coefficient);
            }
        }

        Complex[][] result = phaseOnlyEncode(superposition);

        for (int i = 0; i < function.length; i++, uChange += stepU) {
            for (int j = 0; j < function[0].length; j++, vChange += stepV) {
                function[i][j] = FresnelTransform.
                        transform2DPhaseOnly(uChange, vChange, z, k,
                                yMin, yMax,
                                xMin, xMax,
                                stepX, stepY,
                                N, gauss, width, height, result, coefficient);
            }
            vChange = v;
        }

        draw2DIntensity(function, "intensityPhaseOnly.bmp");
        draw2DPhase(function, "phasePhaseOnly.bmp");
    }*/

}