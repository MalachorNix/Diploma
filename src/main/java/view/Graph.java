package view;

import controller.FresnelTransform;
import controller.Superposition;
import org.apache.commons.math3.complex.Complex;
import org.junit.Test;

import javax.imageio.ImageIO;
import javax.imageio.stream.ImageOutputStream;
import java.awt.*;
import java.awt.color.ColorSpace;
import java.awt.image.*;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Graph {

    public void draw2DIntensity(Complex[][] function) {
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

        writeImage(intensity, maxAmplitude, minAmplitude, "intensity.bmp");

    }

    private void writeImage(double[][] function, double max, double min, String filename) {
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


    public void draw2DPhase(Complex[][] function) {
        double[][] phase = new double[function.length][function[0].length];
        double minPhase = Double.MAX_VALUE;
        double maxPhase = Double.MIN_VALUE;

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                phase[i][j] = function[i][j].getArgument();

                if (phase[i][j] > maxPhase) {
                    maxPhase = phase[i][j];
                }

                if (phase[i][j] < minPhase) {
                    minPhase = phase[i][j];
                }

            }
        }

        writeImage(phase, maxPhase, minPhase, "phase.bmp");
    }

    public void convert() {
        try {
            BufferedImage source = ImageIO.read(new File("intensityExit.bmp"));
            BufferedImage gray = new BufferedImage(source.getWidth(), source.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
            Graphics2D g2d = gray.createGraphics();
            g2d.drawImage(source, 0, 0, null);
            g2d.dispose();

            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            ImageOutputStream ios = ImageIO.createImageOutputStream(baos);
            ImageIO.write(gray, "png", ios);
            ImageIO.write(gray, "bmp", new File("test.bmp"));
            ios.close();

            byte[] array = baos.toByteArray();
        } catch (IOException ex) {
            Logger.getLogger(Test.class.getName()).log(Level.SEVERE, null, ex);
        }
    }


    public Complex[][] phaseOnlyEncode(Complex[][] superposition) {

        Complex[][] result = new Complex[superposition.length][superposition[0].length];

        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[0].length; j++) {
                result[i][j] = Complex.I.multiply(superposition[i][j].getArgument()).exp();
            }
        }


        return result;
    }

    public void elrngl() {
        int N = 10;
        Complex[][] coefficient = new Complex[N][N];

        for (int i = 0; i < coefficient.length; i++) {
            for (int j = 0; j < coefficient[0].length; j++) {
                coefficient[i][j] = new Complex(i, j);
            }
        }

        Complex[][] superposition = new Complex[128][128];

        double x = -64, y = -64;

        for (int i = 0; i < superposition.length; i++, x++) {
            for (int j = 0; j < superposition[0].length; j++, y++) {
                superposition[i][j] = Superposition.superposition2D(10, x, y, 11, coefficient);
            }
            y = -64;
        }

        Complex[][] result = phaseOnlyEncode(superposition);

        Complex[][] function = new Complex[128][128];

        double u = -64, v = -64;

        for (int i = 0; i < function.length; i++, u++) {
            for (int j = 0; j < function[0].length; j++, v++) {
                function[i][j] = FresnelTransform.transform2DSuperposition(u, v, 100, 4, -1, 1, -1, 1, 0.02, 0.02, result[i][j]);
            }
            v = -64;
        }

        draw2DIntensity(function);
    }

    public static void main(String[] args) {
        new Graph().elrngl();
    }
}