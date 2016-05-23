import controller.FresnelTransform;
import org.apache.commons.math3.complex.Complex;
import view.Graph;
import static controller.Superposition.*;

public class Main {

    public static void main(String[] args) {

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
        Complex[][] function = new Complex[width][height];
        double stepU = -2 * u / width;
        double stepV = -2 * v / height;
        double x = -32, y = -32;
        double xStep = -2 * x / width;
        double yStep = -2 * y / height;
        Complex[][] coefficient = new Complex[N + 1][N + 1]; //Набор комплексных коэффициентов.

        // Задаем им значения.
        for (int i = 0; i < coefficient.length; i++) {
            for (int j = 0; j < coefficient[0].length; j++) {
                coefficient[i][j] = new Complex(i, j);
            }
        }

        // Формула 5.
        transform2D(width, height, n, m, z, k, yMin, yMax, xMin, xMax, stepY, stepX, gauss, u, v, function, stepU, stepV);
        // Формула 8.
        transform2DSuperposition(width, height, N, z, k, yMin, yMax, xMin, xMax, stepY, stepX, gauss, u, v, function, stepU, stepV, coefficient);
        // Формула 7.
        superposition(width, height, N, gauss, coefficient, x, y, xStep, yStep);
    }

    private static void transform2D(int width, int height, int n, int m, double z, double k, double yMin, double yMax, double xMin, double xMax, double stepY, double stepX, double gauss, double u, double v, Complex[][] function, double stepU, double stepV) {

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                function[i][j] = FresnelTransform.transform2D(u + i * stepU , v + j * stepV, z, k,
                        yMin, yMax,
                        xMin, xMax,
                        stepX, stepY,
                        n, m, gauss, width, height);
            }
        }

        Graph.draw2DIntensity(function, "intensity5.bmp");
        Graph.draw2DPhase(function, "phase5.bmp");
    }

    private static void superposition(int width, int height, int N, double gauss, Complex[][] coefficient, double x, double y, double xStep, double yStep) {
        Complex[][] superposition = new Complex[width][height];

        for (int i = 0; i < superposition.length; i++) {
            for (int j = 0; j < superposition[0].length; j++) {
                superposition[i][j] = superposition2D(N, x + i * xStep, y + j * yStep, gauss, coefficient);
            }
        }

        Graph.draw2DIntensity(superposition, "intensity7.bmp");
        Graph.draw2DPhase(superposition, "phase7.bmp");
    }

    private static void transform2DSuperposition(int width, int height, int N, double z, double k, double yMin, double yMax, double xMin, double xMax, double stepY, double stepX, double gauss, double u, double v, Complex[][] function, double stepU, double stepV, Complex[][] coefficient) {

        for (int i = 0; i < function.length; i++) {
            for (int j = 0; j < function[0].length; j++) {
                function[i][j] = FresnelTransform.transform2DSuperposition(u + i * stepU, v + j * stepV, z, k,
                                yMin, yMax,
                                xMin, xMax,
                                stepX, stepY,
                                N, gauss, coefficient, width, height);
            }
        }

        Graph.draw2DIntensity(function, "intensity8.bmp");
        Graph.draw2DPhase(function, "phase8.bmp");
    }
}