package model;

public final class HermiteGaussianModes {

    private HermiteGaussianModes() {

    }

    /**
     * Считает одномерную моду Гаусса-Эрмита.
     *
     * @param n     Индекс моды.
     * @param x     Аргумент моды по оси x.
     * @param gauss Гауссовый параметр.
     * @return Значение одномерной моды Гаусса-Эрмита.
     */
    public static double hermiteGauss1D(int n, double x, double gauss) {

        double polynomialArgument = x / gauss;

        return Math.exp(-(x * x) / (2 * gauss * gauss)) * HermitePolynomials.polynomial(n, polynomialArgument);
    }

    /**
     * Считает двумерную моду Гаусса-Эрмита.
     *
     * @param n     Первый индекс моды.
     * @param m     Второй индекс моды.
     * @param x     Аргумент моды по оси x.
     * @param y     Аргумент моды по оси y.
     * @param gauss Гауссовый параметр.
     * @return Значение двумерной моды Гаусса-Эрмита.
     */
    public static double hermiteGauss2D(int n, int m, double x, double y, double gauss) {

        double nPolynomialArg = x / gauss;
        double mPolynomialArg = y / gauss;

        return Math.exp(-(x * x + y * y) / (2 * gauss * gauss))
                * HermitePolynomials.polynomial(n, nPolynomialArg)
                * HermitePolynomials.polynomial(m, mPolynomialArg);
    }
}