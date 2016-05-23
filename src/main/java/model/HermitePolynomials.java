package model;

public final class HermitePolynomials {

    private HermitePolynomials() {

    }

    /**
     * Считает n-полином Эрмита в физическом определении.
     *
     * @param n Номер полинома Эрмита.
     * @param x Аргумент полинома Эрмита.
     * @return Значение n-полинома Эрмита.
     */
    public static double polynomial(int n, double x) {

        if (n == 0) {
            return 1;
        }

        if (n == 1) {
            return 2 * x;
        }

        double result = 0;

        for (int i = 2; i <= n; i++) {
            // result = 2 * x * polynomial(i - 1, x) - 2 * (n - 1) * polynomial(i - 2, x);
            result = 2 * x * polynomial(i - 1, x) - 2 * (i - 1) * polynomial(i - 2, x);
        }

        return result;
    }
}