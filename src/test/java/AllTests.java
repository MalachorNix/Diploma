import controller.FresnelTransform;
import model.GHMode;
import model.HermitePolynomials;
import org.apache.commons.math3.complex.Complex;
import org.junit.Test;

import static org.junit.Assert.*;


public class AllTests {
    @Test
    public void hermiteGauss1D() throws Exception {
        assertEquals(1.06764987, GHMode.hermiteGauss1D(1, 2, 3), 0.01);
    }

    @Test
    public void hermiteGauss2D() throws Exception {
        assertEquals(0.407589, GHMode.hermiteGauss2D(1, 2, 3, 4, 5), 0.01);
    }

    @Test
    public void polynomial() throws Exception {
        assertEquals(200416, HermitePolynomials.polynomial(10, 2), 0.01);
    }

    @Test
    public void transform1D() throws Exception {

    }

    @Test
    public void transform2D() throws Exception {

    }

    @Test
    public void transform2DSuperposition() throws Exception {

    }

    @Test
    public void firstMultiplier1DTransform() throws Exception {

    }

    @Test
    public void firstMultiplier2DTransform() throws Exception {

    }

    @Test
    public void firstExponent() throws Exception {

    }

    @Test
    public void integral1D() throws Exception {

    }

    /*@Test
    public void integral2D() throws Exception {
        Complex testValue = FresnelTransform.integral2D(1, 2, 3, 4, 5, 6, 7, 8, 0.001, 0.001, 4, 3, 11);

        assertEquals(-0.5380, testValue.getReal(), 0.01);
        assertEquals(1.4577, testValue.getImaginary(), 0.01);
    }*/

    @Test
    public void multiTest() throws Exception {
        Complex testValue = FresnelTransform.integrand2D(new Complex(23, 42),1, 2, 3, 4, 5, 6).multiply(0.01).multiply(0.01);

        assertEquals(-0.00477, testValue.getReal(), 0.01);
        assertEquals(0.00034, testValue.getImaginary(), 0.01);
    }

    @Test
    public void integral2DSuperposition() throws Exception {

    }

    @Test
    public void integrand1D() throws Exception {

    }

    @Test
    public void integrand2D() throws Exception {
        Complex function = new Complex(23, 42);
        Complex testValue = FresnelTransform.integrand2D(function, 1, 2, 3, 4, 5, 6);
        assertEquals(-47.7618691, testValue.getReal(), 0.01);
        assertEquals(3.43567368, testValue.getImaginary(), 0.01);
    }

    @Test
    public void transformExponent1D() throws Exception {

    }

    @Test
    public void transformExponent2D() throws Exception {
        Complex testValue = FresnelTransform.transformExponent2D(1, 2, 3, 4, 5, 6);
        assertEquals(-0.4161468, testValue.getReal(), 0.01);
        assertEquals(0.909297426, testValue.getImaginary(), 0.01);
    }

    @Test
    public void getArgumentTest() throws Exception {
        Complex testValue = new Complex(-721, -498);
        assertEquals(-2.537130631, testValue.getArgument(), 0.01);
    }

    @Test
    public void testExp() throws Exception {
        Complex test = Complex.I.multiply(3).exp();
        assertEquals(-0.989992496600, test.getReal(), 0.001);
        assertEquals(0.14112000805, test.getImaginary(), 0.001);
        assertEquals(3, test.getArgument(), 0.001);

        test = Complex.I.multiply(4).exp();
        assertEquals(-2.28318, test.getArgument(), 0.001);
        test = Complex.valueOf(0.5).exp();
        assertEquals(0, test.getArgument(), 0.001);


    }


}