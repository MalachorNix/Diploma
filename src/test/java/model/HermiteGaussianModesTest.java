package model;

import org.junit.Assert;
import org.junit.Test;


public class HermiteGaussianModesTest {
    @Test
    public void hermiteGauss1D() throws Exception {
        Assert.assertEquals(1.06764987, HermiteGaussianModes.hermiteGauss1D(1, 2, 3), 0.01);
    }

    @Test
    public void hermiteGauss2D() throws Exception {
        Assert.assertEquals(0.407589, HermiteGaussianModes.hermiteGauss2D(1, 2, 3, 4, 5), 0.01);
    }

}