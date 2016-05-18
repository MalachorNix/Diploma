package model;

import org.junit.Test;

import static org.junit.Assert.*;

public class HermitePolynomialsTest {
    @Test
    public void polynomial() throws Exception {
        assertEquals(200416, HermitePolynomials.polynomial(10, 2), 0.01);
    }

}