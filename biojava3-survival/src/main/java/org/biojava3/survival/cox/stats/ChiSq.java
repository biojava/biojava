/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.cox.stats;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class ChiSq {

    /**
     *
     * @param z
     * @return
     */
    public static double norm(double z) {
        return ChiSq.chiSq(z * z, 1);
    }
    
    /**
     *
     * @param x
     * @param n
     * @return
     */
    public static double chiSq(double x, int n) {
        double p = Math.exp(-0.5 * x);
        if ((n % 2) == 1) {
            p = p * Math.sqrt(2 * x / Math.PI);
        }
        double k = n;
        while (k >= 2) {
            p = p * x / k;
            k = k - 2;
        }
        double t = p;
        double a = n;
        while (t > 0.000001 * p) {
            a = a + 2;
            t = t * x / a;
            p = p + t;
        }
        return 1 - p;
    }

    
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
    }
}
