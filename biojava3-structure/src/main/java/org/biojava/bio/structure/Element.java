/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 29.04.2010
 *
 */
package org.biojava.bio.structure;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * Element is an enumeration of the elements of the periodic table. In addition,
 * several attributes of each element are accessible.
 * <B>Note:</B> Deuterium and Tritium are treated as separate elements D and T,
 * respectively. Sometimes part of a molecule is represented as an R-group, which
 * is represented as the element R.
 * 
 *  
 * @author Peter Rose
 * @version %I% %G%
 * @since 3.0
 * 
 */

public enum Element implements Serializable {

	// most frequently used elements first
    H(1, 1, 39, 1.10f, 0.32f, 1, 1, 1, 1, 1, 1.008f, 0, 1, 2.20f, ElementType.OTHER_NONMETAL),
    C(6, 2, 0, 1.55f, 0.77f, 4, 4, 4, 4, 4, 12.011f, 2, -4, 2.55f, ElementType.OTHER_NONMETAL),
    N(7, 2, 57, 1.40f, 0.75f, 5, 2, 5, 3, 4, 14.007f, 2, -3, 3.04f, ElementType.OTHER_NONMETAL),
    O(8, 2, 65, 1.35f, 0.73f, 6, 1, 2, 2, 2, 16.000f, 2, -2, 3.44f, ElementType.OTHER_NONMETAL),
    /**
     * Deuterim
     */
    D(1, 1, 27, 1.10f, 0.32f, 1, 1, 1, 1, 1, 1.008f, 0, 1, 2.20f, ElementType.OTHER_NONMETAL), // need to edit properties!
    /**
     * Tritium
     */
    T(1, 1, 90, 1.10f, 0.32f, 1, 1, 1, 1, 1, 1.008f, 0, 1, 2.20f, ElementType.OTHER_NONMETAL), // need to edit properties!
    He(2, 1, 40, 2.20f, 1.60f, 2, 0, 12, 0, 0, 4.003f, 2, 0, 0.0f, ElementType.NOBLE_GAS), // electroneg not reported
    Li(3, 2, 50, 1.22f, 1.34f, 1, 0, 12, 0, 1, 6.940f, 2, 1, 0.98f, ElementType.ALKALI_METAL),
    Be(4, 2, 12, 0.63f, 0.90f, 2, 0, 12, 2, 2, 9.012f, 2, 2, 1.57f, ElementType.TRANSITION_METAL),
    B(5, 2, 10, 1.55f, 0.82f, 3, 3, 5, 3, 4, 10.810f, 2, 3, 2.04f, ElementType.METALLOID),
    F(9, 2, 32, 1.30f, 0.72f, 7, 0, 1, 1, 1, 18.998f, 2, -1, 3.98f, ElementType.HALOGEN),
    Ne(10, 2, 61, 2.02f, 1.12f, 8, 0, 12, 0, 0, 20.170f, 10, 0, 0.00f, ElementType.NOBLE_GAS), // electroneg not reported
    Na(11, 3, 58, 2.20f, 1.54f, 1, 0, 1, 0, 0, 22.990f, 10, 1, 0.93f, ElementType.ALKALI_METAL),
    Mg(12, 3, 54, 1.50f, 1.30f, 2, 0, 2, 0, 2, 24.305f, 10, 2, 1.31f, ElementType.ALKALINE_EARTH_METAL),
    Al(13, 3, 4, 1.50f, 1.18f, 3, 0, 5, 0, 4, 26.982f, 10, 3, 1.61f, ElementType.POST_TRANSITION_METAL),
    Si(14, 3, 86, 2.20f, 1.11f, 4, 4, 4, 4, 4, 28.086f, 10, 4, 1.90f, ElementType.METALLOID),
    P(15, 3, 67, 1.88f, 1.06f, 5, 3, 5, 3, 5, 30.974f, 10, 5, 2.19f, ElementType.OTHER_NONMETAL),
    S(16, 3, 82, 1.81f, 1.02f, 6, 2, 6, 2, 6, 32.060f, 10, -2, 2.58f, ElementType.OTHER_NONMETAL), 
    Cl(17, 3, 21, 1.75f, 0.99f, 7, 0, 1, 1, 1, 35.453f, 10, -1, 3.16f, ElementType.HALOGEN),
    Ar(18, 4, 6, 2.77f, 1.54f, 8, 0, 12, 0, 0, 39.948f, 18, 0, 0.00f, ElementType.NOBLE_GAS), // electroneg not reported
    K(19, 4, 47, 2.39f, 1.96f, 1, 0, 12, 0, 0, 39.102f, 18, 1, 0.82f, ElementType.ALKALI_METAL),
    Ca(20, 4, 17, 1.95f, 1.74f, 2, 0, 2, 0, 0, 40.080f, 18, 2, 1.00f, ElementType.ALKALINE_EARTH_METAL),
    Sc(21, 4, 84, 1.32f, 1.44f, 3, 0, 12, 3, 0, 44.956f, 18, 3, 1.36f, ElementType.TRANSITION_METAL),
    Ti(22, 4, 96, 1.95f, 1.36f, 4, 2, 4, 3, 4, 47.880f, 18, 4, 1.54f, ElementType.TRANSITION_METAL),
    V(23, 4, 100, 1.06f, 1.25f, 5, 0, 12, 3, 0, 50.040f, 18, 5, 1.63f, ElementType.TRANSITION_METAL),
    Cr(24, 4, 24, 1.13f, 1.27f, 6, 0, 12, 2, 0, 51.996f, 18, 3, 1.66f, ElementType.TRANSITION_METAL),
    Mn(25, 4, 55, 1.19f, 1.39f, 7, 0, 12, 0, 0, 54.938f, 18, 2, 1.55f, ElementType.TRANSITION_METAL),
    Fe(26, 4, 33, 1.95f, 1.25f, 3, 0, 8, 0, 0, 55.847f, 18, 3, 1.83f, ElementType.TRANSITION_METAL),
    Co(27, 4, 23, 1.13f, 1.26f, 3, 0, 12, 0, 0, 58.933f, 18, 2, 1.88f, ElementType.TRANSITION_METAL),
    Ni(28, 4, 62, 1.24f, 1.21f, 3, 0, 12, 0, 0, 58.710f, 18, 2, 1.91f, ElementType.TRANSITION_METAL),
    Cu(29, 4, 26, 1.15f, 1.38f, 2, 0, 4, 0, 0, 63.546f, 18, 2, 1.90f, ElementType.TRANSITION_METAL),
    Zn(30, 4, 106, 1.15f, 1.31f, 2, 0, 2, 0, 0, 65.380f, 18, 2, 1.65f, ElementType.TRANSITION_METAL),
    Ga(31, 4, 36, 1.55f, 1.26f, 3, 1, 4, 2, 4, 69.720f, 28, 3, 1.81f, ElementType.POST_TRANSITION_METAL),
    Ge(32, 4, 38, 2.72f, 1.22f, 4, 0, 12, 4, 4, 72.590f, 28, 4, 2.01f, ElementType.METALLOID),
    As(33, 4, 7, 0.83f, 1.19f, 5, 0, 12, 3, 5, 74.922f, 28, -3, 2.18f, ElementType.METALLOID),
    Se(34, 4, 85, 0.90f, 1.16f, 6, 0, 12, 2, 6, 78.960f, 28, 4, 2.55f, ElementType.OTHER_NONMETAL),
    Br(35, 4, 15, 1.95f, 1.14f, 7, 0, 1, 1, 1, 79.904f, 28, -1, 2.96f, ElementType.HALOGEN),
    Kr(36, 4, 48, 1.90f, 1.60f, 8, 0, 12, 0, 0, 83.800f, 28, 0, 3.00f, ElementType.NOBLE_GAS),
    Rb(37, 5, 77, 2.65f, 2.11f, 1, 0, 12, 0, 0, 85.467f, 36, 1, 0.82f, ElementType.ALKALI_METAL),
    Sr(38, 5, 89, 2.02f, 1.92f, 2, 0, 12, 2, 0, 87.620f, 36, 2, 0.95f, ElementType.ALKALINE_EARTH_METAL),
    Y(39, 5, 103, 1.61f, 1.62f, 3, 0, 12, 3, 0, 88.806f, 36, 3, 1.22f, ElementType.TRANSITION_METAL),
    Zr(40, 5, 105, 1.42f, 1.48f, 4, 0, 12, 4, 0, 91.220f, 36, 4, 1.33f, ElementType.TRANSITION_METAL),
    Nb(41, 5, 59, 1.33f, 1.37f, 5, 0, 12, 3, 0, 92.906f, 36, 5, 1.60f, ElementType.TRANSITION_METAL),
    Mo(42, 5, 56, 1.75f, 1.45f, 6, 1, 6, 3, 0, 95.940f, 36, 6, 2.16f, ElementType.TRANSITION_METAL),
    Tc(43, 5, 93, 1.80f, 1.56f, 7, 0, 12, 6, 0, 98.910f, 36, 7, 1.90f, ElementType.TRANSITION_METAL),
    Ru(44, 5, 81, 1.20f, 1.26f, 8, 0, 12, 3, 0, 101.070f, 36, 4, 2.20f, ElementType.TRANSITION_METAL),
    Rh(45, 5, 79, 1.22f, 1.35f, 4, 0, 12, 3, 0, 102.906f, 36, 3, 2.28f, ElementType.TRANSITION_METAL),
    Pd(46, 5, 70, 1.44f, 1.31f, 4, 0, 12, 2, 0, 106.400f, 36, 2, 2.20f, ElementType.TRANSITION_METAL),
    Ag(47, 5, 3, 1.55f, 1.53f, 1, 0, 6, 0, 0, 107.868f, 36, 1, 1.93f, ElementType.TRANSITION_METAL),
    Cd(48, 5, 18, 1.75f, 1.48f, 2, 0, 12, 0, 0, 112.400f, 36, 2, 1.69f, ElementType.TRANSITION_METAL),
    In(49, 5, 45, 1.46f, 1.44f, 3, 0, 12, 3, 0, 114.820f, 46, 3, 1.78f, ElementType.POST_TRANSITION_METAL),
    Sn(50, 5, 88, 1.67f, 1.41f, 4, 0, 12, 2, 4, 118.690f, 46, 4, 1.96f, ElementType.POST_TRANSITION_METAL),
    Sb(51, 5, 83, 1.12f, 1.38f, 5, 0, 12, 4, 5, 121.750f, 46, -3, 2.05f, ElementType.METALLOID),
    Te(52, 5, 94, 1.26f, 1.35f, 6, 0, 12, 2, 6, 127.600f, 46, 4, 2.10f, ElementType.METALLOID),
    I(53, 5, 44, 2.15f, 1.33f, 7, 1, 1, 1, 1, 126.905f, 46, -1, 2.66f, ElementType.HALOGEN),
    Xe(54, 5, 102, 2.10f, 1.70f, 8, 0, 12, 0, 0, 131.300f, 46, 0, 2.60f, ElementType.NOBLE_GAS),
    Cs(55, 6, 25, 3.01f, 2.25f, 1, 0, 12, 0, 0, 132.905f, 54, 1, 0.79f, ElementType.ALKALI_METAL),
    Ba(56, 6, 11, 2.41f, 1.98f, 2, 0, 12, 0, 0, 137.340f, 54, 2, 0.89f, ElementType.ALKALINE_EARTH_METAL),
    La(57, 6, 49, 1.83f, 1.95f, 3, 0, 12, 3, 0, 138.905f, 54, 3, 1.10f, ElementType.LANTHANOID),
    Ce(58, 6, 19, 1.86f, 1.03f, 4, 0, 12, 3, 0, 140.120f, 54, 3, 1.12f, ElementType.LANTHANOID),
    Pr(59, 6, 73, 1.62f, 0.90f, 4, 0, 12, 3, 0, 140.908f, 55, 3, 1.13f, ElementType.LANTHANOID),
    Nd(60, 6, 60, 1.79f, 0.99f, 3, 0, 12, 3, 0, 144.240f, 56, 3, 1.14f, ElementType.LANTHANOID),
    Pm(61, 6, 71, 1.76f, 0.98f, 3, 0, 12, 3, 0, 145.000f, 58, 3, 1.13f, ElementType.LANTHANOID),
    Sm(62, 6, 87, 1.74f, 0.96f, 3, 0, 12, 2, 0, 150.400f, 59, 3, 1.17f, ElementType.LANTHANOID),
    Eu(63, 6, 31, 1.96f, 1.09f, 3, 0, 12, 2, 0, 151.960f, 60, 3, 1.20f, ElementType.LANTHANOID),
    Gd(64, 6, 37, 1.69f, 0.94f, 3, 0, 12, 3, 0, 157.250f, 61, 3, 1.20f, ElementType.LANTHANOID),
    Tb(65, 6, 92, 1.66f, 0.92f, 4, 0, 12, 3, 0, 158.925f, 61, 3, 1.10f, ElementType.LANTHANOID),
    Dy(66, 6, 28, 1.63f, 0.91f, 3, 0, 12, 3, 0, 162.500f, 62, 3, 1.22f, ElementType.LANTHANOID),
    Ho(67, 6, 43, 1.61f, 0.89f, 3, 0, 12, 3, 0, 164.930f, 64, 3, 1.23f, ElementType.LANTHANOID),
    Er(68, 6, 29, 1.59f, 0.88f, 3, 0, 12, 3, 0, 167.260f, 65, 3, 1.24f, ElementType.LANTHANOID),
    Tm(69, 6, 98, 1.57f, 0.87f, 3, 0, 12, 3, 0, 168.934f, 66, 3, 1.25f, ElementType.LANTHANOID),
    Yb(70, 6, 104, 1.54f, 0.86f, 3, 0, 12, 2, 0, 173.040f, 67, 3, 1.10f, ElementType.LANTHANOID),
    Lu(71, 6, 52, 1.53f, 0.85f, 3, 0, 12, 3, 0, 174.970f, 68, 3, 1.27f, ElementType.LANTHANOID),
    Hf(72, 6, 41, 1.40f, 1.58f, 4, 0, 12, 4, 0, 178.490f, 68, 4, 1.30f, ElementType.TRANSITION_METAL),
    Ta(73, 6, 91, 1.22f, 1.38f, 5, 0, 12, 5, 0, 180.850f, 68, 5, 1.50f, ElementType.TRANSITION_METAL),
    W(74, 6, 101, 1.26f, 1.46f, 6, 0, 12, 6, 0, 183.850f, 68, 6, 2.36f, ElementType.TRANSITION_METAL),
    Re(75, 6, 78, 1.30f, 1.59f, 7, 0, 12, 4, 0, 186.200f, 68, 7, 1.90f, ElementType.TRANSITION_METAL),
    Os(76, 6, 66, 1.58f, 1.28f, 8, 0, 12, 2, 0, 190.200f, 68, 4, 2.20f, ElementType.TRANSITION_METAL),
    Ir(77, 6, 46, 1.22f, 1.37f, 6, 0, 12, 3, 0, 192.220f, 68, 4, 2.20f, ElementType.TRANSITION_METAL),
    Pt(78, 6, 74, 1.55f, 1.28f, 4, 0, 6, 0, 0, 195.090f, 68, 4, 2.28f, ElementType.TRANSITION_METAL),
    Au(79, 6, 9, 1.45f, 1.44f, 3, 0, 6, 0, 0, 196.967f, 68, 3, 2.54f, ElementType.TRANSITION_METAL),
    Hg(80, 6, 42, 1.55f, 1.32f, 2, 0, 12, 1, 2, 200.59f, 78, 1, 2.00f, ElementType.TRANSITION_METAL),
    Tl(81, 6, 97, 1.96f, 1.45f, 3, 0, 12, 1, 3, 204.3833f, 78, 1, 1.62f, ElementType.POST_TRANSITION_METAL),
    Pb(82, 6, 69, 2.16f, 1.47f, 4, 0, 12, 2, 4, 207.200f, 78, 2, 2.33f, ElementType.POST_TRANSITION_METAL),
    Bi(83, 6, 13, 1.73f, 1.46f, 5, 0, 12, 3, 3, 208.981f, 78, 3, 2.20f, ElementType.POST_TRANSITION_METAL),
    Po(84, 6, 72, 1.21f, 0.67f, 6, 0, 12, 4, 2, 209.000f, 78, 4, 2.0f, ElementType.METALLOID),
    At(85, 6, 8, 1.12f, 0.62f, 7, 0, 12, 1, 1, 210.000f, 78, -1, 2.20f, ElementType.HALOGEN),
    Rn(86, 6, 80, 2.30f, 1.90f, 8, 0, 12, 0, 0, 222.000f, 78, 0, 0.0f, ElementType.NOBLE_GAS), // electroneg not reported
    Fr(87, 7, 35, 3.24f, 1.80f, 1, 0, 12, 0, 0, 223.000f, -1, 1, 0.70f, ElementType.ALKALI_METAL),
    Ra(88, 7, 76, 2.57f, 1.43f, 2, 0, 12, 2, 0, 226.000f, -1, 2, 0.9f, ElementType.ALKALINE_EARTH_METAL),
    Ac(89, 7, 2, 2.12f, 1.18f, 3, 0, 12, 4, 0, 227.000f, -1, 3, 1.1f, ElementType.ACTINOID),
    Th(90, 7, 95, 1.84f, 1.02f, 4, 0, 12, 1, 0, 232.038f, -1, 4, 1.30f, ElementType.ACTINOID),
    Pa(91, 7, 68, 1.60f, 0.89f, 5, 0, 12, 4, 0, 231.036f, -1, 5, 1.50f, ElementType.ACTINOID),
    U(92, 7, 99, 1.75f, 0.97f, 6, 0, 12, 4, 0, 238.029f, -1, 6, 1.38f, ElementType.ACTINOID),
    Np(93, 7, 64, 1.71f, 0.95f, 6, 0, 12, 4, 0, 237.048f, -1, 5, 1.36f, ElementType.ACTINOID),
    Pu(94, 7, 75, 1.67f, 0.93f, 6, 0, 12, 3, 0, 244.000f, -1, 4, 1.28f, ElementType.ACTINOID),
    Am(95, 7, 5, 1.66f, 0.92f, 6, 0, 12, 3, 0, 243.000f, -1, 3, 1.13f, ElementType.ACTINOID),
    Cm(96, 7, 22, 1.65f, 0.91f, 3, 0, 12, 3, 0, 248.000f, -1, 3, 1.28f, ElementType.ACTINOID),
    Bk(97, 7, 14, 1.64f, 0.90f, 4, 0, 12, 3, 0, 247.000f, -1, 3, 1.30f, ElementType.ACTINOID),
    Cf(98, 7, 20, 1.63f, 0.89f, 3, 0, 12, 4, 0, 251.000f, -1, 3, 1.30f, ElementType.ACTINOID),
    Es(99, 7, 30, 1.62f, 0.88f, -1, 0, 12, 4, 0, 254.000f, -1, 3, 1.30f, ElementType.ACTINOID),
    Fm(100, 7, 34, 1.61f, 0.87f, -1, 0, 12, 4, 0, 257.000f, -1, 3, 1.30f, ElementType.ACTINOID),
    Md(101, 7, 53, 1.60f, 0.86f, -1, 0, 12, 4, 0, 256.000f, -1, 3, 1.30f, ElementType.ACTINOID),
    No(102, 7, 63, 1.59f, 0.85f, -1, 0, 12, 4, 0, 254.000f, -1, 3, 1.30f, ElementType.ACTINOID),
    Lr(103, 7, 51, 1.58f, 0.84f, -1, 0, 12, 4, 0, 257.000f, -1, 3, 0.00f, ElementType.ACTINOID), // electroneg not reported
    /**
     * R-group to represent generic groups that are sometimes present in MDL .sdf
     * files.
     */
    R(104, 0, 105, 0.0f, 0.0f, 0, 0, 4, 1, 0, 0.000f, -1, 3, 0.00f, ElementType.UNKNOWN); // this is an R-group
    // should these be declared final?
    private int atomicNumber;
    private int period;
    //private int hillOrder;
    private float VDWRadius; // in Angstroms
    private float covalentRadius; // in Angstroms
    private int valenceElectronCount;
    private int minimumValence;
    private int maximumValence;
    private int commonValence;
    private int maximumCovalentValence;
    private float atomicMass;
    private int coreElectronCount;
    private int oxidationState;
    // Pauling electronegativity: http://en.wikipedia.org/wiki/Electronegativity
    private float paulingElectronegativity;
    // Element type: http://www.ptable.com/
    private ElementType elementType;
    //private static final Element[] hillOrderIndex;
    
//
//    static {
//        hillOrderIndex = new Element[Element.values().length + 1];
//        for (Element e : Element.values()) {
//            hillOrderIndex[e.getHillOrder()] = e;
//        }
//        hillOrderIndex[Element.H.getHillOrder()] = Element.H; // special case for hydrogen
//    }

    private static final Map<String,Element> allElements ;
   
    static {
    	allElements = new HashMap<String,Element>();
    	for (Element e : Element.values()){
    		allElements.put(e.toString().toLowerCase(), e);
    	}
    }
    private Element(int atomicNumber,
            int period,
            int hillOrder,
            float VDWRadius,
            float covalentRadius,
            int valenceElectronCount,
            int minimumValence,
            int maximumValence,
            int commonValence,
            int maximumCovalentValence,
            float atomicMass,
            int coreElectronCount,
            int oxidationState,
            float paulingElectronegativity,
            ElementType elementType) {

        this.atomicNumber = atomicNumber;
        this.period = period;
        //this.hillOrder = hillOrder;
        this.VDWRadius = VDWRadius;
        this.covalentRadius = covalentRadius;
        this.valenceElectronCount = valenceElectronCount;
        this.minimumValence = minimumValence;
        this.maximumValence = maximumValence;
        this.commonValence = commonValence;
        this.maximumCovalentValence = maximumCovalentValence;
        this.atomicMass = atomicMass;
        this.coreElectronCount = coreElectronCount;
        this.oxidationState = oxidationState;
        this.paulingElectronegativity = paulingElectronegativity;
        this.elementType = elementType;
        
       
    }
    /**
     * Returns the atomic number of this Element.
     * @return the atomic number of this Element.
     */
    public int getAtomicNumber() {
        return atomicNumber;
    }

    /**
     * Returns the period in the periodic table of this Element.
     * @return the period in the periodic table of this Element.
     */
    public int getPeriod() {
        return period;
    }

    /**
     * Returns the Hill Order of this Element. The Hill Order represents the
     * priority by which elements are sorted in molecular formulas.
     * The Hill system is a system of writing chemical formulas such that the
     * number of carbon atoms in a molecule is indicated first, the number of
     * hydrogen atoms next, and then the number of all other chemical elements
     * subsequently, in alphabetical order. When the formula contains no carbon,
     * all the elements, including hydrogen, are listed alphabetically.
     * <p>
     * Edwin A. Hill, "On A System Of Indexing Chemical Literature;
     * Adopted By The Classification Division Of The U. S. Patent Office".
     * J. Am. Chem. Soc. 1900, 22(8), 478-494.
     * <p>
     * <a href="http://en.wikipedia.org/wiki/Hill_system">
     * http://en.wikipedia.org/wiki/Hill_system</a>
     * <p>
     * @return the Hill Order of this Element.
     */
    public int getHillOrder() {
    	throw new RuntimeException("Not implemented, yet!");
    	//throw new NotImplementedYetException();
        //return hillOrder;
    }

    /**
     * Returns the van der Waals radius of this Element.
     * @return the van der Waals radius of this Element, measured in Angstroms.
     */
    public float getVDWRadius() {
        return VDWRadius;
    }

    /**
     * Returns the covalent radius of this Element.
     * @return covalent radius, measured in Angstroms.
     */
    public float getCovalentRadius() {
        return covalentRadius;
    }

    /**
     * Returns the number of valence electrons for this Element.
     * @return the number of valence electrons for this Element.
     */
    public int getValenceElectronCount() {
        return valenceElectronCount;
    }

    /**
     * Returns the minimum valence for this Element.
     * @return the minimum valence of this atom.
     */
    public int getMinimumValence() {
        return minimumValence;
    }

    /**
     * Returns the maximum valence for this Element.
     * @return the maximum valence for this Element.
     */
    public int getMaximumValence() {
        return maximumValence;
    }

    /**
     * Returns the common valence for this Element.
     * @return the common valence for this Element.
     */
    public int getCommonValence() {
        return commonValence;
    }

    /**
     * Returns the maximum valence for this Element.
     * @return the maximum valence of this element.
     */
    public int getMaximumCovalentValence() {
        return maximumCovalentValence;
    }

    /**
     * Returns the atomic mass for this Element.
     * @return the atomic mass for this Element, measured in g/mol.
     */
    public float getAtomicMass() {
        return atomicMass;
    }

    /**
     * Returns the number of core electrons for this Element.
     * @return number of core electrons for this Element.
     */
    public int getCoreElectronCount() {
        return coreElectronCount;
    }

    /**
     * Returns a typical oxidation state for this Element. This information is mostly
     * useful for metals.
     * @return a typical oxidation state for this Element.
     */
    public int getOxidationState() {
        return oxidationState;
    }

    /**
     * Returns the Pauling electronegativity for this Element.
     * @return the Pauling electronegativity for this Element.
     */
    public float getPaulingElectronegativity() {
        return paulingElectronegativity;
    }

    /**
     * Returns the Element Type for this Element.
     * @return the Element Type for this Element.
     */
    public ElementType getElementType() {
        return elementType;
    }
    
    /**
     * Returns the Element that corresponds to the specified element symbol. The case
     * of the element symbol is ignored. Example: FE, fe, Fe represent iron.
     * @param elementSymbol element symbol to specify Element.
     * @return the Element specified by the element symbol.
     */
    public static Element valueOfIgnoreCase(String elementSymbol) throws IllegalArgumentException {

    	Element e = allElements.get(elementSymbol.toLowerCase());
        if ( e != null)
        	return e;
        throw new IllegalArgumentException("Invalid element symbol: " + elementSymbol);
    }

    /**
     * Returns true if this Element is Hydrogen. <B>Note:</B> currently Deuterium (D)
     * and Tritium (T) are not considered Hydrogen.
     * @return <CODE>true</CODE> if the Element is Hydrogen.
     */
    public boolean isHydrogen() {
        return (this == H);
    }

    /**
     * Returns <CODE>true</CODE> is the Element is an not Hydrogen.
     * @return <CODE>true</CODE> is Element is not Hydrogen.
     */
    public boolean isHeavyAtom() {
        return (this != H);
    }

    /**
     * Returns <CODE>true</CODE> if Element is not Hydrogen and not Carbon.
     * @return <CODE>true</CODE> if Element is not Hydrogen and not Carbon.
     */
    public boolean isHeteroAtom() {
        return !(this == C || this == H);
    }

	/**
     * Returns <CODE>true</CODE> if ElementType is a metal.
     * @return <CODE>true</CODE> if ElementType is a metal.
     */
	public boolean isMetal() {
		return elementType.isMetal();
	}
	
	/**
     * Returns <CODE>true</CODE> if ElementType is a metalloid.
     * @return <CODE>true</CODE> if ElementType is a metalloid.
     */
	public boolean isMetalloid() {
		return elementType.isMetalloid();
	}
	
	/**
     * Returns <CODE>true</CODE> if ElementType is a non-metal.
     * @return <CODE>true</CODE> if ElementType is a non-metal.
     */
	public boolean isNonMetal() {
		return elementType.isNonMetal();
	}
	
	/**
     * Returns <CODE>true</CODE> if Element is a halogen (F, Cl, Br, I, At).
     * @return <CODE>true</CODE> if Element is a halogen.
     */
    public boolean isHalogen() {
        return elementType.equals(ElementType.HALOGEN);
    }
    
    /**
     * Returns <CODE>true</CODE> if Element is a chalcogen (O, S, Se, Te, Po).
     * @return <CODE>true</CODE> if Element is a chalcogen.
     */
    public boolean isChalcogen() {
        return (this == O || this == S || this == Se || this == Te ||
                this == Po);
    }

    /**
     * Returns the Element that corresponds to the specified Hill Order.
     * @param index the Hill Order.
     * @return the Element that corresponds to the specified Hill Order.
     * @see #getHillOrder()
     */
    public static Element getElementFromHillIndex(int index) {
    	throw new RuntimeException("Not implemented, yet!");
        //return hillOrderIndex[index];
    }
}
