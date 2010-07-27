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
 * Created on June 9, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment;

import org.biojava3.alignment.template.GapPenalty;

/**
 * Implements a data structure for the gap penalties used during a sequence alignment routine.
 *
 * @author Mark Chapman
 */
public class SimpleGapPenalty implements GapPenalty {

    private static short dgop = 10, dgep = 1;

    /**
     * Sets the default gap extension penalty.
     *
     * @param gep the default gap extension penalty
     */
    public static void setDefaultExtensionPenalty(short gep) {
        dgep = gep;
    }

    /**
     * Sets the default gap open penalty.
     *
     * @param gop the default gap open penalty
     */
    public static void setDefaultOpenPenalty(short gop) {
        dgop = gop;
    }

    private GapPenalty.Type type;
    private short gop, gep;

    /**
     * Creates a new set of gap penalties using the defaults.
     */
    public SimpleGapPenalty() {
        this(dgop, dgep);
    }

    /**
     * Creates a new set of gap penalties.
     *
     * @param gop the gap open penalty
     * @param gep the gap extension penalty
     */
    public SimpleGapPenalty(short gop, short gep) {
        this.gop = (short) -Math.abs(gop);
        this.gep = (short) -Math.abs(gep);
        setType();
    }

    @Override
    public short getExtensionPenalty() {
        return gep;
    }

    @Override
    public short getOpenPenalty() {
        return gop;
    }

    @Override
    public Type getType() {
        return type;
    }

    @Override
    public void setExtensionPenalty(short gep) {
        this.gep = (short) -Math.abs(gep);
        setType();
    }

    @Override
    public void setOpenPenalty(short gop) {
        this.gop = (short) -Math.abs(gop);
        setType();
    }

    // helper method to set the type given the open and extension penalties
    private void setType() {
        type = (gop == 0) ? GapPenalty.Type.LINEAR : ((gep == 0) ? GapPenalty.Type.CONSTANT : GapPenalty.Type.AFFINE);
    }

}
