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

    private GapPenalty.Type type;
    private short gop, gep;

    /**
     * Creates a new set of gap penalties using the defaults.
     */
    public SimpleGapPenalty() {
        this(Default.getOpenPenalty(), Default.getExtensionPenalty());
    }

    /**
     * Creates a new set of gap penalties.
     *
     * @param gop the gap open penalty
     * @param gep the gap extension penalty
     */
    public SimpleGapPenalty(short gop, short gep) {
        this.gop = gop;
        this.gep = gep;
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
        this.gep = gep;
        setType();
    }

    @Override
    public void setOpenPenalty(short gop) {
        this.gop = gop;
        setType();
    }

    // helper method to set the type given the open and extension penalties
    private void setType() {
        type = (gep == 0) ? GapPenalty.Type.CONSTANT :
            (gop == gep) ? GapPenalty.Type.LINEAR :
                GapPenalty.Type.AFFINE;
    }

    /**
     * Stores the default values for the gap penalties.
     */
    public static class Default {

        private static GapPenalty instance = new SimpleGapPenalty((short) 10, (short) 1);

        /**
         * Returns the default gap extension penalty.
         *
         * @return the default gap extension penalty
         */
        public static short getExtensionPenalty() {
            return instance.getExtensionPenalty();
        }

        /**
         * Returns the default gap open penalty.
         *
         * @return the default gap open penalty
         */
        public static short getOpenPenalty() {
            return instance.getOpenPenalty();
        }

        /**
         * Sets a new default set of gap penalties.
         *
         * @param gop the default gap open penalty
         * @param gep the default gap extension penalty
         */
        public static void set(short gop, short gep) {
            instance = new SimpleGapPenalty(gop, gep);
        }

        /**
         * Sets the default gap extension penalty.
         *
         * @param gep the default gap extension penalty
         */
        public static void setExtensionPenalty(short gep) {
            instance.setExtensionPenalty(gep);
        }

        /**
         * Sets the default gap open penalty.
         *
         * @param gop the default gap open penalty
         */
        public static void setOpenPenalty(short gop) {
            instance.setOpenPenalty(gop);
        }

    }

}
