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
 * Created on June 7, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment.template;

/**
 * Defines a data structure for the gap penalties used during a sequence alignment routine.
 *
 * @author Mark Chapman
 */
public interface GapPenalty {

    /**
     * Defines the possible types of gap penalties.  This is:
     * <ul>
     *  <li>CONSTANT, if static and the extension penalty is 0
     *  <li>LINEAR, if static and the open penalty is 0
     *  <li>AFFINE, if static but neither CONSTANT nor LINEAR
     *  <li>DYNAMIC, if penalty values change during alignment
     * </ul>
     */
    enum Type {CONSTANT, LINEAR, AFFINE, DYNAMIC};

    /**
     * Returns penalty given when an already open gap elongates by a single element
     *
     * @return gap extension penalty
     */
    short getExtensionPenalty();

    /**
     * Returns penalty given when a deletion or insertion gap first opens
     *
     * @return gap open penalty
     */
    short getOpenPenalty();

    /**
     * Returns {@link GapPenalty.Type} stored.
     *
     * @return gap penalty type
     */
    Type getType();

    /**
     * Sets penalty given when an already open gap elongates by a single element
     *
     * @param gep gap extension penalty
     */
    void setExtensionPenalty(short gep);

    /**
     * Sets penalty given when a deletion or insertion gap first opens
     *
     * @param gop gap open penalty
     */
    void setOpenPenalty(short gop);

}
