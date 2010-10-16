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
 */

package org.biojava.bio.chromatogram.graphic;

import org.biojava.bio.chromatogram.Chromatogram;

/**
 * Provides the mechanism whereby a ChromatogramGraphic can display
 * a Chromatogram with a non-linear horizontal scale.
 *
 * @author Rhett Sutphin (<a href="http://genome.uiowa.edu/">UI CBCB</a>)
 * @author Matthew Pocock
 * @since 1.3
 */
public interface ChromatogramNonlinearScaler {
    /**
     * Returns the remapped coordinate for the provided trace sample index of the
     * given chromatogram.
     *
     * @param c  the Chromatogram
     * @param sampleIndex the sample index
     * @return the new coordinagte
     */
    public float scale(Chromatogram c, int sampleIndex) 
        throws IndexOutOfBoundsException;

    /**
     * The default scaler that displays the chromatogram 1:1.
     */
    public static class Identity implements ChromatogramNonlinearScaler {
        private static Identity INSTANCE;

      /**
       * Retrieve the singleton instance of this class.
       *
       * @return the Identity instance
       */
        public static Identity getInstance() {
            if (INSTANCE == null)
                INSTANCE = new Identity();
            return INSTANCE;
        }

        private Identity() { }

        public float scale(Chromatogram c, int sampleIndex) {
            return sampleIndex;
        }
    }
}
