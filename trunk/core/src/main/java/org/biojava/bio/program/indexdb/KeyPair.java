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

package org.biojava.bio.program.indexdb;

/**
 * <code>KeyPair</code> represents a mapping of a primary identifier
 * to a secondary identifier for the same <code>Record</code>.
 *
 * @author Matthew Pocock
 * @author Keith James
 */
interface KeyPair {

    /**
     * <code>getPrimary</code> returns the primary identifier of the
     * pair.
     *
     * @return a <code>String</code> primary ID.
     */
    public String getPrimary();

    /**
     * <code>getSecondary</code> returns the secondary identifier of
     * the pair.
     *
     * @return a <code>String</code> secondary ID.
     */
    public String getSecondary();

    /**
     * <code>Impl</code> is the default implementation.
     */
    class Impl implements KeyPair {
        private final String primary;
        private final String secondary;

        public Impl(String primary, String secondary) {
            this.primary = primary;
            this.secondary = secondary;
        }

        public String getPrimary() {
            return primary;
        }

        public String getSecondary() {
            return secondary;
        }

        public String toString() {
            return primary + ":" + secondary;
        }
    }
}
