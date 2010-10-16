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

package org.biojavax.bio.seq;

/**
 * Resolves a position that is fuzzy or covers a range of bases by
 * converting it to a single base.
 * @author Richard Holland
 * @since 1.5
 */
public interface PositionResolver {
    
    /**
     * Resolves the minimum possible base for this position.
     * @param start the position to resolve
     * @return the minimum possible base this resolver can return.
     */
    public int getMin(Position start);   
    
    /**
     * Resolves the maximum possible base for this position. 
     * @param end the position to resolve
     * @return the maximum possible base this resolver can return.
     */
    public int getMax(Position end);
    
    /**
     * The maximal resolver returns the base which provides the
     * largest possible range. 
     */
    public static class MaximalResolver implements PositionResolver {
        
        /**
         * {@inheritDoc}
         * ALWAYS RETURNS s.getStart()
         */
        public int getMin(Position s) {
            return s.getStart();
        }
        
        /**
         * {@inheritDoc}
         * ALWAYS RETURNS e.getEnd()
         */
        public int getMax(Position e) {
            return e.getEnd();
        }
    }
        
    /**
     * The minimal resolver returns the base which provides the
     * smallest possible range. 
     */
    public static class MinimalResolver implements PositionResolver {
        
        /**
         * {@inheritDoc}
         * ALWAYS RETURNS s.getEnd()
         */
        public int getMin(Position s) {
            return s.getEnd();
        }
        
        /**
         * {@inheritDoc}
         * ALWAYS RETURNS e.getStart()
         */
        public int getMax(Position e) {
            return e.getStart();
        }
    }    
        
    /**
     * The minimal resolver returns the base which provides the
     * average range, halfway between maximal and minimal. 
     */
    public static class AverageResolver implements PositionResolver {
        
        /**
         * {@inheritDoc}
         * ALWAYS RETURNS s.getStart()+s.getEnd() / 2
         */
        public int getMin(Position s) {
            int min = s.getStart();
            int max = s.getEnd();
            return (min+max) / 2;
        }
        
        /**
         * {@inheritDoc}
         * ALWAYS RETURNS e.getStart()+e.getEnd() / 2
         */
        public int getMax(Position e) {
            int min = e.getStart();
            int max = e.getEnd();
            return (min+max) / 2;
        }
    }
}
