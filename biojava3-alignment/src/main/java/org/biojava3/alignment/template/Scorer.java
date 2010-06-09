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
 * Defines an algorithm which computes a score.
 *
 * @author Mark Chapman
 */
public interface Scorer {

    /**
     * Returns maximum possible score.
     *
     * @return maximum possible score
     */
    int getMaxScore();

    /**
     * Returns minimum possible score.
     *
     * @return minimum possible score
     */
    int getMinScore();

    /**
     * Returns score resulting from algorithm.  This should normalize between 0 and 1 by calculating
     * ({@link #getScore()} - {@link #getMinScore()}) / ({@link #getMaxScore()} - {@link #getMinScore()}).
     *
     * @return score resulting from algorithm
     */
    int getScore();

}
