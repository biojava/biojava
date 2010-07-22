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
     * Returns score as a distance between 0.0 and 1.0.  This equals ({@link #getMaxScore()} - {@link #getScore()}) /
     * ({@link #getMaxScore()} - {@link #getMinScore()}).
     *
     * @return score as a distance between 0.0 and 1.0
     */
    double getDistance();

    /**
     * Returns score as a distance between 0.0 and scale.  This equals scale * ({@link #getMaxScore()} -
     * {@link #getScore()}) / ({@link #getMaxScore()} - {@link #getMinScore()}).
     *
     * @param scale maximum distance
     * @return score as a distance between 0.0 and scale
     */
    double getDistance(double scale);

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

    /**
     * Returns score as a similarity between 0.0 and 1.0.  This equals ({@link #getScore()} - {@link #getMinScore()}) /
     * ({@link #getMaxScore()} - {@link #getMinScore()}).
     *
     * @return score as a similarity between 0.0 and 1.0
     */
    double getSimilarity();

    /**
     * Returns score as a similarity between 0.0 and scale.  This equals scale * ({@link #getScore()} -
     * {@link #getMinScore()}) / ({@link #getMaxScore()} - {@link #getMinScore()}).
     *
     * @param scale maximum similarity
     * @return score as a similarity between 0.0 and scale
     */
    double getSimilarity(double scale);

}
