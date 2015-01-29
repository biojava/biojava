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
 * Created on July 22, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment.template;

/**
 * Implements common code for algorithms which compute a score.
 *
 * @author Mark Chapman
 */
public abstract class AbstractScorer implements Scorer {

    @Override
    public double getDistance() {
        return getDistance(1.0);
    }

    @Override
    public double getDistance(double scale) {
        return scale * (getMaxScore() - getScore()) / (getMaxScore() - getMinScore());
    }

    @Override
    public double getSimilarity() {
        return getSimilarity(1.0);
    }

    @Override
    public double getSimilarity(double scale) {
        return scale * (getScore() - getMinScore()) / (getMaxScore() - getMinScore());
    }

}
