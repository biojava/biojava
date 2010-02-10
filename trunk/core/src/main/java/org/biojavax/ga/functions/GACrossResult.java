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

package org.biojavax.ga.functions;

import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.SymbolList;

/**
 * <p>Holds the results of a CrossOver event, objects of this type are made by
 * <code>CrossOverFunctions</code></p>
 * 
 * @author Mark Schreiber
 * @version 1.0
 * @since 1.5
 */

public interface GACrossResult {

/**
 * Returns the collection of cross over locations from the last cross
 * @return null if there has been no call to performCrossOver, a zero length
 * array if there was no cross overs or an array of <code>PointLocations</code>
 * describing the cross over points.
 */
  public PointLocation[] getCrossOverPositions();

/**
 * Gets the chromosomes after the cross
 * @return the two Chromosomes, by convention the first will be chromA and
 * the second chromB
 */
  public SymbolList[] getChromosomes();
}