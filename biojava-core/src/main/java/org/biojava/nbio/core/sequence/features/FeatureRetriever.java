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
 * Created on 11-30-2014
 */

package org.biojava.nbio.core.sequence.features;

import java.util.ArrayList;
import java.util.HashMap;
/**
 * If a SequenceProxyReader implements this interface then that external source
 * has a list features
 * @author @author Paolo Pavan
 */
public interface FeatureRetriever {
    HashMap<String, ArrayList<AbstractFeature>> getFeatures();
}
