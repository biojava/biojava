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
 * Created on 01-21-2010
 */
package org.biojava3.core.sequence.location.template;

import org.biojava3.core.sequence.template.Accessioned;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.ProxySequenceReader;

/**
 * A location which is bound to an AccessionID. This is common in 
 * INSDC locations where a location actually points to a remote sequence. This
 * is especially common in records describing Genomic sequence assembly.
 *
 * @author ayates
 */
public interface AccesionedLocation extends Location, Accessioned {

    /**
     * Return the proxy reader used to get sequence for this location. We
     * assume that AccessionID is bound to the instance to facilitate this
     * lookup.
     */
    ProxySequenceReader<? extends Compound> getProxySequenceReader();
}
