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
public interface AccesionedLocation extends Accessioned {

    /**
     * Return the proxy reader used to get sequence for this location. We
     * assume that AccessionID is bound to the instance to facilitate this
     * lookup.
     */
    ProxySequenceReader<? extends Compound> getProxySequenceReader();
}
