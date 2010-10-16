package org.biojava3.core.sequence.template;

import org.biojava3.core.sequence.AccessionID;

/**
 * Indicates an entity is accessioned
 *
 * @author ayates
 */
public interface Accessioned {

    /**
     * Returns the AccessionID this location is currently bound with
     */
    AccessionID getAccession();

}
