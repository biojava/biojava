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
package org.biojava3.core.io;

import java.io.Serializable;
import java.util.Properties;

/**
 * Receives data about things in one format and retransmits them in another.
 * @author Richard Holland
 * @since 3.0
 */
public interface ThingConverter<T extends Serializable> extends ThingParserMember, ThingBuilder<T> {

    /**
     * Conversion rules control how the conversion takes place. For instance,
     * a FASTA to Genbank converter might have a rule specifying a regex to
     * use to convert the FASTA description line, and further rules to tell it
     * which bits of the Genbank option to plug the regex match components into.
     * @param rules the rules for conversion. If null, then any default rules
     * the converter may define will be used instead. If not null, then any 
     * default rules not named in this new set of rules will still be in force.
     */
    public void setConversionRules(Properties rules);

    /**
     * @see {#setConversionRules(Properties)}
     * @return rules the rules for conversion. If null, then no rules apply.
     * Otherwise the returned set is the combination of any defaults, plus any
     * other rules or overriding rules that have been set.
     */
    public Properties getConversionRules();
}