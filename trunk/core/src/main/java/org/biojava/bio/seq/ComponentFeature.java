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

package org.biojava.bio.seq;

import org.biojava.bio.symbol.Location;

/**
 * Feature which represents a component in an assembly (contig).
 * This implies that a portion (possibly all) of the 
 * associated componentSequence is included in this feature's
 * parent sequence.
 *
 * <p>
 * There are important invariants which apply to all ComponentFeatures.
 * The Location returned by getLocation() must
 * contain the same number of unique point locations as that
 * returned by getComponentLocation().
 * </p>
 *
 * <p>
 * In BioJava 1.2, two extra properties were added to <code>ComponentFeature</code>
 * to support the use of these features in environments were it it necessary to
 * represent an assembly built from sequences which are not currently available.
 * Widespread use of such sequences is not encouraged, but it is recognized that
 * they are useful as intermediate objects for data integration applications.
 * </p>
 *
 * @author Thomas Down
 * @since 1.1
 */

public interface ComponentFeature extends StrandedFeature {
    /**
     * Get the sequence object which provides a component of this
     * feature's parent sequence.
     *
     * @return A sequence.
     */

    public Sequence getComponentSequence();

    /**
     * Return a location which identifies a portion of the component
     * sequence which is to be included in the assembly.
     *
     * @return A location within the component sequence.
     */

    public Location getComponentLocation();

    /**
     * Get the name of the component sequence.  In general, this
     * should be equivalent to
     * <code>getComponentSequence().getName()</code>.  However,
     * it may still be defined for un-resolveable components.
     *
     * @since 1.2
     */

    public String getComponentSequenceName();

    /**
     * Determine if the sequence references by this ComponentFeature
     * is available in this context.  If not, calls to getComponentSequence will
     * fail, and getSymbols will return a non-informative 
     * <code>SymbolList</code> (in a DNA context, a list of Ns).
     *
     * @since 1.2
     */

    public boolean isComponentResolvable();

    /**
     * Template for constructing a new ComponentFeature.
     *
     * <p>
     * In BioJava 1.2, we add the <code>componentSequenceName</code>
     * property.  In general, it is only necessary to specify
     * <em>either</em> componentSequenceName or componentSequence
     * </p>
     */

    public static class Template extends StrandedFeature.Template {
	public String   componentSequenceName;
	public Sequence componentSequence;
	public Location componentLocation;
    }
}
