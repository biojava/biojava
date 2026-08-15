/**
 * BioJava development code
 *
 * This code may be freely distributed and modified under the terms of the GNU
 * Lesser General Public Licence. This should be distributed with the code. If
 * you do not have a copy, see:
 *
 * http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual authors. These
 * should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims, or to join the
 * biojava-l mailing list, visit the home page at:
 *
 * http://www.biojava.org/
 */
package org.biojava.nbio.structure.io.density;

import java.io.IOException;

import org.biojava.nbio.core.util.HttpStatusException;

/**
 * Fetches density maps from one particular service.
 * <p>
 * Implementations are combined into an ordered chain by {@link DensityMapCache},
 * which tries each in turn until one produces a map. The contract around
 * exceptions is what makes that chain safe:
 * <ul>
 * <li>Throw {@link HttpStatusException} with {@link HttpStatusException#isNotFound()}
 * &mdash; or return <code>null</code> &mdash; when this service simply has nothing
 * for the entry. The chain moves on to the next source.</li>
 * <li>Throw {@link DensityMapTooLargeException} when the map exists but exceeds
 * the caller's size limit. The chain also moves on, so a smaller representation
 * from another source can be used instead.</li>
 * <li>Throw any other {@link IOException} for a genuine transport failure. The
 * chain stops, so that a network outage is never mistaken for "this entry has no
 * density".</li>
 * </ul>
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public interface DensityMapProvider {

	/**
	 * @return which service this provider talks to
	 */
	DensityMapSource getSource();

	/**
	 * @return the format this provider delivers
	 */
	DensityFileFormat getFormat();

	/**
	 * Whether this provider can serve a given kind of map at all. Used to skip
	 * requests that could not possibly succeed, such as asking an X-ray map service
	 * for a cryo-EM map.
	 *
	 * @param kind the kind of map wanted
	 * @return <code>true</code> if it is worth trying
	 */
	boolean supports(DensityMapKind kind);

	/**
	 * Fetches a map, using the cache if the request's fetch behaviour allows.
	 *
	 * @param request what is wanted. Its kind is always concrete, never
	 *        {@link DensityMapKind#AUTO}.
	 * @return the map, or <code>null</code> if this service has nothing for the entry
	 * @throws DensityMapTooLargeException if the map exceeds the request's size limit
	 * @throws IOException on transport failure; use {@link HttpStatusException} so
	 *         that a missing resource can be told apart from a broken connection
	 */
	DensityMapResult fetch(DensityMapRequest request) throws IOException;
}
