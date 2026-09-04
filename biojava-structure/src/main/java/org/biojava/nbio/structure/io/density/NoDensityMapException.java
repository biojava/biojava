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
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;

import org.biojava.nbio.structure.PdbId;

/**
 * Thrown when no enabled source could supply a density map for an entry.
 * <p>
 * This is a routine outcome rather than an error: not every entry has density.
 * Structures deposited without structure factors have none at all (4HHB, for
 * instance), and a cryo-EM entry has no X-ray maps by construction. The
 * per-source reasons are kept so that a caller &mdash; particularly a user
 * interface &mdash; can explain <i>why</i> rather than just reporting a failure.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class NoDensityMapException extends IOException {

	private static final long serialVersionUID = 1L;

	private final PdbId pdbId;
	private final DensityMapKind kind;
	private final Map<DensityMapSource, String> attempts;

	/**
	 * @param pdbId the entry that was requested, may be <code>null</code> for a
	 *        request made by EMDB identifier
	 * @param kind the kind of map that was requested
	 * @param attempts what happened at each source that was tried, in the order tried
	 */
	public NoDensityMapException(PdbId pdbId, DensityMapKind kind, Map<DensityMapSource, String> attempts) {
		super(buildMessage(pdbId, kind, attempts));
		this.pdbId = pdbId;
		this.kind = kind;
		this.attempts = attempts == null
				? Collections.emptyMap()
				: Collections.unmodifiableMap(new LinkedHashMap<>(attempts));
	}

	private static String buildMessage(PdbId pdbId, DensityMapKind kind, Map<DensityMapSource, String> attempts) {
		StringBuilder sb = new StringBuilder("No ");
		sb.append(kind == null ? "density" : kind).append(" map is available for ");
		sb.append(pdbId == null ? "the requested entry" : pdbId.getId()).append('.');
		if (attempts != null && !attempts.isEmpty()) {
			sb.append(" Tried:");
			for (Map.Entry<DensityMapSource, String> e : attempts.entrySet()) {
				sb.append(' ').append(e.getKey()).append(" (").append(e.getValue()).append(");");
			}
			sb.setLength(sb.length() - 1);
		}
		return sb.toString();
	}

	/**
	 * @return the entry that was requested, or <code>null</code>
	 */
	public PdbId getPdbId() {
		return pdbId;
	}

	/**
	 * @return the kind of map that was requested
	 */
	public DensityMapKind getKind() {
		return kind;
	}

	/**
	 * What happened at each source, in the order they were tried. Suitable for
	 * building a human-readable explanation.
	 *
	 * @return an unmodifiable map from source to reason
	 */
	public Map<DensityMapSource, String> getAttempts() {
		return attempts;
	}
}
