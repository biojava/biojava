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

import java.io.File;
import java.util.Collections;
import java.util.List;

import org.biojava.nbio.structure.PdbId;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;

/**
 * A request for a density map, describing what is wanted and how hard to look
 * for it.
 * <p>
 * Instances are immutable; build them with {@link #builder(PdbId)} or
 * {@link #builder(String)}.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class DensityMapRequest {

	private final PdbId pdbId;
	private final String emdbId;
	private final DensityMapKind kind;
	private final FetchBehavior fetchBehavior;
	private final File cacheDir;
	private final boolean allowNonRenderableFormats;
	private final long maxDownloadBytes;
	private final List<DensityMapSource> sourceChain;

	private DensityMapRequest(Builder b) {
		this.pdbId = b.pdbId;
		this.emdbId = b.emdbId;
		this.kind = b.kind;
		this.fetchBehavior = b.fetchBehavior;
		this.cacheDir = b.cacheDir;
		this.allowNonRenderableFormats = b.allowNonRenderableFormats;
		this.maxDownloadBytes = b.maxDownloadBytes;
		this.sourceChain = b.sourceChain == null ? null : Collections.unmodifiableList(b.sourceChain);
	}

	/**
	 * Starts a request for a PDB entry.
	 *
	 * @param pdbId the entry
	 * @return a new builder
	 */
	public static Builder builder(PdbId pdbId) {
		return new Builder(pdbId, null);
	}

	/**
	 * Starts a request from an identifier string. An identifier beginning with
	 * <code>EMD-</code> (case-insensitively) is taken as an EMDB entry, anything
	 * else as a PDB entry.
	 *
	 * @param id a PDB or EMDB identifier
	 * @return a new builder
	 */
	public static Builder builder(String id) {
		if (id == null) {
			throw new IllegalArgumentException("Identifier must not be null");
		}
		String trimmed = id.trim();
		if (trimmed.toUpperCase().startsWith("EMD-") || trimmed.toUpperCase().startsWith("EMD_")) {
			return new Builder(null, normalizeEmdbId(trimmed)).kind(DensityMapKind.EM);
		}
		return new Builder(new PdbId(trimmed), null);
	}

	/**
	 * Normalises an EMDB identifier to the canonical <code>EMD-1234</code> form.
	 *
	 * @param emdbId an identifier such as <code>emd-1234</code>, <code>EMD_1234</code>
	 *        or a bare number
	 * @return the canonical form
	 */
	public static String normalizeEmdbId(String emdbId) {
		if (emdbId == null) {
			return null;
		}
		String digits = emdbId.trim().toUpperCase().replaceFirst("^EMD[-_]?", "");
		return "EMD-" + digits;
	}

	/**
	 * Extracts the numeric part of an EMDB identifier, as used in file names and
	 * some URL templates.
	 *
	 * @param emdbId an EMDB identifier in any accepted form
	 * @return the digits, without the <code>EMD-</code> prefix
	 */
	public static String emdbNumber(String emdbId) {
		return normalizeEmdbId(emdbId).substring("EMD-".length());
	}

	/** @return the PDB entry requested, or <code>null</code> for an EMDB-only request */
	public PdbId getPdbId() {
		return pdbId;
	}

	/** @return the EMDB entry, if known or explicitly requested; otherwise <code>null</code> */
	public String getEmdbId() {
		return emdbId;
	}

	/** @return the kind of map wanted; never <code>null</code> */
	public DensityMapKind getKind() {
		return kind;
	}

	/** @return how aggressively to re-fetch, or <code>null</code> to use the cache's setting */
	public FetchBehavior getFetchBehavior() {
		return fetchBehavior;
	}

	/** @return an override for the cache directory, or <code>null</code> to use the cache's own */
	public File getCacheDir() {
		return cacheDir;
	}

	/**
	 * Whether sources that deliver something a viewer cannot contour directly
	 * &mdash; currently only map coefficients &mdash; may be used.
	 *
	 * @return <code>true</code> by default; a viewer should set it to
	 *         <code>false</code>
	 */
	public boolean isAllowNonRenderableFormats() {
		return allowNonRenderableFormats;
	}

	/** @return the download size limit in bytes, or 0 for no limit */
	public long getMaxDownloadBytes() {
		return maxDownloadBytes;
	}

	/** @return an explicit source order, or <code>null</code> to let the cache choose */
	public List<DensityMapSource> getSourceChain() {
		return sourceChain;
	}

	/**
	 * Returns a copy of this request with a different map kind, used when expanding
	 * {@link DensityMapKind#AUTO}.
	 *
	 * @param newKind the kind to use
	 * @return a new request
	 */
	public DensityMapRequest withKind(DensityMapKind newKind) {
		return toBuilder().kind(newKind).build();
	}

	/**
	 * Returns a copy of this request with the EMDB entry filled in.
	 *
	 * @param newEmdbId the EMDB identifier
	 * @return a new request
	 */
	public DensityMapRequest withEmdbId(String newEmdbId) {
		return toBuilder().emdbId(newEmdbId).build();
	}

	private Builder toBuilder() {
		Builder b = new Builder(pdbId, emdbId);
		b.kind = kind;
		b.fetchBehavior = fetchBehavior;
		b.cacheDir = cacheDir;
		b.allowNonRenderableFormats = allowNonRenderableFormats;
		b.maxDownloadBytes = maxDownloadBytes;
		b.sourceChain = sourceChain;
		return b;
	}

	@Override
	public String toString() {
		return String.format("DensityMapRequest[%s%s, %s]",
				pdbId == null ? "" : pdbId.getId(),
				emdbId == null ? "" : (pdbId == null ? emdbId : "/" + emdbId),
				kind);
	}

	/**
	 * Builder for {@link DensityMapRequest}.
	 *
	 * @author Amr ALHOSSARY
	 * @since 7.3.0
	 */
	public static final class Builder {

		private final PdbId pdbId;
		private String emdbId;
		private DensityMapKind kind = DensityMapKind.AUTO;
		private FetchBehavior fetchBehavior;
		private File cacheDir;
		private boolean allowNonRenderableFormats = true;
		private long maxDownloadBytes = -1;
		private List<DensityMapSource> sourceChain;

		private Builder(PdbId pdbId, String emdbId) {
			this.pdbId = pdbId;
			this.emdbId = emdbId;
		}

		/** @param kind the kind of map wanted; <code>null</code> means {@link DensityMapKind#AUTO} */
		public Builder kind(DensityMapKind kind) {
			this.kind = kind == null ? DensityMapKind.AUTO : kind;
			return this;
		}

		/** @param emdbId the EMDB entry to use, skipping the PDB-to-EMDB lookup */
		public Builder emdbId(String emdbId) {
			this.emdbId = emdbId == null ? null : normalizeEmdbId(emdbId);
			return this;
		}

		/** @param fetchBehavior how aggressively to re-fetch */
		public Builder fetchBehavior(FetchBehavior fetchBehavior) {
			this.fetchBehavior = fetchBehavior;
			return this;
		}

		/** @param cacheDir an override for the cache directory */
		public Builder cacheDir(File cacheDir) {
			this.cacheDir = cacheDir;
			return this;
		}

		/**
		 * @param allow whether formats that cannot be contoured directly may be used.
		 *        A viewer should pass <code>false</code>.
		 */
		public Builder allowNonRenderableFormats(boolean allow) {
			this.allowNonRenderableFormats = allow;
			return this;
		}

		/** @param bytes the download size limit, or 0 for no limit, or a negative value to use the cache's setting */
		public Builder maxDownloadBytes(long bytes) {
			this.maxDownloadBytes = bytes;
			return this;
		}

		/** @param chain an explicit source order, overriding the cache's choice */
		public Builder sourceChain(List<DensityMapSource> chain) {
			this.sourceChain = chain;
			return this;
		}

		/** @return the finished request */
		public DensityMapRequest build() {
			if (pdbId == null && emdbId == null) {
				throw new IllegalArgumentException("A request needs either a PDB ID or an EMDB ID");
			}
			return new DensityMapRequest(this);
		}
	}
}
