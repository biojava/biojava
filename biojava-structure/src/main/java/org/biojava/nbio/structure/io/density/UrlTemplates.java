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

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Expands the named placeholders used in the configurable density-server URL
 * templates.
 * <p>
 * The recognised placeholders are:
 * <table summary="placeholders">
 * <tr><td><code>{pdbid}</code></td><td>the PDB identifier as given</td></tr>
 * <tr><td><code>{pdbid_lc}</code></td><td>the PDB identifier in lower case</td></tr>
 * <tr><td><code>{pdbid_uc}</code></td><td>the PDB identifier in upper case</td></tr>
 * <tr><td><code>{mid}</code></td><td>the two-character divided-archive directory,
 * e.g. <code>cb</code> for <code>1cbs</code></td></tr>
 * <tr><td><code>{emdb_id}</code></td><td>the EMDB identifier, e.g. <code>EMD-0262</code></td></tr>
 * <tr><td><code>{emdb_num}</code></td><td>the EMDB number alone, e.g. <code>0262</code></td></tr>
 * <tr><td><code>{detail}</code></td><td>the density-server detail level</td></tr>
 * </table>
 * A placeholder with no supplied value is left in place rather than replaced by
 * an empty string, so that a misconfigured template produces an obviously wrong
 * URL instead of a subtly wrong one.
 * <p>
 * This is deliberately separate from the template mechanism in
 * {@code DownloadChemCompProvider}, which resolves a single chemical-component
 * identifier with optional substring indices. The two solve different problems:
 * here several distinct values are substituted by name.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class UrlTemplates {

	private static final Pattern PLACEHOLDER = Pattern.compile("\\{([a-zA-Z_]+)\\}");

	private UrlTemplates() {
	}

	/**
	 * Expands a template against a set of named values.
	 *
	 * @param template the template string
	 * @param values the values, keyed by placeholder name without the braces
	 * @return the expanded string
	 */
	public static String expand(String template, Map<String, String> values) {
		if (template == null) {
			return null;
		}
		Matcher m = PLACEHOLDER.matcher(template);
		StringBuilder out = new StringBuilder();
		int last = 0;
		while (m.find()) {
			String value = values.get(m.group(1));
			out.append(template, last, m.start());
			// An unknown placeholder is kept verbatim: a template that was configured
			// wrongly should fail loudly rather than quietly produce a plausible URL.
			out.append(value == null ? m.group(0) : value);
			last = m.end();
		}
		out.append(template.substring(last));
		return out.toString();
	}

	/**
	 * Builds the standard value map for an entry.
	 *
	 * @param pdbId the PDB identifier, may be <code>null</code>
	 * @param emdbId the EMDB identifier in any accepted form, may be <code>null</code>
	 * @param detail the density-server detail level, or a negative value to omit it
	 * @return a map suitable for {@link #expand(String, Map)}
	 */
	public static Map<String, String> values(String pdbId, String emdbId, int detail) {
		Map<String, String> values = new LinkedHashMap<>();
		if (pdbId != null) {
			values.put("pdbid", pdbId);
			values.put("pdbid_lc", pdbId.toLowerCase());
			values.put("pdbid_uc", pdbId.toUpperCase());
			if (pdbId.length() >= 3) {
				values.put("mid", org.biojava.nbio.structure.io.LocalPDBDirectory.getMiddleHash(pdbId));
			}
		}
		if (emdbId != null) {
			values.put("emdb_id", DensityMapRequest.normalizeEmdbId(emdbId));
			values.put("emdb_num", DensityMapRequest.emdbNumber(emdbId));
		}
		if (detail >= 0) {
			values.put("detail", Integer.toString(detail));
		}
		return values;
	}
}
