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
package org.biojava.nbio.core.util;

import java.io.IOException;

/**
 * Signals that an HTTP request completed but returned a status code outside the
 * 2xx range.
 * <p>
 * This exists so that callers can tell apart the two very different reasons a
 * download can fail:
 * <ul>
 * <li><b>The resource does not exist</b> ({@link #isNotFound()}, i.e. HTTP 404 or
 * 410). For a caller that tries several mirrors or several alternative data
 * sources in turn, this simply means "not here, try the next one".</li>
 * <li><b>Anything else</b> &mdash; a server error, a redirect that was not
 * followed, an authentication failure. These usually mean the whole operation
 * should be abandoned rather than silently treated as "no data available".</li>
 * </ul>
 * Without a distinct exception type the only way to tell these apart is by
 * parsing the message of a plain {@link IOException}, which is brittle.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class HttpStatusException extends IOException {

	private static final long serialVersionUID = 1L;

	private final int statusCode;
	private final String url;

	/**
	 * @param statusCode the HTTP status code returned by the server
	 * @param url the URL that was requested
	 * @param responseMessage the HTTP reason phrase, may be <code>null</code>
	 */
	public HttpStatusException(int statusCode, String url, String responseMessage) {
		super(String.format("HTTP %d%s for %s", statusCode,
				responseMessage == null || responseMessage.isEmpty() ? "" : " " + responseMessage, url));
		this.statusCode = statusCode;
		this.url = url;
	}

	/**
	 * @return the HTTP status code returned by the server
	 */
	public int getStatusCode() {
		return statusCode;
	}

	/**
	 * @return the URL that was requested
	 */
	public String getUrl() {
		return url;
	}

	/**
	 * Whether the status indicates that the resource is simply not there, as
	 * opposed to a transport or server problem.
	 *
	 * @return <code>true</code> for HTTP 404 (Not Found) and 410 (Gone)
	 */
	public boolean isNotFound() {
		return statusCode == 404 || statusCode == 410;
	}
}
