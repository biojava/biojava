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
package org.biojava.nbio.structure.chem;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.net.InetSocketAddress;
import java.nio.charset.StandardCharsets;

import com.sun.net.httpserver.HttpServer;

import org.biojava.nbio.core.util.FlatFileCache;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

/**
 * A server that redirects, or errors, must not have its response body cached as
 * though it were a chemical component definition.
 * <p>
 * This is the failure that took the CATH downloader out when
 * <code>download.cathdb.info</code> moved to https, and the chem comp download had
 * the same shape. A 4xx already failed safely, because
 * <code>getInputStream()</code> throws for those; a redirect did not, because when
 * the JDK declines to follow a 3xx it hands back the redirect's body instead, and
 * that body is short but not empty.
 * <p>
 * The test serves the responses from a local {@link HttpServer} rather than a real
 * service. Pointing it at a third-party server that happens to redirect today would
 * make the test fail on the day they stop, which is precisely the coupling that made
 * the build unreliable in the first place.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
public class TestChemCompRedirectNotCached {

	private HttpServer server;
	private String originalServerUrl;

	@Before
	public void setUp() {
		originalServerUrl = DownloadChemCompProvider.serverBaseUrl;
	}

	@After
	public void tearDown() {
		if (server != null) {
			server.stop(0);
		}
		// Static state: leaving either of these set would corrupt unrelated tests.
		DownloadChemCompProvider.serverBaseUrl = originalServerUrl;
		FlatFileCache.clear();
	}

	/**
	 * Starts a local server that answers every request with the given status and body.
	 *
	 * @return the base URL to point the provider at
	 */
	private String startServer(int status, String location, String body) throws IOException {
		server = HttpServer.create(new InetSocketAddress("127.0.0.1", 0), 0);
		server.createContext("/", exchange -> {
			byte[] bytes = body.getBytes(StandardCharsets.UTF_8);
			if (location != null) {
				exchange.getResponseHeaders().add("Location", location);
			}
			exchange.sendResponseHeaders(status, bytes.length);
			try (OutputStream out = exchange.getResponseBody()) {
				out.write(bytes);
			}
		});
		server.start();
		return "http://127.0.0.1:" + server.getAddress().getPort() + "/";
	}

	private File cacheFileFor(String id) {
		File file = new File(DownloadChemCompProvider.getLocalFileName(id));
		file.delete();
		FlatFileCache.clear();
		return file;
	}

	/**
	 * The case that broke CATH: a redirect the JDK will not follow because it
	 * changes protocol. Its body must not end up on disk under the component's name.
	 */
	@Test
	public void redirectBodyIsNotCached() throws IOException {
		File cached = cacheFileFor("ATP");
		DownloadChemCompProvider.serverBaseUrl =
				startServer(301, "https://example.invalid/ATP.cif", "<html>Moved Permanently</html>");

		ChemComp cc = new DownloadChemCompProvider().getChemComp("ATP");

		assertFalse("the body of a redirect must never be cached as a definition", cached.exists());
		assertNull("nothing parseable was returned, so the component must be empty", cc.getName());
	}

	/**
	 * A 200 is still cached, so the guard has not simply disabled downloading.
	 * <p>
	 * What is under test is the download path, not the CIF parser: the response is
	 * written to the cache before anything tries to parse it, so a parse failure on
	 * this deliberately minimal body says nothing about whether the guard behaved.
	 */
	@Test
	public void aValidResponseIsStillCached() throws IOException {
		File cached = cacheFileFor("ATP");
		DownloadChemCompProvider.serverBaseUrl = startServer(200, null,
				"data_ATP\n#\n_chem_comp.id ATP\n_chem_comp.name \"ADENOSINE-5'-TRIPHOSPHATE\"\n#\n");

		try {
			new DownloadChemCompProvider().getChemComp("ATP");
		} catch (RuntimeException parseFailure) {
			// see the note above
		}

		assertTrue("a 200 response should still be cached", cached.exists());
		cached.delete();
	}

	/** A server error must not be cached either. */
	@Test
	public void serverErrorBodyIsNotCached() throws IOException {
		File cached = cacheFileFor("ATP");
		DownloadChemCompProvider.serverBaseUrl =
				startServer(503, null, "<html>Service Unavailable</html>");

		new DownloadChemCompProvider().getChemComp("ATP");

		assertFalse("the body of a 5xx must never be cached as a definition", cached.exists());
	}
}
