/*
 * BioJava development code
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
 */
package org.biojava.nbio.core.util;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.net.InetSocketAddress;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import com.sun.net.httpserver.HttpServer;

/**
 * Checks that downloads follow the redirects {@link java.net.HttpURLConnection} does
 * not follow by itself.
 * <p>
 * The JDK handles 301, 302 and 303 within a protocol, but never 307 or 308, and never
 * a redirect that changes http to https. Both gaps have broken this project's builds:
 * CATH began answering http with a 301 to https, and ECOD now answers with a 308 to a
 * rewritten path. A browser follows either without comment.
 * <p>
 * The rule tests need no network and no server. The end-to-end tests use a local
 * {@link HttpServer} rather than a real service, so that they cannot fail because a
 * third party is having a bad day.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
class FileDownloadRedirectTest {

	private static final String PAYLOAD = "the file you were looking for\n";

	@Nested
	class RedirectRules {

		private final URL from = url("http://example.org/ecod/distributions/ecod.latest.domains.txt");

		@Test
		void aRelativeLocationIsResolvedAgainstTheRequest() {
			// exactly what ECOD sends: same host, same protocol, relative path
			assertEquals(url("http://example.org/ecod-legacy/distributions/ecod.latest.domains.txt"),
					FileDownloadUtils.redirectTargetFor(308,
							"/ecod-legacy/distributions/ecod.latest.domains.txt", from));
		}

		@Test
		void anAbsoluteLocationIsUsedAsGiven() {
			assertEquals(url("https://example.org/elsewhere.txt"),
					FileDownloadUtils.redirectTargetFor(301, "https://example.org/elsewhere.txt", from));
		}

		@Test
		void everyRedirectStatusWeHandleIsRecognised() {
			for (int code : new int[] { 301, 302, 303, 307, 308 }) {
				assertEquals(url("http://example.org/x"),
						FileDownloadUtils.redirectTargetFor(code, "/x", from),
						"status " + code + " should be followed");
			}
		}

		@Test
		void aSuccessIsNotARedirect() {
			assertNull(FileDownloadUtils.redirectTargetFor(200, null, from));
			assertNull(FileDownloadUtils.redirectTargetFor(404, "/x", from));
		}

		/**
		 * A redirect must never quietly move us onto an unencrypted transport.
		 */
		@Test
		void httpsIsNeverDowngradedToHttp() {
			URL secure = url("https://example.org/file.txt");
			assertNull(FileDownloadUtils.redirectTargetFor(301, "http://example.org/file.txt", secure));
			assertNull(FileDownloadUtils.redirectTargetFor(308, "http://elsewhere.org/file.txt", secure));
		}

		@Test
		void httpToHttpsIsFollowed() {
			// the CATH case
			assertEquals(url("https://example.org/file.txt"),
					FileDownloadUtils.redirectTargetFor(301, "https://example.org/file.txt",
							url("http://example.org/file.txt")));
		}

		@Test
		void anUnusableLocationIsNotFollowed() {
			assertNull(FileDownloadUtils.redirectTargetFor(308, null, from));
			assertNull(FileDownloadUtils.redirectTargetFor(308, "   ", from));
			assertNull(FileDownloadUtils.redirectTargetFor(308, "gopher://example.org/x", from));
		}
	}

	@Nested
	class EndToEnd {

		private HttpServer server;
		private String base;
		private File dir;

		@BeforeEach
		void start() throws IOException {
			server = HttpServer.create(new InetSocketAddress("127.0.0.1", 0), 0);
			base = "http://127.0.0.1:" + server.getAddress().getPort();
			dir = Files.createTempDirectory("redirectTest").toFile();

			serve("/final", 200, null);
			// the ECOD shape: 308 with a relative Location
			serve("/moved", 308, "/final");
			serve("/temp", 307, "/final");
			// a chain that returns to its start
			serve("/loop-a", 308, "/loop-b");
			serve("/loop-b", 308, "/loop-a");
			// longer than the hop limit
			for (int i = 0; i < 9; i++) {
				serve("/hop" + i, 308, "/hop" + (i + 1));
			}
			serve("/hop9", 200, null);
			serve("/nowhere", 308, null);
			server.start();
		}

		private void serve(String path, int status, String location) {
			server.createContext(path, exchange -> {
				byte[] body = PAYLOAD.getBytes(StandardCharsets.UTF_8);
				if (location != null) {
					exchange.getResponseHeaders().add("Location", location);
				}
				exchange.sendResponseHeaders(status, status == 200 ? body.length : -1);
				if (status == 200) {
					try (OutputStream out = exchange.getResponseBody()) {
						out.write(body);
					}
				}
				exchange.close();
			});
		}

		@AfterEach
		void stop() throws IOException {
			server.stop(0);
			FileDownloadUtils.deleteDirectory(dir.getAbsolutePath());
		}

		@Test
		void a308IsFollowed() throws IOException {
			File got = new File(dir, "moved.txt");
			FileDownloadUtils.downloadFile(new URL(base + "/moved"), got);
			assertEquals(PAYLOAD, new String(Files.readAllBytes(got.toPath()), StandardCharsets.UTF_8));
		}

		@Test
		void a307IsFollowed() throws IOException {
			File got = new File(dir, "temp.txt");
			FileDownloadUtils.downloadFile(new URL(base + "/temp"), got);
			assertEquals(PAYLOAD, new String(Files.readAllBytes(got.toPath()), StandardCharsets.UTF_8));
		}

		@Test
		void theRedirectBodyIsNeverWhatWeStore() throws IOException {
			File got = new File(dir, "validated.txt");
			FileDownloadUtils.downloadFileWithValidation(new URL(base + "/moved"), got, null,
					FileDownloadUtils.Hash.UNKNOWN, FileDownloadUtils.ETagPolicy.IGNORE);
			assertEquals(PAYLOAD, new String(Files.readAllBytes(got.toPath()), StandardCharsets.UTF_8));
			assertTrue(FileDownloadUtils.validateFile(got), "the recorded size must describe the real file");
		}

		@Test
		void aLoopIsReportedRatherThanChasedForever() {
			File got = new File(dir, "loop.txt");
			HttpStatusException e = assertThrows(HttpStatusException.class,
					() -> FileDownloadUtils.downloadFile(new URL(base + "/loop-a"), got));
			assertTrue(e.getMessage().contains("loop"), e.getMessage());
		}

		@Test
		void tooManyHopsGivesUp() {
			File got = new File(dir, "hops.txt");
			assertThrows(HttpStatusException.class,
					() -> FileDownloadUtils.downloadFile(new URL(base + "/hop0"), got));
		}

		@Test
		void aRedirectWithNoDestinationIsAnError() {
			File got = new File(dir, "nowhere.txt");
			assertThrows(HttpStatusException.class,
					() -> FileDownloadUtils.downloadFile(new URL(base + "/nowhere"), got));
		}
	}

	private static URL url(String spec) {
		try {
			return new URL(spec);
		} catch (IOException e) {
			throw new IllegalArgumentException(spec, e);
		}
	}
}
