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
 *
 * Created on Feb 23, 2012 Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.nbio.core.util;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.net.HttpURLConnection;
import java.net.SocketTimeoutException;
import java.net.URL;
import java.net.URLConnection;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;
import java.security.DigestInputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class FileDownloadUtils {

	private static final String SIZE_EXT = ".size";
	private static final String HASH_EXT = ".hash";
	private static final Logger logger = LoggerFactory.getLogger(FileDownloadUtils.class);

	/** Buffer used when streaming a file through a {@link MessageDigest}. */
	private static final int DIGEST_BUFFER_SIZE = 64 * 1024;

	/** A bare hex digest, optionally followed by whitespace and a file name (the
	 * layout written by <code>md5sum</code>, <code>sha1sum</code> and friends). */
	private static final Pattern BARE_HEX_HASH = Pattern.compile("^([0-9a-fA-F]{32,128})(?:[\\s*].*)?$");

	/** The BSD layout, e.g. <code>MD5 (file.txt) = d41d8cd9...</code>. */
	private static final Pattern BSD_HASH = Pattern.compile("^\\w+\\s*\\(.*\\)\\s*=\\s*([0-9a-fA-F]{32,128})$");

	public enum Hash{
		MD5, SHA1, SHA256, UNKNOWN
	}

	/**
	 * What to do with the <code>ETag</code> response header when no explicit hash
	 * URL is available.
	 * <p>
	 * Some archives &mdash; notably <code>files.wwpdb.org</code> and
	 * <code>files.rcsb.org</code> &mdash; return the MD5 digest of the content as
	 * the bare <code>ETag</code> value, which lets us record a real checksum
	 * without a second request. Others (the EBI servers, for instance) return a
	 * <code>&lt;modification-time&gt;-&lt;size&gt;</code> form instead; because that
	 * always contains a <code>-</code>, it can never be mistaken for a hex digest.
	 *
	 * @author Amr ALHOSSARY
	 * @since 7.3.0
	 */
	public enum ETagPolicy {
		/** Never look at the <code>ETag</code> header. */
		IGNORE,
		/** Record the <code>ETag</code> as a checksum when it is a bare hex digest
		 * of a length matching one of the supported algorithms. */
		USE_IF_HEX_DIGEST,
		/** As {@link #USE_IF_HEX_DIGEST}, but log a warning when the header is
		 * missing or is not a usable digest. */
		REQUIRE
	}

	/**
	 * Gets the file extension of a file, excluding '.'.
	 * If the file name has no extension the file name is returned.
	 * @param f a File
	 * @return The extension
	 */
	public static String getFileExtension(File f) {
		String fileName = f.getName();
		String ext = "";
		int mid = fileName.lastIndexOf(".");
		ext = fileName.substring(mid + 1);
		return ext;
	}

	/**
	 * Gets the file name up to and excluding the first
	 * '.' character. If there is no extension, the full filename
	 * is returned.
	 * @param f A file
	 * @return A possibly empty but non-null String.
	 */
	public static String getFilePrefix(File f) {
		String fileName = f.getName();
		int mid = fileName.indexOf(".");
		if (mid < 0) {
			return fileName;
		}
		return fileName.substring(0, mid);
	}

	/**
	 * Download the content provided at URL url and store the result to a local
	 * file, using a temp file to cache the content in case something goes wrong
	 * in download. A timeout of 60 seconds is hard-coded and 10 retries are attempted.
	 *
	 * @param url
	 * @param destination
	 * @throws IOException
	 */
	public static void downloadFile(URL url, File destination) throws IOException {
		int count = 0;
		int maxTries = 10;
		int timeout = 60000; //60 sec

		File tempFile = createTempFileFor(destination);

		try {
			while (true) {
				try {
					URLConnection connection = prepareURLConnection(url.toString(), timeout);
					connection.connect();
					checkHttpStatus(connection);
					try (InputStream inputStream = connection.getInputStream()) {
						// Files.copy loops until end of stream. FileChannel.transferFrom(), used
						// here previously, is not guaranteed to drain a socket-backed channel in
						// a single call and could silently truncate a download.
						Files.copy(inputStream, tempFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
					}
					break;
				} catch (SocketTimeoutException e) {
					if (++count == maxTries) throw e;
				}
			}

			logger.debug("Copying temp file [{}] to final location [{}]", tempFile, destination);
			Files.copy(tempFile.toPath(), destination.toPath(), StandardCopyOption.REPLACE_EXISTING);
		} finally {
			// on every path, including failure: the temp file used to leak whenever
			// the download threw.
			deleteQuietly(tempFile);
		}
	}

	/**
	 * Downloads a file and writes its validation metadata in a single pass over a
	 * single connection.
	 * <p>
	 * This is preferable to calling {@link #createValidationFiles(URL, File, URL, Hash)}
	 * followed by {@link #downloadFile(URL, File)}: those open separate connections,
	 * so the <code>Content-Length</code> recorded in the <code>.size</code> file
	 * comes from a different response than the bytes actually written. If the
	 * resource changes between the two requests, the cache entry is left
	 * permanently failing validation.
	 * <p>
	 * The content is streamed to a temporary file and only moved into place once
	 * the declared length and (where available) the checksum have been confirmed,
	 * so a failed download never leaves a partial file at <code>destination</code>.
	 *
	 * @param url the remote file to download
	 * @param destination the local file to download into. Its parent directory must exist.
	 * @param hashURL URL of a file containing the expected hash. May be <code>null</code>.
	 * @param hash the hashing algorithm matching <code>hashURL</code>. Ignored when
	 *             <code>hashURL</code> is <code>null</code>.
	 * @param eTagPolicy what to do with an <code>ETag</code> response header when
	 *             <code>hashURL</code> is <code>null</code>. May be <code>null</code>,
	 *             which is treated as {@link ETagPolicy#IGNORE}.
	 * @throws HttpStatusException if the server answered with a non-2xx status
	 * @throws IOException if the transfer failed, or the transferred content did not
	 *             match the length or checksum the server declared
	 * @author Amr ALHOSSARY
	 * @since 7.3.0
	 */
	public static void downloadFileWithValidation(URL url, File destination, URL hashURL, Hash hash,
			ETagPolicy eTagPolicy) throws IOException {
		int timeout = 60000; //60 sec
		ETagPolicy policy = eTagPolicy == null ? ETagPolicy.IGNORE : eTagPolicy;

		File tempFile = createTempFileFor(destination);
		try {
			URLConnection connection = prepareURLConnection(url.toString(), timeout);
			connection.connect();
			checkHttpStatus(connection);

			long declaredSize = connection.getContentLengthLong();
			String eTag = connection.getHeaderField("ETag");

			// Only digest when we have something to compare against; hashing every
			// download would cost CPU for no benefit.
			Hash eTagHash = policy == ETagPolicy.IGNORE ? Hash.UNKNOWN : hashFromETag(eTag);
			if (policy == ETagPolicy.REQUIRE && eTagHash == Hash.UNKNOWN) {
				logger.warn("ETag [{}] of {} is not a usable hex digest; no checksum will be recorded.", eTag, url);
			}

			MessageDigest digest = eTagHash == Hash.UNKNOWN ? null : newDigest(eTagHash);
			long written;
			try (InputStream raw = connection.getInputStream();
					InputStream in = digest == null ? raw : new DigestInputStream(raw, digest)) {
				written = Files.copy(in, tempFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
			}

			if (declaredSize >= 0 && written != declaredSize) {
				throw new IOException(String.format(
						"Incomplete download of %s: got %d bytes but the server declared %d.",
						url, written, declaredSize));
			}

			String actualDigest = digest == null ? null : toHex(digest.digest());
			if (actualDigest != null && !actualDigest.equalsIgnoreCase(normalizeETag(eTag))) {
				throw new IOException(String.format(
						"Corrupt download of %s: %s of the content is %s but the server's ETag says %s.",
						url, eTagHash, actualDigest, normalizeETag(eTag)));
			}

			moveIntoPlace(tempFile, destination);

			// Sidecars are written only once the content is known good, so a failed
			// download can never leave validation metadata describing a file that is
			// not there.
			writeSizeFile(destination, written);
			if (hashURL != null) {
				if (hash == null || hash == Hash.UNKNOWN) {
					throw new IllegalArgumentException("Hash URL given but algorithm is unknown");
				}
				downloadFile(hashURL, hashFileFor(destination, hash));
			} else if (actualDigest != null) {
				writeHashFile(destination, eTagHash, actualDigest);
			}
		} finally {
			deleteQuietly(tempFile);
		}
	}

	/**
	 * Verifies that an HTTP connection returned a 2xx status. Connections using a
	 * non-HTTP protocol (<code>file:</code>, <code>ftp:</code>, ...) are left alone.
	 * <p>
	 * Without this check a 404 error page is written into the cache as though it
	 * were the requested file &mdash; and because the <code>.size</code> sidecar is
	 * then taken from that same error response, {@link #validateFile(File)} would
	 * subsequently declare it valid.
	 *
	 * @param connection an already connected {@link URLConnection}
	 * @throws HttpStatusException if the status is outside the 2xx range
	 * @throws IOException if the status could not be read
	 * @author Amr ALHOSSARY
	 * @since 7.3.0
	 */
	public static void checkHttpStatus(URLConnection connection) throws IOException {
		if (!(connection instanceof HttpURLConnection)) {
			return;
		}
		HttpURLConnection http = (HttpURLConnection) connection;
		int code = http.getResponseCode();
		if (code >= 200 && code < 300) {
			return;
		}
		if (code == 301 || code == 302 || code == 307 || code == 308) {
			// The JDK follows redirects automatically, but never across protocols, so
			// an http -> https redirect surfaces here and is worth naming explicitly.
			logger.warn("{} returned redirect {} to [{}], which was not followed "
					+ "(the JDK does not follow redirects that change protocol).",
					connection.getURL(), code, http.getHeaderField("Location"));
		}
		throw new HttpStatusException(code, connection.getURL().toString(), http.getResponseMessage());
	}

	/**
	 * Creates validation files beside a file to be downloaded.<br>
	 * Whenever possible, for a <code>file.ext</code> file, it creates 
	 * <code>file.ext.size</code> and <code>file.hash</code> for in the same 
	 * folder where <code>file.ext</code> exists.
	 * If the file connection size could not be deduced from the URL, no size file is created. 
	 * If <code>hashURL</code> is <code>null</code>, no hash file is created.
	 * @param url the remote file URL to download
	 * @param localDestination the local file to download into
	 * @param hashURL the URL of the hash file to download. Can be <code>null</code>.
	 * @param hash The Hashing algorithm. Ignored if <code>hashURL</code> is <code>null</code>.
	 */
	public static void createValidationFiles(URL url, File localDestination, URL hashURL, Hash hash){
		createValidationFiles(url, localDestination, hashURL, hash, ETagPolicy.USE_IF_HEX_DIGEST);
	}

	/**
	 * Creates validation files beside a file to be downloaded, with explicit control
	 * over how the <code>ETag</code> response header is treated.
	 *
	 * @param url the remote file URL to download
	 * @param localDestination the local file to download into
	 * @param hashURL the URL of the hash file to download. Can be <code>null</code>.
	 * @param hash The Hashing algorithm. Ignored if <code>hashURL</code> is <code>null</code>.
	 * @param eTagPolicy how to treat the <code>ETag</code> header when <code>hashURL</code>
	 *        is <code>null</code>. May be <code>null</code>, treated as {@link ETagPolicy#IGNORE}.
	 * @author Amr ALHOSSARY
	 * @since 7.3.0
	 */
	public static void createValidationFiles(URL url, File localDestination, URL hashURL, Hash hash,
			ETagPolicy eTagPolicy){
		try {
			URLConnection resourceConnection = url.openConnection();
			createValidationFiles(resourceConnection, localDestination, hashURL, hash, eTagPolicy);
		} catch (IOException e) {
			logger.warn("could not open connection to resource file due to exception: {}", e.getMessage());
		}
	}
	/**
	 * Creates validation files beside a file to be downloaded.<br>
	 * Whenever possible, for a <code>file.ext</code> file, it creates 
	 * <code>file.ext.size</code> and <code>file.hash_XXXX</code> in the same 
	 * folder where <code>file.ext</code> exists (XXXX may be DM5, SHA1, or SHA256).
	 * If the file connection size could not be deduced from the resourceUrlConnection 
	 * {@link URLConnection}, no size file is created. 
	 * If <code>hashURL</code> is <code>null</code>, no hash file is created.<br>
	 * <b>N.B.</b> None of the hashing algorithms is implemented (yet), because we did not need any of them yet.
	 * @param resourceUrlConnection the remote file URLConnection to download
	 * @param localDestination the local file to download into
	 * @param hashURL the URL of the hash file to download. Can be <code>null</code>.
	 * @param hash The Hashing algorithm. Ignored if <code>hashURL</code> is <code>null</code>.
	 * @since 7.0.0
	 */
	public static void createValidationFiles(URLConnection resourceUrlConnection, File localDestination, URL hashURL, Hash hash){
		createValidationFiles(resourceUrlConnection, localDestination, hashURL, hash, ETagPolicy.USE_IF_HEX_DIGEST);
	}

	/**
	 * Creates validation files beside a file to be downloaded, with explicit control
	 * over how the <code>ETag</code> response header is treated.
	 * <p>
	 * Nothing is written when the connection reports a non-2xx status: previously an
	 * error page's <code>Content-Length</code> would be recorded as the expected
	 * size, so the cached error page then passed validation.
	 *
	 * @param resourceUrlConnection the remote file URLConnection to download
	 * @param localDestination the local file to download into
	 * @param hashURL the URL of the hash file to download. Can be <code>null</code>.
	 * @param hash The Hashing algorithm. Ignored if <code>hashURL</code> is <code>null</code>.
	 * @param eTagPolicy how to treat the <code>ETag</code> header when <code>hashURL</code>
	 *        is <code>null</code>. May be <code>null</code>, treated as {@link ETagPolicy#IGNORE}.
	 * @author Amr ALHOSSARY
	 * @since 7.3.0
	 */
	public static void createValidationFiles(URLConnection resourceUrlConnection, File localDestination, URL hashURL,
			Hash hash, ETagPolicy eTagPolicy){
		try {
			checkHttpStatus(resourceUrlConnection);
		} catch (IOException e) {
			logger.warn("Not writing validation metadata for {}: {}", resourceUrlConnection.getURL(), e.getMessage());
			return;
		}

		long size = resourceUrlConnection.getContentLengthLong();
		if(size == -1) {
			logger.debug("Could not find expected file size for resource {}. Size validation metadata file won't be available for this download.", resourceUrlConnection.getURL());
		} else {
			logger.debug("Content-Length: {}", size);
			writeSizeFile(localDestination, size);
		}

		if(hashURL == null) {
			ETagPolicy policy = eTagPolicy == null ? ETagPolicy.IGNORE : eTagPolicy;
			if (policy != ETagPolicy.IGNORE) {
				String eTag = resourceUrlConnection.getHeaderField("ETag");
				Hash eTagHash = hashFromETag(eTag);
				if (eTagHash == Hash.UNKNOWN) {
					if (policy == ETagPolicy.REQUIRE) {
						logger.warn("ETag [{}] of {} is not a usable hex digest; no checksum recorded.",
								eTag, resourceUrlConnection.getURL());
					}
				} else {
					writeHashFile(localDestination, eTagHash, normalizeETag(eTag));
				}
			}
			return;
		}

		if(hash == null || hash == Hash.UNKNOWN)
			throw new IllegalArgumentException("Hash URL given but algorithm is unknown");
		try {
			downloadFile(hashURL, hashFileFor(localDestination, hash));
		} catch (IOException e) {
			logger.warn("Could not write validation hash file due to exception: {}", e.getMessage());
		}
	}

	/**
	 * Determines which hashing algorithm an <code>ETag</code> header value
	 * represents, based on the length of the hex digest it contains.
	 * <p>
	 * Only a value consisting solely of hex characters is considered. The
	 * <code>&lt;time&gt;-&lt;size&gt;</code> form used by nginx and Apache always
	 * contains a <code>-</code> and therefore never matches.
	 *
	 * @param eTagHeaderValue the raw header value, possibly quoted or weak-prefixed.
	 *        May be <code>null</code>.
	 * @return the matching algorithm, or {@link Hash#UNKNOWN} if the value is not a
	 *         hex digest of a recognised length
	 * @author Amr ALHOSSARY
	 * @since 7.3.0
	 */
	public static Hash hashFromETag(String eTagHeaderValue) {
		String value = normalizeETag(eTagHeaderValue);
		if (value == null || !value.matches("[0-9a-fA-F]+")) {
			return Hash.UNKNOWN;
		}
		switch (value.length()) {
			case 32: return Hash.MD5;
			case 40: return Hash.SHA1;
			case 64: return Hash.SHA256;
			default: return Hash.UNKNOWN;
		}
	}

	/**
	 * Strips the weak-validator prefix and surrounding quotes from an
	 * <code>ETag</code> header value.
	 *
	 * @param eTagHeaderValue the raw header value. May be <code>null</code>.
	 * @return the bare value, or <code>null</code> if the input was <code>null</code>
	 *         or blank
	 * @author Amr ALHOSSARY
	 * @since 7.3.0
	 */
	public static String normalizeETag(String eTagHeaderValue) {
		if (eTagHeaderValue == null) {
			return null;
		}
		String value = eTagHeaderValue.trim();
		if (value.startsWith("W/")) {
			value = value.substring(2).trim();
		}
		if (value.length() >= 2 && value.startsWith("\"") && value.endsWith("\"")) {
			value = value.substring(1, value.length() - 1);
		}
		return value.isEmpty() ? null : value;
	}

	/**
	 * Writes a <code>&lt;name&gt;.hash_&lt;ALGORITHM&gt;</code> sidecar file
	 * containing the given digest as bare lowercase hex.
	 *
	 * @param localDestination the file the digest describes
	 * @param hash the hashing algorithm
	 * @param hexDigest the digest, in hex
	 * @author Amr ALHOSSARY
	 * @since 7.3.0
	 */
	public static void writeHashFile(File localDestination, Hash hash, String hexDigest) {
		if (hash == null || hash == Hash.UNKNOWN || hexDigest == null) {
			return;
		}
		File hashFile = hashFileFor(localDestination, hash);
		try (PrintStream out = new PrintStream(hashFile, StandardCharsets.UTF_8.name())) {
			out.println(hexDigest.toLowerCase());
		} catch (IOException e) {
			logger.warn("Could not write validation hash file due to exception: {}", e.getMessage());
		}
	}

	/**
	 * Computes the digest of a file.
	 *
	 * @param file the file to digest
	 * @param hash the algorithm to use
	 * @return the digest as bare lowercase hex
	 * @throws IOException if the file could not be read
	 * @author Amr ALHOSSARY
	 * @since 7.3.0
	 */
	public static String computeHash(File file, Hash hash) throws IOException {
		try (InputStream in = new BufferedInputStream(Files.newInputStream(file.toPath()), DIGEST_BUFFER_SIZE)) {
			return computeHash(in, hash);
		}
	}

	/**
	 * Computes the digest of a stream. The stream is read to its end but not closed.
	 *
	 * @param in the stream to digest
	 * @param hash the algorithm to use
	 * @return the digest as bare lowercase hex
	 * @throws IOException if the stream could not be read
	 * @author Amr ALHOSSARY
	 * @since 7.3.0
	 */
	public static String computeHash(InputStream in, Hash hash) throws IOException {
		MessageDigest digest = newDigest(hash);
		byte[] buffer = new byte[DIGEST_BUFFER_SIZE];
		int read;
		while ((read = in.read(buffer)) != -1) {
			digest.update(buffer, 0, read);
		}
		return toHex(digest.digest());
	}

	/**
	 * Checks a file against an expected digest.
	 *
	 * @param file the file to check
	 * @param hash the algorithm to use
	 * @param expectedHex the expected digest, in hex; compared case-insensitively
	 * @return <code>true</code> if the digests match
	 * @throws IOException if the file could not be read
	 * @author Amr ALHOSSARY
	 * @since 7.3.0
	 */
	public static boolean verifyHash(File file, Hash hash, String expectedHex) throws IOException {
		return expectedHex != null && expectedHex.trim().equalsIgnoreCase(computeHash(file, hash));
	}

	/**
	 * The JDK name of a hashing algorithm, which differs from the enum constant for
	 * the SHA variants.
	 *
	 * @param hash the algorithm
	 * @return the name to pass to {@link MessageDigest#getInstance(String)}
	 * @author Amr ALHOSSARY
	 * @since 7.3.0
	 */
	public static String getAlgorithmName(Hash hash) {
		switch (hash) {
			case MD5:    return "MD5";
			case SHA1:   return "SHA-1";
			case SHA256: return "SHA-256";
			default: throw new IllegalArgumentException("Hashing algorithm not known: " + hash);
		}
	}

	/**
	 * Reads the expected digest out of a <code>.hash_XXXX</code> sidecar file.
	 * <p>
	 * The file may have been downloaded verbatim from a server, so several common
	 * layouts are accepted: a bare hex digest, the <code>md5sum</code> style
	 * <code>&lt;hex&gt;&nbsp;&nbsp;&lt;filename&gt;</code>, and the BSD style
	 * <code>MD5 (&lt;filename&gt;) = &lt;hex&gt;</code>.
	 *
	 * @param hashFile the sidecar file
	 * @return the digest in hex, or <code>null</code> if nothing recognisable was found
	 */
	static String parseHashFile(File hashFile) {
		try (Scanner scanner = new Scanner(hashFile, StandardCharsets.UTF_8.name())) {
			while (scanner.hasNextLine()) {
				String line = scanner.nextLine().trim();
				if (line.isEmpty()) {
					continue;
				}
				Matcher bare = BARE_HEX_HASH.matcher(line);
				if (bare.matches()) {
					return bare.group(1);
				}
				Matcher bsd = BSD_HASH.matcher(line);
				if (bsd.matches()) {
					return bsd.group(1);
				}
				return null; // first meaningful line was not a digest
			}
		} catch (IOException e) {
			logger.warn("Could not read hash file [{}]: {}", hashFile, e.getMessage());
		}
		return null;
	}

	
	/**
	 * Validate a local file based on pre-existing metadata files for size and hash.<br>
	 * If the passed in <code>localFile</code> parameter is a file named <code>file.ext</code>, the function searches in the same folder for:
	 * <ul>
	 * <li><code>file.ext.size</code>: If found, it compares the size stored in it to the length of <code>localFile</code> (in bytes).</li>
	 * <li><code>file.ext.hash_XXXX (where XXXX is DM5, SHA1, or SHA256)</code>: If found, it compares the size stored in it to the hash code of <code>localFile</code>.</li>
	 * </ul>
	 * If any of these comparisons fail, the function returns <code>false</code>. otherwise it returns true.
	 * <p>
	 * <b>N.B.</b> None of the 3 common verification hashing algorithms are implement yet.
	 * @param localFile The file to validate
	 * @return <code>false</code> if any of the size or hash code metadata files exists but its contents does not match the expected value in the file, <code>true</code> otherwise.
	 * @since 7.0.0
	 */
	public static boolean validateFile(File localFile) {
		// getParentFile() is null for a bare relative name such as new File("x.cif"),
		// which used to make this method throw a NullPointerException.
		File parent = localFile.getAbsoluteFile().getParentFile();
		if (parent == null) {
			logger.debug("Cannot determine the parent directory of [{}]; nothing to validate against.", localFile);
			return true;
		}

		File sizeFile = new File(parent, localFile.getName() + SIZE_EXT);
		if(sizeFile.exists()) {
            try (Scanner scanner = new Scanner(sizeFile)) {
                if (!scanner.hasNextLong()) {
                    // An empty or truncated .size file used to raise an unchecked
                    // NoSuchElementException that escaped the catch below.
                    logger.warn("Size metadata file [{}] is empty or malformed; skipping size validation.", sizeFile);
                } else {
                    long expectedSize = scanner.nextLong();
                    long actualSize = localFile.length();
                    if (expectedSize != actualSize) {
                        logger.warn("File [{}] size ({}) does not match expected size ({}).", localFile, actualSize, expectedSize);
                        return false;
                    }
                }
            } catch (FileNotFoundException e) {
                logger.warn("could not validate size of file [{}] because no size metadata file exists.", localFile);
            }
		}

		File[] hashFiles = parent.listFiles(new FilenameFilter() {
			final String hashPattern = String.format("%s%s_(%s|%s|%s)", Pattern.quote(localFile.getName()), HASH_EXT, Hash.MD5, Hash.SHA1, Hash.SHA256);
			@Override
			public boolean accept(File dir, String name) {
				return name.matches(hashPattern);
			}
		});
		// listFiles() returns null if the parent is not a directory or cannot be read.
		if (hashFiles == null || hashFiles.length == 0) {
			return true;
		}

		// Verify against every sidecar present, not only the first one found.
		for (File hashFile : hashFiles) {
			String name = hashFile.getName();
			String algo = name.substring(name.lastIndexOf('_') + 1);
			Hash hash;
			try {
				hash = Hash.valueOf(algo);
			} catch (IllegalArgumentException e) {
				throw new IllegalArgumentException("Hashing algorithm not known: " + algo, e);
			}
			if (hash == Hash.UNKNOWN) {
				throw new IllegalArgumentException("Hashing algorithm not known: " + algo);
			}

			String expected = parseHashFile(hashFile);
			if (expected == null) {
				// A sidecar we cannot read should not condemn an otherwise good download.
				logger.warn("Could not read a digest from [{}]; skipping {} validation of [{}].", hashFile, hash, localFile);
				continue;
			}
			try {
				if (!verifyHash(localFile, hash, expected)) {
					logger.warn("File [{}] {} does not match the expected digest {}.", localFile, hash, expected);
					return false;
				}
			} catch (IOException e) {
				logger.warn("Could not compute the {} of [{}]: {}", hash, localFile, e.getMessage());
				return false;
			}
		}

		return true;
	}

	/**
	 * The <code>&lt;name&gt;.hash_&lt;ALGORITHM&gt;</code> sidecar file for a
	 * downloaded file.
	 */
	private static File hashFileFor(File localDestination, Hash hash) {
		return new File(localDestination.getAbsoluteFile().getParentFile(),
				String.format("%s%s_%s", localDestination.getName(), HASH_EXT, hash));
	}

	/**
	 * Writes the <code>&lt;name&gt;.size</code> sidecar file.
	 */
	private static void writeSizeFile(File localDestination, long size) {
		File sizeFile = new File(localDestination.getAbsoluteFile().getParentFile(),
				localDestination.getName() + SIZE_EXT);
		try (PrintStream sizePrintStream = new PrintStream(sizeFile, StandardCharsets.UTF_8.name())) {
			sizePrintStream.print(size);
		} catch (IOException e) {
			logger.warn("Could not write size validation metadata file due to exception: {}", e.getMessage());
		}
	}

	private static MessageDigest newDigest(Hash hash) {
		try {
			return MessageDigest.getInstance(getAlgorithmName(hash));
		} catch (NoSuchAlgorithmException e) {
			// MD5, SHA-1 and SHA-256 are required of every Java platform.
			throw new IllegalStateException("Required hashing algorithm is unavailable: " + hash, e);
		}
	}

	private static String toHex(byte[] bytes) {
		StringBuilder sb = new StringBuilder(bytes.length * 2);
		for (byte b : bytes) {
			sb.append(Character.forDigit((b >> 4) & 0xF, 16));
			sb.append(Character.forDigit(b & 0xF, 16));
		}
		return sb.toString();
	}

	/**
	 * Creates a temp file whose name is derived from the destination.
	 * {@link Files#createTempFile} rejects prefixes shorter than 3 characters, so
	 * short names are padded.
	 */
	private static File createTempFileFor(File destination) throws IOException {
		String prefix = getFilePrefix(destination);
		while (prefix.length() < 3) {
			prefix = prefix + "_";
		}
		return Files.createTempFile(prefix, "." + getFileExtension(destination)).toFile();
	}

	private static void moveIntoPlace(File tempFile, File destination) throws IOException {
		try {
			Files.move(tempFile.toPath(), destination.toPath(),
					StandardCopyOption.REPLACE_EXISTING, StandardCopyOption.ATOMIC_MOVE);
		} catch (AtomicMoveNotSupportedException e) {
			// The temp directory is often on a different filesystem than the cache.
			Files.copy(tempFile.toPath(), destination.toPath(), StandardCopyOption.REPLACE_EXISTING);
		}
	}

	private static void deleteQuietly(File file) {
		if (file == null) {
			return;
		}
		try {
			Files.deleteIfExists(file.toPath());
		} catch (IOException e) {
			logger.debug("Could not delete temporary file [{}]: {}", file, e.getMessage());
		}
	}

	/**
	 * Converts path to Unix convention and adds a terminating slash if it was
	 * omitted. 
	 *
	 * @param path original platform dependent path
	 * @return path in Unix convention
	 * @author Peter Rose
	 * @since 3.2
	 */
	public static String toUnixPath(String path) {
		String uPath = path;
		if (uPath.contains("\\")) {
			uPath = uPath.replaceAll("\\\\", "/");
		}
		// this should be removed, it's need since "\" is added AtomCache code
		if (uPath.endsWith("//")) {
			uPath = uPath.substring(0, uPath.length() - 1);
		}
		if (!uPath.endsWith("/")) {
			uPath = uPath + "/";
		}
		return uPath;
	}

	/**
	 * Expands ~ in paths to the user's home directory.
	 *
	 * <p>
	 * This does not work for some special cases for paths: Other users' homes
	 * (~user/...), and Tilde expansion within the path (/.../~/...). In these cases
	 *  the original argument is returned.
	 *
	 * @param file A filepath starting with a tilde
	 * @return An absolute path
	 */
	public static String expandUserHome(String file) {
		// replace any / with the proper separator (/ or \ for Linux and Windows respectively).
		file = file.replaceAll("/", "\\"+File.separator); //The "\\" is to escape the separator if needed.
		if (file.startsWith("~") && (file.length() == 1 || File.separator.equals(file.substring(1, 2)))) {
			file = System.getProperty("user.home") + file.substring(1);
		}
		return file;
	}

	/**
	 * Pings a HTTP URL. This effectively sends a HEAD request and returns
	 * <code>true</code> if the response code is in the 200-399 range.
	 *
	 * @param url The HTTP URL to be pinged.
	 * @param timeout The timeout in millis for both the connection timeout and
	 * the response read timeout. Note that the total timeout is effectively two
	 * times the given timeout.
	 * @return <code>true</code> if the given HTTP URL has returned response
	 * code 200-399 on a HEAD request within the given timeout, otherwise
	 * <code>false</code>.
	 * @author BalusC,
	 * http://stackoverflow.com/questions/3584210/preferred-java-way-to-ping-a-http-url-for-availability
	 */
	public static boolean ping(String url, int timeout) {
		//url = url.replaceFirst("https", "http"); // Otherwise an exception may be thrown on invalid SSL certificates.

		try {
			HttpURLConnection connection = (HttpURLConnection) prepareURLConnection(url, timeout);
			connection.setRequestMethod("HEAD");
			int responseCode = connection.getResponseCode();
			return (200 <= responseCode && responseCode <= 399);
		} catch (IOException exception) {
			return false;
		}
	}

	/**
	 * Prepare {@link URLConnection} with customised timeouts.
	 *
	 * @param url The URL
	 * @param timeout The timeout in millis for both the connection timeout and
	 * the response read timeout. Note that the total timeout is effectively two
	 * times the given timeout.
	 *
	 * <p>
	 * Example of code.      <code>
		 * UrlConnection conn = prepareURLConnection("http://www.google.com/", 20000);
	 * conn.connect();
	 * conn.getInputStream();
	 * </code>
	 * <p>
	 *
	 * <strong>NB. User should execute connect() method before getting input
	 * stream.</strong>
	 * @return
	 * @throws IOException
	 * @author Jacek Grzebyta
	 */
	public static URLConnection prepareURLConnection(String url, int timeout) throws IOException {
		URLConnection connection = new URL(url).openConnection();
		connection.setReadTimeout(timeout);
		connection.setConnectTimeout(timeout);
		return connection;
	}

	/**
	 * Recursively delete a folder &amp; contents
	 *
	 * @param dir directory to delete
	 */
	public static void deleteDirectory(Path dir) throws IOException {
		if(dir == null || !Files.exists(dir))
			return;
		Files.walkFileTree(dir, new SimpleFileVisitor<>() {
	        @Override
	        public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
	            Files.delete(file);
	            return FileVisitResult.CONTINUE;
	        }

	        @Override
	        public FileVisitResult postVisitDirectory(Path dir, IOException e) throws IOException {
	            if (e != null) {
	                throw e;
	            }
	            Files.delete(dir);
	            return FileVisitResult.CONTINUE;
	        }
	    });
	}
	/**
	 * Recursively delete a folder &amp; contents
	 *
	 * @param dir directory to delete
	 */
	public static void deleteDirectory(String dir) throws IOException {
		deleteDirectory(Paths.get(dir));
	}

}
