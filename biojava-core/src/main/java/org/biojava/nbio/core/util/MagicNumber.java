package org.biojava.nbio.core.util;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;

/**
 * The 'magic number' is a sequence of bytes that the beginning of a file that can be used to determine the
 * file type
 *
 * @since 7.0.3
 * @author Gary Murphy
 */
public class MagicNumber {

	/**
	 * The magic number of a gzip file is 0x1F8B. (ref: https://en.wikipedia.org/wiki/Gzip#:~:text=%22gzip%22%20is%20often%20also%20used,and%20the%20operating%20system%20ID.)
	 * @param path the path to the file
	 * @return <code>true</code> if the file has the gzip magic number
	 * @throws IOException if there is an error reading the start of the file
	 */
	public static boolean isGZIP(Path path) throws IOException {
		try (
				InputStream input = Files.newInputStream(path, StandardOpenOption.READ)
		) {
			byte[] magic = new byte[2];
			int count = input.read(magic);
			if (count != 2) {
				return false;
			}
			int id = (int)magic[0] & 0x00ff;
			id <<= 8;
			id += (int)magic[1] & 0x00ff;
			return (id == 0x1f8b);
		}
	}
}
