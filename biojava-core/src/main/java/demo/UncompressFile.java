package demo;

import java.io.FileInputStream;
import java.io.InputStream;

import org.biojava.nbio.core.util.UncompressInputStream;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Uncompresses a single tarred or zipped file, writing output to stdandard out
 */
public class UncompressFile {
    private final static Logger logger
		= LoggerFactory.getLogger(UncompressFile.class);

    /**
	 * Reads a file, uncompresses it, and sends the result to stdout.
	 * Also writes trivial statistics to stderr.
	 * @param args An array with one String element, the name of the file to read.
	 * @throws IOException for any failure
	 */
	public static void main(String[] args) throws Exception {

		if (args.length != 1) {
			logger.info("Usage: UncompressInputStream <file>");
			System.exit(1);
		}
		long beg = System.currentTimeMillis();

		long tot;
		try (InputStream in = new FileInputStream(args[0]); 
         ) {
			tot = UncompressInputStream.uncompress(in, System.out);
		}

		long end = System.currentTimeMillis();
		logger.info("Decompressed {} bytes", tot);
		logger.info("Time: {} seconds", (end - beg) / 1000);
	}
    
}
