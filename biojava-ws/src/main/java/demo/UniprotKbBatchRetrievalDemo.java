package demo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Reader;
import java.util.Collections;
import java.util.Iterator;
import java.util.Scanner;

import org.biojava.nbio.ws.uniprot.BatchRetrieve;
import org.biojava.nbio.ws.uniprot.BatchRetrieve.BatchRetrievalFormat;
import org.biojava.nbio.ws.uniprot.ResponseConsumer;
import org.biojava.nbio.ws.uniprot.UniprotUtilities;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * An example of how to use batch retrieval utility in uniprot
 * @author pbansal
 *
 */
public class UniprotKbBatchRetrievalDemo {

	private final static Logger logger = LoggerFactory.getLogger(UniprotKbBatchRetrievalDemo.class);
	
	public static void main(String[] args) {
		
		Iterable<String> itr = readFile(args[0]);
		BatchRetrieve batchRetrieveUniprotEntries = UniprotUtilities.getBatchUtility();
		
		try {
			batchRetrieveUniprotEntries.retrieve(itr, BatchRetrievalFormat.list, new ResponseConsumer() {
				private Reader reader = null;
				@Override
				public void setReader(Reader reader) {
					this.reader = reader;
				}
				
				@Override
				public void onSucess() {
					try {
						BufferedReader buf = new BufferedReader(reader);
						String line = buf.readLine();
						while(line != null) {
							logger.debug(line);
							line = buf.readLine();
						}
						buf.close();
					} catch(Exception e) {
						logger.error("Error while reading the input file.", e.getMessage());
					}
				}
				
				@Override
				public void onFail(int httpStatusCode) {
					logger.error("HttpStatus received " + httpStatusCode);
				}
			});
		} catch (IOException e) {
			logger.error("Error while retrieving entries", e.getMessage());
		}
	}
	
	static Iterable<String> readFile(final String fileName) {
		return new Iterable<String>() {
			@Override
			public Iterator<String> iterator() {
				File fl = new File(fileName);
				if (!fl.exists()) {
					throw new IllegalArgumentException("File doesn't exist: " + fileName);
				}
				try {
					final Scanner scanner = new Scanner(fl);
					return new Iterator<String>() {
						boolean pastEnd = false;
						@Override
						public boolean hasNext() {
							if (pastEnd) return !pastEnd;
							boolean hasNext = scanner.hasNextLine();
							if (!hasNext) {
								scanner.close();
								pastEnd = true;
							}
							return hasNext;
						}
						public String next() {
							return scanner.nextLine();
						}
					};
				} catch (FileNotFoundException e) {
					logger.error("Error while reading the input file.", e.getMessage());
				}
				return Collections.emptyIterator();
			}
		};
	}
}
