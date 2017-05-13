package org.biojava.nbio.ws.uniprot;

import java.io.IOException;

public interface BatchRetrieve {
	
	/**
	 * List of formats in which uniprot entry set can be 
	 * fetched.
	 * TODO provide tab delimited option too. Once list of columns are listed 
	 * @author pbansal
	 *
	 */
	public enum BatchRetrievalFormat {
		xml("xml"),
		rdf("rdf"),
		list("list"),
		text("text"),
		fastaCannonical("fasta"),
		fastaWithIsoforms("fasta"),
		excel("xls");
		
		private String extension = null;
		private BatchRetrievalFormat(String ext) {
			this.extension = ext;
		}
		/**
		 * @return
		 */
		public String getExtenstion() {
			return extension;
		}
	}
	
	public static final String INVALID_JOB_ID = "invalidJob";
	
	/**
	 * Given uniprotkb accessions, retrieves the ids from 
	 * www.uniprot.org in the said format and consume it using the NetworkResponseConsumer
	 * @param ids
	 * @throws IOException
	 */
	public void retrieve(Iterable<String> ids, BatchRetrievalFormat format, ResponseConsumer out) throws IOException;

} // BatchRetrieve
