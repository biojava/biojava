package org.biojava.nbio.ws.uniprot;

import java.io.Reader;

/**
 * A general purpose content holder to wrap the 
 * network response consumption
 * @author pbansal
 *
 */
public interface ResponseConsumer {
	public void onSucess();
	public void onFail(int httpStatusCode);
	public void setReader(Reader reader);
} // ResponseConsumer
