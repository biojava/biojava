package org.biojava.bio.program.unigene;

import java.net.URL;

import org.biojava.bio.BioException;

/**
 * <p>Objects that can be used to produce a <code>UnigeneDB</code> instance
 * given a URL.</p>
 *
 * <p><em>This class is intended for implementers of the Unigene
 * functionality</em></p>
 *
 * <p>The URL is used to locate the unigene data, as well as potentialy deciding
 * upon the access mechanism. Factory implementations are provided for
 * <code>file</code> and <code>jdbc</code> URLs. You can add another stoorage
 * mechanism to the system by writing an implementation of
 * <code>UnigeneFactory</code> and registering it with UnigeneTools.</p>
 *
 * @author Matthew Pocock
 */
public interface UnigeneFactory {
  public UnigeneDB loadUnigene(URL unigeneURL)
  throws BioException;
  
  public UnigeneDB createUnigene(URL unigeneURL)
  throws BioException;

  public boolean canAccept(URL unigeneURL);
}
