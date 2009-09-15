 
/**
 *                    BioJava development code
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
 *
 */

package org.biojava.bio.seq.io.game;


/**
 * An interface that can be tested for by nested handlers
 * when trying to do a callback.
 * This one handles callbacks from &lt;name&gt; elements.
 *
 * @author David Huen
 * @since 1.8
 */
public interface GAMENameCallbackItf {

/**
 * Callback implemented by nesting class to allow
 * nested class to pass string to nesting class for
 * handling.
 */
  public void NameSetStringValue(String s);

}

