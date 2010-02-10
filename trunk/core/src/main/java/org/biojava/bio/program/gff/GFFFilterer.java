/*
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

package org.biojava.bio.program.gff;


/**
 * An object that filters a stream of GFF, forwarding some
 * <span class="type">GFFRecord</span>s to a
 * listening <span class="type">GFFDocumentHandler</span>, and dropping others.
 * <p>
 * The choice to forward or drop is made by a
 * <span class="type">GFFRecordFilter</span>.
 * <p>
 * Unless otherwise stated, all methods forward to the listening handler without
 * altering the arguments in any way.
 *
 * @author Matthew Pocock
 */
public class GFFFilterer implements GFFDocumentHandler {
  /**
   * The <span class="type">GFFDocumentHandler</span> that will recieve all
   * accepted entries.
   */
  private GFFDocumentHandler handler;
  /**
   * The <span class="type">GFFRecordFilter</span> that will decide what gets
   * forwarded.
   */
  private GFFRecordFilter filter;

  /**
   * Create a new <span class="type">GFFFilterer</span> that will forward
   * to <span class="arg">handler</span> everything that
   * <span class="arg">filter</span> accepts.
   *
   * @param handler the listening <span class="type">GFFDocumentHandler</span>
   * @param filter  the <span class="type">GFFRecordFilter</span> that decides
   *                what is forwarded to <span class="arg">handler</span>
   */
  public GFFFilterer(GFFDocumentHandler handler, GFFRecordFilter filter) {
    this.handler = handler;
    this.filter = filter;
  }
  
  public void startDocument(String locator) {
    handler.startDocument(locator);
  }
  
  public void endDocument() {
    handler.endDocument();
  }
  
  public void commentLine(String comment) {
    handler.commentLine(comment);
  }
  
  /**
   * Only forward the <span class="type">GFFRecord</span>s that match a filter.
   */
  public void recordLine(GFFRecord record) {
    if(filter.accept(record)) {
      handler.recordLine(record);
    }
  }
}
